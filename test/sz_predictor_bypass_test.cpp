#include <compressor/SZ_Predictor_Bypass_Compressor.hpp>
#include <quantizer/IntegerQuantizer.hpp>
#include <predictor/ComposedPredictor.hpp>
#include <predictor/LorenzoPredictor.hpp>
#include <predictor/RegressionPredictor.hpp>
#include <lossless/Lossless_zstd.hpp>
#include <utils/Iterator.hpp>
#include <SZ.hpp>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <memory>
#include <type_traits>

template<typename Type>
std::unique_ptr<Type[]> readfile(const char *file, size_t &num) {
    std::ifstream fin(file, std::ios::binary);
    if (!fin) {
        std::cout << " Error, Couldn't find the file" << "\n";
        return 0;
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    auto data = std::make_unique<Type[]>(num_elements);
    fin.read(reinterpret_cast<char *>(&data[0]), num_elements * sizeof(Type));
    fin.close();
    num = num_elements;
    return data;
}

template<typename Type>
void writefile(const char *file, Type *data, size_t num_elements) {
    std::ofstream fout(file, std::ios::binary);
    fout.write(reinterpret_cast<const char *>(&data[0]), num_elements * sizeof(Type));
    fout.close();
}

template<typename T, uint N>
float compress(const std::unique_ptr<T[]> &data, const std::unique_ptr<T[]> &pred, T eb, std::array<size_t, N> dims) {
    SZ::Config<float, N> conf(eb, dims);
    auto sz = SZ::make_SZ_Predictor_Bypass_Compressor(
            conf,
            SZ::LinearQuantizer<float>(eb),
            SZ::HuffmanEncoder<int>(),
            SZ::Lossless_zstd()
    );
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<unsigned char[]> compressed;
    size_t compressed_size = 0;
    compressed.reset(sz.compress(data.get(), pred.get(), compressed_size));
    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Compression time = "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000 << "s"
              << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;
    float ratio = 1.0 * conf.num * sizeof(float) / compressed_size;
    std::cout << "Compression ratio = " << ratio << std::endl;
    return ratio;

}


int main(int argc, char **argv) {

    std::cout << "Usage .\\predictor_bypass_test Uf48.dat -3 500 100 100 1e-2 Uf48_pred.dat" << std::endl;

    size_t num = 0;
    auto data = SZ::readfile<float>(argv[1], num);
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

    float max = data[0];
    float min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }

    float reb = atof(argv[argp++]);

    auto pred = SZ::readfile<float>(argv[argp], num);

    if (dim == 1) {
        compress<float, 1>(data, pred, reb * (max - min), std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        compress<float, 2>(data, pred, reb * (max - min), std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        compress<float, 3>(data, pred, reb * (max - min), std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        compress<float, 4>(data, pred, reb * (max - min), std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }

    return 0;
}
