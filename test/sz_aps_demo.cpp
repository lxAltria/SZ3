#include <SZ.hpp>
#include <def.hpp>
#include <cstdio>
#include <iostream>
#include <memory>

template<typename T, class Predictor, uint N>
float SZ_APS_Compress(std::unique_ptr<T[]> const &data,
                  SZ::Config<T, N> conf, Predictor predictor) {
    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "eb = " << conf.eb << ", n_dims = " << N << ", quant_bin = " << conf.quant_bin << std::endl;

    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
    std::cout << "****************** Compression ******************" << std::endl;

    auto sz = SZ::make_sz_aps_compressor(conf, predictor, SZ::PastriQuantizer<T>(conf.eb, conf.quant_bin),
                                                 SZ::HuffmanEncoder<int>(),//SZ::FastEncoderType2<int>(), 
                                                 SZ::Lossless_zstd());

    size_t compressed_size = 0;
    std::unique_ptr<SZ::uchar[]> compressed;
    compressed.reset(sz.compress(data.get(), compressed_size));

    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Compression time = "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
              << "s"
              << std::endl;

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    SZ::writefile("compressed.dat", compressed.get(), compressed_size);

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = SZ::readfile<SZ::uchar>("compressed.dat", compressed_size);

    clock_gettime(CLOCK_REALTIME, &start);
    std::unique_ptr<T[]> dec_data;
    dec_data.reset(sz.decompress(compressed.get(), compressed_size));
    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "Decompression time: "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
              << "s"
              << std::endl;

    SZ::verify<T>(data_.data(), dec_data.get(), conf.num);
    SZ::writefile("compressed.dat.out", dec_data.get(), conf.num);
    return ratio;
}

int main(int argc, char **argv) {
    size_t num = 0;
    auto data = SZ::readfile<double>(argv[1], num);
    std::cout << "Read " << num << " elements\n";
    double eb = atof(argv[2]);
    int n_dims = atoi(argv[3]);
    std::vector<size_t> dims(n_dims);
    for(int i=0; i<n_dims; i++){
        dims[i] = atoi(argv[4 + i]);
    }

    using T = double;
    if(eb > 0.5){
        SZ::Config<T, 3> conf(eb, std::array<size_t, 3>{dims[0], dims[1], dims[2]});
        conf.quant_bin = atoi(argv[4 + n_dims]);
        std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, 3>>> predictors;
        predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, 3, 1>>(conf.eb));
        predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, 3>>(conf.block_size, conf.eb));
        SZ_APS_Compress(data, conf, SZ::ComposedPredictor<T, 3>(predictors));
    }
    else{
        SZ::Config<T, 1> conf(eb, std::array<size_t, 1>{dims[0] * dims[1] * dims[2]});
        conf.quant_bin = 2;
        conf.transpose = true;
        auto predictor = SZ::LorenzoPredictor<T, 1, 1>(conf.eb);
        SZ_APS_Compress(data, conf, predictor);
    }
    return 0;
}
