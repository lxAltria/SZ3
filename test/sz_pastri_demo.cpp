#include <SZ.hpp>
#include <def.hpp>
#include <cstdio>
#include <iostream>
#include <memory>

template<typename T>
float SZ_Pastri_Compress(std::unique_ptr<T[]> const &data,
                  const SZ::Config<T, 1> &conf) {

    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "number of patterns = " << conf.num_patterns
                << ", pattern repeated times = " << conf.pattern_repeated_times
                << ", pattern size = " << conf.pattern_size
                << ", eb = " << conf.eb << std::endl;

    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
    std::cout << "****************** Compression ******************" << std::endl;

    auto sz = SZ::make_sz_pastri_compressor(conf, SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());

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
    return ratio;
}

int main(int argc, char **argv) {
    size_t num = 0;
    auto data = SZ::readfile<double>(argv[1], num);
    std::cout << "Read " << num << " elements\n";
    std::array<size_t, 1> dims;
    dims[0] = num;
    double eb = atof(argv[4]);
    double l1 = atof(argv[5]);
    double l2 = atof(argv[6]);
    SZ::Config<double, 1> conf(eb, dims, l1, l2);
    conf.pattern_repeated_times = atoi(argv[2]);
    conf.pattern_size = atoi(argv[3]);
    conf.num_patterns = num / conf.pattern_repeated_times / conf.pattern_size;
    assert(num == conf.num_patterns * conf.pattern_repeated_times * conf.pattern_size);
    SZ_Pastri_Compress(data, conf);
    return 0;
}
