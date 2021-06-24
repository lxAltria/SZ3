#include <SZ.hpp>
#include <def.hpp>
#include <cstdio>
#include <iostream>
#include <memory>

template<class T, uint N>
double SZ_Compress_by_config(int argc, char **argv, int argp, std::unique_ptr<T[]> &data, double eb, std::array<size_t, N> dims) {
    SZ::Config<double, N> conf(eb, dims);
    if (argp < argc) {
        int block_size = atoi(argv[argp++]);
        conf.block_size = block_size;
        conf.stride = block_size;
    }

    if (argp + 1 < argc) {
        int lorenzo_op = atoi(argv[argp++]);
        int regression_op = atoi(argv[argp++]);
        conf.enable_lorenzo = lorenzo_op == 1 || lorenzo_op == 3;
        conf.enable_2ndlorenzo = lorenzo_op == 2 || lorenzo_op == 3;
        conf.enable_regression = regression_op == 1 || regression_op == 3;
        conf.enable_2ndregression = regression_op == 2 || regression_op == 3;
    }

    if (argp < argc) {
        conf.quant_bin = atoi(argv[argp++]);
    }

    return SZ_Compress(data, conf);
}

int main(int argc, char **argv) {
    size_t num = 0;
    auto data = SZ::readfile<double>(argv[1], num);
    std::cout << argv[1] << " " << argv[2] << std::endl;
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2]);
    std::cout << dim << std::endl;
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }
    double reb = atof(argv[argp++]);
    double max = data[0];
    double min = data[0];
    for (int i = 1; i < num; i++) {
        if (max < data[i]) max = data[i];
        if (min > data[i]) min = data[i];
    }
    double eb = reb * (max - min);
    if (argp == argc) {
        if (dim == 1) {
            SZ::SZ_Compress(data, eb, dims[0]);
        } else if (dim == 2) {
            SZ::SZ_Compress(data, eb, dims[0], dims[1]);
        } else if (dim == 3) {
            SZ::SZ_Compress(data, eb, dims[0], dims[1], dims[2]);
        } else if (dim == 4) {
            SZ::SZ_Compress(data, eb, dims[0], dims[1], dims[2], dims[3]);
        }
        return 0;
    }

    if (dim == 1) {
        SZ_Compress_by_config<double, 1>(argc, argv, argp, data, eb, std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        SZ_Compress_by_config<double, 2>(argc, argv, argp, data, eb, std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        SZ_Compress_by_config<double, 3>(argc, argv, argp, data, eb, std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        SZ_Compress_by_config<double, 4>(argc, argv, argp, data, eb, std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }

    return 0;
}
