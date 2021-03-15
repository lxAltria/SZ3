#include "compressor/SZGeneralCompressor.hpp"
#include "compressor/SZIECompressor.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "predictor/RegressionPredictor.hpp"
#include "predictor/ComposedPredictor.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "encoder/HuffmanEncoder.hpp"
#include "lossless/Lossless_zstd.hpp"
#include "utils/FileUtil.h"
#include "utils/Config.hpp"
#include "utils/Verification.hpp"
#include "utils/Timer.hpp"
#include "def.hpp"
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>

template<class T, uint N>
unsigned char * SZ_Compress(std::vector<T>& data, const std::vector<size_t>& dims, const SZ::Config<T, N> &conf, size_t& compressed_size){
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
    predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
    predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
    auto sz = SZ::make_sz_interleaved_encoding_compressor(conf, SZ::ComposedPredictor<T, N>(predictors), SZ::LinearQuantizer<T>(conf.eb, conf.quant_bin),
                                             SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());

    size_t num_elements = conf.num;
    sz.init(num_elements);
    std::vector<T> data_copy(data);
    sz.decorrelate(data_copy.data(), dims, conf.block_size);
    auto compressed = sz.compress(NULL, compressed_size);
    std::cout << "Ratio = " << conf.num * sizeof(float) * 1.0 / compressed_size << std::endl;
    return compressed;
}

template<class T, uint N>
void SZ_Decompress(std::vector<T>& data, const std::vector<size_t>& dims, const SZ::Config<T, N> &conf, const unsigned char* compressed, size_t compressed_size){
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
    predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
    predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
    auto sz = SZ::make_sz_interleaved_encoding_compressor(conf, SZ::ComposedPredictor<T, N>(predictors), SZ::LinearQuantizer<T>(conf.eb, conf.quant_bin),
                                             SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());

    size_t num_elements = conf.num;
    sz.decompress(compressed, compressed_size);
    auto dec_data = sz.recover(dims);
    sz.finalize();

    SZ::verify<T>(data.data(), dec_data, conf.num);
    free(dec_data);
    std::cout << "Ratio = " << conf.num * sizeof(float) * 1.0 / compressed_size << std::endl;
}

template<class T>
std::vector<T> extract_region(T const * data, const std::vector<size_t>& dims, const std::vector<size_t>& offset, const std::vector<size_t>& extract_dims){
    std::vector<T> extracted_data(extract_dims[0] * extract_dims[1] * extract_dims[2]);
    for(int i=0; i<extract_dims[0]; i++){
        for(int j=0; j<extract_dims[1]; j++){
            for(int k=0; k<extract_dims[2]; k++){
                extracted_data[i * extract_dims[1] * extract_dims[2] + j * extract_dims[2] + k]
                    = data[(i + offset[0]) * dims[1] * dims[2] + (j + offset[1]) * dims[2] + (k + offset[2])];
            }
        }
    }
    return extracted_data;
}

size_t compute_size(const std::vector<size_t>& dims){
    size_t size = 1;
    for(int i=0; i<dims.size(); i++){
        size *= dims[i];
    }
    return size;
}

template<class T, uint N>
unsigned char * SZ_Compress_multiple_blocks(std::vector<T>& data, const std::vector<size_t>& dims, const SZ::Config<T, N> &conf, size_t& compressed_size){
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
    predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
    predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
    auto sz = SZ::make_sz_interleaved_encoding_compressor(conf, SZ::ComposedPredictor<T, N>(predictors), SZ::LinearQuantizer<T>(conf.eb, conf.quant_bin),
                                             SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());

    std::vector<size_t> dim1{32, 32, 32};
    std::vector<size_t> dim2{48, 48, 48};
    std::vector<size_t> dim3{64, 64, 64};
    std::vector<size_t> offset1{10, 20, 30};
    std::vector<size_t> offset2{50, 240, 300};
    std::vector<size_t> offset3{24, 100, 100};
    auto data1 = extract_region(data.data(), dims, offset1, dim1);
    auto data2 = extract_region(data.data(), dims, offset2, dim2);
    auto data3 = extract_region(data.data(), dims, offset3, dim3);

    size_t num_elements = compute_size(dim1) + compute_size(dim2) + compute_size(dim3);
    std::cout << compute_size(dim1) << " " << compute_size(dim2) << " " << compute_size(dim3) << std::endl;

    sz.init(num_elements);
    sz.decorrelate(data1.data(), dim1, conf.block_size);
    sz.decorrelate(data2.data(), dim2, conf.block_size);
    sz.decorrelate(data3.data(), dim3, conf.block_size);

    auto compressed = sz.compress(NULL, compressed_size);
    std::cout << "Ratio = " << num_elements * sizeof(float) * 1.0 / compressed_size << std::endl;
    return compressed;
}

template<class T, uint N>
void SZ_Decompress_multiple_blocks(std::vector<T>& data, const std::vector<size_t>& dims, const SZ::Config<T, N> &conf, const unsigned char* compressed, size_t compressed_size){
    std::vector<std::shared_ptr<SZ::concepts::PredictorInterface<T, N>>> predictors;
    predictors.push_back(std::make_shared<SZ::LorenzoPredictor<T, N, 1>>(conf.eb));
    predictors.push_back(std::make_shared<SZ::RegressionPredictor<T, N>>(conf.block_size, conf.eb));
    auto sz = SZ::make_sz_interleaved_encoding_compressor(conf, SZ::ComposedPredictor<T, N>(predictors), SZ::LinearQuantizer<T>(conf.eb, conf.quant_bin),
                                             SZ::HuffmanEncoder<int>(), SZ::Lossless_zstd());

    std::vector<size_t> dim1{32, 32, 32};
    std::vector<size_t> dim2{48, 48, 48};
    std::vector<size_t> dim3{64, 64, 64};
    std::vector<size_t> offset1{10, 20, 30};
    std::vector<size_t> offset2{50, 240, 300};
    std::vector<size_t> offset3{24, 100, 100};
    auto data1 = extract_region(data.data(), dims, offset1, dim1);
    auto data2 = extract_region(data.data(), dims, offset2, dim2);
    auto data3 = extract_region(data.data(), dims, offset3, dim3);

    sz.decompress(compressed, compressed_size);

    auto dec_data1 = sz.recover(dim1);
    SZ::verify<T>(data1.data(), dec_data1, compute_size(dim1));

    auto dec_data2 = sz.recover(dim2);
    SZ::verify<T>(data2.data(), dec_data2, compute_size(dim2));

    auto dec_data3 = sz.recover(dim3);
    SZ::verify<T>(data3.data(), dec_data3, compute_size(dim3));

    sz.finalize();
    free(dec_data1);
    free(dec_data2);
    free(dec_data3);

}

int main(int argc, char **argv) {

    size_t num = 0;
    auto data_pointer = SZ::readfile<float>(argv[1], num);
    // src_file_name = argv[1];
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2]);
    assert(2 <= dim && dim <= 3);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }
    unsigned char * compressed = NULL;
    size_t compressed_size = 0;
    if(dim == 2){
        SZ::Config<float, 2> conf(1e-2, std::array<size_t, 2>{dims[0], dims[1]});
        std::vector<float> data = std::vector<float>(data_pointer.get(), data_pointer.get() + conf.num);
        compressed = SZ_Compress(data, dims, conf, compressed_size);
        SZ_Decompress(data, dims, conf, compressed, compressed_size);
    }
    else if(dim == 3){    
        SZ::Config<float, 3> conf(1e-2, std::array<size_t, 3>{dims[0], dims[1], dims[2]});
        std::vector<float> data = std::vector<float>(data_pointer.get(), data_pointer.get() + conf.num);
        // compressed = SZ_Compress(data, dims, conf, compressed_size);
        // SZ_Decompress(data, dims, conf, compressed, compressed_size);
        // free(compressed);
        std::cout << std::endl;
        compressed = SZ_Compress_multiple_blocks(data, dims, conf, compressed_size); 
        SZ_Decompress_multiple_blocks(data, dims, conf, compressed, compressed_size);
    }
    free(compressed);        
    return 0;
}