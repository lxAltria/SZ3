#ifndef _SZ_SZ_PREDICTOR_BYPASS_COMPRESSOR_HPP
#define _SZ_SZ_PREDICTOR_BYPASS_COMPRESSOR_HPP

#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/Config.hpp"
#include "utils/fileUtil.h"
#include "def.hpp"
#include <cstring>

namespace SZ {
    template<class T, size_t N, class Quantizer, class Encoder, class Lossless>
    class SZ_Predictor_Bypassl_Compressor {
    public:


        SZ_Predictor_Bypassl_Compressor(const Config<T, N> &conf,
                                        Quantizer quantizer, Encoder encoder, Lossless lossless) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                block_size(conf.block_size), stride(conf.stride),
                global_dimensions(conf.dims), num_elements(conf.num) {
            static_assert(std::is_base_of_v<concepts::QuantizerInterface<T>, Quantizer>, "must implement the quatizer interface");
            static_assert(std::is_base_of_v<concepts::EncoderInterface<int>, Encoder>, "must implement the encoder interface");
            static_assert(std::is_base_of_v<concepts::LosslessInterface, Lossless>, "must implement the lossless interface");
        }

//        uchar *compress_demo(T *data, size_t &compressed_size) {
//            for (auto block = inter_block_range->begin(); block != inter_block_range->end(); ++block) {
//                intra_block_range->update(block);
//                predictor.precompress_block(intra_block_range);
//                for (auto element = intra_block_range->begin(); element != intra_block_range->end(); ++element) {
//                    quant_inds[quant_count++] = quantizer.quantize_and_overwrite(*element, predictor.predict(element));
//                }
//            }
//            ...
//        }

        // compress given the error bound
        uchar *compress(T *data, T *prediction, size_t &compressed_size) {
            // TODO: new quantizer if eb does not match

            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), stride, 0);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1, 0);
            std::array<size_t, N> intra_block_dims;
            std::vector<int> quant_inds(num_elements);
            quantizer.precompress_data();
            size_t quant_count = 0;
            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    // std::cout << *block << " " << lp.predict(block) << std::endl;
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1 && global_dimensions[i] - cur_index * stride < block_size) ?
                                              global_dimensions[i] - cur_index * stride : block_size;
                    }

                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());
                    {
                        auto intra_begin = intra_block_range->begin();
                        auto intra_end = intra_block_range->end();
                        for (auto element = intra_begin; element != intra_end; ++element) {
                            quant_inds[quant_count++] = quantizer.quantize_and_overwrite(
                                    *element, prediction[element.get_offset()]);
                        }
                    }
                }
            }

            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = "
                      << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                      << "s" << std::endl;

            quantizer.postcompress_data();
            uchar *compressed_data;
            if (stride > block_size) {
                std::cout << "Sampling Compress Mode is ON" << std::endl;
                quant_inds.resize(quant_count);
                compressed_data = new uchar[3 * quant_count * sizeof(T)];
            } else {
                compressed_data = new uchar[2 * num_elements * sizeof(T)];
            }

            uchar *compressed_data_pos = compressed_data;
            write(global_dimensions.data(), N, compressed_data_pos);
            write(block_size, compressed_data_pos);
            quantizer.save(compressed_data_pos);

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }

    private:
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
    };

    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    SZ_Predictor_Bypassl_Compressor<T, N, Quantizer, Encoder, Lossless>
    make_SZ_Predictor_Bypass_Compressor(const Config<T, N> &conf, Quantizer quantizer, Encoder encoder,
                                         Lossless lossless) {
        return SZ_Predictor_Bypassl_Compressor<T, N, Quantizer, Encoder, Lossless>(conf, quantizer, encoder,
                                                                                   lossless);
    }
}
#endif

