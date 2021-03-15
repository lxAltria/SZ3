#ifndef SZ_IE_COMPRESSOR_HPP
#define SZ_IE_COMPRESSOR_HPP

#include "compressor/Compressor.hpp"
#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/MemoryUtil.hpp"
#include "utils/Config.hpp"
#include "utils/FileUtil.h"
#include "def.hpp"
#include <climits>
#include <cstring>

namespace SZ {
    template<class T, uint N, class Predictor, class Quantizer, class Encoder, class Lossless>
    class SZIECompressor : public concepts::CompressorInterface<T> {
    public:
        SZIECompressor(Predictor predictor, Quantizer quantizer, Encoder encoder, Lossless lossless) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(quantizer.get_eb())),
                predictor(predictor), quantizer(quantizer), encoder(encoder), lossless(lossless) {
            static_assert(std::is_base_of<concepts::PredictorInterface<T, N>, Predictor>::value,
                          "must implement the predictor interface");
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
            compressed_data_dec_pos = compressed_data;
        }

        void init(size_t num_elements){
            quant_inds = std::vector<int>();
            compressed_data = new uchar[num_elements * sizeof(T)];
            // leave blank for encoder offset
            compressed_data_pos = compressed_data + sizeof(size_t);            
        }

        void decorrelate(T * data, const std::vector<size_t>& dims, size_t stride){
            std::cout << compressed_data_pos - compressed_data << std::endl;

            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(dims),
                                                                                         std::end(dims),
                                                                                         stride, 0);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(dims),
                                                                                         std::end(dims), 1,
                                                                                         0);
            int block_size = stride;
            std::array<size_t, N> intra_block_dims;
            predictor.precompress_data(inter_block_range->begin());
            quantizer.precompress_data();
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t cur_dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == cur_dims - 1 &&
                                               dims[i] - cur_index * stride < block_size) ?
                                              dims[i] - cur_index * stride : block_size;
                    }

                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());
                    concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                    if (!predictor.precompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    predictor_withfallback->precompress_block_commit();
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        // quant_inds[quant_count++] = quantizer.quantize_and_overwrite(
                        //         *element, predictor_withfallback->predict(element));
                        quant_inds.push_back(quantizer.quantize_and_overwrite(
                                *element, predictor_withfallback->predict(element)));
                    }
                }
            }

            predictor.postcompress_data(inter_block_range->begin());
            quantizer.postcompress_data();

            write(stride, compressed_data_pos);
            predictor.save(compressed_data_pos);
            quantizer.save(compressed_data_pos);
            predictor.clear();
            quantizer.clear();
        }

        uchar *compress(T *data, size_t& compressed_size) {

            size_t offset = compressed_data_pos - compressed_data;
            *reinterpret_cast<size_t*>(compressed_data) = offset;
            write(quant_inds.size(), compressed_data_pos);
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

        T* decompress(uchar const *lossless_compressed_data, size_t length){
            size_t remaining_length = length;
            compressed_data = lossless.decompress(lossless_compressed_data, remaining_length);
            compressed_data_dec_pos = compressed_data;
            size_t offset = 0;
            read(offset, compressed_data_dec_pos, remaining_length);
            uchar const * encoder_pos = compressed_data + offset;
            size_t total_elements = 0;
            read(total_elements, encoder_pos, remaining_length);
            encoder.load(encoder_pos, remaining_length);
            quant_inds.clear();
            quant_inds = encoder.decode(encoder_pos, total_elements);
            encoder.postprocess_decode();
            // lossless.postdecompress_data(compressed_data);
            quant_count = 0;
            return NULL;
        }

        T * recover(const std::vector<size_t>& dims) {
            std::cout << compressed_data_dec_pos - compressed_data << std::endl;
            size_t remaining_length = INT_MAX;
            size_t stride = 0;
            read(stride, compressed_data_dec_pos, remaining_length);
            int block_size = stride;
            predictor.load(compressed_data_dec_pos, remaining_length);
            quantizer.load(compressed_data_dec_pos, remaining_length);
            std::cout << "load done\n";
            fflush(stdout);

            std::array<size_t, N> intra_block_dims;
            size_t num_elements = 1;
            for(int i=0; i<dims.size(); i++){
                num_elements *= dims[i];
            }
            auto dec_data = std::make_unique<T[]>(num_elements);
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.get(),
                                                                                         std::begin(dims),
                                                                                         std::end(dims),
                                                                                         stride,
                                                                                         0);

            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.get(),
                                                                                         std::begin(dims),
                                                                                         std::end(dims), 1,
                                                                                         0);
            predictor.predecompress_data(inter_block_range->begin());
            quantizer.predecompress_data();
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; block++) {
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t cur_dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == cur_dims - 1) ? dims[i] - cur_index * block_size
                                                                      : block_size;
                    }
                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());

                    concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                    if (!predictor.predecompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        *element = quantizer.recover(predictor_withfallback->predict(element), quant_inds[quant_count ++]);
                    }
                }
            }
            predictor.postdecompress_data(inter_block_range->begin());
            quantizer.postdecompress_data();
            predictor.clear();
            quantizer.clear();
            return dec_data.release();
        }

        void finalize(){
            if(compressed_data){
                lossless.postdecompress_data(compressed_data);
            }
        }

    private:
        Predictor predictor;
        LorenzoPredictor<T, N, 1> fallback_predictor;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        std::vector<int> quant_inds;
        int quant_count = 0;
        uchar * compressed_data = NULL;
        uchar * compressed_data_pos = NULL;
        uchar const * compressed_data_dec_pos = NULL;
    };

    template<class T, uint N, class Predictor, class Quantizer, class Encoder, class Lossless>
    SZIECompressor<T, N, Predictor, Quantizer, Encoder, Lossless>
    make_sz_interleaved_encoding_compressor(const Config<T, N> &conf, Predictor predictor, Quantizer quantizer, Encoder encoder,
                               Lossless lossless) {
        return SZIECompressor<T, N, Predictor, Quantizer, Encoder, Lossless>(predictor, quantizer, encoder, lossless);
    }
}
#endif
