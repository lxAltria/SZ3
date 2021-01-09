#ifndef _SZ_SZ_PASTRI_HPP
#define _SZ_SZ_PASTRI_HPP

#include "predictor/Predictor.hpp"
#include "predictor/PatternPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/Config.hpp"
#include "utils/FileUtil.h"
#include "def.hpp"
#include <cstring>

namespace SZ {
    template<class T, size_t N, class Encoder, class Lossless>
    class GAMESS_Pattern_Based_Compressor {
    public:
        using Quantizer = PastriQuantizer<T>; 
        GAMESS_Pattern_Based_Compressor(const Config<T, N> &conf, Encoder encoder, Lossless lossless) :
                encoder(encoder), lossless(lossless), num_elements(conf.num), 
                eb(conf.eb), quant_bin(conf.quant_bin), pattern_eb(conf.pattern_eb), scale_eb(conf.scale_eb), 
                num_patterns(conf.num_patterns), pattern_repeated_times(conf.pattern_repeated_times), 
                pattern_size(conf.pattern_size){
            // assert if pattern dimensions and data dimensions do not match
            for(int i=0; i<N; i++){            
                assert(num_patterns * pattern_repeated_times * pattern_size == num_elements);
            }
            static_assert(std::is_base_of_v<concepts::QuantizerInterface<T>, Quantizer>, "must implement the quatizer interface");
            static_assert(std::is_base_of_v<concepts::EncoderInterface<int>, Encoder>, "must implement the encoder interface");
            static_assert(std::is_base_of_v<concepts::LosslessInterface, Lossless>, "must implement the lossless interface");
        }

        uchar *compress(T *data, size_t &compressed_size) {
            PatternPredictor<T> predictor = PatternPredictor<T>(pattern_repeated_times, pattern_size);
            Quantizer quantizer(eb, quant_bin);
            // quantizer.precompress_data();
            uchar *compressed_data = new uchar[2 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            write(num_patterns, compressed_data_pos);
            write(pattern_repeated_times, compressed_data_pos);
            write(pattern_size, compressed_data_pos);
            T eb = quantizer.get_eb();
            int radius = quantizer.get_radius();
            write(eb, compressed_data_pos);
            write(radius, compressed_data_pos);
            write(pattern_eb, compressed_data_pos);
            write(scale_eb, compressed_data_pos);
            Quantizer pattern_quantizer(pattern_eb, radius);
            Quantizer scale_quantizer(scale_eb, radius);
            std::vector<int> all_inds(num_elements + num_patterns * pattern_size + num_patterns * pattern_repeated_times);
            int * pattern_quant_inds = all_inds.data();
            int * scale_quant_inds = pattern_quant_inds + num_patterns * pattern_size;
            int * quant_inds = scale_quant_inds + num_patterns * pattern_repeated_times;
            size_t pattern_quant_count = 0;
            size_t scale_quant_count = 0;
            size_t quant_count = 0;
            // leave space for recording compressed predictor size
            uchar * compressed_predictor_pos = compressed_data_pos;
            compressed_data_pos += sizeof(size_t);
            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);
            const T * data_pos = data;
            // deal it in 1D
            int eb_exp = 0;
            frexp(eb, &eb_exp);
            std::cout << "pattern_eb = " << pattern_eb << std::endl;
            std::cout << "scale_eb = " << scale_eb << std::endl;
            std::vector<T> pattern(pattern_size);
            for(int i=0; i<num_patterns; i++){
                // extract pattern
                predictor.extract_pattern(data_pos, pattern.data());
                // quantize pattern by pattern quantizer
                for(int j=0; j<pattern.size(); j++){
                    pattern_quant_inds[pattern_quant_count ++] = pattern_quantizer.quantize_and_overwrite(pattern[j], 0);
                }
                for(int j=0; j<pattern_repeated_times; j++){
                    // compute scaling factor
                    T scale = predictor.compute_scale(data_pos, pattern.data());
                    // quantize scaling factor by scale quantizer
                    scale_quant_inds[scale_quant_count ++] = scale_quantizer.quantize_and_overwrite(scale, 0);
                    for(int k=0; k<pattern_size; k++){
                        T cur_data = *(data_pos++);
                        quant_inds[quant_count] = quantizer.quantize_and_overwrite(cur_data, scale * pattern[k]);
                        quant_count ++;
                    }
                }
                quantizer.save(compressed_data_pos);
                // std::cout << "block encoding size = " << compressed_data_pos - prev_size << std::endl;
            }
            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = "
                      << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                      << "s" << std::endl;
            size_t compressed_predictor_size = compressed_data_pos - compressed_predictor_pos - sizeof(size_t);
            write(compressed_predictor_size, compressed_predictor_pos);

            // quantizer.postcompress_data();
            // std::cout << "quantizer pos = " << compressed_data_pos - compressed_data << std::endl;
            pattern_quantizer.save(compressed_data_pos);
            scale_quantizer.save(compressed_data_pos);

            // std::cout << "encoder pos = " << compressed_data_pos - compressed_data << std::endl;
            auto encode_pos = compressed_data_pos;
            encoder.preprocess_encode(all_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(all_inds, compressed_data_pos);
            encoder.postprocess_encode();

            // std::cout << "before lossless size = " << compressed_data_pos - compressed_data << std::endl;
            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            // std::cout << "lossless size = " << compressed_size << std::endl;
            return lossless_data;
        }


        T *decompress(uchar const *lossless_compressed_data, const size_t length) {
            std::cout << "Pastri decompress, compressed size = " << length << std::endl;
            auto compressed_data = lossless.decompress(lossless_compressed_data, length);
            uchar const *compressed_data_pos = compressed_data;
            size_t remaining_length = length;
            read(num_patterns, compressed_data_pos, remaining_length);
            read(pattern_repeated_times, compressed_data_pos, remaining_length);
            read(pattern_size, compressed_data_pos, remaining_length);
            T eb = 0;
            int radius = 0;
            read(eb, compressed_data_pos, remaining_length);
            read(radius, compressed_data_pos, remaining_length);
            read(pattern_eb, compressed_data_pos, remaining_length);
            read(scale_eb, compressed_data_pos, remaining_length);
            PatternPredictor<T> predictor = PatternPredictor<T>(pattern_repeated_times, pattern_size);
            size_t compressed_predictor_size = 0;
            read(compressed_predictor_size, compressed_data_pos, remaining_length);
            Quantizer quantizer(eb, radius);
            Quantizer pattern_quantizer(pattern_eb, radius);
            Quantizer scale_quantizer(scale_eb, radius);
            // std::cout << "predictor pos = " << compressed_data_pos - compressed_data << std::endl;
            // std::cout << num_patterns << " " << pattern_repeated_times << " " << pattern_size << std::endl;
            uchar const *compressed_predictor_pos = compressed_data_pos;
            compressed_data_pos += compressed_predictor_size;

            // std::cout << "quantizer pos = " << compressed_data_pos - compressed_data << std::endl;
            pattern_quantizer.load(compressed_data_pos, remaining_length);
            scale_quantizer.load(compressed_data_pos, remaining_length);
            // std::cout << "encoder pos = " << compressed_data_pos - compressed_data << std::endl;
            encoder.load(compressed_data_pos, remaining_length);
            // std::cout << num_elements << " loading finished\n";
            // fflush(stdout);

            auto all_inds = encoder.decode(compressed_data_pos, num_elements + num_patterns * pattern_size + num_patterns * pattern_repeated_times);

            int const *pattern_inds = (int const *) all_inds.data();
            int const *scale_inds = pattern_inds + num_patterns * pattern_size;
            int const *quant_inds = scale_inds + num_patterns * pattern_repeated_times;
            encoder.postprocess_decode();
            lossless.postdecompress_data(compressed_data);
            // std::cout << "before lossless size = " << compressed_data_pos - compressed_data << std::endl;

            auto dec_data = std::make_unique<T[]>(num_elements);
            // quantizer.predecompress_data();
            // std::cout << "start decompression" << std::endl;
            size_t pattern_quant_count = 0;
            size_t scale_quant_count = 0;
            size_t quant_count = 0;
            std::vector<T> pattern(pattern_size);
            for(int i=0; i<num_patterns; i++){
                // std::cout << i << " " << compressed_predictor_pos - compressed_data << std::endl;
                // quantizer.load(compressed_predictor_pos, remaining_length);
                quantizer.load(compressed_predictor_pos, remaining_length);
                // recover pattern
                for(int j=0; j<pattern_size; j++){
                    pattern[j] = pattern_quantizer.recover(0, pattern_inds[pattern_quant_count ++]);
                }
                for(int j=0; j<pattern_repeated_times; j++){
                    T scale = scale_quantizer.recover(0, scale_inds[scale_quant_count ++]);
                    for(int k=0; k<pattern_size; k++){
                        dec_data[quant_count] = quantizer.recover(scale * pattern[k], quant_inds[quant_count]);
                        quant_count ++;
                    }
                }
            }
            // std::cout << pattern_quant_count << " " << scale_quant_count << " " << quant_count << std::endl;
            // quantizer.postdecompress_data();
            return dec_data.release();
        }


    private:
        Encoder encoder;
        Lossless lossless;
        size_t num_elements;
        // number of different patterns in the data
        size_t num_patterns;
        // number of repeated times for a pattern
        size_t pattern_repeated_times;
        // number of data points in a pattern
        size_t pattern_size;
        int quant_bin;
        T eb;
        T pattern_eb;
        T scale_eb;
    };

    template<class T, uint N, class Encoder, class Lossless>
    GAMESS_Pattern_Based_Compressor<T, N, Encoder, Lossless>
    make_sz_pastri_compressor(const Config<T, N> &conf, Encoder encoder,
                               Lossless lossless) {
        // Deal with 1D GAMESS data only in this implementation
        // Multidimensional data can be processed with linearization
        assert(N == 1);
        return GAMESS_Pattern_Based_Compressor<T, N, Encoder, Lossless>(conf, encoder, lossless);
    }
}
#endif

