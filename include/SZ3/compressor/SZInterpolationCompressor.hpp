#ifndef _SZ_INTERPOLATION_COMPRESSSOR_HPP
#define _SZ_INTERPOLATION_COMPRESSSOR_HPP

#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"
#include <cstring>
#include <cmath>

namespace SZ {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZInterpolationCompressor {
    public:


        SZInterpolationCompressor(Quantizer quantizer, Encoder encoder, Lossless lossless) :
                quantizer(quantizer), encoder(encoder), lossless(lossless) {

            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) {
            size_t remaining_length = cmpSize;
            uchar *buffer = lossless.decompress(cmpData, remaining_length);
            uchar const *buffer_pos = buffer;

            read(global_dimensions.data(), N, buffer_pos, remaining_length);
            read(blocksize, buffer_pos, remaining_length);
            read(interpolator_id, buffer_pos, remaining_length);
            read(direction_sequence_id, buffer_pos, remaining_length);

            // read auxilliary data
            read(detection_block_size, buffer_pos);
            size_t num_detection_block = 0;
            read(num_detection_block, buffer_pos);
            flushed_block_id = std::vector<uchar>(num_detection_block);
            convertByteArray2IntArray_fast_1b_sz(num_detection_block, buffer_pos, (num_detection_block - 1)/8 + 1, flushed_block_id.data());
            significant_block_id = std::vector<uchar>(num_detection_block);
            convertByteArray2IntArray_fast_1b_sz(num_detection_block, buffer_pos, (num_detection_block - 1)/8 + 1, significant_block_id.data());

            srand(3333);
            init();
            auto num_flushed_elements = compute_auxilliary_data_decompress(decData); 

            quantizer.load(buffer_pos, remaining_length);
            encoder.load(buffer_pos, remaining_length);
            quant_inds = encoder.decode(buffer_pos, num_elements - num_flushed_elements);

            encoder.postprocess_decode();

            lossless.postdecompress_data(buffer);
            double eb = quantizer.get_eb();
            double eb_final = eb / pow(c, interpolation_level - 1);

            std::cout << "start decompression\n";
            // *decData = quantizer.recover(0, quant_inds[quant_index++]);
            recover(0, *decData, 0);

            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
                // if (level >= 3) {
                //     quantizer.set_eb(eb * eb_ratio);
                // } else {
                //     quantizer.set_eb(eb);
                // }
                current_level = level;
                current_base_eb = eb_final;
                quantizer.set_eb(eb_final);
                eb_final *= c;

                size_t stride = 1U << (level - 1);
                auto inter_block_range = std::make_shared<
                        SZ::multi_dimensional_range<T, N>>(decData,
                                                           std::begin(global_dimensions), std::end(global_dimensions),
                                                           stride * blocksize, 0);
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {
                    auto end_idx = block.get_global_index();
                    for (int i = 0; i < N; i++) {
                        end_idx[i] += stride * blocksize;
                        if (end_idx[i] > global_dimensions[i] - 1) {
                            end_idx[i] = global_dimensions[i] - 1;
                        }
                    }
                    block_interpolation(decData, block.get_global_index(), end_idx, PB_recover,
                                        interpolators[interpolator_id], direction_sequence_id, stride);

                }
            }
            quantizer.postdecompress_data();
//            timer.stop("Interpolation Decompress");

            return decData;
        }

        // compress given the error bound
        uchar *compress(const Config &conf, T *data, size_t &compressed_size) {
            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());
            blocksize = conf.interpBlockSize;
            interpolator_id = conf.interpAlgo;
            direction_sequence_id = conf.interpDirection;

            srand(3333);
            init();
            compute_auxilliary_data(conf, data);

            quant_inds.reserve(num_elements);
            size_t interp_compressed_size = 0;

            double eb = quantizer.get_eb();
            double eb_final = eb / pow(c, interpolation_level - 1);

            // quant_inds.push_back(quantizer.quantize_and_overwrite(*data, 0));
            quantize(0, *data, 0);

            Timer timer;
            timer.start();

            for (uint level = interpolation_level; level > 0 && level <= interpolation_level; level--) {
                // if (level >= 3) {
                //     quantizer.set_eb(eb * eb_ratio);
                // } else {
                //     quantizer.set_eb(eb);
                // }
                current_level = level;
                current_base_eb = eb_final;
                quantizer.set_eb(eb_final);
                eb_final *= c;

                size_t stride = 1U << (level - 1);

                auto inter_block_range = std::make_shared<
                        SZ::multi_dimensional_range<T, N>>(data, std::begin(global_dimensions),
                                                           std::end(global_dimensions),
                                                           blocksize * stride, 0);

                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();

                for (auto block = inter_begin; block != inter_end; ++block) {
                    auto end_idx = block.get_global_index();
                    for (int i = 0; i < N; i++) {
                        end_idx[i] += blocksize * stride;
                        if (end_idx[i] > global_dimensions[i] - 1) {
                            end_idx[i] = global_dimensions[i] - 1;
                        }
                    }

                    block_interpolation(data, block.get_global_index(), end_idx, PB_predict_overwrite,
                                        interpolators[interpolator_id], direction_sequence_id, stride);
                }
            }
            assert(quant_inds.size() == num_elements);
//            timer.stop("Prediction & Quantization");

//            writefile("pred.dat", preds.data(), num_elements);
//            writefile("quant.dat", quant_inds.data(), num_elements);
            encoder.preprocess_encode(quant_inds, 0);
            size_t bufferSize = 1.2 * (quantizer.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size());

            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;

            write(global_dimensions.data(), N, buffer_pos);
            write(blocksize, buffer_pos);
            write(interpolator_id, buffer_pos);
            write(direction_sequence_id, buffer_pos);

            // add auxilliary array
            write(detection_block_size, buffer_pos);
            size_t num_detection_block = flushed_block_id.size();
            write(num_detection_block, buffer_pos);
            convertIntArray2ByteArray_fast_1b_to_result_sz(flushed_block_id.data(), flushed_block_id.size(), buffer_pos);
            convertIntArray2ByteArray_fast_1b_to_result_sz(significant_block_id.data(), significant_block_id.size(), buffer_pos);

            quantizer.save(buffer_pos);
            quantizer.postcompress_data();


            timer.start();
            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();
//            timer.stop("Coding");
            assert(buffer_pos - buffer < bufferSize);

            timer.start();
            uchar *lossless_data = lossless.compress(buffer,
                                                     buffer_pos - buffer,
                                                     compressed_size);
            lossless.postcompress_data(buffer);
//            timer.stop("Lossless");

            compressed_size += interp_compressed_size;
            return lossless_data;
        }

    private:

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover
        };

        void init() {
            assert(blocksize % 2 == 0 && "Interpolation block size should be even numbers");
            num_elements = 1;
            interpolation_level = -1;
            for (int i = 0; i < N; i++) {
                if (interpolation_level < ceil(log2(global_dimensions[i]))) {
                    interpolation_level = (uint) ceil(log2(global_dimensions[i]));
                }
                num_elements *= global_dimensions[i];
            }

            dimension_offsets[N - 1] = 1;
            for (int i = N - 2; i >= 0; i--) {
                dimension_offsets[i] = dimension_offsets[i + 1] * global_dimensions[i + 1];
            }

            dimension_sequences = std::vector<std::array<int, N>>();
            auto sequence = std::array<int, N>();
            for (int i = 0; i < N; i++) {
                sequence[i] = i;
            }
            do {
                dimension_sequences.push_back(sequence);
            } while (std::next_permutation(sequence.begin(), sequence.end()));
        }

        inline void quantize(size_t idx, T &d, T pred) {
            if(flushed_block[idx]) {
                d = 0;
            }
            else{
                auto default_eb = quantizer.get_eb();
                if((current_level == 1) && significant_block[idx]){
                    quantizer.set_eb(current_base_eb * c3);
                }
                T noise = 2.0*rand()/RAND_MAX - 1.0;
                if(fabs(pred) > 1e-2) pred += noise * 0.3 * default_eb;
                quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
                quantizer.set_eb(default_eb);        
            }

        }

        inline void recover(size_t idx, T &d, T pred) {
            if(flushed_block[idx]){
                d = 0;
            }
            else{
                auto default_eb = quantizer.get_eb();
                if((current_level == 1) && significant_block[idx]){
                    quantizer.set_eb(current_base_eb * c3);
                }
                T noise = 2.0*rand()/RAND_MAX - 1.0;
                if(fabs(pred) > 1e-2) pred += noise * 0.3 * default_eb;
                d = quantizer.recover(pred, quant_inds[quant_index++]);
                quantizer.set_eb(default_eb); 
            }
        };


        double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride,
                                      const std::string &interp_func,
                                      const PredictorBehavior pb) {
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }
            double predict_error = 0;

            size_t stride3x = 3 * stride;
            size_t stride5x = 5 * stride;
            if (interp_func == "linear" || n < 5) {
                if (pb == PB_predict_overwrite) {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        quantize(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        quantize(d - data, *d, *(d - stride));
                        // if (n < 4) {
                        //     quantize(d - data, *d, *(d - stride));
                        // } else {
                        //     quantize(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                        // }
                    }
                } else {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        recover(d - data, *d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        recover(d - data, *d, *(d - stride));
                        // if (n < 4) {
                        //     recover(d - data, *d, *(d - stride));
                        // } else {
                        //     recover(d - data, *d, interp_linear1(*(d - stride3x), *(d - stride)));
                        // }
                    }
                }
            } else {
                if (pb == PB_predict_overwrite) {

                    T *d;
                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        quantize(d - data, *d,
                                 interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    d = data + begin + stride;
                    quantize(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    d = data + begin + i * stride;
                    quantize(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        quantize(d - data, *d, *(d - stride));
                        // quantize(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    }

                } else {
                    T *d;

                    size_t i;
                    for (i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        recover(d - data, *d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    d = data + begin + stride;

                    recover(d - data, *d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    d = data + begin + i * stride;
                    recover(d - data, *d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));

                    if (n % 2 == 0) {
                        d = data + begin + (n - 1) * stride;
                        recover(d - data, *d, *(d - stride));
                        // recover(d - data, *d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    }
                }
            }

            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 1, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, size_t stride = 1) {
            return block_interpolation_1d(data, begin[0], end[0], stride, interp_func, pb);
        }

        template<uint NN = N>
        typename std::enable_if<NN == 2, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, size_t stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            const std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[0]] - begin[dims[0]]) * dimension_offsets[dims[0]],
                                                        stride * dimension_offsets[dims[0]], interp_func, pb);
            }
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]];
                predict_error += block_interpolation_1d(data, begin_offset,
                                                        begin_offset +
                                                        (end[dims[1]] - begin[dims[1]]) * dimension_offsets[dims[1]],
                                                        stride * dimension_offsets[dims[1]], interp_func, pb);
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 3, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, size_t stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;

            auto default_eb = quantizer.get_eb();
            quantizer.set_eb(default_eb * c2);

            const std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    size_t begin_offset = begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[0]] - begin[dims[0]]) *
                                                            dimension_offsets[dims[0]],
                                                            stride * dimension_offsets[dims[0]], interp_func, pb);
                }
            }

            quantizer.set_eb(default_eb * c1);

            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                          k * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[1]] - begin[dims[1]]) *
                                                            dimension_offsets[dims[1]],
                                                            stride * dimension_offsets[dims[1]], interp_func, pb);
                }
            }

            quantizer.set_eb(default_eb );

            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                          begin[dims[2]] * dimension_offsets[dims[2]];
                    predict_error += block_interpolation_1d(data, begin_offset,
                                                            begin_offset +
                                                            (end[dims[2]] - begin[dims[2]]) *
                                                            dimension_offsets[dims[2]],
                                                            stride * dimension_offsets[dims[2]], interp_func, pb);
                }
            }
            return predict_error;
        }


        template<uint NN = N>
        typename std::enable_if<NN == 4, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, size_t stride = 1) {
            double predict_error = 0;
            size_t stride2x = stride * 2;
            max_error = 0;
            const std::array<int, N> dims = dimension_sequences[direction];
            for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride2x : 0); j <= end[dims[1]]; j += stride2x) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                begin[dims[0]] * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[0]] - begin[dims[0]]) *
                                                                dimension_offsets[dims[0]],
                                                                stride * dimension_offsets[dims[0]], interp_func, pb);
                    }
                }
            }
            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride2x : 0); k <= end[dims[2]]; k += stride2x) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + begin[dims[1]] * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[1]] - begin[dims[1]]) *
                                                                dimension_offsets[dims[1]],
                                                                stride * dimension_offsets[dims[1]], interp_func, pb);
                    }
                }
            }
            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t t = (begin[dims[3]] ? begin[dims[3]] + stride2x : 0);
                         t <= end[dims[3]]; t += stride2x) {
                        size_t begin_offset = i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                              begin[dims[2]] * dimension_offsets[dims[2]] +
                                              t * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[2]] - begin[dims[2]]) *
                                                                dimension_offsets[dims[2]],
                                                                stride * dimension_offsets[dims[2]], interp_func, pb);
                    }
                }
            }

            max_error = 0;
            for (size_t i = (begin[dims[0]] ? begin[dims[0]] + stride : 0); i <= end[dims[0]]; i += stride) {
                for (size_t j = (begin[dims[1]] ? begin[dims[1]] + stride : 0); j <= end[dims[1]]; j += stride) {
                    for (size_t k = (begin[dims[2]] ? begin[dims[2]] + stride : 0); k <= end[dims[2]]; k += stride) {
                        size_t begin_offset =
                                i * dimension_offsets[dims[0]] + j * dimension_offsets[dims[1]] +
                                k * dimension_offsets[dims[2]] +
                                begin[dims[3]] * dimension_offsets[dims[3]];
                        predict_error += block_interpolation_1d(data, begin_offset,
                                                                begin_offset +
                                                                (end[dims[3]] - begin[dims[3]]) *
                                                                dimension_offsets[dims[3]],
                                                                stride * dimension_offsets[dims[3]], interp_func, pb);
                    }
                }
            }
            return predict_error;
        }

        void compute_auxilliary_data(const Config &conf, T *data){
            // TODO: add corner case when dimensions are not divisible by block size
            int block_size = detection_block_size;
            const auto& dims = conf.dims;
            int nx = dims[0] / block_size;
            int ny = dims[1] / block_size;
            int nz = dims[2] / block_size;
            flushed_block = std::vector<uchar>(num_elements, 0);
            flushed_block_id = std::vector<uchar>(nx * ny * nz, 0);
            Timer timer;
            timer.start();
            // compute statistics
            int block_id = 0;
            int flushed_count = 0;
            // using variance for significant blocks
            std::vector<double> var;
            T * x_data_pos = data;
            for(int i=0; i<nx; i++){
                T * y_data_pos = x_data_pos;
                for(int j=0; j<ny; j++){
                    T * z_data_pos = y_data_pos;
                    for(int k=0; k<nz; k++){
                        T max_abs_val = 0;
                        // compute average
                        double sum = 0;
                        T * xx_data_pos = z_data_pos;
                        for(int ii=0; ii<block_size; ii++){
                            T * yy_data_pos = xx_data_pos;
                            for(int jj=0; jj<block_size; jj++){
                                T * zz_data_pos = yy_data_pos;
                                for(int kk=0; kk<block_size; kk++){
                                    sum += *zz_data_pos;
                                    if(fabs(*zz_data_pos) > max_abs_val) max_abs_val = fabs(*zz_data_pos);
                                    zz_data_pos ++;
                                }
                                yy_data_pos += dims[2];
                            }
                            xx_data_pos += dims[1] * dims[2];
                        }
                        if(max_abs_val < quantizer.get_eb()*0.1){
                            T * xx_data_pos = z_data_pos;
                            for(int ii=0; ii<block_size; ii++){
                                T * yy_data_pos = xx_data_pos;
                                for(int jj=0; jj<block_size; jj++){
                                    T * zz_data_pos = yy_data_pos;
                                    for(int kk=0; kk<block_size; kk++){
                                        *zz_data_pos = 0;
                                        flushed_block[zz_data_pos - data] = 1;
                                        zz_data_pos ++;
                                    }
                                    yy_data_pos += dims[2];
                                }
                                xx_data_pos += dims[1] * dims[2];
                            }
                            flushed_block_id[block_id] = 1;
                            flushed_count ++;
                            // skip the variance
                            var.push_back(0);
                        }
                        else{
                            // compute variance 
                            double average = sum / (block_size*block_size*block_size);
                            double variance = 0;
                            T * xx_data_pos = z_data_pos;
                            for(int ii=0; ii<block_size; ii++){
                                T * yy_data_pos = xx_data_pos;
                                for(int jj=0; jj<block_size; jj++){
                                    T * zz_data_pos = yy_data_pos;
                                    for(int kk=0; kk<block_size; kk++){
                                        variance += (*zz_data_pos - average) * (*zz_data_pos - average);
                                        zz_data_pos ++;
                                    }
                                    yy_data_pos += dims[2];
                                }
                                xx_data_pos += dims[1] * dims[2];
                            }
                            variance /= (block_size*block_size*block_size);
                            var.push_back(variance);
                        }
                        block_id ++;
                        z_data_pos += block_size;
                    }
                    y_data_pos += block_size * dims[2];
                }
                x_data_pos += block_size * dims[1] * dims[2];
            }
            std::cout << "flushed_count = " << flushed_count << ", percent = " << flushed_count * 1.0 / (nx * ny * nz) << std::endl;
            auto var_tmp(var);
            std::sort(var_tmp.begin(), var_tmp.end());
            double percent = 0.9;
            double threshold = var_tmp[(int)(percent*var.size())];
            std::cout << percent * 100 << "% threshold = " << threshold << std::endl;
            significant_block = std::vector<uchar>(num_elements, 0);
            significant_block_id = std::vector<uchar>(nx * ny * nz, 0);
            block_id = 0;
            x_data_pos = data;
            for(int i=0; i<nx; i++){
                T * y_data_pos = x_data_pos;
                for(int j=0; j<ny; j++){
                    T * z_data_pos = y_data_pos;
                    for(int k=0; k<nz; k++){
                        if(var[block_id] > threshold){
                            T * xx_data_pos = z_data_pos;
                            for(int ii=0; ii<block_size; ii++){
                                T * yy_data_pos = xx_data_pos;
                                for(int jj=0; jj<block_size; jj++){
                                    T * zz_data_pos = yy_data_pos;
                                    for(int kk=0; kk<block_size; kk++){
                                        significant_block[zz_data_pos - data] = 1;
                                        zz_data_pos ++;
                                    }
                                    yy_data_pos += dims[2];
                                }
                                xx_data_pos += dims[1] * dims[2];
                            }
                            significant_block_id[block_id] = 1;
                        }
                        block_id ++;
                        z_data_pos += block_size;
                    }
                    y_data_pos += block_size * dims[2];
                }
                x_data_pos += block_size * dims[1] * dims[2];
            }
            timer.stop("detect");
        }

        size_t compute_auxilliary_data_decompress(const T *data){
            // TODO: add corner case when dimensions are not divisible by block size
            int block_size = detection_block_size;
            const auto& dims = global_dimensions;
            int nx = dims[0] / block_size;
            int ny = dims[1] / block_size;
            int nz = dims[2] / block_size;
            flushed_block = std::vector<uchar>(num_elements, 0);
            significant_block = std::vector<uchar>(num_elements, 0);
            int block_id = 0;
            size_t num_flushed_elements = 0;
            const T * x_data_pos = data;
            for(int i=0; i<nx; i++){
                const T * y_data_pos = x_data_pos;
                for(int j=0; j<ny; j++){
                    const T * z_data_pos = y_data_pos;
                    for(int k=0; k<nz; k++){
                        if(flushed_block_id[block_id]){
                            num_flushed_elements += block_size * block_size * block_size;
                            const T * xx_data_pos = z_data_pos;
                            for(int ii=0; ii<block_size; ii++){
                                const T * yy_data_pos = xx_data_pos;
                                for(int jj=0; jj<block_size; jj++){
                                    const T * zz_data_pos = yy_data_pos;
                                    for(int kk=0; kk<block_size; kk++){
                                        flushed_block[zz_data_pos - data] = 1;
                                        zz_data_pos ++;
                                    }
                                    yy_data_pos += dims[2];
                                }
                                xx_data_pos += dims[1] * dims[2];
                            }                        
                        } 
                        else if(significant_block_id[block_id]){
                            const T * xx_data_pos = z_data_pos;
                            for(int ii=0; ii<block_size; ii++){
                                const T * yy_data_pos = xx_data_pos;
                                for(int jj=0; jj<block_size; jj++){
                                    const T * zz_data_pos = yy_data_pos;
                                    for(int kk=0; kk<block_size; kk++){
                                        significant_block[zz_data_pos - data] = 1;
                                        zz_data_pos ++;
                                    }
                                    yy_data_pos += dims[2];
                                }
                                xx_data_pos += dims[1] * dims[2];
                            }                        
                        }
                        block_id ++;
                        z_data_pos += block_size;
                    }
                    y_data_pos += block_size * dims[2];
                }
                x_data_pos += block_size * dims[1] * dims[2];
            }
            return num_flushed_elements;
        }

        // *** modified from Sheng's code start ***
        void convertIntArray2ByteArray_fast_1b_to_result_sz(const uchar* intArray, size_t intArrayLength, uchar *& compressed_pos){
            size_t byteLength = 0;
            size_t i, j; 
            if(intArrayLength%8==0)
                byteLength = intArrayLength/8;
            else
                byteLength = intArrayLength/8+1;
                
            size_t n = 0;
            int tmp, type;
            for(i = 0;i<byteLength;i++){
                tmp = 0;
                for(j = 0;j<8&&n<intArrayLength;j++){
                    type = intArray[n];
                    if(type == 1)
                        tmp = (tmp | (1 << (7-j)));
                    n++;
                }
                *(compressed_pos++) = (uchar)tmp;
            }
        }

        void convertByteArray2IntArray_fast_1b_sz(size_t intArrayLength, const uchar*& compressed_pos, size_t byteArrayLength, uchar * intArray){
            if(intArrayLength > byteArrayLength*8){
                printf("Error: intArrayLength > byteArrayLength*8\n");
                printf("intArrayLength=%zu, byteArrayLength = %zu", intArrayLength, byteArrayLength);
                exit(0);
            }
            size_t n = 0, i;
            int tmp;
            for (i = 0; i < byteArrayLength-1; i++) {
                tmp = *(compressed_pos++);
                intArray[n++] = (tmp & 0x80) >> 7;
                intArray[n++] = (tmp & 0x40) >> 6;
                intArray[n++] = (tmp & 0x20) >> 5;
                intArray[n++] = (tmp & 0x10) >> 4;
                intArray[n++] = (tmp & 0x08) >> 3;
                intArray[n++] = (tmp & 0x04) >> 2;
                intArray[n++] = (tmp & 0x02) >> 1;
                intArray[n++] = (tmp & 0x01) >> 0;      
            }
            tmp = *(compressed_pos++);  
            for(int i=0; n < intArrayLength; n++, i++){
                intArray[n] = (tmp & (1 << (7 - i))) >> (7 - i);
            }       
        }
        // *** modified from Sheng's code end ***

        int interpolation_level = -1;
        uint blocksize;
        int interpolator_id;
        double eb_ratio = 0.5;
        std::vector<std::string> interpolators = {"linear", "cubic"};
        std::vector<int> quant_inds;
        size_t quant_index = 0; // for decompress
        double max_error;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
        std::array<size_t, N> dimension_offsets;
        std::vector<std::array<int, N>> dimension_sequences;
        int direction_sequence_id;
        // added for artifact mitigation
        int current_level = 0;
        double c = sqrt(4.4159889);
        double c1 = 1.0 / sqrt(1.640625);
        double c2 = 1.0 / 1.640625;
        double c3 = 1.0 / sqrt(4.4159889);
        int detection_block_size = 16;
        std::vector<uchar> flushed_block;
        std::vector<uchar> flushed_block_id;
        std::vector<uchar> significant_block;
        std::vector<uchar> significant_block_id;
        double current_base_eb;
    };


};


#endif

