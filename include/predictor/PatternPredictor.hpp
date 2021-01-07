#ifndef _SZ_PATTERN_PREDICTOR_HPP
#define _SZ_PATTERN_PREDICTOR_HPP

#include <cassert>
#include "def.hpp"
#include "utils/Iterator.hpp"
#include "quantizer/IntegerQuantizer.hpp"

namespace SZ {

    // Pattern predictor based on Pastri, always 1D *multimendision data should be linearized to 1D before using this predictor*
    // Normalize all the windows to (-1, 1) and solving X that minimizes (X - Y')^2
    // Predict the rest windowes by solving a, b in (a X - Y)^2
    // X: pattern data, Y: current data
    template<class T>
    class PatternPredictor : public concepts::PredictorInterface<T, 1> {
    public:
        static const uint8_t predictor_id = 0b00000011;
        using Range = multi_dimensional_range<T, 1>;
        using iterator = typename multi_dimensional_range<T, 1>::iterator;

        PatternPredictor(int pattern_repeated_times, int pattern_size) 
        : pattern_repeated_times(pattern_repeated_times), pattern_size(pattern_size) {}

        void precompress_data(const iterator &) const {}

        void postcompress_data(const iterator &) const {}

        void predecompress_data(const iterator &) const {}

        void postdecompress_data(const iterator &) const {}

        void extract_pattern(const T * data, T * pattern){
            // directly using max is near optimal
            const T * data_pos = data;
            // identify the max block
            int id = 0;
            T max_val = 0;
            for(int i=0; i<pattern_repeated_times; i++){
                for(int j=0; j<pattern_size; j++){
                    if(max_val < fabs(*data_pos)){
                        max_val = fabs(*data_pos);
                        id = i;
                    }
                    data_pos ++;
                }
            }
            ptrdiff_t offset = id * pattern_size;
            for(int i=0; i<pattern_size; i++){
                pattern[i] = (fabs(data[offset + i]) < 1e-10) ? 0 : data[offset + i];
            }
            // data_pos = data;
            // ptrdiff_t offset = id * pattern_size;
            // int count = 0;
            // compute_max_pattern_id(data + offset);
            // for(int i=0; i<pattern_repeated_times; i++){
            //     T scale = compute_pattern_scale(data_pos, data + offset);
            //     if(scale){
            //         for(int j=0; j<pattern_size; j++){
            //             T y = (fabs(data_pos[j]) < 1e-10) ? 0 : data_pos[j];
            //             patterns[j] += y / scale;
            //         }
            //         count ++;
            //     }
            //     // std::cout << i << " " << scale << std::endl;
            //     data_pos += pattern_size;
            // }
            // if(count){
            //     for(int j=0; j<pattern_size; j++){
            //         patterns[j] /= count;
            //     }
            // }
            // double err1 = 0, err2 = 0;
            // data_pos = data;
            // for(int i=0; i<pattern_repeated_times; i++){
            //     err1 += compute_pattern_error(data_pos, patterns.data());
            //     err2 += compute_pattern_error(data_pos, data + offset);
            //     data_pos += pattern_size;
            // }
            // if(err1 > err2){
            //     for(int i=0; i<pattern_size; i++){
            //         patterns[i] = (fabs(data[offset + i]) < 1e-10) ? 0 : data[offset + i];
            //     }
            // }
        }
        // compute slope and intercept for current window
        T compute_scale(const T * data, const T * pattern){
            // if(!patterns[max_pattern_id]) return 0;
            // return data[max_pattern_id] / patterns[max_pattern_id];
            T sigma_x2 = 0; // sum of x^2
            T sigma_xy = 0; // sum of xy
            for(int i=0; i<pattern_size; i++){
                T x = (fabs(pattern[i]) < 1e-10) ? 0 : pattern[i];
                T y = (fabs(data[i]) < 1e-10) ? 0 : data[i];
                sigma_x2 += x*x;
                sigma_xy += x*y;
            }
            if(sigma_x2 == 0) return 0;
            return sigma_xy / sigma_x2; 
        }

        bool precompress_block(const std::shared_ptr<Range> &range) noexcept { return true;}

        void precompress_block_commit() noexcept {}

        bool predecompress_block(const std::shared_ptr<Range> &) { return true; }

        void save(uchar *&c) const {}

        void load(const uchar *&c, size_t &remaining_length) {}

        void print() const {
            std::cout << "Pattern-based predictor\n";
        }

        inline T estimate_error(const iterator &iter) const noexcept {
            return fabs(*iter - predict(iter));
        }

        inline T predict(const iterator &iter) const noexcept {
            return 0;
        }

        void clear() {}

    private:
        // void compute_max_pattern_id(const T * pattern){
        //     T max_pattern_value = 0;
        //     for(int i=0; i<pattern_size; i++){
        //         if(max_pattern_value < pattern[i]){
        //             max_pattern_value = pattern[i];
        //             max_pattern_id = i;
        //         }
        //     }            
        // }
        // double compute_pattern_error(const T * data, const T * pattern){
        //     T scale = compute_pattern_scale(data, pattern);
        //     double err = 0;
        //     for(int j=0; j<pattern_size; j++){
        //         err += (data[j] - scale * pattern[j])*(data[j] - scale * pattern[j]);
        //     }
        //     return err;            
        // }
        // int max_pattern_id = 0;
        size_t pattern_repeated_times, pattern_size = 0;
    };
}
#endif
