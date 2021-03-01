#ifndef _SZ_INTEGER_QUANTIZER_HPP
#define _SZ_INTEGER_QUANTIZER_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include "def.hpp"
#include "quantizer/Quantizer.hpp"
#include <bitset>

namespace SZ {

    #define UINT8_BITS 8
    #define UINT64_BITS 64
// data with T type
// return int
    template<class T>
    class PredictionBasedQuantizer : public concepts::QuantizerInterface<T> {
    protected:
        T error_bound;
        T error_bound_reciprocal;
        int radius; // quantization interval radius
    public:
        ~PredictionBasedQuantizer() = default;

        PredictionBasedQuantizer() = default;

        PredictionBasedQuantizer(PredictionBasedQuantizer const &) = default;

        PredictionBasedQuantizer(PredictionBasedQuantizer &&) = default;

        PredictionBasedQuantizer &operator=(PredictionBasedQuantizer const &) = default;

        PredictionBasedQuantizer &operator=(PredictionBasedQuantizer &&) = default;

        void precompress_data() const {}

        void postcompress_data() const {}

        void predecompress_data() const {}

        void postdecompress_data() const {}

//        void predecompress_block() {}
//
//        void precompress_block() {}


        PredictionBasedQuantizer(T eb, int r) : error_bound(eb),
                                                error_bound_reciprocal(1.0 / eb),
                                                radius(r) {}

        int get_radius() const { return radius; }

        T get_eb() const { return error_bound; }

    };

    template<class T>
    class LinearQuantizer : public PredictionBasedQuantizer<T> {
    public:
        LinearQuantizer(T eb, int r = 32768) : PredictionBasedQuantizer<T>(eb, r) {
        }

        using value_type = T;
        using reference = T &;

        // quantize the data with a prediction value, and returns the quantization index
        int quantize(T data, T pred);

        // quantize the data with a prediction value, and returns the quantization index and the decompressed data
        // int quantize(T data, T pred, T& dec_data);
        int quantize_and_overwrite(T &data, T pred);

        // recover the data using the quantization index
        T recover(T pred, int quant_index);

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000010;
            c += 1;
            // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            *reinterpret_cast<T *>(c) = this->error_bound;
            c += sizeof(T);
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            *reinterpret_cast<size_t *>(c) = unpred.size();
            c += sizeof(size_t);
            memcpy(c, unpred.data(), unpred.size() * sizeof(T));
            c += unpred.size() * sizeof(T);
        }

        void load(const unsigned char *&c, size_t &remaining_length) {
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            this->error_bound = *reinterpret_cast<const T *>(c);
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            c += sizeof(T);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
            c += unpred_size * sizeof(T);
            // std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            // reset index
            index = 0;
        }

        void clear() {
            unpred.clear();
        }

    protected:
        std::vector<T> unpred;
        size_t index = 0; // used in decompression only
    };

    template<class T>
    int LinearQuantizer<T>::quantize(T data, T pred) {
        int radius = this->radius;
        // compute quantization index
        int quant_index = (int) ((data - pred) * this->error_bound_reciprocal);
        quant_index = (quant_index > 0) ? (quant_index + 1) / 2 : (quant_index - 1) / 2;
        // shift quantization index, set overbound to 0
        return (quant_index > 0) ? (quant_index < radius ? quant_index + radius : 0) : (quant_index > -radius ? quant_index +
                                                                                                                radius : 0);
    }

//    template<class T>
//    inline int LinearQuantizer<T>::quantize_and_overwrite(T &data, T pred) {
//
//        int quant_index = floor((data - pred) * this->error_bound_reciprocal_divided_by_2 + 0.5);
//        int quant = 0;
//        if (abs(quant_index) < this->radius) {
//            T decom = pred + quant_index * this->error_bound_times_2;
//            if (fabs(decom - data) < this->error_bound) {
//                data = decom;
//                quant = quant_index + this->radius;
//            }
//        }
//        if (quant == 0) {
//            unpred.push_back(data);
//        }
//        return quant;
//    }

    template<class T>
    inline int LinearQuantizer<T>::quantize_and_overwrite(T &data, T pred) {

        T diff = data - pred;
        int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
        if ((quant_index >= 0) && (quant_index < this->radius * 2)) {
            quant_index >>= 1;
            int half_index = quant_index;
            quant_index <<= 1;
            int quant_index_shifted;
            if (diff < 0) {
                quant_index = -quant_index;
                quant_index_shifted = this->radius - half_index;
            } else {
                quant_index_shifted = this->radius + half_index;
            }
            T decompressed_data = pred + quant_index * this->error_bound;
            if (fabs(decompressed_data - data) > this->error_bound) {
                unpred.push_back(data);
                return 0;
            } else {
                data = decompressed_data;
                return quant_index_shifted;
            }
        } else {
            // std::cout << data << " " << diff << std::endl;
            // unpred.push_back(data);
            unpred.push_back(data);
            return 0;
        }
    }

    template<class T>
    T LinearQuantizer<T>::recover(T pred, int quant_index) {
        if (quant_index) {
            return pred + 2 * (quant_index - this->radius) * this->error_bound;
        } else {
            return unpred[index++];
        }
    }
}
#endif
