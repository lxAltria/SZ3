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

    template<class T>
    class PastriQuantizer : public LinearQuantizer<T> {
    public:
        PastriQuantizer(T eb, int r = 32768) : LinearQuantizer<T>(eb, r) {}

        inline int quantize_and_overwrite(T &data, T pred) {

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
                    this->unpred.push_back(diff);
                    return 0;
                } else {
                    data = decompressed_data;
                    return quant_index_shifted;
                }
            } else {
                // std::cout << data << " " << diff << std::endl;
                // unpred.push_back(data);
                this->unpred.push_back(diff);
                return 0;
            }
        }
        T recover(T pred, int quant_index) {
            if (quant_index) {
                return pred + 2 * (quant_index - this->radius) * this->error_bound;
            } else {
                return pred + this->unpred[this->index++];
            }
        }

        void save(unsigned char *&c) const {
            uchar* tmp = c;
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000010;
            c += sizeof(uchar);
            *reinterpret_cast<T *>(c) = this->error_bound;
            c += sizeof(T);
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            const std::vector<T>& unpred = this->unpred;
            *reinterpret_cast<size_t *>(c) = unpred.size();
            c += sizeof(size_t);
            if(unpred.size()){
                std::cout << "unpredictable size = " << this->unpred.size() << std::endl;
                std::vector<uchar> sign(unpred.size());
                T max_val = 0;
                for(int i=0; i<unpred.size(); i++){
                    if(max_val < fabs(unpred[i])){
                        max_val = unpred[i];
                    }
                    sign[i] = (unpred[i] < 0);
                }
                uchar buffer = 0;
                int bit_position = 7;
                // record sign
                std::cout << "Recorded size before sign = " << c - tmp << std::endl;
                for(int i=0; i<unpred.size(); i++){
                    buffer += sign[i] << bit_position;
                    bit_position --;
                    if(bit_position < 0){
                        bit_position = 7;
                        *(c++) = buffer;
                        buffer = 0;
                    }
                }
                // flush
                if(bit_position != 7) *(c++) = buffer;
                // record mantissa
                std::cout << "Recorded size after sign = " << c - tmp << std::endl;
                int aligned_exp = 0;
                frexp(max_val, &aligned_exp);
                int eb_exp = 0;
                frexp(this->error_bound, &eb_exp);
                // number of bits recorded per data
                int num_bits = aligned_exp - eb_exp + 1;
                *reinterpret_cast<int*>(c) = aligned_exp;
                c += sizeof(int);
                std::cout << "aligned_exp = " << aligned_exp << std::endl;
                std::cout << "Recorded size after exp = " << c - tmp << std::endl;
                std::cout << "num_bits = " << num_bits << ", num_data = " << unpred.size() << std::endl;
                uint64_t mantissa_buffer = 0;
                std::vector<uint64_t> fix_points(unpred.size());
                for(int i=0; i<unpred.size(); i++){
                    T cur_data = sign[i] ? -unpred[i] : unpred[i];
                    T shifted_data = ldexp(cur_data, num_bits - 1);
                    fix_points[i] = (uint64_t) shifted_data;
                }
                std::cout << "unpred[0] = " << unpred[0] << std::endl;
                size_t compact_size = compact(fix_points, num_bits, reinterpret_cast<uint64_t*>(c));
                c += compact_size;
                std::cout << "Recorded size = " << c - tmp << std::endl;
            }
        }
        void load(const unsigned char *&c, size_t &remaining_length) {
            const uchar* tmp = c;
            c += sizeof(uchar);
            this->error_bound = *reinterpret_cast<const T *>(c);
            std::cout << "eb = " << this->error_bound << std::endl;
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            c += sizeof(T);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            std::vector<T>& unpred = this->unpred;
            unpred = std::vector<T>(unpred_size, 0);
            if(unpred_size){
                // recover sign
                std::cout << "Recorded size before sign = " << c - tmp << std::endl;
                std::vector<bool> sign(unpred_size);
                uchar buffer = *(c++);
                int bit_position = 7;
                for(int i=0; i<unpred_size; i++){
                    sign[i] = (buffer >> bit_position) & 1u;
                    bit_position --;
                    if((bit_position < 0) && (i+1 != unpred_size)){
                        buffer = *(c++);
                        bit_position = 7;
                    }
                }
                std::cout << "Recorded size after sign = " << c - tmp << std::endl;
                // recover data            
                int aligned_exp = *reinterpret_cast<const int*>(c);
                c += sizeof(int);
                std::cout << "aligned_exp = " << aligned_exp << std::endl;
                std::cout << "Recorded size after exp = " << c - tmp << std::endl;
                int eb_exp = 0;
                frexp(this->error_bound, &eb_exp);
                // number of bits recorded per data
                int num_bits = aligned_exp - eb_exp + 1;
                std::vector<uint64_t> fix_points = decompact<uint64_t>(reinterpret_cast<const uint64_t*>(c), unpred_size, num_bits);
                for(int i=0; i<unpred_size; i++){
                    T cur_data = ldexp((T)fix_points[i], -num_bits + 1);
                    unpred[i] = sign[i] ? -cur_data : cur_data;
                }
                std::cout << "unpred[0] = " << unpred[0] << std::endl;
                c += compute_compact_size(unpred_size, num_bits);
                std::cout << "Recorded size = " << c - tmp << std::endl;
            }
        }
    private:
        template<class T1>
        std::vector<T1> get_mask_array() const{
            std::vector<T1> mask = std::vector<T1>(sizeof(T1) * UINT64_BITS + 1, 0);
            T1 m = 1;
            for(int i=1; i<mask.size(); i++){
                mask[i] = m;
                m = (m << 1) + 1;
            }
            return mask;
        }

        inline size_t compute_compact_size(size_t size, size_t index_size) const{
            return ((size * index_size - 1) / UINT64_BITS + 1) * (UINT64_BITS / UINT8_BITS);
        }
        /* 
        @params encoded: addresses of encoded bitplanes
        @params offset: position of encoded bitplanes
        return compact array
        */
        template<class T1>
        size_t compact(const std::vector<T1>& data, uint64_t index_size, uint64_t * compact_data) const{
            static_assert(std::is_unsigned<T1>::value, "codec_utils compact: input array must be unsigned integers.");
            static_assert(std::is_integral<T1>::value, "codec_utils compact: input array must be unsigned integers.");
            auto mask = get_mask_array<T1>();
            uint64_t * compact_data_pos = compact_data;
            uint64_t buffer = 0;
            uint8_t rest_bits = UINT64_BITS;
            for(int i=0; i<data.size(); i++){
                T1 cur_data = data[i];
                if(index_size <= rest_bits){
                    rest_bits -= index_size;
                    buffer <<= index_size;
                    buffer += cur_data;
                }
                else{
                    buffer <<= rest_bits;
                    buffer += cur_data >> (index_size - rest_bits);
                    *(compact_data_pos ++) = buffer;
                    buffer = cur_data & mask[index_size - rest_bits];
                    rest_bits = UINT64_BITS + rest_bits - index_size;
                }
            }
            // flush buffer
            if(rest_bits != UINT64_BITS){
                *(compact_data_pos ++) = buffer << rest_bits;
            }
            return compute_compact_size(data.size(), index_size);
        }

        template<class T1>
        std::vector<T1> decompact(const uint64_t * compact_data, uint32_t size, uint64_t index_size) const{
            static_assert(std::is_unsigned<T1>::value, "codec_utils decompact: input array must be unsigned integers.");
            static_assert(std::is_integral<T1>::value, "codec_utils decompact: input array must be unsigned integers.");
            auto mask = get_mask_array<uint64_t>();
            std::vector<T1> data(size, 0);
            const uint64_t * compact_data_pos = compact_data;
            uint64_t buffer = 0;
            uint8_t rest_bits = 0;
            for(int i=0; i<size; i++){
                T1 cur_data = 0;
                if(index_size <= rest_bits){
                    rest_bits -= index_size;
                    cur_data = buffer >> rest_bits;
                    buffer = buffer & mask[rest_bits];
                }
                else{
                    cur_data = buffer << (index_size - rest_bits);
                    buffer = *(compact_data_pos ++);
                    rest_bits = UINT64_BITS + rest_bits - index_size;
                    cur_data += buffer >> rest_bits;
                    buffer = buffer & mask[rest_bits];
                }
                data[i] = cur_data;
            }
            return data;
        }

    };
}
#endif
