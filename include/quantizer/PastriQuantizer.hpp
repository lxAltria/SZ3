#ifndef _SZ_PASTRI_QUANTIZER_HPP
#define _SZ_PASTRI_QUANTIZER_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include "def.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include <bitset>

namespace SZ {

    template<class T>
    class PastriQuantizer : public LinearQuantizer<T> {
    public:
        class BitEncoder{
        public:
            BitEncoder(uint64_t * stream_begin_pos){
                stream_begin = stream_begin_pos;
                stream_pos = stream_begin;
                buffer = 0;
                position = 0;
            }
            void encode(uint64_t b){
                buffer += b << position;
                position ++;
                if(position == 64){
                    *(stream_pos ++) = buffer;
                    buffer = 0;
                    position = 0;
                }
            }
            void flush(){
                if(position){
                    *(stream_pos ++) = buffer;
                    buffer = 0;
                    position = 0;
                }
            }
            uint32_t size(){
                return (stream_pos - stream_begin);
            }
        private:
            uint64_t buffer = 0;
            uint8_t position = 0;
            uint64_t * stream_pos = NULL;
            uint64_t * stream_begin = NULL;
        };

        class BitDecoder{
        public:
            BitDecoder(uint64_t const * stream_begin_pos){
                stream_begin = stream_begin_pos;
                stream_pos = stream_begin;
                buffer = 0;
                position = 0;
            }
            uint64_t decode(){
                if(position == 0){
                    buffer = *(stream_pos ++);
                    position = 64;
                }
                uint64_t b = buffer & 1u;
                buffer >>= 1;
                position --;
                return b;
            }
            uint32_t size(){
                return (stream_pos - stream_begin);
            }
        private:
            uint64_t buffer = 0;
            uint8_t position = 0;
            uint64_t const * stream_pos = NULL;
            uint64_t const * stream_begin = NULL;
        };
        PastriQuantizer(T eb, int r = 32768) : LinearQuantizer<T>(eb, r) {
            eb_exp = 0;
            frexp(eb, &eb_exp);
            eb_exp -= 1;            
        }

        inline int quantize(T data, T pred) {
            T diff = data - pred;
            if(diff > 0){
                int quant_index = (int) (diff * this->error_bound_reciprocal) + 1;
                return (quant_index >> 1) + this->radius;
            }
            else{
                int quant_index = (int) (diff * this->error_bound_reciprocal) - 1;
                return (quant_index / 2) + this->radius;
            }
        }

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
                data = pred + quant_index * this->error_bound;
                return quant_index_shifted;
            } else {
                // unpred.push_back(data);
                this->unpred.push_back(diff);
                // recover decompressed
                int data_exp = 0;
                frexp(diff, &data_exp);
                int num_bitplanes = data_exp - eb_exp;
                // TODO: add 32 bit
                unsigned char shift_bits = 53 - num_bitplanes;
                uint64_t * tmp = reinterpret_cast<uint64_t*>(&diff);
                (*tmp) = ((*tmp) >> shift_bits) << shift_bits;
                data = pred + diff;
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
        // save unpredictable data in each block
        void save(unsigned char *&c) {
            *reinterpret_cast<size_t *>(c) = this->unpred.size();
            c += sizeof(size_t);
            embedded_encoding(c);
            this->unpred.clear();            
        }
        // load unpredictable data in each block
        void load(const unsigned char *&c, size_t &remaining_length) {
            std::vector<T>& unpred = this->unpred;
            unpred.clear();
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            unpred = std::vector<T>(unpred_size);
            if(unpred_size) embedded_decoding(c, unpred_size);
            this->index = 0;
        }
    private:
        void embedded_encoding(unsigned char* &c){
            auto tmp = c;
            const int block_size = 32;
            size_t n = this->unpred.size();
            if(n == 0) return;
            const T * data = this->unpred.data();
            uint64_t int_data_buffer[block_size];
            unsigned char * encoded_sign_pos = c;
            unsigned char * encoded_data_pos = c + (((n - 1) / 64 + 1) * 8);
            T max_val = 0;
            for(int i=0; i<n; i++){
                if(fabs(data[i]) > max_val){
                    max_val = fabs(data[i]);
                }
            }
            int data_exp = 0;
            frexp(max_val, &data_exp);
            const int num_bitplanes = data_exp - eb_exp;
            // std::cout << data_exp << " " << eb_exp << ", num_bitplanes = " << num_bitplanes << std::endl;
            // encode exponent
            *reinterpret_cast<int *>(encoded_data_pos) = num_bitplanes;
            encoded_data_pos += sizeof(int);
            const T * data_pos = data;
            BitEncoder sign_encoder(reinterpret_cast<uint64_t*>(encoded_sign_pos));
            BitEncoder data_encoder(reinterpret_cast<uint64_t*>(encoded_data_pos));
            int i = 0;
            // std::cout << n << ", sign bytes = " << (((n - 1) / 64 + 1) * 8) << "\n";
            if(n >= block_size){
                for(; i<n - block_size; i+=block_size){
                    for(int j=0; j<block_size; j++){
                        T cur_data = *(data_pos++);
                        T shifted_data = ldexp(cur_data, num_bitplanes - data_exp);
                        int64_t fix_point = (int64_t) shifted_data;
                        bool sign = cur_data < 0;
                        int_data_buffer[j] = sign ? -fix_point : +fix_point;
                        sign_encoder.encode(sign);
                    }
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        for (int j=0; j<block_size; j++){
                            bool bit = (int_data_buffer[j] >> k) & 1u;
                            data_encoder.encode(bit);
                        }
                    }
                }
            }
            // rest
            if(i != n){
                int rest_size = n - i;
                for(int j=0; j<rest_size; j++){
                    T cur_data = *(data_pos++);
                    T shifted_data = ldexp(cur_data, num_bitplanes - data_exp);
                    int64_t fix_point = (int64_t) shifted_data;
                    uint32_t sign = cur_data < 0;
                    int_data_buffer[j] = sign ? -fix_point : +fix_point;
                    sign_encoder.encode(sign);
                }
                for(int k=num_bitplanes - 1; k>=0; k--){
                    for (int j=0; j<rest_size; j++){
                        bool bit = (int_data_buffer[j] >> k) & 1u;
                        data_encoder.encode(bit);
                    }
                }
            }
            sign_encoder.flush();
            data_encoder.flush();
            c = reinterpret_cast<unsigned char*>(encoded_data_pos + data_encoder.size() * 8);
        }
        void embedded_decoding(const unsigned char * &c, size_t unpred_size){
            auto tmp = c;
            const int block_size = 32;
            T * data = this->unpred.data();
            size_t n = unpred_size;
            uint64_t int_data_buffer[block_size];
            const unsigned char * encoded_sign_pos = c;
            const unsigned char * encoded_data_pos = c + (((n - 1) / 64 + 1) * 8);
            int num_bitplanes = *reinterpret_cast<const int*>(encoded_data_pos);
            encoded_data_pos += sizeof(int);
            BitDecoder sign_decoder(reinterpret_cast<const uint64_t*>(encoded_sign_pos));
            BitDecoder data_decoder(reinterpret_cast<const uint64_t*>(encoded_data_pos));
            T * data_pos = data;
            int i = 0;
            if(n >= block_size){
                for(; i<n - block_size; i+=block_size){
                    memset(int_data_buffer, 0, block_size * sizeof(uint64_t));
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        for (int j=0; j<block_size; j++){
                            int_data_buffer[j] += data_decoder.decode() << k;
                        }
                    }
                    for(int j=0; j<block_size; j++){
                        bool sign = sign_decoder.decode();
                        T cur_data = ldexp((T)int_data_buffer[j], eb_exp);
                        *(data_pos++) = sign ? -cur_data : cur_data;
                    }
                }
            }
            // rest
            if(i != n){
                int rest_size = n - i;
                memset(int_data_buffer, 0, rest_size * sizeof(uint64_t));
                for(int k=num_bitplanes - 1; k>=0; k--){
                    for (int j=0; j<rest_size; j++){
                        int_data_buffer[j] += data_decoder.decode() << k;
                    }
                }
                for(int j=0; j<rest_size; j++){
                    bool sign = sign_decoder.decode();
                    T cur_data = ldexp((T)int_data_buffer[j], eb_exp);
                    *(data_pos++) = sign ? -cur_data : cur_data;
                }
            }
            c = encoded_data_pos + data_decoder.size() * 8;
        }

        // exponent for error bound
        int eb_exp;
    };

}
#endif