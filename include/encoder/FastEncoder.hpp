#ifndef _SZ_FAST_ENCODER_HPP
#define _SZ_FAST_ENCODER_HPP

namespace SZ {

    template<class T>
    class FastEncoder : public concepts::EncoderInterface<T> {
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

            FastEncoder(){};
            void preprocess_encode(const std::vector<T> &bins, int stateNum){}

            size_t encode(const std::vector<T> &bins, uchar *&bytes){
                BitEncoder encoder(reinterpret_cast<uint64_t*>(bytes));
                for(int i=0; i<bins.size(); i++){
                    switch(bins[i]){
                        // fixed Huffman encoding:
                        // 0->'10'
                        // 1->'110'
                        // 2->'0'
                        // 3->'111'
                        case 0:{
                            encoder.encode(1);
                            encoder.encode(0);
                            break;                            
                        }
                        case 1:{
                            encoder.encode(1);
                            encoder.encode(1);
                            encoder.encode(0);
                            break;                            
                        }
                        case 2:{
                            encoder.encode(0);
                            break;                            
                        }
                        case 3:{
                            encoder.encode(1);
                            encoder.encode(1);
                            encoder.encode(1);
                            break;                            
                        }
                        default:{
                            printf("bins out of range\n");
                            exit(0);
                        }
                    }
                }
                encoder.flush();
                size_t size = encoder.size() * sizeof(uint64_t);
                bytes += size;
                return size;
            }

            void postprocess_encode(){}

            void preprocess_decode(){}

            std::vector<T> decode(const uchar *&bytes, size_t targetLength){
                std::vector<T> bins(targetLength);
                BitDecoder decorder(reinterpret_cast<const uint64_t*>(bytes));
                for(int i=0; i<targetLength; i++){
                    if(decorder.decode()){
                        if(decorder.decode()){
                            if(decorder.decode()){
                                // '111'
                                bins[i] = 3;
                            }
                            else{
                                // '110'
                                bins[i] = 1;
                            }
                        }
                        else{
                            // '10'
                            bins[i] = 0;
                        }
                    }
                    else{
                        bins[i] = 2;
                    }
                }
                return bins;
            }

            void postprocess_decode(){}

            uint save(uchar *&c){ return 0;}

            void load(const uchar *&c, size_t &remaining_length){}
    };
}
#endif
