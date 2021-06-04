#ifndef _SZ_EST_QUANTIZER_HPP
#define _SZ_EST_QUANTIZER_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include "def.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include <bitset>

namespace SZ {

    template<class T>
    class EstQuantizer : public LinearQuantizer<T> {
    public:
        EstQuantizer(T eb, int r = 32768) : LinearQuantizer<T>(eb, r) {}

        // return unshifted quantization index
        inline int quantize(T data, T pred) {
            T diff = data - pred;
            T quant_index_f = fabs(diff) * this->error_bound_reciprocal + 1;
            if (quant_index_f < this->radius * 2) {
                int quant_index = (int) quant_index_f;
                quant_index >>= 1;
                return quant_index;
            } else {
                return this->radius;
            }
        }
    };

}
#endif