//
// Created by Kai Zhao on 4/28/20.
//

#ifndef SZ_CONFIG_HPP
#define SZ_CONFIG_HPP

#include "def.hpp"

namespace SZ {
    template<class T, uint N>
    class Config {
    public:
        Config(T _eb, std::array<size_t, N> _dims, T l1=1, T l2=1) : eb(_eb), dims(_dims),
            pattern_eb(l1*_eb), scale_eb(l2*_eb) {
            switch (N) {
                case 1:
                    block_size = 128;
                    break;
                case 2:
                    block_size = 16;
                    break;
                default:
                    // >= 3D
                    block_size = 6;
                    break;
            }
            stride = block_size;
            num = 1;
            for (const auto &d:_dims) {
                num *= d;
            }
        }

        std::array<size_t, N> dims = {0};
        size_t num;
        bool enable_lorenzo = true;
        bool enable_2ndlorenzo = false;
        bool enable_regression = true;
        bool enable_2ndregression = false;
        bool enable_lossless = true;
        size_t quant_bin = 2;
        uint block_size, stride;
        T eb;
        // for pastri
        T pattern_eb, scale_eb = 0;
        size_t num_patterns, pattern_repeated_times, pattern_size = 0;
        // for preprocess
        uint pred_dim = 0;
        bool transpose = false;
    };
}

#endif //SZ_CONFIG_HPP
