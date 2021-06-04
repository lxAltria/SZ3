#ifndef _SZ_EST_LORENZO_PREDICTOR_HPP
#define _SZ_EST_LORENZO_PREDICTOR_HPP

#include "def.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "utils/Iterator.hpp"
#include <cassert>

namespace SZ {

    template<class T>
    class EstLorenzoPredictor3D : public LorenzoPredictor<T, 3, 1> {
    public:
        using iterator = typename multi_dimensional_range<T, 3>::iterator;
    	EstLorenzoPredictor3D(T eb){
    		noise_1d = eb * 0.5;
    		noise_2d = eb * 0.81;
    		noise_3d = eb * 1.22;
    	}
    	inline T estimate_error_1d(const iterator &iter) const noexcept {
            return fabs(*iter - predict_1d(iter)) + noise_1d;
        }
    	inline T estimate_error_2d(const iterator &iter) const noexcept {
            return fabs(*iter - predict_2d(iter)) + noise_2d;
        }
    	inline T estimate_error_3d(const iterator &iter) const noexcept {
            return fabs(*iter - predict_3d(iter)) + noise_3d;
        }

    	inline T predict_1d(const iterator &iter) const noexcept {
            return iter.prev(0, 0, 1);
        }
    	inline T predict_2d(const iterator &iter) const noexcept {
            return iter.prev(0, 0, 1) + iter.prev(0, 1, 0) - iter.prev(0, 1, 1);
        }
    	inline T predict_3d(const iterator &iter) const noexcept {
            return iter.prev(0, 0, 1) + iter.prev(0, 1, 0) + iter.prev(1, 0, 0)
                   - iter.prev(0, 1, 1) - iter.prev(1, 0, 1) - iter.prev(1, 1, 0)
                   + iter.prev(1, 1, 1);
        }
    private:
    	T noise_1d = 0;
    	T noise_2d = 0;
    	T noise_3d = 0;    	
    };

}
#endif