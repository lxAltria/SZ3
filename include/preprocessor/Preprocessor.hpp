#ifndef _SZ_PREPROCESSOR_HPP
#define _SZ_PREPROCESSOR_HPP

#include "utils/Concepts.hpp"
#include "utils/Config.hpp"
#include "def.hpp"

namespace SZ {
    namespace concepts {

        template<class T, uint N>
        class PreprocessorInterface {
        public:

            virtual ~PreprocessorInterface() = default;

            virtual void process(T * data, SZ::Config<T, N>& conf) const = 0;

        };

    }
}

#endif
