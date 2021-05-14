#pragma once
#if __cplusplus < 201703L

namespace std {
template<typename T1, typename T2>
using is_convertible_v = typename is_convertible<T1,T2>::value;

}
#endif
