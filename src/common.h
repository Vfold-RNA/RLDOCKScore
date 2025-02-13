#pragma once

#include <limits>
#include <cmath>
#include <vector>

namespace rlscore {

using Float = double;
using Floats = std::vector<Float>;
using Size_Type = std::size_t;
const Size_Type k_max_size_type = (std::numeric_limits<Size_Type>::max)();
// using Size_Types = std::vector<std::size_t>;
// using Float_Pair = std::pair<Float,Float>;
// using Float_Pairs = std::vector<Float_Pair>;
// ############constant related############
const Float k_pi = Float(3.1415926535897931);
const Float k_max_float = (std::numeric_limits<Float>::max)();
const Float k_min_float = (std::numeric_limits<Float>::min)();
const unsigned k_max_unsigned = (std::numeric_limits<unsigned>::max)();
const int k_max_int = (std::numeric_limits<int>::max)();
const Float k_epsilon = std::numeric_limits<Float>::epsilon();
const Float k_tolerance = Float(0.001);
const Float k_not_a_num = std::sqrt(Float(-1)); // FIXME? check


// template<unsigned n>
// inline Float int_pow(Float x) {
// 	return int_pow<n-1>(x)*x; // tests seem to suggest that even for largish n's this optimizes well
// }
// template<>
// inline Float int_pow<0>(Float x) {
// 	return 1;
// }


// inline const Float pK_to_energy(Float pK) { return k_pK_to_energy_factor * pK; } //return the energy from give pK
// Size_Type float_to_size_type(Float x, Size_Type max_size_t) { // return a value in [0, max_sz]
//     if(x <= 0) return 0;
//     if(x >= max_size_t) return max_size_t;
// 	Size_Type tmp = static_cast<Size_Type>(x);
// 	return (std::min)(tmp, max_size_t); // Size_Type -> Float cast loses precision. 'min' is to guard against returning values out of range
// }
// inline const bool not_max(Float x) {
// 	return (x < 0.1 * k_max_float);
// }
inline const bool eq(Float a, Float b) {
	return std::abs(a - b) < k_tolerance;
}
inline const Float square(const Float a) {
	return a*a;
}

}


//define macros
// #define TMD_FOR_IN(i, v) for(Size_Type VINA_MACROS_TMP = (v).size(), (i) = 0; (i) < VINA_MACROS_TMP; ++(i))
// #define TMD_FOR(i, n)    for(Size_Type VINA_MACROS_TMP = (n),        (i) = 0; (i) < VINA_MACROS_TMP; ++(i))
// #define TMD_U_FOR(i, n)  for(unsigned    VINA_MACROS_TMP = (n),        (i) = 0; (i) < VINA_MACROS_TMP; ++(i))

// #define TMD_RANGE(i, a, b)   for(Size_Type VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))
// #define TMD_U_RANGE(i, a, b) for(unsigned    VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))
// #define TMD_I_RANGE(i, a, b) for(int         VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))

// #define TMD_LOOP_CONST(t, i, v) for(t::const_iterator (i) = (v).begin(); (i) != (v).end(); ++(i))
// #define TMD_LOOP(t, i, v) for(t::iterator       (i) = (v).begin(); (i) != (v).end(); ++(i))

// #define TMD_SHOW(x)       do { std::cout << #x << " = " << (x) << std::endl; } while(false)
// #define TMD_SHOW_FAST(x)  do { std::cout << #x << " = " << (x) <<      '\n'; } while(false)
// #define TMD_ESHOW(x)      do { std::cerr << #x << " = " << (x) <<      '\n'; } while(false)