#ifndef EXTRAMATH_H
#define EXTRAMATH_H

// extra misc. math functions used throughout the project

template<typename T>
inline T pow2(T x) { return x*x; }

template<typename T>
inline T pow3(T x) { return x*x*x; }

template<typename T>
inline T pow4(T x) { return pow2( pow2(x) ); }

template<typename T>
inline T pow5(T x) { return pow2(x) * pow3(x); }

template<typename T>
inline T pow6(T x) { return pow3( pow2(x) ); }

template<typename T>
inline T pow7(T x) { return x * pow3( pow2(x) ); }


#endif // EXTRAMATH_H
