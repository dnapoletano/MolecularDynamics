#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace Constants {
  constexpr double NA {6.02214076e+23};
  constexpr double rCutoff {2.}, epsilon {1.0}, sigma {1.0}, mass {1.0};
  constexpr double Delta{1.e-3};
  const double eta {pow(Constants::sigma / Constants::rCutoff, 12) -
    pow(Constants::sigma / Constants::rCutoff, 6)};
}

template <class T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> vec)
{
  os << "{";
  for (auto v{vec.begin()}; v != vec.end(); ++v) {
    if (v != vec.end() - 1)
      os << (*v) << ", ";
    else
      os << (*v);
  }
  os << "}\n";
  return os;
}

template <class T>
inline std::vector<T> linspace(T xmin, T xmax, size_t n,
                               bool include_end = true)
{
  size_t div = include_end ? (n - 1) : n;
  T delta{xmax - xmin};
  std::vector<T> rtn(n);
  for (size_t i{0}; i < n; ++i) {
    if (div > 0) {
      const T step{delta / static_cast<T>(div)};
      rtn[i] = (xmin + i * step);
    } else {
      rtn[i] = (xmin + delta);
    }
  }
  return rtn;
}

template <class T>
inline std::vector<T> logspace(T xmin, T xmax, size_t n, T base,
                               bool include_end = true)
{
  std::vector<T> rtn{linspace<T>(xmin,xmax,n,include_end)};
  for (auto& r: rtn){
    r = static_cast<T>(pow(base,r));
  }

  return rtn;
}

template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, double eps)
{
  return std::fabs(x - y) <=
             eps * std::fabs(x + y)
         // unless the result is subnormal
         || std::fabs(x - y) < eps;
}

inline double WanderVaals(const double d){
  if(d==0.0) return 0.0;
  return 4. * Constants::epsilon *
           (pow(Constants::sigma / d, 12) - pow(Constants::sigma / d, 6) - Constants::eta);
}

inline double WanderVaalsForce(const double d){
  return 4. * Constants::epsilon *
           (12 * pow(Constants::sigma / d, 12) - 6 * pow(Constants::sigma / d, 6)) / d ;
}

#endif