#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace Constants
{
  constexpr double NA{6.02214076e+23};
  constexpr double rCutoff{2.}, epsilon{1.0}, sigma{1.0}, mass{1.0};
  constexpr double Delta{1.e-3};
  const double eta{pow(Constants::sigma / Constants::rCutoff, 12) -
                   pow(Constants::sigma / Constants::rCutoff, 6)};
}
namespace utils
{
 enum om {
      none    ,
      reset   ,
      bold    ,
      blink   ,
      underln ,
      blackbg ,
      red     ,
      green   ,
      blue    ,
      brown   ,
      violet  ,
      lblue   ,
      grey    ,
      redbg   ,
      greenbg ,
      bluebg  ,
      brownbg ,
      violetbg,
      lbluebg ,
      greybg  ,
      yellow,
      pm,
  };

inline std::ostream& operator<<(std::ostream &str,const om modifier)
{
  switch (modifier) {
  case om::reset:
    return str << "\033[0m";
  case om::bold:
    return str << "\033[1m";
  case om::underln:
    return str << "\033[4m";
  case om::blink:
    return str << "\033[5m";
  case om::blackbg:
    return str << "\033[7m";
  case om::red:
    return str << "\033[31m";
  case om::green:
    return str << "\033[32m";
  case om::brown:
    return str << "\033[33m";
  case om::blue:
    return str << "\033[34m";
  case om::violet:
    return str << "\033[35m";
  case om::lblue:
    return str << "\033[36m";
  case om::grey:
    return str << "\033[37m";
  case om::redbg:
    return str << "\033[41m";
  case om::greenbg:
    return str << "\033[42m";
  case om::brownbg:
    return str << "\033[43m";
  case om::bluebg:
    return str << "\033[44m";
  case om::violetbg:
    return str << "\033[45m";
  case om::lbluebg:
    return str << "\033[46m";
  case om::greybg:
    return str << "\033[47m";
  case om::yellow:
    return str << "\e[0;33m";
  case om::pm:
    return str << " \u00B1 ";
  case om::none:
    return str;
  }
  return str;
}


  template <class T>
  inline std::ostream &operator<<(std::ostream &os, const std::vector<T> vec)
  {
    os << "{";
    for (auto v{vec.begin()}; v != vec.end(); ++v)
    {
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
    for (size_t i{0}; i < n; ++i)
    {
      if (div > 0)
      {
        const T step{delta / static_cast<T>(div)};
        rtn[i] = (xmin + i * step);
      }
      else
      {
        rtn[i] = (xmin + delta);
      }
    }
    return rtn;
  }

  template <class T>
  inline std::vector<T> logspace(T xmin, T xmax, size_t n, T base,
                                 bool include_end = true)
  {
    std::vector<T> rtn{linspace<T>(xmin, xmax, n, include_end)};
    for (auto &r : rtn)
    {
      r = static_cast<T>(pow(base, r));
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

  /// Format c++ string using printf like capabilities
  template <typename... Args>
  std::string string_format(const std::string &format, Args... args)
  {
    int size = snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
    if (size <= 0)
    {
      throw std::runtime_error("Error during formatting.");
    }
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
  }
}
#endif