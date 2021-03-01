#ifndef VEC3D_HPP
#define VEC3D_HPP

#include <iostream>
#include <cmath>

#include "utils.hpp"

struct Vec3d
{
  double x, y, z;
  double Mod()  const { return sqrt(x * x + y * y + z * z); }
  double Mod2() const { return x * x + y * y + z * z; }
  static double Distance(const Vec3d &a, const Vec3d &b, const double& lcube);
  static Vec3d  VecDistance(const Vec3d &a, const Vec3d &b, const double& lcube);
  inline Vec3d& operator+=(const Vec3d &b)
  {
    x+=b.x;
    y+=b.y;
    z+=b.z;
    return *this;
  }
  inline Vec3d &operator*=(const double &b)
  {
    x *= b;
    y *= b;
    z *= b;
    return *this;
  }
  inline Vec3d &operator/=(const double &b)
  {
    x /= b;
    y /= b;
    z /= b;
    return *this;
  }
};

inline Vec3d operator+(const Vec3d &a, const Vec3d &b)
{
  return Vec3d{a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vec3d operator-(const Vec3d &a,const  Vec3d &b)
{
  return Vec3d{a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vec3d operator*(const Vec3d &a, const double &b)
{
  return Vec3d{a.x * b, a.y * b, a.z * b};
}

inline Vec3d operator*(const double &b, const Vec3d &a)
{
  return Vec3d{a.x * b, a.y * b, a.z * b};
}

inline Vec3d operator/(const Vec3d &a, const double &b)
{
  return Vec3d{a.x/b, a.y/b, a.z/b};
}

inline std::ostream &operator<<(std::ostream &out, const Vec3d &a)
{
  return out << "{ " << a.x << ", " << a.y << ", " << a.z << "}";
}

inline Vec3d Vec3d::VecDistance(const Vec3d &a, const Vec3d &b, const double& lcube)
{
  Vec3d diff{a-b};
  while (diff.x > lcube / 2.)
    diff.x -= lcube;
  while (diff.y > lcube / 2.)
    diff.y -= lcube;
  while (diff.z > lcube / 2.)
    diff.z -= lcube;
  while (diff.x < -lcube / 2.)
    diff.x += lcube;
  while (diff.y < -lcube / 2.)
    diff.y += lcube;
  while (diff.z < -lcube / 2.)
    diff.z += lcube;
  return diff;
}

inline double Vec3d::Distance(const Vec3d &a, const Vec3d &b, const double& lcube)
{
  return VecDistance(a,b,lcube).Mod();
}

#endif