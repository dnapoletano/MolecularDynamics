#ifndef ATOM_HPP
#define ATOM_HPP

#include "Vec3d.hpp"

class Atom
{
private:
  Vec3d Position, Velocity;

public:
  Atom();
  Atom(const Vec3d& pos, const Vec3d& vel);

  ~Atom();

  inline void SetPosition(const Vec3d& newpos, const double lcube) {
    Position = PeriodicBoundary(newpos,lcube);
  }
  inline void SetVelocity(const Vec3d& newvel) {Velocity = newvel;}
  inline void IncrementPosition(const Vec3d& increment, const double lcube) {
    Position = PeriodicBoundary(Position + increment, lcube);
  }
  inline void IncrementVelocity(const Vec3d& increment) {Velocity = Velocity + increment;}

  inline Vec3d GetPosition() const {return Position;}
  inline Vec3d GetVelocity() const {return Velocity;}

  inline static Vec3d PeriodicBoundary(const Vec3d &x, const double &lcube)
  {
    Vec3d res{x};
    if(res.x > 0.0) while (res.x >= lcube) res.x -= lcube;
    else while(res.x < 0.0) res.x += lcube;
    if(res.y > 0.0) while (res.y >= lcube) res.y -= lcube;
    else while(res.y < 0.0) res.y += lcube;
    if(res.z > 0.0) while (res.z >= lcube) res.z -= lcube;
    else while(res.z < 0.0) res.z += lcube;
    return res;
  }
  inline static double Distance(const Atom& a, const Atom& b, const double& lcube) {
    return Vec3d::Distance(a.GetPosition(),b.GetPosition(), lcube);
  }
  inline friend std::ostream &operator<<(std::ostream &os, const Atom a)
  {
    os << " Atom Position : " << a.GetPosition() << "\n";
    os << " Atom Velocity : " << a.GetVelocity() << "\n";
    return os;
  }
};


#endif