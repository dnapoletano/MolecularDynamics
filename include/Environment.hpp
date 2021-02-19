#ifndef ENVIRONMENT_HPP
#define ENVIRONMENT_HPP

#include <random>

#include "Atom.hpp"

class Environment
{
private:
  std::mt19937_64 re;
  std::uniform_real_distribution<double> urng;
  /// we could use a normal generator for the initial velocities
  /// std::normal_distribution<double> nrng;
  double TargetTemperature, LCube;
  Atoms Particles;
public:
  Environment(const size_t npart, const double Temperature, const double lcube);
  ~Environment();
  void InitAtoms();
  void BoostAndRescale(const double& TotalSquaredVelocity, const Vec3d& TotalVelocity);
  double KineticEnergy() const;
  double PotentialEnergy() const;
  double TotalEnergy() const;
  inline double Temperature() const {return 2. * KineticEnergy() / 3. / Particles.size();}
  double TotalMomentum() const;
};

#endif