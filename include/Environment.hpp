#ifndef ENVIRONMENT_HPP
#define ENVIRONMENT_HPP

#include <functional>
#include <random>

#include "Configuration.hpp"

struct Observables {
  double Pressure{0.0}, PressureError{0.0}, OldPressure{0.0};
};
class Environment
{
private:
  Observables Properties;
  std::mt19937_64 re;
  std::uniform_real_distribution<double> urng;
  double TargetTemperature, LCube;
  double Delta {0.02};
  Configuration Particles;

public:
  Environment(const size_t npart, const double Temperature, const double lcube);
  ~Environment();
  void InitAtoms();
  void InitiPositions();
  void InitVelocities();
  void BoostAndRescale(const double& TotalSquaredVelocity, const Vec3d& TotalVelocity);
  void Evolve(const size_t maxsteps);
  Configuration TrialStep();
  bool Accept(const Configuration& trialconfig);
  inline Observables GetObservables() const {return Properties;}
  inline double GetVolume() const {return Particles.Volume();}
  inline double GetVOverN() const {return Particles.Volume()/Particles.size();}
  inline size_t ActualParticleSize() {return Particles.size();}
};

#endif