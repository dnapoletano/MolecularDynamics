#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <Atom.hpp>

typedef std::vector<class Atom> Atoms;

class Configuration : private Atoms
{
private:
  Vec3d TotMomentum;
  double KinEnergy, PotEnergy, TotEnergy, Temperature, Pressure;
  double TotalSquaredVelocity;
  double LCube;

public:
  Configuration()
  : Atoms{}, TotMomentum{Vec3d{0.0,0.0,0.0}},
    KinEnergy{0.0},
    PotEnergy{0.0},
    TotEnergy{0.0},
    Temperature{0.0},
    Pressure{0.0},
    TotalSquaredVelocity{0.0},
    LCube{30.0}
  {}

  /// copy and move constructors
  Configuration(const Atoms& atoms)
  : Atoms{}, TotMomentum{Vec3d{0.0,0.0,0.0}},
    KinEnergy{0.0},
    PotEnergy{0.0},
    TotEnergy{0.0},
    Temperature{0.0},
    Pressure{0.0},
    TotalSquaredVelocity{0.0},
    LCube{30.0}
  {}
  Configuration(Atoms&& atoms)
  : Atoms{}, TotMomentum{Vec3d{0.0,0.0,0.0}},
    KinEnergy{0.0},
    PotEnergy{0.0},
    TotEnergy{0.0},
    Temperature{0.0},
    Pressure{0.0},
    TotalSquaredVelocity{0.0},
    LCube{30.0}
  {}

  Configuration(const Configuration& config)
  : Atoms{config}, TotMomentum{config.GetTotalMomentum()},
    KinEnergy{config.GetKineticEnergy()},
    PotEnergy{config.GetPotentialEnergy()},
    TotEnergy{config.GetTotalEnergy()},
    Temperature{config.GetTemperature()},
    Pressure{config.GetPressure()},
    TotalSquaredVelocity{config.GetTotalSquaredVelocity()},
    LCube{config.GetLCube()}
  {}

  Configuration(Configuration&& config)
  : Atoms{config}, TotMomentum{config.GetTotalMomentum()},
    KinEnergy{config.GetKineticEnergy()},
    PotEnergy{config.GetPotentialEnergy()},
    TotEnergy{config.GetTotalEnergy()},
    Temperature{config.GetTemperature()},
    Pressure{config.GetPressure()},
    TotalSquaredVelocity{config.GetTotalSquaredVelocity()},
    LCube{config.GetLCube()}
  {}

  Configuration operator=(const Configuration& config);
  Configuration operator=(Configuration&& config);

  using Atoms::at;
  using Atoms::begin;
  using Atoms::clear;
  using Atoms::end;
  using Atoms::size;
  using Atoms::push_back;
  using Atoms::resize;
  using Atoms::const_iterator;
  using Atoms::iterator;
  using Atoms::operator=;
  using Atoms::operator[];

  void CalcEnergy();
  double PartitionFunction() const;
  double LogPartitionFunction() const;

  /// utilities
  inline const Vec3d& GetTotalMomentum()   const {return TotMomentum;}
  inline double GetKineticEnergy()         const {return KinEnergy;}
  inline double GetPotentialEnergy()       const {return PotEnergy;}
  inline double GetTotalEnergy()           const {return TotEnergy;}
  inline double GetTemperature()           const {return Temperature;}
  inline double GetPressure()              const {return Pressure;}
  inline double GetTotalSquaredVelocity()  const {return TotalSquaredVelocity;}
  inline double GetLCube()                 const {return LCube;}
  inline double Volume()                   const {return LCube * LCube * LCube;}

  inline void SetLCube(const double lcube)       { LCube = lcube;}

  inline friend std::ostream &operator<<(std::ostream &os, const Configuration as)
  {
    os<< "Configuration {\n";
    os<< std::setw(20) << "KinEnergy          :" << std::setw(20) << as.GetKineticEnergy() << "\n"
      << std::setw(20) << "PotEnergy          :" << std::setw(20) << as.GetPotentialEnergy() << "\n"
      << std::setw(20) << "TotalEnergy        :" << std::setw(20) << as.GetTotalEnergy() << "\n"
      << std::setw(20) << "Temperature        :" << std::setw(20) << as.GetTemperature() << "\n"
      << std::setw(20) << "Pressure           :" << std::setw(20) << as.GetPressure() << "\n"
      << std::setw(20) << "TotalMomentum      :" << std::setw(20) << as.GetTotalMomentum() << "\n";
    for(const auto& a: as){
      os << std::setw(20) << a << "\n";
    }
    return os;
  }

};
#endif