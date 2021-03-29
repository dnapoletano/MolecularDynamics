#include "Configuration.hpp"

Configuration Configuration::operator=(const Configuration& config)
{
  TotMomentum          = config.GetTotalMomentum();
  KinEnergy            = config.GetKineticEnergy();
  PotEnergy            = config.GetPotentialEnergy();
  TotEnergy            = config.GetTotalEnergy();
  Temperature          = config.GetTemperature();
  Pressure             = config.GetPressure();
  TotalSquaredVelocity = config.GetTotalSquaredVelocity();
  LCube                = config.GetLCube();
  for(size_t i{0}; i < this->size(); ++i){
    (*this)[i].SetPosition(config[i].GetPosition(),LCube);
    (*this)[i].SetVelocity(config[i].GetVelocity());
  }
  return *this;
}

Configuration Configuration::operator=(Configuration&& config)
{
  TotMomentum          = config.GetTotalMomentum();
  KinEnergy            = config.GetKineticEnergy();
  PotEnergy            = config.GetPotentialEnergy();
  TotEnergy            = config.GetTotalEnergy();
  Temperature          = config.GetTemperature();
  Pressure             = config.GetPressure();
  TotalSquaredVelocity = config.GetTotalSquaredVelocity();
  LCube                = config.GetLCube();
  for(size_t i{0}; i < this->size(); ++i){
    (*this)[i].SetPosition(config[i].GetPosition(),LCube);
    (*this)[i].SetVelocity(config[i].GetVelocity());
  }
  return *this;
}

void Configuration::CalcEnergy()
{
  double virial{0.0};
  PotEnergy = 0.0;
  TotalSquaredVelocity = 0.0;
  TotMomentum = Vec3d{0.0, 0.0, 0.0};
  for(Atoms::iterator it1{this->begin()}; it1!=this->end(); ++it1){
    TotMomentum += Constants::mass * it1->GetVelocity();
    TotalSquaredVelocity += it1->GetVelocity().Mod2();
    for(Atoms::iterator it2{it1+1}; it2!=this->end();++it2){
      const Vec3d d {Vec3d::Distance(it1->GetPosition(), it2->GetPosition(), LCube)};
      const double r{d.Mod()};
      if(r < Constants::rCutoff and r != 0){
        PotEnergy += WanderVaals(r);
        virial    += d * WanderVaalsForce(d);
      }
    }
  }
  KinEnergy = Constants::mass * TotalSquaredVelocity / 2.;
  virial  += 2. * KinEnergy;
  TotEnergy = PotEnergy + KinEnergy;
  Temperature =  2. * KinEnergy / 3. / this->size();
  Pressure = virial / 3. / Volume();
}

double Configuration::PartitionFunction() const
{
  return exp(-PotEnergy/Temperature);
}

double Configuration::LogPartitionFunction() const
{
  return PotEnergy/Temperature;
}


double Configuration::WanderVaals(const double d){
  if(d==0.0) return 0.0;
  return 4. * Constants::epsilon *
           (pow(Constants::sigma / d, 12) - pow(Constants::sigma / d, 6) - Constants::eta);
}

Vec3d Configuration::WanderVaalsForce(const Vec3d& d){
  double r {d.Mod()};
  /// f = - U', such that F = - grad U = - rvec * U' / r;
  double f{  - 4. * Constants::epsilon *
           (- 12 * pow(Constants::sigma / r, 12) + 6 * pow(Constants::sigma / r, 6)) / r };
  return Vec3d{d.x * f / r, d.y * f /r, d.z * f /r};
}

