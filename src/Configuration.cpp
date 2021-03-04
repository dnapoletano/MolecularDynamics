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
      const double d{Atom::Distance((*it1), (*it2), LCube)};
      if(d < Constants::rCutoff and d != 0){
        PotEnergy += WanderVaals(d);
        virial  += d * abs(WanderVaalsForce(d));
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

