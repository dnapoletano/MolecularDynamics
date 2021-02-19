#include <functional>
#include "Environment.hpp"

Environment::Environment(const size_t npart, const double Temperature, const double lcube)
: re{0}, urng{0.,1.}, TargetTemperature{Temperature}, LCube{lcube}
{
  /// this is a bit ugly, but let's try:
  /// rounds the number of particles such to the closes perfect cube
  long int roundedpart {lround(pow(lround(pow(npart,1./3)),3))};
  Particles.resize(roundedpart);
  InitAtoms();
}

Environment::~Environment()
{
}

void Environment::InitAtoms() 
{
  /// TODO: I need to implent cells:
  /*
    - create m-cells at a distance of r0 such that only
    interacting cells can talk to each other.
    - each cell can only contain one atom, so need make sure
    that the cell is empty before making a move.
  */
  /// TODO: At the moment this generates v and then adjusts it, maybe worth considering generating it directly the right way
  std::cout << " ---------------- Setting Initial Conditions for " << Particles.size() << " Atoms ...... \n";
  auto urngenerator = std::bind(urng,re);
  Vec3d TotalVelocity{0.0, 0.0, 0.0};
  double TotalSquaredVelocity{0.0};
  long int npart { lround(pow(Particles.size(), 1. / 3)) };
  double lcell {LCube/npart};

  size_t i{0};
  for(size_t ix{0}; ix < npart; ++ix){
    for(size_t iy{0}; iy < npart; ++iy){
      for(size_t iz{0}; iz < npart; ++iz){
        Particles[i].SetPosition(Vec3d{lcell * ix,lcell * iy,lcell * iz}, LCube);
        ++i;
      }
    }
  }

  for(auto& p: Particles){
    /// Init atoms in the cube
    /// give them some initial velocity as well,
    /// which in turn depends on the initial temp
    Vec3d Velocity{-1. + 2. * urngenerator(),
                   -1. + 2. * urngenerator(),
                   -1. + 2. * urngenerator()};
    
    p.SetVelocity(Velocity);
    TotalVelocity += Velocity;
    TotalSquaredVelocity += Velocity.Mod2();
  }
  BoostAndRescale(TotalSquaredVelocity, TotalVelocity/Particles.size());
  std::cout << PotentialEnergy() << " , " <<  KineticEnergy() << std::endl;
  std::cout << " ---------------- Total Momentum : " << TotalMomentum() << "\n";
  std::cout << " ---------------- Done \n";
}
  
void Environment::BoostAndRescale(const double& TotalSquaredVelocity, const Vec3d& AvgVelocity) 
{
  double CoMVelocity{TotalSquaredVelocity - AvgVelocity.Mod2()*Particles.size()};
  double VelocityRescale{sqrt(3. * Particles.size()*TargetTemperature / Constants::mass / CoMVelocity)};
  for(auto& p: Particles){
    p.SetVelocity((p.GetVelocity() - AvgVelocity)*VelocityRescale);
  }
}

double Environment::KineticEnergy() const
{
  double T{0.0};
  for(const auto& p: Particles){
    T += Constants::mass * p.GetVelocity().Mod2() / 2.;
  }
  return T;
}

double Environment::PotentialEnergy() const
{
  double U{0.0};
  double eta{pow(Constants::sigma / Constants::rCutoff, 12) -
             pow(Constants::sigma / Constants::rCutoff, 6)};
  /// for many particles this is a bit too expensive....
  for(size_t i{0}; i < Particles.size() - 1; ++i){
    for(size_t j{i+1}; j < Particles.size(); ++j){
      const double d{Atom::Distance(Particles[i],Particles[j],LCube)};
      if(d > Constants::rCutoff) continue;
      if(d == 0) {
        std::cout << " There is a bit of an issue here!" << std::endl;
        continue;
      }
      U += 4. * Constants::epsilon *
           (pow(Constants::sigma / d, 12) - pow(Constants::sigma / d, 6) - eta);
    }
  }
  return U;
}

double Environment::TotalEnergy() const
{
  double E{0.0}, T{0.0}, U{0.0};
  double eta{pow(Constants::sigma / Constants::rCutoff, 12) -
             pow(Constants::sigma / Constants::rCutoff, 6)};
  for (size_t i{0}; i < Particles.size(); ++i){
    T += Constants::mass * Particles[i].GetVelocity().Mod2() / 2.;
    if(i == Particles.size()-1) continue;
    for (size_t j{i + 1}; j < Particles.size(); ++j) {
      const double d{Atom::Distance(Particles[i], Particles[j], LCube)};
      if (d > Constants::rCutoff) continue;
      if (d == 0) continue;
      U += 4. * Constants::epsilon *
           (pow(Constants::sigma / d, 12) - pow(Constants::sigma / d, 6) - eta);
    }
  }
  return T+U;
}

double Environment::TotalMomentum() const
{
  Vec3d mom{0.0, 0.0, 0.0};
  for (const auto &p : Particles) {
    mom += (Constants::mass * p.GetVelocity());
  }
  return std::abs(mom.x) + std::abs(mom.y) + std::abs(mom.z);
}

int main(){
  Environment env(400, 4.0, 30.0);
  return 0;
}