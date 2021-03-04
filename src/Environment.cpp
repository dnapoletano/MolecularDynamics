/// NOTE: measure if an algorith is better: sum acceptedtrials^2/ computing time
#include <cstdlib>
#include <fstream>

#include "Environment.hpp"

Environment::Environment(const size_t npart, const double Temperature, const double lcube)
: re{0}, urng{0.,1.}, TargetTemperature{Temperature}, LCube{lcube}
{
  /// this is a bit ugly, but let's try:
  /// rounds the number of particles such to the closes perfect cube
  long int roundedpart {lround(pow(lround(pow(npart,1./3)),3))};
  Particles.resize(roundedpart);
  Particles.SetLCube(lcube);
  std::cout << "Density : " << Particles.size()/Particles.Volume() << std::endl;
  InitAtoms();
}

Environment::~Environment()
{
}

void Environment::InitAtoms()
{
  std::cout << " ---------------- Setting Initial Conditions for " << Particles.size() << " Atoms ...... \n";
  InitiPositions();
  InitVelocities();
  Particles.CalcEnergy();
  std::cout << " E, K, U : " <<Particles.GetTotalEnergy() << ", " << Particles.GetKineticEnergy() << ", " << Particles.GetPotentialEnergy() << std::endl;
  std::cout << " ---------------- Total Momentum : " << Particles.GetTotalMomentum() << "\n";
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

void Environment::Evolve(const size_t maxsteps, const double retrialfactor = 1.0)
{
  Properties.Pressure = 0.0;
  Properties.PressureError = 0.0;
  double AvgOfSquares{0.0};
  size_t Accepted{0};
  for(size_t step{0}; step < maxsteps; ++step){
    Configuration TrialConfig{TrialStep(retrialfactor)};
    // std::cout << "#################### Step " << step << " #########################\n";
    // std::cout << "Trial : " << std::endl;
    // std::cout << "Energy : " << TrialConfig.GetTotalEnergy()
    //           << " K : " << TrialConfig.GetKineticEnergy()
    //           << " U : " << TrialConfig.GetPotentialEnergy()
    //           << " P : " << TrialConfig.GetPressure() << "\n";

    // std::cout << "Old : " << std::endl;
    // std::cout << "Energy : " << Particles.GetTotalEnergy()
    //           << " K : " << Particles.GetKineticEnergy()
    //           << " U : " << Particles.GetPotentialEnergy()
    //           << " P : " << Particles.GetPressure() << "\n";
    if(Accept(TrialConfig)) {
      // std::cout << " Accepted ! " << std::endl;
      /// I don't need trial config anymore!
      Particles = std::move(TrialConfig);
      ++Accepted;
      // std::cout << "New : " << std::endl;
      // std::cout << "Energy : " << Particles.GetTotalEnergy()
      //           << " K : " << Particles.GetKineticEnergy()
      //           << " U : " << Particles.GetPotentialEnergy()
      //           << " P : " << Particles.GetPressure() << "\n";
    }
    Properties.Pressure += Particles.GetPressure();
    AvgOfSquares += (Particles.GetPressure() * Particles.GetPressure());
    if(step % 100 == 0){
      std::cout << "  step : " << step << " Pressure : "
                << Properties.Pressure/step << std::flush << "\r";
    }
  }

  Properties.Pressure /= maxsteps;
  Properties.PressureError = sqrt((AvgOfSquares/maxsteps -
    (Properties.Pressure * Properties.Pressure)) / (maxsteps - 1));

  std::cout << "Accepted trials : " << Accepted << ", Pressure: " << Properties.Pressure
            << " +/- " << Properties.PressureError
            <<  ", V/N : " << Particles.Volume() / Particles.size() << "\n";

}

/*
  Copy the actual configuration in a new one.
  select a random particle of the new ensable and modify its position
  as
  r' = r + retrialfactor * Delta * (ran - 0.5)
  with Delta=1.e-3

  At the moment, no modifications of the velocities.
*/
Configuration Environment::TrialStep(const double retrialfactor = 1.0)
{
  Configuration new_config{Particles};
  /// select a random particle
  Configuration::iterator selected {
    (new_config.begin() + static_cast<int>(urng(re) * new_config.size()))
  };


  selected->IncrementPosition(Vec3d{
                                    retrialfactor * Constants::Delta * (urng(re) - 0.5),
                                    retrialfactor * Constants::Delta * (urng(re) - 0.5),
                                    retrialfactor * Constants::Delta * (urng(re) - 0.5)
                                   },LCube);
  new_config.CalcEnergy();
  return new_config;
}

bool Environment::Accept(const Configuration& trialconfig)
{
  if(trialconfig.GetPotentialEnergy() < Particles.GetPotentialEnergy()) return true;
  const double ratioNO{exp(Particles.LogPartitionFunction() - trialconfig.LogPartitionFunction())};
  return (urng(re) <= ratioNO);
}

void Environment::InitiPositions()
{
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
}

/// TODO: At the moment this generates v and then adjusts it, maybe worth considering generating it directly the right way
void Environment::InitVelocities()
{
  Vec3d TotalVelocity{0.0, 0.0, 0.0};
  double TotalSquaredVelocity{0.0};


  for(auto& p: Particles){
    /// Init atoms in the cube
    /// give them some initial velocity as well,
    /// which in turn depends on the initial temp
    Vec3d Velocity{-1. + 2. * urng(re),
                   -1. + 2. * urng(re),
                   -1. + 2. * urng(re)};

    p.SetVelocity(Velocity);
    TotalVelocity += Velocity;
    TotalSquaredVelocity += Velocity.Mod2();
  }
  BoostAndRescale(TotalSquaredVelocity, TotalVelocity/Particles.size());
}

int main(int argc, char* argv[]){
  double refactor {(argc > 1) ? std::atof(argv[1]) : 1.0};
  std::cout << " Actual Delta : " << refactor * Constants::Delta << std::endl;
  std::ofstream outfile;
  std::cout << std::left;
  outfile << "VOverN" << "," << "Pressure" << "," << "PressureError" <<"\n";
  std::vector<double> Temperatures {4.,4.*0.8,4.*pow(0.8,2),
    4.*pow(0.8,3),4.*pow(0.8,4),4.*pow(0.8,5)};
  std::vector<double> VOverN {logspace<double>(log10(0.5),log10(80),23,10.0)};
  //std::vector<double> VOverN {1.58455};
  for(const auto& t: Temperatures){
    outfile.open("results" + std::to_string(t) + ".txt",std::ios::out);
    outfile << std::left;
    /// Keep the number of particles fixed and only change the volume
    for(const auto& vn: VOverN){
      std::cout << " vn " << vn << std::endl;
      auto lcube {pow(vn * 400, 1./3.)};
      std::cout << " LCube : " << lcube << std::endl;
      Environment env(400, t, lcube);
      env.Evolve(1000, refactor);
      outfile << std::scientific << env.GetVOverN() << "," << env.GetObservables().Pressure
              << "," << env.GetObservables().PressureError << "\n";
    }
    outfile.close();
  }
  return 0;
}