/// NOTE: measure if an algorith is better: sum acceptedtrials^2/ computing time
#include <cstdlib>
#include <fstream>

#include "yaml-cpp/yaml.h"

#include "Environment.hpp"

Environment::Environment(const size_t npart, const double Temperature, const double lcube)
: re{0}, urng{0.,1.}, TargetTemperature{Temperature}, LCube{lcube}
{
  /// this is a bit ugly, but let's try:
  /// rounds the number of particles such to the closes perfect cube
  long int roundedpart {lround(pow(lround(pow(npart,1./3)),3))};
  Particles.resize(roundedpart);
  Particles.SetLCube(lcube);
  Delta = LCube / Particles.size();
  std::cout << utils::om::green << "Volume  : " << utils::om::reset
            << Particles.Volume() << std::endl;
  std::cout << utils::om::green << "LCube   : "  << utils::om::reset<< LCube << std::endl;
  std::cout << utils::om::green << "Delta   : "  << utils::om::reset<< Delta << std::endl;
  std::cout << utils::om::green << "Density : "
            << utils::om::reset<< Particles.size()/Particles.Volume() << std::endl;
  InitAtoms();
}

Environment::~Environment()
{
}

void Environment::InitAtoms()
{
  InitiPositions();
  InitVelocities();
  Particles.CalcEnergy();
}


void Environment::BoostAndRescale(const double& TotalSquaredVelocity, const Vec3d& AvgVelocity)
{
  double CoMVelocity{TotalSquaredVelocity - AvgVelocity.Mod2()*Particles.size()};
  double VelocityRescale{sqrt(3. * Particles.size()*TargetTemperature / Constants::mass / CoMVelocity)};
  for(auto& p: Particles){
    p.SetVelocity((p.GetVelocity() - AvgVelocity)*VelocityRescale);
  }
}

void Environment::Evolve(const size_t maxsteps)
{
  size_t multiplier {1};
  Properties.OldPressure = 0.0;
  while (true) {
    Properties.Pressure = 0.0;
    Properties.PressureError = 0.0;
    double AvgOfSquares{0.0};
    for (size_t step{0}; step < multiplier * maxsteps; ++step){
      bool rejected {true};
      while(rejected){
        rejected = (not Accept(TrialStep()));
      }

      Properties.Pressure += Particles.GetPressure();
      AvgOfSquares += (Particles.GetPressure() * Particles.GetPressure());
      if (step % 100 == 0 and step != 0) {
        std::cout << utils::om::red <<"  step : " << step << " Pressure : "
                  << Properties.Pressure / step << utils::om::reset << std::flush << "\r";
      }
    }
    const double totsteps = multiplier * maxsteps;
    Properties.Pressure /= (totsteps);
    AvgOfSquares /= totsteps;
    Properties.PressureError = sqrt(std::abs(AvgOfSquares -
                                             (Properties.Pressure * Properties.Pressure)) /
                                    (totsteps));
    if (Properties.PressureError / Properties.Pressure < 1.e-2 and
        std::abs(Properties.Pressure - Properties.OldPressure)/Properties.Pressure < 1.e-2){
          break;
    }
    std::cout << utils::om::underln << utils::om::yellow << __PRETTY_FUNCTION__
              << "{\n  "<< "Retry with " << maxsteps * multiplier << " steps.\n  "
              << "Pressure = " << Properties.Pressure << utils::om::pm
              << Properties.PressureError << "; OldPressure = " << Properties.OldPressure
              << "\n}" << utils::om::reset << "\n";
    multiplier *= 2;
    Properties.OldPressure = Properties.Pressure;
  }

  std::cout << "------------------------------------------------------------------- \n";
  std::cout << utils::om::red << std::setw(20) << "Pressure: " << Properties.Pressure
            << " +/- " << Properties.PressureError
            <<  ", V/N : " << Particles.Volume() / Particles.size() << utils::om::reset << "\n";
  std::cout << "------------------------------------------------------------------- \n";

}

/*
  Copy the actual configuration in a new one.
  select a random particle of the new ensable and modify its position
  as
  r' = r + Delta * (-1 + 2. * rand)
*/
Configuration Environment::TrialStep()
{
  Configuration new_config{Particles};
  // /// select a random particle
  // Configuration::iterator selected {
  //   (new_config.begin() + static_cast<int>(urng(re) * new_config.size()))
  // };

  // selected->IncrementPosition(Vec3d{
  //                                   Delta * (-1 + 2. * urng(re)),
  //                                   Delta * (-1 + 2. * urng(re)),
  //                                   Delta * (-1 + 2. * urng(re))
  //                                  },LCube);

  for(auto& p: new_config){
    p.IncrementPosition(Vec3d{
                            Delta * (-1 + 2. * urng(re)),
                            Delta * (-1 + 2. * urng(re)),
                            Delta * (-1 + 2. * urng(re))},
                        LCube);
  }

  new_config.CalcEnergy();
  return new_config;
}

bool Environment::Accept(const Configuration& trialconfig)
{
  const double ratioNO{exp( (Particles.GetTotalEnergy() - trialconfig.GetTotalEnergy()) / Particles.GetTemperature() )};
  /// if nothing has changed with the trial, retry until something happens
  if(ratioNO == 1.0) {
    Particles = trialconfig;
    return false;
  } else if (ratioNO > 1) {
    Particles = trialconfig;
    return true;
  } else {
    if (urng(re) < ratioNO) {
      Particles = trialconfig;
      return true;
    }
    return false;
  }
  return false;
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

struct Args {
  int VNSamples {23};
  long int seed {0};
  long int Iterations {100};
  long int NParticles {400};
  double Temperature {4.0};
  std::pair<double,double> VOverNLimits {2.5,80};
  std::string OutFileName {"results"};
};

void save_args(const std::string& filename, Args& args)
{
  YAML::Node file = YAML::LoadFile(filename);
  for(const auto& t : file){
    std::string name = t.first.as<std::string>();
    if(name=="VNSamples")
      args.VNSamples = t.second.as<int>();
    if(name=="Seed" )
      args.seed = t.second.as<long int>();
    if(name=="Iterations" )
      args.Iterations = t.second.as<long int>();
    if(name=="NParticles" )
      args.NParticles = t.second.as<long int>();
    if(name=="Temperature")
      args.Temperature = t.second.as<double>();
    if(name=="VNLimits")
      args.VOverNLimits = t.second.as<std::pair<double, double> >();
    if(name=="OutFileName")
      args.OutFileName = t.second.as<std::string>();
  }
}

int main(int argc, char* argv[]){
  Args args = Args();
  if(argc > 1){
    save_args({argv[1]},args);
  }
  std::ofstream outfile;
  std::vector<double> Temperatures {args.Temperature,4.*0.8,
                                    args.Temperature*pow(0.8,2),
                                    args.Temperature*pow(0.8,3),
                                    args.Temperature*pow(0.8,4),
                                    args.Temperature*pow(0.8,5)
                                   };
  std::vector<double> VOverN {(utils::logspace<double>(log10(args.VOverNLimits.first),
                               log10(args.VOverNLimits.first),args.VNSamples,10.0))};
  for(const auto& t: Temperatures){
    outfile.open(args.OutFileName + std::to_string(t) + ".txt",std::ios::out);
    outfile << std::left;
    outfile << "VOverN" << "," << "Pressure" << "," << "PressureError" <<"\n";
    /// Keep the number of particles fixed and only change the volume
    for(const auto& vn: VOverN){
      auto lcube {pow(vn * args.NParticles, 1./3.)};
      Environment env(args.NParticles, t, lcube);
      long int Iterations = args.Iterations;
      env.Evolve(Iterations);
      outfile << std::scientific << env.GetVOverN() << "," << env.GetObservables().Pressure
              << "," << env.GetObservables().PressureError << "\n";
    }
    outfile.close();
  }
  return 0;
}