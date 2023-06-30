/////////////////////////////////////////////////////////
//  Produces an output ASCII file in the HEPEVT format,
//  which can be used as input to the TextFileGen module
//  of LArSoft to create MCTruth data products that can
//  then be fed through largeant.
//
//  In HEPEVT format, each particle is represented in
//  two lines,
//
//  evtNum 1
//  1 pdgCode 0 0 0 0 Px Py Pz E Mass X Y Z 0
//
//  (Not sure what the 1s and 0s are for... probably
//  something fancy I'm not using.)
//
//  To generate a list, do:
//  > root -b -q 'ElectronGenList.cxx( <num_events>, <number-of-first-event>)'"
//
/////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMath.h"

// ================================================
void neutronGenList(int N = 50, int firstEvent = 1,std::string outFileName="./prodlist_neutrons_11_30_2020.txt")
{
  // Name of output ASCII file
  std::string outputFileName = outFileName;//"./prodlist_gammas_g10_500_a.txt";


  // Particle PDG (11 = electron, 22 = gamma, proton = 2212, neutron = 2112)
  int     pdg       = 2112;
  float   mass      = 0.9395654; //In GeV

  // Upper/lower bounds for generated particle energy [GeV]
  float   energyLo  = 0.005; // 5 MeV
  float   energyHi  = 1.060; // 1 GeV

  //Random smearing spread for points at top of detector
  float lx = 2.5; //x is really just a point on the edge of the detector
  float ly = 2.5; //y point range is roughly 20cm for hotspots, we'll try half. Originally 10
  float lz = 2.5; //z points range is roughly 40 cm for hotspots, we'll try half. Originally 20

  using coords = std::vector<vector<float>>;  //defining a vector of vector of floats for coordinates in uBoone
  //Modified g10 coordinates for now - only including points at top of detector
  coords cg10{{-.5,100.,30.},{-.5,100.,110.},{-.5,100.,230.},{-.5,100.,340.},{-.5,100.,455.},{-.5,100.,580.},{-.5,100.,680.},{-.5,100.,795.},{-.5,100.,910.},{-.5,100.,1000.},{250.5,100.,30.},{250.5,100.,110.},{250.5,100.,230.},{250.5,100.,340.},{250.5,100.,455.},{250.5,100.,580.},{250.5,100.,680.},{250.5,100.,795.},{250.5,100.,910.},{250.5,100.,1000.}};
  int numG10coords = cg10.size();


  // center point of uBoone TPC in G4 coordinates [cm], below, commented: Lariat
  //float x0 = 23.75; //for Lariat
  //float y0 = 0.00; // for Lariat
  //float z0 = 45.00; // for Lariat
  float x0 = 125.0;
  float y0 = 1.00;
  float z0 = 518.25;
  vector<float> randCoord;
  TVector3  tpc_center(x0,y0,z0);

  // span of dimensions of fiducial volume in which
  // to choose neutron start-point
  //float lx = 249.;
  //float ly = 232.;
  //float lz = 1036.;

  // -------------------------------------------------------
  // initialize random number generator
  // and open output file
  TDatime t1;
  int randSeed = t1.GetSecond();
  int y2 = t1.GetMinute();
  if(y2%2 == 1){
    randSeed = -1*randSeed;
  }
  cout << "Using random seed: " << randSeed << "." << endl;
  TRandom rand1(randSeed);
  ofstream output;
  output.open (outputFileName.c_str());

  // ----------------------------------
  // begin loop to create each particle
  for(int i = 0; i < N; i++){

    // assign random momentum, energy
    float p = energyLo + rand1.Rndm()*(energyHi-energyLo);
    float E = sqrt( pow(p,2) + pow(mass,2) );

    // pick starting location
    int randCoordIndex = rand1.Integer(numG10coords); //picks a random int up to the length of the g10 coordinate vector. This will be our starting point.
    randCoord = cg10[randCoordIndex]; //Random coordinate from our array of points
    x0 = randCoord[0];
    y0 = randCoord[1];
    z0 = randCoord[2];
    cout <<"Coord:("<< x0 << ","<<y0<<","<<z0<<")"<<endl;

    TVector3 loc_start(
      x0 + (2.*rand1.Rndm()-1.)*(lx/2.),
      y0 + (2.*rand1.Rndm()-1.)*(ly/2.),
      z0 + (2.*rand1.Rndm()-1.)*(lz/2.));

    // assign random isotropic direction vector
    //rand1.Rndm generates a number between 0,1 --> acos() returns 0 to pi from the vertical
    //double theta  = acos(1.-2.*rand1.Rndm());
    double theta  = acos(-1.*rand1.Rndm()); //trying to limit angle to pi/2 - pi, so downgoing I believe
    double phi    = 2.*TMath::Pi()*rand1.Rndm();
    TVector3 vec;
    vec.SetMagThetaPhi(p,theta,phi);

    // write it out to the file
    output << firstEvent+i << " " << 1 << "\n";
    output << 1 << "  " << pdg
                << "  " << 0
                << "  " << 0
                << "  " << 0
                << "  " << 0
                << "  " << vec.X()
                << "  " << vec.Y()
                << "  " << vec.Z()
                << "  " << E
                << "  " << mass
                << "  " << loc_start.X()
                << "  " << loc_start.Y()
                << "  " << loc_start.Z()
                << "  " << 0
                << "\n";

  }

  output.close();
}
