/*******************************************************************************
*                                                                              *
*  findOrbitMultiPot.cc                                                                *
*                                                                              *
*  C++ code written by Paul McMillan, 2017-                                    *
*  Lund Observatory, Lund University.                                          *
*  address: Box 43, SE-221 00 Lund, Sweden                                     *
*  e-mail:  paul@astro.lu.se                                                   *
*                                                                              *
*******************************************************************************/


#include "GalPot.h"
#include "KeplerPot.h"
#include "MiyamotoNagaiPot.h"
#include "MultiPot.h"
#include "OrbitIntegrator.h"

#include <ctime>

using std::cout;

int main(int argc,char *argv[])
{

 ifstream file;
 ofstream out;
 int nOut = 10; // default just in case
 if(argc<9) {
   cerr << "Input: R z phi v_R v_z v_phi output_file N_points_out\n";
   cerr << "   (distance in kpc, velocity in km/s)\n\n";
   cerr << "This is code which demonstrates the use of MultiPotential\n";
   cerr << "Users are encouraged to write their own code following this example\n";
   return 1;
 }


   Potential **PotList = new Potential* [3];

   string potfile = "pot/PJM16_best.Tpot";

   // Read potential from file
   file.open(potfile.c_str());
   if(!file) {
     cerr << "Input file does not exist. ";
     cerr << "Filename: " << potfile << "\n";
     return 1;
   }
   // GalaxyPotential from parameter file
   PotList[0] = new GalaxyPotential(file);
   file.close();
   // MiyamotoNagaiPotential of mass 8e8, a=0.1, b=0.02 - e.g. a nuclear disc
   PotList[1] = new MiyamotoNagaiPotential(8e8,0.1,0.02);
   // KeplerPotential (i.e. point mass) of mass 4e6, at r=0
   PotList[2] = new KeplerPotential(4e6);

  Potential *Phi = new MultiPotential(PotList,3);




 out.open(argv[7]);
 nOut = atoi(argv[8]);

 // Read position & velocity from input (In galactocentric cylindrical coordinates)
 Vector <double,6> XV=0.;

 // Convert input to code coordinates
 XV[0] = atof(argv[1]) * Units::kpc; // R
 XV[1] = atof(argv[2]) * Units::kpc; // z
 XV[2] = atof(argv[3]) * Units::degree; // phi
 XV[3] = atof(argv[4]) * Units::kms; // v_R
 XV[4] = atof(argv[5]) * Units::kms; // v_z
 XV[5] = atof(argv[6]) * Units::kms; // v_phi

 // Set up Integrator class
 OrbitIntegratorWithStats OI(XV, Phi, 10000.);

 Vector<double,6> *OrbOut = new Vector<double,6>[nOut];

 // Alternatively add :
 // double *tOut = new double[nout];
 //
 // Then run
 //
 //  int IntegrationFail = OI.runWithOutputIncludingTime(OrbOut,tOut,nOut);

 out << "#R    z     phi     v_R      v_z     v_phi\n";

 // run integration
 int IntegrationFail = OI.runWithOutput(OrbOut,nOut);

 // Output results
 if(IntegrationFail == 0) {
   cout << "Guiding Centre radius: " << OI.GuidingRadius << '\n' << std::flush;
   cout << "Minimum, Maximum Cylindrical radius: "<< OI.MinR << ','
	<< OI.MaxR << '\n'
	<< "Maximum z: " << OI.Maxz << '\n'
	<< "Minimum, Maximum Spherical radius: "<< OI.Minr << ','
	<< OI.Maxr << '\n'
	<< "Energy: " << OI.Energy/(Units::kms*Units::kms) << " km^2/s^2\n"
	<< "Angular Momentum (about symmetry axis): "
	<< OI.Lz/(Units::kms*Units::kpc)
	<< " kpc km/s\n"
       	<< "Mean Cylindrical radius: "<< OI.MeanR << '\n' << std::flush;
   for(int i=0;i!=nOut;i++) {
     out << OrbOut[i][0] / Units::kpc  << ' '
	 << OrbOut[i][1] / Units::kpc  << ' '
	 << OrbOut[i][2] / Units::degree  << ' '
	 << OrbOut[i][3] / Units::kms  << ' '
	 << OrbOut[i][4] / Units::kms  << ' '
	 << OrbOut[i][5] / Units::kms  << '\n' << std::flush;
   }
 } else {
   cout << "Input unbound in potential. Energy "
	<< OI.Energy/(Units::kms*Units::kms) << " km^2/s^2\n"
	<< std::flush;
 }

 out.close();
 return 0;

}
