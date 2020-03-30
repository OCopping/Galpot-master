#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdint.h>
 
#include "foo.h"
#include "GalPot.h"
 
using std::cout;
 
bool initialized = false;
GalaxyPotential *Phi = nullptr;
 
/*
gcc -O0 -g -Wall -pedantic  -c -o main.o main.c
gcc -lm -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -c -o one_plummer_with_stream.o one_plummer_with_stream.c
g++ -O3 -ffast-math -Isrc/ -c -o foo.o foo.cc  -Lobj -lPot -lOther -lm
g++ -O3 -ffast-math -Isrc/ -o a.out one_plummer_with_stream.o foo.o -Lobj -lPot -lOther -lm
*/
 
void init() {
  std::cout << "Sanity check : Initializing potential" << std::endl;
  ifstream file;
  char fname[]="pot/PJM17_0.Tpot";
  file.open(fname);
  Phi = new GalaxyPotential(file);
  file.close();
  initialized = true;
}

int fa(int a, int *b)
{
  std::cout << "About to change b\n";
  std::cout << "&b = " << &b << "\n";
  //std::cout << "b = " << b << "\n";
  //std::cout << "*b = " << *b << "\n";
  *b = 2;
  std::cout << "Changed b\n";
  /*double P,dPdR,dPdz;
 
  ifstream file;
  char fname[]="pot/PJM17_best_disk.Tpot";
  file.open(fname);
  GalaxyPotential Phi(file);
  file.close();
 
  P=Phi(10.,1.,dPdR,dPdz);
 
  std::cout << "Hi\n";
  std::cout << P*956044.576449833 << ' ' << dPdR*956044.576449833 << ' ' << dPdz*956044.576449833 << "\n";*/
 
return 2*a;
}
 
void fdisk(double R, double z, double *dPdR_ret, double *dPdz_ret)
{
  double P,dPdR,dPdz;
 
  if (!initialized)
    init();
 
  P=(*Phi)(R,z,dPdR,dPdz);
 
  *dPdR_ret = dPdR*956044.576449833;
  *dPdz_ret = dPdz*956044.576449833;
}
 
void fall(double R, double z, double *dPdR_ret, double *dPdz_ret)
{
  double P,dPdR,dPdz;
 
  if (!initialized)
    init();
 
  P=(*Phi)(R,z,dPdR,dPdz);
 
  *dPdR_ret = dPdR*956044.576449833;
  *dPdz_ret = dPdz*956044.576449833;
}
 
double fR(double R, double z)
{
  double P,dPdR,dPdz;
 
  if (!initialized)
    init();
 
  P=(*Phi)(R,z,dPdR,dPdz);
 
  return dPdR*956044.576449833;
}
 
double fz(double R, double z)
{
  double P,dPdR,dPdz;
  if (!initialized)
    init();
 
  P=(*Phi)(R,z,dPdR,dPdz);
 
  return dPdz*956044.576449833;
}
