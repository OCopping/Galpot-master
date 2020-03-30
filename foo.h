/**
 * gotta use stdint so we can use the same type in C++/C
 * could make the example a bit more complicated with some
 * #defines to fix but it's simpler this way.
 */
#include <stdint.h>

#ifdef __cplusplus
 
#include <string>
 
#endif
 
#ifdef __cplusplus
extern "C" 
{
#endif

  int    fa(int a, int *b);
  void fdisk(double R, double z, double *dPdR, double *dPdz);
  void fall(double R, double z, double *dPdR, double *dPdz);
  double fz(double R, double z);
  double fR(double R, double z);
  
#ifdef __cplusplus
}
#endif
