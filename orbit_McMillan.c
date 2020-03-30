#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_erf.h>
#include "foo.h"
//#include <gsl/gsl_sf_gamma.h>

// Compiling instrucitons
/*
gcc -lm -lgsl -lgslcblas -O3 -I/opt/local/include -L/opt/local/lib -c -o orbit_McMillan.o orbit_McMillan.c
g++ -O3 -ffast-math -Isrc/ -c -o foo.o foo.cc  -Lobj -lPot -lOther -lm
g++ -O3 -ffast-math -Isrc/  -lm -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -o a.out orbit_McMillan.o foo.o -Lobj -lPot -lOther
 */

//note that only first particle is massive
int N = 1; //stores the current number of particles
double M_plum = 0.000002;//0.025;
double M_plum0;
double sigma_plum = 0.01;//0.85;

double sigma_strip = 0.05;

double sigma_vel = 0.5;

double lambda_vlag=1.;

double lambda_vfrac=0.5;

double lambda_fM=0.;

#define NMAX_POS_PROG 150000

//double pos_prog[3*NMAX_POS_PROG]; // for storing pos_prog which will be interpolated later
//double pos_lmc[3*NMAX_POS_PROG];
//double vel_prog[3*NMAX_POS_PROG];
int num_pos_prog=0;

#define NMAX 100000

#define NMAX_r_min 100 // max number of rmins

#define eta 0.01
#define dt_max 1.0
#define dt_min 0.0005

#define N_PER_PERI 200
double SIGMA_T; //0.09 // 0.01

#define G (43007.105731706317) // in kpc (km/s)^2/1e10 Msun

#define POT
#define GENERATE_STREAM

#define GALPOT_ALL

#define INCLUDE_LMC
double M_LMC=15.;
double rs_LMC=15.;
double x_LMC[3] = {-1.,-41.,-28.};
double v_LMC[3] = {-57.,-226.,221.};

/*#define STATICBOWDENROT
double rho0_NFW=0.00326203340141;
double rs_NFW=9.97033771405;
double rho1_y20=0.000278435413775;
double r1_y20=11.680832623;
double l_rot=10.;
double b_rot=10.;
double M_rot[3][3],M_rot_inv[3][3];*/

/*#define STATICLM2010
double vhalo=121.9;
double phi_LM=97.*M_PI/180.;
double q1_LM=1.38;
double qz_LM=1.36;
double rhalo=12.;
double c1;
double c2;
double c3;*/

//#define STATICNFW
double M_NFW=137.;
double rs_NFW=19.6;//0130718954;
double c_NFW=15.4;
double q_NFW=1.;
double theta_NFW=0.;
double phi_NFW=M_PI/2.;

/*#define STATICLOG
double rs_LOG=12.;
double vc_LOG=121.9;*/

//#define STATICMIYAMOTONAGAI
int n_MN = 1;
double a_MN[6]={3., 5.321702747134218, 4.662815777920417,9.021, 9.143, 7.758};
double b_MN[6]={0.28, 0.46986069075175774, 0.46986069075175774, 0.168, 0.168, 0.168};
double M_MN[6]={6.8,-26.436898154682613,30.947757591582164,2.046,2.169,-3.049};//6.8;//10.;

//#define STATICHQDE
double M_HQ=0.5;//0.5;//3.4;
double a_HQ=0.5;//0.5;//0.7;

/*#define STATICBPL
double BPL_alpha=1.8;
double BPL_M=0.5;
double BPL_r=1.9;*/

struct pdata_struct
{
  double x[3];
  double v[3];
  double t_strip;
  double phase;
  double theta;
  double theta_phase;
};

void initialization_steps();
void kick(struct pdata_struct *pdata, double *force, double dt);
void drift(struct pdata_struct *pdata, double dt);
void find_external_force(struct pdata_struct *pdata, double *force, double t);
void find_external_force_forwards(struct pdata_struct *pdata, double *force, double t, double *pos_lmc_t, double *pos_mw_t);
void find_nbody_force(struct pdata_struct *pdata, double *force,double *pos_prog_t);
double find_nbody_force_w_dt(struct pdata_struct *pdata, double *force, double *pos_prog_t);

void update_prog_pos(double *pos_prog_t, double t, double *pos_prog);
void update_prog_vel(double *vel_prog_t, double t, double *vel_prog);
void update_lmc_pos(double *pos_lmc_t, double t, double *pos_lmc);
void update_mw_pos(double *pos_mw_t, double t, double *pos_mw);
void update_mw_vel(double *vel_mw_t, double t, double *vel_mw);
void update_r_prog_lmc(double *r_prog_lmc_t, double t, double *r_prog_lmc);

void generate_stream(struct pdata_struct *pdata, double *force, double t, double *pos_prog_t, double *vel_prog_t, double *pos_mw_t, double *vel_mw_t, int leading_or_trailing);

double dist(double *x, double *y);
double mag(double *x);
double min2(double x, double y);
double max2(double x, double y);

double gen_normal(double mean, double sigma);

double r_lag_calc(double *pos_prog_t, double *vel_prog_t, double *pos_mw_t, double *vel_mw_t);
double d2phidz2_MN(double R, double z);
double d2phidRdz_MN(double R, double z);
double d2phidR2_MN(double R, double z);
double d2phidr2_HQ(double r);
double dphidr_NFW(double r);
double d2phidr2_NFW(double r);

double d2phidR2_NFW(double R, double z);
double d2phidRdz_NFW(double R, double z);
double d2phidz2_NFW(double R, double z);

double d2phidr2(double R, double z);

double pot_HQ(double r);
double pot_NFW(double r);
double pot_MN(double R, double z);
double pot_total(struct pdata_struct *pdata, int counter);

double E_rel_calc(struct pdata_struct *pdata, int counter, double *pos_prog_t, double *vel_prog_t);

int main(int argc, char *argv[])
{

  srand(0746374);

  if(argc!=14)
    {
      printf("Correct usage: ./a.out M_LMC rs_LMC x_LMC y_LMC z_LMC vx_LMC vy_LMC vz_LMC tmax rs_NFW M_NFW c_NFW pid\n");
      exit(1);
    }

  int pid = atoi(argv[13]);

  rs_NFW = atof(argv[10]);
  M_NFW = atof(argv[11]);
  c_NFW = atof(argv[12]);

  M_LMC = atof(argv[1]);
  rs_LMC = atof(argv[2]);
  x_LMC[0] = atof(argv[3]);
  x_LMC[1] = atof(argv[4]);
  x_LMC[2] = atof(argv[5]);
  v_LMC[0] = atof(argv[6]);
  v_LMC[1] = atof(argv[7]);
  v_LMC[2] = atof(argv[8]);
  M_plum0 = 0.001;

  //dt_strip = dt_strip*pow(M_plum/(2.e-6),-1./3.);

  FILE *fileout;
  char fname[999];
  sprintf(fname,"orbit_%d.txt",pid);
  fileout = fopen(fname,"w");

  int steps=0;
  int axes;
  int i,j;
  double t = 0.0;
  double dt;
  struct pdata_struct *pdata = (struct pdata_struct*)malloc(NMAX*sizeof(struct pdata_struct));
  if(pdata == NULL)
    {
      printf("Can't assign memory for pdata\n");
      return -1;
    }

  double *force = malloc(3*NMAX*sizeof(double));

  if(force == NULL)
    {
      printf("Can't assign memory for force\n");
      return -1;
    }

  double *pos_prog = malloc(3*NMAX_POS_PROG*sizeof(double));
  
  if(pos_prog == NULL)
    {
      printf("Can't assign memory for pos_prog\n");
      return -1;
    }

  double *vel_prog = malloc(3*NMAX_POS_PROG*sizeof(double));
  
  if(vel_prog == NULL)
    {
      printf("Can't assign memory for vel_prog\n");
      return -1;
    }

  double *pos_lmc = malloc(3*NMAX_POS_PROG*sizeof(double));
  
  if(pos_lmc == NULL)
    {
      printf("Can't assign memory for pos_lmc\n");
      return -1;
    }

  double *pos_mw = malloc(3*NMAX_POS_PROG*sizeof(double));
  
  if(pos_mw == NULL)
    {
      printf("Can't assign memory for pos_mw\n");
      return -1;
    }

  double *vel_mw = malloc(3*NMAX_POS_PROG*sizeof(double));
  
  if(vel_mw == NULL)
    {
      printf("Can't assign memory for vel_mw\n");
      return -1;
    }

  double *r_prog_lmc = malloc(NMAX_POS_PROG*sizeof(double));
  
  if(r_prog_lmc == NULL)
    {
      printf("Can't assign memory for r_prog_lmc\n");
      return -1;
    }

  double t_r_min[NMAX_r_min];
  int N_r_min=0;
  double r,r_prev,r_prev_2;
  double delta_t_r_min;
  r = r_prev = r_prev_2 = 0.;

  double m = 0.05;
  double sigma = 0.5;
  double xcom[3] = {-8.1196288 ,   0.24419839,  16.92195489}; //{atof(argv[1]), atof(argv[2]), atof(argv[3])};//{-8.1196288 ,   0.24419839,  16.92195489};//-8.1196288 ,  -0.24419839,  16.92195489};
    double vcom[3] = {44.2364676,-109.87367822,-15.846060985767631}; //{atof(argv[4]), atof(argv[5]), atof(argv[6])};
  //double xcom[3] = {-8.15660526,0.23537299,16.672783816552656};//{-8.1196288 ,   0.24419839,  16.92195489};//-8.1196288 ,  -0.24419839,  16.92195489};
  //double vcom[3] = {44.2364676,-109.87367822,-15.846060985767631};//{44.446674611739233, -120.95723526722537, -15.182782258924757};
  double tmax = atof(argv[9])*1.0227; //5.*1.0227*pow(M_plum/(2.e-6),-1./3.);//5.1135;//7.1589;//5.1135;
  //tmax = ceil(tmax/0.00010227)*0.00010227;
  
  double dv_stream = 1.;
  int N_insert=0;
  double N_insert_double = 0.;

  double t_last_insert = 0.;

  double r_part;
  double t_inject;
  int inject_counter,peri_counter;

  double *pos_prog_t = malloc(3*sizeof(double));
  double *vel_prog_t = malloc(3*sizeof(double));
  double *pos_lmc_t = malloc(3*sizeof(double));
  double *pos_mw_t = malloc(3*sizeof(double));
  double *vel_mw_t = malloc(3*sizeof(double));
  double *r_prog_lmc_t = malloc(sizeof(double));
  
  int sign;
  
  double theta,theta_part;

  double E_rel;

  int include_particle;

  double r_min_LMC,t_rmin_LMC;
  double r_max_LMC = 0.;
  double r_LMC;

  
  double r_max,r_min;
  
  initialization_steps();

  /*****************************************************************/
  /*             FIRST INTEGRATE THE ORBIT BACKWARDS               */
  /*****************************************************************/
    
  N+=1; // for LMC
  N+=1; // for MW

  for(i=0;i<3;i++)
    {
      pdata[0].x[i] = xcom[i];
      pdata[0].v[i] = -vcom[i]; //integrating backwards!
      pdata[1].x[i] = x_LMC[i];
      pdata[1].v[i] = -v_LMC[i];
      pdata[2].x[i] = 0.;
      pdata[2].v[i] = 0.;
    }

  find_external_force(pdata,force,t);

  pdata[0].t_strip = 0.;
  pdata[0].phase = 0.;
  pdata[0].theta = 0.;
  pdata[0].theta_phase = 0.;

  for(axes=0;axes<3;axes++)
    {
      pos_prog[axes]= pdata[0].x[axes];
      vel_prog[axes]=-pdata[0].v[axes];
      pos_lmc[axes] = pdata[1].x[axes];
      pos_mw[axes]  = pdata[2].x[axes];
      vel_mw[axes]  =-pdata[2].v[axes];
    }
  r_prog_lmc[0] = dist(pdata[0].x,pdata[1].x);
  
  num_pos_prog++;

  //fprintf(fileout,"%d %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",steps,tmax-t,pdata[0].x[0],pdata[0].x[1],pdata[0].x[2],pdata[0].v[0],pdata[0].v[1],pdata[0].v[2],pdata[1].x[0],pdata[1].x[1],pdata[1].x[2],pdata[1].v[0],pdata[1].v[1],pdata[1].v[2],pdata[2].x[0],pdata[2].x[1],pdata[2].x[2],pdata[2].v[0],pdata[2].v[1],pdata[2].v[2]);

  while(steps < (int)(tmax/0.00010227)) //t < tmax)
    {
      find_external_force(pdata,force,t);
      r_prev = dist(pdata[1].x,pdata[2].x);
            
      dt = 0.0001*1.0227; // time steps are 100000 years. For this to work, the total time must be an integer multiple of 100000 years.

      if(t+dt >= tmax)
	dt = tmax-t;
      
      //KDK
      /******************/
      //kick
      kick(pdata,force,dt/2.);
      //drift
      drift(pdata,dt);
      //kick
      t+=dt;
      r = dist(pdata[1].x,pdata[2].x);
      if(r > r_prev && r_prev < r_prev_2)
	{
	  if(N_r_min + 1 >= NMAX_r_min)
	    {
	      printf("Found %d minimia but only space for %d\n",N_r_min+1,NMAX_r_min);
	      exit(0);
	    }
	  
	  t_r_min[N_r_min]=tmax-t-dt;
	  N_r_min++;
	}


      r_LMC = dist(pdata[1].x,pdata[2].x);
      if(r_LMC > r_max_LMC)
	r_max_LMC = r_LMC;
      
      // kick again!
      find_external_force(pdata,force,t);
      kick(pdata,force,dt/2.);

      r_prev_2 = r_prev;

      

      steps++;

      //fprintf(fileout,"%d %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",steps,tmax-t,pdata[0].x[0],pdata[0].x[1],pdata[0].x[2],pdata[0].v[0],pdata[0].v[1],pdata[0].v[2],pdata[1].x[0],pdata[1].x[1],pdata[1].x[2],pdata[1].v[0],pdata[1].v[1],pdata[1].v[2],pdata[2].x[0],pdata[2].x[1],pdata[2].x[2],pdata[2].v[0],pdata[2].v[1],pdata[2].v[2]);

      
      for(axes=0;axes<3;axes++)
	{
	  pos_prog[3*steps+axes] = pdata[0].x[axes];
	  vel_prog[3*steps+axes] = -pdata[0].v[axes];
	  pos_lmc[3*steps+axes]  = pdata[1].x[axes];
	  pos_mw[3*steps+axes]   = pdata[2].x[axes];
	  vel_mw[3*steps+axes]   = -pdata[2].v[axes];
	}
      r_prog_lmc[steps] = dist(pdata[0].x,pdata[1].x);
      num_pos_prog++;
    }


  fprintf(fileout,"%d %f\n",N_r_min,r_max_LMC);
  
  fclose(fileout);

  if(steps > NMAX_POS_PROG)
    {
      printf("Critical error: took %d steps but NMAX_POS_PROG is %d\n",steps,NMAX_POS_PROG);
      exit(1);
    }

  free(pos_prog_t);
  free(vel_prog_t);
  free(pos_lmc_t);
  free(pos_mw_t);
  free(vel_mw_t);
  free(pdata);
  free(force);
  free(pos_prog);
  free(pos_lmc);
  free(vel_prog);
  free(pos_mw);
  free(vel_mw);
  
  return 0;
}

void initialization_steps()
{
#ifdef STATICLM2010
  c1=pow(cos(phi_LM)/q1_LM,2.)+pow(sin(phi_LM),2.);
  c2=pow(cos(phi_LM),2.)+pow(sin(phi_LM)/q1_LM,2.);
  c3=2.*sin(phi_LM)*cos(phi_LM)*(pow(q1_LM,-2.)-1.);
#endif


#ifdef STATICBOWDENROT

  double theta_rot = l_rot*M_PI/180.;
  double phi_rot = M_PI/2.-b_rot*M_PI/180.;

  

  M_rot[0][0] = cos(theta_rot)*cos(phi_rot);
  M_rot[0][1] = -sin(theta_rot);
  M_rot[0][2] = cos(theta_rot)*sin(phi_rot);
  M_rot[1][0] = sin(theta_rot)*cos(phi_rot);
  M_rot[1][1] = cos(theta_rot);
  M_rot[1][2] = sin(theta_rot)*sin(phi_rot);
  M_rot[2][0] = -sin(phi_rot);
  M_rot[2][1] = 0.;
  M_rot[2][2] = cos(phi_rot);

  M_rot_inv[0][0] = cos(phi_rot)*cos(theta_rot);
  M_rot_inv[0][1] = cos(phi_rot)*sin(theta_rot);
  M_rot_inv[0][2] = -sin(phi_rot);
  M_rot_inv[1][0] = -sin(theta_rot);
  M_rot_inv[1][1] = cos(theta_rot);
  M_rot_inv[1][2] = 0.;
  M_rot_inv[2][0] = sin(phi_rot)*cos(theta_rot);
  M_rot_inv[2][1] = sin(phi_rot)*sin(theta_rot);
  M_rot_inv[2][2] = cos(phi_rot);

#endif //STATICBOWDENROT
}
double dist(double *x, double *y)
{
  return sqrt(pow(x[0]-y[0],2.0)+pow(x[1]-y[1],2.0)+pow(x[2]-y[2],2.0));
}

double mag(double *x)
{
  return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

double min2(double x, double y)
{
  if(x < y)
    return x;
  else
    return y;
}

double max2(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}


void kick(struct pdata_struct *pdata,double *force, double dt)
{
  int i,j;
  for(i=0;i<N;i++)
    for(j=0;j<3;j++)
      pdata[i].v[j]+=force[i*3+j]*dt;

}

void drift(struct pdata_struct *pdata, double dt)
{
  int i,j;
  for(i=0;i<N;i++)
    for(j=0;j<3;j++)
      pdata[i].x[j]+=pdata[i].v[j]*dt;
}

void find_external_force(struct pdata_struct *pdata, double *force, double t)
{
  int i,axes1,axes2;
  double dx,dy,dz,r;
  double mdpsiNFWdr, mdpsiy20dr;
  double a_inertial[3], a_rot[3];
  double denom;

  double dR_mn;
  double rho,rho_0,lambda,v_rel,r_vmax,v_max,sigma,X_vel,rho_mn;
  int MN_counter;
  

#ifdef STATICBPL
  double gamma_store;
  double gamma_inc_store;
#endif
  
#ifdef STATICNFW
  double nhatdotr;
#endif //STATICNFW

#if defined GALPOT_DISK || defined GALPOT_ALL
  double dPdR,dPdz,R;
#endif //GALPOT_DISK

  for(i=0;i<3*N;i++)
    force[i]=0.;

  for(i=0;i<N;i++)
    {
      
      dx = pdata[i].x[0]-pdata[2].x[0];
      dy = pdata[i].x[1]-pdata[2].x[1];
      dz = pdata[i].x[2]-pdata[2].x[2];
      if(i!=2)
	{
#ifdef STATICMIYAMOTONAGAI
	  dx = pdata[i].x[0]-pdata[2].x[0];
	  dy = pdata[i].x[1]-pdata[2].x[1];
	  dz = pdata[i].x[2]-pdata[2].x[2];
	  
	  for(MN_counter=0;MN_counter<n_MN;MN_counter++)
	    {
	      
	      denom = pow(dx*dx+dy*dy+pow(a_MN[MN_counter]+sqrt(dz*dz+b_MN[MN_counter]*b_MN[MN_counter]),2.0),1.5);
	      
	      force[3*i] -= G*M_MN[MN_counter]*dx/denom;
	      force[3*i+1] -= G*M_MN[MN_counter]*dy/denom;
	      force[3*i+2] -= G*M_MN[MN_counter]*dz/denom*(a_MN[MN_counter]+sqrt(dz*dz+b_MN[MN_counter]*b_MN[MN_counter]))/sqrt(dz*dz+b_MN[MN_counter]*b_MN[MN_counter]);
	    }

#endif //STATICMIYAMOTONAGAI

#ifdef GALPOT_ALL

      // dPdR and dPdz are positive
      
      R = sqrt(dx*dx+dy*dy);
      fall(R,dz,&dPdR,&dPdz);

      force[3*i]   += -dPdR*dx/R;
      force[3*i+1] += -dPdR*dy/R;
      force[3*i+2] += -dPdz;
      
#endif //GALPOT_ALL
	  
#ifdef STATICHQDE
	  r = sqrt(dx*dx+dy*dy+dz*dz);
	  force[3*i] -= G*M_HQ/pow(a_HQ+r,2.)*dx/r;
	  force[3*i+1] -= G*M_HQ/pow(a_HQ+r,2.)*dy/r;
	  force[3*i+2] -= G*M_HQ/pow(a_HQ+r,2.)*dz/r;
	  
#endif //STATICHQDE
	  
#ifdef STATICNFW
	  nhatdotr = dx*cos(phi_NFW)*cos(theta_NFW)+dy*cos(phi_NFW)*sin(theta_NFW)+dz*sin(phi_NFW);
	  r = sqrt(dx*dx+dy*dy+dz*dz+(pow(q_NFW,-2.)-1.)*pow(nhatdotr,2.));
	  force[3*i]   += -G*M_NFW*1./pow(r,3.)*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW))*(log(1.+r/rs_NFW)-r/(r+rs_NFW))*(dx+(pow(q_NFW,-2.)-1.)*nhatdotr*cos(phi_NFW)*cos(theta_NFW));
	  force[3*i+1] += -G*M_NFW*1./pow(r,3.)*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW))*(log(1.+r/rs_NFW)-r/(r+rs_NFW))*(dy+(pow(q_NFW,-2.)-1.)*nhatdotr*cos(phi_NFW)*sin(theta_NFW));
	  force[3*i+2] += -G*M_NFW*1./pow(r,3.)*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW))*(log(1.+r/rs_NFW)-r/(r+rs_NFW))*(dz+(pow(q_NFW,-2.)-1.)*nhatdotr*sin(phi_NFW));
#endif
	  
	  
#ifdef STATICLOG
	  r = sqrt(dx*dx+dy*dy+dz*dz);
	  force[3*i] += -pow(vc_LOG,2.)/(r*r+rs_LOG*rs_LOG)*2.*dx;
	  force[3*i+1] += -pow(vc_LOG,2.)/(r*r+rs_LOG*rs_LOG)*2.*dy;
	  force[3*i+2] += -pow(vc_LOG,2.)/(r*r+rs_LOG*rs_LOG)*2.*dz;
#endif
	  
#ifdef STATICLM2010
	  denom = c1*dx*dx+c2*dy*dy+c3*dx*dy+pow(dz/qz_LM,2.)+pow(rhalo,2.);
	  force[3*i] += -pow(vhalo,2.)*(2.*c1*dx+c3*dy)/denom;
	  force[3*i+1] += -pow(vhalo,2.)*(2.*c2*dy+c3*dx)/denom;
	  force[3*i+2] += -pow(vhalo,2.)*(2.*dz/pow(qz_LM,2.))/denom;
	  
#endif //STATICLM2010
	  
	  
#ifdef STATICBPL
	  r = sqrt(dx*dx+dy*dy+dz*dz);
	  gamma_store = gsl_sf_gamma(3-BPL_alpha);
	  gamma_inc_store = gsl_sf_gamma_inc(3-BPL_alpha,r/BPL_r);
	  force[3*i] += -G*BPL_M*(gamma_store-gamma_inc_store)/gamma_store*dx/pow(r,3.);
	  force[3*i+1] += -G*BPL_M*(gamma_store-gamma_inc_store)/gamma_store*dy/pow(r,3.);
	  force[3*i+2] += -G*BPL_M*(gamma_store-gamma_inc_store)/gamma_store*dz/pow(r,3.);
	  //force[3*i+2] += -G*BPL_M*(gsl_sf_gamma(3-BPL_alpha)-gsl_sf_gamma_inc(3-BPL_alpha,r/BPL_r))/gsl_sf_gamma(3-BPL_alpha)*dz/pow(r,3.);
#endif 
	  
#ifdef STATICBOWDENROT
	  dx = 0.;
	  for(axes1=0;axes1<3;axes1++)
	    dx += M_rot_inv[0][axes1]*(pdata[i].x[axes1]-pdata[2].x[axes1]);
	  dy = 0.;
	  for(axes1=0;axes1<3;axes1++)
	    dy += M_rot_inv[1][axes1]*(pdata[i].x[axes1]-pdata[2].x[axes1]);
	  dz = 0.;
	  for(axes1=0;axes1<3;axes1++)
	    dz += M_rot_inv[2][axes1]*(pdata[i].x[axes1]-pdata[2].x[axes1]);      
	  
	  r = sqrt(dx*dx+dy*dy+dz*dz);
	  mdpsiNFWdr = -4.*M_PI*G*rho0_NFW*pow(rs_NFW,3.)/pow(r,2.)*log(1.+r/rs_NFW)+4.*M_PI*G*rho0_NFW*pow(rs_NFW,2.)/r*1./(1.+r/rs_NFW);
	  mdpsiy20dr = -4.*M_PI*G*rho1_y20*pow(r1_y20,3.)/pow(r+r1_y20,2.)*(1.5*pow(dz/r,2.)-0.5)+8.*M_PI*G*rho1_y20*pow(r1_y20,3.)*r/pow(r+r1_y20,3.)*(1.5*pow(dz/r,2.)-0.5)+12.*M_PI*G*rho1_y20*pow(r1_y20,3.)*r/pow(r+r1_y20,2.)*pow(dz,2.)/pow(r,3.);
	  
	  a_inertial[0] = dx/r*(mdpsiNFWdr+mdpsiy20dr);
	  a_inertial[1] = dy/r*(mdpsiNFWdr+mdpsiy20dr);
	  a_inertial[2] = dz/r*(mdpsiNFWdr+mdpsiy20dr) - 12.*M_PI*G*rho1_y20*pow(r1_y20,3.)*r/pow(r+r1_y20,2.)*dz/pow(r,2.);
	  
	  for(axes1=0;axes1<3;axes1++)
	    a_rot[axes1] = 0.;
	  
	  for(axes1=0;axes1<3;axes1++)
	    for(axes2=0;axes2<3;axes2++)
	      a_rot[axes1] += M_rot[axes1][axes2]*a_inertial[axes2];
	  
	  for(axes1=0;axes1<3;axes1++)
	    for(axes2=0;axes2<3;axes2++)
	      force[3*i+axes1] += M_rot[axes1][axes2]*a_inertial[axes2];
#endif //STATICBOWDENROT
	}

#ifdef INCLUDE_LMC
      if(i!=1)
	{
	  dx = pdata[i].x[0]-pdata[1].x[0];
	  dy = pdata[i].x[1]-pdata[1].x[1];
	  dz = pdata[i].x[2]-pdata[1].x[2];
	  r = sqrt(dx*dx+dy*dy+dz*dz);
	  
	  force[3*i] -= G*M_LMC/pow(r+rs_LMC,2.)*dx/r;
	  force[3*i+1] -= G*M_LMC/pow(r+rs_LMC,2.)*dy/r;
	  force[3*i+2] -= G*M_LMC/pow(r+rs_LMC,2.)*dz/r;
	}
    }
      
  dx = pdata[1].x[0]-pdata[2].x[0];
  dy = pdata[1].x[1]-pdata[2].x[1];
  dz = pdata[1].x[2]-pdata[2].x[2];
  r = sqrt(dx*dx+dy*dy+dz*dz);
      
      
      
  rho_0 = M_NFW/(4.*M_PI*pow(rs_NFW,3.))*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW));
  //rho0_NFW;//M_NFW/(4.*M_PI*pow(rs_NFW,3.))*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW));                                                    
  
  rho = rho_0/(r/rs_NFW*pow(1.+r/rs_NFW,2.));
  
  if(rs_LMC > 8.)
    lambda = r/(2.2*rs_LMC-14.);
  else
    lambda = r/(0.45*rs_LMC);
  
  //lambda = r/(0.5*1.4*rs_LMC);
  
  v_rel = dist(pdata[1].v,pdata[2].v);
  r_vmax = 2.16258*rs_NFW;
  v_max = sqrt(G*M_NFW/r_vmax*(log(1.+r_vmax/rs_NFW)-r_vmax/(r_vmax+rs_NFW))/(log(1.+c_NFW)-c_NFW/(1.+c_NFW)));
  //v_max = sqrt(G*rho0_NFW*4.*M_PI*pow(rs_NFW,3.)/r_vmax*(log(1.+r_vmax/rs_NFW)-r_vmax/(r_vmax+rs_NFW)));
  sigma = v_max*1.4393*pow(r/rs_NFW,0.354)/(1.+1.1756*pow(r/rs_NFW,0.725));
  X_vel = v_rel/(sqrt(2.)*sigma);
  
  dR_mn = sqrt(dx*dx+dy*dy);
  rho_mn = 0.;
  for(MN_counter = 0; MN_counter<n_MN; MN_counter++)
    {
      
      rho_mn += (b_MN[MN_counter]*b_MN[MN_counter]*M_MN[MN_counter])/(4.*M_PI)*(a_MN[MN_counter]*dR_mn*dR_mn+(a_MN[MN_counter] + 3.*sqrt(dz*dz+b_MN[MN_counter]*b_MN[MN_counter]))*pow(a_MN[MN_counter] + sqrt(dz*dz+b_MN[MN_counter]*b_MN[MN_counter]),2.))/(pow(dR_mn*dR_mn+pow(a_MN[MN_counter] + sqrt(dz*dz + b_MN[MN_counter]*b_MN[MN_counter]),2.),2.5)*pow(dz*dz+b_MN[MN_counter]*b_MN[MN_counter],1.5));
    }
  //printf("lambda = %f\n",lambda);                                                                                                         
  //printf("gsl_sf_erf(X_vel) = %f\n",gsl_sf_erf(X_vel));          
  
  
  
  if(lambda >= 1.)
    for(i=0;i<3;i++)
      force[3+i]+=4.*M_PI*pow(G,2.)*(M_LMC)*(rho+rho_mn)*log(lambda)/pow(v_rel,3.)*(gsl_sf_erf(X_vel)-2.*X_vel/sqrt(M_PI)*exp(-pow(X_vel,2)))*(pdata[1].v[i]-pdata[2].v[i]);
  


#endif

}

void find_external_force_forwards(struct pdata_struct *pdata, double *force, double t, double *pos_lmc_t, double *pos_mw_t)
{
  int i,axes1,axes2;
  double dx,dy,dz,r;
  double mdpsiNFWdr, mdpsiy20dr;
  double a_inertial[3], a_rot[3];
  double denom;
  int MN_counter;

#ifdef STATICBPL
  double gamma_store;
  double gamma_inc_store;
#endif
  
#ifdef STATICNFW
  double nhatdotr;
#endif //STATICNFW

#if defined GALPOT_DISK || defined GALPOT_ALL
  double dPdR,dPdz,R;
#endif //GALPOT_DISK

  for(i=0;i<3*N;i++)
    force[i]=0.;

  for(i=0;i<N;i++)
    {
      dx = pdata[i].x[0]-pos_mw_t[0];
      dy = pdata[i].x[1]-pos_mw_t[1];
      dz = pdata[i].x[2]-pos_mw_t[2];
#ifdef STATICMIYAMOTONAGAI
      dx = pdata[i].x[0]-pos_mw_t[0];
      dy = pdata[i].x[1]-pos_mw_t[1];
      dz = pdata[i].x[2]-pos_mw_t[2];

      for(MN_counter=0;MN_counter<n_MN;MN_counter++)
	{
	  denom = pow(dx*dx+dy*dy+pow(a_MN[MN_counter]+sqrt(dz*dz+b_MN[MN_counter]*b_MN[MN_counter]),2.0),1.5);
	  
	  force[3*i] -= G*M_MN[MN_counter]*dx/denom;
	  force[3*i+1] -= G*M_MN[MN_counter]*dy/denom;
	  force[3*i+2] -= G*M_MN[MN_counter]*dz/denom*(a_MN[MN_counter]+sqrt(dz*dz+b_MN[MN_counter]*b_MN[MN_counter]))/sqrt(dz*dz+b_MN[MN_counter]*b_MN[MN_counter]);
	}

#endif //STATICMIYAMOTONAGAI

#ifdef GALPOT_ALL

      // dPdR and dPdz are positive
      
      R = sqrt(dx*dx+dy*dy);
      fall(R,dz,&dPdR,&dPdz);

      force[3*i]   += -dPdR*dx/R;
      force[3*i+1] += -dPdR*dy/R;
      force[3*i+2] += -dPdz;
      
#endif //GALPOT_ALL

#ifdef STATICHQDE
      r = sqrt(dx*dx+dy*dy+dz*dz);
      force[3*i] -= G*M_HQ/pow(a_HQ+r,2.)*dx/r;
      force[3*i+1] -= G*M_HQ/pow(a_HQ+r,2.)*dy/r;
      force[3*i+2] -= G*M_HQ/pow(a_HQ+r,2.)*dz/r;

#endif //STATICHQDE

#ifdef STATICNFW
      nhatdotr = dx*cos(phi_NFW)*cos(theta_NFW)+dy*cos(phi_NFW)*sin(theta_NFW)+dz*sin(phi_NFW);
      r = sqrt(dx*dx+dy*dy+dz*dz+(pow(q_NFW,-2.)-1.)*pow(nhatdotr,2.));
      force[3*i]   += -G*M_NFW*1./pow(r,3.)*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW))*(log(1.+r/rs_NFW)-r/(r+rs_NFW))*(dx+(pow(q_NFW,-2.)-1.)*nhatdotr*cos(phi_NFW)*cos(theta_NFW));
      force[3*i+1] += -G*M_NFW*1./pow(r,3.)*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW))*(log(1.+r/rs_NFW)-r/(r+rs_NFW))*(dy+(pow(q_NFW,-2.)-1.)*nhatdotr*cos(phi_NFW)*sin(theta_NFW));
      force[3*i+2] += -G*M_NFW*1./pow(r,3.)*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW))*(log(1.+r/rs_NFW)-r/(r+rs_NFW))*(dz+(pow(q_NFW,-2.)-1.)*nhatdotr*sin(phi_NFW));
#endif

#ifdef STATICLOG
      r = sqrt(dx*dx+dy*dy+dz*dz);
      force[3*i] += -pow(vc_LOG,2.)/(r*r+rs_LOG*rs_LOG)*2.*dx;
      force[3*i+1] += -pow(vc_LOG,2.)/(r*r+rs_LOG*rs_LOG)*2.*dy;
      force[3*i+2] += -pow(vc_LOG,2.)/(r*r+rs_LOG*rs_LOG)*2.*dz;
#endif

#ifdef STATICLM2010
      denom = c1*dx*dx+c2*dy*dy+c3*dx*dy+pow(dz/qz_LM,2.)+pow(rhalo,2.);
      force[3*i] += -pow(vhalo,2.)*(2.*c1*dx+c3*dy)/denom;
      force[3*i+1] += -pow(vhalo,2.)*(2.*c2*dy+c3*dx)/denom;
      force[3*i+2] += -pow(vhalo,2.)*(2.*dz/pow(qz_LM,2.))/denom;

#endif //STATICLM2010


#ifdef STATICBPL
      r = sqrt(dx*dx+dy*dy+dz*dz);
      gamma_store = gsl_sf_gamma(3-BPL_alpha);
      gamma_inc_store = gsl_sf_gamma_inc(3-BPL_alpha,r/BPL_r);
      force[3*i] += -G*BPL_M*(gamma_store-gamma_inc_store)/gamma_store*dx/pow(r,3.);
      force[3*i+1] += -G*BPL_M*(gamma_store-gamma_inc_store)/gamma_store*dy/pow(r,3.);
      force[3*i+2] += -G*BPL_M*(gamma_store-gamma_inc_store)/gamma_store*dz/pow(r,3.);
      //force[3*i+2] += -G*BPL_M*(gsl_sf_gamma(3-BPL_alpha)-gsl_sf_gamma_inc(3-BPL_alpha,r/BPL_r))/gsl_sf_gamma(3-BPL_alpha)*dz/pow(r,3.);
#endif 

#ifdef STATICBOWDENROT
      dx = 0.;
      for(axes1=0;axes1<3;axes1++)
	dx += M_rot_inv[0][axes1]*(pdata[i].x[axes1]-pos_mw_t[axes1]);
      dy = 0.;
      for(axes1=0;axes1<3;axes1++)
	dy += M_rot_inv[1][axes1]*(pdata[i].x[axes1]-pos_mw_t[axes1]);
      dz = 0.;
      for(axes1=0;axes1<3;axes1++)
	dz += M_rot_inv[2][axes1]*(pdata[i].x[axes1]-pos_mw_t[axes1]);      
				   
      r = sqrt(dx*dx+dy*dy+dz*dz);
      mdpsiNFWdr = -4.*M_PI*G*rho0_NFW*pow(rs_NFW,3.)/pow(r,2.)*log(1.+r/rs_NFW)+4.*M_PI*G*rho0_NFW*pow(rs_NFW,2.)/r*1./(1.+r/rs_NFW);
      mdpsiy20dr = -4.*M_PI*G*rho1_y20*pow(r1_y20,3.)/pow(r+r1_y20,2.)*(1.5*pow(dz/r,2.)-0.5)+8.*M_PI*G*rho1_y20*pow(r1_y20,3.)*r/pow(r+r1_y20,3.)*(1.5*pow(dz/r,2.)-0.5)+12.*M_PI*G*rho1_y20*pow(r1_y20,3.)*r/pow(r+r1_y20,2.)*pow(dz,2.)/pow(r,3.);

      a_inertial[0] = dx/r*(mdpsiNFWdr+mdpsiy20dr);
      a_inertial[1] = dy/r*(mdpsiNFWdr+mdpsiy20dr);
      a_inertial[2] = dz/r*(mdpsiNFWdr+mdpsiy20dr) - 12.*M_PI*G*rho1_y20*pow(r1_y20,3.)*r/pow(r+r1_y20,2.)*dz/pow(r,2.);

      for(axes1=0;axes1<3;axes1++)
	a_rot[axes1] = 0.;
      
      for(axes1=0;axes1<3;axes1++)
	for(axes2=0;axes2<3;axes2++)
	  a_rot[axes1] += M_rot[axes1][axes2]*a_inertial[axes2];

      for(axes1=0;axes1<3;axes1++)
	for(axes2=0;axes2<3;axes2++)
	  force[3*i+axes1] += M_rot[axes1][axes2]*a_inertial[axes2];
#endif //STATICBOWDENROT

#ifdef INCLUDE_LMC
      dx = pdata[i].x[0]-pos_lmc_t[0];
      dy = pdata[i].x[1]-pos_lmc_t[1];
      dz = pdata[i].x[2]-pos_lmc_t[2];
      r = sqrt(dx*dx+dy*dy+dz*dz);
      
      force[3*i] -= G*M_LMC/pow(r+rs_LMC,2.)*dx/r;
      force[3*i+1] -= G*M_LMC/pow(r+rs_LMC,2.)*dy/r;
      force[3*i+2] -= G*M_LMC/pow(r+rs_LMC,2.)*dz/r;
#endif


    }

}

double find_nbody_force_w_dt(struct pdata_struct *pdata, double *force, double *pos_prog_t)
{
  int i,k;
  double dt=dt_max;
  double r;
  double accel=0.;

  for(i=0;i<N;i++)
    {
      r = dist(pdata[i].x,pos_prog_t);
      for(k=0;k<3;k++)
        {
          force[3*i+k]+=G*M_plum*(pos_prog_t[k]-pdata[i].x[k])/pow(r*r+sigma_plum*sigma_plum,1.5);
          accel += pow(G*M_plum*(pos_prog_t[k]-pdata[i].x[k])/pow(r*r+sigma_plum*sigma_plum,1.5),2.);
        }
      accel = sqrt(accel);
      dt = min2(dt,eta*sqrt(max2(r,sigma_plum)/accel));
    }
  return dt;
}


void find_nbody_force(struct pdata_struct *pdata, double *force, double *pos_prog_t)
{
  int i,k;
  double r;
  for(i=0;i<N;i++)
    {
      r = dist(pdata[i].x,pos_prog_t);
      for(k=0;k<3;k++)
        force[3*i+k]+=G*M_plum*(pos_prog_t[k]-pdata[i].x[k])/pow(r*r+sigma_plum*sigma_plum,1.5);
    }
}



void generate_stream(struct pdata_struct *pdata, double *force, double t, double *pos_prog_t, double *vel_prog_t, double *pos_mw_t, double *vel_mw_t, int leading_or_trailing)
{
  int i,j;
  int i_prog=0; // remember that the force is actually the one on the progenitor
  //double M_enc = sqrt(pow(force[3*i_prog],2.)+pow(force[3*i_prog+1],2.)+pow(force[3*i_prog+2],2.));
  double r_gal = dist(pos_prog_t,pos_mw_t);

  //M_enc = M_enc*r_gal*r_gal/G;
  double r_lag = r_lag_calc(pos_prog_t,vel_prog_t,pos_mw_t,vel_mw_t); //r_gal*pow(M_plum/(2.*M_enc),1./3.);
  double v1,v2,s;
  
  double omega,r;
  r = dist(pos_prog_t,pos_mw_t);

  double x_rel[3],v_rel[3];

  for(i=0;i<3;i++)
    {
      x_rel[i] = pos_prog_t[i] - pos_mw_t[i];
      v_rel[i] = vel_prog_t[i] - vel_mw_t[i];
    }
  
  omega = sqrt(  pow(x_rel[0]*v_rel[1]-x_rel[1]*v_rel[0],2.) + pow(x_rel[0]*v_rel[2]-x_rel[2]*v_rel[0],2.) + pow(x_rel[1]*v_rel[2]-x_rel[2]*v_rel[1],2.) )/( r*r );
  
  //double sigma_lag = sqrt(2.5)*0.9*pow(M_plum/(1.e-5),0.25)*pow(omega/44.,1./3.); 
  double sigma_lag = lambda_vlag*sqrt(1./3.*G*M_plum/sqrt(r_lag*r_lag+sigma_plum*sigma_plum));
// sqrt(1./3.*G*M_plum/sqrt(r_lag*r_lag+sigma_plum*sigma_plum));//sigma_vel;//0.25*sqrt(G*M_plum/(3.*r_lag));

  double sign; 



  
  if(leading_or_trailing==1)
    sign = 1.;
  else
    sign = -1.;

  if(N+1 <= NMAX)
    {
      for(i=0;i<3;i++)
	{
	  pdata[0].x[i] = x_rel[i]-sign*x_rel[i]/r_gal*r_lag+pos_mw_t[i];
	  pdata[0].t_strip = t;
	  
	  pdata[0].phase = -sign*r_lag;
	  pdata[0].theta = 0.;
	  pdata[0].theta_phase = 0.;
	  
	  do 
	    {
	      v1 = 2.0 * ((float) rand()/RAND_MAX) - 1.0;
	      v2 = 2.0 * ((float) rand()/RAND_MAX) - 1.0;
	      
	      s = v1*v1 + v2*v2;
	    }while ( s >= 1.0 || s == 0);
	  
	  pdata[0].v[i] = v_rel[i]*(1.-lambda_vfrac*sign*r_lag/r_gal)+sigma_lag*v1*sqrt(-2.0*log(s)/s) + vel_mw_t[i];
	  //pdata[0].v[i] = v_rel[i]*(1.-lambda_vfrac*sign*r_lag/r_gal)+sigma_lag*rand_num[rand_counter]+vel_mw_t[i];
	}
    }
}

void update_prog_pos(double *pos_prog_t, double t, double *pos_prog)
{
  int counter = (int)(t/(0.0001*1.0227));
  int axes;
  if(counter < num_pos_prog-1)
    {
      for(axes=0;axes<3;axes++)
	pos_prog_t[axes] = pos_prog[3*(num_pos_prog-counter-1)+axes]+(pos_prog[3*(num_pos_prog-counter-2)+axes]-pos_prog[3*(num_pos_prog-counter-1)+axes])/(0.0001*1.0227)*(t-(double)counter*0.0001*1.0227);
    }
  else
    {
      for(axes=0;axes<3;axes++)
	pos_prog_t[axes] = pos_prog[3*(num_pos_prog-counter-1)+axes];
    }
}

void update_lmc_pos(double *pos_lmc_t, double t, double *pos_lmc)
{
  int counter = (int)(t/(0.0001*1.0227));
  int axes=0;
  if(counter < num_pos_prog - 1)
    {
      for(axes=0;axes<3;axes++)
	pos_lmc_t[axes] = pos_lmc[3*(num_pos_prog-counter-1)+axes]+(pos_lmc[3*(num_pos_prog-counter-2)+axes]-pos_lmc[3*(num_pos_prog-counter-1)+axes])/(0.0001*1.0227)*(t-(double)counter*0.0001*1.0227);
    }
  else
    {
      for(axes=0;axes<3;axes++)
	pos_lmc_t[axes] = pos_lmc[3*(num_pos_prog-counter-1)+axes];
    }
	  
}


void update_prog_vel(double *vel_prog_t, double t, double *vel_prog)
{
  int counter = (int)(t/(0.0001*1.0227));
  int axes;
  if(counter < num_pos_prog - 1)
    {
      for(axes=0;axes<3;axes++)
	vel_prog_t[axes] = vel_prog[3*(num_pos_prog-counter-1)+axes]+(vel_prog[3*(num_pos_prog-counter-2)+axes]-vel_prog[3*(num_pos_prog-counter-1)+axes])/(0.0001*1.0227)*(t-(double)counter*0.0001*1.0227);
    }
  else
    {
      for(axes=0;axes<3;axes++)
	vel_prog_t[axes] = vel_prog[3*(num_pos_prog-counter-1)+axes];
    }
}

void update_mw_pos(double *pos_mw_t, double t, double *pos_mw)
{
  int counter = (int)(t/(0.0001*1.0227));
  int axes=0;
  if(counter < num_pos_prog - 1)
    {
      for(axes=0;axes<3;axes++)
	pos_mw_t[axes] = pos_mw[3*(num_pos_prog-counter-1)+axes]+(pos_mw[3*(num_pos_prog-counter-2)+axes]-pos_mw[3*(num_pos_prog-counter-1)+axes])/(0.0001*1.0227)*(t-(double)counter*0.0001*1.0227);
    }
  else
    {
      for(axes=0;axes<3;axes++)
	pos_mw_t[axes] = pos_mw[3*(num_pos_prog-counter-1)+axes];
    }
	  
}

void update_mw_vel(double *vel_mw_t, double t, double *vel_mw)
{
  int counter = (int)(t/(0.0001*1.0227));
  int axes;
  if(counter < num_pos_prog - 1)
    {
      for(axes=0;axes<3;axes++)
	vel_mw_t[axes] = vel_mw[3*(num_pos_prog-counter-1)+axes]+(vel_mw[3*(num_pos_prog-counter-2)+axes]-vel_mw[3*(num_pos_prog-counter-1)+axes])/(0.0001*1.0227)*(t-(double)counter*0.0001*1.0227);
    }
  else
    {
      for(axes=0;axes<3;axes++)
	vel_mw_t[axes] = vel_mw[3*(num_pos_prog-counter-1)+axes];
    }
}

void update_r_prog_lmc(double *r_prog_lmc_t, double t, double *r_prog_lmc)
{
  int counter = (int)(t/(0.0001*1.0227));
  int axes;
  if(counter < num_pos_prog - 1)
    {
      r_prog_lmc_t[0] = r_prog_lmc[(num_pos_prog-counter-1)]+(r_prog_lmc[(num_pos_prog-counter-2)]-r_prog_lmc[(num_pos_prog-counter-1)])/(0.0001*1.0227)*(t-(double)counter*0.0001*1.0227);
    }
  else
    {
      r_prog_lmc_t[0] = r_prog_lmc[(num_pos_prog-counter-1)];
    }
}

double gen_normal(double mean, double sigma)
{

  double v1,v2,s;

  do
    {
      v1 = 2.0 * ((float) rand()/RAND_MAX) - 1.0;
      v2 = 2.0 * ((float) rand()/RAND_MAX) - 1.0;

      s = v1*v1 + v2*v2;
    }while ( s >= 1.0 || s == 0);

  return mean + sigma*v1*sqrt(-2.0*log(s)/s);

}



double d2phidr2_HQ(double r)
{
  return -2.*G*M_HQ/pow(r+a_HQ,3.);
}

double d2phidR2_MN(double R, double z)
{
  
  double denom;
  double total=0;
  int i_MN;

  for(i_MN=0;i_MN<n_MN;i_MN++)
    {
      denom = pow(R,2.)+pow(a_MN[i_MN]+sqrt(z*z+b_MN[i_MN]*b_MN[i_MN]),2.);
      total += G*M_MN[i_MN]*(pow(a_MN[i_MN]+sqrt(z*z+b_MN[i_MN]*b_MN[i_MN]),2.)-2.*R*R)/pow(denom,2.5);
    }
  return total;
}

double d2phidRdz_MN(double R, double z)
{
  double denom;
  double total=0.;
  int i_MN;
  
  for(i_MN=0;i_MN<n_MN;i_MN++)
    {
      denom = pow(R,2.)+pow(a_MN[i_MN]+sqrt(z*z+b_MN[i_MN]*b_MN[i_MN]),2.);
      total += -3.*G*M_MN[i_MN]*(a_MN[i_MN]+sqrt(z*z+b_MN[i_MN]*b_MN[i_MN]))*z*R/(pow(denom,2.5)*sqrt(z*z+b_MN[i_MN]*b_MN[i_MN]));
    }
  
  return total;

}

double d2phidz2_MN(double R, double z)
{
  double denom1;
  double denom2;
  double numer1;
  double total=0.;
  int i_MN;

  for(i_MN=0;i_MN<n_MN;i_MN++)
    {
      denom1 =  pow(R,2.)+pow(a_MN[i_MN]+sqrt(z*z+b_MN[i_MN]*b_MN[i_MN]),2.);
      denom2 = z*z+b_MN[i_MN]*b_MN[i_MN];
      numer1 = a_MN[i_MN]+sqrt(z*z+b_MN[i_MN]*b_MN[i_MN]);

      total += G*M_MN[i_MN]*numer1/(pow(denom1,1.5)*pow(denom2,0.5)) + G*M_MN[i_MN]*z*z/(pow(denom1,1.5)*denom2)  - G*M_MN[i_MN]*z*z*numer1/(pow(denom1,1.5)*pow(denom2,1.5)) - 3.*G*M_MN[i_MN]*z*z*pow(numer1,2.)/(pow(denom1,2.5)*denom2);
    }

  return total;

}

double dphidr_NFW(double r)
{
  return G*M_NFW/pow(r,2.)*log(1.+r/rs_NFW)/(log(1.+c_NFW)-c_NFW/(1.+c_NFW)) - G*M_NFW/r*1./(r+rs_NFW)*1./(log(1.+c_NFW)-c_NFW/(1.+c_NFW));
}

double d2phidr2_NFW(double r)
{
  return G*M_NFW/pow(r,3.)*(-2.*log(1.+r/rs_NFW)+r*(3.*r+2.*rs_NFW)/pow(r+rs_NFW,2.) )/(log(1.+c_NFW)-c_NFW/(1.+c_NFW));
}

double d2phidR2_NFW(double R, double z)
{
  double r = sqrt(R*R+pow(z/q_NFW,2.));
  return d2phidr2_NFW(r)*R*R/(r*r) + dphidr_NFW(r)*z*z/(q_NFW*q_NFW*pow(r,3.));
}

double d2phidRdz_NFW(double R, double z)
{
  double r = sqrt(R*R+pow(z/q_NFW,2.));
  return d2phidr2_NFW(r)*R/r*z/(q_NFW*q_NFW*r) - dphidr_NFW(r)*R*z/(q_NFW*q_NFW*pow(r,3.));
}

double d2phidz2_NFW(double R, double z)
{
  double r = sqrt(R*R+pow(z/q_NFW,2.));
  return d2phidr2_NFW(r)*z*z/(pow(q_NFW,4.)*r*r) + dphidr_NFW(r)*R*R/(q_NFW*q_NFW*pow(r,3.));
}

double d2phidr2(double R, double z)
{
  double d2phidR2, d2phidRdz, d2phidz2,r;
  r = sqrt(R*R+z*z);
  d2phidR2 = d2phidR2_MN(R,z) + d2phidR2_NFW(R,z);
  d2phidRdz = d2phidRdz_MN(R,z) + d2phidRdz_NFW(R,z);
  d2phidz2 = d2phidz2_MN(R,z) + d2phidz2_NFW(R,z);
  return d2phidr2_HQ(r) + d2phidR2*pow(R/r,2.) + 2.*d2phidRdz*R/r*z/r + d2phidz2*pow(z/r,2.);
}

double r_lag_calc(double *pos_prog_t, double *vel_prog_t, double *pos_mw_t, double *vel_mw_t)
{
  int i;

  double R,z,r_NFW,r;

  R = sqrt( pow(pos_prog_t[0]-pos_mw_t[0],2.)+pow(pos_prog_t[1]-pos_mw_t[1],2.) );
  z = pos_prog_t[2]-pos_mw_t[2];
  
  double omega;

  double x_rel[3],v_rel[3];

  for(i=0;i<3;i++)
    {
      x_rel[i] = pos_prog_t[i] - pos_mw_t[i];
      v_rel[i] = vel_prog_t[i] - vel_mw_t[i];
    }

  r = dist(pos_prog_t,pos_mw_t);
  
  omega = sqrt(  pow(x_rel[0]*v_rel[1]-x_rel[1]*v_rel[0],2.) + pow(x_rel[0]*v_rel[2]-x_rel[2]*v_rel[0],2.) + pow(x_rel[0]*v_rel[2]-x_rel[2]*v_rel[0],2.) )/( r*r );
  
  return pow( G*M_plum/(omega*omega - d2phidr2(R,z)), 1./3.);
}

double pot_NFW(double r)
{
  return -G*M_NFW/r*log(1.+r/rs_NFW)/(log(1.+c_NFW)-c_NFW/(1.+c_NFW));
}

double pot_HQ(double r)
{
  return -G*M_HQ/(r+a_HQ);
}

double pot_MN(double R, double z)
{
  double total_phi,denom;
  int i_MN;
  
  total_phi = 0.;
  
  for(i_MN=0;i_MN<n_MN;i_MN++)
    {
      denom = pow(R,2.)+pow(a_MN[i_MN]+sqrt(z*z+b_MN[i_MN]*b_MN[i_MN]),2.);
      total_phi -= G*M_MN[i_MN]/sqrt(denom);
    }
  return total_phi;
}


double pot_total(struct pdata_struct *pdata, int counter)
{
  double R, z, r;

  R = sqrt( pow(pdata[counter].x[0],2.)+pow(pdata[counter].x[1],2.) );
  z = pdata[counter].x[2];

  r = mag(pdata[counter].x);

  return pot_NFW(r)+pot_HQ(r)+pot_MN(R,z);
}


double E_rel_calc(struct pdata_struct *pdata, int counter, double *pos_prog_t, double *vel_prog_t)
{
  double r_rel;

  r_rel = dist(pdata[counter].x,pos_prog_t);

  
  
  return 0.5*pow(dist(pdata[counter].v,vel_prog_t),2.) - G*M_plum/sqrt(r_rel*r_rel+sigma_plum*sigma_plum);
  

}

