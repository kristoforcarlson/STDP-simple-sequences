//Kris Carlson
//Purdue University

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
//#include "DynArray.h"

//label each one of our channel currents
#define HH_Leak 0
#define HH_K 1
#define HH_Na 2
#define UKB_CA_LTYPE 6

#define NUMSEG 20

//create our compartment structures
//structure for a current from ion channel
//There are currently seven types of currents
struct iion
{
  double gMax,gMaxTot,Erev;
  double n,m,h,i;
  //only used in certain ion channels
  double Na_att,tau;
  double alpha,beta,k,phi;
  double val;//this is where we can be unique to channel
};

//structure for a equipotential compartment
struct cmpt
{
  //variables common to all compartments
  double v,Erest;
  double ra,Ra,RaTot;
  double rm,Rm,RmTot;
  double cm,Cm,CmTot;
  double tau, lambda;
  double length, diameter,distance;
  double Iinj,injectionTime,injectionDuration;//amperes/meters,seconds,seconds
  struct iion Ichan[7];
};

typedef struct iion Iion;
typedef struct cmpt Cmpt;

//HH support functions
double hh_alpha_n(double n),hh_beta_n(double n);
double hh_alpha_m(double n),hh_beta_m(double n);
double hh_alpha_h(double n),hh_beta_h(double n);
//init function
void init(Cmpt [], int n2);
void update_HH_K_n(Cmpt [], int size, double dt);
void update_HH_Na_m(Cmpt [], int size, double dt);
void update_HH_Na_h(Cmpt [], int size, double dt);
void update_UKB_Ca_Ltype_m(Cmpt [], int size, double dt);
void update_UKB_Ca_Ltype_h(Cmpt [], int size, double dt);
void update_leak(Cmpt [], int size);
void update_HH_Na(Cmpt [], int size);
void update_HH_K(Cmpt [], int size);
void update_UKB_Ca_Ltype(Cmpt [], int size);
void update_v_Soma(Cmpt [], int n2, double n3, double n4, double n5);//works
void update_v_Passive(Cmpt Seg[], int size, double dx, double dt, double t);//works
void update_v_HH_Soma(Cmpt Seg[], int size, double dx, double dt, double t);//works
void update_v_HH_Cell(Cmpt [], int n2, double n3, double n4, double n5);

void structTest(Cmpt []);

//new stuff + Urakubo Calcium Current
double i_NMDA_Ca(double v, double t, double gRelease, double gMax);
double i_NMDA_Ca_Two(double v, double t, double gRelease1, double gRelease2, double gMax);
double i_NMDA_Ca_Two_Diff(double v, double t, double gRelease1, double gRelease2, double gMax, double tauCaDiff);
double i_NMDA_Ca_Two_Un_Sat_No_Fast(double v, double t, double gRelease1, double gRelease2, double gMax);
double i_NMDA_Ca_Two_Un_Sat_F_and_S(double v, double t, double gRelease1, double gRelease2, double gMax);
double glutBind(double t, double gRelease);
double glutBindTwo(double t, double gRelease1, double gRelease2);
double glutBindTwoUnSatNoFast(double t, double gRelease1, double gRelease2);
double glutBindTwoUnSatFandS(double t, double gRelease1, double gRelease2);
double glutBindTwoDiff(double t, double gRelease1, double gRelease2, double tauCaDiff);
double i_VDCC_simple(double v, double t, double injectionTime);
double i_VDCC_short(double v, double t, double injectionTime);
double i_VDCC_short_update(double v, double t, double injectionTime, double gMax);
double i_VDCC_short_two_update(double v, double t, double injectionTime1, double injectionTime2, double gMax);
double VDCC_Ca_Diff(double t, double injectionTime, double Ca_VDCC, double Ca_value, double tauCaDiff);
double BPAP(double t, double injectionTime);
double i_UKB_Ca_Ltype(Cmpt [], double t, double dt);

//global vars
double Pi=3.14159;

int main(int argc, char** argv)
{
  double dt,dx;
  double t,totalTime;//seconds 
  double vCoop; //Cooper's voltage for comparison

  //file writing variables
  FILE *fp_v1,*fp_v2,*fp_v3;
  char v1_file[100],v2_file[100],v3_file[100];

  FILE *fp_vCoop;
  char vCoop_file[100];

  FILE *fp_test1,*fp_test2,*fp_test3;
  char test1_file[100],test2_file[100],test3_file[100];

  //Urakubo stuff
  FILE *fp_vend;
  char vend_file[100];

  //Lisman file pointers
  FILE *fp_pK, *fp_P, *fp_Ca;
  char pK_file[100], P_file[100], Ca_file[100];

  //file writing for new different calcium pools
  FILE *fp_Ca_VDCC, *fp_Ca_NMDA;
  char Ca_VDCC_file[100], Ca_NMDA_file[100];

  //file writing for comparing pk and P irrespective of time
  FILE *fp_pK_P;
  char pK_P_file[100];
  
  //Create our compartments
  Cmpt Seg[NUMSEG];

  //Lisman model important variables
  //constants
  double CaConc;
  double Ktot,Ptot,Atot;
  double K0,P0;
  double k1,k2,k3,k4,k11,k12,k13,k14;
  double Km,Km1,Km2,Km11,Km12;
  double c1,c2,c3,c4;
  double Ca;
  //these are the local concentrations 
  double Ca_VDCC, Ca_NMDA;
  double Ca_VDCC_final, Ca_NMDA_final;
  double CaConcAvg, CaAvg;
  double A,Aint;
  double k21,k22;
  double Ca_basal;
  //vars
  double K,pK,P,pP;  
  
  //time of glutamate release
  double gRelease;
  double gRelease1,gRelease2;
  double InjectionOffset;
  double tauCa;
  double tauCaDiff;
  double AvoNum;
  double eleCharge;
  double faradaysConstant;//C/mol
  //features of an individual spine
  //estimate spine as cube with sides 0.46 microns
  double spineArea;
  double activeSpineArea;
  double spineVolume;
  //so we can iterate over gNMDA
  double gNMDA;
  //so we can control gVDCC sooner
  double gVDCC;

  //pre and post repetitive stimulation variables
  double T;
  double f;
  double tp_T_iter1;
  double tp_T_back1;
  double tp_T_iter2;
  double tp_T_back2;  
  double ta_T_iter1;
  double ta_T_back1;
  double ta_T_iter2;
  double ta_T_back2;

  //nano domain variables
  double volFrac;
  
  dt=1e-6;//seconds
  //dt=0.1e-6;//<--original pwalk value
  dx=20e-6;//meters <-Urakubo
  t=0;//seconds
  //totalTime=0.500;//seconds
  //totalTime=0.1;//seconds
  //totalTime=2;
  //eventually our totalTime must be 20 seconds
  //totalTime=0.5;
  totalTime=12;
  //Some Lisman variables
  tauCa=12e-3;
  //New diffusion variable for nano domain calculations
  tauCaDiff=5e-3;
  AvoNum=6.022e23;
  eleCharge=1.602e-19;
  faradaysConstant=96485.339; //C/mol
  //spine attributes
  spineArea=(0.46e-6)*(0.46e-6);
  activeSpineArea=spineArea/6;//one surface
  spineVolume=(0.46e-6)*(0.46e-6)*(0.46e-6);
  
  //original 10ms
  gRelease1=1.010;
  gRelease2=1.020+totalTime;
  
  //for loop to iterate over time intervals
  //for(gRelease=Seg[0].injectionTime-0.060;gRelease<Seg[0].injectionTime+0.060;gRelease=gRelease+0.010)
  //{
  //printf("gRelease=%g\n",gRelease);
  printf("gRelease1=%g\n",gRelease1);
  printf("gRelease2=%g\n",gRelease2);
 
  //Initialize Cmpt and Iion structures
  init(Seg,NUMSEG);
  
  //set the time and duration of injection
  Seg[0].injectionTime=1.000;
  //Seg[0].injectionTime=0.050;
  Seg[0].injectionDuration=0.001;
  printf("post-synaptic injectionTime1=%g\n",Seg[0].injectionTime);
  //Seg[0].injectionDuration=1e-5;
  InjectionOffset=0.020;
  printf("post-synaptic injectionTime2=%g\n",Seg[0].injectionTime+InjectionOffset);
  //InjectionOffset=0.010;
  //quick temporary test
  //InjectionOffset=totalTime+0.420;//pairs
  
  //constants
  Ktot=20;Ptot=20;Atot=1;
  K0=0.5;P0=0.5;
  //careful, check to make sure what values k4 and k14 are
  k1=2;k2=15;k3=1;k4=120;
  k11=2;k12=15;k13=1;k14=80;
  Km=4;Km1=10;Km2=0.3;Km11=10;Km12=1;
  c1=1;c2=1;c3=6;c4=8;
  Ca=0.1;
  Ca_basal=0.1e-6;//must be in SI because it goes through Cooper/Urakubo model first which
  //is in SI
  
  //initialize Lisman vars
  pK=0.3;
  K=Ktot-pK;
  P=0.2;
  pP=Ptot-P;
  A=0.6;
  Aint=Atot-A;
  k21=c1*pK+c3;
  k22=c2*P+c4;
  Ca=0.1;
  CaConc=0.1e-6;
  Ca_NMDA=0;
  Ca_VDCC=0;
  //test on UKB_CA_LTYPE VDCC
  Seg[3].Ichan[UKB_CA_LTYPE].m=0.5;
  Seg[3].Ichan[UKB_CA_LTYPE].h=0.5;
  //nano domain assignments
  volFrac=0.1;
      
  int counter;
  counter=0;
  //Do the calculation
  //for loop for gVDCC conductance (Seg[i].Ichan[UKB_CA_LTYPE].gMax)
  //for(Seg[3].Ichan[UKB_CA_LTYPE].gMax=3;Seg[3].Ichan[UKB_CA_LTYPE].gMax<16;Seg[3].Ichan[UKB_CA_LTYPE].gMax=Seg[3].Ichan[UKB_CA_LTYPE].gMax+3)
  //{
  //for loop for gNMDA
  //for(gNMDA=30;gNMDA<151;gNMDA=gNMDA+30)
  //{
  //gNMDA=175;//original Cooper goes with Pb=0.5
  //gNMDA=175*1.7;//original goes with Pb=0.3
  gNMDA=17.5*1.7;//original goes with Pb=0.3
  //17.5*1.7=29.75
  
  //12/01/2010
  //gVDCC=16;//<- does not actually work with this

  //12/01/2010
  gVDCC=16.5;//<- works with this could probably go down to 16.2?
  //Seg[3].Ichan[UKB_CA_LTYPE].gMax=8; //only used for Urakubo channel versions of VDCCs
  
  //name the files you are opening
  sprintf(v3_file,"v3-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  sprintf(test1_file,"glut-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  sprintf(test2_file,"I_VDCC-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  sprintf(test3_file,"I_NMDA-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  sprintf(pK_file,"pK-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  sprintf(P_file,"P-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  //sprintf(Ca_file,"Ca-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  sprintf(Ca_file,"Ca-gNMDA-%.3f-gVDCC-%.3f",gNMDA,gVDCC);
  //Lee
  //sprintf(Ca_VDCC_file,"Ca_VDCC-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  //sprintf(Ca_NMDA_file,"Ca_NMDA-gNMDA-%.3f-gVDCC-%.3f",gNMDA,Seg[3].Ichan[UKB_CA_LTYPE].gMax);
  sprintf(Ca_VDCC_file,"Ca_VDCC-gNMDA-%.3f-gVDCC-%.3f",gNMDA,gVDCC);
  sprintf(Ca_NMDA_file,"Ca_NMDA-gNMDA-%.3f-gVDCC-%.3f",gNMDA,gVDCC);

  //The pK/P dynamical system
  sprintf(pK_P_file,"pK_P-gNMDA-%.3f-gVDCC-%.3f",gNMDA,gVDCC);

  //open the files for writing before you start
  fp_v3 = fopen(v3_file,"w");
  fp_test1 = fopen(test1_file,"w");
  fp_test2 = fopen(test2_file,"w");
  fp_test3 = fopen(test3_file,"w");
  fp_Ca = fopen(Ca_file,"w");
  fp_pK = fopen(pK_file,"w");
  fp_P = fopen(P_file,"w");
  fp_Ca_NMDA = fopen(Ca_NMDA_file,"w");
  fp_Ca_VDCC = fopen(Ca_VDCC_file,"w"); 
  //The pK/P dynamical system
  fp_pK_P = fopen(pK_P_file,"w");

  //reinitialize everything
  pK=0.3;
  K=Ktot-pK;
  P=0.2;
  pP=Ptot-P;
  A=0.6;
  Aint=Atot-A;
  k21=c1*pK+c3;
  k22=c2*P+c4;
  Ca=0.1;
  CaConc=0.1e-6;
  //test on UKB_CA_LTYPE VDCC
  Seg[3].Ichan[UKB_CA_LTYPE].m=0.5;
  Seg[3].Ichan[UKB_CA_LTYPE].h=0.5;
  
  f=1;
  //T=1/f;
  T=0;//works fine to get single spike pairs
  tp_T_iter1=0;
  tp_T_back1=0;
  tp_T_iter2=0;
  tp_T_back2=0;
  ta_T_iter1=0;
  ta_T_back1=0;
  ta_T_iter2=0;
  ta_T_back2=0;
  
  for(t=0;t<totalTime;t=t+dt)
    {
      
      if(t>(gRelease1+tp_T_iter1))
	{
	  tp_T_back1=tp_T_iter1;
	  tp_T_iter1=tp_T_iter1+T;
	  tp_T_back2=tp_T_iter2;
	  tp_T_iter2=tp_T_iter2+T;
	}
      
      if(t>(Seg[0].injectionTime+ta_T_iter1))
	{
	  ta_T_back1=ta_T_iter1;
	  ta_T_iter1=ta_T_iter1+T;
	  ta_T_back2=ta_T_iter2;
	  ta_T_iter2=ta_T_iter2+T;
	}

      //tell us when to stop the repititive stimulations
      if(t>61)
	{
	  ta_T_back1=0;
	  ta_T_back2=0;
	  tp_T_back1=0;
	  tp_T_back2=0;
	  ta_T_iter1=0;
	  ta_T_iter2=0;
	  tp_T_iter1=0;
	  tp_T_iter2=0;
	}
      
      if(t>=(Seg[0].injectionTime+ta_T_back1) && t<=(Seg[0].injectionTime+ta_T_back1+Seg[0].injectionDuration))
	{
	  //Seg[0].Iinj=20e-9;//had to turn it up
	  Seg[0].Iinj=105e-9;
	  //Seg[0].Iinj=75e-5;
	}//for two postsynaptic spikes
      else if(t>=(Seg[0].injectionTime+ta_T_back2+InjectionOffset) && t<=(Seg[0].injectionTime+ta_T_back2+InjectionOffset+Seg[0].injectionDuration))
	{
	  Seg[0].Iinj=105e-9;
	}
      else
	{
	  Seg[0].Iinj=0;
	}
      //update all channel variables before currents
      //and voltage updates
      update_HH_K_n(Seg,NUMSEG,dt);//update HH-nK
      update_HH_Na_m(Seg,NUMSEG,dt);//update HH-mNa
      update_HH_Na_h(Seg,NUMSEG,dt);//update HH-hNa
      //Urakubo functions
      update_UKB_Ca_Ltype_m(Seg,NUMSEG,dt);//update UKB Ca Ltype m
      update_UKB_Ca_Ltype_h(Seg,NUMSEG,dt);//update UKB Ca Ltype h
      //update current components before voltage compenents
      update_leak(Seg,NUMSEG);
      update_HH_Na(Seg,NUMSEG);
      update_HH_K(Seg,NUMSEG);
      //Urakubo functions
      //update_UKB_Ca_Ltype(Seg,NUMSEG);//<--don't implement this yet.
      //update_v_Soma(Seg,NUMSEG,dx,dt,t);//
      //update_v_Passive(Seg,NUMSEG,dx,dt,t);//
      //update_v_HH_Soma(Seg,NUMSEG,dt,dt,t);//
      update_v_HH_Cell(Seg,NUMSEG,dx,dt,t);
      
      //Cooper's BPAP
      //vCoop=BPAP(t,Seg[0].injectionTime);
      //v=-0.070;
      //hold voltage at 0 mV
      //v=0.0;
      //i_NMDA and i_VDCC are A/m^2 so multiple by active surface area (m^2)
      //with VDCC's
      //CaConc=CaConc+dt*( (i_NMDA_Ca(Seg[3].v,t,gRelease,gNMDA) + i_VDCC_simple(Seg[3].v,t,Seg[0].injectionTime))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((CaConc-Ca_basal)/tauCa));
      //with VDCC's and multiple presynaptic excitations
      //CaConc=CaConc+dt*( (i_NMDA_Ca_Two_Un_Sat_F_and_S(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA) + i_VDCC_short_update(Seg[3].v,t,Seg[0].injectionTime,gVDCC))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((CaConc-Ca_basal)/tauCa));

      //new quadruplet: the addition of the second bpap to pre-post-pre program
      CaConc=CaConc+dt*( (i_NMDA_Ca_Two_Un_Sat_F_and_S(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA) + i_VDCC_short_two_update(Seg[3].v,t,Seg[0].injectionTime+ta_T_back1,Seg[0].injectionTime+InjectionOffset+ta_T_back2,gVDCC))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((CaConc-Ca_basal)/tauCa));

      
      //Urakubo Calcium L-type channel's contribution to calcium current instead of simple version.
      //CaConc=CaConc+dt*( (i_NMDA_Ca(Seg[3].v,t,gRelease,gNMDA) + i_UKB_Ca_Ltype(Seg,t,dt))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((CaConc-Ca_basal)/tauCa));
      //with Urakubo L-type channel contribution and multiple presynaptic excitations
      //original
      //CaConc=CaConc+dt*( (i_NMDA_Ca_Two(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA) + i_UKB_Ca_Ltype(Seg,t,dt))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((CaConc-Ca_basal)/tauCa));
      //just NMDA without delay
      //CaConc=CaConc+dt*( (i_NMDA_Ca_Two(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((CaConc-Ca_basal)/tauCa));

      //Newest global calcium update function for nanodomain calculations
      //CaConc=CaConc+dt*( (i_NMDA_Ca_Two(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA)/(1-volFrac))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((CaConc-Ca_VDCC)/tauCaDiff) - ((CaConc-Ca_basal)/tauCa));
      
      //new Ca2+ avg calculation
      //CaConcAvg=volFrac*Ca_VDCC+(1-volFrac)*Ca_NMDA;

      //original without the delay
      Ca_NMDA=Ca_NMDA+dt*( (i_NMDA_Ca_Two_Un_Sat_F_and_S(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((Ca_NMDA-Ca_basal)/tauCa));

      //Ca_NMDA is calculated so VDCC Diffusion go can down to Ca_NMDA instead of basal
      //Ca_NMDA=Ca_NMDA+dt*( (i_NMDA_Ca_Two_Diff(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA,tauCaDiff)/(1-volFrac))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((Ca_NMDA-CaConcAvg)/tauCaDiff) - ((Ca_NMDA-Ca_basal)/tauCa));
      
      //Ca_VDCC=Ca_VDCC+dt*( (i_UKB_Ca_Ltype(Seg,t,dt))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((Ca_VDCC)/tauCa));

      //Ca_VDCC=Ca_VDCC+dt*((i_VDCC_short_update(Seg[3].v,t,Seg[0].injectionTime+ta_T_back,gVDCC))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((Ca_VDCC)/tauCa));

      Ca_VDCC=Ca_VDCC+dt*( (i_VDCC_short_two_update(Seg[3].v,t,Seg[0].injectionTime+ta_T_back1,Seg[0].injectionTime+InjectionOffset+ta_T_back2,gVDCC))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((Ca_VDCC)/tauCa));

      //list of the terms: 1st VDCC ion channel in small region, 2nd NMDA diffusion Ca2+ current
      //3rd buffering to basal, 4th diffusion to basal 
      //Ca_VDCC=Ca_VDCC+dt*( (i_VDCC_short(Seg[3].v,t,Seg[0].injectionTime+ta_T_back) + i_NMDA_Ca_Two_Diff(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA,tauCaDiff))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((Ca_VDCC-Ca_basal)/tauCa) - VDCC_Ca_Diff(t,Seg[0].injectionTime+ta_T_back,Ca_VDCC,Ca_basal,tauCaDiff));

      //Ca_VDCC=Ca_VDCC+dt*( (i_VDCC_short(Seg[3].v,t,Seg[0].injectionTime+ta_T_back) + i_NMDA_Ca_Two_Diff(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA,tauCaDiff))*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((Ca_VDCC-Ca_basal)/tauCa) - VDCC_Ca_Diff(t,Seg[0].injectionTime+ta_T_back,Ca_VDCC,Ca_NMDA,tauCaDiff));

      //Newest VDCC local calcium update function for nanodomain calculations
      //Ca_VDCC=Ca_VDCC+dt*( (i_VDCC_short_update(Seg[3].v,t,Seg[0].injectionTime+ta_T_back)/volFrac)*(1/(2*faradaysConstant))*(activeSpineArea/spineVolume)*(0.001) - ((Ca_VDCC-CaConcAvg)/tauCaDiff) - ((Ca_VDCC-Ca_basal)/tauCa));


      //change values to micromolar
      Ca=CaConc*1e6;
      Ca_VDCC_final=Ca_VDCC*1e6;
      Ca_NMDA_final=Ca_NMDA*1e6;
      //CaAvg=CaConcAvg*1e6;
      
      //Lisman model active kinase and phosphatase updates
      pK=pK+dt*(k1*pK*(K/(Km1+K))-k2*(P+P0)*(pK/(Km2+pK))+k3*K0+k4*K*((Ca*Ca*Ca*Ca)/((Km*Km*Km*Km)+(Ca*Ca*Ca*Ca))));
      P=P+dt*(k11*P*(pP/(Km11+pP))-k12*(pK+K0)*(P/(Km12+P))+k13*P0+k14*pP*((Ca*Ca*Ca)/((Km*Km*Km)+(Ca*Ca*Ca))));
      
      //conservation of kinases and phosphatases
      K=Ktot-pK;
      pP=Ptot-P;
	      
      //update function for AMPA
      //A=A+dt*(k21*Aint-k22*A);
      //A=A+dt*((c1*pK+c3)*Aint-(c2*P+c4)*A);
      //Aint=Atot-A;
      //k21=c1*pK+c3;
      //k22=c2*P+c4;
      
      //output data to files
      //Do it every 1000 iterations
      if(counter % 100 == 0)
	{
	  //fprintf(fp_v1,"%g\t%g\n",t,(Seg[0].v));
	  //fprintf(fp_v1,"%g\t%g\n",t,(Seg[0].v+0.058));
	  //fprintf(fp_v2,"%g\t%g\n",t,Seg[1].v);
	  
	  //12/1
	  fprintf(fp_v3,"%g\t%g\n",t,Seg[3].v);

	  //fprintf(fp_test1,"%g\t%g\n",t,i_VDCC_simple(Seg[3].v,t,Seg[0].injectionTime));
	  //fprintf(fp_test1,"%g\t%g\n",t,glutBindTwo(t,gRelease1+tp_T_back1,gRelease2+tp_T_back2));
	  
	  //12/1
	  fprintf(fp_test1,"%g\t%g\n",t,glutBindTwoUnSatFandS(t,gRelease1+tp_T_back1,gRelease2+tp_T_back2));

	  //fprintf(fp_test1,"%g\t%g\n",t,glutBind(t,gRelease1));
	  fprintf(fp_test2,"%g\t%g\n",t,i_VDCC_short_two_update(Seg[3].v,t,Seg[0].injectionTime+ta_T_back1,Seg[0].injectionTime+InjectionOffset+ta_T_back2,gVDCC));
	  //fprintf(fp_test2,"%g\t%g\n",t,Seg[3].Ichan[UKB_CA_LTYPE].h);
	  //fprintf(fp_test3,"%g\t%g\n",t,i_NMDA_Ca(Seg[3].v,t,gRelease,gNMDA));
	  //fprintf(fp_test3,"%g\t%g\n",t,i_NMDA_Ca_Two(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA));
	  
	  //12/1
	  fprintf(fp_test3,"%g\t%g\n",t,i_NMDA_Ca_Two_Un_Sat_F_and_S(Seg[3].v,t,gRelease1+tp_T_back1,gRelease2+tp_T_back2,gNMDA));
	  //fprintf(fp_test3,"%g\t%g\n",t,i_UKB_Ca_Ltype(Seg,t,dt));
	  //fprintf(fp_vend,"%g\t%g\n",t,(Seg[NUMSEG-1].v));
	  //fprintf(fp_vCoop,"%g\t%g\n",t,vCoop);
	  //fprintf(fp_vend,"%g\t%g\n",t,(Seg[NUMSEG-1].v+0.058));
	  //output Lisman variables
	 
	 //12/1
	 fprintf(fp_pK,"%g\t%g\n",t,pK);
	 fprintf(fp_P,"%g\t%g\n",t,P);
	 fprintf(fp_Ca,"%g\t%g\n",t,Ca);
	 fprintf(fp_Ca_NMDA,"%g\t\%g\n",t,Ca_NMDA_final);
	 fprintf(fp_Ca_VDCC,"%g\t\%g\n",t,Ca_VDCC_final);
	 
	 //12/11/2011@UCI
	 fprintf(fp_pK_P,"%g\t\%g\n",pK,P);
	 
	}
      counter++;
    }
  //close all files opened for writing
  //fclose(fp_v1);
  //fclose(fp_v2);
  fclose(fp_v3);
  fclose(fp_test1);
  fclose(fp_test2);
  fclose(fp_test3);
  //fclose(fp_test3);
  //fclose(fp_vend);
  //fclose(fp_vCoop);
  fclose(fp_pK);
  fclose(fp_P);
  fclose(fp_Ca);
  fclose(fp_Ca_NMDA);
  fclose(fp_Ca_VDCC);
  //}
  return 0;
}

double hh_alpha_n(double v)
{
        double val;
        v *= 1000;/*convert to mV for now*/
	if(fabs(v+55) > 1e-6) 
	  { 
	    val=(-55-v)/(100.0*(exp((-55-v)/10)-1.0));
	  }
	else 
	  {
	    val=0.10;
	  }
        val *= 1000;
        return(val);
}

double hh_beta_n(double v)
{
        double val;
        v *= 1000;/*convert to mV for now*/
	val=0.125*exp((-65-v)/80.0);
        val *= 1000;
        return(val);
}

double hh_alpha_m(double v)
{
        double val;
        v *= 1000;/*convert to mV for now*/
	val = (-40-v)/(10.0*(exp((-40-v)/10)-1.0));
        val *= 1000;
        return(val);
}

double hh_beta_m(double v)
{
        double val;
        v *= 1000;/*convert to mV for now*/    
	val=4.0*exp((-65-v)/18.0);
        val *= 1000;
        return(val);
}

double hh_alpha_h(double v)
{
        double val;
        v *= 1000;/* convert to mV for now*/
	val=0.07*exp((-65-v)/20); 
        val *= 1000;
        return(val);
}

double hh_beta_h(double v)
{
        double val;
        v *= 1000;/*convert ot mV for now*/
	val=1.0/(exp((-35-v)/10.0)+1.0);
        val *= 1000;//I dont' think this is necessary
        return(val);
}

void init(Cmpt Seg[], int size)
{
  //iterate array to assign correct parameters
  int i;
  //Create the soma compartment (special case)
  //Seg[0].Erest=-0.065;//mV
  Seg[0].Erest=-0.065;
  Seg[0].v=Seg[0].Erest;//mV
  Seg[0].length=20e-6;//meters <-Urakubo soma
  Seg[0].diameter=20e-6;//meters <-Urakubo soma
  Seg[0].distance=0;//distance from soma
  Seg[0].Ra=1.5;//ohms*meter <-Urakubo soma
  Seg[0].RaTot=Seg[0].length*(4*Seg[0].Ra)/(Pi*Seg[0].diameter*Seg[0].diameter);//Ohms
  Seg[0].ra=Seg[0].RaTot/Seg[0].length;
  Seg[0].Rm=2.8;//ohms*meter^2 <- Urakubo soma
  Seg[0].RmTot=Seg[0].Rm/(Pi*Seg[0].diameter*Seg[0].length);//ohms
  Seg[0].rm=Seg[0].RmTot*Seg[0].length;//ohms*meter
  
  Seg[0].Cm=0.01;//farads/meter^2 <- Urakubo soma

  Seg[0].CmTot=Seg[0].Cm*Pi*Seg[0].diameter*Seg[0].length;//farads
  Seg[0].cm=Seg[0].CmTot/Seg[0].length;//farads/meter
  Seg[0].tau=Seg[0].rm*Seg[0].cm;//seconds
  Seg[0].lambda=sqrt(Seg[0].rm/Seg[0].ra);//meters
  Seg[0].injectionTime=1;//seconds
  Seg[0].injectionDuration=0.001;//seconds
  Seg[0].Iinj=105e-9;//amperes/meter//equivalent to 2 pA/meter
  //channel initializations in soma
  //leak channel
  //Seg[0].Ichan[HH_Leak].gMax=0.7;//<--old error
  Seg[0].Ichan[HH_Leak].gMax=1/Seg[0].Rm;
  Seg[0].Ichan[HH_Leak].gMaxTot=Seg[0].Ichan[HH_Leak].gMax*Pi*Seg[0].diameter*Seg[0].length;
  Seg[0].Ichan[HH_Leak].Erev=-0.065;
  Seg[0].Ichan[HH_Leak].n=0;
  Seg[0].Ichan[HH_Leak].m=0;
  Seg[0].Ichan[HH_Leak].h=0;
  Seg[0].Ichan[HH_Leak].i=0;
  Seg[0].Ichan[HH_Leak].Na_att=0;
  Seg[0].Ichan[HH_Leak].alpha=0;
  Seg[0].Ichan[HH_Leak].beta=0;
  Seg[0].Ichan[HH_Leak].k=0;
  Seg[0].Ichan[HH_Leak].tau=0;
  Seg[0].Ichan[HH_Leak].phi=0;
  Seg[0].Ichan[HH_Leak].val=Seg[0].Ichan[HH_Leak].gMaxTot*(Seg[0].Ichan[HH_Leak].Erev-Seg[0].v);
  //HH potassium channel
  Seg[0].Ichan[HH_K].gMax=360;
  Seg[0].Ichan[HH_K].gMaxTot=Seg[0].Ichan[HH_K].gMax*Pi*Seg[0].diameter*Seg[0].length;
  Seg[0].Ichan[HH_K].Erev=-0.077;
  Seg[0].Ichan[HH_K].n=hh_alpha_n(Seg[0].Erest)/(hh_alpha_n(Seg[0].Erest + hh_beta_n(Seg[0].Erest)));
  Seg[0].Ichan[HH_K].m=0;
  Seg[0].Ichan[HH_K].h=0;
  Seg[0].Ichan[HH_K].i=0;
  Seg[0].Ichan[HH_K].Na_att=0;
  Seg[0].Ichan[HH_K].alpha=0;
  Seg[0].Ichan[HH_K].beta=0;
  Seg[0].Ichan[HH_K].k=0;
  Seg[0].Ichan[HH_K].tau=0;
  Seg[0].Ichan[HH_K].phi=0;
  Seg[0].Ichan[HH_K].val=Seg[0].Ichan[HH_K].gMaxTot*pow(Seg[0].Ichan[HH_K].n,4)*(Seg[0].Ichan[HH_K].Erev-Seg[0].v);
  //HH sodium channel
  Seg[0].Ichan[HH_Na].gMax=1200;
  Seg[0].Ichan[HH_Na].gMaxTot=Seg[0].Ichan[HH_Na].gMax*Pi*Seg[0].diameter*Seg[0].length;
  Seg[0].Ichan[HH_Na].Erev=0.050;
  Seg[0].Ichan[HH_Na].n=0;
  Seg[0].Ichan[HH_Na].m=hh_alpha_m(Seg[0].Erest)/(hh_alpha_m(Seg[0].Erest)+hh_beta_m(Seg[0].Erest));
  Seg[0].Ichan[HH_Na].h=hh_alpha_h(Seg[0].Erest)/(hh_alpha_h(Seg[0].Erest)+hh_beta_h(Seg[0].Erest));
  Seg[0].Ichan[HH_Na].i=0;
  Seg[0].Ichan[HH_Na].Na_att=0;
  Seg[0].Ichan[HH_Na].alpha=0;
  Seg[0].Ichan[HH_Na].beta=0;
  Seg[0].Ichan[HH_Na].k=0;
  Seg[0].Ichan[HH_Na].tau=0;
  Seg[0].Ichan[HH_Na].phi=0;
  Seg[0].Ichan[HH_Na].val=Seg[0].Ichan[HH_Na].gMaxTot*pow(Seg[0].Ichan[HH_Na].m,3)*Seg[0].Ichan[HH_Na].h*(Seg[0].Ichan[HH_Na].Erev-Seg[0].v);
  //Urakubo L-type Ca2+ channel
  Seg[0].Ichan[UKB_CA_LTYPE].gMax=930;
  Seg[0].Ichan[UKB_CA_LTYPE].Erev=0;
  Seg[0].Ichan[UKB_CA_LTYPE].m=(1/(1+exp(-1000*(Seg[0].v-0.037))))*(1/(3.6e-3));
  Seg[0].Ichan[UKB_CA_LTYPE].h=(1/(1+exp((Seg[0].v+0.041)/0.0005)))*(1/(2.9e-2));
  Seg[0].Ichan[UKB_CA_LTYPE].n=0;
  Seg[0].Ichan[UKB_CA_LTYPE].i=0;
  Seg[0].Ichan[UKB_CA_LTYPE].Na_att=0;
  Seg[0].Ichan[UKB_CA_LTYPE].k=0;
  Seg[0].Ichan[UKB_CA_LTYPE].alpha=0;
  Seg[0].Ichan[UKB_CA_LTYPE].beta=0;
  Seg[0].Ichan[UKB_CA_LTYPE].tau=0;
  Seg[0].Ichan[UKB_CA_LTYPE].phi=75.6;
  Seg[0].Ichan[UKB_CA_LTYPE].val=Seg[0].Ichan[UKB_CA_LTYPE].gMax*pow(Seg[0].Ichan[UKB_CA_LTYPE].m,3)*Seg[0].Ichan[UKB_CA_LTYPE].h*(-Seg[0].Erest/(1-exp(Seg[0].Ichan[UKB_CA_LTYPE].phi*Seg[0].Erest))); 
  
  //the soma is a special case, let's assign
  //the rest of the dendrite compartments
  for(i=1;i<size;i++)
    {
      //Create the soma compartment (special case)
      //Seg[i].Erest=-0.065;//mV
      Seg[i].Erest=-0.065;
      Seg[i].v=Seg[i].Erest;//mV
      Seg[i].length=20e-6;//meters <-Urakubo dendrite
      Seg[i].diameter=2e-6;//meters <-Urakubo dendrite
      Seg[i].distance=i*Seg[i].length;//distance from soma
      Seg[i].injectionTime=0.0;//seconds
      Seg[i].injectionDuration=0.0;//seconds
      Seg[i].Iinj=0;//amperes/meter
      Seg[i].Ra=1.5;//ohms*meter <-Urakubo dendrite
      Seg[i].RaTot=Seg[i].length*(4*Seg[i].Ra)/(Pi*Seg[i].diameter*Seg[i].diameter);//Ohms
      Seg[i].ra=Seg[i].RaTot/Seg[i].length;
      //Seg[i].Rm=1.4;//ohms*meter^2 <- Urakubo dendrite
      Seg[i].Rm=2.8;//new test
      Seg[i].RmTot=Seg[i].Rm/(Pi*Seg[i].diameter*Seg[i].length);//ohms
      Seg[i].rm=Seg[i].RmTot*Seg[i].length;//ohms*meter
      //Seg[i].Cm=0.02;//farads/meter^2 <- Urakubo dendrite
      
      Seg[i].Cm=0.01;//new test
 
      Seg[i].CmTot=Seg[i].Cm*Pi*Seg[i].diameter*Seg[i].length;//farads
      Seg[i].cm=Seg[i].CmTot/Seg[i].length;//farads/meter
      Seg[i].tau=Seg[i].rm*Seg[i].cm;//seconds
      Seg[i].lambda=sqrt(Seg[i].rm/Seg[i].ra);//meters
      //channel initializations in dendrites
      //leak channel
      //Seg[i].Ichan[HH_Leak].gMax=0.7;//<--- Old error
      Seg[i].Ichan[HH_Leak].gMax=1/Seg[i].Rm;
      Seg[i].Ichan[HH_Leak].gMaxTot=Seg[i].Ichan[HH_Leak].gMax*Pi*Seg[i].diameter*Seg[i].length;
      Seg[i].Ichan[HH_Leak].Erev=-0.065;
      Seg[i].Ichan[HH_Leak].n=0;
      Seg[i].Ichan[HH_Leak].m=0;
      Seg[i].Ichan[HH_Leak].h=0;
      Seg[i].Ichan[HH_Leak].i=0;
      Seg[i].Ichan[HH_Leak].Na_att=0;
      Seg[i].Ichan[HH_Leak].alpha=0;
      Seg[i].Ichan[HH_Leak].beta=0;
      Seg[i].Ichan[HH_Leak].k=0;
      Seg[i].Ichan[HH_Leak].tau=0;
      Seg[i].Ichan[HH_Leak].phi=0;
      Seg[i].Ichan[HH_Leak].val=Seg[i].Ichan[HH_Leak].gMaxTot*(Seg[i].Ichan[HH_Leak].Erev-Seg[i].v);
      //HH potassium channel
      Seg[i].Ichan[HH_K].gMax=360;
      Seg[i].Ichan[HH_K].gMaxTot=Seg[i].Ichan[HH_K].gMax*Pi*Seg[i].diameter*Seg[i].length;
      Seg[i].Ichan[HH_K].Erev=-0.077;
      Seg[i].Ichan[HH_K].n=hh_alpha_n(Seg[i].Erest)/(hh_alpha_n(Seg[i].Erest + hh_beta_n(Seg[i].Erest)));
      Seg[i].Ichan[HH_K].m=0;
      Seg[i].Ichan[HH_K].h=0;
      Seg[i].Ichan[HH_K].i=0;
      Seg[i].Ichan[HH_K].Na_att=0;
      Seg[i].Ichan[HH_K].alpha=0;
      Seg[i].Ichan[HH_K].beta=0;
      Seg[i].Ichan[HH_K].k=0;
      Seg[i].Ichan[HH_K].tau=0;
      Seg[i].Ichan[HH_K].phi=0;
      Seg[i].Ichan[HH_K].val=Seg[i].Ichan[HH_K].gMaxTot*pow(Seg[i].Ichan[HH_K].n,4)*(Seg[i].Ichan[HH_K].Erev-Seg[i].v);
      //HH sodium channel
      //no HH sodium channels in dendrites?
      Seg[i].Ichan[HH_Na].gMax=1200;
      Seg[i].Ichan[HH_Na].gMaxTot=Seg[i].Ichan[HH_Na].gMax*Pi*Seg[i].diameter*Seg[i].length;
      Seg[i].Ichan[HH_Na].Erev=0.050;
      Seg[i].Ichan[HH_Na].n=0;
      Seg[i].Ichan[HH_Na].m=hh_alpha_m(Seg[i].Erest)/(hh_alpha_m(Seg[i].Erest)+hh_beta_m(Seg[i].Erest));
      Seg[i].Ichan[HH_Na].h=hh_alpha_h(Seg[i].Erest)/(hh_alpha_h(Seg[i].Erest)+hh_beta_h(Seg[i].Erest));
      Seg[i].Ichan[HH_Na].i=0;
      Seg[i].Ichan[HH_Na].Na_att=0;
      Seg[i].Ichan[HH_Na].alpha=0;
      Seg[i].Ichan[HH_Na].beta=0;
      Seg[i].Ichan[HH_Na].k=0;
      Seg[i].Ichan[HH_Na].tau=0;
      Seg[i].Ichan[HH_Na].phi=0;
      Seg[i].Ichan[HH_Na].val=Seg[i].Ichan[HH_Na].gMaxTot*pow(Seg[i].Ichan[HH_Na].m,3)*Seg[i].Ichan[HH_Na].h*(Seg[i].Ichan[HH_Na].Erev-Seg[i].v);
       
      //Urakubo L-type Ca2+ channel
      if(i>0 && i<=2)
	{
	  Seg[i].Ichan[UKB_CA_LTYPE].gMax=1460;
	}
      else if(i>2)
	{
	  Seg[i].Ichan[UKB_CA_LTYPE].gMax=32; 
	}
      Seg[i].Ichan[UKB_CA_LTYPE].Erev=0;
      Seg[i].Ichan[UKB_CA_LTYPE].m=(1/(1+exp(-1000*(Seg[i].v-0.037))))*(1/(3.6e-3));
      Seg[i].Ichan[UKB_CA_LTYPE].h=(1/(1+exp((Seg[i].v+0.041)/0.0005)))*(1/(2.9e-2));
      Seg[i].Ichan[UKB_CA_LTYPE].n=0;
      Seg[i].Ichan[UKB_CA_LTYPE].i=0;
      Seg[i].Ichan[UKB_CA_LTYPE].Na_att=0;
      Seg[i].Ichan[UKB_CA_LTYPE].k=0;
      Seg[i].Ichan[UKB_CA_LTYPE].alpha=0;
      Seg[i].Ichan[UKB_CA_LTYPE].beta=0;
      Seg[i].Ichan[UKB_CA_LTYPE].tau=0;
      Seg[i].Ichan[UKB_CA_LTYPE].phi=75.6;
      Seg[i].Ichan[UKB_CA_LTYPE].val=Seg[i].Ichan[UKB_CA_LTYPE].gMax*pow(Seg[i].Ichan[UKB_CA_LTYPE].m,3)*Seg[i].Ichan[UKB_CA_LTYPE].h*(-Seg[i].Erest/(1-exp(Seg[i].Ichan[UKB_CA_LTYPE].phi*Seg[i].Erest))); 
    }
}

void update_HH_K_n(Cmpt Seg[], int size, double dt)
{
  int i;
  for(i=0;i<size;i++)
    {
      Seg[i].Ichan[HH_K].n=Seg[i].Ichan[HH_K].n+dt*hh_alpha_n(Seg[i].v)*(1-Seg[i].Ichan[HH_K].n)-dt*hh_beta_n(Seg[i].v)*Seg[i].Ichan[HH_K].n;
    }
}

void update_HH_Na_m(Cmpt Seg[], int size, double dt)
{
  int i;
  for(i=0;i<size;i++)
    {
      Seg[i].Ichan[HH_Na].m=Seg[i].Ichan[HH_Na].m+dt*hh_alpha_m(Seg[i].v)*(1-Seg[i].Ichan[HH_Na].m)-dt*hh_beta_m(Seg[i].v)*Seg[i].Ichan[HH_Na].m;
    }
}

void update_HH_Na_h(Cmpt Seg[], int size, double dt)
{
  int i;
  for(i=0;i<size;i++)
    {
      Seg[i].Ichan[HH_Na].h=Seg[i].Ichan[HH_Na].h+dt*hh_alpha_h(Seg[i].v)*(1-Seg[i].Ichan[HH_Na].h)-dt*hh_beta_h(Seg[i].v)*Seg[i].Ichan[HH_Na].h;
    }
}

void update_UKB_Ca_Ltype_m(Cmpt Seg[], int size, double dt)
{
  int i;
  double voltage;
  for(i=0;i<size;i++)
    {
      voltage=Seg[i].v*1000;
      Seg[i].Ichan[UKB_CA_LTYPE].m=dt*(1/(3.6e-3))*((1/(1+exp(-voltage-37)))-Seg[i].Ichan[UKB_CA_LTYPE].m)+Seg[i].Ichan[UKB_CA_LTYPE].m;
      //Seg[i].Ichan[UKB_CA_LTYPE].m=dt*((1/(1+exp(-1000*(Seg[i].v-0.037))))*(1/(3.6e-3))-Seg[i].Ichan[UKB_CA_LTYPE].m)+Seg[i].Ichan[UKB_CA_LTYPE].m;
    }
}
void update_UKB_Ca_Ltype_h(Cmpt Seg[], int size, double dt)
{
  int i;
  double voltage;
  for(i=0;i<size;i++)
    {
      voltage=Seg[i].v*1000;
      //Seg[i].Ichan[UKB_CA_LTYPE].h=dt*((1/(1+exp((Seg[i].v+0.041)/0.0005)))*(1/(2.9e-2))-Seg[i].Ichan[UKB_CA_LTYPE].h)+Seg[i].Ichan[UKB_CA_LTYPE].h;
      Seg[i].Ichan[UKB_CA_LTYPE].h=dt*(1/(2.9e-2))*((1/(1+exp((voltage+41)/0.5)))-Seg[i].Ichan[UKB_CA_LTYPE].h)+Seg[i].Ichan[UKB_CA_LTYPE].h;
    }
}


void update_leak(Cmpt Seg[], int size)
{
  int i;
  for(i=0;i<size;i++)
    {
      Seg[i].Ichan[HH_Leak].val=Seg[i].Ichan[HH_Leak].gMaxTot*(Seg[i].Ichan[HH_Leak].Erev-Seg[i].v); 
    }
}

void update_HH_K(Cmpt Seg[], int size)
{
  int i;
  for(i=0;i<size;i++)
    {
      Seg[i].Ichan[HH_K].val=Seg[i].Ichan[HH_K].gMaxTot*pow(Seg[i].Ichan[HH_K].n,4)*(Seg[i].Ichan[HH_K].Erev-Seg[i].v);
    }
}

void update_HH_Na(Cmpt Seg[], int size)
{
  int i;
  for(i=0;i<size;i++)
    {
      Seg[i].Ichan[HH_Na].val=Seg[i].Ichan[HH_Na].gMaxTot*pow(Seg[i].Ichan[HH_Na].m,3)*Seg[i].Ichan[HH_Na].h*(Seg[i].Ichan[HH_Na].Erev-Seg[i].v);
    }
}

void update_UKB_Ca_Ltype(Cmpt Seg[], int size)
{
  int i;
  for(i=0;i<size;i++)
    {
      Seg[i].Ichan[UKB_CA_LTYPE].val=Seg[i].Ichan[UKB_CA_LTYPE].gMax*pow(Seg[i].Ichan[UKB_CA_LTYPE].m,3)*Seg[i].Ichan[UKB_CA_LTYPE].h*(-Seg[i].Erest/(1-exp(Seg[i].Ichan[UKB_CA_LTYPE].phi*Seg[i].Erest)));
    }
}

void structTest(Cmpt Seg[])
{
  Seg[0].v=Seg[0].v+1e-10;
}


void update_v_Soma(Cmpt Seg[], int size, double dx, double dt, double t)
{
  Seg[0].v=(dt/Seg[0].CmTot)*(Seg[0].Ichan[HH_Leak].val+Seg[0].Iinj)+Seg[0].v;
  
  if(t>0.250 && t<0.250+dt)
    {
      /*
      printf("Soma voltage=%f\n",Seg[0].v);
      printf("dendrite cmptmt #1=%f\n",Seg[1].v);
      printf("dendrite cmptmt #2=%f\n",Seg[2].v);
      printf("dt=%f\n",dt);
      printf("dx=%f\n",dx);
      printf("Seg[0].diameter=%f\n",Seg[0].diameter);
      printf("Seg[0].length=%f\n",Seg[0].length);
      printf("Pi=%f\n",Pi);
      printf("Inject current=%.15f\n",Seg[0].Iinj);
      printf("Pi*Seg[0].diameter*Seg[0].length=%.15f\n",Pi*Seg[0].diameter*Seg[0].length);
      printf("Seg[0].Iinj/(Pi*Seg[0].diameter*Seg[0].length)=%.15f\n",Seg[0].Iinj/(Pi*Seg[0].diameter*Seg[0].length));
      printf("Seg[0].Ichan[0].val=%.15f\n",Seg[0].Ichan[0].val);
      */
    }
}

void update_v_Passive(Cmpt Seg[], int size, double dx, double dt, double t)
{
    int i;
  //update the voltage of each compartment
  for(i=0;i<size;i++)
    {
      if(i==0)
	{
	  Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i+1].v-Seg[i].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Iinj)+Seg[i].v;
	}
      else if(i==(size-1))
	{
	  Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i-1].v-Seg[i].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Iinj)+Seg[i].v;
	}
      else
	{
	  Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i-1].v-2*Seg[i].v+Seg[i+1].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Iinj)+Seg[i].v;
	}
    }
}

//works!
void update_v_HH_Soma(Cmpt Seg[], int size, double dx, double dt, double t)
{
  Seg[0].v=(dt/Seg[0].CmTot)*(Seg[0].Ichan[HH_Leak].val+Seg[0].Ichan[HH_Na].val+Seg[0].Ichan[HH_K].val+Seg[0].Iinj)+Seg[0].v;
}

void update_v_HH_Cell(Cmpt Seg[], int size, double dx, double dt, double t)
{
  int i;
  //update the voltage of each compartment
  for(i=0;i<size;i++)
    {
      if(i==0)
	{
	  //with HH channels
	  Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i+1].v-Seg[i].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Ichan[HH_Na].val+Seg[i].Ichan[HH_K].val+Seg[i].Iinj)+Seg[i].v;
	  //without HH channels
	  //Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i+1].v-Seg[i].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Iinj)+Seg[i].v;
	}
      else if(i==(size-1))
	{
	  //with HH channels
	  //Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i-1].v-Seg[i].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Ichan[HH_Na].val+Seg[i].Ichan[HH_K].val+Seg[i].Iinj)+Seg[i].v;
	  //without HH channels
	  Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i-1].v-Seg[i].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Iinj)+Seg[i].v;
	}
      else
	{
	  //with HH channels
	  //Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i-1].v-2*Seg[i].v+Seg[i+1].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Ichan[HH_Na].val+Seg[i].Ichan[HH_K].val+Seg[i].Iinj)+Seg[i].v;
	  //without HH channels
	  Seg[i].v=(dt/Seg[i].CmTot)*((Seg[i-1].v-2*Seg[i].v+Seg[i+1].v)/Seg[i].RaTot+Seg[i].Ichan[HH_Leak].val+Seg[i].Iinj)+Seg[i].v;
	}
    }
  if(t>0.250 && t<0.250+dt)
    {
      /*
      printf("Soma voltage=%f\n",Seg[0].v);
      printf("dendrite cmptmt #1=%f\n",Seg[1].v);
      printf("dendrite cmptmt #2=%f\n",Seg[2].v);
      printf("dt=%f\n",dt);
      printf("dx=%f\n",dx);
      printf("Seg[0].diameter=%f\n",Seg[0].diameter);
      printf("Seg[0].length=%f\n",Seg[0].length);
      printf("Pi=%f\n",Pi);
      printf("Inject current=%.15f\n",Seg[0].Iinj);
      printf("Pi*Seg[0].diameter*Seg[0].length=%.15f\n",Pi*Seg[0].diameter*Seg[0].length);
      printf("Seg[0].Iinj/(Pi*Seg[0].diameter*Seg[0].length)=%.15f\n",Seg[0].Iinj/(Pi*Seg[0].diameter*Seg[0].length));
      printf("Seg[0].Ichan[0].val=%f\n",Seg[0].Ichan[0].val);
      */
    }
}

double i_NMDA_Ca(double v, double t, double gRelease, double gMax)
{
  double val;
  double g_NMDA;
  double v_Cal;
  double NMDA_CaFrac;
  double p0;
  double MgConc;
  
  //g_NMDA=175.0;
  g_NMDA=gMax;
  v_Cal=0.130;//the reversal potential of Calicum (Ca^2+) not v_NMDA which is 0 mV
  MgConc=1.5e-3;
  p0=0.5;
  //NMDA_CaFrac=0.10;

  //val=g_NMDA*NMDA_CaFrac*p0*glutBind(t,gRelease)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  val=g_NMDA*p0*glutBind(t,gRelease)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  return(val);
}

double i_NMDA_Ca_Two(double v, double t, double gRelease1, double gRelease2, double gMax)
{
  double val;
  double g_NMDA;
  double v_Cal;
  double NMDA_CaFrac;
  double p0;
  double MgConc;
  
  //g_NMDA=175.0;
  g_NMDA=gMax;
  v_Cal=0.130;//the reversal potential of Calicum (Ca^2+) not v_NMDA which is 0 mV
  MgConc=1.5e-3;
  p0=0.5;
  //NMDA_CaFrac=0.10;

  //val=g_NMDA*NMDA_CaFrac*p0*glutBindTwo(t,gRelease1,gRelease2)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  val=g_NMDA*p0*glutBindTwo(t,gRelease1,gRelease2)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  return(val);
}

double i_NMDA_Ca_Two_Un_Sat_F_and_S(double v, double t, double gRelease1, double gRelease2, double gMax)
{
  double val;
  double g_NMDA;
  double v_Cal;
  double NMDA_CaFrac;
  double MgConc;
  
  //g_NMDA=175.0;
  g_NMDA=gMax;
  v_Cal=0.130;//the reversal potential of Calicum (Ca^2+) not v_NMDA which is 0 mV
  MgConc=1.5e-3;
  //NMDA_CaFrac=0.10;

  //val=g_NMDA*NMDA_CaFrac*glutBindTwoUnSatFandS(t,gRelease1,gRelease2)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  val=g_NMDA*glutBindTwoUnSatFandS(t,gRelease1,gRelease2)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  return(val);
}

double i_NMDA_Ca_Two_Un_Sat_No_Fast(double v, double t, double gRelease1, double gRelease2, double gMax)
{
  double val;
  double g_NMDA;
  double v_Cal;
  double NMDA_CaFrac;
  double p0;
  double MgConc;
  
  //g_NMDA=175.0;
  g_NMDA=gMax;
  v_Cal=0.130;//the reversal potential of Calicum (Ca^2+) not v_NMDA which is 0 mV
  MgConc=1.5e-3;
  p0=0.5;//take this out
  //NMDA_CaFrac=0.10;

  val=g_NMDA*glutBindTwoUnSatNoFast(t,gRelease1,gRelease2)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  return(val);
}

double i_NMDA_Ca_Two_Diff(double v, double t, double gRelease1, double gRelease2, double gMax, double tauCaDiff)
{
  double val;
  double g_NMDA;
  double v_Cal;
  double NMDA_CaFrac;
  double p0;
  double MgConc;
  
  //g_NMDA=175.0;
  g_NMDA=gMax;
  v_Cal=0.130;//the reversal potential of Calicum (Ca^2+) not v_NMDA which is 0 mV
  MgConc=1.5e-3;
  p0=0.5;
  //NMDA_CaFrac=0.10;

  //val=g_NMDA*NMDA_CaFrac*p0*glutBindTwoDiff(t,gRelease1,gRelease2,tauCaDiff)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  val=g_NMDA*p0*glutBindTwoDiff(t,gRelease1,gRelease2,tauCaDiff)*(v_Cal-v)/(1+(MgConc/3.57e-3)*exp(-62*v));
  return(val);
}

double i_VDCC_simple(double v, double t, double injectionTime)
{
  double val;
  double g_VDCC;
  double v_EVDCC;

  g_VDCC=8;
  v_EVDCC=0.130;

  if(t>injectionTime && t<=(injectionTime+0.002) && (v>=-0.030))
    {
    val=g_VDCC*(v_EVDCC-v);//reversal potential of 130 mV to match NMDA. Could also be 0 (Urakubo2008 et al)
    }
  else
    val=0;
  
  return(val);
}

//this is the same as the i_VDCC_simple function with double the 
//conductance and half the open time
double i_VDCC_short(double v, double t, double injectionTime)
{
  double val;
  double g_VDCC;
  double v_EVDCC;
  double volFrac;
  
  volFrac=0.10;
  g_VDCC=16;
  v_EVDCC=0.130;

  if(t>injectionTime && t<=(injectionTime+0.001) && (v>=-0.030))
    {
      val=g_VDCC*(1/volFrac)*(v_EVDCC-v);
    }
  else
    val=0;
  
  return(val);
}

//this is the same as the i_VDCC_simple function with double the 
//conductance and half the open time
double i_VDCC_short_update(double v, double t, double injectionTime, double gMax)
{
  double val;
  double g_VDCC;
  double v_EVDCC;
  
  //original
  g_VDCC=gMax;
  v_EVDCC=0.130;

  if(t>injectionTime && t<=(injectionTime+0.001) && (v>=-0.030))
    {
      val=g_VDCC*(v_EVDCC-v);
    }
  else
    val=0;
  
  return(val);
}

double i_VDCC_short_two_update(double v, double t, double injectionTime1, double injectionTime2, double gMax)
{
  double val;
  double g_VDCC;
  double v_EVDCC;
  
  g_VDCC=gMax;
  v_EVDCC=0.130;

  if(t>injectionTime1 && t<=(injectionTime1+0.001) && (v>=-0.030))
    {
      val=g_VDCC*(v_EVDCC-v);
    }
  else if(t>injectionTime2 && t<=(injectionTime2+0.001) && (v>=-0.030))
    {
      val=g_VDCC*(v_EVDCC-v);
    }
  else
    val=0;
  
  return(val);
}

double VDCC_Ca_Diff(double t, double injectionTime, double Ca_VDCC, double Ca_value, double tauCaDiff)
{
  double val;
  if(t>injectionTime)
    {
      val=(Ca_VDCC-Ca_value)/tauCaDiff;
    }
  else
    {
      val=0;
    }


  return(val);
}


//Not sure if this works.  Need to test.
//will need this fuction for comparison's sake
double BPAP(double t,double injectionTime)
{
  double val;
  if(t<injectionTime)
    {
      //val=(-0.0578+0.058);
      //original
      val=-0.065;
    }
  //don't need triplet
  /*
  else if (t>=injectionTime && t<=(injectionTime+0.020))
    {
      //val=(-0.0578+0.058)+.100*(0.75*exp((injectionTime-t)/0.003)+0.25*exp((injectionTime-t)/0.025));
      //original
      val=-0.065+.100*(0.75*exp((injectionTime-t)/0.003)+0.25*exp((injectionTime-t)/0.025));
    }
  */
  else
    {
      val=-0.065+.100*(0.75*exp((injectionTime+0.020-t)/0.003)+0.25*exp((injectionTime+0.020-t)/0.025));
    }
  return(val);
}

//think about how to implement this. <--Specifically tp
double glutBind(double t, double gRelease)
{
  double val;
  double INMDAf;
  double INMDAs;
  double tauNMDAf;
  double tauNMDAs;
  double theta;

  INMDAf=0.5;
  INMDAs=0.5;
  tauNMDAf=50e-3;
  tauNMDAs=200e-3;

  if(t<gRelease)
    {
      theta=0;
    }
  else
    {
      theta=1;
    }
  
  val=theta*(INMDAf*exp((gRelease-t)/tauNMDAf)+INMDAs*exp((gRelease-t)/tauNMDAs));
  return(val);
}


double glutBindTwo(double t, double gRelease1, double gRelease2)
{
  double val;
  double INMDAf_total;
  double INMDAs_total;
  double tauNMDAf,tauNMDAs;
  INMDAf_total=0.5;
  INMDAs_total=0.5;
  tauNMDAf=50e-3;
  tauNMDAs=200e-3;

  if (t<gRelease1)
    {
      val=0;
      //printf("val=%f\n",val);
    }
  else if (gRelease1<=t && t<=gRelease2)
    {
      val=INMDAf_total*exp(-(t-gRelease1)/tauNMDAf)+INMDAs_total*exp(-(t-gRelease1)/tauNMDAs);
      //printf("val=%f\n",val);
    }
  else
    {
      val=INMDAf_total*exp(-(t-gRelease2)/tauNMDAf)+INMDAs_total*exp(-(t-gRelease2)/tauNMDAs);
    }
  return(val);
}


double glutBindTwoUnSatFandS(double t, double gRelease1, double gRelease2)
{
  double val;
  double fastFrac,slowFrac;
  double tauNMDAf,tauNMDAs;
  double Pb,Pf,Ps;
  double N1f,N1s,N2f,N2s,N1f_t2,N1s_t2;
  double Navail,Nbound;

  tauNMDAf=50e-3;
  tauNMDAs=200e-3;
  Pb=0.3;
  fastFrac=0.50;//edit this one only
  slowFrac=1-fastFrac;
  Pf=fastFrac*Pb;
  Ps=slowFrac*Pb;
   
  if (t<gRelease1)
    {
      val=0;
    }
  else if (gRelease1<=t && t<=gRelease2)
    {
      N1f=Pf*exp(-(t-gRelease1)/tauNMDAf);
      N1s=Ps*exp(-(t-gRelease1)/tauNMDAs);
      val=N1f+N1s;
    }
  else
    {
      N1f_t2=Pf*exp(-(gRelease2-gRelease1)/tauNMDAf);
      N1s_t2=Ps*exp(-(gRelease2-gRelease1)/tauNMDAs);
      Nbound=N1f_t2+N1s_t2;
      Navail=1-Nbound;
      N1f=Pf*exp(-(t-gRelease1)/tauNMDAf);
      N1s=Ps*exp(-(t-gRelease1)/tauNMDAs);
      N2f=Navail*Pf*exp(-(t-gRelease2)/tauNMDAf);
      N2s=Navail*Ps*exp(-(t-gRelease2)/tauNMDAs);
      val=N1f+N1s+N2f+N2s;
    }
  return(val);
}

double glutBindTwoUnSatNoFast(double t, double gRelease1, double gRelease2)
{
  double val;
  double INMDAf_total;
  double INMDAs_total;
  double tauNMDAf,tauNMDAs;
  double Pb;
  INMDAf_total=0.5;
  INMDAs_total=0.5;
  tauNMDAf=50e-3;
  tauNMDAs=200e-3;
  Pb=0.3;

  if (t<gRelease1)
    {
      val=0;
      //printf("val=%f\n",val);
    }
  else if (gRelease1<=t && t<=gRelease2)
    {
      //val=INMDAf_total*exp(-(t-gRelease1)/tauNMDAf)+INMDAs_total*exp(-(t-gRelease1)/tauNMDAs);
      //with tauNMDA_Slow
      val=Pb*exp(-(t-gRelease1)/tauNMDAs);
      //with tauNMDA_Fast instead
      //val=Pb*exp(-(t-gRelease1)/tauNMDAf);
    }
  else
    {
      //val=INMDAf_total*exp(-(t-gRelease2)/tauNMDAf)+INMDAs_total*exp(-(t-gRelease2)/tauNMDAs);
      //With TauNMDA_Slow
      val=Pb*exp(-(t-gRelease1)/tauNMDAs)+(1-Pb*exp(-(gRelease2-gRelease1)/tauNMDAs))*Pb*exp(-(t-gRelease2)/tauNMDAs);
      //with TauNMDA_Fast instead
      //val=Pb*exp(-(t-gRelease1)/tauNMDAf)+(1-Pb*exp(-(gRelease2-gRelease1)/tauNMDAf))*Pb*exp(-(t-gRelease2)/tauNMDAf);
    }
  return(val);
}


double glutBindTwoDiff(double t, double gRelease1, double gRelease2, double tauCaDiff)
{
  double val;
  double INMDAf_total;
  double INMDAs_total;
  double tauNMDAf,tauNMDAs;
  INMDAf_total=0.5;
  INMDAs_total=0.5;
  tauNMDAf=50e-3;
  tauNMDAs=200e-3;

  if (t<(gRelease1+tauCaDiff))
    {
      val=0;
      //printf("val=%f\n",val);
    }
  else if ((gRelease1+tauCaDiff)<=t && t<=(gRelease2+tauCaDiff))
    {
      //printf("got to glutBindTwoDiff, gRelease1\n");
      //val=(1-exp((tauCaDiff+gRelease1-t)/tauCaDiff))*INMDAf_total*exp(-(t-gRelease1)/tauNMDAf)+INMDAs_total*exp(-(t-gRelease1)/tauNMDAs);
      val=INMDAf_total*exp(-(t-gRelease1-tauCaDiff)/tauNMDAf)+INMDAs_total*exp(-(t-gRelease1-tauCaDiff)/tauNMDAs);
      
      if(val<0)
	{
	  val=0;
	}
    }
  else
    {
      //printf("got to glutBindTwoDiff, gRelease2\n");
      //val=(1-exp((tauCaDiff+gRelease2-t)/tauCaDiff))*INMDAf_total*exp(-(t-gRelease2)/tauNMDAf)+INMDAs_total*exp(-(t-gRelease2)/tauNMDAs);
      
      val=INMDAf_total*exp(-(t-gRelease2-tauCaDiff)/tauNMDAf)+INMDAs_total*exp(-(t-gRelease2-tauCaDiff)/tauNMDAs);
      
      if(val<0)
	val=0;

    }
  return(val);
}

double i_UKB_Ca_Ltype(Cmpt Seg[],double t, double dt)
{
  double voltage;
  voltage=Seg[3].v*1000;
   
  //Reference third segment
  //Check for singularities
  if((1-exp(Seg[3].Ichan[UKB_CA_LTYPE].phi*Seg[3].v))<1e-6)
    {
      return(Seg[3].Ichan[UKB_CA_LTYPE].gMax*pow(Seg[3].Ichan[UKB_CA_LTYPE].m,3.0)*Seg[3].Ichan[UKB_CA_LTYPE].h*(1/(Seg[3].Ichan[UKB_CA_LTYPE].phi)));
    }
  else
    {
      return(Seg[3].Ichan[UKB_CA_LTYPE].val=Seg[3].Ichan[UKB_CA_LTYPE].gMax*pow(Seg[3].Ichan[UKB_CA_LTYPE].m,3)*Seg[3].Ichan[UKB_CA_LTYPE].h*((Seg[3].Ichan[UKB_CA_LTYPE].Erev-Seg[3].v)/(1-exp(Seg[3].Ichan[UKB_CA_LTYPE].phi*Seg[3].v))));
    }
}


