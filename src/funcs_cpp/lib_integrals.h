#ifndef integrals_includes
#define integrals_includes
//C/C++ includes

//Custom includes
#include "../utils_cpp/lib_includes.h"
#include "../scf_main/scf_classes.h"
#include "../classes/classes.h"

// --------------------------------------------------------------------- //
// 			DEFINE MATH FUNCTIONS				 //
// --------------------------------------------------------------------- //

//Calculated normalization factor
double norm(int ax,int ay,int az,double an);
double norm2(int ax,int ay,int az,double an);

//Calculates Factorial of Val
double fact(int val);

//Calculates double factorial of val
double d_fact(int val);

//Calculates n choose k
double nCk(int n, int k);

//Calculates the E_factor from the gaussian product rule
double E_factor(double an, double bm, double AmB2);

//Calculates the individual components of the the overlap integrals
double S_factor(int axyz,int bxyz,double an,double bm,double Axyz,double Bxyz);

//Calculates the error function of x w/ an expansion truncated at the nth term
double erf(double x,int n);

//Calculate Fnu
double Fnu(int nu,double x);

//C factor from the solution of the nuc attract integral
double C_fact(double an,double bm);

//NRI factor from the solution of the electron repulsion integral
double NRI_fact(double gammaP,double gammaQ);

double fi(int in,int l1,int l2,double x1, double x2);

//Transition Dipole Numerical Integrator
double TDNumerInt(int l,double gamma,double P);

// --------------------------------------------------------------------- //
//                      DEFINE INTEGRAL FUNCTIONS                        //
// --------------------------------------------------------------------- //

//Overlap intgral main, fills the array
void overlap_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO);
void kinetic_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO);
void nuclear_main(int ac,MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO);
void electron_repulsion_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO);
void TD_main(MemHandler *data_mem,scf_data *scf,dataOutput *optfile,STOnG_handler *STO);
#endif

