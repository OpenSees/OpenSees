///*
//################################################################################
//# COPY-YES  (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:                                                                     #
//# CLASS:                                                                       #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0                                            #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              May 2004                                                  #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/
#include "time.h"
#include "limits.h"
#include "math.h"
#include "float.h"
// nDarray tools
#include "BJtensor.h"
#include "stresst.h"
#include "straint.h"
#include "BJvector.h"
#include "BJmatrix.h"
//#include "iostream.h"

#include <Tensor.h>


int main(void)
{
::printf("\n\n-------------  MACHINE (CPU) DEPENDENT THINGS  --------------\n\n");

// defining machine epsilon for different built in data types supported by C++
float        float_eps       = f_macheps();  // or use float.h values FLT_EPSILON
double       double_eps      = d_macheps();  // or use float.h values DBL_EPSILON
//long double  long_double_eps = ld_macheps(); // or use float.h values LDBL_EPSILON

::printf("\n float macheps = %.20e \n",float_eps);
::printf(" double macheps = %.20e \n",double_eps);
//::printf(" long double macheps = %.20Le \n\n\n",long_double_eps);


// Usual tolerance is defined as square root of machine epsilon
float        sqrt_f_eps       = sqrt(f_macheps());
double       sqrt_d_eps       = sqrt(d_macheps());
//long double  sqrt_ld_eps      = sqrt(ld_macheps());

::printf("\n sqrt float macheps = %.20e \n",sqrt_f_eps);
::printf(" sqrt double macheps = %.20e \n",sqrt_d_eps);
//::printf(" sqrt long double macheps = %.20Le \n",sqrt_ld_eps);


double one = 1.0;
double onepmc = 1.0+double_eps;
::printf("\n one = %.32e \n",one);
::printf("onepmc = %.32e \n",onepmc);
if (one == onepmc) 
  {
    printf("1.0 == 1.0+double_eps -->>  OK \n");
  }
else 
  {
    printf("1.0 == 1.0+double_eps -->>  NOT OK \n");
  }   



double one2 = 10e20;
double onepmc2 = 10e20*(1.0+double_eps);
::printf("\n one2 = %.32e \n",one2);
::printf("onepmc2 = %.32e \n",onepmc2);
if (one2 == onepmc2) 
  {
    printf("10e20 == 10e20*(1.0+double_eps) -->>  OK \n");
  }
else 
  {
    printf("10e20 == 10e20*(1.0+double_eps) -->>  NOT OK \n");
  }   



double one3 = 1.0;
double onepmc3 = 1.0+double_eps*10000.0;
::printf("\n one3 = %.32e \n",one3);
::printf("onepmc3 = %.32e\n",onepmc3);
if (one3 == onepmc3) 
  {
    printf("1.0 == 1.0+double_eps*10.0 -->>  OK \n");
  }
else 
  {
    printf("1.0 == 1.0+double_eps*10.0 -->>  NOT OK \n");
  }   



::printf("\n\n----------------- TENSOR --------------------\n\n");
  
//  Tensor t0();
//  Tensor *ta = new Tensor;
  

//....................................................................
// There  are several built in tensor types. One of them is:
// Levi-Civita permutation tensor
  tensor e("e",3,def_dim_3);
  e.print("e","\nLevi-Civita permutation tensor e");

// The other one is:
// Kronecker delta tensor
  tensor I2("I", 2, def_dim_2);
  I2.print("I2","\ntensor I2 (2nd order unit tensor: Kronecker Delta -> d_ij)");



	 tensor Ipmnq = I2("pm") * I2("nq");
  tensor Imnpq = I2("mn") * I2("pq");



////  tensor Apqmn = alpha("pq") * I2("mn");
//  tensor X  = (1.0/3.0) * Imnpq ;// - Apqmn * (1.0/3.0);
//	 X.print();

// some of the 4. order tensor used in conituum mechanics:
// the most general representation of the fourth order isotropic tensor
// includes the following fourth orther unit isotropic tensors:
  tensor I_ijkl( 4, def_dim_4, 0.0 );
  I_ijkl = I2("ij")*I2("kl");
//  I_ijkl.print("ijkl","\ntensor I_ijkl = d_ij d_kl)");
  tensor I_ikjl( 4, def_dim_4, 0.0 );
//  I_ikjl = I2("ij")*I2("kl");
  I_ikjl = I_ijkl.transpose0110();
//  I_ikjl.print("ikjl","\ntensor I_ikjl) ");
  tensor I_ikjl_inv_2 = I_ikjl.inverse_2();
//  I_ikjl_inv_2.print("I_ikjl","\ntensor I_ikjl_inverse_2)");
  if ( I_ikjl == I_ikjl_inv_2 ) ::printf("\n\n I_ikjl == I_ikjl_inv_2 !\n");
  else ::printf("\a\n\n           I_ikjl != I_ikjl_inv_2 !\n");

  tensor I_iljk( 4, def_dim_4, 0.0 );
//  I_ikjl = I2("ij")*I2("kl");
//..  I_iljk = I_ikjl.transpose0101();
  I_iljk = I_ijkl.transpose0111();
//  I_iljk.print("iljk","\ntensor I_iljk) ");
// To check out this three fourth-order isotropic tensor please see
// W.Michael Lai, David Rubin, Erhard Krempl
// " Introduction to Continuum Mechanics"
// QA808.2
// ISBN 0-08-022699-X

// tensor additions ans scalar multiplications
// symmetric part of fourth oder unit isotropic tensor
  tensor I4s = (I_ikjl+I_iljk)*0.5;
//  tensor I4s = 0.5*(I_ikjl+I_iljk);
  I4s.print("I4s","\n symmetric tensor I4s = (I_ikjl+I_iljk)*0.5 ");
// skew-symmetric part of fourth oder unit isotropic tensor
  tensor I4sk = (I_ikjl-I_iljk)*0.5;
  I4sk.print("I4sk","\n skew-symmetric tensor I4sk = (I_ikjl-I_iljk)*0.5 ");

//....................................................................
// Linear Isotropic Elasticity Tensor
// building elasticity tensor
  double Ey = 20000;
  double nu = 0.2;
// building stiffness tensor in one command line
  tensor Eel = (I_ijkl*((Ey*nu*2.0)/(2.0*(1.0+nu)*(1-2.0*nu)))) +
              ( (I_ikjl+I_iljk)*(Ey/(2.0*(1.0+nu))));
// building compliance tensor in one command line
  tensor Del= (I_ijkl * (-nu/Ey)) + ( (I_ikjl+I_iljk) * ((1.0+nu)/(2.0*Ey)));

 
double PIo3 = PI/3.0;
//::printf("Pi/3.0 = %f \n",PIo3);
//-------------------------------------------
  double p_start     = 5000.000;
  double q_start     = 2000.0;
  double theta_start = PIo3;
//..//  double theta_start = 0.4;

stresstensor stress = stress.pqtheta2stress(p_start, q_start, theta_start);
stress.report("stress\n");


straintensor strain = D2("ijkl") * stress("kl");
strain.report("strain = D2(\"ijkl\") * stress(\"kl\") \n");






::printf("\n\n\n\n\n");
::printf("\n\n----------------- Derivatives of Yield function --------\n\n");
  


	 static double stressvalues[] = {  1.0, 2.0, 3.0,
                                    2.0, 2.0, 5.0,
                                    3.0, 5.0, 9.0  };

  stresstensor start_sigma(stressvalues);
  start_sigma.report("\n\n starting stress tensor start_sigma )");



	 static double strainvalues[] = {  0.1, 0.0, 0.0,
                                    0.0, 0.1, 0.0,
                                    0.0, 0.0, 0.1  };
  straintensor epsilon(strainvalues);
  epsilon.report("\n\n strain tensor epsilon (2nd-order with values assignment)");


  stresstensor stress_increment = E("ijpq") * strain_incr("pq");
  stress_increment.null_indices();
  stress_increment.report("\n\n stress increment"); 
  
  double f_start = 0.0;
  double f_pred  = 0.0;

  stresstensor elastic_predictor_stress = start_stress + stress_increment;
  elastic_predictor_stress.report("\n\n elastic predictor stress"); 
  

// Yield surface evaluation for von Mises
//    double k = 10.0: //cohesion
//    double k2 = k * k;
//    stresstensor s = elastic_predictor_stress.deviator();     
//    stresstensor temp1 = s("ij") * s("ij");
//    double temp = temp1.trace();
//    temp = temp * 3.0 / 2.0;
//
//    double f   = temp - k2;


//    f_start = 
    //::printf("\n##############  f_start = %.10e  ",f_start);

//    f_pred =  
    //::printf("##############  f_pred = %.10e\n\n",f_pred);


    stresstensor dFods( 2, def_dim_2, 0.0);
    //  stresstensor s;  // deviator
    tensor temp1( 2, def_dim_2, 0.0);
    tensor temp2( 2, def_dim_2, 0.0);
    double lower = 0.0;
    tensor temp3( 2, def_dim_2, 0.0);

    double Delta_lambda = 0.0;


    stresstensor plastic_stress;
    straintensor plastic_strain;
    stresstensor elastic_plastic_stress;
    // ::printf("\n...... felpred = %lf\n",felpred);

    if ( f_pred >= 0 ) 
      {

//        dFods = 

//  tensor upperE1 = Eel("pqkl")*dQods("kl");
//        upperE1.null_indices();
//  tensor upperE2 = dFods("ij")*Eel("ijmn");
//        upperE2.null_indices();
//
//  tensor upperE = upperE1("pq") * upperE2("mn");
//        upperE.null_indices();
//
//
//      temp2 = upperE2("ij")*dQods("ij"); // L_ij * E_ijkl * R_kl
//        temp2.null_indices();
//        lower = temp2.trace();
//
//        tensor Epl = upperE*(1./lower);
//
//        tensor Eep =  Eel - Epl;


    }

::printf("\n\n\n");



// --------------

  exit(1);
//  return 1;
}
