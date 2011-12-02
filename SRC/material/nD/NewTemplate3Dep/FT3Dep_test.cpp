///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, 
// DATE:              Fall 2005
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//


#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

#include "NewTemplate3Dep.h"

#include "ElasticState.h"
#include "Isotropic_Elastic.h"
#include "elnp_Elastic.h"
#include "DM04_Elastic.h"

#include "YieldFunction.h"
#include "DP_YF.h"
#include "VM_YF.h"
#include "CC_YF.h"
#include "DM04_YF.h"

#include "PlasticFlow.h"
#include "DP_PF.h"
#include "VM_PF.h"
#include "CC_PF.h"
#include "DM04_PF.h"

#include "ScalarEvolution.h"
#include "Linear_Eeq.h"
#include "CC_Ev.h"

#include "TensorEvolution.h"
#include "Linear_Eij.h"
#include "AF_Eij.h"
#include "DM04_alpha_Eij.h"
#include "DM04_z_Eij.h"

#include <G3Globals.h>
#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>
#include <OPS_Stream.h>

ErrorHandler *g3ErrorHandler =0;
double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;

OPS_Stream *opserrPtr;

int main() {

cout << "\n\n\n*** EP Finite Deformations: T E S T ***\n\n" << "\n";

ofstream outStress ("Results.txt");



double rho = 0.0;
double e0 = 0.735;
double G0 = 125;
double v = 0.05;
double Pat = 100.0;
double kc = 1.0;
double M_cal = 1.25;
double c = 0.712;
double lambda_c = 0.019;
double xi = 0.7;
double er = 0.934;
double m = 0.01;
double h0 = 7.05;
double ch = 0.968;
double nb = 1.1;
double A0 = 0.704;
double nd = 3.5;
double z_max = 4.0;
double cz = 600.0;
stresstensor zeroT1;
stresstensor zeroT2;
stresstensor initStress;
double p = 1000.0;
initStress.val(1,1) = -p; initStress.val(2,2) = -p; initStress.val(3,3) = -p;
//----------------1    2   3  4  5    6   7      8  9         10  11 12 13   14  15  16  17  18    19
double MC[19] = {rho, e0, G0, v, Pat, kc, M_cal, c, lambda_c, xi, er, m, h0, ch, nb, A0, nd, z_max, cz};
stresstensor TS[2] = {zeroT1, zeroT2};
MaterialParameter matpar(MC,19, TS,2);
DM04_Elastic le(3, 4, 5, 6, 2, initStress);
DM04_YF dpy(0,12, 2,1);
DM04_PF dpf(0, 2, 0,11, 0,9, 0,10, 0,5, 0,12, 0,7, 0,8, 0,16, 0,17, 2,1, 2,2 );
DM04_alpha_Eij Eij1(2, 11, 9, 10, 5, 12, 7, 8, 15, 13, 14, 3, 1, 2 );
DM04_z_Eij Eij2(12, 19, 18, 1, 2 );
TensorEvolution *TE[2] = {&Eij1, &Eij2};

NewTemplate3Dep FTEP(1, &matpar, &le, &dpy, &dpf, TE);

stresstensor thisE;
stresstensor tStress;
straintensor tStrain;

double E11 = 0.0;    double E12 = 0.0;    double E13 = 0.0;
double E21 = 0.0;    double E22 = 0.0;    double E23 = 0.0;
double E31 = 0.0;    double E32 = 0.0;    double E33 = 0.0;

int Num_cyc = 30;
int I_cyc = 0;

//double ee = 0.0;
double pp = 0.0;
double qq = 0.0;

double  q_cut = 300.0;
//double ep_cut = 0.015;

double d_g = 0.00005;
 
 while ( I_cyc < Num_cyc) {
    
    //undrained tri-axial compression 
    E11 += (- d_g);    E22 += (d_g*0.5);    E33 += (d_g*0.5);
    
    thisE.val(1,1) = E11; thisE.val(1,2) = E12; thisE.val(1,3) = E13;
    thisE.val(2,1) = E21; thisE.val(2,2) = E22; thisE.val(2,3) = E23;
    thisE.val(3,1) = E31; thisE.val(3,2) = E32; thisE.val(3,3) = E33;

    FTEP.setTrialStrain(thisE);

    tStress = FTEP.getStressTensor();
    tStrain = FTEP.getStrainTensor();
    
    pp  = tStress.p_hydrostatic();
    qq = - tStress.val(1,1) + tStress.val(2,2);

    outStress << -tStrain.val(1,1)  <<  "  "  <<  pp << "  " << qq << endl;

    //if (fabs(E11) >= ep_cut) {
    if (fabs(qq) >= q_cut) {
        I_cyc++;
        d_g *= (-1.0);
    }

 }

  return 0;

}


/*
double rho = 0.0;
double e0 = 0.735;
double G0 = 125;
double v = 0.25;
double Pat = 100.0;
double kc = 0.1;
double M_cal = 1.25;
double c = 0.712;
double lambda_c = 0.019;
double xi = 0.7;
double er = 0.934;
double m = 0.01;
double h0 = 7.05;
double ch = 0.968;
double nb = 1.1;
double A0 = 0.704;
double nd = 3.5;
double z_max = 4.0;
double cz = 600.0;
stresstensor zeroT1;
stresstensor zeroT2;
stresstensor initStress;
double p = 100;
initStress.val(1,1) = -p; initStress.val(2,2) = -p; initStress.val(3,3) = -p;
//----------------1    2   3  4  5    6   7      8  9         10  11 12 13   14  15  16  17  18    19
double MC[19] = {rho, e0, G0, v, Pat, kc, M_cal, c, lambda_c, xi, er, m, h0, ch, nb, A0, nd, z_max, cz};
stresstensor TS[2] = {zeroT1, zeroT2};
MaterialParameter inpar(MC,19, TS,2);
DM04_Elastic le(3, 4, 5, 6, 2, initStress);
DM04_YF dpy(0,12, 2,1);
DM04_PF dpf(0,2, 0,11, 0,9, 0,10, 0,5, 0,12, 0,7, 0,8, 0,16, 0,17, 2,1, 2,2 );
DM04_alpha_Eij Eij1(2, 11, 9, 10, 5, 12, 7, 8, 15, 13, 14, 3, 1, 2 );
DM04_z_Eij Eij2(12, 19, 18, 1, 2 );
TensorEvolution *TE[2] = {&Eij1, &Eij2};

NewTemplate3Dep FTEP(1, &inpar, &le, &dpy, &dpf, TE);

stresstensor thisE;
stresstensor tStress;
straintensor tStrain;

double E11 = 0.0;    double E12 = 0.0;    double E13 = 0.0;
double E21 = 0.0;    double E22 = 0.0;    double E23 = 0.0;
double E31 = 0.0;    double E32 = 0.0;    double E33 = 0.0;

  //double ee = 0.0;
  double pp = 0.0;
  double qq = 0.0;

int Num_cyc = 20;
int I_cyc = 0;

//double  q_cut = 5.0;
double ep_cut = 0.0005;  
  
double d_g = 0.000001;
    
//for (int i = 0; i < 600; i++) {
while ( I_cyc < Num_cyc) {
 
    //undrained tri-axial compression 
    E11 += (- d_g);
    
    thisE.val(1,1) = E11; thisE.val(1,2) = E12; thisE.val(1,3) = E13;
    thisE.val(2,1) = E21; thisE.val(2,2) = E22; thisE.val(2,3) = E23;
    thisE.val(3,1) = E31; thisE.val(3,2) = E32; thisE.val(3,3) = E33;

    FTEP.setTrialStrain(thisE);

    tStress = FTEP.getStressTensor();
    tStrain = FTEP.getStrainTensor();
    
    //ee = FTEP.getVoid();    
    pp  = tStress.p_hydrostatic();
    //qq  = tStress.q_deviatoric();
    qq = - tStress.val(1,1) + tStress.val(2,2);    

    outStress << -tStrain.val(1,1)  <<  "  "  <<  pp << "  " << qq << endl;
    
    if (fabs(E11) >= ep_cut) {
    //if (fabs(qq) >= q_cut) {
        I_cyc++;
        d_g *= (-1.0);
    }    
    

}
    
  return 0;

}
*/



/*
 // Test for Von_Mises
 double rho = 0.0;
 double E = 1.0e4;
 double v = 0.25;
 double k = 4.0;
 double H = 0.0;
 double Hij = 500.0;
 double ha = 5000.0;
 double Cr = 1000.0;
 stresstensor zero;
 double MC[7] = {rho, E, v, H, Hij, ha, Cr};
 double IS[1] = {k};
 stresstensor TS[1] = {zero};
 MaterialParameter matpar(MC,7, IS,1, TS,1);
 Isotropic_Elastic le(2, 3);
 VM_YF dpy(1,1, 2,1);
 VM_PF dpf(2,1);
 Linear_Eeq Eep(4);
 //Linear_Eij Eij(5);
 AF_Eij Eij(6, 7, 1);
 ScalarEvolution *SE = {&Eep};
 TensorEvolution *TE = {&Eij};

 NewTemplate3Dep FTEP(1, &matpar, &le, &dpy, &dpf, &SE, &TE);

 int Num_step = 200;
 double d_g = 0.00005;

 stresstensor thisE;
 stresstensor tStress;
 straintensor tStrain;
 
 double pp = 0.0;
 double qq = 0.0;

 tensor aStress(2, def_dim_2, 0.0);

 for ( int i = 0; i <= Num_step; i++ )
 {

 //Uniaxial Loading
 double E11 = d_g*float(i);    double E12 = 0.0;    double E13 = 0.0;
 double E21 = 0.0;    double E22 = -d_g*float(i);    double E23 = 0.0;
 double E31 = 0.0;    double E32 = 0.0;    double E33 = -d_g*float(i);
 //

 thisE.val(1,1) = E11; thisE.val(1,2) = E12; thisE.val(1,3) = E13;
 thisE.val(2,1) = E21; thisE.val(2,2) = E22; thisE.val(2,3) = E23;
 thisE.val(3,1) = E31; thisE.val(3,2) = E32; thisE.val(3,3) = E33;
 // thisE.print("F","F:");

 FTEP.setTrialStrain(thisE);
 
 tStress = FTEP.getStressTensor();
 tStrain = FTEP.getStrainTensor();
    
 pp  = tStress.p_hydrostatic();
 qq = tStress.val(1,1) - tStress.val(2,2);

 outStress << -tStrain.val(1,1)  <<  "  "  <<  pp << "  " << qq << endl;

 }

 return 0;
 }
*/


/*

 // Test for Cam Clay 

 double rho = 0.0;
 double e0 = 0.8;
 double M = 0.8;
 double lambda = 0.15;
 double kappa = 0.05;
 double v = 0.25;
 double Kc = 200.0;
 double p0 = 200.0;
 stresstensor initStress;
 initStress.val(1,1) = -p0/1.1; initStress.val(2,2) = -p0/1.1; initStress.val(3,3) = -p0/1.1;

 stresstensor zero;
 double MC[7] = {rho, e0, M, lambda, kappa, v, Kc};
 double IS[1] = {p0};
 MaterialParameter inpar(MC,7, IS,1);
 elnp_Elastic le(5, 6, 7, 2, initStress);
 CC_YF dpy(0,3, 1,1);
 CC_PF dpf(0,3, 1,1);
 CC_Ev Ev(4, 5, 2, 1);
 ScalarEvolution *SE = {&Ev};
                                                                            
 NewTemplate3Dep FTEP(1, &inpar, &le, &dpy, &dpf, &SE);

 int Num_step = 1000;
 double d_g = -0.0001;

 stresstensor thisE;
 stresstensor tStress;
 straintensor tStrain;
 
 double pp = 0.0;
 double qq = 0.0;

 tensor aStress(2, def_dim_2, 0.0);

 for ( int i = 0; i <= Num_step; i++ )
 {

 //Uniaxial Loading
 double E11 = d_g*float(i);    double E12 = 0.0;    double E13 = 0.0;
 double E21 = 0.0;    double E22 = -d_g*float(i)*0.5;    double E23 = 0.0;
 double E31 = 0.0;    double E32 = 0.0;    double E33 = -d_g*float(i)*0.5;
 //

 thisE.val(1,1) = E11; thisE.val(1,2) = E12; thisE.val(1,3) = E13;
 thisE.val(2,1) = E21; thisE.val(2,2) = E22; thisE.val(2,3) = E23;
 thisE.val(3,1) = E31; thisE.val(3,2) = E32; thisE.val(3,3) = E33;

 FTEP.setTrialStrain(thisE);
 
 tStress = FTEP.getStressTensor();
 tStrain = FTEP.getStrainTensor();
    
 pp  = tStress.p_hydrostatic();
 qq = -tStress.val(1,1) + tStress.val(2,2);

 outStress << -tStrain.val(1,1)  <<  "  "  <<  pp << "  " << qq << endl;

 }

 return 0;
 }
*/


//  //*********************   Simple Shear ******************************
//  double E11 = -0.001;    double E12 = 0.0;    double E13 = d_g*float(i);
//  double E21 = 0.0;    double E22 = -0.001;    double E23 = 0.0;
//  double E31 = d_g*float(i);    double E32 = 0.0;    double E33 = -0.001;
//  //*******************************************************************

 
// /***********************************************************************************/
// /*   Uniaxial                                                                      */
// double E11 = d_g*float(i);    double E12 = 0.0;    double E13 = 0.0;
// double E21 = 0.0;    double E22 = d_g*float(i)*(-v);    double E23 = 0.0;
// double E31 = 0.0;    double E32 = 0.0;    double E33 = d_g*float(i)*(-v); 

// /***********************************************************************************/
// /*   Triaxial Compression                                                          */
// double E11 = 0.0;    double E12 = 0.0;    double E13 = 0.0;
// double E21 = 0.0;    double E22 = 0.0;    double E23 = 0.0;
// double E31 = 0.0;    double E32 = 0.0;    double E33 = -d_g*float(i);

////*********************   Uniaxial Loading **************************
//double E11 = d_g*float(i);    double E12 = 0.0;    double E13 = 0.0;
//double E21 = 0.0;    double E22 = d_g*float(i)*(-v);    double E23 = 0.0;
//double E31 = 0.0;    double E32 = 0.0;    double E33 = d_g*float(i)*(-v);
////*******************************************************************
