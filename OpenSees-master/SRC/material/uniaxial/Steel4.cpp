/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2015-02-05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel4.cpp,v $

// Written by Adam Zsarnoczay (zsarnoczay@vbt.bme.hu)
// Created: 2013-04-15 10:24:20 $
//
// Description: This file contains the class implementation for Steel4
// Steel4 is based on Steel02 and its Menegotto-Pinto Steel Material

#include <stdlib.h>
#include <Steel4.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <elementAPI.h>
#include <math.h>

void*
OPS_Steel4(void)
{
  //Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double basicData[2],
    kinData[8],
    isoData[9],
    ultData[4],
    initData[1];
  int    memData[1];
  int argc = 1,
    numBasic = 2,
    numKin = 4,
    numIso = 5,
    numUlt = 2,
    numMem = 1,
    numInit = 1;


  if (OPS_GetIntInput(&argc, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Steel4 tag" << endln;
    return 0;
  }

  //Read the basic parameters
  argc = OPS_GetNumRemainingInputArgs();
  //if the number of arguments is less than the minimum, throw an error
  if (argc < numBasic) {
    opserr << "Invalid #args, want at least two args for Steel4 in the following format:\n" << 
      "uniaxialMaterial Steel4" << iData[0] << " E0? fy?" << endln;
    return 0;
  }
  //if the first two parameters are not doubles, throw an error
  else if (OPS_GetDoubleInput(&numBasic, basicData) != 0) {
    opserr << "Invalid args; E0 and fy for Steel4 (tag: " << iData[0] << 
      ") shall be given as floating point numbers" << endln; 
    return 0;
  }
  //Parsing was successful, the basic data are stored in basicData
  //Load the default values for everything else
  else {
    //kinematic hardening - bilinear perfectly plastic behavior by default
    kinData[0] = 0.0;    //b_k
    kinData[1] = 50.0;   //R_0
    kinData[2] = 0.1;    //r_1
    kinData[3] = 0.15;   //r_2
    for (unsigned int i=0; i<4; i++) { kinData[i+4] = kinData[i]; }
    //isotropic hardening - turned off by default
    isoData[0] = 0.0;  //b_i
    isoData[1] = 1.0;  //rho_i
    isoData[2] = 0.0;  //b_l
    isoData[3] = 50.0; //R_i
    isoData[4] = 0.0;  //l_yp
    for (unsigned int i=0; i<4; i++) { isoData[i+5] = isoData[i]; }
    //ultimate strength limit - turned off by default (i.e. extremely high limit is applied)
    ultData[0] = 100000000.0 * basicData[0];  //f_u
    ultData[1] = 50.0;                        //R_u
    for (unsigned int i=0; i<2; i++) { ultData[i+2] = ultData[i]; }
    //load history memory - turned on by default
    memData[0] = 50;   //cycNum
    //initial stress - zero by default
    initData[0] = 0.0; //sig_0
  }

  //Read the optional string tags to see what else is there to model
  argc = OPS_GetNumRemainingInputArgs();
  //const char *argvLoc = 0;
  while (argc > 1) {
    //char argvLoc[10];
    //if there are no valid string tags throw an error
	const char *argvLoc = OPS_GetString();
  
    //if the material is asymmetric, modify the number of required parameters
    if (strcmp(argvLoc, "-asym") == 0) {
      numKin = 8;
      numIso = 9;
      numUlt = 4;
    } 
    //kinematic hardening
    else if  (strcmp(argvLoc, "-kin") == 0) {
      //check if the right number of parameters are provided
      if (OPS_GetDouble(&numKin, kinData) != 0) {
        opserr << "WARNING invalid -kin args for Steel4 (tag: " << iData[0] << ")\n" << endln;
        return 0;
      }
      //define the unspecified values in case of symmetric hardening
      if (numKin == 4) {
        for (unsigned int i=0; i<4; i++) { kinData[i+4] = kinData[i]; }
      }
    } 
    //isotropic hardening with MP characteristics
    else if  (strcmp(argvLoc, "-iso") == 0) {
      //check if the right number of parameters are provided
      if (OPS_GetDouble(&numIso, isoData) != 0) {
        opserr << "WARNING invalid -iso args for Steel4 (tag: " << iData[0] << ")\n" << endln;
        return 0;
      }
      //define the unspecified values in case of symmetric hardening
      if (numIso == 5) {
        for (unsigned int i=0; i<4; i++) { isoData[i+5] = isoData[i]; }
      }
    }  
    //ultimate strength limit
    else if  (strcmp(argvLoc, "-ult") == 0) {
      //check if the right number of parameters are provided
      if (OPS_GetDouble(&numUlt, ultData) != 0) {
        opserr << "WARNING invalid -ult args for Steel4 (tag: " << iData[0] << ")\n" << endln;
        return 0;
      }
      //define the unspecified values in case of symmetric hardening
      if (numUlt == 2) {
        for (unsigned int i=0; i<2; i++) { ultData[i+2] = ultData[i]; }
      }
    } 
    //load history memory
    else if  (strcmp(argvLoc, "-mem") == 0) {
      //check if the right number of parameters are provided
      if (OPS_GetInt(&numMem, memData) != 0) {
        opserr << "WARNING invalid -mem args for Steel4 (tag: " << iData[0] << ")\n" << endln;
        return 0;
      } 
    } 
    //initial stress
    else if  (strcmp(argvLoc, "-init") == 0) {
      //check if the right number of parameters are provided
      if (OPS_GetDouble(&numInit, initData) != 0) {
        opserr << "WARNING invalid -init args for Steel4 (tag: " << iData[0] << ")\n" << endln;
        return 0;
      }
    }

    argc = OPS_GetNumRemainingInputArgs();
  }

  //Allocate the material
  theMaterial = new Steel4(iData[0], 
                           basicData[0], basicData[1],
                           kinData[0],   kinData[1],   kinData[2],   kinData[3],
                           kinData[4],   kinData[5],   kinData[6],   kinData[7],
                           isoData[0],   isoData[1],   isoData[2],   isoData[3],
                           isoData[4],   isoData[5],   isoData[6],   isoData[7],
                           isoData[8],
                           ultData[0],   ultData[1],   ultData[2],   ultData[3],
                           memData[0],
                           initData[0]
                           );

  //Just in case there was a problem with material creation
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Steel4\n";
    return 0;
  }

  return theMaterial;
}

Steel4::Steel4(  int tag,
              //basics
              double _f_y,   double _E_0, 
              //kinematic hardening
              double _b_k,   double _R_0,    double _r_1,   double _r_2,
              double _b_kc,  double _R_0c,   double _r_1c,  double _r_2c,
              //isotropic hardening
              double _b_i,   double _rho_i,  double _b_l,  double _R_i, double _l_yp,
              double _b_ic,  double _rho_ic, double _b_lc, double _R_ic,
              //ultimate strength limit
              double _f_u,   double _R_u,    double _f_uc,  double _R_uc,
              //load history memory
              int _cycNum,
              //initial stress
              double _sig_init
              ):
    UniaxialMaterial(tag, MAT_TAG_Steel4),
      f_y(_f_y), E_0(_E_0),
      b_k(_b_k), R_0(_R_0), r_1(_r_1), r_2(_r_2), 
      b_kc(_b_kc), R_0c(_R_0c), r_1c(_r_1c), r_2c(_r_2c),
      b_i(_b_i), rho_i(_rho_i), b_l(_b_l), R_i(_R_i), l_yp(_l_yp),
      b_ic(_b_ic), rho_ic(_rho_ic), b_lc(_b_lc), R_ic(_R_ic), 
      f_u(_f_u), R_u(_R_u), f_uc(_f_uc), R_uc(_R_uc),
      cycNum(_cycNum),
      sig_init(_sig_init)
{
  E_t = b_k * E_0;					//kinematic hardening stiffness under tension
  E_c = b_kc * E_0;					//kinematic hardening stiffness under compression
  eps_y0 = f_y / E_0;			  //initial yield strain
  eps_inc = 0.0000001;      //eps increment used for tangent evaluation = 1e-7
  this->revertToStart();
  this->revertToLastCommit();
}

Steel4::Steel4(void):
  UniaxialMaterial(0, MAT_TAG_Steel4)
{
  dir_P = 0;
}

Steel4::~Steel4(void)
{
  // Placeholder, does nothing atm
}

UniaxialMaterial*
Steel4::getCopy(void)
{
  Steel4 *theCopy = new Steel4(this->getTag(),
	  f_y, E_0, 
    b_k, R_0, r_1, r_2,
    b_kc, R_0c, r_1c, r_2c,
	  b_i, rho_i, b_l, R_i, l_yp, 
    b_ic, rho_ic, b_lc, R_ic,
    f_u, R_u, f_uc, R_uc,
    cycNum,
    sig_init);

  //current round
  theCopy-> dir        = dir;      
  theCopy-> eps        = eps;      
  theCopy-> sig        = sig;      
  theCopy-> eps_min    = eps_min;  
  theCopy-> eps_max    = eps_max;  
  theCopy-> eps_l      = eps_l;    
  theCopy-> eps_y      = eps_y;    
  theCopy-> sig_y      = sig_y;    
  theCopy-> eps_0      = eps_0;    
  theCopy-> sig_0      = sig_0;    
  theCopy-> eps_0B     = eps_0B;   
  theCopy-> sig_0B     = sig_0B;   
  theCopy-> eps_0Y     = eps_0Y;   
  theCopy-> eps_plTot  = eps_plTot;
  theCopy-> eps_pl     = eps_pl;   
  theCopy-> E          = E;        
  theCopy-> deps_O     = deps_O;   
  theCopy-> df_yi      = df_yi;    
  theCopy-> df_yk      = df_yk;    

  //previous round
  theCopy-> dir_P      = dir_P;  
  theCopy-> eps_P      = eps_P;     
  theCopy-> sig_P      = sig_P;     
  theCopy-> eps_minP   = eps_minP;  
  theCopy-> eps_maxP   = eps_maxP;  
  theCopy-> eps_lP     = eps_lP;    
  theCopy-> eps_yP     = eps_yP;    
  theCopy-> sig_yP     = sig_yP;    
  theCopy-> eps_0P     = eps_0P;    
  theCopy-> sig_0P     = sig_0P;    
  theCopy-> eps_0BP    = eps_0BP;   
  theCopy-> sig_0BP    = sig_0BP;   
  theCopy-> eps_0YP    = eps_0YP;   
  theCopy-> eps_plTotP = eps_plTotP;
  theCopy-> eps_plP    = eps_plP;   
  theCopy-> E_P        = E_P;       
  theCopy-> deps_OP    = deps_OP;   
  theCopy-> df_yiP     = df_yiP;
  theCopy-> df_ykP     = df_ykP;

  //load history memory
  theCopy-> parentCount = parentCount;
  theCopy-> dir_Par     = dir_Par;
  theCopy-> df_yiPar    = df_yiPar;
  theCopy-> df_ykPar    = df_ykPar;
  theCopy-> eps_01      = eps_01; 
  theCopy-> sig_01      = sig_01; 
  theCopy-> eps_01B     = eps_01B;
  theCopy-> sig_01B     = sig_01B;
  theCopy-> eps_01Par   = eps_01Par;
	theCopy-> sig_01Par   = sig_01Par;
	theCopy-> eps_01BPar  = eps_01BPar;
	theCopy-> sig_01BPar  = sig_01BPar;
  theCopy-> eps_02      = eps_02; 
  theCopy-> sig_02      = sig_02; 
  theCopy-> eps_02B     = eps_02B;
  theCopy-> sig_02B     = sig_02B;
  theCopy-> eps_02Par   = eps_02Par;
	theCopy-> sig_02Par   = sig_02Par;
	theCopy-> eps_02BPar  = eps_02BPar;
	theCopy-> sig_02BPar  = sig_02BPar;
	
  return theCopy;
}

double
Steel4::getInitialTangent(void)
{
  return E_0;
}

double
Steel4::isoHardening(double eps_pl_tot, double b_1, double b_2, double rho, double R)
{
  //calculates the effect of isotropic hardening
  if (eps_pl_tot / eps_y0 > l_yp) {
    double mu = eps_pl_tot / eps_y0 - l_yp;
    if (mu > 10.0*DBL_EPSILON) {
      return mu * (b_2 + (b_1 - b_2) / pow(1.0 + pow(b_1 * mu / rho,R), 1/R) );
    }
  }
  return 0.0;
}

void
Steel4::loadReversal(int loadDir)
{
  //PAY ATTENTION! loadDir here describes the direction of the NEXT cycle!
  //calculate the shift of the origin due to perfectly plastic behavior
  if ((loadDir == 2 && eps_P >= eps_lP) || (loadDir == 1 && eps_P <= eps_lP)) {
    deps_O = deps_O + (eps_P - eps_lP);
  }
  //add the plastic strains from the previous half-cycle to the total stored values
  eps_plTot+=eps_plP;
  eps_pl = 0; 
  //modify the load direction
  dir = loadDir;
  //update the initial sig and eps values
  eps_0 = eps_P;
  sig_0 = sig_P;
  double deltaEpsB = 0.0;
  if (loadDir == 1) {
    if (eps_0>eps_yP) {
      deltaEpsB = (sig_yP-sig_0)/E_0;
    }
    else if (eps_0>eps_lP) {
      deltaEpsB = (sig_yP+(eps_0-eps_yP)*E_0*b_kc-sig_0)/(E_0*(1-b_kc));
    }
    else {
      deltaEpsB = (-f_uc-sig_0)/E_0;
    }
    sig_0B = sig_0 + deltaEpsB*E_0 - (sig_yP - sig_0BP + df_ykP);
    eps_0B = eps_0 + deltaEpsB - (sig_yP - sig_0BP + df_ykP)/E_0;
  }
  else if (loadDir == 2) {
    if (eps_0<eps_yP) {
      deltaEpsB = (sig_yP-sig_0)/E_0;
    }
    else if (eps_0<eps_lP) {
      deltaEpsB = (sig_yP+(eps_0-eps_yP)*E_0*b_k-sig_0)/(E_0*(1-b_k));
    }
    else {
      deltaEpsB = (f_u-sig_0)/E_0;
    }
    sig_0B = sig_0 + deltaEpsB*E_0 - (sig_yP - sig_0BP - df_ykP);
    eps_0B = eps_0 + deltaEpsB - (sig_yP - sig_0BP - df_ykP)/E_0;
  }
  //calculate the influence of asymmetric hardening on the Bauschinger effect
  if (loadDir == 1) {
    df_yi = df_yiP + (isoHardening(eps_plTot,b_ic,b_lc,rho_ic,R_ic) - isoHardening(eps_plTot,b_i,b_l,rho_i,R_i)) * f_y;
    df_yk = E_0 * (eps_0B - sig_0B/E_0 - deps_O) * (b_k - b_kc)/((1-b_k)*(1-b_kc));
  }
  else if (loadDir == 2) {    
    df_yi = df_yiP + (isoHardening(eps_plTot,b_i,b_l,rho_i,R_i) - isoHardening(eps_plTot,b_ic,b_lc,rho_ic,R_ic)) * f_y;
    df_yk = E_0 * (eps_0B - sig_0B/E_0 - deps_O) * (b_k - b_kc)/((1-b_k)*(1-b_kc));
  }
  //update the extreme strains if they have been exceeded
  if (eps_P < eps_min)  eps_min = eps_P;
  if (eps_P > eps_max)  eps_max = eps_P;
  //store the curve origin for parent curve database
  if (loadDir == 1) {
    eps_01  = eps_0;   sig_01  = sig_0;
    eps_01B = eps_0B;  sig_01B = sig_0B;
  }
  else if (loadDir == 2) {
    eps_02  = eps_0;	  sig_02  = sig_0;
    eps_02B = eps_0B;   sig_02B = sig_0B;
  }
}

void
Steel4::calcBreakpoints(int loadDir, double eps_0BC, double sig_0BC, double df_yC, double df_ykC,
                        double eps_pl_tot, double& eps_yC, double& sig_yC, double& eps_lC)
{
  //calculate the total accumulated plastic strain dependent isotropic hardening
  shft = 1.0;
  if (loadDir == 1)      { shft += isoHardening(eps_pl_tot,b_i,b_l,rho_i,R_i);     }
  else if (loadDir == 2) { shft += isoHardening(eps_pl_tot,b_ic,b_lc,rho_ic,R_ic); }
  //calculate the intersection point for elastic and hardening behavior
  sig_D = 0.0;
  if (loadDir == 1) {                             //tension
    sig_D =  shft*f_y + df_ykC + df_yC;
    eps_yC = eps_0BC + sig_D/E_0;
    sig_yC = sig_0BC + sig_D;
  }
  else if (loadDir == 2) {                        //compression
    sig_D =  shft*f_y + df_ykC + df_yC;
    eps_yC = eps_0BC - sig_D/E_0;
    sig_yC = sig_0BC - sig_D;
  }
  //sig_y cannot exceed f_u
  eps_lC = DBL_MAX * (3 - 2 * loadDir) / 2.0;
  if (loadDir == 1) {                              //tension
    if (sig_yC > f_u) {
        eps_yC -= (sig_yC - f_u) / E_0;
        sig_yC = f_u;
    }
    if (E_t > 10.0*DBL_EPSILON) {
      eps_lC =  (f_u - sig_yC) / E_t + eps_yC;
    }
  }
  else if (loadDir == 2) {                        //compression
    if (sig_yC < -f_uc) {
      eps_yC -= (sig_yC + f_uc) / E_0;
      sig_yC = - f_uc;
    }
    if (E_c > 10.0*DBL_EPSILON) {
      eps_lC = (-f_uc - sig_yC) / E_c + eps_yC;        
    }
  }
}

double
Steel4::calcStress(int loadDir, double eps_C, double eps_0C, double sig_0C, double eps_0BC, 
                   double sig_0BC, bool saveProps, double df_yC, double df_ykC)
{
  //calculate the plastic strain at eps_C strain
  if      (dir==1)  eps_plD = std::max(   eps_C-eps_0Y ,0.0);                  
  else if (dir==2)  eps_plD = std::max(- (eps_C-eps_0Y),0.0);
  //calculate the intersection point for elastic and hardening behavior
  calcBreakpoints(dir,eps_0BC,sig_0BC,df_yC,df_ykC,eps_plTot + eps_plD,eps_yD,sig_yD,eps_lD);
  //calculate the value of eps*
  eps_ratY  = (eps_C-eps_0C)/(eps_yD-eps_0C);
  //calculate the value of eps_bar
  eps_ratU = (eps_C-eps_0C)/(eps_lD-eps_0C);
  //calculate the radius for the gradual shift from elastic to hardening to perf. plastic behavior
  R_y = 50.0;
  R_uy = 50.0;
  if (loadDir == 1) {                                //tension
    xi = fabs((eps_max-eps_0Y)/eps_y0);
    R_y = R_0*(1.0 - (r_1*xi)/(r_2+xi));
    R_uy = R_u;
  }
  else if (loadDir == 2) {                           //compression
    xi = fabs((eps_0Y-eps_min)/eps_y0);
    R_y = R_0c*(1.0 - (r_1c*xi)/(r_2c+xi));
    R_uy = R_uc;
  }
  //calculate the stress corresponding to eps_C strain
  sig_D = 0.0;
  if (loadDir == 1) {                              //tension
    sig_D = (sig_yD - sig_0C) * eps_ratY *
            ( b_k        / pow( 1.0 + pow(fabs(eps_ratU), R_uy), (1/R_uy)) +
             (1.0 - b_k) / pow( 1.0 + pow(fabs(eps_ratY), R_y), (1/R_y))  ) +
            sig_0C;
  }
  else if (loadDir == 2) {                        //compression
    sig_D = (sig_yD - sig_0C) * eps_ratY *
			      ( b_kc        / pow( 1.0 + pow(fabs(eps_ratU), R_uy), (1/R_uy)) +
		    	   (1.0 - b_kc) / pow( 1.0 + pow(fabs(eps_ratY), R_y), (1/R_y))  ) +
			      sig_0C;
  }
  //save the calculated properties if needed
  if (saveProps) {
    eps_y = eps_yD;
    sig_y = sig_yD;
    eps_l = eps_lD;
    eps_pl = eps_plD;
  }
  return sig_D;
}

int
Steel4::setTrialStrain(double trialStrain, double strainRate)
{
  //include the effect of initial stress if specified
  if (sig_init != 0.0) eps = trialStrain + sig_init / E_0;  
  else                 eps = trialStrain;            // source: Steel02 mod. by C-P. Lamarche 2006
  deltaEps = eps - eps_P;		//delta epsilon
  //Load the core variables
  eps_max   = eps_maxP;    eps_min   = eps_minP;
  eps_y     = eps_yP;      sig_y     = sig_yP;
  eps_0     = eps_0P;      sig_0     = sig_0P;
  eps_0B    = eps_0BP;     sig_0B    = sig_0BP;
  df_yi     = df_yiP;      df_yk     = df_ykP;
  dir       = dir_P;       eps_plTot = eps_plTotP;
  eps_l     = eps_lP;      eps_0Y    = eps_0YP;
  deps_O    = deps_OP;
  //set the initial values of the controlling variables
  if (dir == 0 || dir == 3) {  // = initial phase
    //if delta epsilon is very small, there is nothing to calculate
    if (fabs(deltaEps) < 10.0*DBL_EPSILON) {
      E = E_0;
      sig = sig_init;
      dir = 3;
      return 0;
    }
    else {  // = delta epsilon is large enough to make a decision on the load direction
      if (deltaEps < 0.0) {    // negative direction -> compression
        dir = 2;
        eps_y   = -eps_y0;   sig_y   = -f_y;
        eps_0B  = 0;         sig_0B  = 0;
        eps_0Y  = -eps_y0;   eps_min = eps_0Y;              
      }
      else {                   // positive direction -> tension
        dir = 1;
        eps_y  = eps_y0;    sig_y   = f_y;
        eps_0B = 0;         sig_0B  = 0;
        eps_0Y = eps_y0;    eps_max = eps_0Y;
      }
    }
  }   
  //in case of load reversal, adjust the corresponding internal variables
  if      (dir == 2 && deltaEps > 0.0) {
    loadReversal(1);     //negative -> positive load reversal
    calcBreakpoints(dir,eps_0B,sig_0B,df_yi,df_yk,eps_plTot,eps_0Y,sig_yD,eps_lD);
  }
  else if (dir == 1 && deltaEps < 0.0) {
    loadReversal(2);     //positive -> negative load reversal
    calcBreakpoints(dir,eps_0B,sig_0B,df_yi,df_yk,eps_plTot,eps_0Y,sig_yD,eps_lD);
  } 
  //calculate the stress corresponding to eps on the actual loading curve
  sig = calcStress(dir,eps,eps_0,sig_0,eps_0B,sig_0B,true,df_yi,df_yk);
  //calculate the stress increase over eps_inc for stiffness assessment
  sig_inc = sig - calcStress(dir,eps - eps_inc,eps_0,sig_0,eps_0B,sig_0B,false,df_yi,df_yk);
  E = sig_inc / eps_inc;
  //check all parent curves and if one of them is exceeded, merge with it
  if (cycNum != 0 && dir == 1) {                                          //tension
	  int i = parentCount;
	  sig_Par = 0.0;                       
	  while (i>0) {
		  if (dir_Par[i]==1) {			
			  if ((eps_01Par[i]-(eps_0-(sig_0-sig_01Par[i])/E_0)) < - 10*DBL_EPSILON) { 
          //calculate the stress corresponding to eps on the parent loading curve
          sig_Par = calcStress(dir,eps,eps_01Par[i],sig_01Par[i],eps_01BPar[i],sig_01BPar[i],false,
                               df_yiPar[i],df_ykPar[i]);
				  if (sig_Par<sig) {            
					  eps_0  = eps_01Par[i]; 		  sig_0  = sig_01Par[i];
            eps_0B = eps_01BPar[i]; 	  sig_0B = sig_01BPar[i];
            df_yi   = df_yiPar[i];      df_yk  = df_ykPar[i];   
            //calculate the stress corresponding to eps on the merged loading curve
            sig = calcStress(dir,eps,eps_0,sig_0,eps_0B,sig_0B,true,df_yi,df_yk);            
            //calculate the stress increase over eps_inc for stiffness assessment
            sig_inc = sig - calcStress(dir,eps - eps_inc,eps_0,sig_0,eps_0B,sig_0B,true,
                                       df_yi,df_yk); 
            E = sig_inc / eps_inc;            
            parentCount = i;
				  }
			  }
			  i = 0;
		  }
		  i --;
	  }
  }
  else if (cycNum != 0 && dir == 2) {                                       //compression
	  int i = parentCount;
	  sig_Par = 0.0;
	  while (i>0) {
		  if (dir_Par[i]==2) {
			  if (((eps_0+(sig_02Par[i]-sig_0)/E_0)-eps_02Par[i]) < - 10 * DBL_EPSILON) {
          //calculate the stress corresponding to eps on the parent loading curve
          sig_Par = calcStress(dir,eps,eps_02Par[i],sig_02Par[i],eps_02BPar[i],sig_02BPar[i],false,
                               df_yiPar[i],df_ykPar[i]);         
				  if (sig_Par>sig) {
					  eps_0  = eps_02Par[i];		  sig_0  = sig_02Par[i];
            eps_0B = eps_02BPar[i];		  sig_0B = sig_02BPar[i];
            df_yi   = df_yiPar[i];      df_yk  = df_ykPar[i];					   
            //calculate the stress corresponding to eps on the parent loading curve
            sig = calcStress(dir,eps,eps_0,sig_0,eps_0B,sig_0B,true,df_yi,df_yk);
            //calculate the stress increase over eps_inc for stiffness assessment          
            sig_inc = sig - calcStress(dir,eps - eps_inc,eps_0,sig_0,eps_0B,sig_0B,true,
                                       df_yi,df_yk);
            E = sig_inc / eps_inc;
            parentCount = i;
				  }
			  }
			  i=0;
		  }
		  i --;
	  }
  }
  return 0;
}

double
Steel4::getStrain(void)
{
  return eps;
}

double
Steel4::getStress(void)
{
  return sig;
}

double
Steel4::getTangent(void)
{
  return E;
}

int
Steel4::commitState(void)
{
  //store the internal variables
  dir_P      = dir;      
  eps_P      = eps;      
  sig_P      = sig;      
  eps_minP   = eps_min;  
  eps_maxP   = eps_max;  
  eps_lP     = eps_l;    
  eps_yP     = eps_y;    
  sig_yP     = sig_y;    
  eps_0P     = eps_0;    
  sig_0P     = sig_0;    
  eps_0BP    = eps_0B;   
  sig_0BP    = sig_0B;   
  eps_0YP    = eps_0Y;   
  eps_plTotP = eps_plTot;
  eps_plP    = eps_pl;   
  E_P        = E;        
  deps_OP    = deps_O;   
  df_yiP     = df_yi;    
  df_ykP     = df_yk;    

  //if load history memory is turned on
  if (cycNum != 0) {
    //and the loading direction changed
    if (fabs(eps_01-eps_01Par[parentCount]) > DBL_EPSILON || 
	      fabs(eps_02-eps_02Par[parentCount]) > DBL_EPSILON || 
	      dir!=dir_Par[parentCount]) {
      //and the previous half-cycle went to inelastic range
      if ( fabs(eps - eps_0) > 2 * f_y / E_0 &&
           ((dir == 1 && eps_0 == eps_01) || (dir == 2 && eps_0 == eps_02))) {
        //then store the data that describes the previous half-cycle
	      parentCount ++;        
        //resize vectors if necessary
        int currentSize = eps_01Par.size();
        if (parentCount >= currentSize) {  
          dir_Par.resize(currentSize+cycNum);
          df_yiPar.resize(currentSize+cycNum);
          df_ykPar.resize(currentSize+cycNum);
          eps_01Par.resize(currentSize+cycNum);
	  sig_01Par.resize(currentSize+cycNum);
	  eps_01BPar.resize(currentSize+cycNum);
	  sig_01BPar.resize(currentSize+cycNum);
          eps_02Par.resize(currentSize+cycNum);
	  sig_02Par.resize(currentSize+cycNum);          
	  eps_02BPar.resize(currentSize+cycNum);
	  sig_02BPar.resize(currentSize+cycNum);                  
        }
        dir_Par[parentCount]   = dir;
	df_yiPar[parentCount]  = df_yi;
        df_ykPar[parentCount]  = df_yk;
        eps_01Par[parentCount] = eps_01;
	sig_01Par[parentCount] = sig_01;
	eps_01BPar[parentCount]= eps_01B;
	sig_01BPar[parentCount]= sig_01B;
        eps_02Par[parentCount] = eps_02;
	sig_02Par[parentCount] = sig_02;         
	eps_02BPar[parentCount]= eps_02B;
	sig_02BPar[parentCount]= sig_02B;               
      }
    }
  }
  return 0;
}

int
Steel4::revertToLastCommit(void)
{
  //load the stored values
  dir       = dir_P;       
  eps       = eps_P;       
  sig       = sig_P;       
  eps_min   = eps_minP;    
  eps_max   = eps_maxP;    
  eps_l     = eps_lP;      
  eps_y     = eps_yP;      
  sig_y     = sig_yP;      
  eps_0     = eps_0P;      
  sig_0     = sig_0P;      
  eps_0B    = eps_0BP;     
  sig_0B    = sig_0BP;     
  eps_0Y    = eps_0YP;     
  eps_plTot = eps_plTotP;  
  eps_pl    = eps_plP;     
  E         = E_P;         
  deps_O    = deps_OP;     
  df_yi     = df_yiP;      
  df_yk     = df_ykP;      

  return 0;
}

int
Steel4::revertToStart(void)
{
  //initialize the stored values
  dir_P      = 0;
  eps_P      = 0.0;
  sig_P      = 0.0;
  eps_minP   = 0.0;
  eps_maxP   = 0.0;
  eps_lP     = 0.0;
  eps_yP     = 0.0;
  sig_yP     = 0.0;
  eps_0P     = 0.0;
  sig_0P     = 0.0;
  eps_0BP    = 0.0;
  sig_0BP    = 0.0;
  eps_0YP    = 0.0;
  eps_plTotP = 0.0;
  eps_plP    = 0.0;
  E_P        = E_0;
  deps_OP    = 0.0;
  df_yiP     = 0.0;
  df_ykP     = 0.0;

  if (sig_init != 0.0) {    //based on Steel02
	  eps_P = sig_init/E_0;   //it might be unnecessary
	  sig_P = sig_init;
   }

  parentCount = 0;
  sig_01 = 0.0;
  eps_01 = 0.0;
  sig_01B = 0.0;
  eps_01B = 0.0;
  eps_02 = 0.0;
  sig_02 = 0.0;
  eps_02B = 0.0;
  sig_02B = 0.0;

  //initialize the vectors that collect load history memory data
  unsigned int initSize = 2;
  if (cycNum != 0) {
    dir_Par   .reserve(cycNum + initSize);
    df_yiPar  .reserve(cycNum + initSize);
    df_ykPar  .reserve(cycNum + initSize);
    eps_01Par .reserve(cycNum + initSize);
    sig_01Par .reserve(cycNum + initSize);
    eps_01BPar.reserve(cycNum + initSize);
    sig_01BPar.reserve(cycNum + initSize);
    eps_02Par .reserve(cycNum + initSize);
    sig_02Par .reserve(cycNum + initSize);    
    eps_02BPar.reserve(cycNum + initSize);
    sig_02BPar.reserve(cycNum + initSize);
    for (unsigned int i=0; i<initSize; i++) {
      dir_Par   .resize(cycNum + initSize);
      df_yiPar  .resize(cycNum + initSize);
      df_ykPar  .resize(cycNum + initSize);
      eps_01Par .resize(cycNum + initSize);
      sig_01Par .resize(cycNum + initSize);
      eps_01BPar.resize(cycNum + initSize);
      sig_01BPar.resize(cycNum + initSize);
      eps_02Par .resize(cycNum + initSize);
      sig_02Par .resize(cycNum + initSize);      
      eps_02BPar.resize(cycNum + initSize);
      sig_02BPar.resize(cycNum + initSize);      
    }
  }

  return 0;
}

int
Steel4::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
Steel4::recvSelf(int commitTag, Channel &theChannel,
	     FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
Steel4::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Steel4 tag: " << this->getTag() << endln;
        s << "  fy: " << f_y << ", ";
        s << "  E0: " << E_0 << ", ";
        s << "  b_kt: " << b_k << ", ";
        s << "  R0_kt: " << R_0 << ", ";
        s << "  cR1_kt: " << r_1 << ", ";
        s << "  cR2_kt: " << r_2 << ", ";
        s << "  b_kc: " << b_kc << ", ";
        s << "  R0_kc: " << R_0c << ", ";
        s << "  cR1_kc: " << r_1c << ", ";
        s << "  cR2_kc: " << r_2c << ", ";
        s << "  lyp: " << l_yp << ", ";
        s << "  b_it: " << b_i << ", ";
        s << "  R_it: " << R_i << ", ";
        s << "  rho_it: " << rho_i << ", ";
        s << "  bl_it: " << b_l << ", ";
        s << "  b_ic: " << b_ic << ", ";
        s << "  R_ic: " << R_ic << ", ";
        s << "  rho_ic: " << rho_ic << ", ";
        s << "  bl_ic: " << b_lc << ", ";
        s << "  fu_t: " << f_u << ", ";
        s << "  Ru_t: " << R_u << ", ";
        s << "  fu_c: " << f_uc << ", ";
        s << "  Ru_c: " << R_uc << ", ";
        s << "  sigini: " << sig_init << ", ";
        s << "  cycNum: " << cycNum;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"Steel4\", ";
        s << "\"E\": " << E_0 << ", ";
        s << "\"fy\": " << f_y << ", ";
        s << "\"b_kt\": " << b_k << ", ";
        s << "\"R0_kt\": " << R_0 << ", ";
        s << "\"cR1_kt\": " << r_1 << ", ";
        s << "\"cR2_kt\": " << r_2 << ", ";
        s << "\"b_kc\": " << b_kc << ", ";
        s << "\"R0_kc\": " << R_0c << ", ";
        s << "\"cR1_kc\": " << r_1c << ", ";
        s << "\"cR2_kc\": " << r_2c << ", ";
        s << "\"lyp\": " << l_yp << ", ";
        s << "\"b_it\": " << b_i << ", ";
        s << "\"R_it\": " << R_i << ", ";
        s << "\"rho_it\": " << rho_i << ", ";
        s << "\"bl_it\": " << b_l << ", ";
        s << "\"b_ic\": " << b_ic << ", ";
        s << "\"R_ic\": " << R_ic << ", ";
        s << "\"rho_ic\": " << rho_ic << ", ";
        s << "\"bl_ic\": " << b_lc << ", ";
        s << "\"fu_t\": " << f_u << ", ";
        s << "\"Ru_t\": " << R_u << ", ";
        s << "\"fu_c\": " << f_uc << ", ";
        s << "\"Ru_c\": " << R_uc << ", ";
        s << "\"sigini\": " << sig_init << ", ";
        s << "\"cycNum\": " << cycNum << "}";
    }
}
