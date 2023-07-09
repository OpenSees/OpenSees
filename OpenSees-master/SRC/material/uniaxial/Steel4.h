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

// $Revision: 1.01 $
// $Date: 2016-03-23 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel4.h,v $

// Written by Adam Zsarnoczay (zsarnoczay@vbt.bme.hu)
// Created: 2013-04-15 10:24:20 $
//
// Description: This file contains the class definition for Steel4
// Steel4 is based on Steel02 and its Menegotto-Pinto Steel Material


#ifndef Steel4_h
#define Steel4_h

#include <UniaxialMaterial.h>
#include <vector>
#include <algorithm>

typedef std::vector<int>    int_vec;
typedef std::vector<double> double_vec;

class Steel4 : public UniaxialMaterial
{
  public:
     Steel4(   int tag,
	       //basics
	       double f_y,   double E_0, 
	       //kinematic hardening
	       double b_k,   double R_0,    double r_1,   double r_2,
	       double b_kc,  double R_0c,   double r_1c,  double r_2c,
	       //isotropic hardening
	       double b_i,   double rho_i,  double b_l,  double R_i,  double l_yp,
	       double b_ic,  double rho_ic, double b_lc, double R_ic,
	       //ultimate strength limit
	       double f_u,   double R_u,    double f_uc,  double R_uc,
	       //load history memory
	       int cycNum,
	       //initial stress
	       double sig_init
	       );
     
     Steel4(void);
     virtual ~Steel4();
     
     const char *getClassType(void) const {return "Steel4";};
     
     double getInitialTangent(void);
     UniaxialMaterial *getCopy(void);
     
     double isoHardening(double mu, double b_i, double b_l, double rho_i, double R_i);
     void loadReversal(int loadDir);
     void calcBreakpoints(int loadDir, double eps_0BC, double sig_0BC, double df_yC, double df_ykC,
			  double eps_pl_tot, double& eps_yC, double& sig_yC, double& eps_uC);
     double calcStress(int loadDir, double eps_C, double eps_0C, double sig_0C, double eps_0BC,
		       double sig_0BC, bool saveProps, double df_yC, double df_ykC);   
     int setTrialStrain(double strain, double strainRate = 0.0);
     
     double getStrain(void);
     double getStress(void);
     double getTangent(void);
     
     int commitState(void);
     int revertToLastCommit(void);
     int revertToStart(void);
     
     int sendSelf(int commitTag, Channel &theChannel);
     int recvSelf(int commitTag, Channel &theChannel,
		  FEM_ObjectBroker &theBroker);
     
     void Print(OPS_Stream &s, int flag =0);
     
 protected:
     
 private:
     // MATERIAL INPUTS -----------------------------------------------------------------------------
     
     // elastic behavior
     double f_y;       //  yield strength
     double E_0;       //  initial stiffness (=Young's modulus)
     
     //kinematic hardening - tension
     double b_k;       //  hardening ratio (E_k/E_0) 
     double R_0;       //  base value for exp transition linear elastic - hardening asymptote 
     double r_1;       //  coefficient for changing R_0 to R_y
     double r_2;       //  coefficient for changing R_0 to R_y
     //                    - compression
     double b_kc;      //  hardening ratio (E_c/E_0)
     double R_0c;      //  base value exp transition linear elastic - hardening asymptote
     double r_1c;      //  coefficient for changing R_0c to R_y
     double r_2c;      //  coefficient for changing R_0c to R_y
     
     //isotropic hardening
     double l_yp;      //  yield plateau length in eps_y0 units
     //                    - tension
     double b_i;       //  initial isotropic hardening ratio
     double R_i;       //  exp transition initial asymptote - saturated asymptote
     double rho_i;     //  isotropic hardening parameter
     double b_l;       //  saturated isotropic hardening ratio
     //                    - compression
     double b_ic;      //  initial isotropic hardening ratio
     double R_ic;      //  exp transition initial asymptote - saturated asymptote
     double rho_ic;    //  isotropic hardening parameter
     double b_lc;      //  saturated isotropic hardening ratio
     
     //ultimate strength limit - tension
     double f_u;       //  ultimate strength
     double R_u;       //  exp transition kinematic hardening - perfectly plastic
     //                        - compression
     double f_uc;      //  ultimate strength
     double R_uc;      //  exp transition kinematic hardening - perfectly plastic
     
     //initial stress
     double sig_init;  //  initial stress value
     
     //size of initial vector for load history memory
     int    cycNum;
     
     
     // FIXED PROPERTIES ----------------------------------------------------------------------------
     double eps_y0;    //  initial yield strain
     double E_t;  	    //  stiffness of the kinematic hardening asymptote under tension
     double E_c;  	    //  stiffness of the kinematic hardening asymptote under compression
     double eps_inc;   //  eps increment used for tangent evaluation = 1e-7
     
     
     // INTERNAL VARIABLES --------------------------------------------------------------------------
     int    dir;       //   loading direction 0-start, 1-tension, 2-compression, 3-prestress
     double eps;       //   eps at current step
     double sig;       //   sig at current step
     double eps_min;   //   minimum strain during the total load history
     double eps_max;   //   maximum strain during the total load history
     double eps_l;     //   eps at the intersection of the hardening and perfectly plastic asymptotes
     double eps_y;     //   eps at intersection of hardening and linear elastic asymptotes
     double sig_y;     //   sig at intersection of hardening and linear elastic asymptotes
     double eps_0;     //   eps at previous load reversal point
     double sig_0;     //   sig at previous load reversal point
     double eps_0B;    //   eps after "unloading" the previous half-cycle
     double sig_0B;    //   sig after "unloading" the previous half-cycle
     double eps_0Y;    //   base yield strain in current half cycle - without in-cycle iso hardening
     double eps_plTot; //   total accumulated plastic strain until the previous load reversal point
     double eps_pl;    //   accumulated plastic strain in the current half-cycle
     double E;         //   stiffness (tangent slope) at the current eps-sig point
     double deps_O;    //   displ. of the origin of the yield surface due to perf. plastic defs.
     double df_yi;     //   yield stress adjustment to consider non-symmetric isotropic hardening
     double df_yk;     //   yield stress adjustment to consider non-symmetric kinematic hardening
     
     
     // STORED VALUES -------------------------------------------------------------------------------
     int    dir_P;  
     double eps_P;       
     double sig_P;       
     double eps_minP;   
     double eps_maxP;   
     double eps_lP;     
     double eps_yP;     
     double sig_yP;     
     double eps_0P;     
     double sig_0P;     
     double eps_0BP;    
     double sig_0BP;    
     double eps_0YP;         
     double eps_plTotP; 
     double eps_plP;    
     double E_P;         
     double deps_OP;    
     double df_yiP;
     double df_ykP;
     
     
     //  VARIABLES FOR LOAD HISTORY MEMORY ----------------------------------------------------------
     int    parentCount;   //  number of half-cycles before the current half-cycle
     int_vec    dir_Par;   //  vector of dir values from previous half-cycles
     double_vec df_yiPar;  //  vector of df_yi values from previous half-cycles
     double_vec df_ykPar;  //  vector of df_yk values from previous half-cycles
     
     //                - tension
     double eps_01;        //  eps_0
     double sig_01;        //  sig_0
     double eps_01B;       //  eps_0B
     double sig_01B;       //  sig_0B
     double_vec eps_01Par; //  vector of eps_0 values from previous half-cycles
     double_vec sig_01Par; //  vector of sig_0 values from previous half-cycles
     double_vec eps_01BPar;//  vector of eps_0B values from previous half-cycles
     double_vec sig_01BPar;//  vector of sig_0B values from previous half-cycles
     
     //                - compression
     double eps_02;        //  eps_0
     double sig_02;        //  sig_0
     double eps_02B;       //  eps_0B
     double sig_02B;       //  sig_0B
     double_vec eps_02Par; //  vector of eps_0 values from previous half-cycles
     double_vec sig_02Par; //  vector of sig_0 values from previous half-cycles
     double_vec eps_02BPar;//  vector of eps_0B values from previous half-cycles
     double_vec sig_02BPar;//  vector of sig_0B values from previous half-cycles


  // TEMPORARY VARIABLES STORED HERE FOR EFFICIENCY ----------------------------------------------
  double deltaEps;  //  delta epsilon during the last load step
  double sig_yD;    //  temp variable for yield stress calculation
  double eps_yD;    //  temp variable for yield strain calculation
  double eps_lD;    //  temp variable for ultimate limit strain calculation
  double eps_plD;   //  temp variable for plastic strain calculation
  double sig_D;     //  temp variable for stress calculation
  double eps_ratY;  //  eps*
  double eps_ratU;  //  eps_bar
  double xi;        //  xi value for R_y calculation
  double R_y;       //  exp transition linear elastic - hardening asymptote 
  double R_uy;      //  exp transition kinematic hardening - perfectly plastic 
  double shft;      //  s for increasing the yield surface size due to isotropic hardening
  double sig_Par;   //  stress in a parent curve at a given strain level
  double sig_inc;   //  sig increment used for tangent evaluation
	
};


#endif

