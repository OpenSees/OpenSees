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
                                                                        
// $Revision: 1.4 $
// $Date: 2010-01-20 21:05:33 $
                                                                        
// File: ~/Bond_SP01.h
//
// Written: 		Jian Zhao, Iowa State University 		04/2004
// Revision:		Jian Zhao, University of Wisconsin, Milwaukee 		04/2006
//
// Description: This file contains the class definition for Uniaxial material Bond_SP01. 
//		Bond_SP01: Strain penetration of rebars at footings w/o bond damage.
/* *********************************************************************** **
** Nonlinear material model with Menegoto-Pinto type functions to define   **
** envelope curve and loading-unloading-reloading curve.  test data of     **
** rebar stress vs. rebar slip is needed to calibrate the model.           **
** *********************************************************************** */

#ifndef Bond_SP01_h
#define Bond_SP01_h

#include <UniaxialMaterial.h>
#include <math.h>
#include <Matrix.h>
#include <Vector.h>

class Bond_SP01 : public UniaxialMaterial
{
  public:
    Bond_SP01(int tag, double fy, double sy, double fu, double su, double Kz, double R, double Cd, double db, double fc, double la); 
    Bond_SP01(int tag, double fy, double sy, double fu, double su, double Kz, double R); 
    Bond_SP01();
    ~Bond_SP01();

    const char *getClassType(void) const {return "Bond_SP01";};  

    int setTrialStrain (double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
    /*** anchorage condition ***/
    double db;				// rebar diameter
    double fc;				// concrete compressive strength (positive ksi)
    double lba;				// bond length (can be shorter than current set -> rigid-body slip)
    double la;				// required minimum anchorage length
    
    /*** strain penetration ***/
    double sy;				// rebar slip at bar yielding
    double su;				// rebar slip at bar failure (assume @ 1.5fy)
    double fy;				// bar yield strength
    double fu;				// bar ultimate strength (assume @ 1.5fy)
    double E0;				// initial slope of the envelope
    double Kz;				// initial hardening ratio (0.3), 0.25<Kz<0.5
    double Cr;				// R factor for the envelope using menegoto-pinto curve
    double Ks;				// reloading stop point (to avoid infinity slope with pinching)
    double slvrg;			// slip corresponding to virgin friction in local bond-slip model 
    
    /*** pinching condition ***/
    double R;				// R factor for reloading using menegoto-pinto type curve
    double Cd;				// bond damage factor (currently no use)
    
    /*** CONVERGED History Variables ***/
    double CRSlip;			// return slip
    double CRLoad;			// return load
    double CRSlope;			// return slope
    double CmaxHSlip;		// Maximum slip in tension
    double CminHSlip;		// Maximum slip in compression
    int Cloading;			// Flag for Stressing/unStressing
							// 1 = loading (positive slip increment)
							// -1 = unloaading (negative slip increment)
							// 0 = initial state
    int CYieldFlag;			// Flag for yielding
							// 1 = yielded
							// 0 = not yielded
    
    /*** CONVERGED State Variables ***/    
    double Cslip;
    double Cload;
    double Ctangent;    

    /*** TRIAL History Variables ***/
    double TRSlip;			// return slip
    double TRLoad;			// return load
    double TRSlope;			// return slope
    double TmaxHSlip;		// Maximum slip in tension
    double TminHSlip;		// Maximum slip in compression
    int Tloading;			// Flag for Stressing/unStressing
    int TYieldFlag;			// Flag for yielding
    
    /*** TRIAL State Variables ***/
    double Tslip;
    double Tload;
    double Ttangent;    
    
    // Calculates the trial state variables based on the trial Strain
    void determineTrialState (double tslip, double dslip);
    
    // Calculates envelope stress for a slip
    double getEnvelopeStress (double slip);
    
    // Determines if a Stress reversal has occurred based on the trial Strain
    void detectStressReversal (double dslip);
};
#endif
