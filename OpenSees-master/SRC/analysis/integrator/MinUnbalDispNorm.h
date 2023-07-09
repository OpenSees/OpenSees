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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:00:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/MinUnbalDispNorm.h,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/MinUnbalDispNorm.h
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the class definition for MinUnbalDispNorm.
// MinUnbalDispNorm is an algorithmic class for performing a static analysis
// using the minimum unbalanced displacement norm (Chan IJNME 26(2657:2669)1988
//
// What: "@(#) MinUnbalDispNorm.h, revA"

#ifndef MinUnbalDispNorm_h
#define MinUnbalDispNorm_h

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class Domain;
#define SIGN_LAST_STEP      1
#define CHANGE_DETERMINANT  2

class MinUnbalDispNorm : public StaticIntegrator
{
  public:
    MinUnbalDispNorm(double lambda1, int specNumIterStep, 
		     double dlambda1min, double dlambda1max,
		     int signFirstStepMethod = SIGN_LAST_STEP);

    ~MinUnbalDispNorm();

    int newStep(void);    
    int update(const Vector &deltaU);
    int domainChanged(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
      //////////////////Sensitivity Begin//////////////////////////////////
      int formEleResidual(FE_Element *theEle);
      int formSensitivityRHS(int gradNum);// it's been modified to compute dLambdadh and dUdh
      int formIndependentSensitivityRHS();
      int saveSensitivity(const Vector &v, int gradNum, int numGrads);
      int saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads);
      int commitSensitivity(int gradNum, int numGrads);
      int computeSensitivities(void);// this function is modified to obtain both dLambdadh and dUdh 

      /////////////////////////////// Abbas //////////////////////////////
      Vector *formTangDispSensitivity(Vector *dUhatdh,int gradNumber); //Obtain *dUhatdh 
     // Vector *formResidualDispSensitivity(Vector *dUIJdh,int gradNumber);// Obtain dKdh*deltaUbar
      double formdLambdaDh(int gradNumber);//calculate dLambdadh for J=1
      double getLambdaSensitivity(int gradNumber);// update the dLambdadh for J>1
      bool computeSensitivityAtEachIteration();// A key that return 1 for loadControl and 2 for DisplacementControl
     // int newStepSens(int gradIndex);
      ////////////////////Sensitivity End/////////////////////////////////////


  protected:
     double dlambdadh; // deltaLambdaI1 for the first iteration J=1
      double Dlambdadh;// deltaLambdaIJ: for J>1
      double dLambda;
     double calldLambda1dh;//Abbas
      int CallParam;

  private:
    double dLambda1LastStep;                  // dLambda1 at step (i-1)
    double specNumIncrStep, numIncrLastStep;    // Jd & J(i-1) 
   double dLambdaj; //for J>1
 
    Vector *deltaUhat, *deltaUbar, *deltaU, *deltaUstep; // vectors for disp measures
    Vector *phat; 	                                 // the reference load vector

    double deltaLambdaStep, currentLambda; // dLambda(i) & current value of lambda  
    int signLastDeltaLambdaStep;           // sign of dLambda(i-1)
    double dLambda1min, dLambda1max;       // min & max values for dlambda1 at step (i) 
    double signLastDeterminant;
    int signFirstStepMethod;

 // int theDofID, theDof; 
      ///////////////////////////////////////Abbas/////////////////////////////////////////
      // Pointers used for sensitivity analysis
      Vector  *dUhatdh,*dUIJdh, *Residual,*Residual2, *sensU,*d_deltaU_dh, *dphatdh, *dLAMBDAdh ;
      // the created pointers shown above are
      // *dUhatdh     : The derivative of the tangent displacement w/r to parameter h
      // *sensU       : Displacement sensitivity using displacement control scheme
      // *d_deltaU_dh : The derivative of the residual displacement
      // *dUIJdh      : The sensitivity of the residual displacement
      // *Residual    : the residual forces that are required to obtain the dLambdadh
      // *Residual2   : the residual forces required to obtain dUdh


      // the reference load vector
        
      double dLambdaStepDh ;//Abbas
      double minIncrement, maxIncrement; // min/max values of deltaU at (i)

      // adding sensitivity
      int gradNumber;
      int sensitivityFlag;
      FE_Element *theEle;

};

#endif

