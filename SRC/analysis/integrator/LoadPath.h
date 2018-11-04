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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/LoadPath.h,v $
                                                                        
                                                                        
#ifndef LoadPath_h
#define LoadPath_h

// File: ~/analysis/integrator/LoadPath.h
// 
// Written: fmk 
// Created: 05/99
// Revision: A
//
// Description: This file contains the class definition for LoadPath.
// LoadPath is an algorithmic class for performing a static analysis
// using a user defined load path (a user specified lambda path)
//
// What: "@(#) LoadPath.h, revA"

#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;

class LoadPath : public StaticIntegrator
{
  public:
    LoadPath(Vector &theLoadPath);
    LoadPath();    

    ~LoadPath();

    int newStep(void);    
    int update(const Vector &deltaU);

    // Public methods for Output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
    Vector *loadPath;
    int currentStep;
};

#endif

