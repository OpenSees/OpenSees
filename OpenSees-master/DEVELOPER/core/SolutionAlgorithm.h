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
// $Date: 2004-11-13 08:08:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/SolutionAlgorithm.h,v $
                                                                        
                                                                        
#ifndef SolutionAlgorithm_h
#define SolutionAlgorithm_h

// File: ~/OOP/analysis/algorithm/SolutionAlgorithm.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for SolutionAlgorithm.
// SolutionAlgorithm is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// What: "@(#) SolutionAlgorithm.h, revA"

#include <MovableObject.h>

class Channel;
class FEM_ObjectBroker;
class Recorder;

extern int SOLUTION_ALGORITHM_tangentFlag;

class SolutionAlgorithm: public MovableObject
{
  public:
    SolutionAlgorithm(int classTag);
    virtual ~SolutionAlgorithm();

    virtual int domainChanged(void);
    
    // methods for monitoring the analysis during an algorithm
    virtual int  addRecorder(Recorder &theRecorder);    	
    virtual int  record(int track);    
    
  protected:
    
  private:
    Recorder **theRecorders;
    int numRecorders;    
};

#endif


