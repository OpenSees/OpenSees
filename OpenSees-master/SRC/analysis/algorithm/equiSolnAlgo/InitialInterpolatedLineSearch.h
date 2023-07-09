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
// $Date: 2003-02-14 23:00:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/InitialInterpolatedLineSearch.h,v $

// Written: fmk 
// Created: 11/01

// Description: This file contains the class definition for LinearInterpolatedSearch.
// This performs the search by using a form of linear interpolation to find the best solution.
// Solution procedure follows the one in Crissfields book.
// (M.A. Crissfield, Nonlinear Finite Element Analysis of Solid and Structures, Wiley. 97).
// NOTE: it is not quite linear interpolation/false-position/regula-falsi as eta(0) = 0.0
// does not change. uses eta(i) = eta(i-1)*s0
//                                -----------
//                                s0 - s(i-1)  to compute eta(i)
//

#ifndef InitialInterpolatedLineSearch_h
#define InitialInterpolatedLineSearch_h

#include <LineSearch.h>
class Vector;
//class OPS_Stream; //Jeremic@ucdavis.edu taken out since there is an include<iOPS_Stream.h> in LineSearch.h 

class InitialInterpolatedLineSearch: public LineSearch
{
  public:
  InitialInterpolatedLineSearch(double tolerance = 0.8, 
				int    maxIter   = 10,
				double minEta    = 0.1,
				double maxEta    = 10.0, 
				int    printFlag = 1);

    ~InitialInterpolatedLineSearch();

    int newStep(LinearSOE &theSOE);
    int search(double s0, 
	       double s1, 
	       LinearSOE &theSOE, 
	       IncrementalIntegrator &theIntegrator);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0) ;    
    
  protected:
    
  private:
    Vector *x;
    double tolerance;
    int    maxIter;
    double minEta;
    double maxEta;
    int    printFlag;
};

#endif


