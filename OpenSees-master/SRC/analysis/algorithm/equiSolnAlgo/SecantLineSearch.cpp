
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
// $Date: 2003-04-02 22:02:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/SecantLineSearch.cpp,v $

// Written: fmk 
// Created: 11/01
// 
// What: "@(#)SecantLineSearch.h, revA"

#include <SecantLineSearch.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <math.h>

SecantLineSearch::SecantLineSearch(double tol, int mIter, double mnEta, double mxEta, int pFlag)
:LineSearch(LINESEARCH_TAGS_SecantLineSearch),
 x(0), tolerance(tol), maxIter(mIter), minEta(mnEta), maxEta(mxEta), printFlag(pFlag)
{   

}

SecantLineSearch::~SecantLineSearch()
{
  if (x != 0)
    delete x;
}


int 
SecantLineSearch::newStep(LinearSOE &theSOE)
{
  const Vector &dU = theSOE.getX();

  if (x == 0)
    x = new Vector(dU);

  if (x->Size() != dU.Size()) {
    delete x;
    x = new Vector(dU);
  }

  return 0;
}

int 
SecantLineSearch::search(double s0, 
				 double s1, 
				 LinearSOE &theSOE, 
				 IncrementalIntegrator &theIntegrator)
{
  double r0 = 0.0;

  if ( s0 != 0.0 ) 
    r0 = fabs( s1 / s0 );
	
  if  (r0 <= tolerance )
    return 0; // Line Search Not Required Residual Decrease Less Than Tolerance

  if (s1 == s0)
    return 0;  // Secant will have a divide-by-zero if continue

  // set some variables
  double eta    = 1.0;
  double s      = s1;
  double etaJ   = 1.0;
  double etaJm1 = 0.0;
  double sJ     = s1;
  double sJm1   = s0;
  double r = r0;

  const Vector &dU = theSOE.getX();

  if (printFlag == 0) {
    opserr << "Secant Line Search - initial: "
	 << "      eta(0) : " << eta << " , Ratio |s/s0| = " << r0 << endln;
  }

  // perform the secant iterations:
  //
  //                eta(j+1) = eta(j) -  s(j) * (eta(j-1)-eta(j))
  //                                     ------------------------
  //                                           s(j-1) - s(j)

  int count = 0; //initial value of iteration counter 
  while ( r > tolerance  &&  count < maxIter ) {
    
    count++;

    eta = etaJ - sJ * (etaJm1-etaJ) / (sJm1 - sJ);

    //-- want to put limits on eta and stop solution blowing up
    if (eta > maxEta)  eta = maxEta;
    if (r   > r0    )  eta =  1.0;
    if (eta < minEta)  eta = minEta;
    
    //update the incremental difference in response and determine new unbalance
    if (eta == etaJ) 
      break; // no change in response

    *x = dU;
    *x *= eta-etaJ;
	    
    if (theIntegrator.update(*x) < 0) {
      opserr << "WARNING SecantLineSearch::search() -";
      opserr << "the Integrator failed in update()\n";	
      return -1;
    }
    
    if (theIntegrator.formUnbalance() < 0) {
      opserr << "WARNING SecantLineSearch::search() -";
      opserr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }	

    //new residual
    const Vector &ResidJ = theSOE.getB();
    
    //new value of s
    s = dU ^ ResidJ;
    
    //new value of r 
    r = fabs( s / s0 ); 

    if (printFlag == 0) {
      opserr << "Secant Line Search - iteration: " << count 
	   << " , eta(j) : " << eta << " , Ratio |sj/s0| = " << r << endln;
    }

    if (etaJ == eta)
      count = maxIter;

    // set variables for next iteration
    etaJm1 = etaJ;
    etaJ = eta;
    sJm1 = sJ;
    sJ = s;

    if (sJm1 == sJ)
      count = maxIter;
    
  } //end while

  // set X in the SOE for the revised dU, needed for convergence tests
  *x = dU;
  if (eta != 0.0)
    *x *= eta;
  theSOE.setX(*x);

  return 0;
}


int
SecantLineSearch::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
SecantLineSearch::recvSelf(int cTag, 
			   Channel &theChannel, 
			   FEM_ObjectBroker &theBroker)
{
  return 0;
}


void
SecantLineSearch::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) {
    s << "SecantLineSearch :: Line Search Tolerance = " << tolerance << endln; 
    s << "                       max num Iterations = " << maxIter << endln;
    s << "                         max value on eta = " << maxEta << endln;
  }
}









