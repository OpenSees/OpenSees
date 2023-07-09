
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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/InitialInterpolatedLineSearch.cpp,v $

// Written: fmk 
// Created: 11/01
// 
// What: "@(#)InitialInterpolatedLineSearch.h, revA"

#include <InitialInterpolatedLineSearch.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <math.h>


InitialInterpolatedLineSearch::InitialInterpolatedLineSearch(double tol, int mIter, double mnEta,
							     double mxEta, int pFlag)
:LineSearch(LINESEARCH_TAGS_InitialInterpolatedLineSearch),
 x(0), tolerance(tol), maxIter(mIter), minEta(mnEta), maxEta(mxEta), printFlag(pFlag)
{   

}

InitialInterpolatedLineSearch::~InitialInterpolatedLineSearch()
{
  if (x != 0)
    delete x;
}

int 
InitialInterpolatedLineSearch::newStep(LinearSOE &theSOE)
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
InitialInterpolatedLineSearch::search(double s0, 
				      double s1, 
				      LinearSOE &theSOE, 
				      IncrementalIntegrator &theIntegrator)
{
  double s = s1;

  //initialize r = ratio of residuals 
  double r0 = 0.0;

  if ( s0 != 0.0 ) 
    r0 = fabs( s / s0 );
	
  if  (r0 <= tolerance )
    return 0; // Line Search Not Required Residual Decrease Less Than Tolerance

  double eta = 1.0;     //initial value of line search parameter
  double etaPrev = 1.0;
  double r = r0;

  const Vector &dU = theSOE.getX();

  int count = 0; //initial value of iteration counter 

  if (printFlag == 0) {
    opserr << "InitialInterpolated Line Search - initial       "
	 << "    eta : " << eta 
	 << " , Ratio |s/s0| = " << r0 << endln;
  }    


  // Solution procedure follows the one in Crissfields book.
  // (M.A. Crissfield, Nonlinear Finite Element Analysis of Solid and Structures, Wiley. 97).
  // NOTE: it is not quite linear interpolation/false-position/regula-falsi as eta(0) = 0.0
  // does not change. uses eta(i) = eta(i-1)*s0
  //                                -----------
  //                                s0 - s(i-1)  to compute eta(i)


  //  opserr << "dU: " << dU;

  while ( r > tolerance  &&  count < maxIter ) {

    count++;
    
    eta *=  s0 / (s0 - s); 

    //-- want to put limits on eta(i)
    if (eta > maxEta)  eta = maxEta;
    if (r   > r0    )  eta =  1.0;
    if (eta < minEta)  eta = minEta;

    if (eta == etaPrev)
      break; // no change in response break

    //dx = ( eta * dx0 ); 
    *x = dU;
    *x *= eta-etaPrev;
	    
    if (theIntegrator.update(*x) < 0) {
      opserr << "WARNInG InitialInterpolatedLineSearch::search() -";
      opserr << "the Integrator failed in update()\n";	
      return -1;
    }
    
    if (theIntegrator.formUnbalance() < 0) {
      opserr << "WARNInG InitialInterpolatedLineSearch::search() -";
      opserr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }	

    //new residual
    const Vector &ResidI = theSOE.getB();
    
    //new value of s
    s = dU ^ ResidI;

    //new value of r 
    r = fabs( s / s0 ); 

    if (printFlag == 0) {
      opserr << "InitialInterpolated Line Search - iteration: " << count 
	   << " , eta(j) : " << eta
	   << " , Ratio |sj/s0| = " << r << endln;
    }    

    // reset the variables, also check not just hitting bounds over and over
    if (eta == etaPrev)
      count = maxIter;
    else
      etaPrev = eta;

  } //end while

  // set X in the SOE for the revised dU, needed for convergence tests
  *x = dU;

  if (eta != 0.0)
    *x *= eta;

  theSOE.setX(*x);

  return 0;
}


int
InitialInterpolatedLineSearch::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
InitialInterpolatedLineSearch::recvSelf(int cTag, 
					Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  return 0;
}


void
InitialInterpolatedLineSearch::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) 
    s << "InitialInterpolatedLineSearch :: Line Search Tolerance = " << tolerance << endln; 
}









