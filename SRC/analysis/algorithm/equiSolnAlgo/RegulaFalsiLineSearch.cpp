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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/RegulaFalsiLineSearch.cpp,v $

// Written: fmk 
// Created: 11/01
// 
// What: "@(#)RegulaFalsiLineSearch.h, revA"

#include <RegulaFalsiLineSearch.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <math.h>

RegulaFalsiLineSearch::RegulaFalsiLineSearch(double tol, int mIter, double mnEta, double mxEta, int pFlag)
:LineSearch(LINESEARCH_TAGS_RegulaFalsiLineSearch),
 x(0), tolerance(tol), maxIter(mIter), minEta(mnEta), maxEta(mxEta), printFlag(pFlag)
{   

}

RegulaFalsiLineSearch::~RegulaFalsiLineSearch()
{
  if (x != 0)
    delete x;
}


int 
RegulaFalsiLineSearch::newStep(LinearSOE &theSOE)
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
RegulaFalsiLineSearch::search(double s0, 
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
    return 0;  // RegulaFalsi will have a divide-by-zero error if continue

  // set some variables
  double eta    = 1.0;
  double s      = s1;
  double etaU   = 1.0;
  double etaL   = 0.0;
  double sU     = s1;
  double sL     = s0;
  double r      = r0;
  double etaJ   = 1.0;
  double compoundFactor = 0.0;

  const Vector &dU = theSOE.getX();

  if (printFlag == 0) {
    opserr << "RegulaFalsi Line Search - initial: "
	 << "      eta(0) : " << eta << " , Ratio |s/s0| = " << r0 << endln;
  }

  // we first search for a bracket to a solution, i.e. we want sU * sL < 0.0
  int count = 0;
  while ((sU * sL > 0.0) && (etaU < maxEta)) {

    count++;

    /*
    if (count == 1)
      etaU = 0.5;
    else
    */
    etaU = etaJ * 4.0;

    //update the incremental difference in response and determine new unbalance
    *x = dU;
    double factor = etaU - etaJ;
    compoundFactor += factor;
    *x *= factor;

    etaJ = etaU;

    if (theIntegrator.update(*x) < 0) {
      opserr << "WARNING BisectionLineSearch::search() -";
      opserr << "the Integrator failed in update()\n";	
      return -1;
    }
    
    if (theIntegrator.formUnbalance() < 0) {
      opserr << "WARNING BisectionLineSearch::search() -";
      opserr << "the Integrator failed in formUnbalance()\n";	
      return -2;
    }	
  
    //new residual
    const Vector &ResidJ = theSOE.getB();
    
    //new value of sU
    sU = dU ^ ResidJ;

    // check if we have a solution we are happy with
    r = fabs( sU / s0 ); 
    if (r < tolerance)
      return 0;

    if (printFlag == 0) {
      opserr << "Bisection Line Search - bracketing: " << count 
	   << " , eta(j) : " << etaU << " , Ratio |sj/s0| = " << r << endln;
    }
  }

  // return if no bracket for a solution found, resetting to initial values
  if (sU * sL > 0.0) {
    *x = dU;
    theSOE.setX(*x);
    *x *= -compoundFactor;
    theIntegrator.update(*x);
    theIntegrator.formUnbalance();
    return 0; 
  }

  // perform the secant iterations:
  //
  //                eta(j+1) = eta(u) -  s(u) * (eta(l) -eta(u))
  //                                     ------------------------
  //                                           s(l) - s(u)

  count = 0; //initial value of iteration counter 
  while ( r > tolerance  &&  count < maxIter ) {
    
    count++;

    eta = etaU - sU * (etaL-etaU) / (sL - sU);


    //-- want to put limits on eta(i)
    if (eta > maxEta)  eta = maxEta;
    if (  r >  r0   )  eta =  1.0;
    if (eta < minEta)  eta = minEta;

    if (eta == etaJ) // break if going to have a zero *x
      break;
    
    //update the incremental difference in response and determine new unbalance
    *x = dU;
    *x *= eta-etaJ;
	    
    if (theIntegrator.update(*x) < 0) {
      opserr << "WARNING RegulaFalsiLineSearch::search() -";
      opserr << "the Integrator failed in update()\n";	
      return -1;
    }
    
    if (theIntegrator.formUnbalance() < 0) {
      opserr << "WARNING RegulaFalsiLineSearch::search() -";
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
      opserr << "RegulaFalsi Line Search - iteration: " << count 
	   << " , eta(j) : " << eta << " , Ratio |sj/s0| = " << r << endln;
    }
    

    if (etaJ == eta)
      count = maxIter;

    // set variables for next iteration
    etaJ = eta;
    
    if (s*sU < 0.0) {
      etaL = eta;
      sL   = s;
    } else if (s*sU == 0.0)
      count = maxIter;
    else {
      etaU = eta;
      sU   = s;
    } 

    if (sL == sU)
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
RegulaFalsiLineSearch::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
RegulaFalsiLineSearch::recvSelf(int cTag, 
				Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
  return 0;
}


void
RegulaFalsiLineSearch::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) {
    s << "RegulaFalsiLineSearch :: Line Search Tolerance = " << tolerance << endln; 
    s << "                         max num Iterations = " << maxIter << endln;
    s << "                         max value on eta = " << maxEta << endln;
  }
}









