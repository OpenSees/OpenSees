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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestNormUnbalance.cpp,v $
                                                                        
                                                                        

#include <CTestNormUnbalance.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>

CTestNormUnbalance::CTestNormUnbalance()	    	
:ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
 theSOE(0), tol(0.0), maxNumIter(0), currentIter(0), printFlag(0)
{
    
}

CTestNormUnbalance::CTestNormUnbalance(double theTol, int maxIter, int printIt)
:ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
 theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0), printFlag(printIt)
{

}

CTestNormUnbalance::~CTestNormUnbalance()
{
    
}

void
CTestNormUnbalance::setTolerance(double newTol)
{
    tol = newTol;
}

int
CTestNormUnbalance::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    theSOE = theAlgo.getLinearSOEptr();
    if (theSOE == 0) {
	cerr << "WARNING: CTestNormUnbalance::setEquiSolnAlgo - no SOE\n";	
	return -1;
    }
    else
	return 0;
}

    
int
CTestNormUnbalance::test(void)
{
    const Vector &x = theSOE->getB();
    double norm = x.Norm();

    // print the data if required
    if (printFlag == 1) {
      cerr << "\t CTestNormUnbalance::test() - iteration: " << currentIter;
      cerr << " current Norm: " << norm << " (max permissable: " << tol << ")\n";
    }

    if (norm <= tol) { // the algorithm converged
      if (printFlag != 0) {
	if (printFlag == 1) 
	  cerr << endl;
	else if (printFlag == 2) {
	  cerr << "\t CTestNormUnbalance::test() - iteration: " << currentIter;
	  cerr << " last Norm: " << norm << " (max permissable: " << tol << ")\n";
	}
      }
      return currentIter;
    } 

    else if (currentIter >= maxNumIter) { // the algorithm failed to converge
      cerr << "WARNING: CTestNormUnbalance::test() - failed to converge \n";
      cerr << "after: " << currentIter << " iterations\n";	
      currentIter++;  // we increment in case analysis does not check for convergence
      return -2;
    } 

    else { // the algorithm has not converged yet
      currentIter++;    
      return -1;
    }
}


int
CTestNormUnbalance::start(void)
{
    currentIter = 1;
    if (theSOE == 0) {
	cerr << "WARNING: CTestNormUnbalance::test() - no SOE returning true\n";
	return -1;
    }
    
    return 0;
}

int 
CTestNormUnbalance::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  Vector x(2);
  x(0) = tol;
  x(1) = maxNumIter;  
  res = theChannel.sendVector(this->getDbTag(), cTag, x);
  if (res < 0) 
    cerr << "CTestNormUnbalance::sendSelf() - failed to send data\n";
    
  return res;
}

int 
CTestNormUnbalance::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  Vector x(2);
  res = theChannel.recvVector(this->getDbTag(), cTag, x);    

  if (res < 0) {
      tol = 1.0e-8;
      maxNumIter = 25;
      cerr << "CTestNormUnbalance::sendSelf() - failed to send data\n";
  }
  else {
      tol = x(0);
      maxNumIter = x(1);
  }
  return res;
}

