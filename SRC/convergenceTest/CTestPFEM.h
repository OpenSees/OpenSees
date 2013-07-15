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

// $Revision: 1.0 $
// $Date: 2013-5-23 11:49:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestPFEM.h,v $

// Written: Minjie
// Date: 05/23

// A CTestPFEM object tests for convergence using the ratio of the 
// current norm to the initial norm (the norm when start is invoked) of the 
// velocity and pressure for fluid nodes and velocity for structural nodes


#ifndef CTestPFEM_h
#define CTestPFEM_h

#include <ConvergenceTest.h>
#include <vector>
class EquiSolnAlgo;
class PFEMLinSOE;


class CTestPFEM: public ConvergenceTest
{
public:
    // constructors
    CTestPFEM();	    	
    CTestPFEM(double tv, double tp, double tv2, double tp2, double tvrel, double tprel,
              int maxNumIter, int maxincr, 
              int printFlag, int normType);

    // destructor
    virtual ~CTestPFEM();
    
    ConvergenceTest *getCopy(int iterations);
    
    void setTolerance(double newTol);
    int setEquiSolnAlgo(EquiSolnAlgo &theAlgo);
    
    int test(void);
    int start(void);
    
    int getNumTests(void);
    int getMaxNumTests(void);        
    double getRatioNumToMax(void);                
    const Vector &getNorms(void);    
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
protected:
    
private:
    PFEMLinSOE *theSOE;
    double tolv;         // the velocity tol on the norm used to test for convergence
    double tolp;         // the pressure tol on the norm used to test for convergence
    double tolv2;        // the velocity tol on the residual used to test for convergence
    double tolp2;        // the pressure tol on the residual used to test for convergence
    double tolvrel;
    double tolprel;
    
    int maxNumIter;     // max number of iterations
    int currentIter;    // number of times test() has been invokes since last start()
    int printFlag;      // a flag indicating if to print on test
    int nType;          // type of norm to use (1-norm, 2-norm, p-norm, max-norm)
    int maxIncr;        // max number of norm increasing
    int numIncr;        // number of norm increasing
    
    std::vector<double> normsv;       // vector to hold the velocity norms
    std::vector<double> normsp;       // vector to hold the pressure norms
    std::vector<double> normsresv;    // velocity norm of residual
    std::vector<double> normsresp;    // pressure norm of residual
    double normv0;
    double normp0;
    Vector norms;
};

#endif
