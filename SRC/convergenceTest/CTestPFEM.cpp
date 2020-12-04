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
// $Source: /usr/local/cvs/OpenSees/SRC/convergenceTest/CTestPFEM.cpp,v $

#include <CTestPFEM.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include "sparseGEN/PFEMLinSOE.h"
#include <typeinfo>
#include <cmath>
#include <elementAPI.h>
#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif

void* OPS_CTestPFEM()
{
    if(OPS_GetNumRemainingInputArgs() < 7) {
	opserr<<"insufficient number of arguments\n";
	return 0;
    }

    // tolerance
    double tol[6];
    int numData = 6;
    if(OPS_GetDoubleInput(&numData,&tol[0]) < 0) {
	opserr << "WARNING PFEM test failed to get tolerance\n";
	return 0;
    }

    // maxIter
    numData = OPS_GetNumRemainingInputArgs();
    if(numData > 4) numData = 4;
    int data[4] = {20,3,0,2};
    if(OPS_GetIntInput(&numData,&data[0]) < 0) {
	opserr << "WARNING PFEM test failed to get int values\n";
	return 0;
    }

    // create test
    return new CTestPFEM(tol[0],tol[1],tol[2],tol[3],tol[4],tol[5],
			 data[0],data[1],data[2],data[3]);
}

CTestPFEM::CTestPFEM()	    	
    : ConvergenceTest(CONVERGENCE_TEST_CTestPFEM),
      theSOE(0), tolv(0), tolp(0), tolv2(0), tolp2(0), tolvrel(0), tolprel(0),
      maxNumIter(0), currentIter(0), printFlag(0), 
      nType(2), maxIncr(-1), numIncr(0), 
      normsv(), normsp(), normsresv(), normsresp(),
      normv0(0), normp0(0), norms()
{
    
}


CTestPFEM::CTestPFEM(double tv, double tp, double tv2, double tp2, double tvrel, double tprel,
                     int maxIter, int maxincr, int printIt, int normType)
    : ConvergenceTest(CONVERGENCE_TEST_CTestPFEM),
      theSOE(0), tolv(tv), tolp(tp), tolv2(tv2), tolp2(tp2), tolvrel(tvrel), tolprel(tprel),
      maxNumIter(maxIter), currentIter(0), printFlag(printIt),
      nType(normType), maxIncr(maxincr), numIncr(0), 
      normsv(), normsp(), normsresv(), normsresp(),
      normv0(0), normp0(0),norms()
{

}


CTestPFEM::~CTestPFEM()
{
    
}


ConvergenceTest* CTestPFEM::getCopy(int iterations)
{
    CTestPFEM *theCopy ;
    theCopy = new CTestPFEM(this->tolv, this->tolp, this->tolv2, this->tolp2, 
                            this->tolvrel, this->tolprel,
                            iterations, this->maxIncr,
                            this->printFlag, this->nType) ;
    
    theCopy->theSOE = this->theSOE ;
    
    return theCopy ;
    
}


void CTestPFEM::setTolerance(double newTol)
{

}


int CTestPFEM::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
    try
    {
        theSOE = dynamic_cast<PFEMLinSOE*>(theAlgo.getLinearSOEptr());
    }
    catch (std::bad_cast& bc)
    {
        opserr << "WARNING: system is not PFEMLinSOE -- CTestPFEM::setEquiSolnAlgo\n";
        return -1;
    }
    if(theSOE == 0) {
        opserr << "WARNING: CTestPFEM::setEquiSolnAlgo() - no SOE\n";	
        return -1;
    }
    else
        return 0;
}


int CTestPFEM::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if(theSOE == 0)
        return -2;
    
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if(currentIter == 0) {
        opserr << "WARNING: CTestPFEM::test() - start() was never invoked.\n";	
        return -2;
    }
    
    // get the X vector & determine it's norm & save the value in norms vector
    const Vector &x = theSOE->getX();
    const Vector &B = theSOE->getB();
    const ID& dofType = theSOE->getDofType();
    int stage = theSOE->getStage();
    if(dofType.Size() != x.Size()) {
        opserr << "WARNING: x and dofType have different size -- CTestPFEM::test()\n";
        return -2;
    }
    Vector velocity(x.Size()), pressure(x.Size()), pi(x.Size()), 
        resv(x.Size()), resp(x.Size()), respi(x.Size());
    for(int i=0; i<x.Size(); i++) {
        if(dofType(i)==0 || dofType(i)==1 || dofType(i)==2) {
            velocity(i) = x(i);
            resv(i) = B(i);
        } else if(dofType(i) == 3) {
            pressure(i) = x(i);
            resp(i) = B(i);
        } else if(dofType(i) == 4) {
            pi(i) = x(i);
            respi(i) = B(i);
        }
    }
    double normv = velocity.pNorm(nType);
    double normp = pressure.pNorm(nType);
    double normpi = pi.pNorm(nType);
    double normresv = resv.pNorm(nType);
    double normresp = resp.pNorm(nType);
    double normrespi = respi.pNorm(nType);

#ifdef _PARALLEL_INTERPRETERS

    // copy norms from host to all processors
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double allnorms[6];
    if(myid == 0) {
	allnorms[0] = normv;
	allnorms[1] = normp;
	allnorms[2] = normpi;
	allnorms[3] = normresv;
	allnorms[4] = normresp;
	allnorms[5] = normrespi;
    }

    if(MPI_Bcast(&allnorms,6,MPI_DOUBLE,0,MPI_COMM_WORLD) != MPI_SUCCESS) {
	opserr<<"WARNING: failed to copy norms to all processors\n";
	return -1;
    }

    if(myid != 0) {
	normv = allnorms[0];
	normp = allnorms[1];
	normpi = allnorms[2];
	normresv = allnorms[3];
	normresp = allnorms[4];
	normrespi = allnorms[5];
    }

    // no print except host
    if(myid != 0) printFlag = 0;
#endif

    // norm at first step
    if(currentIter == 1) {
        normv0 = normv;
        normp0 = normp;
    } 

    double normvrel = normv/normv0;
    double normprel = normp/normp0;
    if (normv0 == 0) {
        normvrel = 0.0;
    }
    if (normp0 == 0) {
        normprel = 0.0;
    }

    if (stage==1 || stage==3) {
        normp = 0;
        normresp = 0;
        normprel = 0;
        normpi = 0;
        normrespi = 0;
    } else if (stage == 2) {
        normv = 0;
        normresv = 0;
        normvrel = 0;
    }

    // record norms
    if(currentIter <= maxNumIter) {
        normsv.push_back(normv);
        normsp.push_back(normp);
        normsresv.push_back(normresv);
        normsresp.push_back(normresp);
    }
    
    // check for norm increase
    if(stage==0 && currentIter > 1 && maxIncr > 0) {
        if(normv>10*normsv[currentIter-2] || normp>10*normsp[currentIter-2] ||
           normresv>10*normsresv[currentIter-2] || normresp>10*normsresp[currentIter-2]) {
            numIncr++;
        }
    }

    // print the data if required
    if(printFlag == 1) {
        opserr << "PFEM: " << currentIter;
        if (stage == 1) {
            opserr << " -- Predictor Stage\n";
        } else if (stage == 2) {
            opserr << " -- Pressure Stage\n";
        } else if (stage == 3) {
            opserr << " -- Corrector Stage\n";
        }
        opserr << " dV(" << normv << "," << normvrel;
        opserr << "), dP(" << normp << "," << normprel;
        opserr << "), dPi(" << normpi;
        opserr << "), resV(" << normresv;
        opserr << "), resP(" << normresp;
        opserr << "), resPi(" << normrespi;
        opserr << "), incr(" << numIncr<<")\n";
    } 
    if(printFlag == 4) {
    } 
    
    //
    // check if the algorithm converged
    //
    
    // if converged - print & return ok
    bool badnorm = normv!=normv || normp!=normp;
    //if((normv<=tolv||normvrel<=tolv2) && normresv && 
    //(normp<=tolp||normresp<=tolp2||normprel<=tolprel)) {
    if((normv<=tolv && normp<=tolp) || (normresv<=tolv2 && normresp<=tolp2) ||
       (normvrel<=tolvrel && normprel<=tolprel)) { 
        
        // do some printing first
        if(printFlag != 0) {
            if(printFlag == 1 || printFlag == 4) 
                opserr << endln;
            else if(printFlag == 2 || printFlag == 6) {
                opserr << "PFEM: " << currentIter;
                opserr << " dV(" << normv << "," << normvrel;
                opserr << "), dP(" << normp << "," << normprel;
                opserr << "), dPi(" << normpi;
                opserr << "), resV(" << normresv;
                opserr << "), resP(" << normresp;
                opserr << "), resPi(" << normrespi;
                opserr << "), incr(" << numIncr<<")\n";
            }
        }

        if (stage == 1) {
            theSOE->setStage(2);
        } else if (stage == 0 || stage == 3) {

            // return the number of times test has been called
            return currentIter;
        }
    }
    
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if((printFlag == 5 || printFlag == 6) 
            && (currentIter >= maxNumIter || numIncr>maxIncr || badnorm)) {
        opserr << "WARNING: CTestPFEM - failed to converge but going on -";
        opserr << " dV(" << normv << "," << normvrel;
        opserr << "), dP(" << normp << "," << normprel;
        opserr << "), dPi(" << normpi;
        opserr << "), resV(" << normresv;
        opserr << "), resP(" << normresp;
        opserr << "), resPi(" << normrespi;
        opserr << "), incr(" << numIncr<<")\n";
        return currentIter;
    }
    
    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if(currentIter >= maxNumIter || numIncr > maxIncr || badnorm) { // failes to converge
        opserr << "WARNING: CTestPFEM - failed to converge \n";
        opserr << "after: " << currentIter << " iterations:";
        opserr << " dV(" << normv << "," << normvrel;
        opserr << "), dP(" << normp << "," << normprel;
        opserr << "), dPi(" << normpi;
        opserr << "), resV(" << normresv;
        opserr << "), resP(" << normresp;
        opserr << "), resPi(" << normrespi;
        opserr << "), incr(" << numIncr<<")\n";
        currentIter++;    
        return -2;
    } 
    
    // algorithm not yet converged - increment counter and return -1
    if (stage == 2) {
        theSOE->setStage(3);
    }
    currentIter++;
    return -1;
}


int CTestPFEM::start(void)
{
    if(theSOE == 0) {
        opserr << "WARNING: CTestPFEM::test() - no SOE returning true\n";
        return -1;
    }
    
    // set iteration count = 1
    normsv.clear();
    normsp.clear();
    normsresv.clear();
    normsresp.clear();
    currentIter = 1;
    numIncr = 0;
    
    return 0;
}


int CTestPFEM::getNumTests()
{
    return currentIter;
}


int CTestPFEM::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestPFEM::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestPFEM::getNorms() 
{
    int size = 0;
    size += normsv.size();
    size += normsp.size();
    size += normsresv.size();
    size += normsresp.size();

    if(size == 0) {
        norms = Vector();
        return norms;
    }
    
    int loc = 0;
    for(int i=0; i<(int)normsv.size(); i++) {
        norms(loc++) = normsv[i];
    }
    for(int i=0; i<(int)normsp.size(); i++) {
        norms(loc++) = normsp[i];
    }
    for(int i=0; i<(int)normsresv.size(); i++) {
        norms(loc++) = normsresv[i];
    }
    for(int i=0; i<(int)normsresp.size(); i++) {
        norms(loc++) = normsresp[i];
    }

    return norms;
}


int CTestPFEM::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    Vector x(10);
    x(0) = tolv;
    x(1) = maxNumIter;
    x(2) = printFlag;
    x(3) = nType;
    x(4) = maxIncr;
    x(5) = tolp;
    x(6) = tolv2;
    x(7) = tolp2;
    x(8) = tolvrel;
    x(9) = tolprel;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if(res < 0) 
        opserr << "CTestPFEM::sendSelf() - failed to send data\n";
    
    return res;
}


int CTestPFEM::recvSelf(int cTag, Channel &theChannel, 
                        FEM_ObjectBroker &theBroker)
{
    int res = 0;
    Vector x(10);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);    
    
    if(res < 0) {
        opserr << "CTestPFEM::sendSelf() - failed to send data\n";
        tolv = 1.0e-8;
        tolp = 1.0e-8;
        maxNumIter = 25;
        printFlag = 0;
        nType = 2;
        maxIncr = 3;
    }
    else {
        tolv = x(0);
        maxNumIter = (int) x(1);
        printFlag = (int) x(2);
        nType = (int) x(3);
        maxIncr = (int) x(4);
        normsv.clear();
        normsp.clear();
        normsresv.clear();
        normsresp.clear();
        tolp = x(5);
        tolv2 = x(6);
        tolp2 = x(7);
        tolvrel = x(8);
        tolprel = x(9);
    }
    return res;
}
