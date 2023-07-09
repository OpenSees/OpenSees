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
// $Date: 2003-02-14 23:02:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/cg/ConjugateGradientSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/cg/ConjugateGradientSolver.C
//
// Written: fmk 
// Created: 06/00
// Revision: A

#include "ConjugateGradientSolver.h"
#include <Vector.h>
#include <OPS_Globals.h>
#include <LinearSOE.h>

ConjugateGradientSolver::ConjugateGradientSolver(int classtag, 
						 LinearSOE *theSOE,
						 double tol)
:LinearSOESolver(classtag),
 r(0),p(0),Ap(0),x(0), 
 theLinearSOE(theSOE), 
 tolerance(tol)
{
    
}

ConjugateGradientSolver::~ConjugateGradientSolver()
{
    if (r != 0)
	delete r;
    if (p != 0)
	delete p;
    if (Ap != 0)
	delete Ap;
    if (x != 0)
	delete x;    
}


int 
ConjugateGradientSolver::setSize(void)
{
    int n = theLinearSOE->getNumEqn();
    if (n <= 0) {
	opserr << "ConjugateGradientSolver::setSize() - n < 0 \n";
	return -1;
    }

    // if old Vector exist and are of incorrect size destroy them
    if (r != 0) {
	if (r->Size() != n) {
	    delete r;
	    delete p;
	    delete Ap;	
	    delete x;		    
	    r = 0;
	    p = 0;
	    Ap = 0;
	    x = 0;	    
	}
    }

    // create new vector of correct size
    if (r == 0) {
	r = new Vector(n);
	p = new Vector(n);
	Ap = new Vector(n);
	x = new Vector(n);	
	if (r == 0 || p == 0 || Ap == 0 || x == 0) {
	    opserr << "ConjugateGradientSolver::setSize() - out of memory\n";
	    if (r != 0)
		delete r;
	    if (p != 0)
		delete p;
	    if (Ap != 0)
		delete Ap;
	    if (x != 0)
		delete x;    	    
	    r = 0;
	    p = 0;
	    Ap = 0;
	    x = 0;	    	    
	    return -2;
	    
	}
    }
    return 0;
}



int
ConjugateGradientSolver::solve(void)
{
    // check for successful setSize
    if (r == 0)
	return -1;
    
    // initialize
    x->Zero();    
    *r = theLinearSOE->getB();
    *p = *r;
    double rdotr = *r ^ *r;
    
    // lopp till convergence
    while (r->Norm() > tolerance) {
	this->formAp(*p, *Ap);

	double alpha = rdotr/(*p ^ *Ap);

	// *x += *p * alpha;
	x->addVector(1.0, *p, alpha);

	// *r -= *Ap * alpha;
	r->addVector(1.0, *Ap, -alpha);

	double oldrdotr = rdotr;

	rdotr = *r ^ *r;

	double beta = rdotr / oldrdotr;

	// *p = *r + *p * beta;
	p->addVector(beta, *r, 1.0);
    }
    return 0;
}











