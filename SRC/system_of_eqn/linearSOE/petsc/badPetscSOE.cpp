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
// $Date: 2003-02-14 23:02:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/badPetscSOE.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/petsc/PetscSOE.C
//
// Written: fmk & om
// Created: 7/98
// Revision: A
//
// Description: This file contains the implementation for BandGenLinSOE


#include "PetscSOE.h"
#include "PetscSolver.h"
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

PetscSOE::PetscSOE(PetscSolver &theSOESolver)
:LinearSOE(theSOESolver, LinSOE_TAGS_PetscSOE),
 isFactored(0),size(0), numProcessors(0), B(0), X(0), 
 vectX(0), vectB(0), A(0), x(0), b(0)
{
    theSOESolver.setLinearSOE(*this);
    MPI_Comm_size(PETSC_COMM_WORLD, &numProcessors);
}


int
PetscSOE::getNumEqn(void) const
{
    return size;
}
    
PetscSOE::~PetscSOE()
{
  /*
  if (vectX != 0) delete vectX;  
  if (vectB != 0) delete vectB;  
  if (B != 0) delete [] B;
  if (X != 0) delete [] X;

  // invoke the petsc destructors
  if (A != 0) MatDestroy(A);
  if (b != 0) VecDestroy(b);
  if (x != 0) VecDestroy(x);
  */
}


int 
PetscSOE::setSize(Graph &theGraph)
{
    int result = 0;
    int oldSize = size;
    size = theGraph.getNumVertex();
    isFactored = 0;

    // delete the old	
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;


    // invoke the petsc destructors
    if (A != 0) MatDestroy(A);
    if (b != 0) VecDestroy(b);
    if (x != 0) VecDestroy(x);

    // create the new
    B = new double[size];
    X = new double[size];
	
    if (B == 0 || X == 0) {
      opserr << "WARNING PetscSOE::PetscSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      size = 0; 
      result = -1;
    }

    // zero the vectors
    for (int j=0; j<size; j++) {
	B[j] = 0;
	X[j] = 0;
    }

    if (vectX != 0) 
      delete vectX;

    if (vectB != 0) 
      delete vectB;
		
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);

    // invoke the petsc constructors
    int ierr = OptionsGetInt(PETSC_NULL,"-n",&size,&flg); CHKERRA(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,size,size,&A); CHKERRA(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,PETSC_DECIDE,size,&x); CHKERRA(ierr); 
    ierr = VecDuplicate(x,&b); CHKERRA(ierr);  
    // invoke setSize() on the Solver
    int solverOK = theSolver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:PetscSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    

    return result;    
}


int 
PetscSOE::setSizeDirectly(int newSize)
{
    int result = 0;
    int oldSize = size;
    size = newSize;
    isFactored = 0;

    // delete the old	
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;

    // invoke the petsc destructors
    if (A != 0) MatDestroy(A);
    if (b != 0) VecDestroy(b);
    if (x != 0) VecDestroy(x);

    // create the new
    B = new double[size];
    X = new double[size];
	
    if (B == 0 || X == 0) {
      opserr << "WARNING PetscSOE::PetscSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      size = 0; 
      result = -1;
    }

    // zero the vectors
    for (int j=0; j<size; j++) {
	B[j] = 0;
	X[j] = 0;
    }

    if (vectX != 0) 
      delete vectX;

    if (vectB != 0) 
      delete vectB;
		
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);

    // invoke the petsc constructors
    int ierr = OptionsGetInt(PETSC_NULL,"-n",&size,&flg); CHKERRA(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,size,size,&A); CHKERRA(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,PETSC_DECIDE,size,&x); CHKERRA(ierr); 
    ierr = VecDuplicate(x,&b); CHKERRA(ierr);  
    // invoke setSize() on the Solver
    int solverOK = theSolver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:PetscSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    

    return result;    
}


int 
PetscSOE::setSizeParallel(int newSize, int N, 
			int d_nz, int *d_nnz,
			int o_nz, int *o_nnz)
{
    int result = 0;
    int oldSize = size;
    size = newSize;
    isFactored = 0;

    // delete the old	
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;

    // invoke the petsc destructors
    if (A != 0) MatDestroy(A);
    if (b != 0) VecDestroy(b);
    if (x != 0) VecDestroy(x);

    // create the new
    B = new double[size];
    X = new double[size];
	
    if (B == 0 || X == 0) {
      opserr << "WARNING PetscSOE::PetscSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      size = 0; 
      result = -1;
    }

    // zero the vectors
    for (int j=0; j<size; j++) {
	B[j] = 0;
	X[j] = 0;
    }

    if (vectX != 0) 
      delete vectX;

    if (vectB != 0) 
      delete vectB;
		
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);

    // invoke the petsc constructors
    int ierr = OptionsGetInt(PETSC_NULL,"-n",&size,&flg); CHKERRA(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,size,size,&A); CHKERRA(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,PETSC_DECIDE,size,&x); CHKERRA(ierr); 
    ierr = VecDuplicate(x,&b); CHKERRA(ierr);  
    // invoke setSize() on the Solver
    int solverOK = theSolver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:PetscSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    

    return result;    
}



int 
PetscSOE::addA(const Matrix &m, const ID &id, double fact)
{
  isFactored = 0;

    // check for a quick return 
    if (fact == 0.0)  return 0;

    
    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "PetscSOE::addA() - Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    int n = id.Size();
    int row;
    int col;
    double value;
    for (int i=0; i<n; i++) {
      row = id(i);
      if (row >= 0) {
	for (int j=0; j<n; j++) {
	  col = id(j);
	  if (col >= 0) {
	    value = m(i,j)*fact;
	    int ierr = MatSetValues(A,1,&row,1,&col,&value,ADD_VALUES); CHKERRA(ierr); 
	  }
	}
      }
    }

    return 0;
}

    
int 
PetscSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "PetscSOE::addB() - Vector and ID not of similar sizes\n";
	return -1;
    }    
    
    int n = id.Size();
    int row;
    double value;
    for (int i=0; i<n; i++) {
      row = id(i);
      if (row >= 0) {
	value = v(i) * fact;
	int ierr = VecSetValues(b,1,&row,&value,ADD_VALUES); CHKERRA(ierr); 
      }
    }

    return 0;
}


int
PetscSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (size != v.Size() ) {
	opserr << "PetscSOE::addB() - Vector not of appropriate size\n";
	return -1;
    }    
    
    int n = size;
    int row;
    double value;
    for (int i=0; i<n; i++) {
      value = v(i) * fact;
      int ierr = VecSetValues(b,1,&i,&value,INSERT_VALUES); CHKERRA(ierr); 
    }
    int ierr = VecAssemblyBegin(b); CHKERRA(ierr); 
    ierr = VecAssemblyEnd(b); CHKERRA(ierr); 

    return 0;
}


void 
PetscSOE::zeroA(void)
{
  isFactored = 0;
  int ierr = MatZeroEntries(A); CHKERRA(ierr);   
}
	
void 
PetscSOE::zeroB(void)
{
  double zero = 0.0;
  VecSet(&zero, b);
}


const Vector &
PetscSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL PetscSOE::getX - vectX == 0!";
	exit(-1);
    }    
    double *theData =0;
    int ierr = VecGetArray(x, &theData); CHKERRA(ierr);       
    for (int i=0; i<size; i++)
      X[i] = theData[i];
    ierr = VecRestoreArray(x, &theData); CHKERRA(ierr);             


    return *vectX;
}


const Vector &
PetscSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL PetscSOE::getB - vectB == 0!";
	exit(-1);
    }    
    double *theData =0;
    int ierr = VecGetArray(b, &theData); CHKERRA(ierr);       
    for (int i=0; i<size; i++)
      B[i] = theData[i];
    ierr = VecRestoreArray(b, &theData); CHKERRA(ierr);             

    return *vectB;
}


double 
PetscSOE::normRHS(void)
{
  this->getB();
  double norm =0.0;
  double *Bptr = B;
  for (int i=0; i<size; i++) {
    double Yi = *Bptr++;
    norm += Yi*Yi;
  }
  return sqrt(norm);
}    


void 
PetscSOE::setX(int loc, double value)
{
  int ierr = VecSetValues(x,1,&loc,&value,INSERT_VALUES); CHKERRA(ierr);   
}

int
PetscSOE::setSolver(PetscSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:PetscSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return solverOK;
	}
    }	
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
PetscSOE::sendSelf(Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    if (size != 0)
	opserr << "WARNING PetscSOE::sendSelf - does not send itself YET\n";
    return 0;
}


int 
PetscSOE::recvSelf(Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}





