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
                                                                        
// $Revision: 1.6 $
// $Date: 2009-05-11 20:56:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 02/06

// Description: This file contains the implementation for MumpsSOE

#include <MumpsSOE.h>
#include <MumpsSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>

#include <stdlib.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

MumpsSOE::MumpsSOE(MumpsSolver &the_Solver, int _matType)
:LinearSOE(the_Solver, LinSOE_TAGS_MumpsSOE),
 size(0), nnz(0), 
 A(0), B(0), X(0), 
 colA(0), rowA(0), rowB(0), colStartA(0),
 vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false), matType(_matType)
{
  the_Solver.setLinearSOE(*this);
}

MumpsSOE::MumpsSOE()
 :LinearSOE(LinSOE_TAGS_MumpsSOE),
  size(0), nnz(0), 
  A(0), B(0), X(0), 
  colA(0), rowA(0), rowB(0), colStartA(0),
  vectX(0), vectB(0),
  Asize(0), Bsize(0),
  factored(false), matType(0)
{

}

MumpsSOE::MumpsSOE(int classTag)
 :LinearSOE(classTag),
  size(0), nnz(0), 
  A(0), B(0), X(0), 
  colA(0), rowA(0), rowB(0), colStartA(0),
  vectX(0), vectB(0),
  Asize(0), Bsize(0),
  factored(false), matType(0)
{

}


MumpsSOE::MumpsSOE(LinearSOESolver &the_Solver, int classTag, int _matType)
  :LinearSOE(the_Solver, classTag),
   size(0), nnz(0), 
   A(0), B(0), X(0), 
   colA(0), rowA(0), rowB(0), colStartA(0),
   vectX(0), vectB(0),
   Asize(0), Bsize(0),
   factored(false), matType(_matType)
{

}


MumpsSOE::~MumpsSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (colStartA != 0) delete [] colStartA;
    if (rowA != 0) delete []rowA;
    if (colA != 0) delete []colA;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;
}


int
MumpsSOE::getNumEqn(void) const
{
    return size;
}

int 
MumpsSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;
  size = theGraph.getNumVertex();
  
  // fist itearte through the vertices of the graph to get nnz
  Vertex *theVertex;
  int newNNZ = 0;
  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    const ID &theAdjacency = theVertex->getAdjacency();
    newNNZ += theAdjacency.Size() +1; // the +1 is for the diag entry
  }

  if (matType !=  0) {
    newNNZ -= size;
    newNNZ /= 2;
    newNNZ += size;
  }

  nnz = newNNZ;
  
  if (newNNZ > Asize) { // we have to get more space for A and rowA
    if (A != 0) delete [] A;
    if (rowA != 0) delete [] rowA;
    
    A = new double[newNNZ];
    rowA = new int[newNNZ];
    colA = new int[newNNZ];
    
    if (A == 0 || rowA == 0 || colA == 0) {
      opserr << "WARNING MumpsSOE::MumpsSOE :";
      opserr << " ran out of memory for A and rowA with nnz = ";
      opserr << newNNZ << " \n";
      size = 0; Asize = 0; nnz = 0;
      result =  -1;
    } 
    
    Asize = newNNZ;
  }
  
  // zero the matrix
  for (int i=0; i<Asize; i++)
    A[i] = 0;
  
  factored = false;
  
  if (size > Bsize) { // we have to get space for the vectors
    
    // delete the old	
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (colStartA != 0) delete [] colStartA;
    
    // create the new
    B = new double[size];
    X = new double[size];
    colStartA = new int[size+1]; 
    
    if (B == 0 || X == 0 || colStartA == 0) {
      opserr << "WARNING MumpsSOE::MumpsSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      size = 0; Bsize = 0;
      result =  -1;
    }
    else
      Bsize = size;
  }
  
  // zero the vectors
  for (int j=0; j<size; j++) {
    B[j] = 0;
    X[j] = 0;
  }
  
  // create new Vectors objects
  if (size != oldSize) {
    if (vectX != 0)
      delete vectX;
    
    if (vectB != 0)
      delete vectB;
    
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);	
  }
  
  // fill in colStartA and rowA
  if (size != 0) {
    colStartA[0] = 0;
    int startLoc = 0;
    int lastLoc = 0;
    for (int a=0; a<size; a++) {
      
      theVertex = theGraph.getVertexPtr(a);
      if (theVertex == 0) {
	opserr << "WARNING:MumpsSOE::setSize :";
	opserr << " vertex " << a << " not in graph! - size set to 0\n";
	size = 0;
	return -1;
      }
      
      int vertexTag = theVertex->getTag();
      rowA[lastLoc++] = vertexTag; // place diag in first
      const ID &theAdjacency = theVertex->getAdjacency();
      int idSize = theAdjacency.Size();
      
      // now we have to place the entries in the ID into order in rowA

      if (matType != 0) {

	for (int i=0; i<idSize; i++) {
	  int row = theAdjacency(i);
	  if (row > vertexTag) {
	    bool foundPlace = false;
	    // find a place in rowA for current col
	    for (int j=startLoc; j<lastLoc; j++)
	      if (rowA[j] > row) { 
		// move the entries already there one further on
		// and place col in current location
		for (int k=lastLoc; k>j; k--)
		  rowA[k] = rowA[k-1];
		rowA[j] = row;
		foundPlace = true;
		j = lastLoc;
	      }
	    
	    if (foundPlace == false) // put in at the end
	      rowA[lastLoc] = row;
	    lastLoc++;
	  }
	}

      } else {

	for (int i=0; i<idSize; i++) {
	  int row = theAdjacency(i);
	  bool foundPlace = false;
	  // find a place in rowA for current col
	  for (int j=startLoc; j<lastLoc; j++)
	    if (rowA[j] > row) { 
	      // move the entries already there one further on
	      // and place col in current location
	      for (int k=lastLoc; k>j; k--)
		rowA[k] = rowA[k-1];
	      rowA[j] = row;
	      foundPlace = true;
	      j = lastLoc;
	    }
	  if (foundPlace == false) // put in at the end
	    rowA[lastLoc] = row;
	  
	  lastLoc++;
	}
      }

      colStartA[a+1] = lastLoc;;	    
      startLoc = lastLoc;
    }
  }

  // fill in colA
  int count = 0;
  for (int i=0; i<size; i++)
    for (int k=colStartA[i]; k<colStartA[i+1]; k++)
      colA[count++] = i;
  
  // invoke setSize() on the Solver    
  LinearSOESolver *the_Solver = this->getSolver();
  int solverOK = the_Solver->setSize();
  if (solverOK < 0) {
    opserr << "WARNING:MumpsSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    
  
  return result;
}

int 
MumpsSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  
	return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "MumpsSOE::addA() ";
	opserr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }

    if (matType != 0) {

      if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	  int col = id(i);
	  if (col < size && col >= 0) {
	    int startColLoc = colStartA[col];
	    int endColLoc = colStartA[col+1];
	    
	    for (int j=0; j<idSize; j++) {
	      int row = id(j);
	      if (row >= col && row < size && row >= 0) {
		// find place in A using rowA
		for (int k=startColLoc; k<endColLoc; k++)
		  if (rowA[k] == row) {
		    A[k] += m(j,i);
		    k = endColLoc;
		  }
	      }
	    }  // for j		
	  } 
	}  // for i
      } else {
	for (int i=0; i<idSize; i++) {
	  int col = id(i);
	  if (col < size && col >= 0) {
	    int startColLoc = colStartA[col];
	    int endColLoc = colStartA[col+1];

	    for (int j=0; j<idSize; j++) {
	      int row = id(j);
	      if (row >= col && row <size && row >= 0) {
		// find place in A using rowA
		for (int k=startColLoc; k<endColLoc; k++)
		  if (rowA[k] == row) {
		    A[k] += fact * m(j,i);
		    k = endColLoc;
		  }
	      }
	    }  // for j		
	  } 
	}  // for i
      }


    } else {
    
      if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	  int col = id(i);
	  if (col < size && col >= 0) {
	    int startColLoc = colStartA[col];
	    int endColLoc = colStartA[col+1];

	    for (int j=0; j<idSize; j++) {
	      int row = id(j);
	      if (row <size && row >= 0) {
		// find place in A using rowA
		for (int k=startColLoc; k<endColLoc; k++)
		  if (rowA[k] == row) {
		    A[k] += m(j,i);
		    k = endColLoc;
		  }
	      }
	    }  // for j		
	  } 
	}  // for i
      } else {
	for (int i=0; i<idSize; i++) {
	  int col = id(i);
	  if (col < size && col >= 0) {
	    int startColLoc = colStartA[col];
	    int endColLoc = colStartA[col+1];
	    for (int j=0; j<idSize; j++) {
	      int row = id(j);
	      if (row <size && row >= 0) {
		// find place in A using rowA
		for (int k=startColLoc; k<endColLoc; k++)
		  if (rowA[k] == row) {
		    A[k] += fact * m(j,i);
		    k = endColLoc;
		  }
	      }
	    }  // for j		
	  } 
	}  // for i
      }
    }

    return 0;
}

    
int 
MumpsSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	opserr << "MumpsSOE::addB() ";
	opserr << " - Vector and ID not of similar sizes\n";
	return -1;
    }    

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	

    return 0;
}


int
MumpsSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != size) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<size; i++) {
	    B[i] = v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<size; i++) {
	    B[i] = -v(i);
	}
    } else {
	for (int i=0; i<size; i++) {
	    B[i] = v(i) * fact;
	}
    }	

    return 0;
}

void 
MumpsSOE::zeroA(void)
{
    double *Aptr = A;
    for (int i=0; i<nnz; i++)
	*Aptr++ = 0;

	factored = false;
}
	
void 
MumpsSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}

void 
MumpsSOE::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

void 
MumpsSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
MumpsSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL MumpsSOE::getX - vectX == 0";
	exit(-1);
    }
    return *vectX;
}

const Vector &
MumpsSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL MumpsSOE::getB - vectB == 0";
	exit(-1);
    }        
    return *vectB;
}

double 
MumpsSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    

int 
MumpsSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int 
MumpsSOE::recvSelf(int cTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)  
{
    return 0;
}

