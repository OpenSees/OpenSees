// File: ~/system_of_eqn/linearSOE/symLinSolver/SymSparseLinSOE.CPP
//
// Written: Jun Peng  (junpeng@stanford.edu)
//          Advisor: Prof. Kincho H. Law
//          Stanford University
// Created: 12/1998
// Revised: 12/2001
//
// Description: This file contains the class definition for 
// SymSparseinSolver. It solves the SymSparseLinSOEobject by calling
// some "C" functions. The solver used here is a symmtric generalized sparse
// solver. The user can choose three different ordering scheme.
//
// What: "@(#) SymSparseLinSOE.C, revA"


#include <fstream>
#include <stdlib.h>

#include <SymSparseLinSOE.h>
#include <SymSparseLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


SymSparseLinSOE::SymSparseLinSOE(SymSparseLinSolver &the_Solver, int lSparse)
:LinearSOE(the_Solver, LinSOE_TAGS_SymSparseLinSOE),
 size(0), nnz(0), B(0), X(0), colA(0), rowStartA(0),
 vectX(0), vectB(0), 
 Bsize(0), factored(false),
 nblks(0), xblk(0), invp(0), diag(0), penv(0), rowblks(0),
 begblk(0), first(0) 
{	
    the_Solver.setLinearSOE(*this);
    this->LSPARSE = lSparse;
}


/* A destructor for cleanning memory.
 * For diag and penv, it is rather straightforward to clean.
 * For row segments, since the memory of nz is allocated for each
 * row, the deallocated needs some special care.
 */
SymSparseLinSOE::~SymSparseLinSOE()
{
    // free the diagonal vector
    if (diag != NULL) free(diag);

    // free the diagonal blocks
    if (penv != NULL) {
	if (penv[0] != NULL) {
	    free(penv[0]);
	}
        free(penv);
    } 

    // free the row segments.
    OFFDBLK *blkPtr = first;
    OFFDBLK *tempBlk;
    int curRow = -1;

    while (1) {
      if (blkPtr->next == blkPtr) {
	if (blkPtr != NULL) {
	  free(blkPtr);
	}     
	break;
      }

      tempBlk = blkPtr->next;
      if (blkPtr->row != curRow) {
	if (blkPtr->nz != NULL) {
	  free(blkPtr->nz);
	}
	curRow = blkPtr->row;
      }
      
      free(blkPtr);

      blkPtr = tempBlk;
    }

    // free the "C" style vectors.
    if (xblk != 0)  free(xblk);
    if (rowblks != 0)   free(rowblks);
    if (invp != 0)  free(invp);
    
    // free the "C++" style vectors.
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;
    if (rowStartA != 0) delete [] rowStartA;
    if (colA != 0) delete [] colA;
}


int SymSparseLinSOE::getNumEqn(void) const
{
    return size;
}


/* the symbolic factorization method defined in <symbolic.c>
 */
extern "C" int symFactorization(int *fxadj, int *adjncy, int neq, int LSPARSE, 
				int **xblkMY, int **invpMY, int **rowblksMY, 
				OFFDBLK ***begblkMY, OFFDBLK **firstMY, 
				double ***penvMY, double **diagMY);


/* Based on the graph (the entries in A), set up the pair (rowStartA, colA).
 * It is the same as the pair (ADJNCY, XADJ).
 * Then perform the symbolic factorization by calling symFactorization().
 */
int SymSparseLinSOE::setSize(Graph &theGraph)
{

    int result = 0;
    int oldSize = size;
    size = theGraph.getNumVertex();

    // first itearte through the vertices of the graph to get nnz
    Vertex *theVertex;
    int newNNZ = 0;
    VertexIter &theVertices = theGraph.getVertices();
    while ((theVertex = theVertices()) != 0) {
        const ID &theAdjacency = theVertex->getAdjacency();
	newNNZ += theAdjacency.Size(); 
    }
    nnz = newNNZ;
 
    colA = new int[newNNZ];	
    if (colA == 0) {
        opserr << "WARNING SymSparseLinSOE::SymSparseLinSOE :";
	opserr << " ran out of memory for colA with nnz = ";
      	opserr << newNNZ << " \n";
       	size = 0; nnz = 0;
       	result =  -1;
    } 
	
    factored = false;
    
    if (size > Bsize) { // we have to get space for the vectors
	
	// delete the old	
	if (B != 0) delete [] B;
	if (X != 0) delete [] X;
	if (rowStartA != 0) delete [] rowStartA;

	// create the new
	B = new double[size];
	X = new double[size];
	rowStartA = new int[size+1]; 
	
	if (B == 0 || X == 0 || rowStartA == 0) {
            opserr << "WARNING SymSparseLinSOE::SymSparseLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    size = 0; Bsize = 0;
	    result =  -1;
	} else {
	    Bsize = size;
	}
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

    // fill in rowStartA and colA
    if (size != 0) {
        rowStartA[0] = 0;
        int startLoc = 0;
	int lastLoc = 0;

	for (int a=0; a<size; a++) {
	   theVertex = theGraph.getVertexPtr(a);
	   if (theVertex == 0) {
	        opserr << "WARNING:SymSparseLinSOE::setSize :";
	        opserr << " vertex " << a << " not in graph! - size set to 0\n";
	        size = 0;
	        return -1;
	   }

	   const ID &theAdjacency = theVertex->getAdjacency();
	   int idSize = theAdjacency.Size();
	
	// now we have to place the entries in the ID into order in colA
	   for (int i=0; i<idSize; i++) {
	      int row = theAdjacency(i);
	      bool foundPlace = false;
	 
	      for (int j=startLoc; j<lastLoc; j++)
	          if (colA[j] > row) { 
	      // move the entries already there one further on
	      // and place col in current location
	              for (int k=lastLoc; k>j; k--)
		          colA[k] = colA[k-1];
                      colA[j] = row;
		      foundPlace = true;
    	              j = lastLoc;
		  }
		  
	      if (foundPlace == false) // put in at the end
	      	   colA[lastLoc] = row;

	      lastLoc++;
	   }
	   rowStartA[a+1] = lastLoc;;	    
	   startLoc = lastLoc;
	}
    }
    
    // call "C" function to form elimination tree and to do the symbolic factorization.
    nblks = symFactorization(rowStartA, colA, size, this->LSPARSE,
			     &xblk, &invp, &rowblks, &begblk, &first, &penv, &diag);

    return result;
}


/* Perform the element stiffness assembly here.
 */
int SymSparseLinSOE::addA(const Matrix &in_m, const ID &in_id, double fact)
{
   // check for a quick return
   if (fact == 0.0)  
       return 0;

   int idSize = in_id.Size();
   if (idSize == 0)  return 0;

   // check that m and id are of similar size
   if (idSize != in_m.noRows() && idSize != in_m.noCols()) {
       opserr << "SymSparseLinSOE::addA() ";
       opserr << " - Matrix and ID not of similiar sizes\n";
       return -1;
   }

   // construct m and id based on non-negative id values.
   int newPt = 0;
   int *id = new int[idSize];
   
   for (int jj = 0; jj < idSize; jj++) {
       if (in_id(jj) >= 0 && in_id(jj) < size) {
	   id[newPt] = in_id(jj);
	   newPt++;
       }
   }

   idSize = newPt;
   if (idSize == 0)  return 0;
   double *m = new double[idSize*idSize];

   int newII = 0;
   for (int ii = 0; ii < in_id.Size(); ii++) {
       if (in_id(ii) >= 0 && in_id(ii) < size) {

	   int newJJ = 0;
	   for (int jj = 0; jj < in_id.Size(); jj++) {
	       if (in_id(jj) >= 0 && in_id(jj) < size) {
		   m[newII*idSize + newJJ] = in_m(ii, jj);
		   newJJ++;
	       }
	   }
	   newII++;
       }
   }

   // forming the new id based on invp.

   int *newID = new int[idSize];
   int *isort = new int[idSize];
   if (newID == 0 || isort ==0) {
       opserr << "WARNING SymSparseLinSOE::SymSparseLinSOE :";
       opserr << " ran out of memory for vectors (newID, isort)";
       return -1;
   }

   for (int kk=0; kk<idSize; kk++) {
       newID[kk] = id[kk];
       if (newID[kk] >= 0)
	   newID[kk] = invp[newID[kk]];
   }
   
   long int  i_eq, j_eq;
   int  i, j, nee, lnee;
   int  k, ipos, jpos;
   int  it, jt;
   int  iblk;
   OFFDBLK  *ptr;
   OFFDBLK  *saveblk;
   double  *fpt, *iloc, *loc;

   nee = idSize;
   lnee = nee;
   
   /* initialize isort */
   for( i = 0, k = 0; i < lnee ; i++ )
   {
       if( newID[i] >= 0 ) {
	   isort[k] = i;
	   k++;
       }
   }
      
   lnee = k;

   /* perform the sorting of isort here */
   i = k - 1;
   do
   {
       k = 0 ;
       for (j = 0 ; j < i ; j++)
       {  
	   if ( newID[isort[j]] > newID[isort[j+1]]) {  
	       isort[j] ^= isort[j+1] ;
	       isort[j+1] ^= isort[j] ;
	       isort[j] ^= isort[j+1] ;
	       k = j ;
	   }
      }
      i = k ;
   }  while ( k > 0) ;

      i = 0 ;
      ipos = isort[i] ;
      k = rowblks[newID[ipos]] ;
      saveblk  = begblk[k] ;

      /* iterate through the element stiffness matrix, assemble each entry */
      for (i=0; i<lnee; i++)
      { 
	 ipos = isort[i] ;
         i_eq = newID[ipos] ;
	 iblk = rowblks[i_eq] ;
	 iloc = penv[i_eq +1] - i_eq ;
	 if (k < iblk)
	    while (saveblk->row != i_eq) saveblk = saveblk->bnext ;
	 
	 ptr = saveblk ;
	 for (j=0; j< i ; j++)
	 {   
	    jpos = isort[j] ;
	    j_eq = newID[jpos] ;

	    if (ipos > jpos) {
	        jt = ipos;
		it = jpos;
	    } else {
	        it = ipos;
		jt = jpos;
	    }

	    if (j_eq >= xblk[iblk]) /* diagonal block (profile) */
	    {  
	        loc = iloc + j_eq ;
		*loc += m[it*idSize + jt] * fact;
            } 
	    else /* row segment */
	    { 
	        while((j_eq >= (ptr->next)->beg) && ((ptr->next)->row == i_eq))
		    ptr = ptr->next ;
		fpt = ptr->nz ;
		fpt[j_eq - ptr->beg] += m[it*idSize + jt] * fact;
            }
         }
	 diag[i_eq] += m[ipos*idSize + ipos] * fact; /* diagonal element */
      }
  	  
    delete [] newID;
    delete [] isort;
    delete [] m;
    delete [] id;

    return 0;
}

    
/* assemble the force vector B (A*X = B).
 */
int SymSparseLinSOE::addB(const Vector &in_v, const ID &in_id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = in_id.Size();    
    // check that m and id are of similar size
    if (idSize != in_v.Size() ) {
	opserr << "SymSparseLinSOE::addB() ";
	opserr << " - Vector and ID not of similar sizes\n";
	return -1;
    } 

   // construct v and id based on non-negative id values.
   int newPt = 0;
   int *id = new int[idSize];
   double *v = new double[idSize];

   for (int ii = 0; ii < idSize; ii++) {
       if (in_id(ii) >= 0 && in_id(ii) < size) {
	   id[newPt] = in_id(ii);
	   v[newPt] = in_v(ii);
	   newPt++;
       }
   }

   idSize = newPt;
   if (idSize == 0)  {
       delete [] id;
       delete [] v;
       return 0;
   }
    int *newID = new int[idSize];
    if (newID == 0) {
       opserr << "WARNING SymSparseLinSOE::SymSparseLinSOE :";
       opserr << " ran out of memory for vectors (newID)";
        return -1;
    }

    for (int i=0; i<idSize; i++) {
       newID[i] = id[i];
	if (newID[i] >= 0)
	    newID[i] = invp[newID[i]];
    }

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = newID[i];
	    if (pos <size && pos >= 0)
		B[pos] += v[i];
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = newID[i];
	    if (pos <size && pos >= 0)
		B[pos] -= v[i];
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = newID[i];
	    if (pos <size && pos >= 0)
		B[pos] += v[i] * fact;  // assemble
	}
    }	

    delete [] newID;
    delete [] v;
    delete [] id;

    return 0;
}


int
SymSparseLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    if (v.Size() != size) {
	opserr << "WARNING SymSparseLinSOE::setB() -";
	opserr << " incomptable sizes " << size << " and " << v.Size() << endln;
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


/* It is used to set all the entries of A to be zero.
 * This method will be called if the structure of A stays the same while the
 * value of A needs to be changed.  e.g. if Newton-Raphson method is used
 * for analysis, after each iteration, the A matrix needs to be updated.
 */
void SymSparseLinSOE::zeroA(void)
{
    memset(diag, 0, size*sizeof(double));

    int profileSize = penv[size] - penv[0];
    memset(penv[0], 0, profileSize*sizeof(double));
    
    OFFDBLK *blkPtr = first;
    int rLen = 0;
    while (1) {
        if (blkPtr->beg == size)  break;
	rLen = xblk[rowblks[blkPtr->beg]+1] - blkPtr->beg;
	memset(blkPtr->nz, 0, rLen*sizeof(double));

	blkPtr = blkPtr->next;
    }

    factored = false;
}
	
void 
SymSparseLinSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}

void 
SymSparseLinSOE::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

void
SymSparseLinSOE::setX(const Vector &x)
{
    if (x.Size() == size && vectX != 0) 
        *vectX = x;
}


const Vector &
SymSparseLinSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL SymSparseLinSOE::getX - vectX == 0";
	exit(-1);
    }
    return *vectX;
}

const Vector &
SymSparseLinSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL SymSparseLinSOE::getB - vectB == 0";
	exit(-1);
    }        
    return *vectB;
}

double 
SymSparseLinSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
}    


/* Create a linkage between SOE and Solver.
 */
int SymSparseLinSOE::setSymSparseLinSolver(SymSparseLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
        int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:SymSparseLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
SymSparseLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    // not implemented.
    return 0;
}


int 
SymSparseLinSOE::recvSelf(int cTag, 
			  Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // not implemented.
    return 0;
}

