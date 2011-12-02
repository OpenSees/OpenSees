// File: ~/system_of_eqn/linearSOE/symLinSolver/SymSparseLinSOE.C
//
// Written: Jun Peng
// Created: 12/98
// Revision: A
//
// Description: This file contains the class definition for 
// LawSparseinSolver. It solves the SymSparseLinSOEobject by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) SymSparseLinSOE.C, revA"


#include <fstream.h>

#include "SymSparseLinSOE.h"
#include "SymSparseLinSolver.h"
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

extern "C" {
   #include "FeStructs.h"
}

SymSparseLinSOE::SymSparseLinSOE(SymSparseLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_SymSparseLinSOE),
 size(0), nnz(0), B(0), X(0), colA(0), rowStartA(0),
 vectX(0), vectB(0), 
 Bsize(0), factored(false),
 nblks(0), xblk(0), invp(0), diag(0), penv(0), rowblks(0),
 begblk(0), first(0) 
{	
    the_Solver.setLinearSOE(*this);
}


SymSparseLinSOE::~SymSparseLinSOE()
{
	
    if (invp != 0) delete [] invp;
    if (xblk != 0) delete [] xblk;
    if (penv != 0) {
        delete [] penv[0];
        delete [] penv;
    }
    if (diag != 0) delete [] diag;
    if (rowblks != 0) delete [] rowblks;
    if (begblk != 0) delete [] begblk;
    if (first != 0) delete [] first;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;
    if (rowStartA != 0) delete [] rowStartA;
    if (colA != 0) delete [] colA;
}


int
SymSparseLinSOE::getNumEqn(void) const
{
    return size;
}

// extern "C" int symFactorization(int *fxadj, int *adjncy, int neq, int LSPARSE);

extern "C" int symFactorization(int *fxadj, int *adjncy, int neq, int LSPARSE, 
				int **xblkMY, int **invpMY, int **rowblksMY, 
				OFFDBLK ***begblkMY, OFFDBLK **firstMY, 
				double ***penvMY, double **diagMY);


int 
SymSparseLinSOE::setSize(Graph &theGraph)
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
	newNNZ += theAdjacency.Size(); 
    }
    nnz = newNNZ;
 
    colA = new int[newNNZ];	
    if (colA == 0) {
        cerr << "WARNING SymSparseLinSOE::SymSparseLinSOE :";
	cerr << " ran out of memory for colA with nnz = ";
      	cerr << newNNZ << " \n";
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
            cerr << "WARNING SymSparseLinSOE::SymSparseLinSOE :";
	    cerr << " ran out of memory for vectors (size) (";
	    cerr << size << ") \n";
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
	        cerr << "WARNING:SymSparseLinSOE::setSize :";
	        cerr << " vertex " << a << " not in graph! - size set to 0\n";
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
    
    // begin to choose different ordering schema.
    cout<< "Enter DOF Numberer Type: \n";
    cout<< "[1] Minimum Degree, [2] Nested Dissection, [3] RCM: ";
    int LSPARSE;
    cin >>LSPARSE;

    // call "C" function to form elimination tree and do the symbolic factorization.
    nblks = symFactorization(rowStartA, colA, size, LSPARSE,
			     &xblk, &invp, &rowblks, &begblk, &first, &penv, &diag);
    return result;
}


int 
SymSparseLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return
    if (fact == 0.0)  
    return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	cerr << "SymSparseLinSOE::addA() ";
	cerr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    int *newID = new int[idSize];
    int *isort = new int[idSize];
    if (newID == 0 || isort ==0) {
        cerr << "WARNING SymSparseLinSOE::SymSparseLinSOE :";
        cerr << " ran out of memory for vectors (newID, isort)";
        return -1;
    }

    for (int ii=0; ii<idSize; ii++) {
        newID[ii] = id(ii);
        if (newID[ii] >= 0)
             newID[ii] = invp[newID[ii]];
    }
       
   long int  i_eq, j_eq, iadd;
   int  i, j, nee, lnee;
   int  k, ipos, jpos;
   int  it, jt;
   int  iblk;
   int  jblk;
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

	    if (j_eq >= xblk[iblk]) /* diagonal block */
	    {  
	        loc = iloc + j_eq ;
		*loc += m(it,jt) * fact;
            } 
	    else /* row segment *//**/
	    { 
	        while((j_eq >= (ptr->next)->beg) && ((ptr->next)->row == i_eq))
		    ptr = ptr->next ;
		fpt = ptr->nz ;
		fpt[j_eq - ptr->beg] += m(it,jt) * fact;
            }
         }
	 diag[i_eq] += m(ipos,ipos) * fact;
      }
  	  
    delete [] newID;
    delete [] isort;

    return 0;
}

    
int 
SymSparseLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	cerr << "UmfpackGenLinSOE::addB() ";
	cerr << " - Vector and ID not of similar sizes\n";
	return -1;
    }    

    int *newID = new int[idSize];
    if (newID == 0) {
       cerr << "WARNING SymSparseLinSOE::SymSparseLinSOE :";
       cerr << " ran out of memory for vectors (newID)";
        return -1;
    }

    for (int i=0; i<idSize; i++) {
       newID[i] = id(i);
	if (newID[i] >= 0)
	    newID[i] = invp[newID[i]];
    }

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = newID[i];
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = newID[i];
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = newID[i];
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	

    delete [] newID;
    return 0;
}


int
SymSparseLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    if (v.Size() != size) {
	cerr << "WARNING BandGenLinSOE::setB() -";
	cerr << " incomptable sizes " << size << " and " << v.Size() << endl;
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
SymSparseLinSOE::zeroA(void)
{
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

const Vector &
SymSparseLinSOE::getX(void)
{
    if (vectX == 0) {
	cerr << "FATAL UmfpackGenLinSOE::getX - vectX == 0";
	exit(-1);
    }
    return *vectX;
}

const Vector &
SymSparseLinSOE::getB(void)
{
    if (vectB == 0) {
	cerr << "FATAL UmfpackGenLinSOE::getB - vectB == 0";
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


int
SymSparseLinSOE::setSymSparseLinSolver(SymSparseLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
        int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    cerr << "WARNING:SymSparseLinSOE::setSolver :";
	    cerr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
SymSparseLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}


int 
SymSparseLinSOE::recvSelf(int cTag, 
			  Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

