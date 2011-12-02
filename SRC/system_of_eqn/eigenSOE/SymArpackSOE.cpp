// File: ~/system_of_eqn/linearSOE/LawSolver/SymArpackSOE.C
//
// Written: Jun Peng
// Created: 12/98
// Revision: A
//
// Description: This file contains the class definition for 
// LawSparseinSolver. It solves the SymArpackSOEobject by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) SymArpackSOE.C, revA"


#include <fstream.h>

#include "SymArpackSOE.h"
#include "SymArpackSolver.h"
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>
#include <Vector.h>

SymArpackSOE::SymArpackSOE(SymArpackSolver &the_Solver, AnalysisModel &aModel,
			 double theShift)
:EigenSOE(the_Solver, EigenSOE_TAGS_SymArpackSOE),
 size(0), nnz(0), colA(0), rowStartA(0), 
 factored(false), shift(theShift), theModel(&aModel),
 nblks(0), xblk(0), invp(0), diag(0), penv(0), rowblks(0),
 begblk(0), first(0)
{	
    the_Solver.setEigenSOE(*this);
}


SymArpackSOE::~SymArpackSOE()
{
  /*
    if (invp != 0) free((void *)invp);
    if (xblk != 0) free((void *)xblk); 
    if (penv != 0) {
      free((void *)penv[0]); 
      free((void *)penv); 
    }
    if (diag != 0) free((void *)diag); 
    if (rowblks != 0) free((void *)rowblks); 
    if (begblk != 0) free((void *)begblk); 
    if (first  != 0) free((void *)first); 

    if (rowStartA != 0) delete [] rowStartA;
    if (colA != 0) delete [] colA;
  */
}


int
SymArpackSOE::getNumEqn(void) const
{
    return size;
}

// extern "C" int symFactorization(int *fxadj, int *adjncy, int neq, int LSPARSE);
extern "C" int symFactorization(int *fxadj, int *adjncy, int neq, int LSPARSE, 
				int **xblkMY, int **invpMY, int **rowblksMY, 
				OFFDBLK ***begblkMY, OFFDBLK **firstMY, 
				double ***penvMY, double **diagMY);
				
int 
SymArpackSOE::setSize(Graph &theGraph)
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

    if (colA != 0)
      delete [] colA;

    colA = new int[newNNZ];	
    if (colA == 0) {
	cerr << "WARNING SymArpackSOE::SymArpackSOE :";
	cerr << " ran out of memory for colA with nnz = ";
      	cerr << newNNZ << " \n";
       	size = 0; nnz = 0;
       	result =  -1;
    } 
	
    factored = false;

    if (rowStartA != 0) 
      delete [] rowStartA;
    rowStartA = new int[size+1]; 
    if (rowStartA == 0) {
        cerr << "SymArpackSOE::ran out of memory for rowStartA." << endl;
        result = -1;
    }

    // fill in rowStartA and colA
    if (size != 0) {
        rowStartA[0] = 0;
        int startLoc = 0;
	int lastLoc = 0;

	for (int a=0; a<size; a++) {
	   theVertex = theGraph.getVertexPtr(a);
	   if (theVertex == 0) {
	        cerr << "WARNING:SymArpackSOE::setSize :";
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
    //   cout << "Enter DOF Numberer Type: \n";
    //   cout << "[1] Minimum Degree, [2] Nested Dissection, [3] RCM: ";
    int LSPARSE = 1;
    //   cin >> LSPARSE;

// call "C" function to form elimination tree and do the symbolic factorization.
//    nblks = symFactorization(rowStartA, colA, size, LSPARSE);
    nblks = symFactorization(rowStartA, colA, size, LSPARSE,
			     &xblk, &invp, &rowblks, &begblk, &first, &penv, &diag);

    // invoke setSize() on the Solver
    EigenSolver *theSolvr = this->getSolver();
    int solverOK = theSolvr->setSize();
    if (solverOK < 0) {
	cerr << "WARNING:BandArpackSOE::setSize :";
	cerr << " solver failed setSize()\n";
	return solverOK;
    } 


    return result;
}


int 
SymArpackSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return
	if (fact == 0.0)  
	return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	cerr << "SymArpackSOE::addA() ";
	cerr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    int *newID = new int [idSize];
    int *isort = new int[idSize];
    if (newID == 0 || isort ==0) {
        cerr << "WARNING SymArpackSOE::addA() :";
        cerr << " ran out of memory for vectors (newID, isort)";
        return -1;
    }

    int i;    
    for (i=0; i<idSize; i++) {
        newID[i] = id(i);
	 if (newID[i] >= 0)
	       newID[i] = invp[newID[i]];
    }
       
   long int  i_eq, j_eq, iadd;
   int  j, nee, lnee;
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
   for(i = 0, k = 0; i < lnee ; i++ )
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
	       int temp = isort[j];
	       isort[j] = isort[j+1];
	       isort[j+1] = temp;

/*	       isort[j] ^= isort[j+1] ;
	       isort[j+1] ^= isort[j] ;
	       isort[j] ^= isort[j+1] ;
*/
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
	    while (saveblk->row != i_eq) 
	         saveblk = saveblk->bnext ;
	 
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
	    {  loc = iloc + j_eq ;
	       *loc += m(it, jt) * fact;
            } else /* row segment */
	    {  while((j_eq >= (ptr->next)->beg) && ((ptr->next)->row == i_eq))
		  ptr = ptr->next ;
	       fpt = ptr->nz ;
	       fpt[j_eq - ptr->beg] += m(it,jt) * fact;
            }
         }
	 diag[i_eq] += m(ipos, ipos) * fact;
      }
   
    delete [] newID;
    delete [] isort;

    return 0;
}

    
int 
SymArpackSOE::addM(const Matrix &m, const ID &id, double fact)
{
    return this->addA(m, id, -shift);
}


double 
SymArpackSOE::getShift(void)
{
    return shift;
}

void 
SymArpackSOE::zeroA(void)
{
    factored = false;
}
	
void 
SymArpackSOE::zeroM(void)
{
  
}

int 
SymArpackSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int 
SymArpackSOE::recvSelf(int cTag, Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
    return 0;
}












