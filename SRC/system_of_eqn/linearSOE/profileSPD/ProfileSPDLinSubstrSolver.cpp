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
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSubstrSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/ProfileSPD/ProfileSPDLinSubstrSolver.C
//
// Written: fmk 
// Created: Febuary 1997
// Revision: A
//
// Description: This file contains the implementation for ProfileSPDLinSubstrSolver


#include <ProfileSPDLinSubstrSolver.h>
#include <ProfileSPDLinSOE.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

ProfileSPDLinSubstrSolver::ProfileSPDLinSubstrSolver(double tol)
:ProfileSPDLinDirectSolver(tol),
 DomainSolver(SOLVER_TAGS_ProfileSPDLinSubstrSolver),
 dSize(0),DU(0),Aext(0),Yext(0)
{

}

    
ProfileSPDLinSubstrSolver::~ProfileSPDLinSubstrSolver()
{
}

int
ProfileSPDLinSubstrSolver::solve()
{
    return this->ProfileSPDLinDirectSolver::solve();
}

int
ProfileSPDLinSubstrSolver::setSize(void)
{
    return this->ProfileSPDLinDirectSolver::setSize();
}    

/* ProfileSPDLinSubstrSolver::condenseA(int numInt)
**
** purpose: A procedure which takes the stifness matrix A = | A11 A12 |
**							    | A21 A22 |
**	    and does the following:
**
**		1) replaces A11 with D1 & U11, where A11 = U11'*D1*U11
**
**		2) replaces A12 with M, where M = inv(U11')*A12
**
**		3) replaces A22 with Kdash, where Kdash = A22 - A21*inv(A11)*A12
**							= A22 - M'*inv(D1)*M
**
** inputs: K[]          - the nxn matrix symmetric matrix  
**         Profile[]    - a vector outlining the profile of K
**         n         	- the size of the system
**	   in		- the no of internal dof (<n)
**
** outputs: int    0 if no error
**                -1 if diag term < MINDIAG
**                -2 if diag tem becomes <= 0
**
*/


int 
ProfileSPDLinSubstrSolver::condenseA(int numInt)
{
  /*
for (int iii=0; iii<theSOE->size; iii++) {
  int rowiiitop = RowTop[iii];
  double *ajiptr = topRowPtr[iii];
  cerr << "\n COLUMN " << iii << " TopRow " << rowiiitop << " -> ";
  for (int jjj = rowiiitop; jjj <=iii; jjj++)
    cerr << *ajiptr++ << " ";
}
cerr << endl;
*/

    
    // check for quick return
    if (theSOE == 0)
	return -1;

    if (numInt == 0) {
	theSOE->numInt = numInt;	
	return 0;
    }

    if (dSize != size) {
	if (DU != 0) delete [] DU;
	DU = new double[numInt];
	if (DU == 0) {
	    cerr << "ProfileSPDLinSubstrSolver::condenseA()";
	    cerr << " - ran out of memory for work space\n";	    
	    return -1;
	}
	dSize = numInt;
    }

		
    //
    //  form D1 & U11, store in A11
    //    - where A11 = U11'*D1*U11
    //    - done using Crout decomposition
    //





    this->factor(numInt);
    
    /*
     *  form M, leave in A12
     *   - where M = inv(U11')*A12
     */
    int i;

    for (i=numInt; i<size; i++) {
	
	int rowitop = RowTop[i];
	double *ajiPtr = topRowPtr[i];

	int jstrt = rowitop;
	if (rowitop == 0) {
	    jstrt = 1;
	    ajiPtr++;
	} else {
	    jstrt = rowitop;
	}

	
	for (int j=jstrt; j<numInt; j++) {
	    double tmp = *ajiPtr;
	  
	    int rowjtop = RowTop[j];

	    double *akiPtr, *akjPtr;
	   
	    if (rowitop > rowjtop) {
		akiPtr = topRowPtr[i];
		akjPtr = topRowPtr[j] + (rowitop-rowjtop);
		
		for (int k=rowitop; k<j; k++) 
		    tmp -= *akjPtr++ * *akiPtr++ ;
	      
	    } else {
		akiPtr = topRowPtr[i] + (rowjtop-rowitop);
		akjPtr = topRowPtr[j];
	       
		for (int k=rowjtop; k<j; k++) 
		    tmp -= *akjPtr++ * *akiPtr++ ;	    
	    }

	    *ajiPtr++ = tmp;
	}
    }


    /*
     * Now form K*, leave in A22
     *  - where K* = A22 - M'*inv(D11)*M
     */

    for (i=numInt; i<size; i++) {

	int rowitop = RowTop[i];
	double *ajiPtr =  topRowPtr[i];;

	int jstrt;

	if (rowitop < numInt) {
	    ajiPtr += (numInt-rowitop);
	    jstrt = numInt;
	}
	else
	  jstrt = rowitop;

	double *DUPtr = DU; 
	double *akiPtr =  topRowPtr[i];

	for (int k=rowitop; k<numInt; k++)
	  *DUPtr++ = *akiPtr++ * invD[k];

       
	for (int j=jstrt; j<=i; j++) {
	   
	    double tmp = *ajiPtr;
	    int rowjtop = RowTop[j];
	    double *akiPtr, *akjPtr;
	    
	    if (rowitop > rowjtop) {
	        akiPtr = DU;
		akjPtr = topRowPtr[j] + (rowitop-rowjtop);
		
		for (int k=rowitop; k<numInt; k++) 
		  tmp -= *akjPtr++ * *akiPtr++;

	    } else {
		akiPtr = &DU[rowjtop-rowitop];
		akjPtr = topRowPtr[j];
		
		for (int k=rowjtop; k<numInt; k++) 
		    tmp -= *akjPtr++ * *akiPtr++;
	    }

	    *ajiPtr++ = tmp;
	}
    }      


    theSOE->isAcondensed = true;
    theSOE->numInt = numInt;

cout << "ProfileSPDLinSubstrSolver::condenseA  numDOF: " << size << "  numInt: " << numInt << "  numExt: " << size-numInt << endl;

    return 0;

}




/* name: condenseRHS
**
** purpose: A procedure which takes the stifness matrix A = | D1\U11   M  | 
**							    | A21    A22* |
**	    
**	    and load vector B = | B1 |
**				| B2 |
**	    does the following:
**
**		1) replaces R1 with R1*, where R1* = inv(D11)inv(U11')R1
**
**		2) replaces R2 with R2*, where R2* = R2 - A21*inv(A11)*R1
**						   = R2 - M'*R1*
**
** inputs: A[nxn] - the nxn matrix symmetric matrix  STORED AS ROW VECTOR IN COL MAJOR ORDER
**	   R[n]   - the nX1 load vector
**         n         	- the size of the system
**	   in		- the no of internal dof (<n)
**
** outputs: int    0 if no error
**                -1 if diag term < MINDIAG
**                -2 if diag tem becomes <= 0
**
*/

int 
ProfileSPDLinSubstrSolver::condenseRHS(int numInt, Vector *v)
{

    // check for quick return
    if (theSOE == 0)
	return -1;

    if (numInt == 0) {
	theSOE->numInt = numInt;
	return 0;
    }

    // check A has been condensed & numInt was numInt given
    if (theSOE->isAcondensed != true) {
	int ok = this->condenseA(numInt);
	if (ok < 0) {
	    cerr << "ProfileSPDLinSubstrSolver::condenseRHS()";
	    cerr << " - failed to condenseA\n";
	    return ok;
	}
    }

    if (theSOE->numInt != numInt) {	
	cerr << "ProfileSPDLinSubstrSolver::condenseRHS()";
	cerr << " - numInt " << numInt << "does not agree with condensedA";
	cerr << " numInt " << theSOE->numInt << endl;
	return -1;
    }


    // set some pointers
    double *B = theSOE->B;

    //
    // form Y1*, leaving in Y1
    // - simply a triangular solve
    //

    // do forward substitution 
    for (int i=1; i<numInt; i++) {
	
	int rowitop = RowTop[i];	    
	double *ajiPtr = topRowPtr[i];
	double *bjPtr  = &B[rowitop];  
	double tmp = 0;	    
	
	for (int j=rowitop; j<i; j++) 
	    tmp -= *ajiPtr++ * *bjPtr++; 
	
	B[i] += tmp;
    }

    // divide by diag term 
    double *bjPtr = B; 
    double *aiiPtr = invD;
    for (int j=0; j<numInt; j++) 
	*bjPtr++ = *aiiPtr++ * B[j];

    //
    // Now form Y2*, leaving in Y2
    //

    for (int ii=numInt; ii<size; ii++) {
	int rowitop = RowTop[ii];	    
	double *ajiPtr = topRowPtr[ii];
	double *bjPtr  = &B[rowitop];  
	double tmp = 0;	    
	
	for (int j=rowitop; j<numInt; j++)	
	    tmp -= *ajiPtr++ * *bjPtr++;

	B[ii] += tmp;
    }
    /*
      for (int iii=0; iii<theSOE->size; iii++) {
      cerr << "COLUMN " << iii << " Biii -> " << B[iii] << endl;
      }
      cerr << endl;
      */

    return 0;
}


int 
ProfileSPDLinSubstrSolver::computeCondensedMatVect(
			       int numInt, const Vector &theVect) 
{
    cerr << "ProfileSPDLinSubstrSolver::computeCondensedMatVect() -";
    cerr << " not implemented yet\n";
    return -1;
}


const Matrix &
ProfileSPDLinSubstrSolver::getCondensedA(void)
{
    int numInt = theSOE->numInt;
    int matSize = size - numInt;

    // check Aext exists, if not create it
    if (Aext == 0) {
	Aext = new Matrix(matSize,matSize);
	
	if (Aext == 0 || Aext->noRows() == 0) {
	    cerr << "ProfileSPDLinSubstrSolver::getCondensedA";
	    cerr << "- ran out of memory for matSize " << matSize << " \n";
	    exit(-1);
	}
    }
    
    // check that current Aext is big enough if not enllarge
    if (Aext->noRows() != matSize) {
	delete Aext;
	Aext = new Matrix(matSize,matSize);
	
	if (Aext == 0 || Aext->noRows() == 0) {
	    cerr << "ProfileSPDLinSubstrSolver::getCondensedA";
	    cerr << "- ran out of memory for matSize " << matSize << " \n";
	    exit(-1);
	}
    }
    
    // set the components of Aee to be the matrix
    Aext->Zero();
    for (int i=numInt; i<size; i++) {
	int col = i - numInt;

	int rowitop = RowTop[i];
	double *akiPtr = topRowPtr[i];
	
	int row =0;
	if (rowitop < numInt) 
	    akiPtr += (numInt - rowitop);
	else
	    row = rowitop - numInt; // numInt FORTRAN index, rowitop C index
	    
	while (row < col){
	    double aki = *akiPtr++;
	    (*Aext)(row,col) = aki;
	    (*Aext)(col,row) = aki;	    
	    row++;
	}
	(*Aext)(col,row) = *akiPtr;
    }
    return *Aext;
}


const Vector &
ProfileSPDLinSubstrSolver::getCondensedRHS(void)
{
    int numInt = theSOE->numInt;
    int matSize = size - numInt;

    double *Y = &(theSOE->B[numInt]);

    // check Yext exists, if not create it
    if (Yext == 0) {
	Yext = new Vector(Y,matSize);
	
	if (Yext == 0 || Yext->Size() == 0) {
	    cerr << "ProfileSPDLinSubstrSolver::getCondensedRHS";
	    cerr << "- ran out of memory for vector Size " << matSize << " \n";
	    exit(-1);
	}
    }
    
    // check that current Yext is big enough if not enllarge
    if (Yext->Size() != matSize) {
	delete Yext;
	Yext = new Vector(Y,matSize);
	
	if (Yext == 0 || Yext->Size() == 0) {
	    cerr << "ProfileSPDLinSubstrSolver::getCondensedRHS";
	    cerr << "- ran out of memory for vect Size " << matSize << " \n";
	    exit(-1);
	}
    }
    
    return *Yext;
}

const Vector &
ProfileSPDLinSubstrSolver::getCondensedMatVect(void)
{
    cerr << "ProfileSPDLinSubstrSolver::computeCondensedMatVect() -";
    cerr << " not implemented yet\n";
    exit(-1);

    // needed for a strict compiler, program will exit before this
    Vector *errVect = new Vector(0); 
    return *errVect;
      
}


int 
ProfileSPDLinSubstrSolver::setComputedXext(const Vector &xExt)
{
    if (xExt.Size() != (size - theSOE->numInt)) {
	cerr << "ProfileSPDLinSubstrSolver::setComputedxExt() :";
	cerr << " - size mismatch " << xExt.Size() << " and ";
	cerr << size-theSOE->numInt << endl;
	return -1;
    }

    double *xPtr = &(theSOE->X)[theSOE->numInt];
    for (int i=0; i<xExt.Size(); i++)
	*xPtr++ = xExt(i);

    return 0;
}    


/* name: S O L V r I
**
** purpose: A procedure which takes the stifness matrix K = | D1\U11   M  | 
**							    | K21    K22* |
**	    
**	    and load/displ vector R = | R1* |
**				      | r2  |
**	    does the following:
**
**		1) replaces R1* with r1, where r1 = inv(K11){R1 - K12*r2}
**						  = inv(D11)inv(U11'){(D1)R1* - Mr2}
**
** inputs: K[nxn] - the nxn matrix symmetric matrix  STORED AS ROW VECTOR IN COL MAJ ORD
**	   R[n]   - the nX1 load vector
**         n         	- the size of the system
**	   in		- the no of internal dof (<n)
**
**
*/
int 
ProfileSPDLinSubstrSolver::solveXint(void)
{

  /*
   * form X1 = (D1)Y1* - M *X2, store in X1*
   */

    int numInt = theSOE->numInt;
    double *X = theSOE->X;
    double *B = theSOE->B;
    
    for (int j=0; j<numInt; j++) 
	X[j] = B[j]/invD[j];

    for (int i=numInt; i<size; i++) {
	
	int rowitop = RowTop[i];
	double *ajiPtr = topRowPtr[i];
	double *XjPtr  = &X[rowitop];
	double Xi = X[i];

	for (int j=rowitop; j<numInt; j++)
	    *XjPtr++ -= *ajiPtr++ * Xi;
    }

    for (int l=0; l<numInt; l++) 
       X[l] = X[l]*invD[l];

  /*
   * form inv(U11)*A, store in A
   * - simply a triangular solve (back sub)
   */
    for (int k=(numInt-1); k>0; k--) {

	int rowktop = RowTop[k];
	double Xk = X[k];
	double *ajiPtr = topRowPtr[k]; 		

	for (int j=rowktop; j<k; j++) 
	    X[j] -= *ajiPtr++ * Xk;
    }   	     
    return 0;
}

int
ProfileSPDLinSubstrSolver::getClassTag(void) const
{
    return SOLVER_TAGS_ProfileSPDLinSubstrSolver;
}

int
ProfileSPDLinSubstrSolver::sendSelf(int cTag,
				    Channel &theChannel)
{
    if (size != 0)
	cerr << "ProfileSPDLinSubstrSolver::sendSelf - does not send itself YET\n"; 
    return 0;
}


int 
ProfileSPDLinSubstrSolver::recvSelf(int cTag,
				    Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
    return 0;
}
