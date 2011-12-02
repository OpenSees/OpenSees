#include <stdio.h>
#include <string.h>
#include "snopt.h"
#include "toyOptfunction.h"


extern SNOPTAnalysis * theSNOPTAnalysis;// = new testaa();

int toyOptusrf_(integer    *Status, integer *n,    doublereal x[],
	     integer    *needF,  integer *neF,  doublereal F[],
	     integer    *needG,  integer *neG,  doublereal G[],
	     char       *cu,     integer *lencu,
	     integer    iu[],    integer *leniu,
	     doublereal ru[],    integer *lenru )
{

	int i;

	int result = theSNOPTAnalysis->updateXFG(x);

	if (result < 0) {
		opserr << "toyOptFunction:SNOPTAnalysis->update wrong!" << endln;
	
		*Status = -1;
	//	return -1;
	}

    double * tempF = theSNOPTAnalysis->getScaledFunctionPtr();
	for( i=0; i< *neF; i++)
		F[i]=tempF[i];


/*	static jjj=0;
	jjj++;
	opserr<<"Fata::The "<<jjj<<" th iteration, F is:" <<F[0]<<endln; //<<","<<F[1]<<","<<F[2]<<endln;
*/


	return 0;
}

int toyOptusrfg_( integer    *Status, integer *n,    doublereal x[],
	       integer    *needF,  integer *neF,  doublereal F[],
	       integer    *needG,  integer *neG,  doublereal G[],
	       char       *cu,     integer *lencu,
	       integer    iu[],    integer *leniu,
	       doublereal ru[],    integer *lenru )
{
  //==================================================================
  // Computes the nonlinear objective and constraint terms for the toy
  // problem featured in the SnoptA users guide.
  // neF = 3, n = 2.
  //
  //   Minimize     1/2*x'*x
  //
  //   subject to   g(x)=0 i.e., (x1-3)^2+(x2-6)^2+x3^2=0
  //
  // The triples (g(k),iGfun(k),jGvar(k)), k = 1:neG, define
  // the sparsity pattern and values of the nonlinear elements
  // of the Jacobian.
  //==================================================================




if( *needF > 0 ) {

	int result = theSNOPTAnalysis->updateXFG(x);

	if (result < 0) {
		opserr << "toyOptFunction:SNOPTAnalysis->update wrong!" << endln;
		*Status = -1;
	//	return -1;
	}
    double * tempF = theSNOPTAnalysis->getScaledFunctionPtr();
	for(int i=0; i< *neF; i++)
		F[i]=tempF[i];

// 	static iii=0;
//	iii++;
//	opserr<<"The "<<iii<<" th iteration, F is:" <<F[0]<<","<<F[1]<<","<<F[2]<<endln;
//	opserr<<"toyOptFunction, F is :"<<*F<<endln;
} //if *needF

//------------------ grad of [F;G] ---------------------------------------------
 
  if( *needG > 0 ){
	  if (*needF<=0) {
		  opserr<<"needG but no needF !!   in snopt userfunction calling! "<<endln; 
		  //exit(-1);
	  }
	double * tempG = theSNOPTAnalysis->getScaledGradientPtr(); 
	for(int i=0; i< (*n)*(*neF); i++){
		G[i] = tempG[i];
		if (G[i]>1.0e10) {
			opserr<<"warning: G is large the 1e10 in touoptfunction. something maybe wrong"<<endln;
		}
	}

//	opserr<<"toyOptFunction, G is :"<<G[0]<<" "<<G[1]<<" "<<G[2]<<" "<<G[3]<<" "<<G[4]<<" "<<G[5]<<" "endln;

/*    // iGfun[*neG] = 1
    // jGvar[*neG] = 1
    *neG = 0;
    G[*neG] = x[0];

    // iGfun[*neG] = 1
    // jGvar[*neG] = 2
    *neG = *neG + 1;
    G[*neG] = x[1];

    // iGfun[*neG] = 1
    // jGvar[*neG] = 3
    *neG = *neG + 1;
    G[*neG] = x[2];

    // iGfun[*neG] = 2
    // jGvar[*neG] = 1
    *neG = *neG + 1;
	//---------------
    G[*neG] = 2*x[0]-6;

    // iGfun[*neG] = 2
    // jGvar[*neG] = 2
    *neG = *neG + 1;
    G[*neG] = 2*x[1]-12;

    // iGfun[*neG] = 2
    // jGvar[*neG] = 3
    *neG = *neG + 1;
    G[*neG] = 2*x[2];
    *neG = *neG + 1;  */
  }




  return 0;
}
