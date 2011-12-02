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
                                                                        
// $Revision: 1.5 $
// $Date: 2005-11-29 22:42:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/BFGS.cpp,v $
                                                                        
// Written: Ed Love
// Created: 06/01

// What: "@(#)BFGS.cpp, revA"

#include <BFGS.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>

// Constructor
BFGS::BFGS(int theTangentToUse, int n )
:EquiSolnAlgo(EquiALGORITHM_TAGS_BFGS),
 theTest(0), tangent(theTangentToUse), numberLoops(n) 
{
  s  = new Vector*[numberLoops+3];

  z  = new Vector*[numberLoops+3];

  //  r  = new (Vector*)[numberLoops+3];

  residOld = 0;
  residNew = 0;
  du = 0;
  b  = 0;

  temp = 0;

  rdotz = 0;
  sdotr = 0;

  for ( int i =0; i < numberLoops+3; i++ ) {
    s[i] = 0;
    z[i] = 0;
  }
 
  localTest = 0;

}

//Constructor
BFGS::BFGS(ConvergenceTest &theT, int theTangentToUse, int n)
:EquiSolnAlgo(EquiALGORITHM_TAGS_BFGS),
 theTest(&theT), tangent(theTangentToUse), numberLoops(n) 
{
  s  = new Vector*[numberLoops+3];

  z  = new Vector*[numberLoops+3];

  residOld = 0;
  residNew = 0;
  du = 0;
  b  = 0;
  temp = 0;

  rdotz = 0;
  sdotr = 0;

  for ( int i =0; i < numberLoops+3; i++ ) {
    s[i] = 0;
    z[i] = 0;
    //r[i] = 0;
  }

  localTest = theTest->getCopy( this->numberLoops );

}

// Destructor
BFGS::~BFGS()
{

  if (temp != 0) delete temp;
  temp = 0;

  if (residOld != 0 ) delete residOld;  
  residOld = 0;

  if (residNew != 0) delete residNew;
  residNew = 0;

  if (du != 0) delete du;
  du = 0;

  if (b != 0 ) delete b;
  b = 0;

  if (rdotz != 0 ) delete [] rdotz;
  rdotz = 0;
  
  if (sdotr != 0 ) delete [] sdotr;
  sdotr = 0;

  for ( int i =0; i < numberLoops+3; i++ ) {
    if ( s[i] != 0 ) delete s[i];
    if ( z[i] != 0 ) delete z[i];
    //delete r[i];
    s[i] = 0;
    z[i] = 0;
    //r[i] = 0;
  } //end for i

  if ( s != 0 ) delete[] s; 
  if ( z != 0 ) delete[] z;
  s = 0;  
  z = 0;

  if ( localTest != 0 )
     delete localTest;
  localTest = 0;

  
}


int
BFGS::setConvergenceTest(ConvergenceTest *newTest)
{
    theTest = newTest;

    if ( localTest != 0 )  
      delete localTest;

    localTest = theTest->getCopy( this->numberLoops );
    if (localTest != 0) {
      opserr << "BFGS::setConvergenceTest() - could not get copy for local test\n";
      return -1;
    } else
      return 0;
}



int 
BFGS::solveCurrentStep(void)
{
 
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass

    AnalysisModel   *theAnaModel = this->getAnalysisModelPtr();

    IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();

    LinearSOE  *theSOE = this->getLinearSOEptr();

    if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
	|| (theTest == 0)){
	opserr << "WARNING BFGS::solveCurrentStep() - setLinks() has";
	opserr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "BFGS::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    localTest->setEquiSolnAlgo(*this);

    if (rdotz == 0)
       rdotz = new double[numberLoops+3];

    if (sdotr == 0)
	sdotr = new double[numberLoops+3];


    int result = -1;
    int count = 0;
    do {

      // opserr << "      BFGS -- Forming New Tangent" << endln;

      //form the initial tangent
      if (theIntegrator->formTangent(tangent) < 0){
         opserr << "WARNING BFGS::solveCurrentStep() -";
         opserr << "the Integrator failed in formTangent()\n";
         return -1; 
      }

      //form the initial residual 
      if (theIntegrator->formUnbalance() < 0) {
        opserr << "WARNING BFGS::solveCurrentStep() -";
        opserr << "the Integrator failed in formUnbalance()\n";	
      }	    

      //solve
      if (theSOE->solve() < 0) {
	  opserr << "WARNING BFGS::solveCurrentStep() -";
	  opserr << "the LinearSysOfEqn failed in solve()\n";	
	  return -3;
	}	    

      //update
      if ( theIntegrator->update(theSOE->getX() ) < 0) {
	opserr << "WARNING BFGS::solveCurrentStep() -";
	opserr << "the Integrator failed in update()\n";	
	return -4;
      }	        


      //    int systemSize = ( theSOE->getB() ).Size();
      int systemSize = theSOE->getNumEqn( );

      //temporary vector
      if (temp == 0 )
	temp = new Vector(systemSize);

      //initial displacement increment
      if ( s[1] == 0 ) 
	s[1] = new Vector(systemSize);

      *s[1] = theSOE->getX( );

      if ( residOld == 0 ) 
	residOld = new Vector(systemSize);

      *residOld = theSOE->getB( ) ;
      *residOld *= (-1.0 );

      //form the residual again
      if (theIntegrator->formUnbalance() < 0) {
        opserr << "WARNING BFGS::solveCurrentStep() -";
        opserr << "the Integrator failed in formUnbalance()\n";	
      }	    

      if ( residNew == 0 ) 
	residNew = new Vector(systemSize);
 
      if ( du == 0 ) 
	du = new Vector(systemSize);

      if ( b == 0 )
	b = new Vector(systemSize);

      localTest->start();

      int nBFGS = 1;
      do {

        //save residual
        *residNew =  theSOE->getB( ); 
        *residNew *= (-1.0 );

      
        //solve
        if (theSOE->solve() < 0) {
	    opserr << "WARNING BFGS::solveCurrentStep() -";
	    opserr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
        }	    

	//save right hand side
        *b = theSOE->getB( );

        //save displacement increment
        *du = theSOE->getX( );

        //BFGS modifications to du
        BFGSUpdate( theIntegrator, theSOE, *du, *b, nBFGS ) ;

        if ( theIntegrator->update( *du ) < 0 ) {
	   opserr << "WARNING BFGS::solveCurrentStep() -";
	   opserr << "the Integrator failed in update()\n";	
	   return -4;
        }	        

	/* opserr << "        BFGS Iteration " << nBFGS 
            << " Residual Norm = " 
            << sqrt( (*residNew) ^ (*residNew) ) << endln;
	*/
        
        //increment broyden counter
        nBFGS += 1;

        //save displacement increment
        if ( s[nBFGS] == 0 ) 
	  s[nBFGS] = new Vector(systemSize);

        *s[nBFGS] = *du;

        //swap residuals
	*residOld = *residNew;

        //form the residual again
        if (theIntegrator->formUnbalance() < 0) {
          opserr << "WARNING BFGS::solveCurrentStep() -";
          opserr << "the Integrator failed in formUnbalance()\n";	
        }	    

        result = localTest->test();
 
        
      } while ( result == -1 && nBFGS <= numberLoops );


      result = theTest->test();
      this->record(count++);

    }  while (result == -1);


    if (result == -2) {
      opserr << "BFGS::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }

    // note - if postive result we are returning what the convergence test returned
    // which should be the number of iterations
    return result;
}



void  BFGS::BFGSUpdate(IncrementalIntegrator *theIntegrator, 
		       LinearSOE *theSOE, 
		       Vector &du, 
		       Vector &b,
		       int nBFGS) 
{

  static const double eps = 1.0e-16;

  //  int systemSize = ( theSOE->getB() ).Size();
      int systemSize = theSOE->getNumEqn( );


  //compute z
  //  theSOE->setB( (*r[nBFGS]) - (*r[nBFGS-1]) );
  //    theSOE->setB( (*residNew) - (*residOld) );
  *temp = *residNew;
  *temp -= *residOld;
  theSOE->setB(*temp);


  if (theSOE->solve() < 0) {
       opserr << "WARNING BFGS::solveCurrentStep() -";
       opserr << "the LinearSysOfEqn failed in solve()\n";	
   }	    
  
  if ( z[nBFGS] == 0 ) 
    z[nBFGS] = new Vector(systemSize);

  *z[nBFGS] = theSOE->getX(); 
  //  *z[nBFGS] *= (-1.0);

    int i;
  for ( i=1; i<=(nBFGS-1); i++ ) {

    if ( sdotr[i] < eps ) 
      break; 

    double fact1 = 1.0 + ( rdotz[i] / sdotr[i] );

    fact1 /= sdotr[i];

    double pdotb = (*s[i]) ^ ( theSOE->getB() );


    fact1 *= pdotb;

    //    *z[nBFGS] +=  fact1 * ( *s[i] );
    *temp = *s[i];
    *temp *= fact1;
    *z[nBFGS] += *temp;


    double bdotz = (*z[i]) ^ ( theSOE->getB() );  

    //    *z[nBFGS] -= (1.0/sdotr[i]) * 
    //             ( bdotz * (*s[i])   +  pdotb * (*z[i]) );   
    *temp = *s[i];
    *temp *= bdotz;
    *temp /= sdotr[i];
    *z[nBFGS] -= *temp;

    *temp = *z[i];
    *temp *= pdotb;
    *temp /= sdotr[i];
    *z[nBFGS] -= *temp;
 
  } //end for i


  //sdotr[nBFGS] = *s[nBFGS] ^ ( *residNew - *residOld );

  //rdotz[nBFGS] = *z[nBFGS] ^ ( *residNew - *residOld );   

  *temp = *residNew;
  *temp -= *residOld;

  sdotr[nBFGS] = *s[nBFGS] ^ (*temp);

  rdotz[nBFGS] = *z[nBFGS] ^ (*temp);


  //BFGS modifications to du
  for ( i=1; i<=nBFGS; i++ ) {

    if ( sdotr[i] < eps )
      break;

    double fact1 = 1.0 + ( rdotz[i] / sdotr[i] );

    fact1 /= sdotr[i];

    double sdotb = (*s[i]) ^ b;

    fact1 *= sdotb;

    //du +=  fact1 * ( *s[i] );
    *temp = *s[i];
    *temp *= fact1;
    du += *temp;


    double bdotz = (*z[i]) ^ b;  

    //du -= (1.0/sdotr[i]) * 
    //             ( bdotz * (*s[i])   +  sdotb * (*z[i]) );   
    *temp = *s[i];
    *temp *= bdotz;
    *temp /= sdotr[i];
    du -= *temp;

    *temp = *z[i];
    *temp *= sdotb;
    *temp /= sdotr[i];
    du -= *temp;


  } //end for i

}


ConvergenceTest *
BFGS::getConvergenceTest(void)
{
  return theTest;
}

int
BFGS::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
BFGS::recvSelf(int cTag, 
	       Channel &theChannel, 
	       FEM_ObjectBroker &theBroker)
{
    return -1;
}


void
BFGS::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) {
      s << "BFGS" << endln;
      s << "  Number of Iterations = " << numberLoops << endln;
    }
}
