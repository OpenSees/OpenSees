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
                                                                        
// $Revision: 1.7 $
// $Date: 2010-04-22 23:42:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/Broyden.cpp,v $
                                                                        
                                                                        
// Written: Ed C++ Love
// Created: 04/01

// Description: This file contains the class definition implementation of
// Broyden.  
// 
// What: "@(#)Broyden.h, revA"


#include <Broyden.h>
#include <AnalysisModel.h>
#include <StaticAnalysis.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>
#include <math.h>
#include <elementAPI.h>

void* OPS_Broyden()
{
    int formTangent = CURRENT_TANGENT;
    int count = -1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();
	
	if (strcmp(flag,"-secant") == 0) {
	    formTangent = CURRENT_SECANT;
	    
	} else if (strcmp(flag,"-initial") == 0) {
	    formTangent = INITIAL_TANGENT;
	    
	} else if (strcmp(flag,"-count") == 0 && OPS_GetNumRemainingInputArgs()>0) {
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &count) < 0) {
		opserr << "WARNING Broyden failed to read count\n";
		return 0;
	    }
	}
    }

    if (count == -1)
	return new Broyden(formTangent); 
    else
	return new Broyden(formTangent, count); 
}

// Constructor
Broyden::Broyden(int theTangentToUse, int n )
:EquiSolnAlgo(EquiALGORITHM_TAGS_Broyden),
 tangent(theTangentToUse), numberLoops(n) 
{
  s  = new Vector*[numberLoops+3] ;

  z  = new Vector*[numberLoops+3] ;

  residOld = 0 ;
  residNew = 0 ;
  du = 0 ;
  temp = 0 ;

  for ( int i =0; i < numberLoops+3; i++ ) {
    s[i] = 0 ;
    z[i] = 0 ;
    //r[i] = 0 ;
  }

  localTest = 0;
}

//Constructor
Broyden::Broyden(ConvergenceTest &theT, int theTangentToUse, int n)
:EquiSolnAlgo(EquiALGORITHM_TAGS_Broyden),
 tangent(theTangentToUse), numberLoops(n) 
{
  s  = new Vector*[numberLoops+3] ;
  z  = new Vector*[numberLoops+3] ;

  residOld = 0 ;
  residNew = 0 ;
  du = 0 ;
  temp = 0 ;

  for ( int i =0; i < numberLoops+3; i++ ) {
    s[i] = 0 ;
    z[i] = 0 ;
  }

  localTest= 0;
}

// Destructor
Broyden::~Broyden()
{
  if ( residOld != 0 ) delete residOld ;  
  residOld = 0 ;

  if ( residNew != 0 ) delete residNew ;
  residNew = 0 ;

  if ( du != 0 ) delete du ;
  du = 0 ;

  if ( temp != 0 ) delete temp ;
  temp = 0 ;

  for ( int i =0; i < numberLoops+3; i++ ) {
    if ( s[i] != 0 ) delete s[i] ;
    if ( z[i] != 0 ) delete z[i] ;
    s[i] = 0 ;
    z[i] = 0 ;
  } //end for i

  if ( s != 0 ) delete[] s ; 

  if ( z != 0 ) delete[] z ;

  s = 0 ;  
  z = 0 ;

  if ( localTest != 0 )
     delete localTest ;
  localTest = 0;
}


void 
Broyden::setLinks(AnalysisModel &theModel, 
		  IncrementalIntegrator &theIntegrator,
		  LinearSOE &theSOE,
		  ConvergenceTest *theTest)
{
  this->EquiSolnAlgo::setLinks(theModel, theIntegrator, theSOE, theTest);

  if (theTest == 0)
    return;

  if ( localTest != 0 ) 
    delete localTest ;

  localTest = theTest->getCopy(this->numberLoops);
  if (localTest == 0) {
    opserr << "Broyden::setTest() - could not get a copy\n";
  } 
}


int
Broyden::setConvergenceTest(ConvergenceTest *newTest)
{
  this->EquiSolnAlgo::setConvergenceTest(newTest);

  if (theTest == 0)
    return 0;
  
  if ( localTest != 0 )  
    delete localTest ;
  
  localTest = theTest->getCopy( this->numberLoops ) ;
  if (localTest == 0) {
    opserr << "Broyden::setTest() - could not get a copy\n";
    return -1;
  } 
  
  return 0;
}



int 
Broyden::solveCurrentStep(void)
{
 
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel   *theAnaModel = this->getAnalysisModelPtr();

    IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();

    LinearSOE  *theSOE = this->getLinearSOEptr();


    if ((theAnaModel == 0) || (theIntegrator == 0) || (theSOE == 0)
	|| (theTest == 0)){
	opserr << "WARNING Broyden::solveCurrentStep() - setLinks() has";
	opserr << " not been called - or no ConvergenceTest has been set\n";
	return -5;
    }	

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "Broyden::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return -3;
    }

    localTest->setEquiSolnAlgo(*this);

    int result = -1 ;
    int count = 0 ;
    do {

      // opserr << "      Broyden -- Forming New Tangent" << endln ;

      //form the initial tangent
      if (theIntegrator->formTangent(tangent) < 0){
         opserr << "WARNING Broyden::solveCurrentStep() -";
         opserr << "the Integrator failed in formTangent()\n";
         return -1; 
      }

      //form the initial residual 
      if (theIntegrator->formUnbalance() < 0) {
        opserr << "WARNING Broyden::solveCurrentStep() -";
        opserr << "the Integrator failed in formUnbalance()\n";	
      }	    

      //solve
      if (theSOE->solve() < 0) {
	  opserr << "WARNING Broyden::solveCurrentStep() -";
	  opserr << "the LinearSysOfEqn failed in solve()\n";	
	  return -3;
	}	    

      //update
      if ( theIntegrator->update(theSOE->getX() ) < 0) {
	opserr << "WARNING Broyden::solveCurrentStep() -";
	opserr << "the Integrator failed in update()\n";	
	return -4;
      }	        

      //    int systemSize = ( theSOE->getB() ).Size();
      int systemSize = theSOE->getNumEqn( ) ;

      //temporary vector
      if ( temp == 0 ) 
	temp = new Vector(systemSize) ;

      //initial displacement increment
      if ( s[1] == 0 ) 
	s[1] = new Vector(systemSize) ;

      *s[1] = theSOE->getX( ) ;

      //initial residual
      /* if ( r[0] == 0 ) r[0] = new Vector(systemSize) ;
      *r[0] = theSOE->getB( )  ;
      *r[0] *= (-1.0 ) ;
      */

      if ( residOld == 0 ) 
	residOld = new Vector(systemSize) ;

      *residOld = theSOE->getB( )  ;
      *residOld *= (-1.0 ) ;

      //form the residual again
      if (theIntegrator->formUnbalance() < 0) {
        opserr << "WARNING Broyden::solveCurrentStep() -";
        opserr << "the Integrator failed in formUnbalance()\n";	
      }	    

      if ( residNew == 0 ) 
	residNew = new Vector(systemSize) ;
 
      if ( du == 0 ) 
	du = new Vector(systemSize) ;


      localTest->start() ;

      int nBroyden = 1 ;
      do {

        //save residual
	/*        if ( r[nBroyden] == 0 ) r[nBroyden] = new Vector(systemSize) ;
        *r[nBroyden] =  theSOE->getB( ) ; 
        *r[nBroyden] *= (-1.0 ) ; 
        */

        *residNew =  theSOE->getB( ) ; 
        *residNew *= (-1.0 ) ;
      
        //solve
        if (theSOE->solve() < 0) {
	    opserr << "WARNING Broyden::solveCurrentStep() -";
	    opserr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
        }	    

        //save displacement increment
        *du = theSOE->getX( ) ;

        //broyden modifications to du
        BroydenUpdate( theIntegrator, theSOE, *du, nBroyden )  ;

        if ( theIntegrator->update( *du ) < 0 ) {
	   opserr << "WARNING Broyden::solveCurrentStep() -";
	   opserr << "the Integrator failed in update()\n";	
	   return -4;
        }	        


	/*	opserr << "        Broyden Iteration " << nBroyden 
            << " Residual Norm = " 
            << sqrt( (*residNew) ^ (*residNew) ) << endln ;
	*/
        
        //increment broyden counter
        nBroyden += 1 ;

        //save displacement increment
        if ( s[nBroyden] == 0 ) 
	  s[nBroyden] = new Vector(systemSize) ;

        *s[nBroyden] = *du ;

        //swap residuals
	*residOld = *residNew ;

        //form the residual again
        if (theIntegrator->formUnbalance() < 0) {
          opserr << "WARNING Broyden::solveCurrentStep() -";
          opserr << "the Integrator failed in formUnbalance()\n";	
        }	    
	
	result = localTest->test() ;
        
      } while ( result == -1 && nBroyden <= numberLoops );


      result = theTest->test();
      this->record(count++);

    }  while (result == -1);


    if (result == -2) {
      opserr << "Broyden::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in test()\n";
      return -3;
    }


    // note - if positive result we are returning what the convergence test returned
    // which should be the number of iterations
    return result;
}



void  Broyden::BroydenUpdate( IncrementalIntegrator *theIntegrator, 
			       LinearSOE *theSOE, 
                               Vector &du, 
			       int nBroyden ) 
{

  static const double eps = 1.0e-16 ;

  //  int systemSize = ( theSOE->getB() ).Size();
      int systemSize = theSOE->getNumEqn( ) ;


  //compute z
  //  theSOE->setB( (*r[nBroyden]) - (*r[nBroyden-1]) ) ;
  //    theSOE->setB( (*residNew) - (*residOld) ) ;
  *temp  = (*residNew) ;
  *temp -= (*residOld) ;
  theSOE->setB( *temp ) ;

  if (theSOE->solve() < 0) {
       opserr << "WARNING Broyden::solveCurrentStep() -";
       opserr << "the LinearSysOfEqn failed in solve()\n";	
   }	    
  
  if ( z[nBroyden] == 0 ) 
    z[nBroyden] = new Vector(systemSize) ;

  *z[nBroyden] = theSOE->getX() ; 
  *z[nBroyden] *= (-1.0) ;

  int i;
  for ( i=1; i<=(nBroyden-1); i++ ) {

    double p = - ( (*s[i]) ^ (*z[i]) ) ;

    if ( fabs(p) < eps ) break ;

    double sdotz = (*s[i]) ^ (*z[nBroyden]) ;

    //*z[nBroyden] += (1.0/p) * sdotz * ( *s[i] + *z[i] ) ;
    *temp  = (*s[i]) ;
    *temp += (*z[i]) ;
    *temp *= ( (1.0/p) * sdotz ) ;
    *z[nBroyden] += (*temp) ;


  } //end for i


  //broyden modifications to du
  for ( i=1; i<=nBroyden; i++ ) {

    double p = - ( (*s[i]) ^ (*z[i]) ) ;

    if ( fabs(p) < eps ) break ;

    double sdotdu = (*s[i]) ^ du ;

    //du += (1.0/p) * sdotdu * ( *s[i] + *z[i] ) ;
    *temp  = (*s[i]) ;
    *temp += (*z[i]) ;
    *temp *= ( (1.0/p) * sdotdu ) ;
    du += (*temp) ;


  } //end for i

}


ConvergenceTest *
Broyden::getConvergenceTest(void)
{
  return theTest;
}

int
Broyden::sendSelf(int cTag, Channel &theChannel)
{
  static ID data(2);
  data(0) = tangent;
  data(1) = numberLoops;
  if (theChannel.sendID(0, cTag, data) < 0) {
    opserr << "Broyden::sendSelf() - failed to send data\n";
    return -1;
  }
  return 0;
}

int
Broyden::recvSelf(int cTag, 
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
  static ID data(2);
  if (theChannel.recvID(0, cTag, data) < 0) {
    opserr << "Broyden::recvSelf() - failed to recv data\n";
    return -1;
  }
  tangent = data(0);

  if (numberLoops != data(1)) {

    // remove old
    if (s != 0 && z != 0) {
      for ( int i =0; i < numberLoops+3; i++ ) {
	if ( s[i] != 0 ) delete s[i] ;
	if ( z[i] != 0 ) delete z[i] ;
      }
      delete [] s;
      delete [] z;
    }

    numberLoops = data(1);

    // create new
    s  = new Vector*[numberLoops+3] ;
    z  = new Vector*[numberLoops+3] ;
    for ( int i =0; i < numberLoops+3; i++ ) {
      s[i] = 0 ;
      z[i] = 0 ;
    }
  }

  return 0;
}


void
Broyden::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) {
      s << "Broyden" << endln ;
      s << "  Number of Iterations = " << numberLoops << endln ;
    }
}


//const Vector &pi = *p[i] ;

//    //residual at this iteration before next solve 
//    Resid0 = theSOE->getB() ;

//    //line search direction 
//    dx0 = theSOE->getX() ;





/*

        //first solve step

 	if (theIntegrator->formTangent(tangent) < 0){
	    opserr << "WARNING Broyden::solveCurrentStep() -";
	    opserr << "the Integrator failed in formTangent()\n";
	    return -1;
	}		    
	
	if (theSOE->solve() < 0) {
	    opserr << "WARNING Broyden::solveCurrentStep() -";
	    opserr << "the LinearSysOfEqn failed in solve()\n";	
	    return -3;
	}	    


	if (theIntegrator->update(theSOE->getX()) < 0) {
	    opserr << "WARNING Broyden::solveCurrentStep() -";
	    opserr << "the Integrator failed in update()\n";	
	    return -4;
	}	        


	if (theIntegrator->formUnbalance() < 0) {
	    opserr << "WARNING Broyden::solveCurrentStep() -";
	    opserr << "the Integrator failed in formUnbalance()\n";	
	    return -2;
	}	
	
	result = theTest->test();
	this->record(nBroyden++);

      const Vector &du = BroydengetX( theIntegrator, theSOE, nBroyden )  ;

*/
