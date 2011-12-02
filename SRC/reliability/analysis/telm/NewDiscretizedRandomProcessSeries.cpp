/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.4 $
// $Date: 2010-02-04 18:31:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewDiscretizedRandomProcessSeries.cpp,v $

#include <NewDiscretizedRandomProcessSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <ModulatingFunction.h>
#include <Filter.h>
#include <classTags.h>
#include <Parameter.h>

NewDiscretizedRandomProcessSeries::NewDiscretizedRandomProcessSeries(int numMod, 
								     ModulatingFunction **theModFuncs,
								     double p_mean,
								     double p_maxStdv)
:TimeSeries(TSERIES_TAG_DiscretizedRandomProcessSeries)
{
  numRandVar=0;
  randomVariables = 0;
  kickInTimes = 0;
  theModulatingFunctions = theModFuncs;
  numModFuncs = numMod;
  mean = p_mean;
  maxStdv = p_maxStdv;
  parameterID = -1;
  c = maxStdv;
  active = 0;
  
  output.open("excitation.txt", ios::out);
}


NewDiscretizedRandomProcessSeries::~NewDiscretizedRandomProcessSeries()
{
  if (randomVariables != 0)  {
    delete randomVariables;
    randomVariables = 0;
  }
  
  if (kickInTimes != 0) {
    delete kickInTimes;
    kickInTimes = 0;
  }

  if (active !=0 ) {
    delete [] active;
    active = 0;
  }
  if ( theModulatingFunctions != 0 ){
    delete [] theModulatingFunctions;
    theModulatingFunctions = 0;
  }
}

TimeSeries *
NewDiscretizedRandomProcessSeries::getCopy(void)
{
  opserr << "NewDiscretizedRandomProcessSeries::getCopy() - not yet implemented\n";
  return 0;
}

int
NewDiscretizedRandomProcessSeries::setParameter (const char **argv, int argc, Parameter &param)
{
  if (argc < 2)
    return -1;
 
  if (strcmp(argv[0], "-randomProcessDiscretizer") != 0)
    return -1;
 
  int i;
  
  if (param.addObject(numRandVar, this) < 0) {
    opserr << "NewDiscretizedRandomProcessSeries::setParameter() - failed to add this to parameter\n";
    return -1;
  }
  
  int rvNumber = numRandVar;
  
  // The second argument tells when the random variable "kicks in".
  // Store this in a table...
  // In case the vector doesn't exist
  numRandVar++;

  if( kickInTimes != 0 ){
    delete kickInTimes;
    kickInTimes=0;
  }

  kickInTimes = new Vector(numRandVar);
  //  double ktime=(double)atof(argv[1])+1.0e-8;
  (*kickInTimes)(rvNumber)=(double)atof(argv[1])+1.0e-8;
  
  if( active != 0 ) {
    delete [] active;
    active = 0;
  }
  active = new bool[numRandVar];
  active[rvNumber]=false;
  
  if( randomVariables != 0 ){
    delete randomVariables;
    randomVariables = 0;
  }
  randomVariables = new Vector(numRandVar);
  for( i=0; i<numRandVar; i++) (*randomVariables)(i)=0.0; 
  
  // The random variable number is returned as a parameter ID
  return rvNumber;
}

double
NewDiscretizedRandomProcessSeries::getFactor(double time)
{
  if (time == 0.0) {
    return 0.0;
  }
  else if (randomVariables == 0 || kickInTimes == 0) {
    opserr << "ERROR in DiscretizedRandomProcessSeries::getFactor(): " << endln
	   << " random variables or kick-in times vector(s) do not exist. " << endln;
    return 0.0;
  }
  else if (kickInTimes->Size() != randomVariables->Size() ) {
    opserr << "ERROR in DiscretizedRandomProcessSeries::getFactor(): " << endln
	   << " number of random variables is not the same as kick-in times. " << endln;
    return 0.0;
  }
  else {
    double sum1;
    double sum2;
    //		int nrv = 0;
    double modFuncAmplitude, filterAmplitude;
    Filter *theFilter;
    
    //		nrv = randomVariables->Size();
    for (int i=0; i<numRandVar; i++) {
      active[i]=false;
      if(time>(*kickInTimes)(i)-1.0e-7){
	active[i]=true;
      }
    }
    
    // Loop over all modulating functions
    sum1 = 0.0;
    double dtime;
    for (int k=0; k<numModFuncs; k++) {
      
      // Get value of modulating function number k at time t
      modFuncAmplitude = theModulatingFunctions[k]->getAmplitude(time);
      theFilter = theModulatingFunctions[k]->getFilter();
      
      // Number of discretizing random variables
      //			nrv = randomVariables->Size();
      
      // Loop over all active rv's 
      sum2 = 0.0;
      for (int i=0; i<numRandVar; i++) {
	
	// Get value of filter for argument (t-ti)
	dtime=time-(*kickInTimes)(i);
	theFilter->setKickTime((*kickInTimes)(i));
	filterAmplitude = theFilter->getAmplitude(dtime, 0.0);
	
	// Add contribution 'ui * hi'
	sum2 += (*randomVariables)(i) * filterAmplitude;
	
	// Break when we get to inactive rv's
	if (dtime <= -1.0e-7) {
	  break;
	}
      }
      
      sum1 += sum2*modFuncAmplitude;
    }
    
    double result = mean + c*sum1;
    
    output << "time... ," << time << ","<< "value.. ," << sum1 << "\n"; 
    output.flush();
    
    return result;
  }
}


double
NewDiscretizedRandomProcessSeries::getFactorSensitivity(double time)
{
	// The parameterID has been set to the number of 
	// the random variable in question

	// So, do the same thing as above, just set x(i-1) equal to 1.0
	// for i==parameterID

	if (time == 0.0 || parameterID<0) {
		return 0.0;
	}
	else if (randomVariables == 0 || kickInTimes == 0) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " random variables or kick-in times vector(s) do not exist. " << endln;
		return 0.0;
	}
	else if (kickInTimes->Size() != randomVariables->Size() ) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " number of random variables is not the same as kick-in times. " << endln;
		return 0.0;
	}
	else {

		double sum1;
		double sum2;
//		int nrv = 0;
		double modFuncAmplitude;
		Filter *theFilter;

		// Loop over all modulating functions
		double dtime;
		sum1 = 0.0;
		for (int k=0; k<numModFuncs; k++) {

			// Get value of modulating function number k at time t
			modFuncAmplitude = theModulatingFunctions[k]->getAmplitude(time);
			theFilter = theModulatingFunctions[k]->getFilter();

			// Number of discretizing random variables
			//nrv = randomVariables->Size();

			// Loop over all rv's (even though some may be zero at this time)
			dtime=time-(*kickInTimes)(parameterID);
//			if(fabs(dtime)<=1.0e-7) dtime=0.0;
			sum2 = theFilter->getAmplitude(dtime, 0.0);
			sum1 += sum2*modFuncAmplitude;
		}

//		double result=0.0;
//		if(fabs(time-(*kickInTimes)(parameterID-1))<= 1.0e-8) result=1.0;
		double result = mean + c*sum1;
		return result;
	}
}


double
NewDiscretizedRandomProcessSeries::getFactorSensitivity(double time, double ktime)
{
	// The parameterID has been set to the number of 
	// the random variable in question

	// So, do the same thing as above, just set x(i-1) equal to 1.0
	// for i==parameterID

	if (time == 0.0) {
		return 0.0;
	}
	else if (randomVariables == 0 || kickInTimes == 0) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " random variables or kick-in times vector(s) do not exist. " << endln;
		return 0.0;
	}
	else if (kickInTimes->Size() != randomVariables->Size() ) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " number of random variables is not the same as kick-in times. " << endln;
		return 0.0;
	}
	else {

		double sum1;
		double sum2;
//		int nrv = 0;
		double modFuncAmplitude;
		Filter *theFilter;

		// Loop over all modulating functions
		double dtime;
		sum1 = 0.0;
		for (int k=0; k<numModFuncs; k++) {

			// Get value of modulating function number k at time t
			modFuncAmplitude = theModulatingFunctions[k]->getAmplitude(time);
			theFilter = theModulatingFunctions[k]->getFilter();

			// Number of discretizing random variables
			//nrv = randomVariables->Size();

			// Loop over all rv's (even though some may be zero at this time)
			dtime=time-ktime;
			sum2 = theFilter->getAmplitude(dtime);
			sum1 += sum2*modFuncAmplitude;
		}

//		double result=0.0;
//		if(fabs(time-(*kickInTimes)(parameterID-1))<= 1.0e-8) result=1.0;
		double result = mean + c*sum1;
		return result;
	}
}

int
NewDiscretizedRandomProcessSeries::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}


int 
NewDiscretizedRandomProcessSeries::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;    
}


void
NewDiscretizedRandomProcessSeries::Print(OPS_Stream &s, int flag)
{
}


int
NewDiscretizedRandomProcessSeries::updateParameter (int parameterID, Information &info)
{
  // In case the vector doesn't exist
  if (parameterID < 0 || parameterID >= numRandVar) {
    opserr << "DiscretizedRandomProcessSeries::updateParameter - out of range error\n";
    return -1;
  }

  (*randomVariables)(parameterID) = info.theDouble;
  
  return 0;
}

int
NewDiscretizedRandomProcessSeries::activateParameter(int passedParameterID)
{
  opserr << "ERROR FMK NewDiscretizedRandomProcessSeries::activateParameter(int passedParameterID) - " << passedParameterID << endln; 
  return -1;

  /*
  if( passedParameterID == 0 ){
    //      ------ inactivate ---------
    parameterID=-1;
    for( int i=0; i<numRandVar; i++) active[i]=false; 
    return 0;
  }else if( passedParameterID < 0 ){
    //      ------ inactivate ---------
    parameterID=-1;
    int itemp=arrayID[-passedParameterID];
    active[itemp]=false; 
    return 0;
  }else{
    //      ------ activate ---------
    parameterID=arrayID[passedParameterID];
    if(active[parameterID]) return 1;
    else return 0;
  }
  */
}

double NewDiscretizedRandomProcessSeries::getkickInTimes(int rvNum) const
{ 
  if( rvNum < 0 && rvNum < numRandVar){
    return (*kickInTimes)(rvNum);
  }
}
int NewDiscretizedRandomProcessSeries::getPulseSequentialID( int rvNum ) const
{
  return rvNum;
}


int
NewDiscretizedRandomProcessSeries::updateRV (int nrv, double value)
{
  if (nrv < 0 || nrv >= numRandVar) {
    opserr << "DiscretizedRandomProcessSeries::updateParameter - out of range error\n";
    return -1;
  }

  (*randomVariables)(nrv) = value;
  
  return 0;
}

