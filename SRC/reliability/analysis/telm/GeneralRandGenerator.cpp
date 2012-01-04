/*
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
                                                                        
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/GeneralRandGenerator.cpp,v $


#include <RandomNumberGenerator.h>
#include <GeneralRandGenerator.h>
#include <NormalRV.h>
#include <Vector.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

GeneralRandGenerator::GeneralRandGenerator(int passedType, int passedSeed)
:RandomNumberGenerator()
{
	generatedNumbers = 0;

	Type=passedType;

	uint32 seed32;
	if(passedSeed!=0){
		seed=passedSeed;
		seed32=seed;
	}else{
		seed = time(NULL);
		seed32=seed;
	}

	if(Type==1){
	  theUniformGenerator= new TRandomMersenne(seed32);
	}else if(Type==2){
	  theUniformGenerator= new TRanrotWGenerator(seed32);
	}else{
	  theUniformGenerator= new TRandomMotherOfAll(seed32);
	}
	if (theUniformGenerator==0) {
		opserr << "GeneralRandGenerator - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}
}

GeneralRandGenerator::~GeneralRandGenerator()
{
	if (generatedNumbers != 0){
		delete generatedNumbers;
		generatedNumbers=0;
	}
	if (theUniformGenerator != 0){
		delete theUniformGenerator;
		theUniformGenerator=0;
	}
}

int
GeneralRandGenerator::generate_nIndependentUniformNumbers(int n, double lower, double upper, int seedIn)
{
	// Initial declarations
	int j;
	Vector randomArray(n);


	// Create array of standard normal random numbers
	if (seedIn != 0) {
		uint32 seed32=seedIn;
		theUniformGenerator->RandomInit(seed32);
		seed=seedIn;
	}
	for ( j=0; j<n; j++)
	{
		// Modify it so that the value lies between 0 and 1
		randomNumberBetween0And1 = theUniformGenerator->Random();
		// Transform according to uniform distribution
		randomArray(j) = (upper-lower)*randomNumberBetween0And1 + lower;
	}
	if (generatedNumbers == 0) {
		generatedNumbers = new Vector(n);
	}
	else if (generatedNumbers->Size() != n) {
		delete generatedNumbers;
		generatedNumbers=0;
		generatedNumbers = new Vector(n);
	}
	(*generatedNumbers) = randomArray;


	return 0;
}

int
GeneralRandGenerator::generate_nIndependentStdNormalNumbers(int n, int seedIn)
{
	// Initial declarations
	int j;
	Vector randomArray(n);

	static NormalRV aStdNormRV(1, 0.0, 1.0);

	// Create array of standard normal random numbers
	if (seedIn != 0) {
		uint32 seed32=seedIn;
		theUniformGenerator->RandomInit(seed32);
		seed=seedIn;
	}
	for ( j=0; j<n; j++)
	{
		randomNumberBetween0And1 = theUniformGenerator->Random();
		// Treat two special cases
		if (randomNumberBetween0And1 == 0.0) {
			randomNumberBetween0And1 = 1.0e-10;
		}
		// Transform that number into a standard normal variable
		//    Phi(z) = F(x)
		//    z = invPhi( F(x) )
		//       where F(x) for the uniform distribution 
		//       from 0 to 1 in fact is equal to x itself.
		randomArray(j) = aStdNormRV.getInverseCDFvalue(randomNumberBetween0And1); 
	}
	if (generatedNumbers == 0) {
		generatedNumbers = new Vector(n);
	}
	else if (generatedNumbers->Size() != n) {
		delete generatedNumbers;
		generatedNumbers=0;
		generatedNumbers = new Vector(n);
	}
	(*generatedNumbers) = randomArray;

	return 0;
}
const Vector &
GeneralRandGenerator::getGeneratedNumbers()
{
	return (*generatedNumbers);
}
int
GeneralRandGenerator::getSeed()
{
	return seed;
}
void
GeneralRandGenerator::setSeed(int passedSeed)
{
	uint32 seed32;
	if(passedSeed!=0){
		seed=passedSeed;
		seed32=seed;
	}else{
		seed = time(NULL);
		seed32=seed;
	}
	theUniformGenerator->RandomInit(seed32);
}
double 
GeneralRandGenerator::generate_singleUniformNumber(double lower, double upper)
{
	// Initial declarations
	if(seed==0) {
		uint32 seed32;
		seed = time(NULL);
		seed32=seed;
		theUniformGenerator->RandomInit(seed32);
	}

	randomNumberBetween0And1 = theUniformGenerator->Random();
	randomNumber=randomNumberBetween0And1;
	if(lower!=0.0||upper!=1.0) 
		randomNumber = (upper-lower)*randomNumberBetween0And1 + lower;
	return randomNumber;
}
double
GeneralRandGenerator::generate_singleStdNormalNumber(void)
{
	if(seed==0) {
		uint32 seed32;
		seed = time(NULL);
		seed32=seed;
		theUniformGenerator->RandomInit(seed32);
	}

	static NormalRV aStdNormRV(1, 0.0, 1.0);

	randomNumberBetween0And1 = theUniformGenerator->Random();
	// Treat two special cases
	if (randomNumberBetween0And1 == 0.0) {
		randomNumberBetween0And1 = 1.0e-10;
	}
	// Transform that number into a standard normal variable
	//    Phi(z) = F(x)
	//    z = invPhi( F(x) )
	//       where F(x) for the uniform distribution 
	//       from 0 to 1 in fact is equal to x itself.
	randomNumber=aStdNormRV.getInverseCDFvalue(randomNumberBetween0And1); 
	return randomNumber;
}

