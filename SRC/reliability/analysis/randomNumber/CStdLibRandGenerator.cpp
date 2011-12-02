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
                                                                        
// $Revision: 1.11 $
// $Date: 2008-05-12 16:05:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/randomNumber/CStdLibRandGenerator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <RandomNumberGenerator.h>
#include <CStdLibRandGenerator.h>
#include <NormalRV.h>
#include <Vector.h>


CStdLibRandGenerator::CStdLibRandGenerator()
:RandomNumberGenerator()
{
	generatedNumbers = 0;
}

CStdLibRandGenerator::~CStdLibRandGenerator()
{
	if (generatedNumbers != 0)
		delete generatedNumbers;
}





int
CStdLibRandGenerator::generate_nIndependentUniformNumbers(int n, double lower, double upper, int seedIn)
{
	// Initial declarations
	int j;
	int randomNumberBetween0AndRAND_MAX;
	double randomNumberBetween0And1;

	if (generatedNumbers == 0) {
		generatedNumbers = new Vector(n);
	}
	else if (generatedNumbers->Size() != n) {
		delete generatedNumbers;
		generatedNumbers = new Vector(n);
	}
	Vector &randomArray = *generatedNumbers;


	// Create array of standard normal random numbers
	if (seedIn != 0) {
		srand(seedIn);
	}
	for ( j=0; j<n; j++)
	{
		// Generate a number between 0 and RAND_MAX
		randomNumberBetween0AndRAND_MAX = rand();

		// Modify it so that the value lies between 0 and 1
		randomNumberBetween0And1 = (double)randomNumberBetween0AndRAND_MAX/RAND_MAX;

		// Transform according to uniform distribution
		randomArray(j) = (upper-lower)*randomNumberBetween0And1 + lower;
	}

	seed = randomNumberBetween0AndRAND_MAX;
	
	return 0;
}




int
CStdLibRandGenerator::generate_nIndependentStdNormalNumbers(int n, int seedIn)
{
	// Initial declarations
	int j;
	int randomNumberBetween0AndRAND_MAX;
	double randomNumberBetween0And1;
	static NormalRV aStdNormRV(1,0.0,1.0,0.0);

	if (generatedNumbers == 0) {
		generatedNumbers = new Vector(n);
	}
	else if (generatedNumbers->Size() != n) {
		delete generatedNumbers;
		generatedNumbers = new Vector(n);
	}
	Vector &randomArray = *generatedNumbers;

	// Create array of standard normal random numbers
	if (seedIn != 0) {
		srand(seedIn);
	}
	for ( j=0; j<n; j++)
	{
		// Generate a number between 0 and RAND_MAX
		randomNumberBetween0AndRAND_MAX = rand();

		// Modify it so that the value lies between 0 and 1
		randomNumberBetween0And1 = (double)randomNumberBetween0AndRAND_MAX/RAND_MAX;


		// Treat two special cases
		if (randomNumberBetween0And1 == 0.0) {
			randomNumberBetween0And1 = 0.0000001;
		}
		if (randomNumberBetween0And1 == 1.0) {
			randomNumberBetween0And1 = 0.9999999;
		}


		// Transform that number into a standard normal variable
		//    Phi(z) = F(x)
		//    z = invPhi( F(x) )
		//       where F(x) for the uniform distribution 
		//       from 0 to 1 in fact is equal to x itself.
		randomArray(j) = aStdNormRV.getInverseCDFvalue(randomNumberBetween0And1); 
	}
	seed = randomNumberBetween0AndRAND_MAX;

	return 0;
}



const Vector&
CStdLibRandGenerator::getGeneratedNumbers()
{
	return (*generatedNumbers);
}


int
CStdLibRandGenerator::getSeed()
{
	return seed;
}

void
CStdLibRandGenerator::setSeed(int passedSeed)
{
	if(passedSeed!=0){
		srand(passedSeed);
		seed=passedSeed;
	}else{
		seed=time(NULL);
		srand(seed);
	}
}
double 
CStdLibRandGenerator::generate_singleUniformNumber(double lower, double upper)
{
	// Initial declarations
	if(seed==0) {
		seed=time(NULL);
		srand(seed);
	}

	randomNumberBetween0AndRAND_MAX = rand();
	randomNumberBetween0And1 = (double)randomNumberBetween0AndRAND_MAX/RAND_MAX;
	randomNumber=randomNumberBetween0And1;
 	if(lower!=0.0||upper!=1.0) 
 		randomNumber = (upper-lower)*randomNumberBetween0And1 + lower;
 	return randomNumber;
}
double
CStdLibRandGenerator::generate_singleStdNormalNumber(void)
{
	static NormalRV aStdNormRV(1,0.0,1.0,0.0);

	if(seed==0) {
		seed=time(NULL);
		srand(seed);
	}
	randomNumberBetween0AndRAND_MAX = rand();
	// Modify it so that the value lies between 0 and 1
	randomNumberBetween0And1 = (double)randomNumberBetween0AndRAND_MAX/RAND_MAX;
 	// Treat two special cases
 	if (randomNumberBetween0And1 == 0.0) {
 		randomNumberBetween0And1 = 0.0000001;
 	}
 	if (randomNumberBetween0And1 == 1.0) {
 		randomNumberBetween0And1 = 0.9999999;
 	}
 	// Transform that number into a standard normal variable
 	//    Phi(z) = F(x)
 	//    z = invPhi( F(x) )
 	//       where F(x) for the uniform distribution 
 	//       from 0 to 1 in fact is equal to x itself.
 	randomNumber=aStdNormRV.getInverseCDFvalue(randomNumberBetween0And1); 
 	return randomNumber;
}

