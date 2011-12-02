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
                                                                        
// $Revision: 1.7 $
// $Date: 2005-03-21 21:35:56 $
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
	Vector randomArray(n);


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
	
	if (generatedNumbers == 0) {
		generatedNumbers = new Vector(n);
	}
	else if (generatedNumbers->Size() != n) {
		delete generatedNumbers;
		generatedNumbers = new Vector(n);
	}
	(*generatedNumbers) = randomArray;


	return 0;
}




int
CStdLibRandGenerator::generate_nIndependentStdNormalNumbers(int n, int seedIn)
{
	// Initial declarations
	int j;
	int randomNumberBetween0AndRAND_MAX;
	double randomNumberBetween0And1;
	Vector randomArray(n);
	NormalRV *aStdNormRV = 0;
	aStdNormRV = new NormalRV(1,0.0,1.0,0.0);


	// Check if computer ran out of memory
	if (aStdNormRV==0) {
		opserr << "CStdLibRandGenerator::generate_nIndependentStdNormalNumbers() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		return -1;
	}


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
		randomArray(j) = aStdNormRV->getInverseCDFvalue(randomNumberBetween0And1); 
	}
	seed = randomNumberBetween0AndRAND_MAX;

	if (generatedNumbers == 0) {
		generatedNumbers = new Vector(n);
	}
	else if (generatedNumbers->Size() != n) {
		delete generatedNumbers;
		generatedNumbers = new Vector(n);
	}
	(*generatedNumbers) = randomArray;

	delete aStdNormRV;

	return 0;
}



Vector
CStdLibRandGenerator::getGeneratedNumbers()
{
	return (*generatedNumbers);
}


int
CStdLibRandGenerator::getSeed()
{
	return seed;
}
