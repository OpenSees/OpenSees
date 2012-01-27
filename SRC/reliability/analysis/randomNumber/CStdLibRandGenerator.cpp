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
#include <math.h>
#include <stdlib.h>
#include <time.h>


CStdLibRandGenerator::CStdLibRandGenerator()
:RandomNumberGenerator()
{
	generatedNumbers = 0;
    setSeed(0);
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
	int randomNumberBetween0AndRAND_MAX;
	double randomNumberBetween0And1;
    
    // set RNG seed if necessary
    if (seedIn != 0)
		setSeed(seedIn);

    // size output vector
	if (generatedNumbers == 0) {
		generatedNumbers = new Vector(n);
	}
	else if (generatedNumbers->Size() != n) {
		delete generatedNumbers;
		generatedNumbers = new Vector(n);
	}
	Vector &randomArray = *generatedNumbers;

	// Create array of standard normal random numbers
	for ( int j=0; j<n; j++) {
		// Generate a number between 0 and RAND_MAX
		randomNumberBetween0AndRAND_MAX = rand();

		// Modify it so that the value lies between 0 and 1
		randomNumberBetween0And1 = (double)randomNumberBetween0AndRAND_MAX/RAND_MAX;

		// Transform according to uniform distribution
		randomArray(j) = (upper-lower)*randomNumberBetween0And1 + lower;
	}

    // KRM - not sure what this is meant to do, should probably use setSeed(0)
	seed = randomNumberBetween0AndRAND_MAX;
	
	return 0;
}


int
CStdLibRandGenerator::generate_nIndependentStdNormalNumbers(int n, int seedIn)
{
	// Initial declarations
	int randomNumberBetween0AndRAND_MAX;
	double randomNumberBetween0And1;

    // set RNG seed if necessary
    if (seedIn != 0)
		setSeed(seedIn);
    
    // size output vector
	if (generatedNumbers == 0) {
		generatedNumbers = new Vector(n);
	}
	else if (generatedNumbers->Size() != n) {
		delete generatedNumbers;
		generatedNumbers = new Vector(n);
	}
	Vector &randomArray = *generatedNumbers;
    
	// Create array of standard normal random numbers
    static NormalRV uRV(1, 0.0, 1.0);
	for ( int j=0; j<n; j++) {
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
		//    where F(x) for the uniform distribution from 0 to 1 in fact is equal to x itself.
		randomArray(j) = uRV.getInverseCDFvalue(randomNumberBetween0And1); 
	}
    
    // KRM - not sure what this is meant to do, should probably use setSeed(0)
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
	if (passedSeed!=0) {
		srand(passedSeed);
		seed = passedSeed;
	} else {
		seed = time(NULL);
		srand(seed);
	}
}


double 
CStdLibRandGenerator::generate_singleUniformNumber(double lower, double upper)
{
    // this is silly, shouldn't repeat code with _single and _n functions - KRM
    // modifying to use previous methods
    generate_nIndependentUniformNumbers(1, lower, upper);
    return (*generatedNumbers)(0);
}


double
CStdLibRandGenerator::generate_singleStdNormalNumber(void)
{
    // this is silly, shouldn't repeat code with _single and _n functions - KRM
    // modifying to use previous methods
    generate_nIndependentStdNormalNumbers(1);
    return (*generatedNumbers)(0);
}

