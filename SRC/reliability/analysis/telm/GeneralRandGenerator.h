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
                                                                        
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/GeneralRandGenerator.h,v $

#ifndef GeneralRandGenerator_h
#define GeneralRandGenerator_h

#include <RandomNumberGenerator.h>
#include <randomc.h>
#include <NormalRV.h>

class GeneralRandGenerator : public RandomNumberGenerator
{

public:
	GeneralRandGenerator(int passedType=1,int passedSeed=0);
	~GeneralRandGenerator();

	int		generate_nIndependentStdNormalNumbers(int n, int seed=0);
	int     generate_nIndependentUniformNumbers(int n, double lower, double upper, int seed=0);
	const Vector	&getGeneratedNumbers();
	int     getSeed();
	double  generate_singleStdNormalNumber();		
	double  generate_singleUniformNumber(double lower=0.0, double upper=1.0);		
	void setSeed(int passedSeed=0);

protected:

private:
	UniformGenerator* theUniformGenerator;
	Vector *generatedNumbers;
	int seed;
	int	randomNumberBetween0And32767;
	double randomNumberBetween0And1;
	double randomNumber;

	int Type;


};

#endif

