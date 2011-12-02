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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-02-29 19:47:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/randomNumber/CStdLibRandGenerator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef CStdLibRandGenerator_h
#define CStdLibRandGenerator_h

#include <NormalRV.h>

class CStdLibRandGenerator : public RandomNumberGenerator
{

public:
	CStdLibRandGenerator();
	~CStdLibRandGenerator();

	int		generate_nIndependentStdNormalNumbers(int n, int seed=0);
	int     generate_nIndependentUniformNumbers(int n, double lower, double upper, int seed=0);
	const Vector& getGeneratedNumbers();//added by K.F.
	//const Vector& getGeneratedNumbers();
	int     getSeed();
	/////S added By K Fujimura /////
 	double  generate_singleStdNormalNumber();		
 	double  generate_singleUniformNumber(double lower=0.0, double upper=1.0);		
 	void setSeed(int passedSeed=0);
    /////E added By K Fujimura /////

protected:

private:
	Vector *generatedNumbers;
	int seed;

    /////S added by K Fujimrua /////
	NormalRV *aStdNormRV;
// 	int	randomNumberBetween0And32767;
 	int	randomNumberBetween0AndRAND_MAX;
 	double randomNumberBetween0And1;
 	double randomNumber;
   /////E added by K Fujimrua /////
};

#endif

