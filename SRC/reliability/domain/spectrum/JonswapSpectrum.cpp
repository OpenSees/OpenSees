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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-10-27 23:04:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/spectrum/JonswapSpectrum.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <JonswapSpectrum.h>
#include <Spectrum.h>
#include <Vector.h>
#include <classTags.h>
#include <math.h>

JonswapSpectrum::JonswapSpectrum(int tag, double min, double max,
								 double alp, double w, double gam)
:Spectrum(tag,SPECTRUM_jonswap)
{
	minFreq = min;
	maxFreq = max;
	alpha = alp;
	wp = w;
	gamma = gam;
}

JonswapSpectrum::~JonswapSpectrum()
{
}

void
JonswapSpectrum::Print(OPS_Stream &s, int flag)  
{
}


double
JonswapSpectrum::getMinFrequency()
{
	return minFreq;
}
double
JonswapSpectrum::getMaxFrequency()
{
	return maxFreq;
}
double
JonswapSpectrum::getAmplitude(double frequency)
{
	if (frequency < minFreq  ||  frequency > maxFreq) {
		return 0.0;
	}

	double sigma;
	if (frequency < wp) {
		sigma = 0.07;
	}
	else {
		sigma = 0.09;
	}

	double power = exp(-(frequency-wp)/(2.0*sigma*sigma*wp*wp));
	double GAMMA = pow(gamma,power);

	return GAMMA*alpha*pow(frequency,-5.0)*exp(-5.0/4.0*pow((frequency/wp),-4.0));
}
