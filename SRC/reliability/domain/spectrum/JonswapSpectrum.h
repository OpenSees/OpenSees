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
// $Date: 2003-03-04 00:44:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/spectrum/JonswapSpectrum.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef JonswapSpectrum_h
#define JonswapSpectrum_h

#include <Spectrum.h>

class JonswapSpectrum : public Spectrum
{

public:
	JonswapSpectrum(int tag, double minFreq, double maxFreq,
							 double alpha, double wp, double gamma);
	~JonswapSpectrum();

	void Print(OPS_Stream &s, int flag =0);

	double getMinFrequency();
	double getMaxFrequency();
	double getAmplitude(double frequency);


protected:

private:
	double minFreq;
	double maxFreq;
	double alpha;
	double wp;
	double gamma;
};

#endif
