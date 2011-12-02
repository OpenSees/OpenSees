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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-02-24 01:38:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/spectrum/PointsSpectrum.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef PointsSpectrum_h
#define PointsSpectrum_h

#include <Spectrum.h>
#include <Vector.h>

class PointsSpectrum : public Spectrum
{

public:
	PointsSpectrum(int tag, const Vector &frequencies, const Vector &amplitudes);
	~PointsSpectrum();

	void Print(OPS_Stream &s, int flag =0);

	double getMinFrequency();
	double getMaxFrequency();
	double getAmplitude(double frequency);


protected:

private:
	Vector frequencies;
	Vector amplitudes;
};

#endif
