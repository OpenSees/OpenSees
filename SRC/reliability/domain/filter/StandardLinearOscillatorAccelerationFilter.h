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
// $Date: 2008-03-13 22:36:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/filter/StandardLinearOscillatorAccelerationFilter.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef StandardLinearOscillatorAccelerationFilter_h
#define StandardLinearOscillatorAccelerationFilter_h

#include <Filter.h>

class StandardLinearOscillatorAccelerationFilter : public Filter
{

public:
	StandardLinearOscillatorAccelerationFilter(int tag, double period, double dampingRatio);
	~StandardLinearOscillatorAccelerationFilter();
	double getAmplitude(double time, double dT = 0.0);
	double getMaxAmplitude();
	double getTimeOfMaxAmplitude();

	void Print(OPS_Stream &s, int flag =0);

protected:

private:
	double wn;
	double xi;
};

#endif
