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
// $Date: 2008-03-13 22:36:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/filter/DeltaFilter.h,v $

//  Written by Quan and Michele  Dec. 2005 


#ifndef DeltaFilter_h
#define DeltaFilter_h

#include <Filter.h>

class DeltaFilter : public Filter
{

public:
	DeltaFilter(int tag);
	~DeltaFilter();
	double getAmplitude(double time,double dT =0.0);
	double getMaxAmplitude();
	double getTimeOfMaxAmplitude();

	void Print(OPS_Stream &s, int flag =0);

protected:

private:

};

#endif



