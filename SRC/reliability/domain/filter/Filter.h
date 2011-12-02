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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-03-13 22:36:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/filter/Filter.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef Filter_h
#define Filter_h

#include <Vector.h>
#include <ReliabilityDomainComponent.h>

class Filter : public ReliabilityDomainComponent
{

public:
	Filter(int tag, int classtag);
	virtual ~Filter();
//	virtual double getAmplitude(double time) = 0;
	virtual double getAmplitude(double time,double dT =0.0) = 0;
	virtual double getMaxAmplitude() = 0;
	virtual double getTimeOfMaxAmplitude() = 0;
	virtual void setKickTime(double time);	///// added by K Fujimura /////

protected:

private:

};

#endif
