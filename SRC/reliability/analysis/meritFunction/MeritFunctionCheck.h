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
// $Date: 2003-03-04 00:39:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/meritFunction/MeritFunctionCheck.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef MeritFunctionCheck_h
#define MeritFunctionCheck_h

#include <Vector.h>

class MeritFunctionCheck
{

public:
	MeritFunctionCheck();
	virtual ~MeritFunctionCheck();

	virtual int	check(Vector u_old, 
					  double g_old, 
					  Vector grad_G_old, 
					  double stepSize,
					  Vector stepDirection,
					  double g_new, 
					  Vector grad_G_new) = 0;
	virtual double getMeritFunctionValue(Vector u, double g, Vector grad_G) = 0;
	virtual int updateMeritParameters(Vector u, double g, Vector grad_G) = 0;

protected:

private:

};

#endif
