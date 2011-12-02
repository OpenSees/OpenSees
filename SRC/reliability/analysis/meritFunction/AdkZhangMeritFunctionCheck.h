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
// $Date: 2008-02-29 19:47:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/meritFunction/AdkZhangMeritFunctionCheck.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef AdkZhangMeritFunctionCheck_h
#define AdkZhangMeritFunctionCheck_h

#include <MeritFunctionCheck.h>

class AdkZhangMeritFunctionCheck : public MeritFunctionCheck
{

public:
	AdkZhangMeritFunctionCheck(double multi, double add, double a=0.5);
	~AdkZhangMeritFunctionCheck();

	int	check(const Vector &u_old, 
		      double g_old, 
		      const Vector &grad_G_old, 
		      double stepSize,
		      const Vector &stepDirection,
		      double g_new,
/////S added byt K Fujimura
			  int reschk=0);
	double getMeritFunctionValue(const Vector &u, double g, const Vector &grad_G);
	int updateMeritParameters(const Vector &u, double g, const Vector &grad_G,
		int reschk);
 
	/*double getMeritFunctionValue(const Vector &u, double g,
				     const Vector &grad_G);
	int updateMeritParameters(const Vector &u, double g,
				  const Vector &grad_G);*/
 /////E added byt K Fujimura

protected:

private:
	double multi, add;
	double a;
	double c;

};

#endif
