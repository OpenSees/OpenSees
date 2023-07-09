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
// $Date: 2008-02-29 19:47:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/meritFunction/AdkZhangMeritFunctionCheck.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <AdkZhangMeritFunctionCheck.h>
#include <MeritFunctionCheck.h>
#include <math.h>
#include <Vector.h>


AdkZhangMeritFunctionCheck::AdkZhangMeritFunctionCheck(double pmulti, double padd, double pa)
:MeritFunctionCheck()
{
	multi = pmulti;
	add = padd;
	a = pa;
}

AdkZhangMeritFunctionCheck::~AdkZhangMeritFunctionCheck()
{
}







int
AdkZhangMeritFunctionCheck::check(const Vector &u_old, 
				  double g_old, 
				  const Vector &grad_G_old, 
				  double stepSize,
				  const Vector &stepDirection,
				  double g_new,
			      int reschk  ///// added by K Fujimura /////							 				  
				  )

{

	// Update penalty parameter 'c' (should remain constant along the search direction)
		/////S added by K Fujimura /////
	//this->updateMeritParameters(u_old,g_old,grad_G_old);
	this->updateMeritParameters(u_old,g_old,grad_G_old,reschk);
		/////E added by K Fujimura /////


	// New point in standard normal space    (//different from K.F.)
	//Vector u_new = u_old + stepSize*stepDirection;
	Vector u_new(u_old);
	u_new.addVector(1.0, stepDirection, stepSize);


	// Compute value of merit functions
	static Vector dummy(1);
	double merit_old = this->getMeritFunctionValue(u_old,g_old,dummy);
	double merit_new = this->getMeritFunctionValue(u_new,g_new,dummy);

	
	// Gradient of the merit function
	double signumG;
	if (g_old != 0.0) {
		signumG = g_old/fabs(g_old);
	}
	else {
		signumG = 1.0;
	}
	//Vector gradM_old = u_old + c * signumG * grad_G_old;
	//Vector gradM_old(u_old);
	//gradM_old.addVector(1.0, grad_G_old, c*signumG);


	// Without forming temporary Vector object
	double checkValue = u_old^stepDirection;
	checkValue += c*signumG*(grad_G_old^stepDirection);
	checkValue *= a*stepSize;

	// Do the check
	//if (  (merit_new-merit_old)  <=  a*stepSize*(gradM_old^stepDirection)  ) {
	if (  (merit_new-merit_old)  <=  checkValue  ) {
		return 0;  // ok
	}
	else {
		return -1; // not ok
	}
	
}





int
AdkZhangMeritFunctionCheck::updateMeritParameters(const Vector &u, 
						  double g,
							/////S added by K Fujimura /////
						  /*const Vector &grad_G,*/ 
						  const Vector &grad_G,
						  int reschk 
						  /////E added by K Fujimura /////
						  )
{
	   /////S modificed by K Fujimura /////
	// Update penalty factor 'c'
	if(reschk==-2){
		c = (u.Norm() / grad_G.Norm());
	}else{
		c = (u.Norm() / grad_G.Norm()) * multi + add;
	}
	return 0;

//	// Update penalty factor 'c'
//	c = (u.Norm() / grad_G.Norm()) * multi + add;

	return 0;
	/////E modificed by K Fujimura /////

}


double
AdkZhangMeritFunctionCheck::getMeritFunctionValue(const Vector &u, 
						  double g,
						  const Vector &grad_G)
{
	// Note that it is correct to keep 'c' constant

	// Compute merit function
	double merit = 0.5 * (u ^ u) + c * fabs(g);


	// Return the result
	return merit;
}
