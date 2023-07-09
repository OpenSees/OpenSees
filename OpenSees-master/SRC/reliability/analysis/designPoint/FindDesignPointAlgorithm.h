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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-05-13 16:30:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/designPoint/FindDesignPointAlgorithm.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef FindDesignPointAlgorithm_h
#define FindDesignPointAlgorithm_h

#include <Vector.h>

class FindDesignPointAlgorithm
{

public:
	FindDesignPointAlgorithm();
	virtual ~FindDesignPointAlgorithm();

    virtual int gradientStandardNormal(double) =0;
	virtual int findDesignPoint() =0;

	virtual const Vector &get_x() =0;
	virtual const Vector &get_u() =0;
	virtual const Vector &get_alpha() =0;
	virtual const Vector &get_gamma() =0;
	virtual int getNumberOfSteps() =0;
    virtual double getFirstCurvature() =0;
	virtual const Vector &getLastSearchDirection() =0;
	virtual double getFirstGFunValue() =0;
	virtual double getLastGFunValue() =0;
	virtual const Vector &getGradientInStandardNormalSpace() =0;
    virtual const Vector &getGradientInOriginalSpace() =0;
	virtual int getNumberOfEvaluations() = 0;
	
	/////S added by K Fujimura /////
	virtual int    getNumberOfSensAna();
	virtual double get_check1_init();
	virtual double get_check2_init();
	virtual double get_check1_conv();
	virtual double get_check2_conv();
	
	virtual void set_u(Vector&);
	virtual void set_x(Vector&);
	virtual double get_beta();
	virtual Matrix getJacobian_x_u();
	/////E added by K Fujimura /////

protected:

private:

};

#endif
