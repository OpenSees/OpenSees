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
// $Date: 2008-05-11 19:52:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/misc/CorrelatedStandardNormal.h,v $


//
// Written by Kevin Mackie (kmackie@mail.ucf.edu)
//

#ifndef CorrelatedStandardNormal_h
#define CorrelatedStandardNormal_h

class CorrelatedStandardNormal
{

public:
	CorrelatedStandardNormal(double);
	~CorrelatedStandardNormal();
	
	int setCorrelation(double);
	double bivariatePDF(double b1, double b2, double r);
	double getPDF(double b1, double b2);
	double getCDF(double b1, double b2);

protected:

private:
	double exponentialForm(double, double, double);
	double SimpsonOwen(double, double, double, double);
	double SimpsonSheppard(double, double, double);
	double getCDFowen(double b1, double b2, int popt);
	double getCDFsheppard(double b1, double b2, int popt);
	double getCDFadaptive(double b1, double b2, int popt);
	double getAdaptiveIntegralValue(double tol, double lowerBound, double upperBound, 
		double fa, double fb, double fc, double beta1, double beta2);
	void testCDF();
	
	double rho;
    

};

#endif

