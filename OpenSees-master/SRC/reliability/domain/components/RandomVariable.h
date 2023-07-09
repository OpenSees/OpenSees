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
                                                                        
// $Revision: 1.17 $
// $Date: 2008-04-10 16:22:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariable.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef RandomVariable_h
#define RandomVariable_h

#include <ReliabilityDomainComponent.h>
#include <Vector.h>

class RandomVariable : public ReliabilityDomainComponent
{

public:
	RandomVariable(int tag, int classTag);
	virtual ~RandomVariable();
	
	// pure virtual defining variable type and properties
	virtual const char* getType() = 0;
	virtual double getMean() = 0;
	virtual double getStdv() = 0;
	virtual const Vector &getParameters() = 0;
	virtual int setParameters(double mean, double stdv) {/*MHS 9/28/2011*/return 0;}
	
	// RV functionality
	virtual double getPDFvalue(double rvValue) = 0;
	virtual double getCDFvalue(double rvValue) = 0;
	virtual double getInverseCDFvalue(double rvValue) = 0; 
	
	// starting realization
	virtual int setStartValue(double newVal) {startValue = newVal; return 0;}
	virtual double getStartValue() {return startValue;}

	// current realization
	virtual int setCurrentValue(double newVal) {rvValue = newVal; return 0;}
	virtual double getCurrentValue() {return rvValue;}
    
    // standardization of random variables
    virtual double transform_x_to_u(void);
    virtual double transform_u_to_x(double uVal);
    virtual double gradient_x_to_u(double uVal);
	
	// sensitivity of CDF with respect to distribution parameters
	double getCDFMeanSensitivity(void);
	double getCDFStdvSensitivity(void);
    virtual int getCDFparameterSensitivity(Vector &dFdP) {return 0;}
    virtual int getParameterMeanSensitivity(Vector &dPdmu) {return 0;}
    virtual int getParameterStdvSensitivity(Vector &dPdstdv) {return 0;}
	
	// other public functions
	virtual void Print(OPS_Stream &s, int flag = 0);
	int setNewTag(int tag);
	
    // utility functions for gamma and beta
	double gammaFunction(double x);
	double incompleteGammaFunction(double a, double x);
	double betaFunction(double passed_q, double passed_r);
	double incompleteBetaFunction(double passed_q, double passed_r, double passed_x);
    
    // standard normal utility functions
	double errorFunction(double x);
	double inverseErrorFunction(double y);
    double standardNormalPhi(double u);
    double standardNormalInversePhi(double p);
    double harmonicNumber(double x);
    

protected:
	static const double pi;
	static const double euler;
    static const double zeta3;
    static const double zeta5;
	
private:
	double startValue;
	double rvValue;
};

#endif

