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
// $Date: 2003-03-04 00:39:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/transformation/NatafProbabilityTransformation.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef NatafProbabilityTransformation_h
#define NatafProbabilityTransformation_h

#include <Vector.h>
#include <Matrix.h>
#include <ReliabilityDomain.h>
#include <MatrixOperations.h>

class NatafProbabilityTransformation : public ProbabilityTransformation
{

public:
	NatafProbabilityTransformation(ReliabilityDomain *passedReliabilityDomain,
						  int printFlag);
	~NatafProbabilityTransformation();

	int set_x(Vector x);
	int set_u(Vector u);

	int transform_x_to_u();
	int transform_u_to_x();
	int transform_u_to_x_andComputeJacobian();

	Vector get_x();
	Vector get_u();
	Matrix getJacobian_x_u();
	Matrix getJacobian_u_x();

	Vector meanSensitivityOf_x_to_u(Vector &x, int gradNumber);
	Vector stdvSensitivityOf_x_to_u(Vector &x, int gradNumber);

protected:

private:

	// Private data members
	Vector *x;
	Vector *u;
	Matrix *jacobian_x_u;
	Matrix *jacobian_u_x;
	ReliabilityDomain *theReliabilityDomain;
	MatrixOperations *theMatrixOperations;
	int nrv;
	Matrix *correlationMatrix;
	Matrix *lowerCholesky;
	Matrix *inverseLowerCholesky;
	int printFlag;

	// Private member functions
	void setCorrelationMatrix(int pertMeanOfThisRV, int pertStdvOfThisRV, double h);
	Matrix getJacobian_z_x(Vector x, Vector z);
	Vector z_to_x(Vector z);
	Vector x_to_z(Vector x);


};

#endif
