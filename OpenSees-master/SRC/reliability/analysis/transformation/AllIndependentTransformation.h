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
// $Date: 2008-02-29 19:43:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/AllIndependentTransformation.h,v $

#ifndef AllIndependentTransformation_h
#define AllIndependentTransformation_h

#include <Vector.h>
#include <Matrix.h>
#include <ReliabilityDomain.h>
#include <MatrixOperations.h>
#include <ProbabilityTransformation.h>

class AllIndependentTransformation : public ProbabilityTransformation
{
protected:

public:
	AllIndependentTransformation(ReliabilityDomain *passedReliabilityDomain,
						  int printFlag);
	~AllIndependentTransformation();

	int transform_x_to_u(Vector &u);
	int transform_u_to_x(const Vector &u, Vector &x);
	int getJacobian_x_to_u(Matrix &Jxu);
	int getJacobian_u_to_x(const Vector &u, Matrix &Jux);

	Vector meanSensitivityOf_x_to_u(const Vector &x, int rvTag);
	Vector stdvSensitivityOf_x_to_u(const Vector &x, int rvTag);

    void setReliabilityDomain(ReliabilityDomain* theRelDom);

protected:

private:

	// Private data members
	Vector *x;
	Vector *xtemp;
	Vector *u;
	Vector *z;
	Vector *DzDmean;
	Vector *DzDstdv;
	Matrix *jacobian_x_u;
	Matrix *jacobian_u_x;
	Matrix *jacobian_z_x;
	ReliabilityDomain *theReliabilityDomain;
	int nrv;
	Matrix *correlationMatrix;
	Matrix *lowerCholesky;
	Matrix *inverseLowerCholesky;
	int printFlag;

	// Private member functions
	// Representing diag matrix as vector
	int getJacobian_z_x(const Vector &z, Vector &Jzx);
	int z_to_x(const Vector &z, Vector &x);
	int x_to_z(Vector &z);

};

#endif
