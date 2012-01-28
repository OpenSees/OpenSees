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
// $Date: 2007-11-06 19:32:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/transformation/ProbabilityTransformation.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef ProbabilityTransformation_h
#define ProbabilityTransformation_h

#include <Vector.h>
#include <Matrix.h>

class ProbabilityTransformation
{

public:
	ProbabilityTransformation();
	virtual ~ProbabilityTransformation();

	virtual int transform_x_to_u(Vector &u) = 0;
	virtual int transform_u_to_x(const Vector &u, Vector &x) = 0;
	virtual int getJacobian_x_to_u(Matrix &Jxu) = 0;
	virtual int getJacobian_u_to_x(const Vector &u, Matrix &Jux) = 0;

	virtual Vector meanSensitivityOf_x_to_u(const Vector &x, int rvTag) = 0;
	virtual Vector stdvSensitivityOf_x_to_u(const Vector &x, int rvTag) = 0;

protected:

private:

};

#endif
