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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-03-13 22:31:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/hessianApproximation/HessianApproximation.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef HessianApproximation_h
#define HessianApproximation_h

#include <Matrix.h>
#include <Vector.h>

class HessianApproximation
{

public:
	HessianApproximation();
	virtual ~HessianApproximation();

	virtual const Matrix &getHessianApproximation() = 0;
	virtual int setHessianToIdentity(int size) = 0;
	virtual int updateHessianApproximation(const Vector &u_old,
					       double g_old,
					       const Vector &gradG_old,
					       double stepSize,
					       const Vector &searchDirection,
					       double g_new,
					       const Vector &gradG_new) = 0;
protected:

private:

};

#endif
