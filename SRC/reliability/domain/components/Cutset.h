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
// $Date: 2008-05-08 15:34:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/Cutset.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef Cutset_h
#define Cutset_h

#include <ReliabilityDomainComponent.h>
#include <Vector.h>
#include <Matrix.h>

class Cutset : public ReliabilityDomainComponent
{

public:
	Cutset(int passedTag, Vector passedComponents);
	~Cutset();

	int setBetaCutset(const Vector &);
	int setRhoCutset(const Matrix &);

	int getNumberOfComponents();
	const Vector& getComponents();
	const Vector& getBetaCutset();
	const Matrix& getRhoCutset();
	
	void Print(OPS_Stream &s, int flag =0);

protected:

private:
	int numComponents;
	Vector *components;
	Vector *betas;
	Matrix *rhos;
};

#endif
