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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/Cutset.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <Cutset.h>
#include <Vector.h>
#include <Matrix.h>
#include <classTags.h>
#include <OPS_Globals.h>

Cutset::Cutset(int passedTag, Vector passedComponents)
  :ReliabilityDomainComponent(passedTag, CUTSET)
{

	numComponents = passedComponents.Size();
	
	components = new Vector(numComponents);
	betas = new Vector(numComponents);
	rhos = new Matrix(numComponents,numComponents);

	for (int i = 0; i < numComponents; i++)
		(*components)(i) = passedComponents(i);

}


Cutset::~Cutset()
{
	if (components != 0)
		delete components;
	if (betas != 0)
		delete betas;
	if (rhos != 0 )
		delete rhos;
}


void
Cutset::Print(OPS_Stream &s, int flag)  
{
  s << "Cutset, tag: " << this->getTag() << endln;
  s << "\tcomponents: " << *components << endln;
}


int
Cutset::setBetaCutset(const Vector &bin)
{
	*betas = bin;
	return 0;
}

int 
Cutset::setRhoCutset(const Matrix &rin)
{
	*rhos = rin;
	return 0;
}


int
Cutset::getNumberOfComponents(void)
{
	return numComponents;
}

const Vector &
Cutset::getComponents(void)
{
	return *components;
}

const Vector &
Cutset::getBetaCutset(void)
{
	return *betas;
}

const Matrix &
Cutset::getRhoCutset(void)
{
	return *rhos;
}

