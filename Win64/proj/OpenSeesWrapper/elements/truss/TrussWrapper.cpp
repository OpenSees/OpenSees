#include "stdafx.h"

#include "TrussWrapper.h"
using namespace OpenSees::Elements;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::Sections;

TrussWrapper::TrussWrapper(int tag, int dimension, int Nd1, int Nd2, UniaxialMaterialWrapper^ theMaterial, double A, double rho, int doRayleighDamping,
	int cMass)
{
	_Element = new Truss(tag, dimension, Nd1, Nd2, *theMaterial->_UniaxialMaterial, A, rho, doRayleighDamping, cMass);
}

CorotTrussWrapper::CorotTrussWrapper(int tag, int dim,
	int Nd1, int Nd2,
	UniaxialMaterialWrapper^ theMaterial,
	double A, double rho,
	int doRayleighDamping,
	int cMass)
{
	_Element = new CorotTruss(tag, dim,
		Nd1, Nd2,
		*theMaterial->_UniaxialMaterial,
		A, rho,
		doRayleighDamping,
		cMass);
}

CorotTruss2Wrapper::CorotTruss2Wrapper(int tag, int dim,
	int Nd1, int Nd2, int oNd1, int oNd2,
	UniaxialMaterialWrapper^ theMaterial,
	double A, double rho)
{
	_Element = new CorotTruss2(tag, dim,
		Nd1, Nd2, oNd1, oNd2,
		*theMaterial->_UniaxialMaterial,
		A, rho);
}

CorotTrussSectionWrapper::CorotTrussSectionWrapper(int tag, int dim,
	int Nd1, int Nd2,
	SectionForceDeformationWrapper^ theMaterial,
	double rho,
	int doRayleighDamping,
	int cMass)
{
	_Element = new CorotTrussSection(tag, dim,
		Nd1, Nd2,
		*theMaterial->_SectionForceDeformation,
		rho,
		doRayleighDamping,
		cMass);
}

N4BiaxialTrussWrapper::N4BiaxialTrussWrapper(int tag,
	int dimension,
	int Nd1, int Nd2,
	int GNd1, int GNd2,
	UniaxialMaterialWrapper^ theMaterial,
	double A,
	double rho,
	int doRayleighDamping)
{
	_Element = new N4BiaxialTruss(tag,
		dimension,
		Nd1, Nd2,
		GNd1, GNd2,
		*theMaterial->_UniaxialMaterial,
		A,
		rho,
		doRayleighDamping);
}


Truss2Wrapper::Truss2Wrapper(int tag,
	int dimension,
	int Nd1, int Nd2, int oNd1, int oNd2,
	UniaxialMaterialWrapper^ theMaterial,
	double A,
	double rho,
	int doRayleighDamping)
{
	_Element = new Truss2(tag,
		dimension,
		Nd1, Nd2, oNd1, oNd2,
		*theMaterial->_UniaxialMaterial,
		A,
		rho,
		doRayleighDamping);
}

TrussSectionWrapper::TrussSectionWrapper(int tag, int dimension,
	int Nd1, int Nd2,
	SectionForceDeformationWrapper^ theSection,
	double rho,
	int doRayleighDamping,
	int cMass)
{
	_Element = new TrussSection(tag, dimension,
		Nd1, Nd2,
		*theSection->_SectionForceDeformation,
		rho,
		doRayleighDamping,
		cMass);
}


