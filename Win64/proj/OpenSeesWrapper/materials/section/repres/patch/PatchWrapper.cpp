#include "stdafx.h"
#include "PatchWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials::Sections::Repres;
PatchWrapper::PatchWrapper()
{

}

QuadPatchWrapper::QuadPatchWrapper(int materialID, int numSubdivIJ,
	int numSubdivJK, MatrixWrapper^ vertexCoords)
{
	_Patch = new QuadPatch(materialID, numSubdivIJ, numSubdivJK, *vertexCoords->_Matrix);
}

CircPatchWrapper::CircPatchWrapper(int materialID, int numSubdivCircunf, int numSubdivRadial,
	VectorWrapper^ centerPosition, double internRadius,
	double externRadius, double initialAngle, double finalAngle)
{
	_Patch = new CircPatch(materialID, numSubdivCircunf, numSubdivRadial, *centerPosition->_Vector,
		internRadius, externRadius, initialAngle, finalAngle);
}

