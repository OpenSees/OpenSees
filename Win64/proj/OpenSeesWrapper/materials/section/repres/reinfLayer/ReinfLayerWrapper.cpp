#include "stdafx.h"
#include "ReinfLayerWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials::Sections::Repres;
ReinfLayerWrapper::ReinfLayerWrapper()
{

}

StraightReinfLayerWrapper::StraightReinfLayerWrapper(int materialID, int numReinfBars, double reinfBarArea, VectorWrapper ^ InitialPosition, VectorWrapper ^ FinalPosition)
{
	_ReinfLayer = new StraightReinfLayer(materialID, numReinfBars, reinfBarArea, *InitialPosition->_Vector, *FinalPosition->_Vector);
}

CircReinfLayerWrapper::CircReinfLayerWrapper(int materialID, int numReinfBars, double  reinfBarArea,
	VectorWrapper^ centerPosition, double arcRadius, double
	initialAngle, double finalAngle)
{
	_ReinfLayer = new CircReinfLayer(materialID, numReinfBars, reinfBarArea, *centerPosition->_Vector,
		arcRadius, initialAngle, finalAngle);
}
