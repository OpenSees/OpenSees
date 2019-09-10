#pragma once
#include <ReinfLayer.h>
#include <StraightReinfLayer.h>
#include <CircReinfLayer.h>

#include "../../../../matrix/VectorWrapper.h"

using namespace OpenSees;
namespace OpenSees {
	namespace Materials {
		namespace Sections {
			namespace Repres {
				public ref class ReinfLayerWrapper
				{
				public:
					ReinfLayerWrapper();
					~ReinfLayerWrapper() {
						if (_ReinfLayer != 0)
							delete _ReinfLayer;
					};
				internal:
					ReinfLayer * _ReinfLayer;
				};

				public ref class StraightReinfLayerWrapper : ReinfLayerWrapper
				{
				public:
					StraightReinfLayerWrapper(int materialID, int numReinfBars,
						double reinfBarArea,
						VectorWrapper^ InitialPosition,
						VectorWrapper^ FinalPosition);
					~StraightReinfLayerWrapper() {
						if (_ReinfLayer != 0)
							delete _ReinfLayer;
					};
				internal:

				};

				public ref class CircReinfLayerWrapper : ReinfLayerWrapper
				{
				public:
					CircReinfLayerWrapper(int materialID, int numReinfBars, double  reinfBarArea,
						VectorWrapper^ centerPosition, double arcRadius, double
						initialAngle, double finalAngle);
					~CircReinfLayerWrapper() {
						if (_ReinfLayer != 0)
							delete _ReinfLayer;
					};
				internal:

				};
			}
		}
	}
}
