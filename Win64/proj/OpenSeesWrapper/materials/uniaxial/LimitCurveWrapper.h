#pragma once



// fedeas
#include <limitstate\limitCurve\LimitCurve.h>
#include <limitstate\limitCurve\AxialCurve.h>
#include <limitstate\limitCurve\RotationShearCurve.h>
#include <limitstate\limitCurve\ShearCurve.h>
#include <limitstate\limitCurve\ThreePointCurve.h>
//#include <limitstate\limitCurve\WrapperLimitCurve.h>




#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../../matrix/IDWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../../matrix/VectorWrapper.h"
#include "../../domains/domain/BaseDomainWrapper.h"
#include "../../domains/nodes/NodeWrapper.h"
#include "../../elements/BaseElementWrapper.h"
//#include "../../materials/uniaxial/UniaxialMaterialWrapper.h"



using namespace System;
using namespace OpenSees;
using namespace OpenSees::Components;
using namespace OpenSees::Elements;
//using namespace OpenSees::Materials::Uniaxials;


namespace OpenSees {
	namespace Materials {
		namespace Uniaxials {

#pragma region LimitState
			public ref class LimitCurveWrapper abstract : TaggedObjectWrapper, IMovableObjectWrapper
			{
			public:
				LimitCurveWrapper() {};

				~LimitCurveWrapper() {
					if (_LimitCurve != 0)
						delete _LimitCurve;
				};

			internal:
				LimitCurve * _LimitCurve;
			};

			public ref class AxialCurveWrapper : LimitCurveWrapper
			{
			public:
				AxialCurveWrapper(int tag, int eleTag, BaseDomainWrapper^ theDomain,
					double Fsw, double Kdeg, double Fres,
					int defType, int forType,
					int ndI, int ndJ, int dof, int perpDirn,
					double delta, int eleRemove);

				~AxialCurveWrapper() {
					if (_LimitCurve != 0)
						delete _LimitCurve;
				};
			};
			
			public ref class RotationShearCurveWrapper : LimitCurveWrapper
			{
			public:
				RotationShearCurveWrapper(int crvTag, int eleTag,
					int ndI, int ndj, int rotAxis, double Vn, double Vr, double Kdeg, double rotLim, int defTyp,
					double b, double d, double h, double L, double st, double As, double Acc, double ld, double db, double rhot, double fc,
					double fy, double fyt, double delta, BaseDomainWrapper^ theDom, BaseElementWrapper^ theEle, NodeWrapper^ theNdI, NodeWrapper^ theNdJ);

				~RotationShearCurveWrapper() {
					if (_LimitCurve != 0)
						delete _LimitCurve;
				};
			};
			
			public ref class ThreePointCurveWrapper : LimitCurveWrapper
			{
			public:
				ThreePointCurveWrapper(int tag, int eleTag, BaseDomainWrapper^ theDomain,
					double x1, double y1,
					double x2, double y2,
					double x3, double y3,
					double Kdeg, double Fres,
					int defType, int forType,
					int ndI, int ndJ, int dof, int perpDirn);

				~ThreePointCurveWrapper() {
					if (_LimitCurve != 0)
						delete _LimitCurve;
				};
			};

			


#pragma endregion

			

			

			

			

			
		}
	}
	
}
