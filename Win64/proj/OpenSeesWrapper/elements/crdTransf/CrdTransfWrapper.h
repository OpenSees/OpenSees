#pragma once
#include <CrdTransf.h>
#include <LinearCrdTransf3d.h>
#include <CorotCrdTransf3d.h>
#include <PDeltaCrdTransf3d.h>
#include <LinearCrdTransf2d.h>
#include <CorotCrdTransf2d.h>
#include <PDeltaCrdTransf2d.h>
#include <CorotCrdTransfWarping2d.h>
#include <dispBeamColumnInt/LinearCrdTransf2dInt.h>


#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../../matrix/VectorWrapper.h"

using namespace OpenSees;

namespace OpenSees {
	namespace Elements {
		namespace CrdTransfs {
			public ref class CrdTransfWrapper : TaggedObjectWrapper, IMovableObjectWrapper
			{
			public:
				CrdTransfWrapper();
				double GetInitialLength() {
					return _CrdTransf->getInitialLength();
				}

				double GetDeformedLength() {
					return _CrdTransf->getDeformedLength();
				}

				int GetLocalAxes(VectorWrapper^ xAxis, VectorWrapper^ yAxis, VectorWrapper^ zAxis) {
					Vector x = Vector(3);
					Vector y = Vector(3);
					Vector z = Vector(3);
					int ret = _CrdTransf->getLocalAxes(x, y, z);
					xAxis->_Vector = &x;
					yAxis->_Vector = &y;
					zAxis->_Vector = &z;
					return ret;
				}

				VectorWrapper^  getBasicTrialDisp() {
					Vector mat = _CrdTransf->getBasicTrialDisp();
					VectorWrapper^ wmat = gcnew VectorWrapper();
					wmat->_Vector = &mat;
					return wmat;
				};

				VectorWrapper^  getBasicIncrDisp() {
					Vector mat = _CrdTransf->getBasicIncrDisp();
					VectorWrapper^ wmat = gcnew VectorWrapper();
					wmat->_Vector = &mat;
					return wmat;
				};

				VectorWrapper^  getBasicIncrDeltaDisp() {
					Vector mat = _CrdTransf->getBasicIncrDeltaDisp();
					VectorWrapper^ wmat = gcnew VectorWrapper();
					wmat->_Vector = &mat;
					return wmat;
				};

				VectorWrapper^  getBasicTrialVel() {
					Vector mat = _CrdTransf->getBasicTrialVel();
					VectorWrapper^ wmat = gcnew VectorWrapper();
					wmat->_Vector = &mat;
					return wmat;
				};

				VectorWrapper^  getBasicTrialAccel() {
					Vector mat = _CrdTransf->getBasicTrialAccel();
					VectorWrapper^ wmat = gcnew VectorWrapper();
					wmat->_Vector = &mat;
					return wmat;
				};

				

				~CrdTransfWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:
				CrdTransf* _CrdTransf;
			private:

			};

			public ref class PDeltaCrdTransf2dWrapper : CrdTransfWrapper
			{
			public:
				PDeltaCrdTransf2dWrapper();
				~PDeltaCrdTransf2dWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:

			private:

			};

			public ref class LinearCrdTransf2dWrapper : CrdTransfWrapper
			{
			public:
				LinearCrdTransf2dWrapper();
				~LinearCrdTransf2dWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:

			private:

			};

			public ref class CorotCrdTransf2dWrapper : CrdTransfWrapper
			{
			public:
				CorotCrdTransf2dWrapper();
				~CorotCrdTransf2dWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:

			private:

			};

			public ref class LinearCrdTransf3dWrapper : CrdTransfWrapper
			{
			public:
				LinearCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane);
				LinearCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane, VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ);
				~LinearCrdTransf3dWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:

			private:

			};

			public ref class CorotCrdTransf3dWrapper : CrdTransfWrapper
			{
			public:
				CorotCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane, VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ);
				~CorotCrdTransf3dWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:

			private:

			};

			public ref class PDeltaCrdTransf3dWrapper : CrdTransfWrapper
			{
			public:
				PDeltaCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane);
				PDeltaCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane, VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ);
				~PDeltaCrdTransf3dWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:

			private:

			};

			public ref class CorotCrdTransfWarping2dWrapper : CrdTransfWrapper
			{
			public:
				CorotCrdTransfWarping2dWrapper(VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ);
				CorotCrdTransfWarping2dWrapper();
				~CorotCrdTransfWarping2dWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:

			private:

			};

			public ref class LinearCrdTransf2dIntWrapper : CrdTransfWrapper
			{
			public:
				LinearCrdTransf2dIntWrapper();
				LinearCrdTransf2dIntWrapper(VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ);
				~LinearCrdTransf2dIntWrapper() {
					if (_CrdTransf != 0)
						delete _CrdTransf;
				};

			internal:

			private:

			};
		}
	}
}
