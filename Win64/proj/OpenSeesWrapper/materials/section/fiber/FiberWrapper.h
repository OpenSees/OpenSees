#pragma once
#include <Fiber.h>
#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>
#include <NDFiber2d.h>
#include <NDFiber3d.h>


#include "../../../taggeds/TaggedObjectWrapper.h"
#include "../../../actors/IMovableObjectWrapper.h"
#include "../../uniaxial/UniaxialMaterialWrapper.h"
#include "../../ndmaterial/NDMaterialWrapper.h"
#include "../../section/NDFiberSection2d.h"
#include "../../../matrix/VectorWrapper.h"
#include "../../../matrix/MatrixWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;
namespace OpenSees {
	namespace Materials {
		namespace Sections {
			namespace Repres {
				public ref class FiberWrapper : TaggedObjectWrapper, IMovableObjectWrapper
				{
				public:
					FiberWrapper();
					VectorWrapper^ GetFiberStressResultants() {
						Vector vec = _Fiber->getFiberStressResultants();
						return VectorWrapper::GetVectorWrapper(&vec);
					}

					MatrixWrapper^ GetFiberTangentStiffContr() {
						Matrix mat = _Fiber->getFiberTangentStiffContr();
						return MatrixWrapper::GetMatrixWrapper(&mat);
					}

					double GetArea() {
						return _Fiber->getArea();
					}

					VectorWrapper^ GetFiberLocation() {
						double y, z;
						VectorWrapper^ vec = gcnew VectorWrapper(2);
						_Fiber->getFiberLocation(y, z);
						vec[0] = y;
						vec[1] = z;
						return vec;
					}

					UniaxialMaterialWrapper^ GetMaterial() {
						UniaxialMaterialWrapper^ mat = gcnew UniaxialMaterialWrapper();
						mat->_UniaxialMaterial = _Fiber->getMaterial();
						mat->_TaggedObject = _Fiber->getMaterial();
						return mat;
					}

					NDMaterialWrapper^ GetNDMaterial() {
						NDMaterialWrapper^ mat = gcnew NDMaterialWrapper();
						mat->_NDMaterial = _Fiber->getNDMaterial();
						mat->_TaggedObject = _Fiber->getNDMaterial();
						return mat;
					}


					~FiberWrapper() {
						if (_Fiber != 0)
							delete _Fiber;
					};
				internal:
					Fiber * _Fiber;
				};

				public ref class UniaxialFiber2dWrapper : FiberWrapper
				{
				public:
					UniaxialFiber2dWrapper(int tag, UniaxialMaterialWrapper^ theMat, double area, double position);
					~UniaxialFiber2dWrapper() {
						if (_Fiber != 0)
							delete _Fiber;
					};
				private:

				};

				public ref class NDFiber2dWrapper : FiberWrapper
				{
				public:
					NDFiber2dWrapper(int tag, NDMaterialWrapper^ theMat, double Area, double position);
					~NDFiber2dWrapper() {
						if (_Fiber != 0)
							delete _Fiber;
					};
				private:

				};

				public ref class UniaxialFiber3dWrapper : FiberWrapper
				{
				public:
					UniaxialFiber3dWrapper(int tag, UniaxialMaterialWrapper^ theMat, double area, VectorWrapper^ position);
					~UniaxialFiber3dWrapper() {
						if (_Fiber != 0)
							delete _Fiber;
					};
				private:

				};

				public ref class NDFiber3dWrapper : FiberWrapper
				{
				public:
					NDFiber3dWrapper(int tag, NDMaterialWrapper^ theMat, double Area, double y, double z);
					~NDFiber3dWrapper() {
						if (_Fiber != 0)
							delete _Fiber;
					};
				private:

				};
			}
		}
	}
}
