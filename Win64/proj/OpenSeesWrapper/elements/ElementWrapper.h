#pragma once
#include <Element.h>
#include <InformationProxy.h>

#include "BaseElementWrapper.h"

#include "../matrix/IDWrapper.h"
#include "../matrix/VectorWrapper.h"
#include "../matrix/MatrixWrapper.h"
#include "../materials/uniaxial/UniaxialMaterialWrapper.h"
#include "../domains/components/DomainComponentWrapper.h"
#include "../domains/nodes/NodeWrapper.h"
//#include "../handlers/HandlerWrapper.h"
//#include "../recorders/ResponseWrapper.h"
#include "../recorders/InformationWrapper.h"

//#include "InformationWrapper.h"
using namespace OpenSees;
using namespace OpenSees::Recorders;
//using namespace OpenSees::Handlers;
using namespace OpenSees::Components;
using namespace OpenSees::Materials::Uniaxials;
using namespace System;
using namespace System::Runtime::InteropServices;


namespace OpenSees {
	namespace Elements {

		public ref class ElementWrapper : BaseElementWrapper
		{
		public:
			ElementWrapper();

			IDWrapper^ GetNodeTagsIDWrapper() {
				return IDWrapper::GetIDWrapper(_Element->getExternalNodes());
			}

			array<int>^ GetNodeTags() {
				return IDWrapper::GetArray(_Element->getExternalNodes());
			}

			MatrixWrapper^ GetTangentStiffMatrixWrapper() {
				return MatrixWrapper::GetMatrixWrapper(_Element->getTangentStiff());
			}

			array<double, 2>^ GetTangentStiff() {
				return MatrixWrapper::GetArray(_Element->getTangentStiff());
			}

			MatrixWrapper^ GetMassMatrixWrapper() {
				return MatrixWrapper::GetMatrixWrapper(_Element->getMass());
			}

			array<double, 2>^ GetMass() {
				return MatrixWrapper::GetArray(_Element->getMass());
			}

			MatrixWrapper^ GetDampMatrixWrapper() {
				return MatrixWrapper::GetMatrixWrapper(_Element->getDamp());
			}

			array<double, 2>^ GetDamp() {
				return MatrixWrapper::GetArray(_Element->getDamp());;
			}

			VectorWrapper^ GetResistingForceVectorWrapper() {
				return VectorWrapper::GetVectorWrapper(_Element->getResistingForce());
			}

			array<double>^ GetResistingForce() {
				return VectorWrapper::GetArray(_Element->getResistingForce());
			}

			VectorWrapper^ GetResistingForceIncInertiaVectorWrapper() {
				return VectorWrapper::GetVectorWrapper(_Element->getResistingForceIncInertia());
			}

			array<double>^ GetResistingForceIncInertia() {
				return VectorWrapper::GetArray(_Element->getResistingForceIncInertia());
			}

			VectorWrapper^ GetRayleighDampingForcesVectorWrapper() {
				return VectorWrapper::GetVectorWrapper(_Element->getRayleighDampingForces());
			}

			array<double>^ GetRayleighDampingForces() {
				return VectorWrapper::GetArray(_Element->getRayleighDampingForces());
			}

			virtual int GetTag() override {
				return _Element->getTag();
			}

			int GetNumExternalNodes() {
				return _Element->getNumExternalNodes();
			}

			int GetNumDOF() {
				return _Element->getNumDOF();
			}

			array<NodeWrapper^>^ GetNodes() {
				int num = _Element->getNumExternalNodes();
				Node** nodes = _Element->getNodePtrs();
				array<NodeWrapper^>^ _nodes = gcnew array<NodeWrapper^>(num);
				for (int i = 0; i < num; i++)
				{
					_nodes[i] = gcnew NodeWrapper(nodes[i]);
				}
				return _nodes;
			}

			String^ GetClassType() {
				const char* type = _Element->getClassType();
				return gcnew String(type);
			}

			/*ResponseWrapper^ SetResponse(array<String^>^ argv, OpenSees::Handlers::OPS_StreamWrapper^ output) {
				char ** _argv = new char*[argv->Length];
				for (int i = 0; i < argv->Length; i++)
				{
					_argv[i] = OPS::StringToChar(argv[i]);
				}

				Response* response = this->_Element->setResponse(_argv, argv->Length, *output->_OPS_StreamPtr);
				ResponseWrapper^ _response = gcnew ResponseWrapper();
				_response->_Response = response;
				return _response;
			}*/

			int GetResponse(int responseID, InformationWrapper^ info) {
				return _Element->getResponse(responseID, *info->_Information);
			}

			~ElementWrapper();


			String^ GetStringResponse(int id) {
				if (_InformationProxy == 0)
					_InformationProxy = new InformationProxy(this->_Element);
				char* vec = _InformationProxy->getString(id);
				String^ vecw = gcnew String(vec);
				return vecw;
			};

			double GetDoubleResponse(int id) {
				if (_InformationProxy == 0)
					_InformationProxy = new InformationProxy(this->_Element);

				return _InformationProxy->getDouble(id);
			};

			void GetVectorResponse(int id, VectorWrapper^% rvec) {
				if (_InformationProxy == 0)
					_InformationProxy = new InformationProxy(this->_Element);
				if (rvec->_Vector != 0)
					rvec->~VectorWrapper();

				Vector* vec = _InformationProxy->getVector(id);
				rvec->_Vector = vec;
				return;
			};

			void GetMatrixResponse(int id, MatrixWrapper^% rmat) {
				if (_InformationProxy == 0)
					_InformationProxy = new InformationProxy(this->_Element);

				if (rmat->_Matrix != 0)
					rmat->~MatrixWrapper();

				Matrix* mat = _InformationProxy->getMatrix(id);
				rmat->_Matrix = mat;
				return;
			};

			void GetIDResponse(int id, IDWrapper^% rid) {
				if (_InformationProxy == 0)
					_InformationProxy = new InformationProxy(this->_Element);

				if (rid->_ID != 0)
					rid->~IDWrapper();

				ID* vec = _InformationProxy->getID(id);
				rid->_ID = vec;
				return;
			};

			int GetIntResponse(int id) {
				if (_InformationProxy == 0)
					_InformationProxy = new InformationProxy(this->_Element);
				return _InformationProxy->getInt(id);
			};


		internal:
			//Element* _Element;
			InformationProxy* _InformationProxy;
		private:

		};

	}
}
