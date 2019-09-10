#pragma once

#include <UniaxialMaterial.h>



#include "../MaterialWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
//#include "../../taggeds/TaggedObjectWrapper.h"
//#include "../../domains/domain/BaseDomainWrapper.h"
//#include "../../domains/timeSeries/TimeSeriesWrapper.h"
//#include "../../elements/BaseElementWrapper.h"
//#include "../../domains/nodes/NodeWrapper.h"
//#include "../../matrix/VectorWrapper.h"

#include "../../handlers/HandlerWrapper.h"
#include "../../recorders/ResponseWrapper.h"
#include "../../recorders/InformationWrapper.h"


using namespace System;
using namespace System::Runtime::InteropServices;

using namespace OpenSees;
using namespace OpenSees::Recorders;
using namespace OpenSees::Handlers;

//using namespace OpenSees::Components;
//using namespace OpenSees::Components::Timeseries;
//using namespace OpenSees::Elements;
//using namespace OpenSees::DamageModels;
using namespace OpenSees::Materials;


namespace OpenSees {
	namespace Materials {
		namespace Uniaxials {

			public ref class UniaxialMaterialWrapper : MaterialWrapper, IMovableObjectWrapper
			{
			public:
				UniaxialMaterialWrapper();
				UniaxialMaterialWrapper(int tag);
				~UniaxialMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

				UniaxialMaterialWrapper^ GetCopy() {
					UniaxialMaterial* mat = _UniaxialMaterial->getCopy();
					UniaxialMaterialWrapper^ theMat = gcnew UniaxialMaterialWrapper();
					theMat->_TaggedObject = mat;
					theMat->_UniaxialMaterial = mat;
					return theMat;
				}

				double GetStrain() { return _UniaxialMaterial->getStrain(); };
				double GetStrainRate() { return _UniaxialMaterial->getStrainRate(); };
				double GetStress() { return _UniaxialMaterial->getStress(); };
				double GetTangent() { return _UniaxialMaterial->getTangent(); };
				double GetInitialTangent() { return _UniaxialMaterial->getInitialTangent(); };
				double GetDampTangent() { return _UniaxialMaterial->getDampTangent(); };
				double GetRho() { return _UniaxialMaterial->getRho(); };
				int CommitState() {
					return _UniaxialMaterial->commitState();
				}
				int RevertToLastCommit() {
					return _UniaxialMaterial->revertToLastCommit();
				}
				int RevertToStart() {
					return _UniaxialMaterial->revertToStart();
				}
				String^ GetType() {
					return gcnew String(_UniaxialMaterial->getClassType());
				}
				int GetTypeTag() {
					return _UniaxialMaterial->getClassTag();
				}

				ResponseWrapper^ SetResponse(array<String^>^ argv, OpenSees::Handlers::OPS_StreamWrapper^ output) {
					const char ** _argv = new const char*[argv->Length];
					for (int i = 0; i < argv->Length; i++)
					{
						_argv[i] = OPS::StringToChar(argv[i]);
					}
					
					Response* response = this->_UniaxialMaterial->setResponse(_argv, argv->Length, *output->_OPS_StreamPtr);
					ResponseWrapper^ _response = gcnew ResponseWrapper();
					_response->_Response = response;
					return _response;
				}

				int GetResponse(int responseID, InformationWrapper^ info) {
					return _UniaxialMaterial->getResponse(responseID, *info->_Information);
				}

			internal:
				UniaxialMaterial * _UniaxialMaterial;
			};	
		}
	}
	
}
