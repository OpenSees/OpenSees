#pragma once
#include <ExternalElement.h>
#include "../ElementWrapper.h"
#include "../../domains/domain/DomainWrapper.h"
#include "../../domains/loads/LoadWrapper.h"
#include "../../recorders/InformationWrapper.h"
#include "../../recorders/ResponseWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Components;
using namespace OpenSees::Components::Loads;
using namespace OpenSees::Recorders;

namespace OpenSees {
	namespace Elements {
		public delegate String^ DEle_GetClassType();
		delegate const char* UDEle_GetClassType();

		public delegate int DEle_GetNumExternalNodes(void);
		delegate int UDEle_GetNumExternalNodes(void);

		public delegate IDWrapper^ DEle_GetExternalNodes(void);
		delegate ID* UDEle_GetExternalNodes(void);

		public delegate array<NodeWrapper^>^ DEle_GetNodePtrs(void);
		delegate Node** UDEle_GetNodePtrs(void);

		public delegate int DEle_GetNumDOF(void);
		delegate int UDEle_GetNumDOF(void);

		public delegate void DEle_SetDomain(DomainWrapper^ theDomain);
		delegate void UDEle_SetDomain(Domain* theDomain);

		public delegate int DEle_CommitState(void);
		delegate int UDEle_CommitState(void);

		public delegate int DEle_RevertToLastCommit(void);
		delegate int UDEle_RevertToLastCommit(void);

		public delegate int DEle_RevertToStart(void);
		delegate int UDEle_RevertToStart(void);

		public delegate int DEle_Update(void);
		delegate int UDEle_Update(void);

		public delegate MatrixWrapper^ DEle_GetTangentStiff(void);
		delegate Matrix * UDEle_GetTangentStiff(void);

		public delegate MatrixWrapper^ DEle_GetInitialStiff(void);
		delegate Matrix * UDEle_GetInitialStiff(void);

		public delegate MatrixWrapper^ DEle_GetDamp(void);
		delegate Matrix * UDEle_GetDamp(void);

		public delegate MatrixWrapper^ DEle_GetMass(void);
		delegate Matrix * UDEle_GetMass(void);

		public delegate void DEle_ZeroLoad(void);
		delegate void UDEle_ZeroLoad(void);

		public delegate int DEle_AddLoad(ElementalLoadWrapper^ theLoad, double loadFactor);
		delegate int UDEle_AddLoad(ElementalLoad *theLoad, double loadFactor);

		public delegate int DEle_AddInertiaLoadToUnbalance(VectorWrapper^ accel);
		delegate int UDEle_AddInertiaLoadToUnbalance(const Vector &accel);

		public delegate VectorWrapper^ DEle_GetResistingForce(void);
		delegate Vector * UDEle_GetResistingForce(void);

		public delegate VectorWrapper^ DEle_GetResistingForceIncInertia(void);
		delegate Vector * UDEle_GetResistingForceIncInertia(void);

		public delegate void DEle_Print(int);
		delegate void UDEle_Print(int);

		public delegate ResponseWrapper^ DEle_SetResponse(array<String^>^ argv, OPS_StreamWrapper^ output);
		delegate Response* UDEle_SetResponse(const char **argv, int argc, OPS_Stream &output);

		public delegate int DEle_GetResponse(int responseID, InformationWrapper^ eleInfo);
		delegate int UDEle_GetResponse(int responseID, Information &eleInfo);

		public ref class ExternalElementWrapper : ElementWrapper
		{
		public:
			ExternalElementWrapper(int tag);
			~ExternalElementWrapper();
			void SetLinks(
				DEle_GetClassType^ _Ele_GetClassType,
				DEle_GetNumExternalNodes^ _Ele_GetNumExternalNodes,
				DEle_GetNodePtrs^ _Ele_GetNodePtrs,
				DEle_GetExternalNodes^ _Ele_GetExternalNodes,
				DEle_GetNumDOF^ _Ele_GetNumDOF,
				DEle_SetDomain^ _Ele_SetDomain,
				DEle_CommitState^ _Ele_CommitState,
				DEle_RevertToLastCommit^ _Ele_RevertToLastCommit,
				DEle_RevertToStart^ _Ele_RevertToStart,
				DEle_Update^ _Ele_Update,
				DEle_GetTangentStiff^ _Ele_GetTangentStiff,
				DEle_GetInitialStiff^ _Ele_GetInitialStiff,
				DEle_GetDamp^ _Ele_GetDamp,
				DEle_GetMass^ _Ele_GetMass,
				DEle_ZeroLoad^ _Ele_ZeroLoad,
				DEle_AddLoad^ _Ele_AddLoad,
				DEle_AddInertiaLoadToUnbalance^ _Ele_AddInertiaLoadToUnbalance,
				DEle_GetResistingForce^ _Ele_GetResistingForce,
				DEle_GetResistingForceIncInertia^ _Ele_GetResistingForceIncInertia,
				DEle_Print^ _Ele_Print,
				DEle_SetResponse^ _Ele_SetResponse,
				DEle_GetResponse^ _Ele_GetResponse
			)
			{
				this->_Ele_GetClassType = _Ele_GetClassType;
				UDEle_GetClassType^ udEle_GetClassType = gcnew UDEle_GetClassType(this, &ExternalElementWrapper::GetClassType_Method);
				udEle_GetClassType_GC = GCHandle::Alloc(udEle_GetClassType);
				IntPtr ip = Marshal::GetFunctionPointerForDelegate(udEle_GetClassType);
				Ele_GetClassType GetClassType = static_cast<Ele_GetClassType>(ip.ToPointer());

				this->_Ele_GetNumExternalNodes = _Ele_GetNumExternalNodes;
				UDEle_GetNumExternalNodes^ udEle_GetNumExternalNodes = gcnew UDEle_GetNumExternalNodes(this, &ExternalElementWrapper::GetNumExternalNodes_Method);
				udEle_GetNumExternalNodes_GC = GCHandle::Alloc(udEle_GetNumExternalNodes);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetNumExternalNodes);
				Ele_GetNumExternalNodes GetNumExternalNodes = static_cast<Ele_GetNumExternalNodes>(ip.ToPointer());

				this->_Ele_GetExternalNodes = _Ele_GetExternalNodes;
				UDEle_GetExternalNodes^ udEle_GetExternalNodes = gcnew UDEle_GetExternalNodes(this, &ExternalElementWrapper::GetExternalNodes_Method);
				udEle_GetExternalNodes_GC = GCHandle::Alloc(udEle_GetExternalNodes);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetExternalNodes);
				Ele_GetExternalNodes GetExternalNodes = static_cast<Ele_GetExternalNodes>(ip.ToPointer());

				this->_Ele_GetNodePtrs = _Ele_GetNodePtrs;
				UDEle_GetNodePtrs^ udEle_GetNodePtrs = gcnew UDEle_GetNodePtrs(this, &ExternalElementWrapper::GetNodePtrs_Method);
				udEle_GetNodePtrs_GC = GCHandle::Alloc(udEle_GetNodePtrs);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetNodePtrs);
				Ele_GetNodePtrs GetNodePtrs = static_cast<Ele_GetNodePtrs>(ip.ToPointer());

				this->_Ele_GetNumDOF = _Ele_GetNumDOF;
				UDEle_GetNumDOF^ udEle_GetNumDOF = gcnew UDEle_GetNumDOF(this, &ExternalElementWrapper::GetNumDOF_Method);
				udEle_GetNumDOF_GC = GCHandle::Alloc(udEle_GetNumDOF);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetNumDOF);
				Ele_GetNumDOF GetNumDOF = static_cast<Ele_GetNumDOF>(ip.ToPointer());
				
				this->_Ele_SetDomain = _Ele_SetDomain;
				UDEle_SetDomain^ udEle_SetDomain = gcnew UDEle_SetDomain(this, &ExternalElementWrapper::SetDomain_Method);
				udEle_SetDomain_GC = GCHandle::Alloc(udEle_SetDomain);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_SetDomain);
				Ele_SetDomain SetDomain = static_cast<Ele_SetDomain>(ip.ToPointer());

				this->_Ele_CommitState = _Ele_CommitState;
				UDEle_CommitState^ udEle_CommitState = gcnew UDEle_CommitState(this, &ExternalElementWrapper::CommitState_Method);
				udEle_CommitState_GC = GCHandle::Alloc(udEle_CommitState);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_CommitState);
				Ele_CommitState CommitState = static_cast<Ele_CommitState>(ip.ToPointer());

				this->_Ele_RevertToLastCommit = _Ele_RevertToLastCommit;
				UDEle_RevertToLastCommit^ udEle_RevertToLastCommit = gcnew UDEle_RevertToLastCommit(this, &ExternalElementWrapper::RevertToLastCommit_Method);
				udEle_RevertToLastCommit_GC = GCHandle::Alloc(udEle_RevertToLastCommit);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_RevertToLastCommit);
				Ele_RevertToLastCommit RevertToLastCommit = static_cast<Ele_RevertToLastCommit>(ip.ToPointer());

				this->_Ele_RevertToStart = _Ele_RevertToStart;
				UDEle_RevertToStart^ udEle_RevertToStart = gcnew UDEle_RevertToStart(this, &ExternalElementWrapper::RevertToStart_Method);
				udEle_RevertToStart_GC = GCHandle::Alloc(udEle_RevertToStart);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_RevertToStart);
				Ele_RevertToStart RevertToStart = static_cast<Ele_RevertToStart>(ip.ToPointer());

				this->_Ele_Update = _Ele_Update;
				UDEle_Update^ udEle_Update = gcnew UDEle_Update(this, &ExternalElementWrapper::Update_Method);
				udEle_Update_GC = GCHandle::Alloc(udEle_Update);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_Update);
				Ele_Update Update = static_cast<Ele_Update>(ip.ToPointer());

				this->_Ele_GetTangentStiff = _Ele_GetTangentStiff;
				UDEle_GetTangentStiff^ udEle_GetTangentStiff = gcnew UDEle_GetTangentStiff(this, &ExternalElementWrapper::GetTangentStiff_Method);
				udEle_GetTangentStiff_GC = GCHandle::Alloc(udEle_GetTangentStiff);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetTangentStiff);
				Ele_GetTangentStiff GetTangentStiff = static_cast<Ele_GetTangentStiff>(ip.ToPointer());

				this->_Ele_GetInitialStiff = _Ele_GetInitialStiff;
				UDEle_GetInitialStiff^ udEle_GetInitialStiff = gcnew UDEle_GetInitialStiff(this, &ExternalElementWrapper::GetInitialStiff_Method);
				udEle_GetInitialStiff_GC = GCHandle::Alloc(udEle_GetInitialStiff);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetInitialStiff);
				Ele_GetInitialStiff GetInitialStiff = static_cast<Ele_GetInitialStiff>(ip.ToPointer());

				this->_Ele_GetDamp = _Ele_GetDamp;
				UDEle_GetDamp^ udEle_GetDamp = gcnew UDEle_GetDamp(this, &ExternalElementWrapper::GetDamp_Method);
				udEle_GetDamp_GC = GCHandle::Alloc(udEle_GetDamp);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetDamp);
				Ele_GetDamp GetDamp = static_cast<Ele_GetDamp>(ip.ToPointer());

				this->_Ele_GetMass = _Ele_GetMass;
				UDEle_GetMass^ udEle_GetMass = gcnew UDEle_GetMass(this, &ExternalElementWrapper::GetMass_Method);
				udEle_GetMass_GC = GCHandle::Alloc(udEle_GetMass);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetMass);
				Ele_GetMass GetMass = static_cast<Ele_GetMass>(ip.ToPointer());

				this->_Ele_ZeroLoad = _Ele_ZeroLoad;
				UDEle_ZeroLoad^ udEle_ZeroLoad = gcnew UDEle_ZeroLoad(this, &ExternalElementWrapper::ZeroLoad_Method);
				udEle_ZeroLoad_GC = GCHandle::Alloc(udEle_ZeroLoad);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_ZeroLoad);
				Ele_ZeroLoad ZeroLoad = static_cast<Ele_ZeroLoad>(ip.ToPointer());

				this->_Ele_AddLoad = _Ele_AddLoad;
				UDEle_AddLoad^ udEle_AddLoad = gcnew UDEle_AddLoad(this, &ExternalElementWrapper::AddLoad_Method);
				udEle_AddLoad_GC = GCHandle::Alloc(udEle_AddLoad);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_AddLoad);
				Ele_AddLoad AddLoad = static_cast<Ele_AddLoad>(ip.ToPointer());

				this->_Ele_AddInertiaLoadToUnbalance = _Ele_AddInertiaLoadToUnbalance;
				UDEle_AddInertiaLoadToUnbalance^ udEle_AddInertiaLoadToUnbalance = gcnew UDEle_AddInertiaLoadToUnbalance(this, &ExternalElementWrapper::AddInertiaLoadToUnbalance_Method);
				udEle_AddInertiaLoadToUnbalance_GC = GCHandle::Alloc(udEle_AddInertiaLoadToUnbalance);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_AddInertiaLoadToUnbalance);
				Ele_AddInertiaLoadToUnbalance AddInertiaLoadToUnbalance = static_cast<Ele_AddInertiaLoadToUnbalance>(ip.ToPointer());

				this->_Ele_GetResistingForce = _Ele_GetResistingForce;
				UDEle_GetResistingForce^ udEle_GetResistingForce = gcnew UDEle_GetResistingForce(this, &ExternalElementWrapper::GetResistingForce_Method);
				udEle_GetResistingForce_GC = GCHandle::Alloc(udEle_GetResistingForce);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetResistingForce);
				Ele_GetResistingForce GetResistingForce = static_cast<Ele_GetResistingForce>(ip.ToPointer());

				this->_Ele_GetResistingForceIncInertia = _Ele_GetResistingForceIncInertia;
				UDEle_GetResistingForceIncInertia^ udEle_GetResistingForceIncInertia = gcnew UDEle_GetResistingForceIncInertia(this, &ExternalElementWrapper::GetResistingForceIncInertia_Method);
				udEle_GetResistingForceIncInertia_GC = GCHandle::Alloc(udEle_GetResistingForceIncInertia);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetResistingForceIncInertia);
				Ele_GetResistingForceIncInertia GetResistingForceIncInertia = static_cast<Ele_GetResistingForceIncInertia>(ip.ToPointer());

				this->_Ele_Print = _Ele_Print;
				UDEle_Print^ udEle_Print = gcnew UDEle_Print(this, &ExternalElementWrapper::Print_Method);
				udEle_Print_GC = GCHandle::Alloc(udEle_Print);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_Print);
				Ele_Print Print = static_cast<Ele_Print>(ip.ToPointer());

				this->_Ele_SetResponse = _Ele_SetResponse;
				UDEle_SetResponse^ udEle_SetResponse = gcnew UDEle_SetResponse(this, &ExternalElementWrapper::SetResponse_Method);
				udEle_SetResponse_GC = GCHandle::Alloc(udEle_SetResponse);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_SetResponse);
				Ele_SetResponse setResponse = static_cast<Ele_SetResponse>(ip.ToPointer());

				this->_Ele_GetResponse = _Ele_GetResponse;
				UDEle_GetResponse^ udEle_GetResponse = gcnew UDEle_GetResponse(this, &ExternalElementWrapper::GetResponse_Method);
				udEle_GetResponse_GC = GCHandle::Alloc(udEle_GetResponse);
				ip = Marshal::GetFunctionPointerForDelegate(udEle_GetResponse);
				Ele_GetResponse getResponse = static_cast<Ele_GetResponse>(ip.ToPointer());

				((ExternalElement*)_Element)->SetLinks(GetClassType,
					GetNumExternalNodes,
					GetNodePtrs,
					GetExternalNodes,
					GetNumDOF,
					SetDomain,
					CommitState,
					RevertToLastCommit,
					RevertToStart,
					Update,
					GetTangentStiff,
					GetInitialStiff,
					GetDamp,
					GetMass,
					ZeroLoad,
					AddLoad,
					AddInertiaLoadToUnbalance,
					GetResistingForce,
					GetResistingForceIncInertia,
					Print, getResponse, setResponse
				);
			}

			MatrixWrapper^ GetBaseElementDamp() {
				return MatrixWrapper::GetMatrixWrapper(_Element->Element::getDamp());
			}
			 
		internal:

			const char* GetClassType_Method() {
				String^ str = _Ele_GetClassType();
				char* str2 = (char*)(void*)Marshal::StringToHGlobalAnsi(str);
				return str2;
			}

			int GetNumExternalNodes_Method() {
				return _Ele_GetNumExternalNodes();
			}

			ID *GetExternalNodes_Method() {
				IDWrapper^ id = _Ele_GetExternalNodes();
				if (ExternalNodes == 0)
				{
					ExternalNodes = new ID(id->Size());
					for (int i = 0; i < id->Size(); i++)
					{
						ExternalNodes->operator()(i) = id[i];
					}
				}

				//int x = id[0];
				return id->_ID;
			}

			Node** GetNodePtrs_Method() {
				array<NodeWrapper^>^ nodes = _Ele_GetNodePtrs();
				Node** _nodes = new Node*[nodes->Length];
				for (int i = 0; i < nodes->Length; i++)
					_nodes[i] = nodes[i]->_Node;
				return _nodes;
			}

			int GetNumDOF_Method() {
				return _Ele_GetNumDOF();
			}

			void SetDomain_Method(Domain *theDomain) {
				DomainWrapper^ domain = gcnew DomainWrapper();
				domain->SetDomain(theDomain);
				this->_Element->DomainComponent::setDomain(theDomain);
				return _Ele_SetDomain(domain);
			}

			int CommitState_Method() {
				this->_Element->Element::commitState();
				return _Ele_CommitState();
			}

			int RevertToLastCommit_Method() {
				return _Ele_RevertToLastCommit();
			}

			int RevertToStart_Method() {
				return _Ele_RevertToStart();
			}

			int Update_Method() {
				return _Ele_Update();
			}

			Matrix * GetTangentStiff_Method() {
				return _Ele_GetTangentStiff()->_Matrix;
			}

			Matrix * GetInitialStiff_Method() {
				return _Ele_GetInitialStiff()->_Matrix;
			}

			Matrix * GetDamp_Method() {
				return _Ele_GetDamp()->_Matrix;
			}

			

			Matrix * GetMass_Method() {
				return _Ele_GetMass()->_Matrix;
			}

			void ZeroLoad_Method() {
				_Ele_ZeroLoad();
			}

			int AddLoad_Method(ElementalLoad *theLoad, double loadFactor) {
				ElementalLoadWrapper^ eleload = gcnew ElementalLoadWrapper(theLoad->getTag());
				eleload->_ElementalLoad = theLoad;
				eleload->_TaggedObject = theLoad;
				return _Ele_AddLoad(eleload, loadFactor);
			}

			int AddInertiaLoadToUnbalance_Method(const Vector &accel) {
				VectorWrapper^ vec = VectorWrapper::GetVectorWrapper(accel);
				return _Ele_AddInertiaLoadToUnbalance(vec);
			}

			Vector * GetResistingForce_Method()
			{
				return _Ele_GetResistingForce()->_Vector;
			}

			Vector * GetResistingForceIncInertia_Method()
			{
				return _Ele_GetResistingForceIncInertia()->_Vector;
			}

			void Print_Method(int flag) {
				_Ele_Print(flag);
			}

			Response* SetResponse_Method(const char **argv, int argc, OPS_Stream &output) {
				OPS_StreamWrapper^ ops_stream = gcnew OPS_StreamWrapper();
				ops_stream->_OPS_StreamPtr = &output;
				array<String^>^ _argv = gcnew array<String^>(argc);
				for (int i = 0; i < argc; i++)
				{
					_argv[i] = gcnew String(argv[i]);
				}
				ResponseWrapper^ response = _Ele_SetResponse(_argv, ops_stream);
				return response->_Response;
			}

			int GetResponse_Method(int id, Information &info) {
				InformationWrapper^ _info = gcnew InformationWrapper();
				_info->_Information = &info;
				return _Ele_GetResponse(id, _info);
			}

		private:
			GCHandle udEle_GetClassType_GC;
			GCHandle udEle_GetNumExternalNodes_GC;
			GCHandle udEle_GetExternalNodes_GC;
			GCHandle udEle_GetNodePtrs_GC;
			GCHandle udEle_GetNumDOF_GC;
			GCHandle udEle_SetDomain_GC;
			GCHandle udEle_CommitState_GC;
			GCHandle udEle_RevertToLastCommit_GC;
			GCHandle udEle_RevertToStart_GC;
			GCHandle udEle_Update_GC;
			GCHandle udEle_GetTangentStiff_GC;
			GCHandle udEle_GetInitialStiff_GC;
			GCHandle udEle_GetDamp_GC;
			GCHandle udEle_GetMass_GC;
			GCHandle udEle_ZeroLoad_GC;
			GCHandle udEle_AddLoad_GC;
			GCHandle udEle_AddInertiaLoadToUnbalance_GC;
			GCHandle udEle_GetResistingForce_GC;
			GCHandle udEle_GetResistingForceIncInertia_GC;
			GCHandle udEle_Print_GC;
			GCHandle udEle_SetResponse_GC;
			GCHandle udEle_GetResponse_GC;




			DEle_GetClassType^ _Ele_GetClassType;
			DEle_GetNumExternalNodes^ _Ele_GetNumExternalNodes;
			DEle_GetNodePtrs^ _Ele_GetNodePtrs;
			DEle_GetExternalNodes^ _Ele_GetExternalNodes;
			DEle_GetNumDOF^ _Ele_GetNumDOF;
			DEle_SetDomain^ _Ele_SetDomain;
			DEle_CommitState^ _Ele_CommitState;
			DEle_RevertToLastCommit^ _Ele_RevertToLastCommit;
			DEle_RevertToStart^ _Ele_RevertToStart;
			DEle_Update^ _Ele_Update;
			DEle_GetTangentStiff^ _Ele_GetTangentStiff;
			DEle_GetInitialStiff^ _Ele_GetInitialStiff;
			DEle_GetDamp^ _Ele_GetDamp;
			DEle_GetMass^ _Ele_GetMass;
			DEle_ZeroLoad^ _Ele_ZeroLoad;
			DEle_AddLoad^ _Ele_AddLoad;
			DEle_AddInertiaLoadToUnbalance^ _Ele_AddInertiaLoadToUnbalance;
			DEle_GetResistingForce^ _Ele_GetResistingForce;
			DEle_GetResistingForceIncInertia^ _Ele_GetResistingForceIncInertia;
			DEle_Print^ _Ele_Print;
			DEle_SetResponse^ _Ele_SetResponse;
			DEle_GetResponse^ _Ele_GetResponse;
			ID* ExternalNodes;
		};
	}
}
