#include "stdafx.h"
#include "ExternalElementWrapper.h"

using namespace System::Runtime::InteropServices;
using namespace OpenSees;
using namespace OpenSees::Elements;

ExternalElementWrapper::ExternalElementWrapper(int tag) {
	_Element = new ExternalElement(tag);
}

OpenSees::Elements::ExternalElementWrapper::~ExternalElementWrapper()
{
	if (_Element != 0)
		delete _Element;

	udEle_GetClassType_GC.Free();
	udEle_GetNumExternalNodes_GC.Free();
	udEle_GetExternalNodes_GC.Free();
	udEle_GetNodePtrs_GC.Free();
	udEle_GetNumDOF_GC.Free();
	udEle_SetDomain_GC.Free();
	udEle_CommitState_GC.Free();
	udEle_RevertToLastCommit_GC.Free();
	udEle_RevertToStart_GC.Free();
	udEle_Update_GC.Free();
	udEle_GetTangentStiff_GC.Free();
	udEle_GetInitialStiff_GC.Free();
	udEle_GetDamp_GC.Free();
	udEle_GetMass_GC.Free();
	udEle_ZeroLoad_GC.Free();
	udEle_AddLoad_GC.Free();
	udEle_AddInertiaLoadToUnbalance_GC.Free();
	udEle_GetResistingForce_GC.Free();
	udEle_GetResistingForceIncInertia_GC.Free();
	udEle_Print_GC.Free();
	udEle_SetResponse_GC.Free();
	udEle_GetResponse_GC.Free();

}

