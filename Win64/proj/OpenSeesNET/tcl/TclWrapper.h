#pragma once
#include <tcl.h>
#include <tclDecls.h>
#include <commands.h>
#include <Domain.h>

#include "../domains/domain/DomainWrapper.h"
using namespace System;
namespace OpenSees {
	namespace Tcl {
			public ref class TclWrapper 
			{
			public:
				TclWrapper();
				int Init();
				void TclEvalFile(String^ filename);
				void TclEval(String^ command);
				DomainWrapper^ GetActiveDomain();
				~TclWrapper();

			internal:
				
			private:
				Tcl_Interp* interp;
			};
	}
}
