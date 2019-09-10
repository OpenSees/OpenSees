#include "stdafx.h"
#include "TclWrapper.h"

using namespace OpenSees::Tcl;
using namespace System::Runtime::InteropServices;
TclWrapper::TclWrapper()
{
}

int
TclWrapper::Init() {
	this->interp = Tcl_CreateInterp();
	int ret = OpenSeesAppInit(this->interp);
	if (ret != 0)
		return ret;
	else
		return 0;
}

void
TclWrapper::TclEvalFile(String^ filename)
{
	char* _filename = (char*)(void*)Marshal::StringToHGlobalAnsi(filename);
	Tcl_EvalFile(this->interp, _filename);
	return;
}

void
TclWrapper::TclEval(String^ command)
{
	char* _command = (char*)(void*)Marshal::StringToHGlobalAnsi(command);
	Tcl_Eval(this->interp, _command);
	return;
}

DomainWrapper^
TclWrapper::GetActiveDomain() {
	//Domain* domain = ops_TheActiveDomain;
	DomainWrapper^ dw = gcnew DomainWrapper();
	dw->_Domain = ops_TheActiveDomain;
	return dw;
}

TclWrapper::~TclWrapper()
{
	if (this->interp != 0)
	{
		Tcl_Eval(this->interp, "quit");
		delete this->interp;
	}
}



