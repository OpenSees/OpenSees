#include "stdafx.h"

#include "RecorderWrapper.h"

using namespace System::Runtime::InteropServices;
using namespace System;
using namespace OpenSees;
using namespace OpenSees::DamageModels;
using namespace OpenSees::Recorders;
using namespace OpenSees::Components;
using namespace OpenSees::Handlers;


DamageRecorderWrapper::DamageRecorderWrapper(int elemid, IDWrapper^ secIDs, int dofid, DamageModelWrapper^ dmgPtr,
	BaseDomainWrapper^ theDomainPtr,
	bool echotimeflag, double deltat, OPS_StreamWrapper^ theOutputStream) {
	_Recorder = new DamageRecorder(elemid, *secIDs->_ID, dofid, dmgPtr->_DamageModel, *theDomainPtr->_Domain, echotimeflag, deltat, *theOutputStream->_OPS_StreamPtr);
}

NodeRecorderWrapper::NodeRecorderWrapper(IDWrapper^ theDof,
	IDWrapper^ theNodes,
	int sensitivity,
	String^ dataToStore,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theOutputHandler) {
	_Recorder = new NodeRecorder(*theDof->_ID, theNodes->_ID, sensitivity, (char*)(void*)Marshal::StringToHGlobalAnsi(dataToStore), *theDomain->_Domain, *theOutputHandler->_OPS_StreamPtr);
}

DriftRecorderWrapper::DriftRecorderWrapper(int ndI, int ndJ, int dof, int perpDirn,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler,
	bool echoTime,
	double dT) {
	_Recorder = new DriftRecorder(ndI, ndJ, dof, perpDirn, *theDomain->_Domain, *theHandler->_OPS_StreamPtr, echoTime, dT);
}


DriftRecorderWrapper::DriftRecorderWrapper(int ndI, int ndJ, int dof, int perpDirn,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler) {
	_Recorder = new

		DriftRecorder(ndI, ndJ, dof, perpDirn, *theDomain->_Domain, *theHandler->_OPS_StreamPtr);
}


DriftRecorderWrapper::DriftRecorderWrapper(IDWrapper^ ndI, IDWrapper^ ndJ, int dof, int perpDirn,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler,
	bool echoTime,
	double dT) {
	_Recorder = new

		DriftRecorder(*ndI->_ID, *ndJ->_ID, dof, perpDirn, *theDomain->_Domain, *theHandler->_OPS_StreamPtr, echoTime, dT);
}


DriftRecorderWrapper::DriftRecorderWrapper(IDWrapper^ ndI, IDWrapper^ ndJ, int dof, int perpDirn,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler) {
	_Recorder = new

		DriftRecorder(*ndI->_ID, *ndJ->_ID, dof, perpDirn, *theDomain->_Domain, *theHandler->_OPS_StreamPtr);
}

ElementRecorderWrapper::ElementRecorderWrapper(IDWrapper^ eleID,
	array<String^>^ argv,
	bool echoTime,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theOutputHandler,
	double deltaT,
	IDWrapper^ dof) {


	const char** cargv = new const char*[argv->Length];
	for (int i = 0; i < argv->Length; i++)
	{
		cargv[i] = (char*)(void*)Marshal::StringToHGlobalAnsi(argv[i]);
	}

	_Recorder = new ElementRecorder(eleID->_ID, cargv,
		argv->Length, echoTime,
		*theDomain->_Domain, *theOutputHandler->_OPS_StreamPtr, deltaT, dof->_ID);
}

ElementRecorderWrapper::ElementRecorderWrapper(IDWrapper^ eleID,
	array<String^>^ argv,
	bool echoTime,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theOutputHandler) {

	const char** cargv = new const char*[argv->Length];
	for (int i = 0; i < argv->Length; i++)
	{
		cargv[i] = (char*)(void*)Marshal::StringToHGlobalAnsi(argv[i]);
	}

	_Recorder = new		ElementRecorder(eleID->_ID, cargv, argv->Length, echoTime, *theDomain->_Domain, *theOutputHandler->_OPS_StreamPtr);
}

EnvelopeDriftRecorderWrapper::EnvelopeDriftRecorderWrapper(int ndI, int ndJ, int dof, int perpDirn,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler,
	bool echoTime) {
	_Recorder = new EnvelopeDriftRecorder(ndI, ndJ, dof, perpDirn, *theDomain->_Domain, *theHandler->_OPS_StreamPtr, echoTime);
}

EnvelopeDriftRecorderWrapper::
EnvelopeDriftRecorderWrapper(IDWrapper^ ndI, IDWrapper^ ndJ, int dof, int perpDirn,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler,
	bool echoTime) {
	_Recorder = new
		EnvelopeDriftRecorder(*ndI->_ID, *ndJ->_ID, dof, perpDirn, *theDomain->_Domain, *theHandler->_OPS_StreamPtr, echoTime);
}

EnvelopeElementRecorderWrapper::EnvelopeElementRecorderWrapper(IDWrapper^ eleID,
	array<String^>^ argv,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler,
	double deltaT,
	bool echoTimeFlag,
	IDWrapper^ dof) {
	const char** cargv = new const char*[argv->Length];
	for (int i = 0; i < argv->Length; i++)
	{
		cargv[i] = (char*)(void*)Marshal::StringToHGlobalAnsi(argv[i]);
	}
	_Recorder = new EnvelopeElementRecorder(eleID->_ID, cargv, argv->Length, *theDomain->_Domain, *theHandler->_OPS_StreamPtr, deltaT, echoTimeFlag, dof->_ID);
}


EnvelopeElementRecorderWrapper::EnvelopeElementRecorderWrapper(IDWrapper^ eleID,
	array<String^>^ argv,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler) {

	const char** cargv = new const char*[argv->Length];
	for (int i = 0; i < argv->Length; i++)
	{
		cargv[i] = (char*)(void*)Marshal::StringToHGlobalAnsi(argv[i]);
	}

	_Recorder = new		EnvelopeElementRecorder(eleID->_ID, cargv, argv->Length, *theDomain->_Domain, *theHandler->_OPS_StreamPtr);
}

EnvelopeNodeRecorderWrapper::EnvelopeNodeRecorderWrapper(IDWrapper^ theDof,
	IDWrapper^ theNodes,
	String^ dataToStore,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler,
	double deltaT,
	bool echoTimeFlag,
	array<TimeSeriesWrapper^>^ theTimeSeries) {

	TimeSeries** timeseties = new TimeSeries*[theTimeSeries->Length];
	for (int i = 0; i < theTimeSeries->Length; i++)
	{
		timeseties[i] = theTimeSeries[i]->_TimeSeries;
	}

	_Recorder = new EnvelopeNodeRecorder(*theDof->_ID, theNodes->_ID, (char*)(void*)Marshal::StringToHGlobalAnsi(dataToStore), *theDomain->_Domain, *theHandler->_OPS_StreamPtr, deltaT, echoTimeFlag, timeseties);
}

EnvelopeNodeRecorderWrapper::
EnvelopeNodeRecorderWrapper(IDWrapper^ theDof,
	IDWrapper^ theNodes,
	String^ dataToStore,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler) {
	_Recorder = new
		EnvelopeNodeRecorder(*theDof->_ID, theNodes->_ID, (char*)(void*)Marshal::StringToHGlobalAnsi(dataToStore), *theDomain->_Domain, *theHandler->_OPS_StreamPtr);
}

MaxNodeDispRecorderWrapper::MaxNodeDispRecorderWrapper(int dof, IDWrapper^ theNodes, BaseDomainWrapper^ theDomain) {
	_Recorder = new MaxNodeDispRecorder(dof, *theNodes->_ID, *theDomain->_Domain);
}

NormElementRecorderWrapper::NormElementRecorderWrapper(IDWrapper^ eleID,
	array<String^>^ argv,
	bool echoTimeFlag,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler,
	double deltaT,
	IDWrapper^ dof) {

	const char** cargv = new const char*[argv->Length];
	for (int i = 0; i < argv->Length; i++)
	{
		cargv[i] = (char*)(void*)Marshal::StringToHGlobalAnsi(argv[i]);
	}

	_Recorder = new NormElementRecorder(eleID->_ID, cargv, argv->Length, echoTimeFlag, *theDomain->_Domain, *theHandler->_OPS_StreamPtr, deltaT, dof->_ID);
}


NormElementRecorderWrapper::NormElementRecorderWrapper(IDWrapper^ eleID,
	array<String^>^ argv,
	bool echoTimeFlag,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler) {

	const char** cargv = new const char*[argv->Length];
	for (int i = 0; i < argv->Length; i++)
	{
		cargv[i] = (char*)(void*)Marshal::StringToHGlobalAnsi(argv[i]);
	}

	_Recorder = new		NormElementRecorder(eleID->_ID, cargv, argv->Length, echoTimeFlag, *theDomain->_Domain, *theHandler->_OPS_StreamPtr);
}


NormEnvelopeElementRecorderWrapper::NormEnvelopeElementRecorderWrapper(IDWrapper^ eleID,
	array<String^>^ argv,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler,
	double deltaT,
	bool echoTimeFlag,
	IDWrapper^ dof) {

	const char** cargv = new const char*[argv->Length];
	for (int i = 0; i < argv->Length; i++)
	{
		cargv[i] = (char*)(void*)Marshal::StringToHGlobalAnsi(argv[i]);
	}

	_Recorder = new NormEnvelopeElementRecorder(eleID->_ID, cargv, argv->Length, *theDomain->_Domain, *theHandler->_OPS_StreamPtr, deltaT, echoTimeFlag, dof->_ID);
}


NormEnvelopeElementRecorderWrapper::NormEnvelopeElementRecorderWrapper(IDWrapper^ eleID,
	array<String^>^ argv,
	BaseDomainWrapper^ theDomain,
	OPS_StreamWrapper^ theHandler) {

	const char** cargv = new const char*[argv->Length];
	for (int i = 0; i < argv->Length; i++)
	{
		cargv[i] = (char*)(void*)Marshal::StringToHGlobalAnsi(argv[i]);
	}
	_Recorder = new		NormEnvelopeElementRecorder(eleID->_ID, cargv, argv->Length, *theDomain->_Domain, *theHandler->_OPS_StreamPtr);
}

PatternRecorderWrapper::PatternRecorderWrapper(int thePattern,
	BaseDomainWrapper^ theDomain,
	String^ argv,
	double deltaT,
	int startFlag) {
	_Recorder = new PatternRecorder(thePattern, *theDomain->_Domain, (char*)(void*)Marshal::StringToHGlobalAnsi(argv), deltaT, startFlag);
}








