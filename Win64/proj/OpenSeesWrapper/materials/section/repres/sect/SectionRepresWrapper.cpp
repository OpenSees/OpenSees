#include "stdafx.h"
#include "SectionRepresWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Sections::Repres;
SectionRepresWrapper::SectionRepresWrapper()
{

}

SectionRepresWrapper::~SectionRepresWrapper()
{
	if (_FiberSectionRepr != 0)
		delete _FiberSectionRepr;
}

FiberSectionReprWrapper::FiberSectionReprWrapper(int sectionID, array<PatchWrapper^>^ patches,
	array<ReinfLayerWrapper^>^ reinfLayers)
{
	int numPatches = patches->Length;
	int numReinfLayers = reinfLayers->Length;
	Patch** _patchs = new Patch*[numPatches];
	ReinfLayer** _reinfLayers = new ReinfLayer*[numReinfLayers];


	_FiberSectionRepr = new FiberSectionRepr(sectionID, numPatches, numReinfLayers);
	for (int i = 0; i < numPatches; i++)
		_FiberSectionRepr->addPatch(*patches[i]->_Patch);
	for (int i = 0; i < numReinfLayers; i++)
		_FiberSectionRepr->addReinfLayer(*reinfLayers[i]->_ReinfLayer);

}

FiberSectionReprWrapper::~FiberSectionReprWrapper()
{
	if (_FiberSectionRepr != 0)
		delete _FiberSectionRepr;
}
