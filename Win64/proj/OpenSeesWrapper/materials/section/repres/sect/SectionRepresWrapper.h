#pragma once
#include "../../../../taggeds/TaggedObjectWrapper.h"
#include <FiberSectionRepr.h>


#include "../patch/PatchWrapper.h"
#include "../reinfLayer/ReinfLayerWrapper.h"


namespace OpenSees {
	namespace Materials {
		namespace Sections {
			namespace Repres {
				public ref class SectionRepresWrapper abstract : TaggedObjectWrapper
				{
				public:
					SectionRepresWrapper();
					~SectionRepresWrapper();
				internal:
					FiberSectionRepr * _FiberSectionRepr;
				};

				public ref class FiberSectionReprWrapper : SectionRepresWrapper
				{
				public:
					FiberSectionReprWrapper(int sectionID, array<PatchWrapper^>^ patches,
						array<ReinfLayerWrapper^>^ reinfLayers);
					~FiberSectionReprWrapper();

				internal:

				};
			}
		}
	}
}
