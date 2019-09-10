#pragma once
#include <Patch.h>
#include <QuadPatch.h>
#include <CircPatch.h>

#include "../../../../matrix/VectorWrapper.h"
#include "../../../../matrix/MatrixWrapper.h"

using namespace OpenSees;
namespace OpenSees {
	namespace Materials {
		namespace Sections {
			namespace Repres {
				public ref class  PatchWrapper abstract
				{
				public:
					PatchWrapper();
					~PatchWrapper() {
						if (_Patch != 0)
							delete _Patch;
					};
				internal:
					Patch * _Patch;
				};

				public ref class QuadPatchWrapper : PatchWrapper
				{
				public:
					QuadPatchWrapper(int materialID, int numSubdivIJ, int numSubdivJK,
						MatrixWrapper^ vertexCoords);
					~QuadPatchWrapper() {
						if (_Patch != 0)
							delete _Patch;
					};
				internal:

				};
				public ref class CircPatchWrapper : PatchWrapper
				{
				public:
					CircPatchWrapper(int materialID, int numSubdivCircunf, int numSubdivRadial,
						VectorWrapper^ centerPosition, double internRadius,
						double externRadius, double initialAngle, double finalAngle);
					~CircPatchWrapper() {
						if (_Patch != 0)
							delete _Patch;
					};
				internal:

				};
			}
		}
	}
}
