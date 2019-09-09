#pragma once

#include <ProfileSPDLinSolver.h>
#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <DistributedProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <ProfileSPDLinSubstrSolver.h>
#include <SProfileSPDLinSOE.h>
#include <SProfileSPDLinSolver.h>

#include "../actors/IMovableObjectWrapper.h"
#include "../analysis/models/AnalysisModelWrapper.h"

#include "LinearSOESolverWrapper.h"
#include "LinearSOEWrapper.h"

using namespace OpenSees;
namespace OpenSees {
	namespace Systems {
		namespace Linears {

			public ref class ProfileSPDLinSolverWrapper abstract : LinearSOESolverWrapper
			{
			public:
				ProfileSPDLinSolverWrapper() {};
				~ProfileSPDLinSolverWrapper() {
					if (_ProfileSPDLinSolver != 0)
						delete _ProfileSPDLinSolver;
				};
			internal:
				ProfileSPDLinSolver * _ProfileSPDLinSolver;
			private:

			};

			public ref class SProfileSPDLinSolverWrapper : LinearSOESolverWrapper
			{
			public:
				SProfileSPDLinSolverWrapper();
				~SProfileSPDLinSolverWrapper() {
					if (_SProfileSPDLinSolver != 0)
						delete _SProfileSPDLinSolver;
				};
			internal:
				SProfileSPDLinSolver* _SProfileSPDLinSolver;
			private:

			};

			public ref class ProfileSPDLinDirectSolverWrapper : ProfileSPDLinSolverWrapper
			{
			public:
				ProfileSPDLinDirectSolverWrapper();
				~ProfileSPDLinDirectSolverWrapper() {
					if (_ProfileSPDLinSolver != 0)
						delete _ProfileSPDLinSolver;
				};
			internal:

			private:

			};

			public ref class ProfileSPDLinSubstrSolverWrapper : ProfileSPDLinSolverWrapper
			{
			public:
				ProfileSPDLinSubstrSolverWrapper();
				~ProfileSPDLinSubstrSolverWrapper() {
					if (_ProfileSPDLinSolver != 0)
						delete _ProfileSPDLinSolver;
				};
			internal:

			private:

			};

			public ref class ProfileSPDLinSOEWrapper : LinearSOEWrapper
			{
			public:
				ProfileSPDLinSOEWrapper() {};
				ProfileSPDLinSOEWrapper(ProfileSPDLinSolverWrapper^ solver);
				~ProfileSPDLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};

			public ref class SProfileSPDLinSOEWrapper : LinearSOEWrapper
			{
			public:
				SProfileSPDLinSOEWrapper() {};
				SProfileSPDLinSOEWrapper(SProfileSPDLinSolverWrapper^ solver);
				~SProfileSPDLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};
		
			public ref class DistributedProfileSPDLinSOEWrapper : ProfileSPDLinSOEWrapper
			{
			public:
				DistributedProfileSPDLinSOEWrapper(ProfileSPDLinSolverWrapper^ solver);
				~DistributedProfileSPDLinSOEWrapper() {
					if (_LinearSOE != 0)
						delete _LinearSOE;
				};
			internal:

			private:

			};
		}
	}
}
