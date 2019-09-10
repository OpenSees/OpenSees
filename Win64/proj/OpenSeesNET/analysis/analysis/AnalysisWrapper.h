#pragma once

#include <StaticAnalysis.h>
#include <TransientAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <PFEMAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#include "../../domains/domain/DomainWrapper.h"
#include "../../analysis/handlers/ConstraintHandlerWrapper.h"
#include "../../analysis/numberers/DOF_NumbererWrapper.h"
#include "../../analysis/models/AnalysisModelWrapper.h"
#include "../../analysis/algorithms/AlgorithmsWrapper.h"
#include "../../systems/LinearSOEWrapper.h"
#include "../../systems/EigenSOEWrapper.h"
#include "../../analysis/integrators/StaticIntegratorWrapper.h"
#include "../../convergances/ConvergenceTestWrapper.h"
#include "../../analysis/integrators/TransientIntegratorWrapper.h"
using namespace OpenSees::Algorithms;
using namespace OpenSees::Systems;
using namespace OpenSees::Systems::Eigens;
using namespace OpenSees::Systems::Linears;
using namespace OpenSees::Handlers;
using namespace OpenSees::Integrators::Static;
using namespace OpenSees::Integrators::Transient;
using namespace OpenSees::Numberers;
using namespace OpenSees::ConvergenceTests;
namespace OpenSees {
	namespace Analysis {
		public ref class AnalysisWrapper abstract
		{
		public:
			AnalysisWrapper();
			~AnalysisWrapper();
		internal:

		private:

		};

		public ref class StaticAnalysisWrapper : AnalysisWrapper
		{
		public:
			StaticAnalysisWrapper(DomainWrapper^ theDomain,
				ConstraintHandlerWrapper^ theHandler,
				DOF_NumbererWrapper^ theNumberer,
				AnalysisModelWrapper^ theModel,
				EquiSolnAlgoWrapper^ theSolnAlgo,
				LinearSOEWrapper^ theSOE,
				StaticIntegratorWrapper^ theIntegrator,
				ConvergenceTestWrapper^ theTest);
			int Analyze(int steps);
			int SetEigenSOE(EigenSOEWrapper^ theEigenSOE) {
				return _StaticAnalysis->setEigenSOE(*theEigenSOE->_EigenSOE);
			};
			~StaticAnalysisWrapper();
			int SetIntegrator(StaticIntegratorWrapper^ theIntegrator) {
				return _StaticAnalysis->setIntegrator(*theIntegrator->_StaticIntegrator);
			}

			int SetAlgorithm(EquiSolnAlgoWrapper^ theSolnAlgo) {
				return _StaticAnalysis->setAlgorithm(*theSolnAlgo->_EquiSolnAlgo);
			}

			int SetConvergenceTest(ConvergenceTestWrapper^ theTest) {
				return _StaticAnalysis->setConvergenceTest(*theTest->_ConvergenceTest);
			}

			int SetLinearSOE(LinearSOEWrapper^ theSOE) {
				return _StaticAnalysis->setLinearSOE(*theSOE->_LinearSOE);
			}

			int SeNumberer(DOF_NumbererWrapper^ theNumberer) {
				return _StaticAnalysis->setNumberer(*theNumberer->_DOFNumberer);
			}


			int Eigen(int numMode, bool generalized, bool findSmallest) {
				return _StaticAnalysis->eigen(numMode, generalized, findSmallest);
			}

			[ObsoleteAttribute("Use ClearAll method")]
			void Wipe() {
				if (_StaticAnalysis != 0)
					_StaticAnalysis->clearAll();
			}

			void ClearAll() {
				if (_StaticAnalysis != 0)
					_StaticAnalysis->clearAll();
			}

		internal:
			StaticAnalysis * _StaticAnalysis;
		private:

		};

		public ref class TransientAnalysisWrapper abstract : AnalysisWrapper
		{
		public:
			TransientAnalysisWrapper();
			virtual int Analyze(int steps, double dt) {
				throw gcnew System::NotImplementedException();
			};
			~TransientAnalysisWrapper();
		internal:

		private:

		};

		public ref class DirectIntegrationAnalysisWrapper : TransientAnalysisWrapper
		{
		public:
			DirectIntegrationAnalysisWrapper(DomainWrapper^ theDomain,
				ConstraintHandlerWrapper^ theHandler,
				DOF_NumbererWrapper^ theNumberer,
				AnalysisModelWrapper^ theModel,
				EquiSolnAlgoWrapper^ theSolnAlgo,
				LinearSOEWrapper^ theSOE,
				TransientIntegratorWrapper^ theIntegrator,
				ConvergenceTestWrapper^ theTest);
			int SetEigenSOE(EigenSOEWrapper^ theEigenSOE) {
				return _DirectIntegrationAnalysis->setEigenSOE(*theEigenSOE->_EigenSOE);
			};
			virtual int Analyze(int steps, double dt) override;
			int Eigen(int numMode, bool generalized, bool findSmallest) {
				return _DirectIntegrationAnalysis->eigen(numMode, generalized, findSmallest);
			}
			~DirectIntegrationAnalysisWrapper();
			int SetIntegrator(TransientIntegratorWrapper^ theIntegrator) {
				return _DirectIntegrationAnalysis->setIntegrator(*theIntegrator->_TransientIntegrator);
			}

			int SetAlgorithm(EquiSolnAlgoWrapper^ theSolnAlgo) {
				return _DirectIntegrationAnalysis->setAlgorithm(*theSolnAlgo->_EquiSolnAlgo);
			}

			int SetConvergenceTest(ConvergenceTestWrapper^ theTest) {
				return _DirectIntegrationAnalysis->setConvergenceTest(*theTest->_ConvergenceTest);
			}

			int SetLinearSOE(LinearSOEWrapper^ theSOE) {
				return _DirectIntegrationAnalysis->setLinearSOE(*theSOE->_LinearSOE);
			}

			int SeNumberer(DOF_NumbererWrapper^ theNumberer) {
				return _DirectIntegrationAnalysis->setNumberer(*theNumberer->_DOFNumberer);
			}

			[ObsoleteAttribute("Use ClearAll method")]
			void Wipe() {
				if (_DirectIntegrationAnalysis != 0)
					_DirectIntegrationAnalysis->clearAll();
			}

			void ClearAll() {
				if (_DirectIntegrationAnalysis != 0)
					_DirectIntegrationAnalysis->clearAll();
			}

		internal:
			DirectIntegrationAnalysis * _DirectIntegrationAnalysis;
		private:

		};

		public ref class PFEMAnalysisAnalysisWrapper : TransientAnalysisWrapper
		{
		public:
			PFEMAnalysisAnalysisWrapper(DomainWrapper^ theDomain,
				ConstraintHandlerWrapper^ theHandler,
				DOF_NumbererWrapper^ theNumberer,
				AnalysisModelWrapper^ theModel,
				EquiSolnAlgoWrapper^ theSolnAlgo,
				LinearSOEWrapper^ theSOE,
				TransientIntegratorWrapper^ theIntegrator,
				ConvergenceTestWrapper^ theTest,
				double max, double min, double g, double r);

			virtual int Analyze() override;
			void Wipe() {
				if (_PFEMAnalysis != 0)
					_PFEMAnalysis->clearAll();
			}

			~PFEMAnalysisAnalysisWrapper();
		internal:
			PFEMAnalysis * _PFEMAnalysis;
		private:

		};

		public ref class VariableTimeStepDirectIntegrationAnalysisWrapper : TransientAnalysisWrapper
		{
		public:
			VariableTimeStepDirectIntegrationAnalysisWrapper(DomainWrapper^ theDomain,
				ConstraintHandlerWrapper^ theHandler,
				DOF_NumbererWrapper^ theNumberer,
				AnalysisModelWrapper^ theModel,
				EquiSolnAlgoWrapper^ theSolnAlgo,
				LinearSOEWrapper^ theSOE,
				TransientIntegratorWrapper^ theIntegrator,
				ConvergenceTestWrapper^ theTest);

			virtual int Analyze(int numSteps, double dT, double dtMin, double dtMax, int Jd) override;
			void Wipe() {
				if (_VariableTimeStepDirectIntegrationAnalysis != 0)
					_VariableTimeStepDirectIntegrationAnalysis->clearAll();
			}
			~VariableTimeStepDirectIntegrationAnalysisWrapper();
		internal:
			VariableTimeStepDirectIntegrationAnalysis * _VariableTimeStepDirectIntegrationAnalysis;
		private:

		};
	}
}