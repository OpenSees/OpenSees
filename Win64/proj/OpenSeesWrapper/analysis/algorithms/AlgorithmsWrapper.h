#pragma once
#include <EquiSolnAlgo.h>
#include "AcceleratorWrapper.h"
#include "LineSearchWrapper.h"
#include "../../convergances/ConvergenceTestWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"

using namespace System;
using namespace System::ComponentModel;
using namespace OpenSees::ConvergenceTests;
using namespace OpenSees::Algorithms::Accelerators;
using namespace OpenSees::Algorithms::LineSearchs;
namespace OpenSees {
	namespace Algorithms {

		public ref class SolutionAlgorithmWrapper abstract : IMovableObjectWrapper
		{
		public:
			SolutionAlgorithmWrapper();
		};

		public ref class EquiSolnAlgoWrapper abstract :SolutionAlgorithmWrapper
		{
		public:
			EquiSolnAlgoWrapper();
			void SetConvergenceTest(ConvergenceTestWrapper^ test) {
				_EquiSolnAlgo->setConvergenceTest(test->_ConvergenceTest);
			}

			double GetTotalTimeCPU(){
				return _EquiSolnAlgo->getAccelTimeCPU();
			}

			double GetTotalTimeReal() {
				return _EquiSolnAlgo->getTotalTimeReal();
			}

			double GetSolveTimeCPU() {
				return _EquiSolnAlgo->getSolveTimeCPU();
			}

			double GetSolveTimeReal() {
				return _EquiSolnAlgo->getSolveTimeReal();
			}

			double GetAccelTimeReal() {
				return _EquiSolnAlgo->getAccelTimeReal();
			}

			int GetNumIterations() {
				return _EquiSolnAlgo->getNumIterations();
			}

			int GetNumFactorizations() {
				return _EquiSolnAlgo->getNumFactorizations();
			}


			~EquiSolnAlgoWrapper();
		internal:

			EquiSolnAlgo * _EquiSolnAlgo;
		};

		public ref class LinearWrapper :EquiSolnAlgoWrapper
		{
		public:
			LinearWrapper();
			LinearWrapper(int theTangent, int factorOnce);
			~LinearWrapper();

		internal:
		};

		public ref class NewtonRaphsonWrapper :EquiSolnAlgoWrapper
		{
		public:
			NewtonRaphsonWrapper();
			NewtonRaphsonWrapper(int theTangentToUse, double iFact, double cFact);
			~NewtonRaphsonWrapper();
		internal:

		private:

		};

		public ref class AcceleratedNewtonWrapper :EquiSolnAlgoWrapper
		{
		public:
			AcceleratedNewtonWrapper();
			AcceleratedNewtonWrapper(int tangent);
			AcceleratedNewtonWrapper(ConvergenceTestWrapper^ theTest, AcceleratorWrapper^ theAccel, int tangent);
			~AcceleratedNewtonWrapper();
		internal:

		};

		public ref class BFGSWrapper :EquiSolnAlgoWrapper
		{
		public:
			BFGSWrapper();
			BFGSWrapper(int tangent, int n);
			BFGSWrapper(ConvergenceTestWrapper^ theTest, int tangent, int n);
			~BFGSWrapper();
		internal:

		};

		public ref class BroydenWrapper :EquiSolnAlgoWrapper
		{
		public:
			BroydenWrapper();
			BroydenWrapper(int tangent, int n);
			BroydenWrapper(ConvergenceTestWrapper^ theTest, int tangent, int n);
			~BroydenWrapper();
		internal:

		};

		public ref class KrylovNewtonWrapper :EquiSolnAlgoWrapper
		{
		public:
			KrylovNewtonWrapper();
			KrylovNewtonWrapper(int tangent, int maxDim);
			KrylovNewtonWrapper(ConvergenceTestWrapper^ theTest, int tangent, int maxDim);
			~KrylovNewtonWrapper();
		internal:

		};

		public ref class ModifiedNewtonWrapper :EquiSolnAlgoWrapper
		{
		public:
			ModifiedNewtonWrapper(int tangent);
			ModifiedNewtonWrapper(int tangent, double iFact, double cFact);
			ModifiedNewtonWrapper(ConvergenceTestWrapper^ theTest, int tangent, double iFact, double cFact);
			~ModifiedNewtonWrapper();
		internal:

		};

		public ref class NewtonHallMWrapper :EquiSolnAlgoWrapper
		{
		public:
			NewtonHallMWrapper();
			NewtonHallMWrapper(double initFactor, int method, double alpha, double c);
			~NewtonHallMWrapper();
		internal:

		};

		public ref class NewtonLineSearchWrapper :EquiSolnAlgoWrapper
		{
		public:
			NewtonLineSearchWrapper();
			NewtonLineSearchWrapper(ConvergenceTestWrapper^ theTest, LineSearchWrapper^ theLineSearch);
			~NewtonLineSearchWrapper();
		internal:

		};


		public ref class PeriodicNewtonWrapper :EquiSolnAlgoWrapper
		{
		public:
			PeriodicNewtonWrapper();
			PeriodicNewtonWrapper(int tangent, int maxCount);
			[DescriptionAttribute("tangent = 0, maxCount =3")]
			PeriodicNewtonWrapper(ConvergenceTestWrapper^ theTest, int tangent, int maxCount);
			~PeriodicNewtonWrapper();
		internal:

		};
	}
}