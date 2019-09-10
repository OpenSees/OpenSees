#pragma once
#include <ConvergenceTest.h>
#include <CTestEnergyIncr.h>
#include <CTestNormDispIncr.h>
#include <CTestFixedNumIter.h>
#include <CTestNormUnbalance.h>
#include <CTestPFEM.h>
#include <CTestRelativeEnergyIncr.h>
#include <CTestRelativeNormDispIncr.h>
#include <CTestRelativeNormUnbalance.h>
#include <CTestRelativeTotalNormDispIncr.h>
#include <NormDispAndUnbalance.h>
#include <NormDispOrUnbalance.h>

#include "../actors/IMovableObjectWrapper.h"


namespace OpenSees {
	namespace ConvergenceTests {
		public ref class ConvergenceTestWrapper abstract :IMovableObjectWrapper
		{
		public:
			ConvergenceTestWrapper();
			~ConvergenceTestWrapper();
			array<double>^ GetNorms()
			{
				Vector vec = _ConvergenceTest->getNorms();
				int size = vec.Size();
				array<double>^ _vec = gcnew array<double>(size);
				for (int i = 0; i < size; i++)
					_vec[i] = vec[i];

				return _vec;
			}

			int GetNumTests() {
				return _ConvergenceTest->getNumTests();
			}

		internal:
			ConvergenceTest * _ConvergenceTest;
		};



		public ref class CTestEnergyIncrWrapper :ConvergenceTestWrapper
		{
		public:

			CTestEnergyIncrWrapper(double tol, int maxNumIter, int printFlag, int normType);
			~CTestEnergyIncrWrapper();
		private:

		};

		public ref class CTestFixedNumIterWrapper :ConvergenceTestWrapper
		{
		public:

			CTestFixedNumIterWrapper(int maxNumIter, int printFlag, int normType);
			~CTestFixedNumIterWrapper();
		private:

		};

		public ref class CTestNormDispIncrWrapper :ConvergenceTestWrapper
		{
		public:
			//CTestNormDispIncrWrapper();
			CTestNormDispIncrWrapper(double tol, int maxNumIter, int printFlag, int normType, double maxTol);
			~CTestNormDispIncrWrapper();
		private:

		};

		public ref class CTestNormUnbalanceWrapper :ConvergenceTestWrapper
		{
		public:
			//CTestNormUnbalanceWrapper();
			CTestNormUnbalanceWrapper(double tol, int maxNumIter, int printFlag, int normType, int maxincr, double maxTol);
			~CTestNormUnbalanceWrapper();
		private:

		};

		public ref class CTestPFEMWrapper :ConvergenceTestWrapper
		{
		public:
			//CTestPFEMWrapper();
			CTestPFEMWrapper(double tv, double tp, double tv2, double tp2, double tvrel, double tprel,
				int maxNumIter, int maxincr,
				int printFlag, int normType);
			~CTestPFEMWrapper();
		private:

		};

		public ref class CTestRelativeEnergyIncrWrapper :ConvergenceTestWrapper
		{
		public:
			//CTestRelativeEnergyIncrWrapper();
			CTestRelativeEnergyIncrWrapper(double tol, int maxNumIter, int printFlag, int normType);
			~CTestRelativeEnergyIncrWrapper();
		private:

		};

		public ref class CTestRelativeNormDispIncrWrapper :ConvergenceTestWrapper
		{
		public:
			//CTestRelativeNormDispIncrWrapper();
			CTestRelativeNormDispIncrWrapper(double tol, int maxNumIter, int printFlag, int normType);
			~CTestRelativeNormDispIncrWrapper();
		private:

		};

		public ref class CTestRelativeNormUnbalanceWrapper :ConvergenceTestWrapper
		{
		public:
			//CTestRelativeNormUnbalanceWrapper();
			CTestRelativeNormUnbalanceWrapper(double tol, int maxNumIter, int printFlag, int normType);
			~CTestRelativeNormUnbalanceWrapper();
		private:

		};

		public ref class CTestRelativeTotalNormDispIncrWrapper :ConvergenceTestWrapper
		{
		public:
			//CTestRelativeTotalNormDispIncrWrapper();
			CTestRelativeTotalNormDispIncrWrapper(double tol, int maxNumIter, int printFlag, int normType);
			~CTestRelativeTotalNormDispIncrWrapper();
		private:

		};

		public ref class NormDispAndUnbalanceWrapper :ConvergenceTestWrapper
		{
		public:
			//NormDispAndUnbalanceWrapper();
			NormDispAndUnbalanceWrapper(double tolDisp,
				double tolUnbalance,
				int maxNumIter,
				int printFlag,
				int normType, int maxincr);
			~NormDispAndUnbalanceWrapper();
		private:

		};

		public ref class NormDispOrUnbalanceWrapper :ConvergenceTestWrapper
		{
		public:
			//NormDispOrUnbalanceWrapper();
			NormDispOrUnbalanceWrapper(double tolDisp,
				double tolUnbalance,
				int maxNumIter,
				int printFlag,
				int normType, int maxincr);
			~NormDispOrUnbalanceWrapper();
		private:

		};
	}
}
