#pragma once

#include <ConstraintHandler.h>
#include <PlainHandler.h>
#include <TransformationConstraintHandler.h>
#include <LagrangeConstraintHandler.h>
#include <PenaltyConstraintHandler.h>
#include "../../actors/IMovableObjectWrapper.h"




namespace OpenSees {
	namespace Handlers {
		public ref class ConstraintHandlerWrapper abstract : IMovableObjectWrapper
		{
		public:
			ConstraintHandlerWrapper();
			~ConstraintHandlerWrapper();
		internal:
			ConstraintHandler * _ConstraintHandler;
		private:

		};

		public ref class PlainHandlerWrapper : ConstraintHandlerWrapper
		{
		public:
			PlainHandlerWrapper();
			~PlainHandlerWrapper();
		internal:

		private:

		};

		public ref class TransformationConstraintHandlerWrapper : ConstraintHandlerWrapper
		{
		public:
			TransformationConstraintHandlerWrapper();
			~TransformationConstraintHandlerWrapper();
		internal:

		private:

		};

		public ref class LagrangeConstraintHandlerWrapper : ConstraintHandlerWrapper
		{
		public:
			LagrangeConstraintHandlerWrapper(double alphaSP, double alphaMP);
			~LagrangeConstraintHandlerWrapper();
		internal:

		private:

		};

		public ref class PenaltyConstraintHandlerWrapper : ConstraintHandlerWrapper
		{
		public:
			PenaltyConstraintHandlerWrapper(double alphaSP, double alphaMP);
			~PenaltyConstraintHandlerWrapper();
		internal:

		private:

		};
	}
}