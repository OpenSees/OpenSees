#pragma once

#include <BeamIntegration.h>
#include <CompositeSimpsonBeamIntegration.h>
#include <DistHingeIntegration.h>
#include <FixedLocationBeamIntegration.h>
#include <HingeEndpointBeamIntegration.h>
#include <HingeMidpointBeamIntegration.h>
#include <HingeRadauBeamIntegration.h>
#include <HingeRadauTwoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <LobattoBeamIntegration.h>
#include <LowOrderBeamIntegration.h>
#include <MidDistanceBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <RegularizedHingeIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <UserDefinedBeamIntegration.h>
#include <UserDefinedHingeIntegration.h>


#include "../../../actors//IMovableObjectWrapper.h"
#include "../../../matrix/VectorWrapper.h"

using namespace OpenSees;

namespace OpenSees {
	namespace Elements {
		namespace BeamIntegrations {
			public enum class  BeamIntegrationType { 
				Lobatto, 
				Legendre, 
				Radau, 
				CompositeSimpson, 
				NewtonCotes,
				Trapezoidal,
			};

			public ref class BeamIntegrationWrapper : IMovableObjectWrapper
			{
			public:
				BeamIntegrationWrapper() {};
				BeamIntegrationWrapper(BeamIntegrationType^ type);

				VectorWrapper^ GetSectionPositions(int numSections,double l) {
					double* xi = new double[numSections];
					_BeamIntegration->getSectionLocations(numSections, l, xi);
					return VectorWrapper::GetVectorWrapper(new Vector(xi, numSections));
				};

				VectorWrapper^ GetSectionWeights(int numSections, double l) {
					double* wi = new double[numSections];
					_BeamIntegration->getSectionWeights(numSections, l, wi);
					return VectorWrapper::GetVectorWrapper(new Vector(wi, numSections));
				};

				~BeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:
				BeamIntegration * _BeamIntegration;
			private:
				BeamIntegrationType _type;
			};

			public ref class DistHingeIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				DistHingeIntegrationWrapper() {};
				DistHingeIntegrationWrapper(double lpI, double lpJ, BeamIntegrationWrapper^ bi);
				~DistHingeIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:
				
			private:
				
			};

			public ref class FixedLocationBeamIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				FixedLocationBeamIntegrationWrapper() {};
				FixedLocationBeamIntegrationWrapper(int nIP, VectorWrapper^ pt);
				~FixedLocationBeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class HingeEndpointBeamIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				HingeEndpointBeamIntegrationWrapper() {};
				HingeEndpointBeamIntegrationWrapper(double lpI, double lpJ);
				~HingeEndpointBeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class HingeMidpointBeamIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				HingeMidpointBeamIntegrationWrapper() {};
				HingeMidpointBeamIntegrationWrapper(double lpI, double lpJ);
				~HingeMidpointBeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class HingeRadauBeamIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				HingeRadauBeamIntegrationWrapper() {};
				HingeRadauBeamIntegrationWrapper(double lpI, double lpJ);
				~HingeRadauBeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class HingeRadauTwoBeamIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				HingeRadauTwoBeamIntegrationWrapper() {};
				HingeRadauTwoBeamIntegrationWrapper(double lpI, double lpJ);
				~HingeRadauTwoBeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class LowOrderBeamIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				LowOrderBeamIntegrationWrapper() {};
				LowOrderBeamIntegrationWrapper(int nIP, VectorWrapper^ pt, int nc, VectorWrapper^ wt);
				~LowOrderBeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class MidDistanceBeamIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				MidDistanceBeamIntegrationWrapper() {};
				MidDistanceBeamIntegrationWrapper(int nIP, VectorWrapper^ pt);
				~MidDistanceBeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class RegularizedHingeIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				RegularizedHingeIntegrationWrapper() {};
				RegularizedHingeIntegrationWrapper(BeamIntegrationWrapper^ bi,
					double lpI, double lpJ,
					double epsI, double epsJ);
				~RegularizedHingeIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class UserDefinedBeamIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				UserDefinedBeamIntegrationWrapper() {};
				UserDefinedBeamIntegrationWrapper(int nIP, VectorWrapper^ pt, VectorWrapper^ wt);
				~UserDefinedBeamIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};

			public ref class UserDefinedHingeIntegrationWrapper : BeamIntegrationWrapper
			{
			public:
				UserDefinedHingeIntegrationWrapper() {};
				UserDefinedHingeIntegrationWrapper(int npL, VectorWrapper^ ptL, VectorWrapper^ wtL,
					int npR, VectorWrapper^ ptR, VectorWrapper^ wtR);
				~UserDefinedHingeIntegrationWrapper() {
					if (_BeamIntegration != 0)
						delete _BeamIntegration;
				};

			internal:

			private:

			};
		}
	}
}
