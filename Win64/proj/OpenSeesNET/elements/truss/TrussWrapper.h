#pragma once

#include <CorotTruss.h>
#include <CorotTruss2.h>
#include <CorotTrussSection.h>
#include <N4BiaxialTruss.h>
#include <Truss.h>
#include <Truss2.h>
#include <TrussSection.h>

#include "../../materials/uniaxial/UniaxialMaterialWrapper.h"
#include "../../materials/section/SectionForceDeformationWrapper.h"
#include "../ElementWrapper.h"

using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::Sections;
namespace OpenSees {
	namespace Elements {
		public ref class CorotTrussWrapper : ElementWrapper
		{
		public:
			CorotTrussWrapper(int tag, int dim,
				int Nd1, int Nd2,
				UniaxialMaterialWrapper^ theMaterial,
				double A, double rho,
				int doRayleighDamping,
				int cMass);
			CorotTrussWrapper() {};
			~CorotTrussWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		private:

		};

		public ref class CorotTruss2Wrapper : ElementWrapper
		{
		public:
			CorotTruss2Wrapper(int tag, int dim,
				int Nd1, int Nd2, int oNd1, int oNd2,
				UniaxialMaterialWrapper^ theMaterial,
				double A, double rho);
			CorotTruss2Wrapper() {};
			~CorotTruss2Wrapper() {
				if (_Element != 0)
					delete _Element;
			};
		private:

		};

		public ref class CorotTrussSectionWrapper : ElementWrapper
		{
		public:
			CorotTrussSectionWrapper(int tag, int dim,
				int Nd1, int Nd2,
				SectionForceDeformationWrapper^ theMaterial,
				double rho,
				int doRayleighDamping,
				int cMass);
			CorotTrussSectionWrapper() {};
			~CorotTrussSectionWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		private:

		};

		public ref class N4BiaxialTrussWrapper : ElementWrapper
		{
		public:
			N4BiaxialTrussWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				int GNd1, int GNd2,
				UniaxialMaterialWrapper^ theMaterial,
				double A,
				double rho ,
				int doRayleighDamping);
			N4BiaxialTrussWrapper() {};
			~N4BiaxialTrussWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		private:

		};

		public ref class TrussWrapper : ElementWrapper
		{
		public:
			TrussWrapper(int tag, int dimension, int Nd1, int Nd2, UniaxialMaterialWrapper^ theMaterial, double A, double rho, int doRayleighDamping,
				int cMass);
			~TrussWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		private:

		};

		public ref class Truss2Wrapper : ElementWrapper
		{
		public:
			Truss2Wrapper(int tag,
				int dimension,
				int Nd1, int Nd2, int oNd1, int oNd2,
				UniaxialMaterialWrapper^ theMaterial,
				double A,
				double rho ,
				int doRayleighDamping );
			~Truss2Wrapper() {
				if (_Element != 0)
					delete _Element;
			};
		private:

		};

		public ref class TrussSectionWrapper : ElementWrapper
		{
		public:
			TrussSectionWrapper(int tag, int dimension,
				int Nd1, int Nd2,
				SectionForceDeformationWrapper^ theSection,
				double rho ,
				int doRayleighDamping ,
				int cMass );
			TrussSectionWrapper() {};
			~TrussSectionWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		private:

		};
	}
}
