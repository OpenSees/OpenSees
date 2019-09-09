#pragma once

#include <shell/ShellANDeS.h>
#include <shell/ShellDKGQ.h>
#include <shell/ShellDKGT.h>
#include <shell/ShellMITC4.h>
#include <shell/ShellMITC4Thermal.h>
#include <shell/ShellMITC9.h>
#include <shell/ShellNLDKGQ.h>
#include <shell/ShellNLDKGQThermal.h>
#include <shell/ShellNLDKGT.h>


#include "../BaseElementWrapper.h"
#include "../ElementWrapper.h"
#include "../../materials/ndmaterial/NDMaterialWrapper.h"
#include "../../materials/section/SectionForceDeformationWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials::NDMaterials;
namespace OpenSees {
	namespace Elements {
		public ref class ShellANDeSWrapper : ElementWrapper
		{
		public:
			ShellANDeSWrapper(int element_number,
				int node_numb_1, int node_numb_2, int node_numb_3, double t,
				double E, double nu, double rho);
			~ShellANDeSWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ShellDKGQWrapper : ElementWrapper
		{
		public:
			ShellDKGQWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				SectionForceDeformationWrapper^ theMaterial);
			~ShellDKGQWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		};

		public ref class ShellDKGTWrapper : ElementWrapper
		{
		public:
			ShellDKGTWrapper(int tag,
				int node1,
				int node2,
				int node3,
				SectionForceDeformationWrapper^ theMaterial);
			~ShellDKGTWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		};

		public ref class ShellMITC4Wrapper : ElementWrapper
		{
		public:
			ShellMITC4Wrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				SectionForceDeformationWrapper^ theMaterial,
				bool updateBasis);
			~ShellMITC4Wrapper() {
				if (_Element != 0)
					delete _Element;
			};

		};

		public ref class ShellMITC4ThermalWrapper : ElementWrapper
		{
		public:
			ShellMITC4ThermalWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				SectionForceDeformationWrapper^ theMaterial);
			~ShellMITC4ThermalWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		};

		public ref class ShellMITC9Wrapper : ElementWrapper
		{
		public:
			ShellMITC9Wrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				int node5,
				int node6,
				int node7,
				int node8,
				int node9,
				SectionForceDeformationWrapper^ theMaterial);
			~ShellMITC9Wrapper() {
				if (_Element != 0)
					delete _Element;
			};

		};

		public ref class ShellNLDKGQWrapper : ElementWrapper
		{
		public:
			ShellNLDKGQWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				SectionForceDeformationWrapper^ theMaterial);
			~ShellNLDKGQWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		};

		public ref class ShellNLDKGQThermalWrapper : ElementWrapper
		{
		public:
			ShellNLDKGQThermalWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				SectionForceDeformationWrapper^ theMaterial);
			~ShellNLDKGQThermalWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		};

		public ref class ShellNLDKGTWrapper : ElementWrapper
		{
		public:
			ShellNLDKGTWrapper(int tag,
				int node1,
				int node2,
				int node3,
				SectionForceDeformationWrapper^ theMaterial);
			~ShellNLDKGTWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		};
	}
}
