#pragma once
#include <Load.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>

#include <Beam2dPartialUniformLoad.h>
#include <Beam2dPointLoad.h>
#include <ElementalLoad.h>
#include <Beam2dTempLoad.h>
#include <Beam2dThermalAction.h>
#include <Beam2dUniformLoad.h>
#include <Beam3dPointLoad.h>
#include <Beam3dThermalAction.h>
#include <Beam3dUniformLoad.h>
#include <BrickSelfWeight.h>
#include <ShellThermalAction.h>
#include <NodalThermalAction.h>
#include <SelfWeight.h>
#include <SurfaceLoader.h>
#include <ThermalActionWrapper.h>

#include "../../matrix/VectorWrapper.h"
#include "../components/DomainComponentWrapper.h"
#include "../timeSeries/TimeSeriesWrapper.h"

using namespace OpenSees::Components::Timeseries;
namespace OpenSees {
	namespace Components {
		namespace Loads {
			public ref class LoadWrapper abstract : DomainComponentWrapper
			{
			public:
				LoadWrapper(int tag);
				~LoadWrapper();
			internal:
				Load * _Load;
			private:

			};

			public ref class NodalLoadWrapper : LoadWrapper
			{
			public:
				NodalLoadWrapper() :LoadWrapper(0) {};
				NodalLoadWrapper(int tag) : LoadWrapper(tag) {};
				NodalLoadWrapper(int tag, int node, VectorWrapper^ values, bool isLoadConstant);
				int GetNodeTag() {
					return _NodalLoad->getNodeTag();
				}
				~NodalLoadWrapper() {
					if (_NodalLoad != 0)
						delete _NodalLoad;
				};
			internal:
				NodalLoad * _NodalLoad;
			private:

			};

			public ref class ElementalLoadWrapper : LoadWrapper
			{
			public:
				ElementalLoadWrapper() :ElementalLoadWrapper(0) {};
				ElementalLoadWrapper(int tag) : LoadWrapper(tag) {};
				~ElementalLoadWrapper() {};
			internal:
				ElementalLoad * _ElementalLoad;
			private:

			};

			public ref class Beam2dPartialUniformLoadWrapper : ElementalLoadWrapper
			{
			public:
				Beam2dPartialUniformLoadWrapper(int tag,double wTrans, double wAxial, double aL, double bL, int eleTag);
				~Beam2dPartialUniformLoadWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class Beam2dPointLoadWrapper : ElementalLoadWrapper
			{
			public:
				Beam2dPointLoadWrapper(int tag, double Pt, double x, int eleTag, double Pa);
				~Beam2dPointLoadWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class Beam2dTempLoadWrapper : ElementalLoadWrapper
			{
			public:
				Beam2dTempLoadWrapper(int tag, double Ttop1, double Tbot1, double Ttop2,
					double Tbot2,
					int theElementTag);
				Beam2dTempLoadWrapper(int tag, double Tuniform,
					int theElementTag);
				Beam2dTempLoadWrapper(int tag, double Ttop, double Tbot,
					int theElementTag);
				Beam2dTempLoadWrapper(int tag, int theElementTag);
				~Beam2dTempLoadWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class Beam2dThermalActionWrapper : ElementalLoadWrapper
			{
			public:
				Beam2dThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
					double t3, double locY3, double t4, double locY4,
					double t5, double locY5, double t6, double locY6,
					double t7, double locY7, double t8, double locY8,
					double t9, double locY9,
					int theElementTag);
				Beam2dThermalActionWrapper(int tag, double locY1, double locY2,
					TimeSeriesWrapper^ theSeries, int theElementTag);
				Beam2dThermalActionWrapper(int tag, VectorWrapper^ locs,
					TimeSeriesWrapper^ theSeries, int theElementTag);
				Beam2dThermalActionWrapper(int tag, int theElementTag);
				~Beam2dThermalActionWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class Beam2dUniformLoadWrapper : ElementalLoadWrapper
			{
			public:
				Beam2dUniformLoadWrapper(int tag, double wTrans, double wAxial, int eleTag);
				~Beam2dUniformLoadWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class Beam3dPointLoadWrapper : ElementalLoadWrapper
			{
			public:
				Beam3dPointLoadWrapper(int tag, double Py, double Pz, double x, int eleTag, double Pa);
				~Beam3dPointLoadWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class Beam3dThermalActionWrapper : ElementalLoadWrapper
			{
			public:
				Beam3dThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
					double t3, double locY3, double t4, double locY4,
					double t5, double locY5, double t6, double t7, double locZ1,
					double t8, double t9, double locZ2, double t10, double t11, double locZ3,
					double t12, double t13, double locZ4, double t14, double t15, double locZ5,
					int theElementTag);
				Beam3dThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
					double t3, double locY3, double t4, double locY4,
					double t5, double locY5, double t6, double locY6,
					double t7, double locY7, double t8, double locY8,
					double t9, double locY9,
					int theElementTag);
				Beam3dThermalActionWrapper(int tag, double locY1, double locY2, double locZ1, double Z2,
					TimeSeriesWrapper^ theSeries,
					int theElementTag);

				Beam3dThermalActionWrapper(int tag, VectorWrapper^ locs, TimeSeriesWrapper^ theSeries,
					int theElementTag);

				Beam3dThermalActionWrapper(int tag, int theElementTag);

				~Beam3dThermalActionWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class Beam3dUniformLoadWrapper : ElementalLoadWrapper
			{
			public:
				Beam3dUniformLoadWrapper(int tag, double wy, double wz, double wx, int eleTag);
				~Beam3dUniformLoadWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class BrickSelfWeightWrapper : ElementalLoadWrapper
			{
			public:
				BrickSelfWeightWrapper(int tag, int eleTag);
				~BrickSelfWeightWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class ShellThermalActionWrapper : ElementalLoadWrapper
			{
			public:
				ShellThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
					double t3, double locY3, double t4, double locY4,
					double t5, double locY5, double t6, double locY6,
					double t7, double locY7, double t8, double locY8,
					double t9, double locY9,
					int theElementTag);
				ShellThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
					double t3, double locY3, double t4, double locY4,
					double t5, double locY5, int theElementTag);
				ShellThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
					int theElementTag);
				ShellThermalActionWrapper(int tag, double locY1, double locY2,
					TimeSeriesWrapper^ theSeries, int theElementTag);
				ShellThermalActionWrapper(int tag, int theElementTag);
				~ShellThermalActionWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class NodalThermalActionWrapper : NodalLoadWrapper
			{
			public:
				NodalThermalActionWrapper(int tag) :NodalLoadWrapper(tag) {};
				NodalThermalActionWrapper(int tag, int theNodeTag,
					VectorWrapper^ locy, TimeSeriesWrapper^ theSeries, VectorWrapper^ crds);
				NodalThermalActionWrapper(int tag, int theNodeTag,
					double locY1, double locY2, double locZ1, double locZ2,
					TimeSeriesWrapper^ theSeries, VectorWrapper^ crds);
				NodalThermalActionWrapper(int tag, int theNodeTag,
					double t1, double locY1, double t2, double locY2, VectorWrapper^ crds);
				~NodalThermalActionWrapper() {
					if (_NodalLoad != 0)
						delete _NodalLoad;
				};
			};

			public ref class SelfWeightWrapper : ElementalLoadWrapper
			{
			public:
				SelfWeightWrapper(int tag, double xFact, double yFact, double zFact, int eleTag);
				~SelfWeightWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class SurfaceLoaderWrapper : ElementalLoadWrapper
			{
			public:
				SurfaceLoaderWrapper(int tag, int eleTag);
				~SurfaceLoaderWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};

			public ref class ThermalActionWrapperWrapper : ElementalLoadWrapper
			{
			public:
				ThermalActionWrapperWrapper(int tag, int EleTag,
					NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2);
				ThermalActionWrapperWrapper(int tag, int EleTag,
					NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2,
					NodalThermalActionWrapper^ theNodalTA3);
				ThermalActionWrapperWrapper(int tag, int EleTag,
					NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2,
					NodalThermalActionWrapper^ theNodalTA3, NodalThermalActionWrapper^ theNodalTA4);
				ThermalActionWrapperWrapper(int tag, int EleTag,
					NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2,
					NodalThermalActionWrapper^ theNodalTA3, NodalThermalActionWrapper^ theNodalTA4,
					NodalThermalActionWrapper^ theNodalTA5);
				ThermalActionWrapperWrapper(int tag, int EleTag,
					NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2,
					NodalThermalActionWrapper^ theNodalTA3, NodalThermalActionWrapper^ theNodalTA4,
					NodalThermalActionWrapper^ theNodalTA5, NodalThermalActionWrapper^ theNodalTA6);
				~ThermalActionWrapperWrapper() {
					if (_ElementalLoad != 0)
						delete _ElementalLoad;
				};
			};
		}
	}
}
