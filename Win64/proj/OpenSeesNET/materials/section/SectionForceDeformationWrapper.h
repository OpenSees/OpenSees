#pragma once
#include <ID.h>
#include <Fiber.h>
#include <repres\cell\Cell.h>
#include <repres\reinfBar\ReinfBar.h>

#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>
#include <UniaxialMaterial.h>

#include <SectionForceDeformation.h>

#include <dispBeamColumnInt/FiberSection2dInt.h>

#include <yieldSurface\YieldSurfaceSection2d.h>
#include <yieldSurface\YS_Section2D01.h>
#include <yieldSurface\YS_Section2D02.h>

#include <Bidirectional.h>
#include <ElasticMembranePlateSection.h>
#include <ElasticPlateSection.h>
#include <ElasticSection2d.h>
#include <ElasticSection3d.h>
#include <ElasticShearSection2d.h>
#include <ElasticShearSection3d.h>
#include <ElasticTubeSection3d.h>
#include <ElasticWarpingShearSection2d.h>
#include <Elliptical2.h>
#include <ElasticShearSection3d.h>
#include <FiberSection2d.h>
#include <FiberSection2dThermal.h>
#include <FiberSection3d.h>
#include <FiberSection3dThermal.h>
#include <FiberSectionGJ.h>
#include <FiberSectionGJThermal.h>
#include <GenericSection1d.h>
#include <Isolator2spring.h>
#include <LayeredShellFiberSection.h>
#include <LayeredShellFiberSectionThermal.h>
#include <MembranePlateFiberSection.h>
#include <MembranePlateFiberSectionThermal.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>
#include <NDFiberSectionWarping2d.h>
#include <ParallelSection.h>
#include <SectionAggregator.h>




#include "../uniaxial/UniaxialMaterialWrapper.h"
#include "../ndmaterial/NDMaterialWrapper.h"
#include "../yieldSurface/YieldSurfaceWrapper.h"
#include "fiber/FiberWrapper.h"
#include "repres/sect/SectionRepresWrapper.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace OpenSees::Materials;
using namespace OpenSees::Materials::YieldSurfaces;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;
using namespace OpenSees::Materials::Sections::Repres;
using namespace OpenSees::Materials::Sections;

namespace OpenSees {
	namespace Materials {
		namespace Sections {
			public ref class SectionForceDeformationWrapper : MaterialWrapper
			{
			public:
				SectionForceDeformationWrapper();

				VectorWrapper^ GetSectionDeformation() {
					Vector vec = _SectionForceDeformation->getSectionDeformation();
					return VectorWrapper::GetVectorWrapper(&vec);
				}

				VectorWrapper^ GetStressResultant() {
					Vector vec = _SectionForceDeformation->getStressResultant();
					return VectorWrapper::GetVectorWrapper(&vec);
				}

				MatrixWrapper^ GetSectionTangent() {
					Matrix mat = _SectionForceDeformation->getSectionTangent();
					return MatrixWrapper::GetMatrixWrapper(&mat);
				}

				MatrixWrapper^ GetInitialTangent() {
					Matrix mat = _SectionForceDeformation->getInitialTangent();
					return MatrixWrapper::GetMatrixWrapper(&mat);
				}

				MatrixWrapper^ GetSectionFlexibility() {
					Matrix mat = _SectionForceDeformation->getSectionFlexibility();
					return MatrixWrapper::GetMatrixWrapper(&mat);
				}

				MatrixWrapper^ GetInitialFlexibility() {
					Matrix mat = _SectionForceDeformation->getInitialFlexibility();
					return MatrixWrapper::GetMatrixWrapper(&mat);
				}

				double GetRho() {
					return _SectionForceDeformation->getRho();
				}

				~SectionForceDeformationWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};
			internal:
				SectionForceDeformation * _SectionForceDeformation;
			};

#pragma region ysSection
			public ref class YieldSurfaceSection2dWrapper abstract : SectionForceDeformationWrapper
			{
			public:
				YieldSurfaceSection2dWrapper(int tag, int classtag,
					YieldSurface_BCWrapper^ ptrys,
					bool use_kr) {};
				YieldSurfaceSection2dWrapper() {};
				~YieldSurfaceSection2dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:
			};

			public ref class YS_Section2D01Wrapper : YieldSurfaceSection2dWrapper
			{
			public:
				YS_Section2D01Wrapper(int tag, double E, double A, double I,
					YieldSurface_BCWrapper^ ptrys,
					bool use_kr);
				~YS_Section2D01Wrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:
			};

			public ref class YS_Section2D02Wrapper : YieldSurfaceSection2dWrapper
			{
			public:
				YS_Section2D02Wrapper(int tag,
					double E, double A, double I,
					double theta_p_max,
					YieldSurface_BCWrapper^ ptrys,
					bool use_kr);
				~YS_Section2D02Wrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:
			};

#pragma endregion

			public ref class BidirectionalWrapper : SectionForceDeformationWrapper
			{
			public:
				BidirectionalWrapper(int tag, double E, double sigY, double Hiso, double Hkin,
					int code1,
					int code2);

				BidirectionalWrapper() {};
				~BidirectionalWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:
			};

			public ref class ElasticMembranePlateSectionWrapper : SectionForceDeformationWrapper
			{
			public:
				ElasticMembranePlateSectionWrapper(int    tag,
					double E,
					double nu,
					double h ,
					double rho);

				ElasticMembranePlateSectionWrapper() {};
				~ElasticMembranePlateSectionWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:
			};

			public ref class ElasticPlateSectionWrapper : SectionForceDeformationWrapper
			{
			public:
				ElasticPlateSectionWrapper(int    tag,
					double E,
					double nu,
					double h);

				ElasticPlateSectionWrapper() {};
				~ElasticPlateSectionWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:
			};

			public ref class ElasticSection2dWrapper : SectionForceDeformationWrapper
			{
			public:
				ElasticSection2dWrapper(int tag, double E, double A, double I);

				ElasticSection2dWrapper() {};
				~ElasticSection2dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:
			};

			public ref class ElasticSection3dWrapper : SectionForceDeformationWrapper
			{
			public:
				ElasticSection3dWrapper(int tag, double E, double A, double Iz,
					double Iy, double G, double J);
				~ElasticSection3dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class ElasticShearSection2dWrapper : SectionForceDeformationWrapper
			{
			public:
				ElasticShearSection2dWrapper(int tag, double E, double A, double I,
					double G, double alpha);
				~ElasticShearSection2dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class ElasticShearSection3dWrapper : SectionForceDeformationWrapper
			{
			public:
				ElasticShearSection3dWrapper(int tag, double E, double A, double Iz,
					double Iy, double G, double J, double alphaY, double alphaZ);
				~ElasticShearSection3dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class ElasticTubeSection3dWrapper : SectionForceDeformationWrapper
			{
			public:
				ElasticTubeSection3dWrapper(int tag, double E, double d, double tw, double G);
				~ElasticTubeSection3dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class ElasticWarpingShearSection2dWrapper : SectionForceDeformationWrapper
			{
			public:
				ElasticWarpingShearSection2dWrapper(int tag, double E, double A, double I,
					double G, double alpha, double J, double B, double C);
				~ElasticWarpingShearSection2dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class Elliptical2Wrapper : SectionForceDeformationWrapper
			{
			public:
				Elliptical2Wrapper(int tag, double E1, double E2, double sigY1, double sigY2,
					double Hiso, double Hkin1, double Hkin2,
					int c1 , int c2 );
				~Elliptical2Wrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class FiberSection2dThermalWrapper : SectionForceDeformationWrapper
			{
			public:
				FiberSection2dThermalWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers);
				~FiberSection2dThermalWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class FiberSection3dThermalWrapper : SectionForceDeformationWrapper
			{
			public:
				FiberSection3dThermalWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers);
				~FiberSection3dThermalWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class FiberSectionGJWrapper : SectionForceDeformationWrapper
			{
			public:
				FiberSectionGJWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double GJ);
				~FiberSectionGJWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class FiberSectionGJThermalWrapper : SectionForceDeformationWrapper
			{
			public:
				FiberSectionGJThermalWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double GJ);
				~FiberSectionGJThermalWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class GenericSection1dWrapper : SectionForceDeformationWrapper
			{
			public:
				GenericSection1dWrapper(int tag, UniaxialMaterialWrapper^ m, int code);
				~GenericSection1dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class Isolator2springWrapper : SectionForceDeformationWrapper
			{
			public:
				Isolator2springWrapper(int tag, double tol_in, double k1_in, double Fy_in, double kb_in, double kvo_in,
					double hb_in, double Pe_in, double po_in);
				~Isolator2springWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class LayeredShellFiberSectionWrapper : SectionForceDeformationWrapper
			{
			public:
				LayeredShellFiberSectionWrapper(int tag,
					int iLayers,
					array<double>^ thickness,
					array<NDMaterialWrapper^>^ fibers);
				~LayeredShellFiberSectionWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class LayeredShellFiberSectionThermalWrapper : SectionForceDeformationWrapper
			{
			public:
				LayeredShellFiberSectionThermalWrapper(int tag,
					int iLayers,
					array<double>^ thickness,
					array<NDMaterialWrapper^>^ fibers);
				~LayeredShellFiberSectionThermalWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class MembranePlateFiberSectionWrapper : SectionForceDeformationWrapper
			{
			public:
				MembranePlateFiberSectionWrapper(int tag,
					double thickness,
					NDMaterialWrapper^ fibers);
				~MembranePlateFiberSectionWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class MembranePlateFiberSectionThermalWrapper : SectionForceDeformationWrapper
			{
			public:
				MembranePlateFiberSectionThermalWrapper(int tag,
					double thickness,
					NDMaterialWrapper^ fibers);
				~MembranePlateFiberSectionThermalWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			
			public ref class NDFiberSection2dWrapper : SectionForceDeformationWrapper
			{
			public:
				NDFiberSection2dWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double a);
				~NDFiberSection2dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class NDFiberSection3dWrapper : SectionForceDeformationWrapper
			{
			public:
				NDFiberSection3dWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double a);
				~NDFiberSection3dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class NDFiberSectionWarping2dWrapper : SectionForceDeformationWrapper
			{
			public:
				NDFiberSectionWarping2dWrapper(int tag, int numFibers, array<FiberWrapper^>^ fibers, double a);
				~NDFiberSectionWarping2dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class ParallelSectionWrapper : SectionForceDeformationWrapper
			{
			public:
				ParallelSectionWrapper(int tag, int numSections, array<SectionForceDeformationWrapper^>^ theSections);
				~ParallelSectionWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class SectionAggregatorWrapper : SectionForceDeformationWrapper
			{
			public:
				SectionAggregatorWrapper(int tag, SectionForceDeformationWrapper^ theSection,
					int numAdditions, array<UniaxialMaterialWrapper^>^ theAdditions,
					IDWrapper^ code);

				SectionAggregatorWrapper(int tag, int numAdditions, array<UniaxialMaterialWrapper^>^ theAdditions,
					IDWrapper^ code);

				SectionAggregatorWrapper(int tag, SectionForceDeformationWrapper^ theSection, UniaxialMaterialWrapper^ theAdditions,
					int code);

				~SectionAggregatorWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class FiberSection3dWrapper : SectionForceDeformationWrapper
			{
			public:
				FiberSection3dWrapper(int tag, array<FiberWrapper^>^ fibers, UniaxialMaterialWrapper^ GJ);
				FiberSection3dWrapper(int tag, FiberSectionReprWrapper ^ fiberSectionRepr,
					Dictionary<int, UniaxialMaterialWrapper^>^ materials, array<UniaxialFiber3dWrapper^>^ fibers, UniaxialMaterialWrapper^ GJ);

				/*FiberSection3dWrapper(int tag, array<FiberWrapper^>^ fibers, array<QuadPatchFiber3dWrapper^>^ quadPatchFiber3ds, UniaxialMaterialWrapper^ GJ);
				FiberSection3dWrapper(int tag, FiberSectionReprWrapper ^ fiberSectionRepr,
					Dictionary<int, UniaxialMaterialWrapper^>^ materials, array<QuadPatchFiber3dWrapper^>^ quadPatchFiber3ds, UniaxialMaterialWrapper^ GJ);*/
				~FiberSection3dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;

					/*if (_QuadPatchFiber3ds != 0)
					{
						for (int i = 0; i < numQuadFibers; i++)
							delete _QuadPatchFiber3ds[i];
						delete[] _QuadPatchFiber3ds;
					}*/

					if (_GJ != 0)
						delete _GJ;
				};

			internal:
				//QuadPatchFiber3d **_QuadPatchFiber3ds;
				int numQuadFibers;
				UniaxialMaterial *_GJ;
			};

			public ref class FiberSection2dWrapper : SectionForceDeformationWrapper
			{
			public:
				FiberSection2dWrapper(int tag, array<FiberWrapper^>^ fibers);
				FiberSection2dWrapper(int tag, FiberSectionReprWrapper ^ fiberSectionRepr,
					Dictionary<int, UniaxialMaterialWrapper^>^ materials, array<UniaxialFiber2dWrapper^>^ fibers);
				~FiberSection2dWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			public ref class FiberSection2dIntWrapper : SectionForceDeformationWrapper
			{
			public:
				FiberSection2dIntWrapper(int tag,
					array<FiberWrapper^>^ fibers,
					array<FiberWrapper^>^ Hfibers,
					int NStrip1,
					double tavg1,
					int NStrip2,
					double tavg2,
					int NStrip3,
					double tavg3);
				~FiberSection2dIntWrapper() {
					if (_SectionForceDeformation != 0)
						delete _SectionForceDeformation;
				};

			internal:

			};

			
		}
	}
}
