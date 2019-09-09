#pragma once
#include <PlasticHardeningMaterial.h>
#include <ExponReducing.h>
#include <MultiLinearKp.h>
#include <NullPlasticMaterial.h>

#include <YS_Evolution.h>
#include <BkStressLimSurface2D.h>
#include <BoundingSurface2D.h>
#include <CombinedIsoKin2D01.h>
#include <CombinedIsoKin2D02.h>
#include <Kinematic2D01.h>
#include <Kinematic2D02.h>
#include <NullEvolution.h>
#include <PeakOriented2D01.h>
#include <PeakOriented2D02.h>
#include <PlasticHardening2D.h>
#include <YS_Evolution2D.h>
#include <Isotropic2D01.h>

#include <Attalla2D.h>
#include <ElTawil2D.h>
#include <ElTawil2DUnSym.h>
#include <Hajjar2D.h>
#include <NullYS2D.h>
#include <Orbison2D.h>
#include <YieldSurface_BC.h>
#include <YieldSurface_BC2D.h>


#include "../MaterialWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../../matrix/VectorWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../../matrix/IDWrapper.h"

#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"

namespace OpenSees {
	namespace Materials {
		namespace YieldSurfaces {

			public ref class YS_EvolutionWrapper abstract : TaggedObjectWrapper, IMovableObjectWrapper
			{
			public:
				YS_EvolutionWrapper() {};
				~YS_EvolutionWrapper() {
				};
			internal:
				YS_Evolution * _YS_Evolution;
			};

			public ref class YS_Evolution2DWrapper abstract : YS_EvolutionWrapper
			{
			public:
				YS_Evolution2DWrapper() {};
				~YS_Evolution2DWrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:
				YS_Evolution2D * _YS_Evolution2D;
			};

			public ref class PlasticHardeningMaterialWrapper abstract : MaterialWrapper
			{
			public:
				PlasticHardeningMaterialWrapper() {};
				~PlasticHardeningMaterialWrapper() {
					if (_PlasticHardeningMaterial != 0)
						delete _PlasticHardeningMaterial;
				};

			internal:
				PlasticHardeningMaterial * _PlasticHardeningMaterial;
			};

			public ref class YieldSurface_BCWrapper abstract : TaggedObjectWrapper, IMovableObjectWrapper
			{
			public:
				YieldSurface_BCWrapper() {};
				~YieldSurface_BCWrapper() {
					if (_YieldSurface_BC != 0)
						delete _YieldSurface_BC;
				};
			internal:
				YieldSurface_BC * _YieldSurface_BC;
			};

			public ref class YieldSurface_BC2DWrapper abstract : YieldSurface_BCWrapper
			{
			public:
				YieldSurface_BC2DWrapper() {};
				YieldSurface_BC2DWrapper(int tag, int classTag, double xmax, double ymax,
					YS_EvolutionWrapper^ model) {};
				~YieldSurface_BC2DWrapper() {
					if (_YieldSurface_BC2D != 0)
						delete _YieldSurface_BC2D;
				};
			internal:
				YieldSurface_BC2D * _YieldSurface_BC2D;
			};

#pragma region plasticHardening
			

			public ref class ExponReducingWrapper : PlasticHardeningMaterialWrapper
			{
			public:
				ExponReducingWrapper(int tag, double kp0, double alfa);
				ExponReducingWrapper(int tag, double kp0, double alfa, double res_fact);
				~ExponReducingWrapper() {
					if (_PlasticHardeningMaterial != 0)
						delete _PlasticHardeningMaterial;
				};

			internal:

			};

			public ref class MultiLinearKpWrapper : PlasticHardeningMaterialWrapper
			{
			public:
				MultiLinearKpWrapper(int tag, VectorWrapper^ sum_plas_defo, VectorWrapper^ kp);
				~MultiLinearKpWrapper() {
					if (_PlasticHardeningMaterial != 0)
						delete _PlasticHardeningMaterial;
				};

			internal:

			};

			public ref class NullPlasticMaterialWrapper : PlasticHardeningMaterialWrapper
			{
			public:
				NullPlasticMaterialWrapper(int tag);
				~NullPlasticMaterialWrapper() {
					if (_PlasticHardeningMaterial != 0)
						delete _PlasticHardeningMaterial;
				};

			internal:

			};

#pragma endregion

#pragma region yieldSurfaceBC
			

			public ref class Attalla2DWrapper : YieldSurface_BC2DWrapper
			{
			public:
				Attalla2DWrapper(int tag, double xmax, double ymax,
					YS_EvolutionWrapper^ model);
				Attalla2DWrapper(int tag, double xmax, double ymax,
					YS_EvolutionWrapper^ model, double a01, double a02, double a03,
					double a04, double a05, double a06);
				~Attalla2DWrapper() {
					if (_YieldSurface_BC != 0)
						delete _YieldSurface_BC;
				};
			
			};

			public ref class ElTawil2DWrapper : YieldSurface_BC2DWrapper
			{
			public:
				ElTawil2DWrapper(int tag, double xbal, double ybal, double ypos, double yneg,
					YS_EvolutionWrapper^ model);
				ElTawil2DWrapper(int tag, double xbal, double ybal, double ypos, double yneg,
					YS_EvolutionWrapper^ model, double cz, double ty);
				~ElTawil2DWrapper() {
					if (_YieldSurface_BC != 0)
						delete _YieldSurface_BC;
				};
			
			};

			public ref class ElTawil2DUnSymWrapper : YieldSurface_BC2DWrapper
			{
			public:
				ElTawil2DUnSymWrapper(int tag, double xPosBal, double yPosBal,
					double xNegBal, double yNegBal,
					double ypos, double yneg,
					YS_EvolutionWrapper^ model);
				ElTawil2DUnSymWrapper(int tag, double xPosBal, double yPosBal,
					double xNegBal, double yNegBal,
					double ypos, double yneg,
					YS_EvolutionWrapper^ model, double czPos, double tyPos,
					double czNeg, double tyNeg);
				~ElTawil2DUnSymWrapper() {
					if (_YieldSurface_BC != 0)
						delete _YieldSurface_BC;
				};
			internal:
				YieldSurface_BC * _YieldSurface_BC;
			};


			public ref class Hajjar2DWrapper : YieldSurface_BC2DWrapper
			{
			public:
				Hajjar2DWrapper(int tag, double xmax, double ymax,
					YS_EvolutionWrapper^ model,
					double centroid_y, double c1, double c2, double c3);
				Hajjar2DWrapper(int tag, YS_EvolutionWrapper^ model,
					double D, double b, double t, double fc_, double fy_);
				~Hajjar2DWrapper() {
					if (_YieldSurface_BC != 0)
						delete _YieldSurface_BC;
				};
			internal:
				YieldSurface_BC * _YieldSurface_BC;
			};


			public ref class NullYS2DWrapper : YieldSurface_BC2DWrapper
			{
			public:
				NullYS2DWrapper(int tag);
				~NullYS2DWrapper() {
					if (_YieldSurface_BC != 0)
						delete _YieldSurface_BC;
				};
			internal:
				YieldSurface_BC * _YieldSurface_BC;
			};

			public ref class Orbison2DWrapper : YieldSurface_BC2DWrapper
			{
			public:
				Orbison2DWrapper(int tag, double xmax, double ymax, YS_EvolutionWrapper^ model);
				~Orbison2DWrapper() {
					if (_YieldSurface_BC != 0)
						delete _YieldSurface_BC;
				};
			internal:
				YieldSurface_BC * _YieldSurface_BC;
			};

#pragma endregion

#pragma region evolution

			public ref class BkStressLimSurface2DWrapper abstract : YS_Evolution2DWrapper
			{
			public:
				BkStressLimSurface2DWrapper() {};
				BkStressLimSurface2DWrapper(int tag, int classTag, double min_iso_factor,
					double iso_ratio, double kin_ratio,
					YieldSurface_BCWrapper^  lim_surface,
					PlasticHardeningMaterialWrapper^ kinX,
					PlasticHardeningMaterialWrapper^  kinY,
					PlasticHardeningMaterialWrapper^  isoXPos,
					PlasticHardeningMaterialWrapper^  isoXNeg,
					PlasticHardeningMaterialWrapper^  isoYPos,
					PlasticHardeningMaterialWrapper^  isoYNeg,
					int restype, double res_Fact, double app_Fact, double dir) {};
				~BkStressLimSurface2DWrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class BoundingSurface2DWrapper abstract : YS_Evolution2DWrapper
			{
			public:
				BoundingSurface2DWrapper(int tag, int classTag, double min_iso_factor,
					double iso_ratio, double kin_ratio,
					PlasticHardeningMaterialWrapper^ kpx,
					PlasticHardeningMaterialWrapper^  kpy,
					YieldSurface_BCWrapper^  lim_surface) {};
				~BoundingSurface2DWrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class PlasticHardening2DWrapper abstract : YS_Evolution2DWrapper
			{
			public:
				PlasticHardening2DWrapper() {};
				PlasticHardening2DWrapper(int tag, int classTag, double min_iso_factor,
					double iso_ratio, double kin_ratio,
					PlasticHardeningMaterialWrapper^ kpx_pos,
					PlasticHardeningMaterialWrapper^  kpx_neg,
					PlasticHardeningMaterialWrapper^ kpy_pos,
					PlasticHardeningMaterialWrapper^  kpy_neg,
					double dir) {};
				~PlasticHardening2DWrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class CombinedIsoKin2D01Wrapper :PlasticHardening2DWrapper
			{
			public:
				CombinedIsoKin2D01Wrapper(int tag,
					double iso_ratio, double kin_ratio,
					double shr_iso_ratio, double shr_kin_ratio,
					double min_iso_factor,
					PlasticHardeningMaterialWrapper^ kpx_pos,
					PlasticHardeningMaterialWrapper^  kpx_neg,
					PlasticHardeningMaterialWrapper^ kpy_pos,
					PlasticHardeningMaterialWrapper^  kpy_neg,
					bool isDeformable, double dir);
				~CombinedIsoKin2D01Wrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class CombinedIsoKin2D02Wrapper :BkStressLimSurface2DWrapper
			{
			public:
				CombinedIsoKin2D02Wrapper(int tag, double min_iso_factor,
					double iso_ratio, double kin_ratio,
					YieldSurface_BCWrapper^  lim_surface,
					PlasticHardeningMaterialWrapper^ kinX,
					PlasticHardeningMaterialWrapper^  kinY,
					PlasticHardeningMaterialWrapper^ isoXPos,
					PlasticHardeningMaterialWrapper^  isoXNeg,
					PlasticHardeningMaterialWrapper^ isoYPos,
					PlasticHardeningMaterialWrapper^  isoYNeg,
					bool isDeformable,
					int  algo, double resfact, double appfact, double dir);
				~CombinedIsoKin2D02Wrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class Isotropic2D01Wrapper :PlasticHardening2DWrapper
			{
			public:
				Isotropic2D01Wrapper(int tag, double min_iso_factor,
					PlasticHardeningMaterialWrapper^ kpx,
					PlasticHardeningMaterialWrapper^  kpy);
				~Isotropic2D01Wrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class Kinematic2D01Wrapper :PlasticHardening2DWrapper
			{
			public:
				Kinematic2D01Wrapper(int tag, double min_iso_factor,
					PlasticHardeningMaterialWrapper^ kpx,
					PlasticHardeningMaterialWrapper^  kpy, double dir);
				~Kinematic2D01Wrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class Kinematic2D02Wrapper :BkStressLimSurface2DWrapper
			{
			public:
				Kinematic2D02Wrapper(int tag, double min_iso_factor,
					YieldSurface_BCWrapper^ lim_surface,
					PlasticHardeningMaterialWrapper^ kpx,
					PlasticHardeningMaterialWrapper^  kpy,
					int algo, double resfact, double appfact, double dir);
				~Kinematic2D02Wrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class NullEvolutionWrapper :YS_EvolutionWrapper
			{
			public:
				NullEvolutionWrapper(int tag, double isox);
				NullEvolutionWrapper(int tag, double isox, double isoy);
				NullEvolutionWrapper(int tag, double isox, double isoy, double isoz);
				~NullEvolutionWrapper() {
					if (_YS_Evolution != 0)
						delete _YS_Evolution;
				};
			internal:

			};

			public ref class PeakOriented2D01Wrapper :PlasticHardening2DWrapper
			{
			public:
				PeakOriented2D01Wrapper(int tag, double min_iso_factor,
					PlasticHardeningMaterialWrapper^ kpx,
					PlasticHardeningMaterialWrapper^  kpy);
				
				~PeakOriented2D01Wrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			public ref class PeakOriented2D02Wrapper :BkStressLimSurface2DWrapper
			{
			public:
				PeakOriented2D02Wrapper(int tag, double min_iso_factor,
					YieldSurface_BCWrapper^  lim_surface,
					PlasticHardeningMaterialWrapper^ kinX,
					PlasticHardeningMaterialWrapper^  kinY,
					PlasticHardeningMaterialWrapper^ isoX,
					PlasticHardeningMaterialWrapper^  isoY, int algo);

				~PeakOriented2D02Wrapper() {
					if (_YS_Evolution2D != 0)
						delete _YS_Evolution2D;
				};
			internal:

			};

			


#pragma endregion



		}
	}
}
