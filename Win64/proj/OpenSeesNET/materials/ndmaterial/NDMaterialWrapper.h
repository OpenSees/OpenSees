#pragma once
#include <NDMaterial.h>

#include <ExternalNDMaterial.h>

#include <cyclicSoil/MultiaxialCyclicPlasticity.h>
#include <cyclicSoil/MultiaxialCyclicPlasticity3D.h>
#include <cyclicSoil/MultiaxialCyclicPlasticityAxiSymm.h>
#include <cyclicSoil/MultiaxialCyclicPlasticityPlaneStrain.h>

#include <ElasticIsotropic3DThermal.h>
#include <ElasticIsotropicAxiSymm.h>
#include <ElasticIsotropicBeamFiber.h>
#include <ElasticIsotropicBeamFiber2d.h>
#include <ElasticIsotropicMaterial.h>
#include <ElasticIsotropicMaterialThermal.h>
#include <ElasticIsotropicPlaneStrain2D.h>
#include <ElasticIsotropicPlaneStress2D.h>
#include <ElasticIsotropicPlateFiber.h>
#include <ElasticIsotropicThreeDimensional.h>
#include <PressureDependentElastic3D.h>

#include <FeapMaterial.h>
#include <feap/FeapMaterial01.h>
#include <feap/FeapMaterial02.h>
#include <feap/FeapMaterial03.h>

#include <J2AxiSymm.h>
#include <J2BeamFiber2d.h>
#include <J2BeamFiber3d.h>
#include <uwmaterials/J2CyclicBoundingSurface.h>
#include <J2PlaneStrain.h>
#include <J2PlaneStress.h>
#include <J2Plasticity.h>
#include <J2PlasticityThermal.h>
#include <J2PlateFiber.h>
#include <J2PlateFibre.h>
#include <J2ThreeDimensional.h>
#include <J2ThreeDimensionalThermal.h>


#include <matcmm/MaterialCMM.h>



#include <reinforcedConcretePlaneStress\FAFourSteelPCPlaneStress.h>
#include <reinforcedConcretePlaneStress\FAFourSteelRCPlaneStress.h>
#include <reinforcedConcretePlaneStress\FAPrestressedConcretePlaneStress.h>
#include <reinforcedConcretePlaneStress\FAReinforcedConcretePlaneStress.h>
#include <reinforcedConcretePlaneStress\PrestressedConcretePlaneStress.h>
#include <reinforcedConcretePlaneStress\RAFourSteelPCPlaneStress.h>
#include <reinforcedConcretePlaneStress\FAPrestressedConcretePlaneStress.h>
#include <reinforcedConcretePlaneStress\RAFourSteelRCPlaneStress.h>
#include <reinforcedConcretePlaneStress\ReinforcedConcretePlaneStress.h>



#include <soil\FluidSolidPorousMaterial.h>
#include <soil\MultiYieldSurface.h>
#include <soil\MultiYieldSurfaceClay.h>
#include <soil\PressureDependMultiYield.h>
#include <soil\PressureDependMultiYield02.h>
#include <soil\PressureIndependMultiYield.h>


//#include <stressdensitymodel\StressDensityModel.h>
//#include <stressdensitymodel\StressDensityModel2D.h>
//#include <stressdensitymodel\StressDensityModel3D.h>


#include <BeamFiberMaterial.h>
#include <BeamFiberMaterial2d.h>
#include <PlaneStrainMaterial.h>
#include <PlaneStressLayeredMaterial.h>
#include <PlaneStressMaterial.h>
#include <PlaneStressRebarMaterial.h>
#include <PlaneStressSimplifiedJ2.h>
#include <PlaneStressUserMaterial.h>
#include <PlateFiberMaterial.h>
#include <PlateFiberMaterialThermal.h>
#include <PlateFromPlaneStressMaterial.h>
#include <PlateFromPlaneStressMaterialThermal.h>
#include <PlateRebarMaterial.h>
#include <PlateRebarMaterialThermal.h>
#include <AcousticMedium.h>

#include <uwmaterials/BoundingCamClay.h>
#include <uwmaterials/BoundingCamClay3D.h>
#include <uwmaterials/BoundingCamClayPlaneStrain.h>

#include <CapPlasticity.h>
#include <ConcreteS.h>
#include <uwmaterials/ContactMaterial2D.h>
#include <uwmaterials/ContactMaterial3D.h>
#include <CycLiqCP.h>
#include <CycLiqCP3D.h>
#include <CycLiqCPPlaneStrain.h>
#include <CycLiqCPSP.h>
#include <CycLiqCPSP3D.h>
#include <CycLiqCPSPPlaneStrain.h>
#include <uwmaterials/DruckerPrager.h>
#include <uwmaterials/DruckerPrager3D.h>
#include <DruckerPrager3DThermal.h>
#include <uwmaterials/DruckerPragerPlaneStrain.h>
#include <DruckerPragerThermal.h>
#include <ElasticOrthotropicMaterial.h>
#include <ElasticOrthotropicThreeDimensional.h>
#include <FSAM.h>
#include <uwmaterials/InitialStateAnalysisWrapper.h>
#include <InitStressNDMaterial.h>
#include <LinearCap.h>
#include <uwmaterials/ManzariDafalias.h>
#include <uwmaterials/ManzariDafalias3D.h>
#include <uwmaterials/ManzariDafalias3DRO.h>
#include <uwmaterials/ManzariDafaliasPlaneStrain.h>
#include <uwmaterials/ManzariDafaliasPlaneStrainRO.h>
#include <uwmaterials/ManzariDafaliasRO.h>
#include <PlasticDamageConcrete3d.h>
#include <PlasticDamageConcretePlaneStress.h>
#include <uwmaterials/PM4Sand.h>
#include <uwmaterials/PM4Silt.h>
#include <SimplifiedJ2.h>
#include <ElasticIsotropicMaterial.h>

#include <classTags.h>

#include "../MaterialWrapper.h"
#include "../uniaxial/UniaxialMaterialWrapper.h"
#include "../uniaxial/UniaxialMaterialWrapper_all.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../../matrix/VectorWrapper.h"
#include "../../handlers/HandlerWrapper.h"
#include "../../recorders/ResponseWrapper.h"
#include "../../recorders/InformationWrapper.h"


using namespace OpenSees;
using namespace OpenSees::Recorders;
using namespace OpenSees::Handlers;

using namespace System;
using namespace System::Runtime::InteropServices;
using namespace OpenSees::Materials::Uniaxials;

namespace OpenSees {
	namespace Materials {
		namespace NDMaterials {
			public ref class NDMaterialWrapper : MaterialWrapper, IMovableObjectWrapper
			{
			public:
				NDMaterialWrapper();
				
				array<double, 2>^ GetTangent() {
					return MatrixWrapper::GetArray(_NDMaterial->getTangent());
				}

				MatrixWrapper^ GetTangentMatrixWrapper() {
					return MatrixWrapper::GetMatrixWrapper(_NDMaterial->getTangent());
				}

				array<double, 2>^ GetInitialTangent() {
					return MatrixWrapper::GetArray(_NDMaterial->getInitialTangent());
				}

				MatrixWrapper^ GetInitialTangentMatrixWrapper() {
					return MatrixWrapper::GetMatrixWrapper(_NDMaterial->getInitialTangent());
				}

				array<double>^ GetStress() {
					return VectorWrapper::GetArray(_NDMaterial->getStress());
				}

				VectorWrapper^ GetStressVectorWrapper() {
					return VectorWrapper::GetVectorWrapper(_NDMaterial->getStress());
				}

				array<double>^ GetStrain() {
					return VectorWrapper::GetArray(_NDMaterial->getStrain());
				}

				VectorWrapper^ GetStrainVectorWrapper() {
					return VectorWrapper::GetVectorWrapper(_NDMaterial->getStrain());
				}

				String^ GetType() {
					return gcnew String(_NDMaterial->getClassType());
				}

				int GetTypeTag() {
					return _NDMaterial->getClassTag();
				}

				int GetOrder() {
					return _NDMaterial->getOrder();
				}

				int CommitState() {
					return _NDMaterial->commitState();
				}
				int RevertToLastCommit() {
					return _NDMaterial->revertToLastCommit();
				}
				int RevertToStart() {
					return _NDMaterial->revertToStart();
				}

				ResponseWrapper^ SetResponse(array<String^>^ argv, OpenSees::Handlers::OPS_StreamWrapper^ output) {
					const char ** _argv = new const char*[argv->Length];
					for (int i = 0; i < argv->Length; i++)
					{
						_argv[i] = OPS::StringToChar(argv[i]);
					}

					Response* response = this->_NDMaterial->setResponse(_argv, argv->Length, *output->_OPS_StreamPtr);
					ResponseWrapper^ _response = gcnew ResponseWrapper();
					_response->_Response = response;
					return _response;
				}

				int GetResponse(int responseID, InformationWrapper^ info) {
					return _NDMaterial->getResponse(responseID, *info->_Information);
				}

				~NDMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:
				NDMaterial * _NDMaterial;
			};

#pragma region cyclicsoil
			public ref class MultiaxialCyclicPlasticityWrapper abstract : NDMaterialWrapper
			{
			public:
				MultiaxialCyclicPlasticityWrapper(int    tag,
					double rho,
					double K,
					double G,
					double Su,
					double Ho_kin,
					double Parameterh,
					double Parameter_m,
					double Parameter_beta,
					double Kcoeff,
					double viscosity);

				MultiaxialCyclicPlasticityWrapper(int tag, double rho, double K, double G);

				MultiaxialCyclicPlasticityWrapper() {};
				~MultiaxialCyclicPlasticityWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class MultiaxialCyclicPlasticity3DWrapper abstract : MultiaxialCyclicPlasticityWrapper
			{
			public:
				MultiaxialCyclicPlasticity3DWrapper(int    tag,
					double rho,
					double K,
					double G,
					double Su,
					double Ho_kin,
					double Parameterh,
					double Parameter_m,
					double Parameter_beta,
					double Kcoeff,
					double viscosity);

				MultiaxialCyclicPlasticity3DWrapper(int tag, double rho, double K, double G);

				MultiaxialCyclicPlasticity3DWrapper() {};
				~MultiaxialCyclicPlasticity3DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class MultiaxialCyclicPlasticityAxiSymmWrapper abstract : MultiaxialCyclicPlasticityWrapper
			{
			public:
				MultiaxialCyclicPlasticityAxiSymmWrapper(int    tag,
					double rho,
					double K,
					double G,
					double Su,
					double Ho_kin,
					double Parameterh,
					double Parameter_m,
					double Parameter_beta,
					double Kcoeff,
					double viscosity);

				MultiaxialCyclicPlasticityAxiSymmWrapper(int tag, double rho, double K, double G);

				MultiaxialCyclicPlasticityAxiSymmWrapper() {};
				~MultiaxialCyclicPlasticityAxiSymmWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class MultiaxialCyclicPlasticityPlaneStrainWrapper abstract : MultiaxialCyclicPlasticityWrapper
			{
			public:
				MultiaxialCyclicPlasticityPlaneStrainWrapper(int    tag,
					double rho,
					double K,
					double G,
					double Su,
					double Ho_kin,
					double Parameterh,
					double Parameter_m,
					double Parameter_beta,
					double Kcoeff,
					double viscosity);

				MultiaxialCyclicPlasticityPlaneStrainWrapper(int tag, double rho, double K, double G);

				MultiaxialCyclicPlasticityPlaneStrainWrapper() {};
				~MultiaxialCyclicPlasticityPlaneStrainWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

#pragma endregion

#pragma region elasticIsotropic
			/*public ref class ElasticIsotropicMaterialWrapper abstract : NDMaterialWrapper
			{
			public:
				ElasticIsotropicMaterialWrapper() {};
				ElasticIsotropicMaterialWrapper(int tag, int classTag, double E, double nu, double rho);
				ElasticIsotropicMaterialWrapper(int tag, double E, double nu, double rho);

				~ElasticIsotropicMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};*/

			public ref class ElasticIsotropicMaterialThermalWrapper abstract : NDMaterialWrapper
			{
			public:
				ElasticIsotropicMaterialThermalWrapper(int tag, double E, double nu, double rho, double alpha, int softindex);
				ElasticIsotropicMaterialThermalWrapper() {};
				~ElasticIsotropicMaterialThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class ElasticIsotropic3DThermalWrapper abstract : ElasticIsotropicMaterialThermalWrapper
			{
			public:
				ElasticIsotropic3DThermalWrapper(int tag, double e, double nu, double rho, double alpha, int softindex);
				ElasticIsotropic3DThermalWrapper() {};
				~ElasticIsotropic3DThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class ElasticIsotropicAxiSymmWrapper : NDMaterialWrapper
			{
			public:
				ElasticIsotropicAxiSymmWrapper(int tag, double E, double nu, double rho);
				~ElasticIsotropicAxiSymmWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class ElasticIsotropicBeamFiberWrapper : NDMaterialWrapper
			{
			public:
				ElasticIsotropicBeamFiberWrapper(int tag, double E, double nu, double rho);
				~ElasticIsotropicBeamFiberWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class ElasticIsotropicBeamFiber2dWrapper : NDMaterialWrapper
			{
			public:
				ElasticIsotropicBeamFiber2dWrapper(int tag, double E, double nu, double rho);
				~ElasticIsotropicBeamFiber2dWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class ElasticIsotropicPlaneStress2DWrapper :NDMaterialWrapper
			{
			public:
				ElasticIsotropicPlaneStress2DWrapper(int tag, double E, double nu, double rho);
				~ElasticIsotropicPlaneStress2DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class ElasticIsotropicPlaneStrain2DWrapper : NDMaterialWrapper
			{
			public:
				ElasticIsotropicPlaneStrain2DWrapper(int tag, double E, double nu, double rho);
				~ElasticIsotropicPlaneStrain2DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class ElasticIsotropicPlateFiberWrapper :NDMaterialWrapper
			{
			public:
				ElasticIsotropicPlateFiberWrapper(int tag, double E, double nu, double rho);
				~ElasticIsotropicPlateFiberWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class ElasticIsotropicThreeDimensionalWrapper :NDMaterialWrapper
			{
			public:
				ElasticIsotropicThreeDimensionalWrapper(int tag, double E, double nu, double rho);
				~ElasticIsotropicThreeDimensionalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class PressureDependentElastic3DWrapper :NDMaterialWrapper
			{
			public:
				PressureDependentElastic3DWrapper(int tag,
					double E,
					double nu,
					double rhop,
					double expp,
					double pr,
					double pop);
				~PressureDependentElastic3DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

#pragma endregion

#pragma region Feap
			public ref class FeapMaterialWrapper :NDMaterialWrapper
			{
			public:
				FeapMaterialWrapper(int tag, int classTag, int numHV, int numData,
					double rho);
				~FeapMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class FeapMaterial01Wrapper :NDMaterialWrapper
			{
			public:
				FeapMaterial01Wrapper(int tag, double E, double nu, double rho);
				~FeapMaterial01Wrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class FeapMaterial02Wrapper :NDMaterialWrapper
			{
			public:
				FeapMaterial02Wrapper(int tag, double K, double G, double muK, double muG,
					double lamK, double lamG, double theta);
				~FeapMaterial02Wrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class FeapMaterial03Wrapper :NDMaterialWrapper
			{
			public:
				FeapMaterial03Wrapper(int tag, double K, double G, double sigY, double Hiso);
				~FeapMaterial03Wrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

#pragma endregion

#pragma region j2Plasticity
			public ref class J2PlasticityWrapper :NDMaterialWrapper
			{
			public:
				J2PlasticityWrapper(int    tag,
					int    classTag,
					double K,
					double G,
					double yield0,
					double yield_infty,
					double d,
					double H,
					double viscosity,
					double rho);
				J2PlasticityWrapper(int tag, int classTag, double K, double G);
				J2PlasticityWrapper() {};
				~J2PlasticityWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2AxiSymmWrapper :J2PlasticityWrapper
			{
			public:
				J2AxiSymmWrapper(int    tag,
					double K,
					double G,
					double yield0,
					double yield_infty,
					double d,
					double H,
					double viscosity,
					double rho);
				J2AxiSymmWrapper(int tag, double K, double G);
				J2AxiSymmWrapper() {};
				~J2AxiSymmWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2BeamFiber2dWrapper :NDMaterialWrapper
			{
			public:
				J2BeamFiber2dWrapper(int tag, double E, double G, double sigY, double Hi, double Hk);
				J2BeamFiber2dWrapper() {};
				~J2BeamFiber2dWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2BeamFiber3dWrapper :NDMaterialWrapper
			{
			public:
				J2BeamFiber3dWrapper(int tag, double E, double G, double sigY, double Hi, double Hk);
				J2BeamFiber3dWrapper() {};
				~J2BeamFiber3dWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2CyclicBoundingSurfaceWrapper :NDMaterialWrapper
			{
			public:
				J2CyclicBoundingSurfaceWrapper(int    tag,
					double G,
					double K,
					double su,
					double rho,
					double h,
					double m,
					double k_in,
					double beta);
				J2CyclicBoundingSurfaceWrapper() {};
				~J2CyclicBoundingSurfaceWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2PlaneStrainWrapper :J2PlasticityWrapper
			{
			public:
				J2PlaneStrainWrapper(int    tag,
					double K,
					double G,
					double yield0,
					double yield_infty,
					double d,
					double H,
					double viscosity);
				J2PlaneStrainWrapper() {};
				~J2PlaneStrainWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2PlaneStressWrapper :J2PlasticityWrapper
			{
			public:
				J2PlaneStressWrapper(int    tag,
					double K,
					double G,
					double yield0,
					double yield_infty,
					double d,
					double H,
					double viscosity,
					double rho);
				J2PlaneStressWrapper(int tag, double K, double G);
				J2PlaneStressWrapper() {};
				~J2PlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2PlasticityThermalWrapper :NDMaterialWrapper
			{
			public:
				J2PlasticityThermalWrapper(int    tag,
					int    classTag,
					double K,
					double G,
					double yield0,
					double yield_infty,
					double d,
					double H,
					double viscosity,
					double rho);
				J2PlasticityThermalWrapper(int tag, int classTag, double K, double G);
				J2PlasticityThermalWrapper() {};
				~J2PlasticityThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2PlateFiberWrapper :J2PlasticityWrapper
			{
			public:
				J2PlateFiberWrapper(int    tag,
					double K,
					double G,
					double yield0,
					double yield_infty,
					double d,
					double H,
					double viscosity ,
					double rho);
				J2PlateFiberWrapper(int tag, double K, double G);
				J2PlateFiberWrapper() {};
				~J2PlateFiberWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2PlateFibreWrapper :NDMaterialWrapper
			{
			public:
				J2PlateFibreWrapper(int tag, double E, double G, double sigY, double Hi, double Hk);
				J2PlateFibreWrapper() {};
				~J2PlateFibreWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2ThreeDimensionalWrapper :J2PlasticityWrapper
			{
			public:
				J2ThreeDimensionalWrapper(int    tag,
					double K,
					double G,
					double yield0,
					double yield_infty,
					double d,
					double H,
					double viscosity,
					double rho);
				J2ThreeDimensionalWrapper() {};
				~J2ThreeDimensionalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

			public ref class J2ThreeDimensionalThermalWrapper :J2PlasticityThermalWrapper
			{
			public:
				J2ThreeDimensionalThermalWrapper(int    tag,
					double K,
					double G,
					double yield0,
					double yield_infty,
					double d,
					double H,
					double viscosity,
					double rho);
				J2ThreeDimensionalThermalWrapper(int tag, double K, double G);
				J2ThreeDimensionalThermalWrapper() {};
				~J2ThreeDimensionalThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};

			internal:

			};

#pragma endregion

#pragma region matCMM
			public ref class MaterialCMMWrapper : NDMaterialWrapper
			{
			public:
				MaterialCMMWrapper(int tag, int layer, array<double>^ props);
				~MaterialCMMWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};
#pragma endregion

#pragma region reinforcedConcretePlaneStress
			

			public ref class FAFourSteelPCPlaneStressWrapper : NDMaterialWrapper
			{
			public:
				FAFourSteelPCPlaneStressWrapper(int      tag,
					double   RHO,
					UniaxialMaterialWrapper^ t1,
					UniaxialMaterialWrapper^ t2,
					UniaxialMaterialWrapper^ s1,
					UniaxialMaterialWrapper^ s2,
					UniaxialMaterialWrapper^ c1,
					UniaxialMaterialWrapper^ c2,
					double   ANGLE1,
					double   ANGLE2,
					double   ANGLE3,
					double   ANGLE4,
					double   ROU1,
					double   ROU2,
					double   ROU3,
					double   ROU4,
					double	PSTRAIN1,
					double	PSTRAIN2,
					double   FPC,
					double   FY1,
					double	FY2,
					double   E,
					double   EPSC0);
				~FAFourSteelPCPlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};

			public ref class FAFourSteelRCPlaneStressWrapper : NDMaterialWrapper
			{
			public:
				FAFourSteelRCPlaneStressWrapper(int      tag,
					double   RHO,
					UniaxialMaterialWrapper^ t1,
					UniaxialMaterialWrapper^ t2,
					UniaxialMaterialWrapper^ s1,
					UniaxialMaterialWrapper^ s2,
					UniaxialMaterialWrapper^ c1,
					UniaxialMaterialWrapper^ c2,
					double   ANGLE1,
					double   ANGLE2,
					double   ANGLE3,
					double   ANGLE4,
					double   ROU1,
					double   ROU2,
					double   ROU3,
					double   ROU4,
					double   FPC,
					double   FY,
					double   E,
					double   EPSC0);
				~FAFourSteelRCPlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};

			public ref class FAPrestressedConcretePlaneStressWrapper : NDMaterialWrapper
			{
			public:
				FAPrestressedConcretePlaneStressWrapper(int      tag,
					double   RHO,
					UniaxialMaterialWrapper^ t1,
					UniaxialMaterialWrapper^ s1,
					UniaxialMaterialWrapper^ c1,
					UniaxialMaterialWrapper^ c2,
					double   ANGLE1,
					double   ANGLE2,
					double   ROU1,
					double   ROU2,
					double	PSTRAIN,
					double   FPC,
					double   FY1,
					double	FY2,
					double   E,
					double   EPSC0);
				~FAPrestressedConcretePlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};

			public ref class FAReinforcedConcretePlaneStressWrapper : NDMaterialWrapper
			{
			public:
				FAReinforcedConcretePlaneStressWrapper(int      tag,
					double   RHO,
					UniaxialMaterialWrapper^ t1,
					UniaxialMaterialWrapper^ s1,
					UniaxialMaterialWrapper^ c1,
					UniaxialMaterialWrapper^ c2,
					double   ANGLE1,
					double   ANGLE2,
					double   ROU1,
					double   ROU2,
					double   FPC,
					double   FY,
					double   E,
					double   EPSC0);
				~FAReinforcedConcretePlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};

			public ref class PrestressedConcretePlaneStressWrapper : NDMaterialWrapper
			{
			public:
				PrestressedConcretePlaneStressWrapper(int      tag,
					double   RHO,
					UniaxialMaterialWrapper^ t1,
					UniaxialMaterialWrapper^ s1,
					UniaxialMaterialWrapper^ c1,
					UniaxialMaterialWrapper^ c2,
					double   ANGLE1,
					double   ANGLE2,
					double   ROU1,
					double   ROU2,
					double	PSTRAIN,
					double   FPC,
					double   FY1,
					double	FY2,
					double   E,
					double   EPSC0);
				~PrestressedConcretePlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};

			public ref class RAFourSteelPCPlaneStressWrapper : NDMaterialWrapper
			{
			public:
				RAFourSteelPCPlaneStressWrapper(int      tag,
					double   RHO,
					UniaxialMaterialWrapper^ t1,
					UniaxialMaterialWrapper^ t2,
					UniaxialMaterialWrapper^ s1,
					UniaxialMaterialWrapper^ s2,
					UniaxialMaterialWrapper^ c1,
					UniaxialMaterialWrapper^ c2,
					double   ANGLE1,
					double   ANGLE2,
					double   ANGLE3,
					double   ANGLE4,
					double   ROU1,
					double   ROU2,
					double   ROU3,
					double   ROU4,
					double	PSTRAIN1,
					double	PSTRAIN2,
					double   FPC,
					double   FY1,
					double	FY2,
					double   E,
					double   EPSC0);
				~RAFourSteelPCPlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};

			public ref class RAFourSteelRCPlaneStressWrapper : NDMaterialWrapper
			{
			public:
				RAFourSteelRCPlaneStressWrapper(int      tag,
					double   RHO,
					UniaxialMaterialWrapper^ s1,
					UniaxialMaterialWrapper^ s2,
					UniaxialMaterialWrapper^ s3,
					UniaxialMaterialWrapper^ s4,
					UniaxialMaterialWrapper^ c1,
					UniaxialMaterialWrapper^ c2,
					double   ANGLE1,
					double   ANGLE2,
					double   ANGLE3,
					double   ANGLE4,
					double   ROU1,
					double   ROU2,
					double   ROU3,
					double   ROU4,
					double   FPC,
					double   FY,
					double   E,
					double   EPSC0);
				~RAFourSteelRCPlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};

			public ref class ReinforcedConcretePlaneStressWrapper : NDMaterialWrapper
			{
			public:
				ReinforcedConcretePlaneStressWrapper(int      tag,
					double   RHO,
					UniaxialMaterialWrapper^ s1,
					UniaxialMaterialWrapper^ s2,
					UniaxialMaterialWrapper^ c1,
					UniaxialMaterialWrapper^ c2,
					double   ANGLE1,
					double   ANGLE2,
					double   ROU1,
					double   ROU2,
					double   FPC,
					double   FY,
					double   E,
					double   EPSC0);
				~ReinforcedConcretePlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			};


#pragma endregion

#pragma region soilModels

			public ref class FluidSolidPorousMaterialWrapper : NDMaterialWrapper
			{
			public:
				FluidSolidPorousMaterialWrapper(int tag, int nd, NDMaterialWrapper^ soilMat,
					double combinedBulkModul, double atm);
				~FluidSolidPorousMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class MultiYieldSurfaceClayWrapper : NDMaterialWrapper
			{
			public:
				MultiYieldSurfaceClayWrapper(int tag,
					int nd,
					double rho,
					double refShearModul,
					double refBulkModul,
					double cohesi,
					double peakShearStra,
					double frictionAng,
					double refPress,
					double pressDependCoe,
					int   numberOfYieldSurf,
					array<double>^ gredu);
				~MultiYieldSurfaceClayWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PressureDependMultiYieldWrapper : NDMaterialWrapper
			{
			public:
				PressureDependMultiYieldWrapper(int tag,
					int nd,
					double rho,
					double refShearModul,
					double refBulkModul,
					double frictionAng,
					double peakShearStra,
					double refPress,
					double pressDependCoe,
					double phaseTransformAngle,
					double contractionParam1,
					double dilationParam1,
					double dilationParam2,
					double liquefactionParam1,
					double liquefactionParam2,
					double liquefactionParam4,
					int   numberOfYieldSurf,
					array<double>^ gredu,
					double e,
					double volLimit1,
					double volLimit2,
					double volLimit3,
					double atm,
					double cohesi,
					double hv,
					double pv);
				~PressureDependMultiYieldWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PressureDependMultiYield02Wrapper : NDMaterialWrapper
			{
			public:
				PressureDependMultiYield02Wrapper(int tag,
					int nd,
					double rho,
					double refShearModul,
					double refBulkModul,
					double frictionAng,
					double peakShearStra,
					double refPress,
					double pressDependCoe,
					double phaseTransformAngle,
					double contractionParam1,
					double contractionParam3,
					double dilationParam1,
					double dilationParam3,
					int   numberOfYieldSurf,
					array<double>^ gredu,
					double contractionParam2,
					double dilationParam2,
					double liquefactionParam1,
					double liquefactionParam2,
					double e,
					double volLimit1,
					double volLimit2,
					double volLimit3,
					double atm,
					double cohesi,
					double hv,
					double pv);
				~PressureDependMultiYield02Wrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PressureIndependMultiYieldWrapper : NDMaterialWrapper
			{
			public:
				PressureIndependMultiYieldWrapper(int tag,
					int nd,
					double rho,
					double refShearModul,
					double refBulkModul,
					double cohesi,
					double peakShearStra,
					double frictionAng,
					double refPress,
					double pressDependCoe,
					int   numberOfYieldSurf,
					array<double>^ gredu);
				~PressureIndependMultiYieldWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

#pragma endregion

#pragma region stressDensityModel

			//public ref class StressDensityModelWrapper : NDMaterialWrapper
			//{
			//public:
			//	StressDensityModelWrapper(int tag, int classTag, double constDensity,
			//		// SD model  parameters		
			//		double initialVoidRatio, double constA, double exponentN,
			//		double poissonRatio, double constAlpha1, double constBeta1,
			//		double constAlpha2, double constBeta2, double constAlpha3,
			//		double constBeta3, double constDegradation, double constMumin,
			//		double constMucyclic, double constDilatancyStrain,
			//		double constMumax, double constPatm);

			//	StressDensityModelWrapper(int tag, int classTag, double constDensity,
			//		// SD model  parameters		
			//		double initialVoidRatio, double constA, double exponentN,
			//		double poissonRatio, double constAlpha1, double constBeta1,
			//		double constAlpha2, double constBeta2, double constAlpha3,
			//		double constBeta3, double constDegradation, double constMumin,
			//		double constMucyclic, double constDilatancyStrain,
			//		double constMumax, double constPatm, double constsslvoidatP1, double constsslvoidatP2, double constsslvoidatP3,
			//		double constsslvoidatP4, double constsslvoidatP5, double constsslvoidatP6,
			//		double constsslvoidatP7, double constsslvoidatP8, double constsslvoidatP9,
			//		double constsslvoidatP10,
			//		// hydrostatic state line void ratio
			//		double consthslvoid,
			//		// reference pressures
			//		double constP1, double constP2, double constP3, double constP4,
			//		double constP5, double constP6, double constP7, double constP8,
			//		double constP9, double constP10,
			//		// offset of the failure surface
			//		double constRxx, double constRyy, double constRzz,
			//		double constRxy, double constRyz, double constRzx);
			//	~StressDensityModelWrapper() {
			//		if (_NDMaterial != 0)
			//			delete _NDMaterial;
			//	};
			//private:

			//};

			//public ref class StressDensityModel2DWrapper : NDMaterialWrapper
			//{
			//public:
			//	StressDensityModel2DWrapper(int tag, double constDensity,
			//		// SD model  parameters
			//		double initialVoidRatio, double constA, double exponentN, double poissonRatio,
			//		double constAlpha1, double constBeta1, double constAlpha2, double constBeta2,
			//		double constAlpha3, double constBeta3, double constDegradation, double constMumin,
			//		double constMucyclic, double constDilatancyStrain, double constMumax, double constPatm,
			//		// steady state line void ratio
			//		double constsslvoidatP1, double constsslvoidatP2, double constsslvoidatP3,
			//		double constsslvoidatP4, double constsslvoidatP5, double constsslvoidatP6,
			//		double constsslvoidatP7, double constsslvoidatP8, double constsslvoidatP9,
			//		double constsslvoidatP10,
			//		// hydrostatic state line void ratio
			//		double consthslvoid,
			//		// reference pressures 
			//		double constP1, double constP2, double constP3, double constP4, double constP5,
			//		double constP6, double constP7, double constP8, double constP9, double constP10,
			//		// offset of the failure surface
			//		double constRxx, double constRyy, double constRxy);
			//	~StressDensityModel2DWrapper() {
			//		if (_NDMaterial != 0)
			//			delete _NDMaterial;
			//	};
			//private:

			//};

			//public ref class StressDensityModel3DWrapper : NDMaterialWrapper
			//{
			//public:
			//	StressDensityModel3DWrapper(int tag, double constDensity,
			//		// SD model  parameters
			//		double initialVoidRatio, double constA, double exponentN, double poissonRatio,
			//		double constAlpha1, double constBeta1, double constAlpha2, double constBeta2,
			//		double constAlpha3, double constBeta3, double constDegradation, double constMumin,
			//		double constMucyclic, double constDilatancyStrain, double constMumax, double constPatm,
			//		// steady state line void ratio
			//		double constsslvoidatP1, double constsslvoidatP2, double constsslvoidatP3,
			//		double constsslvoidatP4, double constsslvoidatP5, double constsslvoidatP6,
			//		double constsslvoidatP7, double constsslvoidatP8, double constsslvoidatP9,
			//		double constsslvoidatP10,
			//		// hydrostatic state line void ratio
			//		double consthslvoid,
			//		// reference pressures 
			//		double constP1, double constP2, double constP3, double constP4, double constP5,
			//		double constP6, double constP7, double constP8, double constP9, double constP10,
			//		// offset of the failure surface
			//		double constRxx, double constRyy, double constRzz, double constRxy,
			//		double constRyz, double constRzx);
			//	~StressDensityModel3DWrapper() {
			//		if (_NDMaterial != 0)
			//			delete _NDMaterial;
			//	};
			//private:

			//};

#pragma endregion

#pragma region wrappers
			public ref class BeamFiberMaterialWrapper : NDMaterialWrapper
			{
			public:
				BeamFiberMaterialWrapper(int tag, NDMaterialWrapper^ theMat);
				~BeamFiberMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class BeamFiberMaterial2dWrapper : NDMaterialWrapper
			{
			public:
				BeamFiberMaterial2dWrapper(int tag, NDMaterialWrapper^ theMat);
				~BeamFiberMaterial2dWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlaneStrainMaterialWrapper : NDMaterialWrapper
			{
			public:
				PlaneStrainMaterialWrapper(int tag, NDMaterialWrapper^ theMat);
				~PlaneStrainMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlaneStressLayeredMaterialWrapper : NDMaterialWrapper
			{
			public:
				PlaneStressLayeredMaterialWrapper(int tag,
					int iLayers,
					array<double>^ thickness,
					array<NDMaterialWrapper^>^ fibers);

				~PlaneStressLayeredMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlaneStressMaterialWrapper : NDMaterialWrapper
			{
			public:
				PlaneStressMaterialWrapper(int tag, NDMaterialWrapper^ the3DMaterial);
				~PlaneStressMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlaneStressRebarMaterialWrapper : NDMaterialWrapper
			{
			public:
				PlaneStressRebarMaterialWrapper(int tag, UniaxialMaterialWrapper^ uniMat,
					double ang);
				~PlaneStressRebarMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlaneStressSimplifiedJ2Wrapper : NDMaterialWrapper
			{
			public:
				PlaneStressSimplifiedJ2Wrapper(int tag, int nd, NDMaterialWrapper^ uniMat);
				~PlaneStressSimplifiedJ2Wrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlaneStressUserMaterialWrapper : NDMaterialWrapper
			{
			public:
				PlaneStressUserMaterialWrapper(int tag, int istatevs, int iprops, array<double>^ props);
				~PlaneStressUserMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlateFiberMaterialWrapper : NDMaterialWrapper
			{
			public:
				PlateFiberMaterialWrapper(int    tag,
					NDMaterialWrapper^ the3DMaterial);
				~PlateFiberMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlateFiberMaterialThermalWrapper : NDMaterialWrapper
			{
			public:
				PlateFiberMaterialThermalWrapper(int    tag,
					NDMaterialWrapper^ the3DMaterial);
				~PlateFiberMaterialThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlateFromPlaneStressMaterialWrapper : NDMaterialWrapper
			{
			public:
				PlateFromPlaneStressMaterialWrapper(int    tag,
					NDMaterialWrapper^ ndMat, double g);
				~PlateFromPlaneStressMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlateFromPlaneStressMaterialThermalWrapper : NDMaterialWrapper
			{
			public:
				PlateFromPlaneStressMaterialThermalWrapper(int    tag,
					NDMaterialWrapper^ ndMat, double g);
				~PlateFromPlaneStressMaterialThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlateRebarMaterialWrapper : NDMaterialWrapper
			{
			public:
				PlateRebarMaterialWrapper(int    tag,
					UniaxialMaterialWrapper^ uniMat, double ang);
				~PlateRebarMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlateRebarMaterialThermalWrapper : NDMaterialWrapper
			{
			public:
				PlateRebarMaterialThermalWrapper(int    tag,
					UniaxialMaterialWrapper^ uniMat, double ang);
				~PlateRebarMaterialThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};
#pragma endregion

			public ref class ElasticIsotropicMaterialWrapper : NDMaterialWrapper
			{
			public:
				ElasticIsotropicMaterialWrapper(int tag, double E, double nu, double rho);
				~ElasticIsotropicMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class AcousticMediumWrapper : NDMaterialWrapper
			{
			public:
				AcousticMediumWrapper(int tag, int classTag, double k, double rho, double gamma);
				~AcousticMediumWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class BoundingCamClayWrapper : NDMaterialWrapper
			{
			public:
				BoundingCamClayWrapper(int tag, int classTag, double massDen, double C, double bulk, double OCR,
					double mu_o, double Alpha, double lambda, double h, double m);
				~BoundingCamClayWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class BoundingCamClay3DWrapper : NDMaterialWrapper
			{
			public:
				BoundingCamClay3DWrapper(int tag, double mDen, double c, double bulk, double OCR, double mu_o,
					double alpha, double lambda, double h, double m);
				~BoundingCamClay3DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class BoundingCamClayPlaneStrainWrapper : NDMaterialWrapper
			{
			public:
				BoundingCamClayPlaneStrainWrapper(int tag, double mDen, double c, double bulk, double OCR, double mu_o,
					double alpha, double lambda, double h, double m);
				~BoundingCamClayPlaneStrainWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class DruckerPragerThermalWrapper : NDMaterialWrapper
			{
			public:
				DruckerPragerThermalWrapper(int tag, int classTag, double bulk, double shear,
					double s_y, double r, double r_bar, double Kinfinity, double Kinit,
					double d1, double d2, double H, double t, double massDen, double atm);
				~DruckerPragerThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class DruckerPragerWrapper : NDMaterialWrapper
			{
			public:
				DruckerPragerWrapper(int tag, int classTag, double bulk, double shear,
					double s_y, double r, double r_bar, double Kinfinity, double Kinit,
					double d1, double d2, double H, double t, double massDen, double atm);
				DruckerPragerWrapper() {};
				~DruckerPragerWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class DruckerPrager3DWrapper : DruckerPragerWrapper
			{
			public:
				DruckerPrager3DWrapper(int tag, double bulk, double shear,
					double s_y, double r, double r_bar, double Kinfinity, double Kinit,
					double d1, double d2, double H, double t, double massDen, double atm);
				~DruckerPrager3DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class CycLiqCPSPWrapper : NDMaterialWrapper
			{
			public:
				CycLiqCPSPWrapper(int    tag,
					int classTag,
					double G01,
					double kappa1,
					double h1,
					double Mfc1,       //critical state
					double dre11,
					double dre21,
					double rdr1,
					double eta1,
					double dir1,
					double lamdac1,
					double ksi1,
					double e01,
					double nb1,
					double nd1,
					double ein1,      //initial void ratio
					double rho1);
				CycLiqCPSPWrapper() {};
				~CycLiqCPSPWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class CycLiqCPWrapper : NDMaterialWrapper
			{
			public:
				CycLiqCPWrapper(int    tag,
					int classTag,
					double G01,
					double kappa1,
					double h1,
					double Mfc1,       //critical state
					double dre11,
					double Mdc1,
					double dre21,
					double rdr1,
					double eta1,
					double dir1,
					double ein1,      //initial void ratio
					double rho1);
				CycLiqCPWrapper() {};
				~CycLiqCPWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class CycLiqCPSP3DWrapper : CycLiqCPSPWrapper
			{
			public:
				CycLiqCPSP3DWrapper(int    tag,
					int classTag,
					double G01,
					double kappa1,
					double h1,
					double Mfc1,       //critical state
					double dre11,
					double dre21,
					double rdr1,
					double eta1,
					double dir1,
					double lamdac1,
					double ksi1,
					double e01,
					double nb1,
					double nd1,
					double ein1,      //initial void ratio
					double rho1);
				~CycLiqCPSP3DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class CycLiqCPPlaneStrainWrapper : CycLiqCPWrapper
			{
			public:
				CycLiqCPPlaneStrainWrapper(int    tag,
					double G01,
					double kappa1,
					double h1,
					double Mfc1,       //critical state
					double dre11,
					double Mdc1,
					double dre21,
					double rdr1,
					double eta1,
					double dir1,
					double ein1,      //initial void ratio
					double rho1);
				~CycLiqCPPlaneStrainWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class CycLiqCP3DWrapper : CycLiqCPWrapper
			{
			public:
				CycLiqCP3DWrapper(int    tag,
					double G01,
					double kappa1,
					double h1,
					double Mfc1,       //critical state
					double dre11,
					double Mdc1,
					double dre21,
					double rdr1,
					double eta1,
					double dir1,
					double ein1,      //initial void ratio
					double rho1);
				~CycLiqCP3DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ContactMaterial3DWrapper : NDMaterialWrapper
			{
			public:
				ContactMaterial3DWrapper(int tag, double mu, double G, double c, double t);
				~ContactMaterial3DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ContactMaterial2DWrapper : NDMaterialWrapper
			{
			public:
				ContactMaterial2DWrapper(int tag, double mu, double G, double c, double t);
				~ContactMaterial2DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ConcreteSWrapper : NDMaterialWrapper
			{
			public:
				ConcreteSWrapper(int tag, double rE, double rnu, double rfc, double rft, double rEs);
				~ConcreteSWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class CapPlasticityWrapper : NDMaterialWrapper
			{
			public:
				CapPlasticityWrapper(int    tag,
					double G,
					double K,
					double rho,
					double X,
					double D,
					double W,
					double R,
					double lambda,
					double theta,
					double beta,
					double alpha,
					double T,
					int ndm,
					double pTol_k);
				~CapPlasticityWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class CycLiqCPSPPlaneStrainWrapper : CycLiqCPSPWrapper
			{
			public:
				CycLiqCPSPPlaneStrainWrapper(int    tag,
					double G01,
					double kappa1,
					double h1,
					double Mfc1,       //critical state
					double dre11,
					double Mdc1,
					double dre21,
					double rdr1,
					double eta1,
					double dir1,
					double lamdac1,
					double ksi1,
					double e01,
					double nb1,
					double nd1,
					double ein1,      //initial void ratio
					double rho1);
				~CycLiqCPSPPlaneStrainWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			

			public ref class DruckerPrager3DThermalWrapper : NDMaterialWrapper
			{
			public:
				DruckerPrager3DThermalWrapper(int tag, double bulk, double shear,
					double s_y, double r, double r_bar, double Kinfinity, double Kinit,
					double d1, double d2, double H, double t, double massDen, double atm);
				~DruckerPrager3DThermalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class DruckerPragerPlaneStrainWrapper : DruckerPragerWrapper
			{
			public:
				DruckerPragerPlaneStrainWrapper(int tag, double bulk, double shear,
					double s_y, double r, double r_bar, double Kinfinity, double Kinit,
					double d1, double d2, double H, double t, double massDens, double atm);
				~DruckerPragerPlaneStrainWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ElasticOrthotropicMaterialWrapper : DruckerPragerWrapper
			{
			public:
				ElasticOrthotropicMaterialWrapper(int tag, int classTag,
					double Ex, double Ey, double Ez,
					double vxy, double vyz, double vzx,
					double Gxy, double Gyz, double Gz, double rho);
				ElasticOrthotropicMaterialWrapper() {};
				~ElasticOrthotropicMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ElasticOrthotropicThreeDimensionalWrapper : ElasticOrthotropicMaterialWrapper
			{
			public:
				ElasticOrthotropicThreeDimensionalWrapper(int tag, double Ex, double Ey,
					double Ez, double vxy, double vyz, double vzx, double Gxy, double Gyz,
					double Gzx, double rho);
				~ElasticOrthotropicThreeDimensionalWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class FSAMWrapper : NDMaterialWrapper
			{
			public:
				FSAMWrapper(int tag,             
					double RHO,             
					UniaxialMaterialWrapper^ s1,   
					UniaxialMaterialWrapper^ s2,   
					UniaxialMaterialWrapper^ c1,   
					UniaxialMaterialWrapper^ c2,   
					UniaxialMaterialWrapper^ cA1,  
					UniaxialMaterialWrapper^ cA2,  
					UniaxialMaterialWrapper^ cB1,  
					UniaxialMaterialWrapper^ cB2,  
					double ROUX,            
					double ROUY,			
					double NU,				
					double ALFADOW);
				~FSAMWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class InitialStateAnalysisWrapperWrapper : NDMaterialWrapper
			{
			public:
				InitialStateAnalysisWrapperWrapper(int tag, NDMaterialWrapper^ mainMat, int ndim);
				~InitialStateAnalysisWrapperWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class InitStressNDMaterialWrapper : NDMaterialWrapper
			{
			public:
				InitStressNDMaterialWrapper(int tag, NDMaterialWrapper^ material,  VectorWrapper^ sigInit, int ndim);
				~InitStressNDMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class LinearCapWrapper : NDMaterialWrapper
			{
			public:
				LinearCapWrapper(int    tag,
					double G,
					double K,
					double rho,
					double theta,
					double alpha,
					double T,
					int ndm,
					double pTol_k);
				~LinearCapWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ManzariDafaliasWrapper : NDMaterialWrapper
			{
			public:
				ManzariDafaliasWrapper(int tag, int classTag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
					double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen,
					int integrationScheme , int tangentType, int JacoType, double TolF, double TolR);
				ManzariDafaliasWrapper() {};
				~ManzariDafaliasWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ManzariDafalias3DWrapper : ManzariDafaliasWrapper
			{
			public:
				ManzariDafalias3DWrapper(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
					double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme ,
					int tangentType, int JacoType , double TolF, double TolR );
				~ManzariDafalias3DWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ManzariDafalias3DROWrapper : ManzariDafaliasWrapper
			{
			public:
				ManzariDafalias3DROWrapper(int tag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c,
					double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
					double z_max, double cz, double mDen, double kappa , int integrationScheme , int tangentType , int JacoType ,
					double TolF , double TolR );
				~ManzariDafalias3DROWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ManzariDafaliasPlaneStrainWrapper : ManzariDafaliasWrapper
			{
			public:
				ManzariDafaliasPlaneStrainWrapper(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
					double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme ,
					int tangentType , int JacoType , double TolF , double TolR);
				~ManzariDafaliasPlaneStrainWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ManzariDafaliasPlaneStrainROWrapper : ManzariDafaliasWrapper
			{
			public:
				ManzariDafaliasPlaneStrainROWrapper(int tag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c,
					double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
					double z_max, double cz, double mDen, double kappa, int integrationScheme, int tangentType , int JacoType ,
					double TolF , double TolR );
				~ManzariDafaliasPlaneStrainROWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class ManzariDafaliasROWrapper : ManzariDafaliasWrapper
			{
			public:
				ManzariDafaliasROWrapper(int tag, int classTag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c,
					double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
					double z_max, double cz, double mDen, double kappa, int integrationScheme, int tangentType, int JacoType,
					double TolF, double TolR);
				~ManzariDafaliasROWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlasticDamageConcrete3dWrapper : NDMaterialWrapper
			{
			public:
				PlasticDamageConcrete3dWrapper(int tag,
					double E,
					double nu,
					double ft,
					double fc,
					double beta ,
					double Ap ,
					double An ,
					double Bn );
				~PlasticDamageConcrete3dWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PlasticDamageConcretePlaneStressWrapper : NDMaterialWrapper
			{
			public:
				PlasticDamageConcretePlaneStressWrapper(int tag,
					double E,
					double nu,
					double ft,
					double fc,
					double beta ,
					double Ap ,
					double An ,
					double Bn );
				~PlasticDamageConcretePlaneStressWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PM4SandWrapper : NDMaterialWrapper
			{
			public:
				PM4SandWrapper(int tag, int classTag, double Dr, double G0, double hp0, double mDen, double P_atm, double h0, double emax,
					double emin, double nb, double nd, double Ado, double z_max, double cz,
					double ce, double phi_cv, double nu, double Cgd, double Cdr, double Ckaf, double Q,
					double R, double m, double Fsed_min, double p_sdeo, int integrationScheme, int tangentType,
					double TolF, double TolR);
				PM4SandWrapper(int tag, double Dr, double G0, double hp0, double mDen, double P_atm, double h0, double emax,
					double emin, double nb, double nd, double Ado, double z_max, double cz,
					double ce, double phi_cv, double nu, double Cgd, double Cdr, double Ckaf, double Q,
					double R, double m, double Fsed_min, double p_sdeo, int integrationScheme, int tangentType,
					double TolF, double TolR);
				~PM4SandWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class PM4SiltWrapper : NDMaterialWrapper
			{
			public:
				PM4SiltWrapper(int tag, int classTag, double Su, double Su_rate, double G0, double hpo, double mDen, double Fsu, double P_atm, double nu,
					double nG, double h0, double einit, double lambda, double phi_cv, double nbwet, double nbdry, double nd,
					double Ado, double ru_max, double z_max, double cz, double ce, double Cgd, double Ckaf, double m,
					double CG_consol, int integrationScheme, int tangentType, double TolF, double TolR);

				PM4SiltWrapper(int tag, double Su, double Su_rate, double G0, double hpo, double mDen, double Fsu, double P_atm, double nu,
					double nG, double h0, double einit, double lambda, double phi_cv, double nbwet, double nbdry, double nd,
					double Ado, double ru_max, double z_max, double cz, double ce, double Cgd, double Ckaf, double m,
					double CG_consol, int integrationScheme, int tangentType, double TolF, double TolR);

				~PM4SiltWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public ref class SimplifiedJ2Wrapper : NDMaterialWrapper
			{
			public:
				SimplifiedJ2Wrapper(int tag,
					int nd,
					double G,
					double K,
					double sigmaY0,
					double H_kin,
					double H_iso);

				~SimplifiedJ2Wrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;
				};
			private:

			};

			public delegate NDMaterialWrapper^ NDMD_GetCopy();
			delegate NDMaterial* UNDMD_GetCopy();

			public delegate NDMaterialWrapper^ NDMD_GetCopy_Type(String^ type);
			delegate NDMaterial* UNDMD_GetCopy_Type(const char *type);

			public delegate void NDMD_Print(int);
			delegate void UNDMD_Print(int);

			public delegate double NDMD_GetRho(void);
			delegate double UNDMD_GetRho(void);

			public delegate int NDMD_SetTrialStrain_V(VectorWrapper^ v);
			delegate int UNDMD_SetTrialStrain_V(const Vector &v);

			public delegate int NDMD_SetTrialStrain_VR(VectorWrapper^ v, VectorWrapper^ r);
			delegate int UNDMD_SetTrialStrain_VR(const Vector &v, const Vector &r);

			public delegate int NDMD_SetTrialStrainIncr_V(VectorWrapper^ v);
			delegate int UNDMD_SetTrialStrainIncr_V(const Vector &v);

			public delegate int NDMD_SetTrialStrainIncr_VR(VectorWrapper^ v, VectorWrapper^ r);
			delegate int UNDMD_SetTrialStrainIncr_VR(const Vector &v, const Vector &r);

			public delegate MatrixWrapper^ NDMD_GetTangent(void);
			delegate Matrix * UNDMD_GetTangent(void);

			public delegate MatrixWrapper^ NDMD_GetInitialTangent(void);
			delegate Matrix * UNDMD_GetInitialTangent(void);

			public delegate VectorWrapper^ NDMD_GetStress(void);
			delegate Vector * UNDMD_GetStress(void);

			public delegate VectorWrapper^ NDMD_GetStrain(void);
			delegate Vector * UNDMD_GetStrain(void);

			public delegate int NDMD_RevertToStart(void);
			delegate int UNDMD_RevertToStart(void);

			public delegate int NDMD_CommitState(void);
			delegate int UNDMD_CommitState(void);

			public delegate int NDMD_RevertToLastCommit(void);
			delegate int UNDMD_RevertToLastCommit(void);

			public delegate String^ NDMD_GetType(void);
			delegate const char* UNDMD_GetType(void);

			public delegate int NDMD_GetOrder(void);
			delegate int UNDMD_GetOrder(void);

			public ref class ExternalNDMaterialWrapper : NDMaterialWrapper
			{
				
			public:
				ExternalNDMaterialWrapper(int tag);

				~ExternalNDMaterialWrapper() {
					if (_NDMaterial != 0)
						delete _NDMaterial;

					undmd_SetTrialStrain_V_GC.Free();
					undmd_SetTrialStrain_VR_GC.Free();
					undmd_SetTrialStrainIncr_V_GC.Free();
					undmd_SetTrialStrainIncr_VR_GC.Free();
					undmd_GetTangent_GC.Free();
					undmd_GetInitialTangent_GC.Free();
					undmd_GetStress_GC.Free();
					undmd_GetStrain_GC.Free();
					undmd_RevertToStart_GC.Free();
					undmd_CommitState_GC.Free();
					undmd_RevertToLastCommit_GC.Free();
					undmd_GetType_GC.Free();
					undmd_GetOrder_GC.Free();
					undmd_GetCopy_GC.Free();
					undmd_GetCopy_Type_GC.Free();
					undmd_Print_GC.Free();
					undmd_GetRho_GC.Free();
				};

				void SetLinks(
					NDMD_SetTrialStrain_V^ _NDMD_SetTrialStrain_V,
					NDMD_SetTrialStrain_VR^ _NDMD_SetTrialStrain_VR,
					NDMD_SetTrialStrainIncr_V^ _NDMD_SetTrialStrainIncr_V,
					NDMD_SetTrialStrainIncr_VR^ _NDMD_SetTrialStrainIncr_VR,
					NDMD_GetTangent^ _NDMD_GetTangent,
					NDMD_GetInitialTangent^ _NDMD_GetInitialTangent,
					NDMD_GetStress^ _NDMD_GetStress,
					NDMD_GetStrain^ _NDMD_GetStrain,
					NDMD_RevertToStart^ _NDMD_RevertToStart,
					NDMD_CommitState^ _NDMD_CommitState,
					NDMD_RevertToLastCommit^ _NDMD_RevertToLastCommit,
					NDMD_GetType^ _NDMD_GetType,
					NDMD_GetOrder^ _NDMD_GetOrder,
					NDMD_GetCopy^ _NDMD_GetCopy,
					NDMD_GetCopy_Type^ _NDMD_GetCopy_Type,
					NDMD_Print^ _NDMD_Print,
					NDMD_GetRho^ _NDMD_GetRho
				) {
					this->_NDMD_SetTrialStrain_V = _NDMD_SetTrialStrain_V;
					UNDMD_SetTrialStrain_V^ undmd_SetTrialStrain_V = gcnew UNDMD_SetTrialStrain_V(this, &ExternalNDMaterialWrapper::SetTrialStrain_V_Method);
					undmd_SetTrialStrain_V_GC = GCHandle::Alloc(undmd_SetTrialStrain_V);
					IntPtr ip = Marshal::GetFunctionPointerForDelegate(undmd_SetTrialStrain_V);
					NDM_SetTrialStrain_V setTrialStrain_V = static_cast<NDM_SetTrialStrain_V>(ip.ToPointer());

					this->_NDMD_SetTrialStrain_VR = _NDMD_SetTrialStrain_VR;
					UNDMD_SetTrialStrain_VR^ undmd_SetTrialStrain_VR = gcnew UNDMD_SetTrialStrain_VR(this, &ExternalNDMaterialWrapper::SetTrialStrain_VR_Method);
					undmd_SetTrialStrain_VR_GC = GCHandle::Alloc(undmd_SetTrialStrain_VR);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_SetTrialStrain_VR);
					NDM_SetTrialStrain_VR SetTrialStrain_VR = static_cast<NDM_SetTrialStrain_VR>(ip.ToPointer());

					this->_NDMD_SetTrialStrainIncr_V = _NDMD_SetTrialStrainIncr_V;
					UNDMD_SetTrialStrainIncr_V^ undmd_SetTrialStrainIncr_V = gcnew UNDMD_SetTrialStrainIncr_V(this, &ExternalNDMaterialWrapper::SetTrialStrainIncr_V_Method);
					undmd_SetTrialStrainIncr_V_GC = GCHandle::Alloc(undmd_SetTrialStrainIncr_V);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_SetTrialStrainIncr_V);
					NDM_SetTrialStrainIncr_V SetTrialStrainIncr_V = static_cast<NDM_SetTrialStrainIncr_V>(ip.ToPointer());

					this->_NDMD_SetTrialStrainIncr_VR = _NDMD_SetTrialStrainIncr_VR;
					UNDMD_SetTrialStrainIncr_VR^ undmd_SetTrialStrainIncr_VR = gcnew UNDMD_SetTrialStrainIncr_VR(this, &ExternalNDMaterialWrapper::SetTrialStrainIncr_VR_Method);
					undmd_SetTrialStrainIncr_VR_GC = GCHandle::Alloc(undmd_SetTrialStrainIncr_VR);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_SetTrialStrainIncr_VR);
					NDM_SetTrialStrainIncr_VR SetTrialStrainIncr_VR = static_cast<NDM_SetTrialStrainIncr_VR>(ip.ToPointer());

					this->_NDMD_GetTangent = _NDMD_GetTangent;
					UNDMD_GetTangent^ undmd_GetTangent = gcnew UNDMD_GetTangent(this, &ExternalNDMaterialWrapper::GetTangent_Method);
					undmd_GetTangent_GC = GCHandle::Alloc(undmd_GetTangent);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetTangent);
					NDM_GetTangent GetTangent = static_cast<NDM_GetTangent>(ip.ToPointer());

					this->_NDMD_GetInitialTangent = _NDMD_GetInitialTangent;
					UNDMD_GetInitialTangent^ undmd_GetInitialTangent = gcnew UNDMD_GetInitialTangent(this, &ExternalNDMaterialWrapper::GetInitialTangent_Method);
					undmd_GetInitialTangent_GC = GCHandle::Alloc(undmd_GetInitialTangent);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetInitialTangent);
					NDM_GetInitialTangent GetInitialTangent = static_cast<NDM_GetInitialTangent>(ip.ToPointer());

					this->_NDMD_GetStress = _NDMD_GetStress;
					UNDMD_GetStress^ undmd_GetStress = gcnew UNDMD_GetStress(this, &ExternalNDMaterialWrapper::GetStress_Method);
					undmd_GetStress_GC = GCHandle::Alloc(undmd_GetStress);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetStress);
					NDM_GetStress GetStress = static_cast<NDM_GetStress>(ip.ToPointer());

					this->_NDMD_GetStrain = _NDMD_GetStrain;
					UNDMD_GetStrain^ undmd_GetStrain = gcnew UNDMD_GetStrain(this, &ExternalNDMaterialWrapper::GetStrain_Method);
					undmd_GetStrain_GC = GCHandle::Alloc(undmd_GetStrain);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetStrain);
					NDM_GetStrain GetStrain = static_cast<NDM_GetStrain>(ip.ToPointer());

					this->_NDMD_RevertToStart = _NDMD_RevertToStart;
					UNDMD_RevertToStart^ undmd_RevertToStart = gcnew UNDMD_RevertToStart(this, &ExternalNDMaterialWrapper::RevertToStart_Method);
					undmd_RevertToStart_GC = GCHandle::Alloc(undmd_RevertToStart);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_RevertToStart);
					NDM_RevertToStart RevertToStart = static_cast<NDM_RevertToStart>(ip.ToPointer());

					this->_NDMD_CommitState = _NDMD_CommitState;
					UNDMD_CommitState^ undmd_CommitState = gcnew UNDMD_CommitState(this, &ExternalNDMaterialWrapper::CommitState_Method);
					undmd_CommitState_GC = GCHandle::Alloc(undmd_CommitState);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_CommitState);
					NDM_CommitState CommitState = static_cast<NDM_CommitState>(ip.ToPointer());

					this->_NDMD_RevertToLastCommit = _NDMD_RevertToLastCommit;
					UNDMD_RevertToLastCommit^ undmd_RevertToLastCommit = gcnew UNDMD_RevertToLastCommit(this, &ExternalNDMaterialWrapper::RevertToLastCommit_Method);
					undmd_RevertToLastCommit_GC = GCHandle::Alloc(undmd_RevertToLastCommit);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_RevertToLastCommit);
					NDM_RevertToLastCommit RevertToLastCommit = static_cast<NDM_RevertToLastCommit>(ip.ToPointer());

					this->_NDMD_GetType = _NDMD_GetType;
					UNDMD_GetType^ undmd_GetType = gcnew UNDMD_GetType(this, &ExternalNDMaterialWrapper::GetType_Method);
					undmd_GetType_GC = GCHandle::Alloc(undmd_GetType);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetType);
					NDM_GetType GetType = static_cast<NDM_GetType>(ip.ToPointer());

					this->_NDMD_GetOrder = _NDMD_GetOrder;
					UNDMD_GetOrder^ undmd_GetOrder = gcnew UNDMD_GetOrder(this, &ExternalNDMaterialWrapper::GetOrder_Method);
					undmd_GetOrder_GC = GCHandle::Alloc(undmd_GetOrder);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetOrder);
					NDM_GetOrder GetOrder = static_cast<NDM_GetOrder>(ip.ToPointer());

					this->_NDMD_GetCopy = _NDMD_GetCopy;
					UNDMD_GetCopy^ undmd_GetCopy = gcnew UNDMD_GetCopy(this, &ExternalNDMaterialWrapper::GetCopy_Method);
					undmd_GetCopy_GC = GCHandle::Alloc(undmd_GetCopy);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetCopy);
					NDM_GetCopy GetCopy = static_cast<NDM_GetCopy>(ip.ToPointer());

					this->_NDMD_GetCopy_Type = _NDMD_GetCopy_Type;
					UNDMD_GetCopy_Type^ undmd_GetCopy_Type = gcnew UNDMD_GetCopy_Type(this, &ExternalNDMaterialWrapper::GetCopy_Type_Method);
					undmd_GetCopy_Type_GC = GCHandle::Alloc(undmd_GetCopy_Type);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetCopy_Type);
					NDM_GetCopy_Type GetCopy_Type = static_cast<NDM_GetCopy_Type>(ip.ToPointer());

					this->_NDMD_Print = _NDMD_Print;
					UNDMD_Print^ undmd_Print = gcnew UNDMD_Print(this, &ExternalNDMaterialWrapper::Print_Method);
					undmd_Print_GC = GCHandle::Alloc(undmd_Print);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_Print);
					NDM_Print Print = static_cast<NDM_Print>(ip.ToPointer());

					this->_NDMD_GetRho = _NDMD_GetRho;
					UNDMD_GetRho^ undmd_GetRho = gcnew UNDMD_GetRho(this, &ExternalNDMaterialWrapper::GetRho_Method);
					undmd_GetRho_GC = GCHandle::Alloc(undmd_GetRho);
					ip = Marshal::GetFunctionPointerForDelegate(undmd_GetRho);
					NDM_GetRho GetRho = static_cast<NDM_GetRho>(ip.ToPointer());

					((ExternalNDMaterial*)_NDMaterial)->SetLinks(GetCopy, GetCopy_Type, Print, GetRho, setTrialStrain_V, SetTrialStrain_VR,
						SetTrialStrainIncr_V, SetTrialStrainIncr_VR, GetTangent, GetInitialTangent, GetStress, GetStrain, RevertToStart, CommitState,
						RevertToLastCommit, GetType, GetOrder);

				}

			internal :
				int SetTrialStrain_V_Method(const Vector &v) {
					VectorWrapper^ V = VectorWrapper::GetVectorWrapper(v);
					return this->_NDMD_SetTrialStrain_V(V);
				}

				int SetTrialStrain_VR_Method(const Vector &v, const Vector &r) {
					VectorWrapper^ V = VectorWrapper::GetVectorWrapper(v);
					VectorWrapper^ R = VectorWrapper::GetVectorWrapper(r);
					return this->_NDMD_SetTrialStrain_VR(V,R);
				}

				int SetTrialStrainIncr_V_Method(const Vector &v) {
					VectorWrapper^ V = VectorWrapper::GetVectorWrapper(v);
					return this->_NDMD_SetTrialStrainIncr_V(V);
				}

				int SetTrialStrainIncr_VR_Method(const Vector &v, const Vector &r) {
					VectorWrapper^ V = VectorWrapper::GetVectorWrapper(v);
					VectorWrapper^ R = VectorWrapper::GetVectorWrapper(r);
					return this->_NDMD_SetTrialStrainIncr_VR(V, R);
				}

				Matrix * GetTangent_Method() {
					return this->_NDMD_GetTangent()->_Matrix;
				}

				Matrix * GetInitialTangent_Method() {
					return this->_NDMD_GetInitialTangent()->_Matrix;
				}

				Vector * GetStress_Method() {
					return this->_NDMD_GetStress()->_Vector;
				}

				Vector * GetStrain_Method() {
					return this->_NDMD_GetStrain()->_Vector;
				}

				int RevertToStart_Method() {
					return this->_NDMD_RevertToStart();
				}

				int CommitState_Method() {
					return this->_NDMD_CommitState();
				}

				int RevertToLastCommit_Method() {
					return this->_NDMD_RevertToLastCommit();
				}

				const char* GetType_Method() {
					String^ strType = this->_NDMD_GetType();
					char* str2 = (char*)(void*)Marshal::StringToHGlobalAnsi(strType);
					return str2;
				}

				NDMaterial* GetCopy_Method() {
					NDMaterialWrapper^ mat = this->_NDMD_GetCopy();
					return mat->_NDMaterial;
				}

				NDMaterial* GetCopy_Type_Method(const char* type) {
					String^ str = gcnew String(type);
					NDMaterialWrapper^ mat = this->_NDMD_GetCopy_Type(str);
					return mat->_NDMaterial;
				}

				void Print_Method(int flag) {
					this->_NDMD_Print(flag);
				}

				double GetRho_Method() {
					return this->_NDMD_GetRho();
				}

				int GetOrder_Method() {
					return this->_NDMD_GetOrder();
				}

			private:
				GCHandle undmd_SetTrialStrain_V_GC;
				GCHandle undmd_SetTrialStrain_VR_GC;
				GCHandle undmd_SetTrialStrainIncr_V_GC;
				GCHandle undmd_SetTrialStrainIncr_VR_GC;
				GCHandle undmd_GetTangent_GC;
				GCHandle undmd_GetInitialTangent_GC;
				GCHandle undmd_GetStress_GC;
				GCHandle undmd_GetStrain_GC;
				GCHandle undmd_RevertToStart_GC;
				GCHandle undmd_CommitState_GC;
				GCHandle undmd_RevertToLastCommit_GC;
				GCHandle undmd_GetType_GC;
				GCHandle undmd_GetOrder_GC;
				GCHandle undmd_GetCopy_GC;
				GCHandle undmd_GetCopy_Type_GC;
				GCHandle undmd_Print_GC;
				GCHandle undmd_GetRho_GC;

				NDMD_GetRho^ _NDMD_GetRho;
				NDMD_SetTrialStrain_V^ _NDMD_SetTrialStrain_V;
				NDMD_SetTrialStrain_VR^ _NDMD_SetTrialStrain_VR;
				NDMD_SetTrialStrainIncr_V^ _NDMD_SetTrialStrainIncr_V;
				NDMD_SetTrialStrainIncr_VR^ _NDMD_SetTrialStrainIncr_VR;
				NDMD_GetTangent^ _NDMD_GetTangent;
				NDMD_GetInitialTangent^ _NDMD_GetInitialTangent;
				NDMD_GetStress^ _NDMD_GetStress;
				NDMD_GetStrain^ _NDMD_GetStrain;
				NDMD_RevertToStart^ _NDMD_RevertToStart;
				NDMD_CommitState^ _NDMD_CommitState;
				NDMD_RevertToLastCommit^ _NDMD_RevertToLastCommit;
				NDMD_GetType^ _NDMD_GetType;
				NDMD_GetOrder^ _NDMD_GetOrder;
				NDMD_GetCopy^ _NDMD_GetCopy;
				NDMD_GetCopy_Type^ _NDMD_GetCopy_Type;
				NDMD_Print^ _NDMD_Print;

			};

		}
	}
}
