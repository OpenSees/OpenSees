#pragma once


#include <generic/GenericClient.h>
#include <generic/GenericCopy.h>
#include <joint/BeamColumnJoint2d.h>
#include <joint/BeamColumnJoint3d.h>
#include <joint/ElasticTubularJoint.h>
#include <joint/Joint2D.h>
#include <joint/Joint3D.h>
#include <joint/MP_Joint2D.h>
#include <joint/MP_Joint3D.h>

#include <pyMacro/PY_Macro2D.h>

#include <tetrahedron\FourNodeTetrahedron.h>

#include <triangle\Tri31.h>

#include <twoNodeLink\TwoNodeLink.h>

#include <up-UCSD/BBarBrickUP.h>
#include <up-UCSD/BBarFourNodeQuadUP.h>
#include <up-UCSD/BrickUP.h>
#include <up-UCSD/FourNodeQuadUP.h>
#include <up-UCSD/Nine_Four_Node_QuadUP.h>

#include <UWelements/BeamContact2D.h>
#include <UWelements/BeamContact2Dp.h>
#include <UWelements/BeamContact3D.h>
#include <UWelements/BeamContact3Dp.h>
#include <UWelements/BeamEndContact3D.h>
#include <UWelements/BeamEndContact3Dp.h>
#include <UWelements/Brick8FiberOverlay.h>
#include <UWelements/PileToe3D.h>

#include <UWelements/QuadBeamEmbedContact.h>
#include <UWelements/SimpleContact2D.h>
#include <UWelements/SimpleContact3D.h>
#include <UWelements/SSPbrick.h>
#include <UWelements/SSPbrickUP.h>
#include <UWelements/SSPquad.h>
#include <UWelements/SSPquadUP.h>

#include <XMUelements/AC3D8HexWithSensitivity.h>
#include <XMUelements/ASI3D8QuadWithSensitivity.h>
#include <XMUelements/AV3D4QuadWithSensitivity.h>
#include <XMUelements/VS3D4QuadWithSensitivity.h>


#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../ElementWrapper.h"
#include "../crdTransf/CrdTransfWrapper.h"
#include "../../materials/section/SectionForceDeformationWrapper.h"
#include "../../materials/uniaxial/UniaxialMaterialWrapper.h"
#include "../../materials/ndmaterial/NDMaterialWrapper.h"
#include "../beamColumn/beamIntegration/BeamIntegrationWrapper.h"
#include "../../matrix/IDWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../../matrix/VectorWrapper.h"
#include "../../domains/domain/DomainWrapper.h"
#include "../../damage/DamageModelWrapper.h"


using namespace System;
using namespace OpenSees;
using namespace OpenSees::DamageModels;
using namespace OpenSees::Components;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Elements::CrdTransfs;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;

namespace OpenSees {
	namespace Elements {
		public enum class  PlaneElementType { PlaneStrain, PlaneStress };
#pragma region generic 
		public ref class GenericClientWrapper : ElementWrapper
		{

		public:
			GenericClientWrapper(int tag, IDWrapper^ nodes, IDWrapper^ dof,
				int port);
			GenericClientWrapper(int tag, IDWrapper^ nodes, IDWrapper^ dof,
				int port, String^ machineInetAddr,
				int ssl, int udp, int dataSize,
				int addRayleigh);

			GenericClientWrapper() {};
			~GenericClientWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class GenericCopyWrapper : ElementWrapper
		{

		public:
			GenericCopyWrapper(int tag, IDWrapper^ nodes, int srcTag);
			GenericCopyWrapper() {};
			~GenericCopyWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};
#pragma endregion

#pragma region joint

		public ref class MP_Joint2DWrapper : MP_ConstraintWrapper
		{
		public:

			MP_Joint2DWrapper(DomainWrapper^ theDomain, int nodeRetain, int nodeConstr,
				int Maindof, int fixedend, int LrgDsp);

			MP_Joint2DWrapper() {};
			~MP_Joint2DWrapper() {
				if (_MP_Constraint != 0)
					delete _MP_Constraint;
			};
		internal:
			MP_Constraint * _MP_Constraint;
		};

		public ref class MP_Joint3DWrapper : MP_ConstraintWrapper
		{
		public:

			MP_Joint3DWrapper(DomainWrapper^ theDomain, int nodeRetain, int nodeConstr,
				int nodeRot, int Rotdof, int nodeDisp, int Dispdof, int LrgDsp);

			MP_Joint3DWrapper() {};
			~MP_Joint3DWrapper() {
				if (_MP_Constraint != 0)
					delete _MP_Constraint;
			};
		internal:
			MP_Constraint * _MP_Constraint;
		};

		public ref class BeamColumnJoint2dWrapper : ElementWrapper
		{

		public:
			BeamColumnJoint2dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
				UniaxialMaterialWrapper^ theMat1, UniaxialMaterialWrapper^ theMat2,
				UniaxialMaterialWrapper^ theMat3, UniaxialMaterialWrapper^ theMat4,
				UniaxialMaterialWrapper^ theMat5, UniaxialMaterialWrapper^ theMat6,
				UniaxialMaterialWrapper^ theMat7, UniaxialMaterialWrapper^ theMat8,
				UniaxialMaterialWrapper^ theMat9, UniaxialMaterialWrapper^ theMat10,
				UniaxialMaterialWrapper^ theMat11, UniaxialMaterialWrapper^ theMat12,
				UniaxialMaterialWrapper^ theMat13);
			BeamColumnJoint2dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
				UniaxialMaterialWrapper^ theMat1, UniaxialMaterialWrapper^ theMat2,
				UniaxialMaterialWrapper^ theMat3, UniaxialMaterialWrapper^ theMat4,
				UniaxialMaterialWrapper^ theMat5, UniaxialMaterialWrapper^ theMat6,
				UniaxialMaterialWrapper^ theMat7, UniaxialMaterialWrapper^ theMat8,
				UniaxialMaterialWrapper^ theMat9, UniaxialMaterialWrapper^ theMat10,
				UniaxialMaterialWrapper^ theMat11, UniaxialMaterialWrapper^ theMat12,
				UniaxialMaterialWrapper^ theMat13, double Hgtfac, double Wdtfac);
			BeamColumnJoint2dWrapper() {};
			~BeamColumnJoint2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class BeamColumnJoint3dWrapper : ElementWrapper
		{

		public:
			BeamColumnJoint3dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
				UniaxialMaterialWrapper^ theMat1, UniaxialMaterialWrapper^ theMat2,
				UniaxialMaterialWrapper^ theMat3, UniaxialMaterialWrapper^ theMat4,
				UniaxialMaterialWrapper^ theMat5, UniaxialMaterialWrapper^ theMat6,
				UniaxialMaterialWrapper^ theMat7, UniaxialMaterialWrapper^ theMat8,
				UniaxialMaterialWrapper^ theMat9, UniaxialMaterialWrapper^ theMat10,
				UniaxialMaterialWrapper^ theMat11, UniaxialMaterialWrapper^ theMat12,
				UniaxialMaterialWrapper^ theMat13);
			BeamColumnJoint3dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
				UniaxialMaterialWrapper^ theMat1, UniaxialMaterialWrapper^ theMat2,
				UniaxialMaterialWrapper^ theMat3, UniaxialMaterialWrapper^ theMat4,
				UniaxialMaterialWrapper^ theMat5, UniaxialMaterialWrapper^ theMat6,
				UniaxialMaterialWrapper^ theMat7, UniaxialMaterialWrapper^ theMat8,
				UniaxialMaterialWrapper^ theMat9, UniaxialMaterialWrapper^ theMat10,
				UniaxialMaterialWrapper^ theMat11, UniaxialMaterialWrapper^ theMat12,
				UniaxialMaterialWrapper^ theMat13, double Hgtfac, double Wdtfac);
			BeamColumnJoint3dWrapper() {};
			~BeamColumnJoint3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ElasticTubularJointWrapper : ElementWrapper
		{

		public:
			ElasticTubularJointWrapper(int tag, int iNode, int jNode,
				double Brace_Diameter,
				double Brace_Angle,
				double e,
				double Chord_Diameter,
				double Chord_Thickness,
				double Chord_Angle);
			ElasticTubularJointWrapper() {};
			~ElasticTubularJointWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class Joint2DWrapper : ElementWrapper
		{

		public:
			Joint2DWrapper(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
				UniaxialMaterialWrapper^ spring1,
				UniaxialMaterialWrapper^ spring2,
				UniaxialMaterialWrapper^ spring3,
				UniaxialMaterialWrapper^ spring4,
				UniaxialMaterialWrapper^ springC,
				DomainWrapper^ theDomain,
				int LrgDisp);

			Joint2DWrapper(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
				UniaxialMaterialWrapper^ spring1,
				UniaxialMaterialWrapper^ spring2,
				UniaxialMaterialWrapper^ spring3,
				UniaxialMaterialWrapper^ spring4,
				UniaxialMaterialWrapper^ springC,
				DomainWrapper^ theDomain,
				int LrgDisp,
				DamageModelWrapper^ dmg1,
				DamageModelWrapper^ dmg2,
				DamageModelWrapper^ dmg3,
				DamageModelWrapper^ dmg4,
				DamageModelWrapper^ dmgC);

			Joint2DWrapper() {};
			~Joint2DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class Joint3DWrapper : ElementWrapper
		{

		public:
			Joint3DWrapper(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int IntNodeTag,
				UniaxialMaterialWrapper^ springx,
				UniaxialMaterialWrapper^ springy,
				UniaxialMaterialWrapper^ springz,
				DomainWrapper^ theDomain,
				int LrgDisp);

			Joint3DWrapper() {};
			~Joint3DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

#pragma endregion

#pragma region pyMacro
		public ref class PY_Macro2DWrapper : ElementWrapper
		{

		public:
			PY_Macro2DWrapper(int tag,
				int Nd1,
				int Nd2,
				double K,
				double py,
				double a,
				double b,
				double g,
				double m1,
				double m2,
				double w1,
				double p1,
				double S1,
				double beta,
				double s1,
				double tolerance,
				int maxNumIter);

			PY_Macro2DWrapper() {};
			~PY_Macro2DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};
#pragma endregion


#pragma region tetrahedron
		public ref class FourNodeTetrahedronWrapper : ElementWrapper
		{

		public:
			FourNodeTetrahedronWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				NDMaterialWrapper^ theMaterial,
				double b1, double b2, double b3);

			FourNodeTetrahedronWrapper() {};
			~FourNodeTetrahedronWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};
#pragma endregion

#pragma region Tri31
		

		public ref class Tri31Wrapper : ElementWrapper
		{

		public:
			Tri31Wrapper(int tag, int nd1, int nd2, int nd3,
				NDMaterialWrapper^ m, PlaneElementType^ type,
				double t);

			Tri31Wrapper(int tag, int nd1, int nd2, int nd3,
				NDMaterialWrapper^ m, PlaneElementType^ type,
				double t, double pressure,
				double rho,
				double b1, double b2);

			Tri31Wrapper() {};
			~Tri31Wrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};
#pragma endregion

#pragma region twoNodeLink
		public ref class TwoNodeLinkWrapper : ElementWrapper
		{

		public:
			TwoNodeLinkWrapper(int tag, int dimension, int Nd1, int Nd2,
				IDWrapper^ direction, array<UniaxialMaterialWrapper^>^ theMaterials);

			TwoNodeLinkWrapper(int tag, int dimension, int Nd1, int Nd2,
				IDWrapper^ direction, array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^ y, VectorWrapper^ x,
				VectorWrapper^ Mratio, VectorWrapper^ shearDistI,
				int addRayleigh, double mass);

			TwoNodeLinkWrapper() {};
			~TwoNodeLinkWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};
#pragma endregion

#pragma region up
		public ref class BBarBrickUPWrapper : ElementWrapper
		{

		public:
			BBarBrickUPWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				int node5,
				int node6,
				int node7,
				int node8,
				NDMaterialWrapper^ theMaterial, double bulk, double rhof,
				double perm1, double perm2, double perm3,
				double b1, double b2, double b3);

			BBarBrickUPWrapper() {};
			~BBarBrickUPWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class BBarFourNodeQuadUPWrapper : ElementWrapper
		{

		public:
			BBarFourNodeQuadUPWrapper(int tag, int nd1, int nd2, int nd3, int nd4,
				NDMaterialWrapper^ m, PlaneElementType^ type,
				double t, double bulk, double rhof, double perm1, double perm2,
				double b1, double b2, double p);

			BBarFourNodeQuadUPWrapper() {};
			~BBarFourNodeQuadUPWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class BrickUPWrapper : ElementWrapper
		{

		public:
			BrickUPWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				int node5,
				int node6,
				int node7,
				int node8,
				NDMaterialWrapper^ theMaterial, double bulk, double rhof,
				double perm1, double perm2, double perm3,
				double b1, double b2, double b3);

			BrickUPWrapper() {};
			~BrickUPWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class FourNodeQuadUPWrapper : ElementWrapper
		{

		public:
			FourNodeQuadUPWrapper(int tag, int nd1, int nd2, int nd3, int nd4,
				NDMaterialWrapper^ m, PlaneElementType^ type,
				double t, double bulk, double rhof, double perm1, double perm2,
				double b1, double b2, double p);

			FourNodeQuadUPWrapper() {};
			~FourNodeQuadUPWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class NineFourNodeQuadUPWrapper : ElementWrapper
		{

		public:
			NineFourNodeQuadUPWrapper(int tag, int nd1, int nd2, int nd3, int nd4,
				int nd5, int nd6, int nd7, int nd8, int nd9,
				NDMaterialWrapper^ m, PlaneElementType^ type,
				double t, double bulk, double rhof, double perm1, double perm2,
				double b1, double b2);

			NineFourNodeQuadUPWrapper() {};
			~NineFourNodeQuadUPWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

#pragma endregion

#pragma region UWelements
		public ref class BeamContact2DWrapper : ElementWrapper
		{

		public:
			BeamContact2DWrapper(int tag, int Nd1, int Nd2, int NdS, int NdL, NDMaterialWrapper^ theMat,
				double width, double tolG, double tolF, int cSwitch);

			BeamContact2DWrapper() {};
			~BeamContact2DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class BeamContact2DpWrapper : ElementWrapper
		{

		public:
			BeamContact2DpWrapper(int tag, int Nd1, int Nd2, int NdS, NDMaterialWrapper^ theMat,
				double width, double pen, int cSwitch);

			BeamContact2DpWrapper() {};
			~BeamContact2DpWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class BeamContact3DWrapper : ElementWrapper
		{

		public:
			BeamContact3DWrapper(int tag, int Nd1, int Nd2,
				int NdS, int NdL, double rad, CrdTransfWrapper^ coordTransf,
				NDMaterialWrapper^ theMat, double tolG, double tolF, int cSwitch);

			BeamContact3DWrapper() {};
			~BeamContact3DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class BeamContact3DpWrapper : ElementWrapper
		{

		public:
			BeamContact3DpWrapper(int tag, int Nd1, int Nd2, int NdS, double rad, CrdTransfWrapper^ coordTransf,
				NDMaterialWrapper^ theMat, double pen, int cSwitch);

			BeamContact3DpWrapper() {};
			~BeamContact3DpWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class BeamEndContact3DWrapper : ElementWrapper
		{

		public:
			BeamEndContact3DWrapper(int tag, int Nd1, int Nd2, int NdS, int NdL,
				double rad, double tolG, double tolF, int cSwitch);

			BeamEndContact3DWrapper() {};
			~BeamEndContact3DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class BeamEndContact3DpWrapper : ElementWrapper
		{

		public:
			BeamEndContact3DpWrapper(int tag, int Nd1, int Nd2, int NdS, double rad, double pen, int cSwitch);

			BeamEndContact3DpWrapper() {};
			~BeamEndContact3DpWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class Brick8FiberOverlayWrapper : ElementWrapper
		{

		public:
			Brick8FiberOverlayWrapper(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int nd7, int nd8,
				UniaxialMaterialWrapper^ m, double Af, double beta1, double beta2, double beta3, double beta4);

			Brick8FiberOverlayWrapper() {};
			~Brick8FiberOverlayWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class PileToe3DWrapper : ElementWrapper
		{

		public:
			PileToe3DWrapper(int tag, int Nd1, int BNd1, int BNd2, double rad, double k, CrdTransfWrapper^ coordTransf);

			PileToe3DWrapper() {};
			~PileToe3DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		

		public ref class QuadBeamEmbedContactWrapper : ElementWrapper
		{

		public:
			QuadBeamEmbedContactWrapper(int tag, int Qnd1, int Qnd2, int Qnd3, int Qnd4,
				int Bnd1, int Bnd2, double r, double mu, double ep, double et);

			QuadBeamEmbedContactWrapper() {};
			~QuadBeamEmbedContactWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class SimpleContact2DWrapper : ElementWrapper
		{

		public:
			SimpleContact2DWrapper(int tag, int Nd1, int Nd2, int NdS, int NdL,
				NDMaterialWrapper^ theMat, double tolG, double tolF);

			SimpleContact2DWrapper() {};
			~SimpleContact2DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class SimpleContact3DWrapper : ElementWrapper
		{

		public:
			SimpleContact3DWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
				int NdS, int NdL, NDMaterialWrapper^ theMat, double tolG, double tolF);

			SimpleContact3DWrapper() {};
			~SimpleContact3DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class SSPbrickWrapper : ElementWrapper
		{

		public:
			SSPbrickWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
				NDMaterialWrapper^ theMat, double b1, double b2, double b3);

			SSPbrickWrapper() {};
			~SSPbrickWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class SSPbrickUPWrapper : ElementWrapper
		{

		public:
			SSPbrickUPWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
				NDMaterialWrapper^ theMat, double Kf, double Rf, double k1, double k2, double k3,
				double eVoid, double alpha, double b1, double b2, double b3);

			SSPbrickUPWrapper() {};
			~SSPbrickUPWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class SSPquadWrapper : ElementWrapper
		{

		public:
			SSPquadWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterialWrapper^ theMat,
				PlaneElementType^ type, double thick, double b1, double b2);

			SSPquadWrapper() {};
			~SSPquadWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class SSPquadUPWrapper : ElementWrapper
		{

		public:
			SSPquadUPWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterialWrapper^ theMat,
				double thick, double Kf, double Rf, double k1, double k2,
				double eVoid, double alpha, double b1, double b2,
				double Pup, double Plow, double Pleft, double Pright);

			SSPquadUPWrapper() {};
			~SSPquadUPWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};


#pragma endregion

#pragma region XMUelements
		public ref class AC3D8HexWithSensitivityWrapper : ElementWrapper
		{

		public:
			AC3D8HexWithSensitivityWrapper(int element_number,
				int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
				int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
				NDMaterialWrapper^ Globalmmodel);

			AC3D8HexWithSensitivityWrapper(int element_number,
				int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
				int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8);

			AC3D8HexWithSensitivityWrapper() {};
			~AC3D8HexWithSensitivityWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ASI3D8QuadWithSensitivityWrapper : ElementWrapper
		{

		public:
			ASI3D8QuadWithSensitivityWrapper(int element_number,
				int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
				int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8);

			ASI3D8QuadWithSensitivityWrapper() {};
			~ASI3D8QuadWithSensitivityWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class AV3D4QuadWithSensitivityWrapper : ElementWrapper
		{

		public:
			AV3D4QuadWithSensitivityWrapper(int element_number,
				int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4);

			AV3D4QuadWithSensitivityWrapper(int element_number,
				int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
				NDMaterialWrapper^ Globalmmodel);

			AV3D4QuadWithSensitivityWrapper() {};
			~AV3D4QuadWithSensitivityWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class VS3D4QuadWithSensitivityWrapper : ElementWrapper
		{

		public:
			VS3D4QuadWithSensitivityWrapper(int element_number,
				int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
				double _E, double _G, double _rho, double _R, double _alphaN,
				double _alphaT);

			VS3D4QuadWithSensitivityWrapper() {};
			~VS3D4QuadWithSensitivityWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};
#pragma endregion
	}

}
