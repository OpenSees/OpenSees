#include "stdafx.h"
#include "OtherElementWrapper.h"


using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;

GenericClientWrapper::GenericClientWrapper(int tag, IDWrapper^ nodes, IDWrapper^ dof,
	int port) {
	_Element = new GenericClient(tag, *nodes->_ID, dof->_ID, port);
}

GenericClientWrapper::GenericClientWrapper(int tag, IDWrapper^ nodes, IDWrapper^ dof,
	int port, String^ machineInetAddr,
	int ssl, int udp, int dataSize,
	int addRayleigh) {
	_Element = new GenericClient(tag, *nodes->_ID, dof->_ID,
		port, (char*)(void*)Marshal::StringToHGlobalAnsi(machineInetAddr),
		ssl, udp, dataSize,
		addRayleigh);
}

GenericCopyWrapper::GenericCopyWrapper(int tag, IDWrapper^ nodes, int srcTag) {
	_Element = new GenericCopy( tag, *nodes->_ID,  srcTag);
}

BeamColumnJoint2dWrapper::BeamColumnJoint2dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
	UniaxialMaterialWrapper^ theMat1, UniaxialMaterialWrapper^ theMat2,
	UniaxialMaterialWrapper^ theMat3, UniaxialMaterialWrapper^ theMat4,
	UniaxialMaterialWrapper^ theMat5, UniaxialMaterialWrapper^ theMat6,
	UniaxialMaterialWrapper^ theMat7, UniaxialMaterialWrapper^ theMat8,
	UniaxialMaterialWrapper^ theMat9, UniaxialMaterialWrapper^ theMat10,
	UniaxialMaterialWrapper^ theMat11, UniaxialMaterialWrapper^ theMat12,
	UniaxialMaterialWrapper^ theMat13) {
	
	_Element = new BeamColumnJoint2d(tag, Nd1, Nd2, Nd3, Nd4,
		*theMat1->_UniaxialMaterial, *theMat2->_UniaxialMaterial,
		*theMat3->_UniaxialMaterial, *theMat4->_UniaxialMaterial,
		*theMat5->_UniaxialMaterial, *theMat6->_UniaxialMaterial,
		*theMat7->_UniaxialMaterial, *theMat8->_UniaxialMaterial,
		*theMat9->_UniaxialMaterial, *theMat10->_UniaxialMaterial,
		*theMat11->_UniaxialMaterial, *theMat12->_UniaxialMaterial,
		*theMat13->_UniaxialMaterial);
}

BeamColumnJoint2dWrapper::BeamColumnJoint2dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
	UniaxialMaterialWrapper^ theMat1, UniaxialMaterialWrapper^ theMat2,
	UniaxialMaterialWrapper^ theMat3, UniaxialMaterialWrapper^ theMat4,
	UniaxialMaterialWrapper^ theMat5, UniaxialMaterialWrapper^ theMat6,
	UniaxialMaterialWrapper^ theMat7, UniaxialMaterialWrapper^ theMat8,
	UniaxialMaterialWrapper^ theMat9, UniaxialMaterialWrapper^ theMat10,
	UniaxialMaterialWrapper^ theMat11, UniaxialMaterialWrapper^ theMat12,
	UniaxialMaterialWrapper^ theMat13, double Hgtfac, double Wdtfac) {

	_Element = new BeamColumnJoint2d(tag, Nd1, Nd2, Nd3, Nd4,
		*theMat1->_UniaxialMaterial, *theMat2->_UniaxialMaterial,
		*theMat3->_UniaxialMaterial, *theMat4->_UniaxialMaterial,
		*theMat5->_UniaxialMaterial, *theMat6->_UniaxialMaterial,
		*theMat7->_UniaxialMaterial, *theMat8->_UniaxialMaterial,
		*theMat9->_UniaxialMaterial, *theMat10->_UniaxialMaterial,
		*theMat11->_UniaxialMaterial, *theMat12->_UniaxialMaterial,
		*theMat13->_UniaxialMaterial, Hgtfac, Wdtfac);
}

BeamColumnJoint3dWrapper::BeamColumnJoint3dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
	UniaxialMaterialWrapper^ theMat1, UniaxialMaterialWrapper^ theMat2,
	UniaxialMaterialWrapper^ theMat3, UniaxialMaterialWrapper^ theMat4,
	UniaxialMaterialWrapper^ theMat5, UniaxialMaterialWrapper^ theMat6,
	UniaxialMaterialWrapper^ theMat7, UniaxialMaterialWrapper^ theMat8,
	UniaxialMaterialWrapper^ theMat9, UniaxialMaterialWrapper^ theMat10,
	UniaxialMaterialWrapper^ theMat11, UniaxialMaterialWrapper^ theMat12,
	UniaxialMaterialWrapper^ theMat13) {

	_Element = new BeamColumnJoint3d(tag, Nd1, Nd2, Nd3, Nd4,
		*theMat1->_UniaxialMaterial, *theMat2->_UniaxialMaterial,
		*theMat3->_UniaxialMaterial, *theMat4->_UniaxialMaterial,
		*theMat5->_UniaxialMaterial, *theMat6->_UniaxialMaterial,
		*theMat7->_UniaxialMaterial, *theMat8->_UniaxialMaterial,
		*theMat9->_UniaxialMaterial, *theMat10->_UniaxialMaterial,
		*theMat11->_UniaxialMaterial, *theMat12->_UniaxialMaterial,
		*theMat13->_UniaxialMaterial);
}

BeamColumnJoint3dWrapper::BeamColumnJoint3dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
	UniaxialMaterialWrapper^ theMat1, UniaxialMaterialWrapper^ theMat2,
	UniaxialMaterialWrapper^ theMat3, UniaxialMaterialWrapper^ theMat4,
	UniaxialMaterialWrapper^ theMat5, UniaxialMaterialWrapper^ theMat6,
	UniaxialMaterialWrapper^ theMat7, UniaxialMaterialWrapper^ theMat8,
	UniaxialMaterialWrapper^ theMat9, UniaxialMaterialWrapper^ theMat10,
	UniaxialMaterialWrapper^ theMat11, UniaxialMaterialWrapper^ theMat12,
	UniaxialMaterialWrapper^ theMat13, double Hgtfac, double Wdtfac) {

	_Element = new BeamColumnJoint3d(tag, Nd1, Nd2, Nd3, Nd4,
		*theMat1->_UniaxialMaterial, *theMat2->_UniaxialMaterial,
		*theMat3->_UniaxialMaterial, *theMat4->_UniaxialMaterial,
		*theMat5->_UniaxialMaterial, *theMat6->_UniaxialMaterial,
		*theMat7->_UniaxialMaterial, *theMat8->_UniaxialMaterial,
		*theMat9->_UniaxialMaterial, *theMat10->_UniaxialMaterial,
		*theMat11->_UniaxialMaterial, *theMat12->_UniaxialMaterial,
		*theMat13->_UniaxialMaterial, Hgtfac, Wdtfac);
}

ElasticTubularJointWrapper::ElasticTubularJointWrapper(int tag, int iNode, int jNode,
	double Brace_Diameter,
	double Brace_Angle,
	double e,
	double Chord_Diameter,
	double Chord_Thickness,
	double Chord_Angle) {

	_Element = new ElasticTubularJoint(tag, iNode, jNode,
		Brace_Diameter,
		Brace_Angle,
		e,
		Chord_Diameter,
		Chord_Thickness,
		Chord_Angle);
}

Joint2DWrapper::Joint2DWrapper(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
	UniaxialMaterialWrapper^ spring1,
	UniaxialMaterialWrapper^ spring2,
	UniaxialMaterialWrapper^ spring3,
	UniaxialMaterialWrapper^ spring4,
	UniaxialMaterialWrapper^ springC,
	DomainWrapper^ theDomain,
	int LrgDisp) {

	_Element = new Joint2D(tag, nd1, nd2, nd3, nd4, IntNodeTag,
		*spring1->_UniaxialMaterial,
		*spring2->_UniaxialMaterial,
		*spring3->_UniaxialMaterial,
		*spring4->_UniaxialMaterial,
		*springC->_UniaxialMaterial,
		theDomain->_Domain,
		LrgDisp);
}

MP_Joint2DWrapper::MP_Joint2DWrapper(DomainWrapper^ theDomain, int nodeRetain, int nodeConstr,
	int Maindof, int fixedend, int LrgDsp) {
	_MP_Constraint = new MP_Joint2D(theDomain->_Domain, nodeRetain, nodeConstr,
		Maindof, fixedend, LrgDsp);
}


MP_Joint3DWrapper::MP_Joint3DWrapper(DomainWrapper^ theDomain, int nodeRetain, int nodeConstr,
	int nodeRot, int Rotdof, int nodeDisp, int Dispdof, int LrgDsp) {
	_MP_Constraint = new MP_Joint3D(theDomain->_Domain, nodeRetain, nodeConstr,
		nodeRot, Rotdof, nodeDisp, Dispdof, LrgDsp);
}

Joint2DWrapper::Joint2DWrapper(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
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
	DamageModelWrapper^ dmgC) {

	_Element = new Joint2D(tag, nd1, nd2, nd3, nd4, IntNodeTag,
		*spring1->_UniaxialMaterial,
		*spring2->_UniaxialMaterial,
		*spring3->_UniaxialMaterial,
		*spring4->_UniaxialMaterial,
		*springC->_UniaxialMaterial,
		theDomain->_Domain,
		LrgDisp,
		*dmg1->_DamageModel,
		*dmg2->_DamageModel,
		*dmg3->_DamageModel,
		*dmg4->_DamageModel,
		*dmgC->_DamageModel);
}

Joint3DWrapper::Joint3DWrapper(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int IntNodeTag,
	UniaxialMaterialWrapper^ springx,
	UniaxialMaterialWrapper^ springy,
	UniaxialMaterialWrapper^ springz,
	DomainWrapper^ theDomain,
	int LrgDisp) {

	_Element = new Joint3D( tag,  nd1,  nd2,  nd3,  nd4,  nd5,  nd6,  IntNodeTag,
		*springx->_UniaxialMaterial,
		*springy->_UniaxialMaterial,
		*springz->_UniaxialMaterial,
		theDomain->_Domain,
		LrgDisp);
}

PY_Macro2DWrapper::PY_Macro2DWrapper(int tag,
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
	int maxNumIter) {

	_Element = new PY_Macro2D( tag,
		 Nd1,
		 Nd2,
		 K,
		 py,
		 a,
		 b,
		 g,
		 m1,
		 m2,
		 w1,
		 p1,
		 S1,
		 beta,
		 s1,
		 tolerance,
		 maxNumIter);
}

FourNodeTetrahedronWrapper::FourNodeTetrahedronWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	NDMaterialWrapper^ theMaterial,
	double b1, double b2, double b3) {

	_Element = new FourNodeTetrahedron(tag,
		node1,
		node2,
		node3,
		node4,
		*theMaterial->_NDMaterial,
		b1, b2, b3);
}

Tri31Wrapper::Tri31Wrapper(int tag, int nd1, int nd2, int nd3,
	NDMaterialWrapper^ m, PlaneElementType^ type,
	double t) {

	_Element = new Tri31(tag, nd1, nd2, nd3,
		*m->_NDMaterial, (char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()),
		t);
}

Tri31Wrapper::Tri31Wrapper(int tag, int nd1, int nd2, int nd3,
	NDMaterialWrapper^ m, PlaneElementType^ type, double t, double pressure,
	double rho,
	double b1, double b2) {

	_Element = new Tri31(tag, nd1, nd2, nd3,
		*m->_NDMaterial, (char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()),
		t, pressure,
		rho,
		b1, b2);
}

UniaxialMaterial** array2pointer2(array<UniaxialMaterialWrapper^>^ theMaterials) {
	UniaxialMaterial** mats = new UniaxialMaterial*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_UniaxialMaterial;
	}

	return mats;
}

TwoNodeLinkWrapper::TwoNodeLinkWrapper(int tag, int dimension, int Nd1, int Nd2,
	IDWrapper^ direction, array<UniaxialMaterialWrapper^>^ theMaterials) {

	_Element = new TwoNodeLink(tag, dimension, Nd1, Nd2,
		*direction->_ID, array2pointer2(theMaterials));
}

TwoNodeLinkWrapper::TwoNodeLinkWrapper(int tag, int dimension, int Nd1, int Nd2,
	IDWrapper^ direction, array<UniaxialMaterialWrapper^>^ theMaterials, VectorWrapper^ y, VectorWrapper^ x,
	VectorWrapper^ Mratio, VectorWrapper^ shearDistI,
	int addRayleigh, double mass) {

	_Element = new TwoNodeLink(tag, dimension, Nd1, Nd2,
		*direction->_ID, array2pointer2(theMaterials), *y->_Vector, *x->_Vector,
		*Mratio->_Vector, *shearDistI->_Vector,
		addRayleigh, mass);
}

BBarBrickUPWrapper::BBarBrickUPWrapper(int tag,
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
	double b1, double b2, double b3) {

	_Element = new BBarBrickUP(tag,
		node1,
		node2,
		node3,
		node4,
		node5,
		node6,
		node7,
		node8,
		*theMaterial->_NDMaterial, bulk, rhof,
		perm1, perm2, perm3,
		b1, b2, b3);
}

BBarFourNodeQuadUPWrapper::BBarFourNodeQuadUPWrapper(int tag, int nd1, int nd2, int nd3, int nd4,
	NDMaterialWrapper^ m, PlaneElementType^ type,
	double t, double bulk, double rhof, double perm1, double perm2,
	double b1, double b2, double p) {

	_Element = new BBarFourNodeQuadUP(tag, nd1, nd2, nd3, nd4,
		*m->_NDMaterial, (char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()),
		t, bulk, rhof, perm1, perm2,
		b1, b2, p);
}

BrickUPWrapper::BrickUPWrapper(int tag,
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
	double b1, double b2, double b3) {

	_Element = new BrickUP(tag,
		node1,
		node2,
		node3,
		node4,
		node5,
		node6,
		node7,
		node8,
		*theMaterial->_NDMaterial, bulk, rhof,
		perm1, perm2, perm3,
		b1, b2, b3);
}

FourNodeQuadUPWrapper::FourNodeQuadUPWrapper(int tag, int nd1, int nd2, int nd3, int nd4,
	NDMaterialWrapper^ m, PlaneElementType^ type,
	double t, double bulk, double rhof, double perm1, double perm2,
	double b1, double b2, double p) {

	_Element = new FourNodeQuadUP(tag, nd1, nd2, nd3, nd4,
		*m->_NDMaterial, (char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()),
		t, bulk, rhof, perm1, perm2,
		b1, b2, p);
}

NineFourNodeQuadUPWrapper::NineFourNodeQuadUPWrapper(int tag, int nd1, int nd2, int nd3, int nd4,
	int nd5, int nd6, int nd7, int nd8, int nd9,
	NDMaterialWrapper^ m, PlaneElementType^ type,
	double t, double bulk, double rhof, double perm1, double perm2,
	double b1, double b2) {

	_Element = new NineFourNodeQuadUP(tag, nd1, nd2, nd3, nd4,
		nd5, nd6, nd7, nd8, nd9,
		*m->_NDMaterial, (char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()),
		t, bulk, rhof, perm1, perm2,
		b1, b2);
}

BeamContact2DWrapper::BeamContact2DWrapper(int tag, int Nd1, int Nd2, int NdS, int NdL, NDMaterialWrapper^ theMat,
	double width, double tolG, double tolF, int cSwitch) {

	_Element = new BeamContact2D(tag, Nd1, Nd2, NdS, NdL, *theMat->_NDMaterial,
		width, tolG, tolF, cSwitch);
}

BeamContact2DpWrapper::BeamContact2DpWrapper(int tag, int Nd1, int Nd2, int NdS, NDMaterialWrapper^ theMat,
	double width, double pen, int cSwitch) {

	_Element = new BeamContact2Dp( tag,  Nd1,  Nd2,  NdS, *theMat->_NDMaterial,
		 width,  pen,  cSwitch);
}

BeamContact3DWrapper::BeamContact3DWrapper(int tag, int Nd1, int Nd2,
	int NdS, int NdL, double rad, CrdTransfWrapper^ coordTransf,
	NDMaterialWrapper^ theMat, double tolG, double tolF, int cSwitch) {

	_Element = new BeamContact3D(tag, Nd1, Nd2,
		NdS, NdL, rad, *coordTransf->_CrdTransf,
		*theMat->_NDMaterial, tolG, tolF, cSwitch);
}

BeamContact3DpWrapper::BeamContact3DpWrapper(int tag, int Nd1, int Nd2,
	int NdS, double rad, CrdTransfWrapper^ coordTransf,
	NDMaterialWrapper^ theMat, double pen, int cSwitch) {

	_Element = new BeamContact3Dp(tag, Nd1, Nd2,
		NdS, rad, *coordTransf->_CrdTransf,
		*theMat->_NDMaterial, pen, cSwitch);
}

BeamEndContact3DWrapper::BeamEndContact3DWrapper(int tag, int Nd1, int Nd2, int NdS, int NdL,
	double rad, double tolG, double tolF, int cSwitch) {

	_Element = new BeamEndContact3D(tag, Nd1, Nd2, NdS, NdL,
		rad, tolG, tolF, cSwitch);
}

BeamEndContact3DpWrapper::BeamEndContact3DpWrapper(int tag, int Nd1, int Nd2, int NdS, double rad, double pen, int cSwitch) {

	_Element = new BeamEndContact3Dp(tag, Nd1, Nd2, NdS, rad, pen, cSwitch);
}

Brick8FiberOverlayWrapper::Brick8FiberOverlayWrapper(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int nd7, int nd8,
	UniaxialMaterialWrapper^ m, double Af, double beta1, double beta2, double beta3, double beta4) {

	_Element = new Brick8FiberOverlay(tag, nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8,
		*m->_UniaxialMaterial, Af, beta1, beta2, beta3, beta4);
}

PileToe3DWrapper::PileToe3DWrapper(int tag, int Nd1, int BNd1, int BNd2, double rad, double k, CrdTransfWrapper^ coordTransf) {

	_Element = new PileToe3D(tag, Nd1, BNd1, BNd2, rad, k,
		*coordTransf->_CrdTransf);
}



QuadBeamEmbedContactWrapper::QuadBeamEmbedContactWrapper(int tag, int Qnd1, int Qnd2, int Qnd3, int Qnd4,
	int Bnd1, int Bnd2, double r, double mu, double ep, double et) {

	_Element = new QuadBeamEmbedContact(tag, Qnd1, Qnd2, Qnd3, Qnd4,
		Bnd1, Bnd2, r, mu, ep, et);
}

SimpleContact2DWrapper::SimpleContact2DWrapper(int tag, int Nd1, int Nd2, int NdS, int NdL,
	NDMaterialWrapper^ theMat, double tolG, double tolF) {

	_Element = new SimpleContact2D(tag, Nd1, Nd2, NdS, NdL,
		*theMat->_NDMaterial, tolG, tolF);
}

SimpleContact3DWrapper::SimpleContact3DWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
	int NdS, int NdL, NDMaterialWrapper^ theMat, double tolG, double tolF) {

	_Element = new SimpleContact3D(tag, Nd1, Nd2, Nd3, Nd4,
		NdS, NdL, *theMat->_NDMaterial, tolG, tolF);
}

SSPbrickWrapper::SSPbrickWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
	NDMaterialWrapper^ theMat, double b1, double b2, double b3) {

	_Element = new SSPbrick(tag, Nd1, Nd2, Nd3, Nd4, Nd5, Nd6, Nd7, Nd8,
		*theMat->_NDMaterial, b1, b2, b3);
}

SSPbrickUPWrapper::SSPbrickUPWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
	NDMaterialWrapper^ theMat, double Kf, double Rf, double k1, double k2, double k3,
	double eVoid, double alpha, double b1, double b2, double b3) {

	_Element = new SSPbrickUP(tag, Nd1, Nd2, Nd3, Nd4, Nd5, Nd6, Nd7, Nd8,
		*theMat->_NDMaterial, Kf, Rf, k1, k2, k3,
		eVoid, alpha, b1, b2, b3);
}

SSPquadWrapper::SSPquadWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterialWrapper^ theMat,
	PlaneElementType^ type, double thick, double b1, double b2) {

	_Element = new SSPquad(tag, Nd1, Nd2, Nd3, Nd4, *theMat->_NDMaterial,
		(char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()), thick, b1, b2);
}

SSPquadUPWrapper::SSPquadUPWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, NDMaterialWrapper^ theMat,
	double thick, double Kf, double Rf, double k1, double k2,
	double eVoid, double alpha, double b1, double b2,
	double Pup, double Plow, double Pleft, double Pright) {

	_Element = new SSPquadUP( tag,  Nd1,  Nd2,  Nd3,  Nd4, *theMat->_NDMaterial,
		 thick,  Kf,  Rf,  k1,  k2,
		 eVoid,  alpha,  b1,  b2,
		 Pup,  Plow,  Pleft,  Pright);
}

AC3D8HexWithSensitivityWrapper::AC3D8HexWithSensitivityWrapper(int element_number,
	int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
	int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
	NDMaterialWrapper^ Globalmmodel) {

	_Element = new AC3D8HexWithSensitivity(element_number,
		node_numb_1, node_numb_2, node_numb_3, node_numb_4,
		node_numb_5, node_numb_6, node_numb_7, node_numb_8,
		Globalmmodel->_NDMaterial);
}

AC3D8HexWithSensitivityWrapper::AC3D8HexWithSensitivityWrapper(int element_number,
	int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
	int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8) {

	_Element = new AC3D8HexWithSensitivity(element_number,
		node_numb_1, node_numb_2, node_numb_3, node_numb_4,
		node_numb_5, node_numb_6, node_numb_7, node_numb_8);
}

ASI3D8QuadWithSensitivityWrapper::ASI3D8QuadWithSensitivityWrapper(int element_number,
	int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
	int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8) {

	_Element = new ASI3D8QuadWithSensitivity(element_number,
		node_numb_1, node_numb_2, node_numb_3, node_numb_4,
		node_numb_5, node_numb_6, node_numb_7, node_numb_8);
}

AV3D4QuadWithSensitivityWrapper::AV3D4QuadWithSensitivityWrapper(int element_number,
	int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4) {

	_Element = new AV3D4QuadWithSensitivity(element_number,
		node_numb_1, node_numb_2, node_numb_3, node_numb_4);
}

AV3D4QuadWithSensitivityWrapper::AV3D4QuadWithSensitivityWrapper(int element_number,
	int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,NDMaterialWrapper^ Globalmmodel) {

	_Element = new AV3D4QuadWithSensitivity(element_number,
		node_numb_1, node_numb_2, node_numb_3, node_numb_4, Globalmmodel->_NDMaterial);
}

VS3D4QuadWithSensitivityWrapper::VS3D4QuadWithSensitivityWrapper(int element_number,
	int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
	double _E, double _G, double _rho, double _R, double _alphaN,
	double _alphaT) {

	_Element = new VS3D4QuadWithSensitivity(element_number,
		node_numb_1, node_numb_2, node_numb_3, node_numb_4,
		_E, _G, _rho, _R, _alphaN,
		_alphaT);
}






