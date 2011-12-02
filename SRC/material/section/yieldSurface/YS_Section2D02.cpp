// @ rkaul@stanford.edu
//   ggd@stanford.edu

#include <math.h>
#include <YS_Section2D02.h>


YS_Section2D02::YS_Section2D02(void)
:YieldSurfaceSection2d(0, SEC_TAG_YieldSurface2D02, 0, true),
 E(0), A(0), I(0), maxPlstkRot(0.0), peakPlstkRot(0.0), iFactor(1.0)
{
}

YS_Section2D02::YS_Section2D02
(int tag, double E_in, double A_in, double I_in, double theta_p_max, YieldSurface_BC *ptrys, bool use_kr)
:YieldSurfaceSection2d(tag, SEC_TAG_YieldSurface2D02, ptrys, use_kr),
 E(E_in), A(A_in), I(I_in), maxPlstkRot(theta_p_max),
 peakPlstkRot(0.0), iFactor(1.0)
{
    if (E <= 0.0)
    {
      opserr << "YS_Section2D02::YS_Section2D02s -- Input E <= 0.0 ... setting E to 1.0\n";
      E = 1.0;
    }

    if (A <= 0.0)
    {
      opserr << "YS_Section2D02::YS_Section2D02s -- Input A <= 0.0 ... setting A to 1.0\n";
      A = 1.0;
    }

    if (I <= 0.0)
    {
      opserr << "YS_Section2D02::YS_Section2D02s -- Input I <= 0.0 ... setting I to 1.0\n";
      I = 1.0;
    }

    if(maxPlstkRot <= 0.0)
    {
      opserr << "YS_Section2D02::YS_Section2D02s -- Input maxPlstkRot <= 0.0 ... setting to 0.0\n";
      maxPlstkRot = 0.0;
    }

}


YS_Section2D02::~YS_Section2D02(void)
{

}

int YS_Section2D02::commitState(void)
{
	double pRot = fabs(ys->hModel->getTrialPlasticStrains(0));
	if(pRot > peakPlstkRot)
		peakPlstkRot = pRot;

	if(fabs(maxPlstkRot) <= 1e-10)
		iFactor = 1.0;
	else
		iFactor =  1 - (peakPlstkRot/maxPlstkRot);

	if(iFactor < 0.02)
		iFactor = 0.02;

	opserr << peakPlstkRot << "\t" << iFactor << endln;
	return this->YieldSurfaceSection2d::commitState();
}

void YS_Section2D02::getSectionStiffness(Matrix &Ks)
{
//	cout << "iFactor: " << iFactor << "\tpeak rotP: " <<  peakPlstkRot
//	     << "\tmaxP: " << maxPlstkRot << endln;

	Ks(0,0) = E*A; Ks(0,1) = 0.0;
	Ks(1,0) = 0.0; Ks(1,1) = E*I*iFactor;
}

const Matrix &
YS_Section2D02::getInitialTangent(void)
{
	ks(0,0) = E*A; ks(0,1) = 0.0;
	ks(1,0) = 0.0; ks(1,1) = E*I;
	
	return ks;
}


SectionForceDeformation* YS_Section2D02::getCopy ()
{
    // Make a copy of the hinge
    YS_Section2D02 *theCopy =
    new YS_Section2D02 (this->getTag(), E, A, I, maxPlstkRot, ys, use_Kr_orig);
    theCopy->eCommit = eCommit;
    theCopy->sCommit = sCommit;
    theCopy->peakPlstkRot =  peakPlstkRot;

    return theCopy;
}

void YS_Section2D02::Print (OPS_Stream &s, int flag)
{
    s << "YS_Section2D02, tag: " << this->getTag() << endln;
	s << "\tE: " << E << endln;
	s << "\tA: " << A << endln;
	s << "\tI: " << I << endln;
	s << "\tDegradation Factor      :" << iFactor << endln;
	s << "\tPeak plastic-rotation   :" << peakPlstkRot << endln;
	s << "\tMaximum plastic-rotation:" << maxPlstkRot << endln;
	this->YieldSurfaceSection2d::Print(s, flag);
}
