// @ rkaul@stanford.edu
// @ ggd@stanford.edu

#include <YS_Section2D01.h>

YS_Section2D01::YS_Section2D01(void)
:YieldSurfaceSection2d(0, SEC_TAG_YieldSurface2D01, 0, true),
 E(0), A(0), I(0)
{
}

YS_Section2D01::YS_Section2D01
(int tag, double E_in, double A_in, double I_in, YieldSurface_BC *ptrys, bool use_kr)
:YieldSurfaceSection2d(tag, SEC_TAG_YieldSurface2D01, ptrys, use_kr),
 E(E_in), A(A_in), I(I_in)
{
    if (E <= 0.0)
    {
      opserr << "YS_Section2D01::YS_Section2D01 -- Input E <= 0.0 ... setting E to 1.0\n";
      E = 1.0;
    }

    if (A <= 0.0)
    {
      opserr << "YS_Section2D01::YS_Section2D01 -- Input A <= 0.0 ... setting A to 1.0\n";
      A = 1.0;
    }

    if (I <= 0.0)
    {
      opserr << "YS_Section2D01::YS_Section2D01 -- Input I <= 0.0 ... setting I to 1.0",
	I = 1.0;
    }

}


YS_Section2D01::~YS_Section2D01(void)
{

}

void YS_Section2D01::getSectionStiffness(Matrix &Ks)
{
	Ks(0,0) = E*A; Ks(0,1) = 0.0;
	Ks(1,0) = 0.0; Ks(1,1) = E*I;
}

const Matrix &
YS_Section2D01::getInitialTangent(void)
{
	ks(0,0) = E*A; ks(0,1) = 0.0;
	ks(1,0) = 0.0; ks(1,1) = E*I;
	
	return ks;
}

SectionForceDeformation* YS_Section2D01::getCopy ()
{
    // Make a copy of the hinge
    YS_Section2D01 *theCopy =
    new YS_Section2D01 (this->getTag(), E, A, I, ys, use_Kr_orig);
    theCopy->eCommit = eCommit;
    theCopy->sCommit = sCommit;

    return theCopy;
}

void YS_Section2D01::Print (OPS_Stream &s, int flag)
{
    s << "YS_Section2D01, tag: " << this->getTag() << endln;
	s << "\tE: " << E << endln;
	s << "\tA: " << A << endln;
	s << "\tI: " << I << endln;
	this->YieldSurfaceSection2d::Print(s, flag);
}
