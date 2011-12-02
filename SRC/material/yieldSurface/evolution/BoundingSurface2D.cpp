//BoundingSurface2D
//
//////////////////////////////////////////////////////////////////////

#include "BoundingSurface2D.h"
#include <math.h>
#include <YieldSurface_BC.h>

#define evolDebug 0
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BoundingSurface2D::BoundingSurface2D(int tag, int classTag, double min_iso_factor,
						double iso_ratio, double kin_ratio,
                        PlasticHardeningMaterial &kpx,
						PlasticHardeningMaterial &kpy,
						YieldSurface_BC &bound_surface)
:YS_Evolution2D(tag, classTag, min_iso_factor, iso_ratio, kin_ratio)
{
	kpMatX = kpx.getCopy();
	kpMatY = kpy.getCopy();
	boundSurface = bound_surface.getCopy();
}

BoundingSurface2D::~BoundingSurface2D()
{
  if (kpMatX != 0)
    delete kpMatX;

  if (kpMatY != 0)
    delete kpMatY;

  if (boundSurface != 0)
    delete boundSurface;
}

int BoundingSurface2D::commitState()
{
	this->YS_Evolution2D::commitState();

    int res  = kpMatX->commitState();
	    res += kpMatY->commitState();

	return res;
}

int BoundingSurface2D::revertToLastCommit(void)
{
	this->YS_Evolution2D::revertToLastCommit();

	kpMatX->revertToLastCommit();
	kpMatY->revertToLastCommit();

	return 0;
}

const Vector &BoundingSurface2D::getEquiPlasticStiffness(void)
{
// find the conjugate point
// find "dR"
// find Kp
	  {
		v2(0) = 0;
		v2(1) = 0;
	  }
	return v2;
}

void    BoundingSurface2D::setTrialPlasticStrains(double ep, const Vector &f, const Vector &g)
{
// set for isotropic materials	
	
}

double  BoundingSurface2D::getIsoPlasticStiffness(int dir)
{
	return 0;
}

double  BoundingSurface2D::getKinPlasticStiffness(int dir)
{
	return 0;	
}

Vector& BoundingSurface2D::getEvolDirection(Vector &f_new)
{
// find the conjugate point
// return that as dir	
	return v2;
}


void BoundingSurface2D::Print(OPS_Stream &s, int flag)
{
	s << "BoundingSurface2D \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
	
