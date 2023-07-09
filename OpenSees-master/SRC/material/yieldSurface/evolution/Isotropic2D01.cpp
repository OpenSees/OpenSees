//Isotropic2D01.cpp: implementation of the YS_HardeningModel class.
//
//////////////////////////////////////////////////////////////////////

#include "Isotropic2D01.h"
#include <math.h>

#define evolDebug 0
#define ISOTROPIC2D01_CLASSTAG 1
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Isotropic2D01::Isotropic2D01(int tag, double min_iso_factor,
                        PlasticHardeningMaterial &kpx,
						PlasticHardeningMaterial &kpy)
:PlasticHardening2D(tag, ISOTROPIC2D01_CLASSTAG, min_iso_factor,1, 0,
                    kpx, kpx, kpy, kpy, 0.0)
{

}

Isotropic2D01::~Isotropic2D01()
{

}

YS_Evolution *Isotropic2D01::getCopy(void)
{
	Isotropic2D01 *theCopy = new  Isotropic2D01(this->getTag(), minIsoFactor, *kpMatXPos, *kpMatYPos);
	if(theCopy==0)
	{
		opserr << "WARNING - Isotropic2D01, unable to get copy\n";
	}
	
	return theCopy;
}

void Isotropic2D01::Print(OPS_Stream &s, int flag)
{
	s << "Isotropic2D01 \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
	
