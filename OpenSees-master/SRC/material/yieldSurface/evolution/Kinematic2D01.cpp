//Kinematic2D01.cpp: implementation of the YS_HardeningModel class.
//
//////////////////////////////////////////////////////////////////////

#include "Kinematic2D01.h"
#include <math.h>

#define evolDebug 0
#define KINEMATIC2D01_CLASSTAG -1
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Kinematic2D01::Kinematic2D01(int tag, double min_iso_factor,
                        PlasticHardeningMaterial &kpx,
						PlasticHardeningMaterial &kpy, double dir)
:PlasticHardening2D(tag, KINEMATIC2D01_CLASSTAG, min_iso_factor,0, 1,
                    kpx, kpx, kpy, kpy, dir)
{

}

Kinematic2D01::~Kinematic2D01()
{

}

YS_Evolution *Kinematic2D01::getCopy(void)
{
	Kinematic2D01 *theCopy = new  Kinematic2D01(this->getTag(), minIsoFactor, *kpMatXPos, *kpMatYPos, direction);
	if(theCopy==0)
	{
		opserr << "WARNING - Kinematic2D01, unable to get copy\n";
	}
	
	return theCopy;
}

void Kinematic2D01::Print(OPS_Stream &s, int flag)
{
	s << "Kinematic2D01 \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
	
