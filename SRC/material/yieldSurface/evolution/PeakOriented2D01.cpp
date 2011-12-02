//PeakOriented2D01.cpp: implementation of the YS_HardeningModel class.
//
//////////////////////////////////////////////////////////////////////

#include "PeakOriented2D01.h"
#include <math.h>

#define evolDebug 0
#define PEAK_ORIENTED2D01_CLASSTAG -1
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PeakOriented2D01::PeakOriented2D01(int tag, double min_iso_factor,
                                   PlasticHardeningMaterial &kpx,
                                   PlasticHardeningMaterial &kpy)
:PlasticHardening2D(tag, PEAK_ORIENTED2D01_CLASSTAG, min_iso_factor, 0.5, 0.5,
                    kpx, kpx, kpy, kpy, 0.0)
{

}

PeakOriented2D01::~PeakOriented2D01()
{

}
	
YS_Evolution *PeakOriented2D01::getCopy(void)
{
	PeakOriented2D01 *theCopy = new  PeakOriented2D01(this->getTag(), minIsoFactor, *kpMatXPos, *kpMatYPos);
	if(theCopy==0)
	{
		opserr << "WARNING - PeakOriented2D, unable to get copy\n";
	}
	
	return theCopy;
}

void PeakOriented2D01::Print(OPS_Stream &s, int flag)
{
	s << "PeakOriented2D \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
	
