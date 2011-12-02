//CombinedIsoKin2D01.cpp: implementation of the YS_HardeningModel class.
//
//////////////////////////////////////////////////////////////////////

#include "CombinedIsoKin2D01.h"
#include <math.h>

#define evolDebug 0
#define COMBINEDISOKIN2D01_CLASSTAG -1
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CombinedIsoKin2D01::CombinedIsoKin2D01(int tag,
				  double iso_ratio, double kin_ratio,
                  double shr_iso_ratio, double shr_kin_ratio,
				  double min_iso_factor,
                  PlasticHardeningMaterial &kpx_pos,
                  PlasticHardeningMaterial &kpx_neg,
				  PlasticHardeningMaterial &kpy_pos,
				  PlasticHardeningMaterial &kpy_neg,				  
				  bool isDeformable, double dir)
:PlasticHardening2D(tag, COMBINEDISOKIN2D01_CLASSTAG, min_iso_factor,
                    iso_ratio, kin_ratio, kpx_pos, kpx_neg, kpy_pos, kpy_neg, dir)
{
	deformable = isDeformable;
	isotropicRatio_shrink = shr_iso_ratio;
	kinematicRatio_shrink = shr_kin_ratio;
}

CombinedIsoKin2D01::~CombinedIsoKin2D01()
{

}

YS_Evolution *CombinedIsoKin2D01::getCopy(void)
{
	CombinedIsoKin2D01 *theCopy = new  
	CombinedIsoKin2D01(	this->getTag(),
						isotropicRatio, kinematicRatio,
						isotropicRatio_shrink, kinematicRatio_shrink, 
						minIsoFactor,
						*kpMatXPos, *kpMatXNeg, *kpMatYPos, *kpMatYNeg,
						deformable, direction);
	if(theCopy==0)
	{
		opserr << "WARNING - CombinedIsoKin2D01, unable to get copy\n";
	}
	
	return theCopy;
}

void CombinedIsoKin2D01::Print(OPS_Stream &s, int flag)
{
	s << "CombinedIsoKin2D01 \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
	
