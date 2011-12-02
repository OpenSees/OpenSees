// CombinedIsoKin2D01.h
//
//////////////////////////////////////////////////////////////////////

#if !defined USERDEFINED2D_H
#define USERDEFINED2D_H

#include "PlasticHardening2D.h"

class  CombinedIsoKin2D01: public PlasticHardening2D
{
public:
    CombinedIsoKin2D01(int tag,
				  double iso_ratio, double kin_ratio,
                  double shr_iso_ratio, double shr_kin_ratio,
				  double min_iso_factor,
                PlasticHardeningMaterial &kpx_pos,
                PlasticHardeningMaterial &kpx_neg,                  
				  PlasticHardeningMaterial &kpy_pos,
				  PlasticHardeningMaterial &kpy_neg,
				  bool isDeformable, double dir);
	virtual ~CombinedIsoKin2D01();
	virtual void	Print(OPS_Stream &s, int flag =0);
	virtual YS_Evolution *getCopy(void);

};

#endif
