// 
//
//////////////////////////////////////////////////////////////////////

#if !defined ISOTROPIC2D01_H
#define ISOTROPIC2D01_H

#include "PlasticHardening2D.h"

class Isotropic2D01 : public PlasticHardening2D
{
public:
    Isotropic2D01(int tag, double min_iso_factor,
				PlasticHardeningMaterial  &kpx,
				PlasticHardeningMaterial  &kpy);

	virtual ~Isotropic2D01();

	virtual void	Print(OPS_Stream &s, int flag =0);
	virtual YS_Evolution *getCopy(void);
	//virtual int 	evolveSurface(YieldSurfaceBC *ys, double magPlasticDefo);
	
private:

};

#endif
