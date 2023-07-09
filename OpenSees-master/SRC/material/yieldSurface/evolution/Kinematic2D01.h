// Kinematic2D.h
//
//////////////////////////////////////////////////////////////////////

#if !defined KINEMATIC2D01_H
#define KINEMATIC2D01_H

#include "PlasticHardening2D.h"

class Kinematic2D01 : public PlasticHardening2D
{
public:
    Kinematic2D01(int tag, double min_iso_factor,
				PlasticHardeningMaterial &kpx,
				PlasticHardeningMaterial &kpy, double dir);

	virtual ~Kinematic2D01();

	virtual void	Print(OPS_Stream &s, int flag =0);
	virtual YS_Evolution *getCopy(void);
	//virtual int 	evolveSurface(YieldSurfaceBC *ys, double magPlasticDefo);
	
private:
//	double minIsoFactor;
};

#endif
