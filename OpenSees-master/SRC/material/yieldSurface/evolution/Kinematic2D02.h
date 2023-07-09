// Kinematic2D02.h
//
//////////////////////////////////////////////////////////////////////

#if !defined KINEMATIC2D02_H
#define KINEMATIC2D02_H

#include "BkStressLimSurface2D.h"
#include "NullPlasticMaterial.h"

//class YieldSurface_BC;
class Kinematic2D02 : public BkStressLimSurface2D
{
public:
    Kinematic2D02(int tag, double min_iso_factor, 
				YieldSurface_BC &lim_surface,
				PlasticHardeningMaterial &kpx,
				PlasticHardeningMaterial &kpy,
				int algo, double resfact, double appfact, double dir);

	virtual ~Kinematic2D02();
//	virtual int  displaySelf(Renderer &theViewer, int displayMode, float fact);
	virtual void Print(OPS_Stream &s, int flag =0);
	virtual YS_Evolution *getCopy(void);

static NullPlasticMaterial nullMat;
};

#endif
