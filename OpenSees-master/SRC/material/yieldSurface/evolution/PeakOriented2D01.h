//
//
//////////////////////////////////////////////////////////////////////

#if !defined PEAKORIENTED2D01_H
#define PEAKORIENTED2D01_H

#include "PlasticHardening2D.h"

class PeakOriented2D01 : public PlasticHardening2D
{
public:
    PeakOriented2D01(int tag, double min_iso_factor,
                        PlasticHardeningMaterial &kpx,
                        PlasticHardeningMaterial &kpy );
	virtual ~PeakOriented2D01();
	
	virtual void	Print(OPS_Stream &s, int flag =0);
	virtual YS_Evolution *getCopy(void);
	
private:

};

#endif
