/***************************************************************************
                          PeakOriented2D02.h  -  description
                             -------------------
    begin                : Fri Jul 12 2002
    email                : rkaul@ce-blume215-pent-2.stanford.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef PEAKORIENTED2D02_H
#define PEAKORIENTED2D02_H

#include <BkStressLimSurface2D.h>
#include <NullPlasticMaterial.h>

/** Grows/ shrinks one side while keeping the diametrically opposite end 
the same.  Kp_isoXPos = Kp_isoXNeg = (should be) Kp_kinX
Since kinematic is function of back-stress and
isotropic is a function of plastic-strains - calibration
will be difficult.

Therefore, this is a pseudo PeakOriented model, for true
PeakOriented, use PeakOriented01

isotropic:kinematic ratio = 0.5:0.5
  *@author rkaul
  */

class PeakOriented2D02 : public BkStressLimSurface2D  {
public: 
	PeakOriented2D02(int tag, double min_iso_factor,
					YieldSurface_BC &lim_surface,
					PlasticHardeningMaterial &kinX,
					PlasticHardeningMaterial &kinY,
					PlasticHardeningMaterial &isoX,
					PlasticHardeningMaterial &isoY,
					int algo
					);
	~PeakOriented2D02();
	YS_Evolution * getCopy();
	void Print(OPS_Stream & s, int flag);

private: // Private attributes
  /**  */
  static NullPlasticMaterial nullMat;
};

#endif
