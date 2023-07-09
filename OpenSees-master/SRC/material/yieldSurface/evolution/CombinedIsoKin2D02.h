/***************************************************************************
                          CombinedIsoKin2D02.h  -  description
                             -------------------
    begin                : Fri Jul 12 2002
    email                : rkaul@ce-blume215-pent-2.stanford.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef COMBINEDISOKIN2D02_H
#define COMBINEDISOKIN2D02_H

#include <BkStressLimSurface2D.h>

/**Kinematic model is based on back-stress
Isotropic model on plastic-deformations
This one is also deformable
  *@author rkaul
  */

class CombinedIsoKin2D02 : public BkStressLimSurface2D  {
public: 
	CombinedIsoKin2D02( int tag,  double min_iso_factor,
				double iso_ratio, double kin_ratio,
				YieldSurface_BC  &lim_surface,
				PlasticHardeningMaterial &kinX,
				PlasticHardeningMaterial &kinY,
				PlasticHardeningMaterial &isoXPos,
				PlasticHardeningMaterial &isoXNeg,
				PlasticHardeningMaterial &isoYPos,
				PlasticHardeningMaterial &isoYNeg,
				bool isDeformable,
				int  algo, double resfact, double appfact, double dir
	);
	
	~CombinedIsoKin2D02();
	
  void Print(OPS_Stream & s, int flag);
  YS_Evolution * getCopy();


};

#endif
