
#ifndef NULLPLASTICMATERIAL_H
#define NULLPLASTICMATERIAL_H

#include "PlasticHardeningMaterial.h"

/**
  *@author rkaul
  */

class NullPlasticMaterial : public PlasticHardeningMaterial  {
public: 
	NullPlasticMaterial(int tag);
	NullPlasticMaterial();
	
	~NullPlasticMaterial();

  double getTrialPlasticStiffness();
  PlasticHardeningMaterial * getCopy();

};
 
#endif
