
#include "NullPlasticMaterial.h"
#define CLASSTAG_NULL_PLASTIC_MAT -1

NullPlasticMaterial::NullPlasticMaterial()
: PlasticHardeningMaterial(-1, CLASSTAG_NULL_PLASTIC_MAT)
{
}

NullPlasticMaterial::NullPlasticMaterial(int tag)
: PlasticHardeningMaterial(tag, CLASSTAG_NULL_PLASTIC_MAT)
{
}

NullPlasticMaterial::~NullPlasticMaterial()
{

}


double NullPlasticMaterial::getTrialPlasticStiffness()
{
  return 0;
}

PlasticHardeningMaterial * NullPlasticMaterial::getCopy()
{
     NullPlasticMaterial *nullPM = new  NullPlasticMaterial(this->getTag());
     return nullPM;
}
