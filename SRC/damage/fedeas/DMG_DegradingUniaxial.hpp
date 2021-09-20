/*
 * libdmg
 * TODO: Copyright, authors
 */

// Claudio Perez
#ifndef DMG_DEGRADING_UNIAXIAL_H
#define DMG_DEGRADING_UNIAXIAL_H
#include <functional>
#include <string>
#include <unordered_map>

#include <tcl.h>


struct UniaxialState {
/* This struct defines the interface between
 * the public wrapper and it's external 
 * implementation. */
 double  e, ep, De, se, kt, ke;
};

template<typename BaseMaterial>
// class DMG_DegradingUniaxial:
class Damage1D:
  public BaseMaterial
{
public:
  // inherit constructors
  using BaseMaterial::BaseMaterial;

  Damage1D<BaseMaterial>*
  getCopy(void) override {
    return new Damage1D<BaseMaterial>(*this);
  }

  UniaxialState past,pres;
  typedef std::function<int(void*, void*)> degrade_f;

  int setDamageWrapper(Tcl_Interp *interp, std::string tag) {

    typedef std::function<int(void*,void*)> degrade_t;
    // std::string tag = Tcl_GetString(objv[1]);
    
    auto map = (std::unordered_map<std::string, degrade_t(*)(Tcl_Interp*, std::string)>*)
      Tcl_GetAssocData(interp, "elle::libdmg::DamageGenerators", NULL);
    if (!map){
      printf("Failed to get constructor map");
      return -1;
    }

    this->degrade = (degrade_t) ((*map)["Uniaxial"](interp, tag));

    return 1;
  }
  

/*
  double getStrain();
  double getStress();
  double getTangent();
  double getInitialTangent();
*/

/* int commitState(void){
    puts("commit from wrapper");
    double prevStrain = BaseMaterial::getStrain();*/

  virtual int setTrialStrain(double trialStrain, double strainRate=0.0) override {
    
    double pastStrain = BaseMaterial::getStrain();
    BaseMaterial::setTrialStrain(trialStrain, strainRate);

    if (degrade){
      // double currStrain = BaseMaterial::getStrain();

      UniaxialState pres, past = {
        .e   = pastStrain,
        .ep  = trialStrain,
        .De  = trialStrain - pastStrain,
        .se  = BaseMaterial::getStress(),
        .kt  = BaseMaterial::getTangent(),
        .ke  = BaseMaterial::getInitialTangent()
      };

      degrade((void*)&past, (void*)&pres);

      this->sig = pres.se;
      this->e   = pres.kt;
    }
    return 0;
  }

private:
  using BaseMaterial::setTrialStrain;
  degrade_f degrade = NULL;
}; // class DMG_DegradingBeamColumn

#endif // DMG_DEGRADING_UNIAXIAL_H

