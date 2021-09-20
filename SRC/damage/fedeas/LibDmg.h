#ifndef DMGAPI_H
#define DMGAPI_H
#include <functional>
#include <fedeas/DMG_DegradingUniaxial.hpp>

typedef std::function<int(void*, void*)> degrade_f;

degrade_f* getNMM_Degradation(Tcl_Interp*, const char *);
degrade_f* getUniaxialDegradation(Tcl_Interp*, const char *);

/*
static int
Damage_Cmd([[maybe_unused]] ClientData cdata, Tcl_Interp *interp, [[maybe_unused]] int objc, [[maybe_unused]] Tcl_Obj *const objv[])
{
  typedef std::function<int(void*,void*)> degrade_t;
  std::string tag = Tcl_GetString(objv[1]);
  puts("tst::material called");
  auto map = (std::unordered_map<std::string, degrade_t(*)(Tcl_Interp*, std::string)>*)
    Tcl_GetAssocData(interp, "elle::libdmg::DamageGenerators", NULL);
  if (!map)
    printf("Failed to get constructor map");
  else
    printf("Retrieved uniaxial constructor");
  auto degrade = (degrade_t) ((*map)["Uniaxial"](interp, tag));


  // std::function<int(void*,void*)> degrade = (*get)(interp, tag);

  return TCL_OK;
}
*/

#endif // DMGAPI_H
