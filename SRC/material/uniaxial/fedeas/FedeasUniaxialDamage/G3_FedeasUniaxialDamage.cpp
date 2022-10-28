#include <stdio.h>
#include <InputAPI.h>
#include "DegradingUniaxialWrapper.h"

#define WRAPPER_CMD "FedeasUniaxialDamage"
// #define WRAPPER_CMD "FedeasDamage"

UniaxialMaterial*
G3Parse_newFedeasUniaxialDamage(G3_Runtime* rt, int argc, TCL_Char **argv)
{
  // Pointer to a uniaxial material that will be returned
  DegradingUniaxialWrapper *theMaterial = 0;
  UniaxialMaterial *theWrappedMaterial = 0;
  int tags[2];

  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial " WRAPPER_CMD " $tag "
              "$wrapTag <-damage $damageTag>"
           << endln;
    return nullptr;
  }

  // Get wrapper tag
  if (G3Parse_getInt(rt, argv[2], &tags[0]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return nullptr;
  }
  // Get base tag
  if (G3Parse_getInt(rt, argv[3], &tags[1]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return nullptr;
  }

  // Get base material
  theWrappedMaterial = G3_getUniaxialMaterialInstance(rt, tags[1]);
  if (theWrappedMaterial == 0) {
    opserr << "WARNING unable to retrieve uniaxialMaterial with tag" WRAPPER_CMD " tag: "
           << tags[1] << endln;
    return nullptr;
  }

  int argn = 4;
  const char *dmgtag = 0;
  double Ccd = 0.5;
  StateOperator *damage = new StateOperator;
  while (argn < argc) {
    const char *param = argv[argn];

    if ((strcmp(param, "-damage") == 0) || 
        (strcmp(param, "-dmg") == 0)    ||
        (strcmp(param, "-DMG") == 0))   {
      *damage = *(StateOperator*)Tcl_GetAssocData(G3_getInterpreter(rt), 
                                                  "fedeas::damage::UniaxialDamage", NULL);
      Tcl_Interp* interp = G3_getInterpreter(rt);
      
      damage->call(damage, interp, ISW_CREATE, argc - argn, &argv[++argn], 0, 0, 0, 0, 0);
      damage->call(damage, interp, ISW_MALLOC, 0, 0, 0, 0, 0, 0, 0);

    } else if ((strcmp(param, "-couple") == 0) || 
               (strcmp(param, "-ccd") == 0)    ||
               (strcmp(param, "-Ccd") == 0))   {
      Ccd = std::stod(argv[++argn]);
    } else {
      break;
    }
    argn++;
  }

  // Parsing was successful, allocate the material
  theMaterial = new DegradingUniaxialWrapper(tags[0], *theWrappedMaterial, damage);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type " WRAPPER_CMD << endln;
    return nullptr;
  }
  theMaterial->setCoupling(Ccd);

  // if (dmgtag){
  //   if (theMaterial->setDamageWrapper(G3_getInterpreter(rt), dmgtag) > 0)
  //     opserr << "#Set damage wrapper '" << dmgtag << "'\n";
  // }

  // return G3_addUniaxialMaterial(rt, theMaterial);
  return theMaterial;
}

