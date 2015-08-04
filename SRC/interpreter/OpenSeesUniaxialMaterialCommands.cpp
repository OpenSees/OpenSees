
#include <UniaxialMaterial.h>
#include <string.h>

extern void *OPS_ElasticMaterial(void);
extern void *OPS_ElasticPPMaterial(void);
extern void *OPS_ParallelMaterial(void);

UniaxialMaterial *
OPS_ParseUniaxialMaterialCommand(const char *matType) {
  UniaxialMaterial  *theMaterial = 0;

  if (strcmp(matType,"Elastic") == 0) {
    
    void *theMat = OPS_ElasticMaterial();
    if (theMat != 0) 
      theMaterial = (UniaxialMaterial *)theMat;
    else 
      return 0;
    
  } else if (strcmp(matType,"ElasticPP") == 0) {
    
    void *theMat = OPS_ElasticPPMaterial();
    if (theMat != 0) 
      theMaterial = (UniaxialMaterial *)theMat;
    else 
      return 0;
    
  } else if (strcmp(matType,"Parallel") == 0) {
    void *theMat = OPS_ParallelMaterial();
    if (theMat != 0) 
      theMaterial = (UniaxialMaterial *)theMat;
    else 
      return 0;
  }

  return theMaterial;
}
