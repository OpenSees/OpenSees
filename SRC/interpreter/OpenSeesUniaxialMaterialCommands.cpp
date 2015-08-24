#include <UniaxialMaterial.h>
#include <string.h>
#include <map>


struct char_cmp { 
  bool operator () (const char *a,const char *b) const 
  {
    return strcmp(a,b)<0;
  } 
};

typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;

static OPS_ParsingFunctionMap uniaxialMaterialsMap;
static bool initDone = false;

extern void *OPS_ElasticMaterial(void);
extern void *OPS_ElasticPPMaterial(void);
extern void *OPS_ParallelMaterial(void);

int OPS_SetUpUniaxialMaterials(void) {
  uniaxialMaterialsMap.insert(std::make_pair("Elastic", &OPS_ElasticMaterial));
  uniaxialMaterialsMap.insert(std::make_pair("ElasticPP", &OPS_ElasticPPMaterial));
  uniaxialMaterialsMap.insert(std::make_pair("Parallel", &OPS_ParallelMaterial));

  return 0;
}

UniaxialMaterial *
OPS_ParseUniaxialMaterialCommand(const char *matType) {

  if (initDone == false) {
    OPS_SetUpUniaxialMaterials();
    initDone = true;
  }

  UniaxialMaterial  *theMaterial = 0;
  OPS_ParsingFunctionMap::const_iterator iter = uniaxialMaterialsMap.find(matType);
  if (iter == uniaxialMaterialsMap.end()) {
    return 0;
  }

  void *theMat = (*iter->second)();
  if (theMat != 0) 
    theMaterial = (UniaxialMaterial *)theMat;
  else 
    return 0;

  return theMaterial;

  /* OLD CODE .. IF ELSE
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
  */
}
