/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

/*                                                                        
** $Revision: 1.12 $
** $Date: 2010-03-05 22:32:36 $
** $Source: /usr/local/cvs/OpenSees/SRC/api/elementAPI.cpp,v $
                                                                        
** Written: fmk 
*/

#include "elementAPI.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//#include <packages.h>
#include <Domain.h>
#include <Node.h>
//#include <TclModelBuilder.h>
//#include <WrapperElement.h>

#include <map>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
//#include <WrapperUniaxialMaterial.h>
//#include <WrapperNDMaterial.h>

#include <OPS_Globals.h>

#include <CrdTransf.h>

typedef struct elementFunction {
  char *funcName;
  eleFunct theFunct; 
  struct elementFunction *next;
} ElementFunction;

typedef struct materialFunction {
  char *funcName;
  matFunct theFunct; 
  struct materialFunction *next;
} MaterialFunction;

static ElementFunction *theElementFunctions = NULL;
static MaterialFunction *theMaterialFunctions = NULL;
static Domain *theDomain = 0;
static int currentArg = 0;
static int maxArg = 0;
extern FE_Datastore *theDatabase;
//static int uniaxialMaterialObjectCount =0;

modelState theModelState;

struct cmp_str {
  bool operator()(const char *a, const char *b) {
    return strcmp(a,b)<0;
  }
};

std::map<char *, eleFunct, cmp_str>theEleFunctions;              // map of user added ele functions
std::map<char *, eleFunct, cmp_str>theUniaxialMaterialFunctions; // map of user added material functions

//std::map<int, UniaxialMaterial *>theUniaxialMaterials;           // map for UniaxialMaterial objects needed by user added ele functions'

static 
void OPS_InvokeMaterialObject(struct matObject *theMat, modelState *theModel,double *strain, double *tang, double *stress, int *isw, int *result)
{
  int matType = theMat->theParam[0];

  if (matType == 1) {
    //  UniaxialMaterial *theMaterial = theUniaxialMaterials[matCount];
    UniaxialMaterial *theMaterial = (UniaxialMaterial *)theMat->matObjectPtr;
    if (theMaterial == 0) {
      *result = -1;
      return;
    }
    
    if (*isw == ISW_COMMIT) {
      *result =  theMaterial->commitState();
      return;
    } else if (*isw == ISW_REVERT) {
      *result =  theMaterial->revertToLastCommit();
      return;
    } else if (*isw == ISW_REVERT_TO_START) {
      *result =  theMaterial->revertToStart();
      return;
    } else if (*isw == ISW_FORM_TANG_AND_RESID) {
      double matStress = 0.0;
      double matTangent = 0.0;
      int res = theMaterial->setTrial(strain[0], matStress, matTangent);
      stress[0] = matStress;
      tang[0] = matTangent;
      *result = res;
      return;
    }
  }
  
  return;
}

extern "C" 
int OPS_Error(char *errorMessage, int length)
{
  opserr << errorMessage;
  opserr << endln;
  return 0;
}

extern "C"   
int OPS_GetNumRemainingInputArgs()
{
  return maxArg-currentArg;
}

extern "C"   
int OPS_GetIntInput(int *numData, int*data)
{
  return 0;
}

extern "C" 
int OPS_GetDoubleInput(int *numData, double *data)
{
  return 0;  
}



extern "C" 
const char *OPS_GetString(void)
{
  return 0;  
}


int OPS_GetStringCopy(char **arrayData)
{
  return 0;  
}


extern "C" 
matObj *OPS_GetMaterial(int *matTag, int *matType)
{

  if (*matType == OPS_UNIAXIAL_MATERIAL_TYPE) {
    UniaxialMaterial *theUniaxialMaterial = OPS_getUniaxialMaterial(*matTag);
    
    if (theUniaxialMaterial != 0) {
      
      UniaxialMaterial *theCopy = theUniaxialMaterial->getCopy();
      //  uniaxialMaterialObjectCount++;
      // theUniaxialMaterials[uniaxialMaterialObjectCount] = theCopy;
      
      matObject *theMatObject = new matObject;
      theMatObject->tag = *matTag;
      theMatObject->nParam = 1;
      theMatObject->nState = 0;
      
      theMatObject->theParam = new double[1];
      //  theMatObject->theParam[0] = uniaxialMaterialObjectCount;
      theMatObject->theParam[0] = 1; // code for uniaxial material
      
      theMatObject->tState = 0;
      theMatObject->cState = 0;
      theMatObject->matFunctPtr = OPS_InvokeMaterialObject;

      theMatObject->matObjectPtr = theCopy;
      
      return theMatObject;
    }
    
    fprintf(stderr,"getMaterial - no uniaxial material exists with tag %d\n", *matTag);    
    return 0;

  } else if (*matType == OPS_SECTION_TYPE) {
    fprintf(stderr,"getMaterial - not yet implemented for Section\n");    
    return 0;
  } else {

    //    NDMaterial *theNDMaterial = theModelBuilder->getNDMaterial(*matTag);

    //    if (theNDMaterial != 0) 
      //      theNDMaterial = theNDMaterial->getCopy(matType);
      //    else {
      //      fprintf(stderr,"getMaterial - no nd material exists with tag %d\n", *matTag);          
      //      return 0;
      //    }

      //    if (theNDMaterial == 0) {
    //      fprintf(stderr,"getMaterial - material with tag %d cannot deal with %d\n", *matTag, matType);          
    //      return 0;
    //    }

    fprintf(stderr,"getMaterial - not yet implemented for nDMaterial\n");    
    return 0;
  }

  fprintf(stderr,"getMaterial - unknown material type\n");    
  return 0;

}

/*
extern "C" 
void OPS_GetMaterialPtr(int *matTag, matObj *theRes)
{
  UniaxialMaterial *theUniaxialMaterial = theModelBuilder->getUniaxialMaterial(*matTag);

  if (theUniaxialMaterial != 0) {

    UniaxialMaterial *theCopy = theUniaxialMaterial->getCopy();
    if (theCopy  == 0) {
      fprintf(stderr,"OPS_GetMaterialPtr() failed - no material of type %d \n", *matTag);      
      theRes = 0;
      return;
    }

    uniaxialMaterialObjectCount++;
    theUniaxialMaterials[uniaxialMaterialObjectCount] = theCopy;

    matObject *theMatObject = new matObject;
    theMatObject->tag = *matTag;
    theMatObject->nParam = 1;
    theMatObject->nState = 0;

    theMatObject->theParam = new double[1];
    theMatObject->theParam[0] = uniaxialMaterialObjectCount;

    theMatObject->tState = 0;
    theMatObject->cState = 0;
    theMatObject->matFunctPtr = OPS_UniaxialMaterialFunction;

    theRes = theMatObject;
  }

  theRes = 0;
}
*/


extern "C" 
eleObj *OPS_GetElement(int *eleTag) {
  return 0;
}

extern "C" 
eleObj *OPS_GetElementType(char *type, int sizeType) {

  // try existing loaded routines

  ElementFunction *eleFunction = theElementFunctions;
  bool found = false;
  while (eleFunction != NULL && found == false) {
    if (strcmp(type, eleFunction->funcName) == 0) {
      
      // create a new eleObject, set the function ptr &  return it
      
      eleObj *theEleObject = new eleObj;
      theEleObject->eleFunctPtr = eleFunction->theFunct;
      return theEleObject;
    }
    else
      eleFunction = eleFunction->next;
  }

  // ty to load new routine from dynamic library in load path
  
  eleFunct eleFunctPtr;
  void *libHandle;

  return 0;
}

extern "C" 
matObj *OPS_GetMaterialType(char *type, int sizeType) {

  // try existing loaded routines
  MaterialFunction *matFunction = theMaterialFunctions;
  bool found = false;
  while (matFunction != NULL && found == false) {
    if (strcmp(type, matFunction->funcName) == 0) {
      
      // create a new eleObject, set the function ptr &  return it
      
      matObj *theMatObject = new matObj;
      theMatObject->matFunctPtr = matFunction->theFunct;
      /* opserr << "matObj *OPS_GetMaterialType() - FOUND " << endln;  */
      return theMatObject;
    }
    else
      matFunction = matFunction->next;
  }

  // ty to load new routine from dynamic library in load path
  matFunct matFunctPtr;
  void *libHandle;
  
  

  return 0;
}

extern "C" 
int OPS_AllocateMaterial(matObject *theMat){

  /*fprintf(stderr,"allocateMaterial Address %p\n",theMat);*/

  if (theMat->nParam > 0)
    theMat->theParam = new double[theMat->nParam];

  int nState = theMat->nState;

  if (nState > 0) {
    theMat->cState = new double[nState];
    theMat->tState = new double[nState];
    for (int i=0; i<nState; i++) {
      theMat->cState[i] = 0;
      theMat->tState[i] = 0;
    }
  } else {
    theMat->cState = 0;
    theMat->tState = 0;
  }

  return 0;
}  

extern "C" 
int OPS_AllocateElement(eleObject *theEle, int *matTags, int *matType){
  if (theEle->nNode > 0)
    theEle->node = new int[theEle->nNode];

  if (theEle->nParam > 0)
    theEle->param = new double[theEle->nParam];

  if (theEle->nState > 0) {
    theEle->cState = new double[theEle->nState];
    theEle->tState = new double[theEle->nState];
  }

  int numMat = theEle->nMat;
  if (numMat > 0)
    theEle->mats = new matObject *[numMat];

  
  for (int i=0; i< numMat; i++) {
  /*  opserr << "AllocateElement - matTag " << matTags[i] << "\n"; */

    matObject *theMat = OPS_GetMaterial(&(matTags[i]), matType);
    //    matObject *theMat = OPS_GetMaterial(&(matTags[i]));

    theEle->mats[i] = theMat;
  }

  return 0;
}  



extern "C" 
int OPS_GetNodeCrd(int *nodeTag, int *sizeCrd, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);
  if (theNode == 0) {
    opserr << "OPS_GetNodeCrd - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeCrd;
  const Vector &crd = theNode->getCrds();
  if (crd.Size() != size) {
    opserr << "OPS_GetNodeCrd - crd size mismatch\n";
    return -1;
  }
  for (int i=0; i < size; i++) 
    data[i] = crd(i);
    
  return 0;
}

extern "C" 
int OPS_GetNodeDisp(int *nodeTag, int *sizeData, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeDisp - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &disp = theNode->getTrialDisp();

  if (disp.Size() != size) {
    opserr << "OPS_GetNodeDisp - crd size mismatch\n";
    return -1;
  }
  for (int i=0; i < size; i++) 
    data[i] = disp(i);
    
  return 0;
}

extern "C" 
int OPS_GetNodeVel(int *nodeTag, int *sizeData, double *data)
{
	  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeVel - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &vel = theNode->getTrialVel();

  if (vel.Size() != size) {
    opserr << "OPS_GetNodeVel - crd size mismatch\n";
    return -1;
  }
  for (int i=0; i < size; i++) 
    data[i] = vel(i);
    
  return 0;
}

extern "C" 
int OPS_GetNodeAccel(int *nodeTag, int *sizeData, double *data)
{
   Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeAccel - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &accel = theNode->getTrialAccel();

  if (accel.Size() != size) {
    opserr << "OPS_GetNodeAccel - accel size mismatch\n";
    return -1;
  }
  for (int i=0; i < size; i++) 
    data[i] = accel(i);
    
  return 0;
}

extern "C" 
int OPS_GetNodeIncrDisp(int *nodeTag, int *sizeData, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeIncrDisp - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &disp = theNode->getIncrDisp();

  if (disp.Size() != size) {
    opserr << "OPS_GetNodeIncrDis - crd size mismatch\n";
    return -1;
  }
  for (int i=0; i < size; i++) 
    data[i] = disp(i);
    
  return 0;
}


extern "C" 
int OPS_GetNodeIncrDeltaDisp(int *nodeTag, int *sizeData, double *data)
{
  Node *theNode = theDomain->getNode(*nodeTag);

  if (theNode == 0) {
    opserr << "OPS_GetNodeIncrDisp - no node with tag " << *nodeTag << endln;
    return -1;
  }
  int size = *sizeData;
  const Vector &disp = theNode->getIncrDeltaDisp();

  if (disp.Size() != size) {
    opserr << "OPS_GetNodeIncrDis - crd size mismatch\n";
    return -1;
  }
  for (int i=0; i < size; i++) 
    data[i] = disp(i);
    
  return 0;
}




extern "C" int        
OPS_InvokeMaterial(eleObject *theEle, int *mat, modelState *model, double *strain, double *stress, double *tang, int *isw)
{
  return 0;
}

extern "C" int        
OPS_InvokeMaterialDirectly(matObject **theMat, modelState *model, double *strain, double *stress, double *tang, int *isw)
{
  return 0;
}


extern "C" int        
OPS_InvokeMaterialDirectly2(matObject *theMat, modelState *model, double *strain, double *stress, double *tang, int *isw)
{
		return -1;
}


UniaxialMaterial *
OPS_GetUniaxialMaterial(int matTag) {
  return 0;
}

NDMaterial *
OPS_GetNDMaterialPointer(int matTag)
{
  return 0;
}

CrdTransf * 
OPS_GetCrdTransfPtr(int tag)
{
  return 0;
}

SectionForceDeformation *
OPS_GetSectionForceDeformation(int matTag)
{
  return 0;
}


int OPS_ResetCurrentInputArg(int cArg)
{
    return 0;
}

int
OPS_ResetInput(Domain *domain)
{
  return 0;
}

int
OPS_ResetInputNoBuilder(Domain *domain)
{
  return 0;
}
	       
			       
int     
OPS_GetNDF()
{
  return 0;
}

int     
OPS_GetNDM()
{
  return 0;
}

FE_Datastore *OPS_GetFEDatastore()
{
	return 0;
}

const char *OPS_GetInterpPWD()
{
  return 0;
}

LimitCurve *OPS_GetLimitCurve(int LimCrvTag)
{
  return 0;
}
Domain *OPS_GetDomain(void) {
  return 0;
}
