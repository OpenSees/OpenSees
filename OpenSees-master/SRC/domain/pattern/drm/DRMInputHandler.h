/*
 *  DRMInputHandler.h
 *  
 *
 *  Created by george  petropoulos on 2/10/06.
 *  <gnp> <petropoulos@gmail.com>
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/* 
* @author: gnp <petropoulos@gmail.com>
 *
 * @Description: Main DRM Functionality
 *
 * @Date: 2/10/06
 */


#ifndef DRMInputHandler_h
#define DRMInputHandler_h


#include <LoadPattern.h>
#include <Domain.h>
#include "DRMBoundaryLayerDecorator.h"
#include "Mesh3DSubdomain.h"
#include <Vector.h>
#include <Element.h>
#include <Matrix.h>
#include <map>

class DRMInputHandler
{
  

 public:

  DRMInputHandler(int tag, char** in_files, int files,  double dt, double* time_array,
		  int num_steps, int* file_data, int fileData_size, double* domain_crds, double* drm_box_crds,
		  Mesh3DSubdomain* my_mesher);
  
  virtual ~DRMInputHandler();
  
  void seteNodeMap(std::map<int,int>& nodes);
  
  void seteleMap(std::map<int,Element*>& ele, std::map<int,Vector*>& stor, std::map<int,int>& stor2); 
   
  virtual void getMotions(Element* eleTag, double time, Vector& U, Vector& Ud, Vector& Udd);
  
 protected:
  
  Domain* myDomain;
  
  
  double* timeArray;
  int numSteps;
  double deltaT;
  
  double startTime;
  
  Mesh3DSubdomain* myMesher;
  
  double* domain_Crds;
  double* drm_box_Crds;
  
  char** filePtrs;
  int numFiles;
  int* fileData;

  std::map<int, Vector*> ele_str;
  std::map<int, int> ele_str2;
};
#endif
