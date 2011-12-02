/*
 *  DRMInputHandler.cpp
 *  
 *
 *  Created by george  petropoulos on 2/10/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "DRMInputHandler.h"

DRMInputHandler::DRMInputHandler(int tag, char** in_files, int files,  double dt, double* time_array,
				 int num_steps, int* file_data, int fileData_size, double* domain_crds, double* drm_box_crds,
				 Mesh3DSubdomain* my_mesher)
{
    
  this->myMesher = my_mesher;
  this->timeArray = time_array;
  this->deltaT = dt;
  this->numSteps = num_steps;
  
  this->domain_Crds = domain_crds;
  this->drm_box_Crds = drm_box_crds;
  
  this->numFiles = files;
  this->filePtrs = new char*[this->numFiles];
  for (int i=0; i<this->numFiles; i++)
    this->filePtrs[i] = in_files[i];
  this->fileData = new int[fileData_size];
  for (int i=0; i<fileData_size; i++)
    fileData[i] = file_data[i];
  
}


DRMInputHandler::~DRMInputHandler()
{
  delete [] filePtrs;
  delete [] fileData;
  delete [] domain_Crds;
  delete [] drm_box_Crds;
  if (timeArray != 0)
    delete [] timeArray;
}

void DRMInputHandler::seteNodeMap(std::map<int,int>& enodes)
{
  double xMin = this->drm_box_Crds[0]; double xMax = this->drm_box_Crds[1];
  double yMin = this->drm_box_Crds[2]; double yMax = this->drm_box_Crds[3];
  double zMin = this->drm_box_Crds[4]; double zMax = this->drm_box_Crds[5];
  
  this->myMesher->allocate_e_Nodes(xMin, xMax,
				   yMin, yMax,
				   zMin, zMax,
				   enodes);
}

void DRMInputHandler::seteleMap(std::map<int,Element* >& ele, std::map<int,Vector*>& stor, std::map<int,int>& stor2)
{
  double xMin = this->drm_box_Crds[0]; double xMax = this->drm_box_Crds[1];
  double yMin = this->drm_box_Crds[2]; double yMax = this->drm_box_Crds[3];
  double zMin = this->drm_box_Crds[4]; double zMax = this->drm_box_Crds[5];
  
  this->myMesher->allocateBoundaryLayerElements(xMin, xMax,
						yMin, yMax,
						zMin, zMax,
						ele,
						stor,
						stor2);
  this->ele_str = stor;
  this->ele_str2 = stor2;
}


void DRMInputHandler::getMotions(Element*  eleTag, double time, Vector& U, Vector& Ud, Vector& Udd)
{
  
}
