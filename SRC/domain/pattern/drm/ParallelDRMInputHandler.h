/*
 *  PlaneDRMInputHandler.h
 *  
 *
 *  Created by george  petropoulos on 3/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PlaneDRMInputHandler_h
#define PlaneDRMInputHandler_h

#include <fstream>
#include <cstdlib>
#include <Domain.h>
#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include "DRMInputHandler.h"
#include "GeometricBrickDecorator.h"
#include "Mesh3DSubdomain.h"
#include <math.h>

class PlaneDRMInputHandler : public DRMInputHandler {
  
 public:
    
  PlaneDRMInputHandler(int tag, char** in_files, int files,  double dt, double* time_array, 
			  int num_steps, int* file_data, int fileData_size, double* domain_crds, double* drm_box_crds, double* eleD,
			  Mesh3DSubdomain* my_mesher, int steps_to_cache, Domain* domain);
  
  virtual ~PlaneDRMInputHandler();
  
  void populateBuffers();
  void getMotions(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd);
  void computeHistory(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd, bool updateDm1); 
  void handle_elementAtface5(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd); 
  void handle_elementAtface1(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd); 
  void handle_elementAtface2(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd); 
  void handle_elementAtface3(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd); 
  void handle_elementAtface4(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd); 
  void getLocations(double x, double y, double dx, double dy, int* xloc, int* yloc);
  void getTemporal(double time, int* tloc);
  int getIndex(double time);
  void getf5pointer(Node* node_tag, int local_tag, int index);
  void getf1pointer(Node* node_tag, int local_tag, int index);
  void getf2pointer(Node* node_tag, int local_tag, int index);
  void getf3pointer(Node* node_tag, int local_tag, int index);
  void getf4pointer(Node* node_tag, int local_tag, int index);
  void pointerCopy(int node_from, int node_to);
  void populateTempBuffers(int index, int fileptr, double ksi, double eta);

  
  private :
    
    bool initial;
  
  double* f1buffer;
  double* f2buffer;
  double* f3buffer;
  double* f4buffer;
  double* f5buffer;

  double** buffers;
  double* eleD;
  
  double* timeBuf;
  
  Domain* myDomain;

  std::ifstream ifile1;
  std::ifstream ifile2;
  std::ifstream ifile3;
  std::ifstream ifile4;
  std::ifstream ifile5a;
  std::ifstream ifile5b;
  

  static Vector Vtm2;
  static Vector Vtm1;
  static Vector Vtp1;
  static Vector Vtp2;

  static Vector Vtm2_n1;
  static Vector Vtm2_n2;
  static Vector Vtm2_n3;
  static Vector Vtm2_n4;
  
  static Vector Vtm1_n1;
  static Vector Vtm1_n2;
  static Vector Vtm1_n3;
  static Vector Vtm1_n4;

  static Vector Vtp1_n1;
  static Vector Vtp1_n2;
  static Vector Vtp1_n3;
  static Vector Vtp1_n4;

  static Vector Vtp2_n1;
  static Vector Vtp2_n2;
  static Vector Vtp2_n3;
  static Vector Vtp2_n4;

  static Vector Vtempm2;
  static Vector Vtempm1;
  static Vector Vtempp1;
  static Vector Vtempp2;


  
  
  int* which;
  int cacheValue;
  int globalCounter;
  int localCounter;
  GeometricBrickDecorator* myDecorator;
  
  // some helpers perhaps
};
#endif
