/*
 *  PlaneDRMInputHandler.cpp
 *  
 *
 *  Created by george  petropoulos on 3/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "PlaneDRMInputHandler.h"

Vector PlaneDRMInputHandler::Vtm2(24);
Vector PlaneDRMInputHandler::Vtm1(24);
Vector PlaneDRMInputHandler::Vtp1(24);
Vector PlaneDRMInputHandler::Vtp2(24);

Vector PlaneDRMInputHandler::Vtm2_n1(3);
Vector PlaneDRMInputHandler::Vtm2_n2(3);
Vector PlaneDRMInputHandler::Vtm2_n3(3);
Vector PlaneDRMInputHandler::Vtm2_n4(3);
  
Vector PlaneDRMInputHandler::Vtm1_n1(3);
Vector PlaneDRMInputHandler::Vtm1_n2(3);
Vector PlaneDRMInputHandler::Vtm1_n3(3);
Vector PlaneDRMInputHandler::Vtm1_n4(3);

Vector PlaneDRMInputHandler::Vtp1_n1(3);
Vector PlaneDRMInputHandler::Vtp1_n2(3);
Vector PlaneDRMInputHandler::Vtp1_n3(3);
Vector PlaneDRMInputHandler::Vtp1_n4(3);

Vector PlaneDRMInputHandler::Vtp2_n1(3);
Vector PlaneDRMInputHandler::Vtp2_n2(3);
Vector PlaneDRMInputHandler::Vtp2_n3(3);
Vector PlaneDRMInputHandler::Vtp2_n4(3);

Vector PlaneDRMInputHandler::Vtempm2(3);
Vector PlaneDRMInputHandler::Vtempm1(3);
Vector PlaneDRMInputHandler::Vtempp1(3);
Vector PlaneDRMInputHandler::Vtempp2(3);





PlaneDRMInputHandler::PlaneDRMInputHandler(int tag, char** in_files, int files,  double dt, double* time_array,
					   int num_steps, int* file_data, int fileData_size,
					   int _nd1, int _nd2,
					   double* domain_crds, double* drm_box_crds, double* eleD,
					   Mesh3DSubdomain* my_mesher, int steps_to_cache, Domain* domain)
  
  : DRMInputHandler(tag, in_files, files,  dt, time_array,
		    num_steps, file_data, fileData_size, domain_crds,drm_box_crds,
		    my_mesher)
{
  /*
   *  Current implementations assumes num_steps + 2 = total actual time steps,
   *  included in the files where the first one ist he zero time step
   *  ie time = 0.0 at the first time step in the files
   *
   *  eleD refers to the dx,dy,dz of the cmu mesh
   *
   *  Berkeley coordinate system, ie actual models mesh assumes crd system origin to be lower left corner of DRM box
   *
   *  always need to have node at origin of cmu mesh
   *
   *  file data always refers to cmu mesh: ie: fileData[0] is num nodes in face 1 of cmu, 
   *                                           fileData[1] is num nodes along horiz traversal of face 1
   *                                           fileData[2] is num nodes along vert traversal of face 1
   *                                           and so on for fileData[...]
   *
   *
   */
  

  if (steps_to_cache < 1) {
    opserr << " need steps_to_cache at least >=1 ";
    steps_to_cache = 1;
  }
  this->cacheValue = steps_to_cache;
  
  this->f1buffer = new double[3*(this->cacheValue + 4)*this->fileData[0]];		
  this->f2buffer = new double[3*(this->cacheValue + 4)*this->fileData[3]];
  this->f3buffer = new double[3*(this->cacheValue + 4)*this->fileData[6]];
  this->f4buffer = new double[3*(this->cacheValue + 4)*this->fileData[9]];
  this->f5buffer = new double[3*(this->cacheValue + 4)*this->fileData[12]];
 

  this->buffers = new double*[5];

  this->buffers[0] = f1buffer;
  this->buffers[1] = f2buffer;
  this->buffers[2] = f3buffer;
  this->buffers[3] = f4buffer;
  this->buffers[4] = f5buffer;
 

  if (f1buffer==0 || f2buffer==0 || f3buffer==0 || f4buffer==0 || f5buffer==0) {
    opserr << "Error in memory allocations for DRM Load pattern, try smaller caching " << endln;
    if (f1buffer !=0)
      delete [] f1buffer;
    if (f2buffer !=0)
      delete [] f2buffer;
    if (f3buffer !=0)
      delete [] f3buffer;
    if (f4buffer !=0)
      delete [] f4buffer;
    if (f5buffer !=0)
      delete [] f5buffer;
     opserr << " Need abort ";
    delete [] buffers;
    exit(-1);
  }
	
  // Open File streams
  
  this->ifile1.open(filePtrs[0]);
  this->ifile2.open(filePtrs[1]);
  this->ifile3.open(filePtrs[2]);
  this->ifile4.open(filePtrs[3]);
  this->ifile5a.open(filePtrs[4]);
  this->ifile5b.open(filePtrs[5]);

  if (ifile1.bad() ) {
    opserr << " Bad file 1 " << endln;
    exit(-1);
  }
  if (ifile2.bad() ) {
    opserr << " Bad file 2 " << endln;
    exit(-1);
  }
  if (ifile3.bad() ) {
    opserr << " Bad file 3 " << endln;
    exit(-1);
  }
  if (ifile4.bad() ) {
    opserr << " Bad file 4 " << endln;
    exit(-1);
  }
  if (ifile5a.bad() ) {
    opserr << " Bad file 5a " << endln;
    exit(-1);
  }
  if (ifile5b.bad() ) {
    opserr << " Bad file 5b " << endln;
    exit(-1);
  }

  this->cacheValue = cacheValue;
  this->initial = true;
  
  this->timeBuf = new double[this->cacheValue+4];

  // set to ONE to be comparable to cacheValue
  this->localCounter = 1;
  this->globalCounter = 0;
  
  this->myDecorator = new GeometricBrickDecorator();
  this->myDecorator->setDomain(domain);
  this->eleD = eleD;
  this->myDomain = domain;
  
  this->which = new int[4];

  Vtm2.Zero();
  Vtm1.Zero();
  Vtp1.Zero();
  Vtp2.Zero();

  nd1 = _nd1;
  nd2 = _nd2;

  populateBuffers();
}

PlaneDRMInputHandler::~PlaneDRMInputHandler()
{
  delete [] f1buffer;
  delete [] f2buffer;
  delete [] f3buffer;
  delete [] f4buffer;
  delete [] f5buffer;
  
  delete [] buffers;
  delete [] which;
  delete myDecorator;
}



void PlaneDRMInputHandler::populateBuffers()
{
//  int pid;
//  bool debug = false;
  
  if (this->initial) {
    
    int temp = this->fileData[0];
    double dataIn;
    for (int i=0; i<3*temp; i++) {
      f1buffer[i] = 0.0;
    }
    for (int i=3*temp; i<3*(cacheValue + 4)*temp; i++) {
      ifile1 >> dataIn;
      f1buffer[i] = dataIn;
    }
    temp = this->fileData[3];
    for (int i=0; i<3*temp; i++) {
      f2buffer[i] = 0.0;
    }
    for (int i=3*temp; i<3*(cacheValue + 4)*temp; i++) {
      ifile2 >> dataIn;
      f2buffer[i] = dataIn;
    }
    temp = this->fileData[6];
    for (int i=0; i<3*temp; i++) {
      f3buffer[i] = 0.0;
    }
    for (int i=3*temp; i<3*(cacheValue + 4)*temp; i++) {
      ifile3 >> dataIn;
      f3buffer[i] = dataIn;
    }
    temp = this->fileData[9];
    for (int i=0; i<3*temp; i++) {
      f4buffer[i] = 0.0;
    }
    for (int i=3*temp; i<3*(cacheValue + 4)*temp; i++) {
      ifile4 >> dataIn;
      f4buffer[i] = dataIn;
    }
    temp = this->fileData[12];
    for (int i=0; i<3*temp; i++) {
      f5buffer[i] = 0.0;
    }
    int index = 3*temp;
    for (int k=1; k<(cacheValue +4); k++) {
      for (int i=0; i<3*nd1; i++) {
	ifile5a >> dataIn;
	f5buffer[index++] = dataIn;
      }
      for (int i=0; i<3*nd2; i++) {
	ifile5b >> dataIn;
	f5buffer[index++] = dataIn;
      }
    } 
    initial = false;
    globalCounter += cacheValue + 1;
    
    //update time buff
    timeBuf[0] = -deltaT;
    timeBuf[1] = 0.0;
    timeBuf[2] = deltaT;
    for (int i=0; i<this->cacheValue+1; i++)
      timeBuf[3+i] = timeBuf[2+i]+this->deltaT;
    
  }
  else {
    int rem = numSteps - globalCounter;
    if ( rem < cacheValue ) {
      if (rem < 0)
	return;
    }
    else {
      rem = cacheValue;
    }
    int temp = this->fileData[0];
    for (int i=0; i<3*temp; i++) {
      f1buffer[i] = f1buffer[i+ 3*(cacheValue+1)*temp];
      f1buffer[i+3*temp] = f1buffer[i+3*(cacheValue+2)*temp];			
      f1buffer[i+6*temp] = f1buffer[i+3*(cacheValue+3)*temp];			
    }
    double dataIn=0.0;
    for (int i=9*temp; i<3*(rem + 4)*temp; i++) {
      ifile1 >> dataIn;
      f1buffer[i] = dataIn;
    }
    temp = this->fileData[3];
    for (int i=0; i<3*temp; i++) {
      f2buffer[i] = f2buffer[i+ 3*(cacheValue+1)*temp];
      f2buffer[i+3*temp] = f2buffer[i+3*(cacheValue+2)*temp];			
      f2buffer[i+6*temp] = f2buffer[i+3*(cacheValue+3)*temp];			
    }
    dataIn=0.0;
    for (int i=9*temp; i<3*(rem + 4)*temp; i++) {
      ifile2 >> dataIn;
      f2buffer[i] = dataIn;
    }
    temp = this->fileData[6];
    for (int i=0; i<3*temp; i++) {
      f3buffer[i] = f3buffer[i+ 3*(cacheValue+1)*temp];
      f3buffer[i+3*temp] = f3buffer[i+3*(cacheValue+2)*temp];			
      f3buffer[i+6*temp] = f3buffer[i+3*(cacheValue+3)*temp];			
    }
    dataIn=0.0;
    for (int i=9*temp; i<3*(rem + 4)*temp; i++) {
      ifile3 >> dataIn;
      f3buffer[i] = dataIn;
    }
    temp = this->fileData[9];
    for (int i=0; i<3*temp; i++) {
      f4buffer[i] = f4buffer[i+ 3*(cacheValue+1)*temp];
      f4buffer[i+3*temp] = f4buffer[i+3*(cacheValue+2)*temp];			
      f4buffer[i+6*temp] = f4buffer[i+3*(cacheValue+3)*temp];			
    }
    dataIn=0.0;
    for (int i=9*temp; i<3*(rem + 4)*temp; i++) {
      ifile4 >> dataIn;
      f4buffer[i] = dataIn;
    }
    temp = this->fileData[12];
    for (int i=0; i<3*temp; i++) {
      f5buffer[i] = f5buffer[i+ 3*(cacheValue+1)*temp];
      f5buffer[i+3*temp] = f5buffer[i+3*(cacheValue+2)*temp];			
      f5buffer[i+6*temp] = f5buffer[i+3*(cacheValue+3)*temp];			
    }
    dataIn=0.0;
    int index = 9*temp;
    //need some change here for two files
    //for (int i=6*temp; i<3*(rem + 3)*temp; i++) {
    for (int k=3; k<(cacheValue +4); k++) {
      for (int i=0; i<3*nd1; i++) {
	ifile5a >> dataIn;
	f5buffer[index++] = dataIn;
      }
      for (int i=0; i<3*nd2; i++) {
	ifile5b >> dataIn;
	f5buffer[index++] = dataIn;
      }
    }
    globalCounter += cacheValue+1;

    //update time buff
    timeBuf[0] = timeBuf[this->cacheValue+1];
    timeBuf[1] = timeBuf[this->cacheValue+2];
    timeBuf[2] = timeBuf[this->cacheValue+3];
    for (int i=0; i<this->cacheValue+1; i++)
      timeBuf[3+i] = timeBuf[2+i]+this->deltaT;
    
  }
}

void PlaneDRMInputHandler::getMotions(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd)
{
  
//  int tg = eletag->getTag();

  if (time > numSteps*deltaT)
    return;
  ///Start by finding the face to which the element belongs;
  // check for bottom face i.e f5
  // check for 1,2,3,4
  double xMin,xMax,yMin,yMax,zMin,zMax;
  xMin = this->drm_box_Crds[0];
  xMax = this->drm_box_Crds[1];
  yMin = this->drm_box_Crds[2];
  yMax = this->drm_box_Crds[3];
  zMin = this->drm_box_Crds[4];
  zMax = this->drm_box_Crds[5];
  
  this->myDecorator->setBrick(eletag);
  int result = 0;
  if (this->myDecorator->isLeftBoundary(xMin, xMax,  yMin,  yMax,  zMin,  zMax)) 
    result =3;
  if (this->myDecorator->isRightBoundary(xMin,xMax,  yMin,  yMax,  zMin,  zMax)) 
    result =5;
  if (this->myDecorator->isFrontBoundary(xMin,xMax,  yMin,  yMax,  zMin,  zMax)) 
    result =7;
  if (this->myDecorator->isRearBoundary(xMin, xMax,  yMin,  yMax,  zMin,  zMax)) 
    result =11;
  if (this->myDecorator->isBottomBoundary(xMin,xMax, yMin,  yMax,  zMin,  zMax)) 
    result =1;

  switch (result) {
  case 1:
    // face 5 of cmu that is bottom face of ucb
    handle_elementAtface5(eletag, time, U, Ud, Udd);
    break;
  case 5:
    // face 1 of cmu that is right face of ucb
    handle_elementAtface1(eletag, time, U, Ud, Udd);
    break;
  case 3:
    // face 2 of cmu thas is left face of ucb
    handle_elementAtface2(eletag, time, U, Ud, Udd);
    break;
  case 7:
    // face 3 of cmu that is front face of ucb
    handle_elementAtface3(eletag, time, U, Ud, Udd);
    break;
  case 11:
    // face 4 of cmu that is rear face of ucb
    handle_elementAtface4(eletag, time, U, Ud, Udd);
    break;
  default:
    opserr << " SHOULDN'T SEE THIS \n";
    break;
  }
}


void PlaneDRMInputHandler::computeHistory(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd, bool updateDm1) 
{
  double oo2dt = 0.5/deltaT;

  int per = floor(time/deltaT);
  double mt = (time -((double) per)*deltaT)/deltaT;


  // linear interpolation for the velocities
  Ud = (1.0-mt)*Vtm1 + mt*Vtp1;
  
  // central differences for the accelerations
  Udd = (1.0 -mt)*oo2dt*(Vtp1-Vtm2) + mt*oo2dt*(Vtp2-Vtm1);

  // trapezoidal rule for the displacements
  // old way with ptrs kept in element.cpp
  //Vector& Dm1 = eletag->getDtm1();

  // get ele* tag
  int tag = eletag->getTag();
  Vector& Dm1 = *((this->ele_str)[tag]);


  if (updateDm1)
    Dm1 += 0.5*deltaT*(Vtm2 + Vtm1);
  U = Dm1 + mt*0.5*deltaT*(Vtp1+Vtm1);
}


void PlaneDRMInputHandler::handle_elementAtface5(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd) 
{
  bool updateDm1 =  false;
  Node** nodes = eletag->getNodePtrs();
  int index = getIndex(time);
//   int lastIndex = eletag->getLastDRMIndex();
//   if (index != lastIndex) {
//     updateDm1 = true;
//     eletag->setLastDRMIndex(index);
//   }

  int tag = eletag->getTag();
  int lastIndex = ele_str2[tag];
  if (index != lastIndex) {
    updateDm1 = true;
    ele_str2[tag] = index;
  }


  getf5pointer(nodes[0],0,index);
  getf5pointer(nodes[1],1,index);
  getf5pointer(nodes[2],2,index);
  getf5pointer(nodes[3],3,index);
  pointerCopy(0, 4);
  pointerCopy(1, 5);
  pointerCopy(2, 6);
  pointerCopy(3, 7);
  computeHistory(eletag, time, U, Ud, Udd, updateDm1); 
}

void PlaneDRMInputHandler::handle_elementAtface1(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd) 
{
  bool updateDm1 =  false;
  Node** nodes = eletag->getNodePtrs();
  int index = getIndex(time);
  //  int lastIndex = eletag->getLastDRMIndex();
  //  if (index != lastIndex) {
  //    updateDm1 = true;
  //    eletag->setLastDRMIndex(index);
  //  }

  int tag = eletag->getTag();
  int lastIndex = ele_str2[tag];
  if (index != lastIndex) {
    updateDm1 = true;
    ele_str2[tag] = index;
  }

  getf1pointer(nodes[1],1,index);
  getf1pointer(nodes[2],2,index);
  getf1pointer(nodes[5],5,index);
  getf1pointer(nodes[6],6,index);
  pointerCopy(1, 0);
  pointerCopy(2, 3);
  pointerCopy(5, 4);
  pointerCopy(6, 7);
  computeHistory(eletag, time, U, Ud, Udd, updateDm1); 
}

void PlaneDRMInputHandler::handle_elementAtface2(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd)
{
  bool updateDm1 =  false;
  Node** nodes = eletag->getNodePtrs();
  int index = getIndex(time);	
//   int lastIndex = eletag->getLastDRMIndex();
//   if (index != lastIndex) {
//     updateDm1 = true;
//     eletag->setLastDRMIndex(index);
//   }

  int tag = eletag->getTag();
  int lastIndex = ele_str2[tag];
  if (index != lastIndex) {
    updateDm1 = true;
    ele_str2[tag] = index;
  }

  getf2pointer(nodes[0],0,index);
  getf2pointer(nodes[3],3,index);
  getf2pointer(nodes[4],4,index);
  getf2pointer(nodes[7],7,index);
  pointerCopy(0, 1);
  pointerCopy(3, 2);
  pointerCopy(4, 5);
  pointerCopy(7, 6);
  computeHistory(eletag, time, U, Ud, Udd, updateDm1); 
}

void PlaneDRMInputHandler::handle_elementAtface3(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd)
{
  bool updateDm1 =  false;
  Node** nodes = eletag->getNodePtrs();
  int index = getIndex(time);
//   int lastIndex = eletag->getLastDRMIndex();
//   if (index != lastIndex) {
//     updateDm1 = true;
//     eletag->setLastDRMIndex(index);
//   }

  int tag = eletag->getTag();
  int lastIndex = ele_str2[tag];
  if (index != lastIndex) {
    updateDm1 = true;
    ele_str2[tag] = index;
  }

  getf3pointer(nodes[0],0,index);
  getf3pointer(nodes[1],1,index);
  getf3pointer(nodes[4],4,index);
  getf3pointer(nodes[5],5,index);
  pointerCopy(0, 3);
  pointerCopy(1, 2);
  pointerCopy(4, 7);
  pointerCopy(5, 6);
  computeHistory(eletag, time, U, Ud, Udd, updateDm1); 
}

void PlaneDRMInputHandler::handle_elementAtface4(Element* eletag, double time, Vector& U, Vector& Ud, Vector& Udd)
{
  bool updateDm1 =  false;
  Node** nodes = eletag->getNodePtrs();
  int index = getIndex(time);
//   int lastIndex = eletag->getLastDRMIndex();
//   if (index != lastIndex) {
//     updateDm1 = true;
//     eletag->setLastDRMIndex(index);
//   }

  int tag = eletag->getTag();
  int lastIndex = ele_str2[tag];
  if (tag ==1)
    opserr << " Index " << index << endln;
  if (index != lastIndex) {
    updateDm1 = true;
    ele_str2[tag] = index;
  }

  getf4pointer(nodes[3],3,index);
  getf4pointer(nodes[2],2,index);
  getf4pointer(nodes[7],7,index);
  getf4pointer(nodes[6],6,index);
  pointerCopy(3, 0);
  pointerCopy(2, 1);
  pointerCopy(7, 4);
  pointerCopy(6, 5);
  computeHistory(eletag, time, U, Ud, Udd, updateDm1); 
}


void PlaneDRMInputHandler::getLocations(double x, double y, double dx, double dy, int* xloc, int* yloc)
{
  *xloc = floor(x/dx);
  *yloc = floor(y/dy);
}

void PlaneDRMInputHandler::getTemporal(double time, int* tloc)
{
  *tloc = floor(time/deltaT);
}

int PlaneDRMInputHandler::getIndex(double time)
{
  // kill roundoffs
  double tmptime = time + 0.0000000001;

  if (tmptime >= timeBuf[this->cacheValue+2]) {
    populateBuffers();
  }

  //now we are in memory for sure so get location in local buffers:
  int index;
  for (index=1; index<this->cacheValue+1; index++) {
    double tmp =  timeBuf[index];
    if ( (tmptime >= tmp) && (tmptime <= tmp+deltaT) )
      break;
  }
  if ( (index < 1) || (index > cacheValue+1)) {
    opserr << " Severe error aborting tasks index is: "<< index << " gc is: " << globalCounter<< "\n";
    return -1;
  }
  return index-1;
}
			
void PlaneDRMInputHandler::getf5pointer(Node* node_tag, int local_tag, int index)
{
  const Vector& crd = node_tag->getCrds();
  double x = crd(0);
  double y = crd(1);
//  double z = crd(2);
  double dx = this->eleD[0];
  double dy = this->eleD[1];
//  double dz = this->eleD[2];
  int xloc,yloc;
  int numx,numy;
  numx = this->fileData[13];
  numy = this->fileData[14];
  int temp = fileData[12];
  //  the num_ele below should be coming from eleData cause they refer to the cmu mesh
  x = -x +numx*dx; //change of crds for cmu crd system
  y = -y +numy*dy;
  // in there goes cmu ele data
  getLocations(x,y, dx,dy, &xloc, &yloc);

  //now we have the index in the buffers populate the vectors accordingly:
  index *= 3*temp;
  index += 3*(numx+1)*yloc;
  index += 3*xloc;
  bool no_interpolation = false;

//  int tttag = node_tag->getTag();

  if ( (xloc*dx == x) && ( yloc*dy == y) )
    no_interpolation = true;

  if (no_interpolation) {
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = f5buffer[index + i];
      Vtm1(3*local_tag+i) = f5buffer[index + i + 3*temp];
      Vtp1(3*local_tag+i) = f5buffer[index + i + 6*temp];
      Vtp2(3*local_tag+i) = f5buffer[index + i + 9*temp];
    }
  }
  else {
    double ksi = -1.0+2.0*(x - xloc*dx)/dx;
    double eta = -1.0+2.0*(y - yloc*dy)/dy;
    populateTempBuffers(index,5,ksi,eta);
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = Vtempm2[i];
      Vtm1(3*local_tag+i) = Vtempm1[i];
      Vtp1(3*local_tag+i) = Vtempp1[i];
      Vtp2(3*local_tag+i) = Vtempp2[i];
    }
  }
}
			
void PlaneDRMInputHandler::getf1pointer(Node* node_tag, int local_tag, int index)
{
//  int tag = node_tag->getTag();
  const Vector& crd = node_tag->getCrds();
//  double x = crd(0);
  double y = crd(1);
  double z = crd(2);
 // double dx = this->eleD[0];
  double dy = this->eleD[1];
  double dz = this->eleD[2];
  int yloc,zloc;
  int numy,numz;
  numy = this->fileData[1];
  numz = this->fileData[2];
  int temp = fileData[0];
  z = -z +numz*dz; //change of crds for cmu crd system
  getLocations(y,z, dy,dz, &yloc, &zloc);

  //now we have the index in the buffers populate the vectors accordingly:
  index *= 3*temp;
  index += 3*(numy+1)*zloc;
  index += 3*yloc;
  

  bool no_interpolation = false;

  if ( (zloc*dz == z) && ( yloc*dy == y) )
    no_interpolation = true;

  if (no_interpolation) {
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = f1buffer[index + i];
      Vtm1(3*local_tag+i) = f1buffer[index + i + 3*temp];
      Vtp1(3*local_tag+i) = f1buffer[index + i + 6*temp];
      Vtp2(3*local_tag+i) = f1buffer[index + i + 9*temp];
    }
  }
  else {
    double ksi = -1.0+2.0*(y -yloc*dy)/dy;
    double eta = -1.0+ 2.0*(z-zloc*dz)/dz;
    populateTempBuffers(index,1,ksi,eta);
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = Vtempm2[i];
      Vtm1(3*local_tag+i) = Vtempm1[i];
      Vtp1(3*local_tag+i) = Vtempp1[i];
      Vtp2(3*local_tag+i) = Vtempp2[i];
    }
  }
}

void PlaneDRMInputHandler::getf2pointer(Node* node_tag, int local_tag, int index)
{
  const Vector& crd = node_tag->getCrds();
//  double x = crd(0);
  double y = crd(1);
  double z = crd(2);
  //double dx = this->eleD[0];
  double dy = this->eleD[1];
  double dz = this->eleD[2];
  int yloc,zloc;
  int numy,numz;
  numy = this->fileData[4];
  numz = this->fileData[5];
  int temp = fileData[3];
  z = -z +numz*dz; //change of crds for cmu crd system
  getLocations(y,z, dy,dz, &yloc, &zloc);
		
  //now we have the index in the buffers populate the vectors accordingly:
  index *= 3*temp;
  index += 3*(numy+1)*zloc;
  index += 3*yloc;
  
  bool no_interpolation = false;
  
  if ( (zloc*dz == z) && ( yloc*dy == y) )
    no_interpolation = true;

  if (no_interpolation) {
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = f2buffer[index + i];
      Vtm1(3*local_tag+i) = f2buffer[index + i + 3*temp];
      Vtp1(3*local_tag+i) = f2buffer[index + i + 6*temp];
      Vtp2(3*local_tag+i) = f2buffer[index + i + 9*temp];
    }
  }
  else {
    double ksi = -1.0+2.0*(y -yloc*dy)/dy;
    double eta = -1.0+2.0*(z-zloc*dz)/dz;
    populateTempBuffers(index, 2,ksi,eta);
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = Vtempm2[i];
      Vtm1(3*local_tag+i) = Vtempm1[i];
      Vtp1(3*local_tag+i) = Vtempp1[i];
      Vtp2(3*local_tag+i) = Vtempp2[i];
    }
  }
}

void PlaneDRMInputHandler::getf3pointer(Node* node_tag, int local_tag, int index)
{
  const Vector& crd = node_tag->getCrds();
  double x = crd(0);
//  double y = crd(1);
  double z = crd(2);
  double dx = this->eleD[0];
//  double dy = this->eleD[1];
  double dz = this->eleD[2];
  int xloc,zloc;
  int numx,numz;
  numx = this->fileData[7];
  numz = this->fileData[8];
  int temp = fileData[6];
  x = -x +numx*dx; //change of crds for cmu crd system
  z = -z +numz*dz; //change of crds for cmu crd system
  getLocations(x,z, dx,dz, &xloc, &zloc);
    
  //now we have the index in the buffers populate the vectors accordingly:
  index *= 3*temp;
  index += 3*(numx+1)*zloc;
  index += 3*xloc;
  
  bool no_interpolation = false;

  if ( (xloc*dx == x) && ( zloc*dz == z) )
    no_interpolation = true;

  if (no_interpolation) {
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = f3buffer[index + i];
      Vtm1(3*local_tag+i) = f3buffer[index + i + 3*temp];
      Vtp1(3*local_tag+i) = f3buffer[index + i + 6*temp];
      Vtp2(3*local_tag+i) = f3buffer[index + i + 9*temp];
    }
  }
  else {
    double ksi = -1.0+2.0*(x -xloc*dx)/dx;
    double eta = -1.0+2.0*(z-zloc*dz)/dz;
    populateTempBuffers(index, 3,ksi,eta);
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = Vtempm2[i];
      Vtm1(3*local_tag+i) = Vtempm1[i];
      Vtp1(3*local_tag+i) = Vtempp1[i];
      Vtp2(3*local_tag+i) = Vtempp2[i];
    }
  }
}

void PlaneDRMInputHandler::getf4pointer(Node* node_tag, int local_tag, int index)
{
//  bool debug = false;
  const Vector& crd = node_tag->getCrds();
  double x = crd(0);
  //double y = crd(1);
  double z = crd(2);
  double dx = this->eleD[0];
//  double dy = this->eleD[1];
  double dz = this->eleD[2];
  int xloc,zloc;
  int numx,numz;
  numx = this->fileData[10];
  numz = this->fileData[11];
  int temp = fileData[9];
  x = -x +numx*dx; //change of crds for cmu crd system
  z = -z +numz*dz; //change of crds for cmu crd system
  getLocations(x,z, dx,dz, &xloc, &zloc);
    
  //now we have the index in the buffers populate the vectors accordingly:
  index *= 3*temp;
  index += 3*(numx+1)*zloc;
  index += 3*xloc;
  bool no_interpolation = false;

  if ( (xloc*dx == x) && ( zloc*dz == z) )
    no_interpolation = true;

  if (no_interpolation) {
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = f4buffer[index + i];
      Vtm1(3*local_tag+i) = f4buffer[index + i + 3*temp];
      Vtp1(3*local_tag+i) = f4buffer[index + i + 6*temp];
      Vtp2(3*local_tag+i) = f4buffer[index + i + 9*temp];
    }
  }
  else {
    double ksi = -1.0+2.0*(x -xloc*dx)/dx;
    double eta = -1.0+2.0*(z-zloc*dz)/dz;
    populateTempBuffers(index,4,ksi,eta);
    for (int i=0; i<3; i++) {
      Vtm2(3*local_tag+i) = Vtempm2[i];
      Vtm1(3*local_tag+i) = Vtempm1[i];
      Vtp1(3*local_tag+i) = Vtempp1[i];
      Vtp2(3*local_tag+i) = Vtempp2[i];
    }
  }
}
			
void PlaneDRMInputHandler::pointerCopy(int node_from, int node_to)
{
  int indf = 3*node_from;
  int indt = 3*node_to;
  for (int i=0; i<3; i++) {
    Vtm2(indt+i) = Vtm2(indf+i);
    Vtm1(indt+i) = Vtm1(indf+i);
    Vtp1(indt+i) = Vtp1(indf+i);
    Vtp2(indt+i) = Vtp2(indf+i);
  }
}

void PlaneDRMInputHandler::populateTempBuffers(int index, int fileptr, double ksi, double eta)
{
  //num nodes in file
  int temp = fileData[(fileptr-1)*3];
  //num nodes along horizontal traversal of file
  int fileD = fileData[(fileptr-1)*3+1];
  fileD++;
  fileD*= 3;
  Vtm2_n1(0) = buffers[fileptr-1][index];
  Vtm2_n1(1) = buffers[fileptr-1][index+1];
  Vtm2_n1(2) = buffers[fileptr-1][index+2];
  
  Vtm2_n2(0) = buffers[fileptr-1][index+3];
  Vtm2_n2(1) = buffers[fileptr-1][index+3+1];
  Vtm2_n2(2) = buffers[fileptr-1][index+3+2];
  
  Vtm2_n3(0) = buffers[fileptr-1][index+3+fileD];
  Vtm2_n3(1) = buffers[fileptr-1][index+3+1+fileD];
  Vtm2_n3(2) = buffers[fileptr-1][index+3+2+fileD];
  
  Vtm2_n4(0) = buffers[fileptr-1][index+fileD];
  Vtm2_n4(1) = buffers[fileptr-1][index+1+fileD];
  Vtm2_n4(2) = buffers[fileptr-1][index+2+fileD];
  
  Vtm1_n1(0) = buffers[fileptr-1][index+3*temp];
  Vtm1_n1(1) = buffers[fileptr-1][index+1+3*temp];
  Vtm1_n1(2) = buffers[fileptr-1][index+2+3*temp];
  
  Vtm1_n2(0) = buffers[fileptr-1][index+3+3*temp];
  Vtm1_n2(1) = buffers[fileptr-1][index+3+1+3*temp];
  Vtm1_n2(2) = buffers[fileptr-1][index+3+2+3*temp];

  Vtm1_n3(0) = buffers[fileptr-1][index+3+fileD+3*temp];
  Vtm1_n3(1) = buffers[fileptr-1][index+3+1+fileD+3*temp];
  Vtm1_n3(2) = buffers[fileptr-1][index+3+2+fileD+3*temp];
  
  Vtm1_n4(0) = buffers[fileptr-1][index+fileD+3*temp];
  Vtm1_n4(1) = buffers[fileptr-1][index+1+fileD+3*temp];
  Vtm1_n4(2) = buffers[fileptr-1][index+2+fileD+3*temp];
  
  Vtp1_n1(0) = buffers[fileptr-1][index+6*temp];
  Vtp1_n1(1) = buffers[fileptr-1][index+1+6*temp];
  Vtp1_n1(2) = buffers[fileptr-1][index+2+6*temp];
  
  Vtp1_n2(0) = buffers[fileptr-1][index+3+6*temp];
  Vtp1_n2(1) = buffers[fileptr-1][index+3+1+6*temp];
  Vtp1_n2(2) = buffers[fileptr-1][index+3+2+6*temp];
  
  Vtp1_n3(0) = buffers[fileptr-1][index+3+fileD+6*temp];
  Vtp1_n3(1) = buffers[fileptr-1][index+3+1+fileD+6*temp];
  Vtp1_n3(2) = buffers[fileptr-1][index+3+2+fileD+6*temp];
  
  Vtp1_n4(0) = buffers[fileptr-1][index+fileD+6*temp];
  Vtp1_n4(1) = buffers[fileptr-1][index+1+fileD+6*temp];
  Vtp1_n4(2) = buffers[fileptr-1][index+2+fileD+6*temp];

  Vtp2_n1(0) = buffers[fileptr-1][index+9*temp];
  Vtp2_n1(1) = buffers[fileptr-1][index+1+9*temp];
  Vtp2_n1(2) = buffers[fileptr-1][index+2+9*temp];
  
  Vtp2_n2(0) = buffers[fileptr-1][index+3+9*temp];
  Vtp2_n2(1) = buffers[fileptr-1][index+3+1+9*temp];
  Vtp2_n2(2) = buffers[fileptr-1][index+3+2+9*temp];
  
  Vtp2_n3(0) = buffers[fileptr-1][index+3+fileD+9*temp];
  Vtp2_n3(1) = buffers[fileptr-1][index+3+1+fileD+9*temp];
  Vtp2_n3(2) = buffers[fileptr-1][index+3+2+fileD+9*temp];
  
  Vtp2_n4(0) = buffers[fileptr-1][index+fileD+9*temp];
  Vtp2_n4(1) = buffers[fileptr-1][index+1+fileD+9*temp];
  Vtp2_n4(2) = buffers[fileptr-1][index+2+fileD+9*temp];

  //finally perform spatial interpolation @ the node for the four time steps
  Vtempm2 = 0.25*(1-ksi)*(1-eta)*Vtm2_n1 + 0.25*(1+ksi)*(1-eta)*Vtm2_n2 + 0.25*(1+ksi)*(1+eta)*Vtm2_n3 + 0.25*(1-ksi)*(1+eta)*Vtm2_n4;
  Vtempm1 = 0.25*(1-ksi)*(1-eta)*Vtm1_n1 + 0.25*(1+ksi)*(1-eta)*Vtm1_n2 + 0.25*(1+ksi)*(1+eta)*Vtm1_n3 + 0.25*(1-ksi)*(1+eta)*Vtm1_n4;
  Vtempp1 = 0.25*(1-ksi)*(1-eta)*Vtp1_n1 + 0.25*(1+ksi)*(1-eta)*Vtp1_n2 + 0.25*(1+ksi)*(1+eta)*Vtp1_n3 + 0.25*(1-ksi)*(1+eta)*Vtp1_n4;
  Vtempp2 = 0.25*(1-ksi)*(1-eta)*Vtp2_n1 + 0.25*(1+ksi)*(1-eta)*Vtp2_n2 + 0.25*(1+ksi)*(1+eta)*Vtp2_n3 + 0.25*(1-ksi)*(1+eta)*Vtp2_n4;

}
