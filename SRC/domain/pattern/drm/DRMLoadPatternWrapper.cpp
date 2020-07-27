/*
 *
 * @author: gnp <petropoulos@gmail.com>
 *
 * @Description: DRM Load Pattern Wrapper Implementation
 *               This class exists in order to allow 
 *               for easy integration in the (*****) opensees software
 *               of the drm load pattern related classes.
 *               Essentially creates the required classes 
 *               and implements the send/recv Self methods.
 *           
 * @Date: 7/7/08
 *
 */
#include "DRMLoadPatternWrapper.h"
#include "PlaneDRMInputHandler.h"

DRMLoadPatternWrapper::DRMLoadPatternWrapper(int ttag, double cfact, char** on_files, int sfiles,  double ddt,
					     int numsteps, int* filedata, int fileDatasize,
					     int nd_1, int nd_2,
					     double* drmboxcrds, double* ele_D,
					     int stepstocache)
  :LoadPattern(ttag, PATTERN_TAG_DRMLoadPattern)
{
  this->infiles = on_files;
  this->fileData_size = fileDatasize;
  this->file_data = filedata;
  this->factor = cfact;
  this->files = sfiles;
  this->dt = ddt;
  this->num_steps = numsteps;  
  this->nd1 = nd_1;
  this->nd2 = nd_2;
  this->drm_box_crds = drmboxcrds;
  this->eleD = ele_D;
  this->steps_to_cache = stepstocache;
  
  this->myPattern = 0;
  this->initialized = false;
  this->cleanUpAfterMySelf = false;
}

DRMLoadPatternWrapper::DRMLoadPatternWrapper()
  :LoadPattern(0,PATTERN_TAG_DRMLoadPattern)
{
  this->infiles = 0;
  this->fileData_size = 0;
  this->file_data = 0;
  this->factor = 0;
  this->files = 0;
  this->dt = 0;
  this->num_steps = 0;
  this->nd1 = 0;
  this->nd2 = 0;
  this->drm_box_crds = 0;
  this->eleD = 0;
  this->steps_to_cache = 0;
  
  this->myPattern = 0;
  this->initialized = false;
  this->cleanUpAfterMySelf = false;
}



DRMLoadPatternWrapper::~DRMLoadPatternWrapper()
{
  //  nothing to deallocate here.
  if (this->cleanUpAfterMySelf) {
    delete [] this->eleD;
    delete [] this->drm_box_crds;
    delete [] this->file_data;
    for (int i=0; i<6; i++)
      delete [] infiles[i];
    delete [] infiles;
  }

}

void DRMLoadPatternWrapper::initialize() 
{
  Mesh3DSubdomain* my_mesher = new Mesh3DSubdomain(this->getDomain());
  
  PlaneDRMInputHandler* patternhandler = new PlaneDRMInputHandler(1, infiles, files,  dt, 0,
 								  num_steps, file_data, fileData_size,
 								  nd1, nd2,
 								  drm_box_crds, drm_box_crds, eleD,
 								  my_mesher, steps_to_cache, this->getDomain());
  
  DRMLoadPattern* ptr = new DRMLoadPattern(1,factor,patternhandler,this->getDomain());
  
  this->myPattern = ptr;
  initialized = true;
}

void DRMLoadPatternWrapper::applyLoad(double time) 
{
  if (!initialized) 
    this->initialize();
  
  this->myPattern->applyLoad(time); 

}

int DRMLoadPatternWrapper::sendSelf(int commitTag, Channel& theChannel)
{
  int dbTag = this->getDbTag();
  
  static ID i_Data(22);
  
  i_Data(0) = this->getTag();
  i_Data(1) = this->fileData_size;
  for (int i=0; i<this->fileData_size; i++)
    i_Data(2+i) = file_data[i];
  i_Data(17) = this->files;
  i_Data(18) = this->num_steps;
  i_Data(19) = this->nd1;
  i_Data(20) = this->nd2;
  i_Data(21) = this->steps_to_cache;
 
  if ( theChannel.sendID(dbTag, commitTag, i_Data) < 0 ) {
    opserr << "DRMLoadPatternWrapper::sendSelf L.121 failed to sendID \n";
    return -1;
  }

  static Vector d_Data(11);
  
  d_Data(0) = this->dt;
  for (int i=0; i<6; i++)
    d_Data(1+i) = this->drm_box_crds[i];
  for (int i=0; i<3; i++)
    d_Data(7+i) = this->eleD[i];
  d_Data(10) = this->factor;

  if ( theChannel.sendVector(dbTag, commitTag, d_Data) < 0 ) {
    opserr << "DRMLoadPatternWrapper::sendSelf L.135 failed to sendVector \n";
    return -1;
  } 
  
  static ID c_Data_sz(this->files+1);
//  int pos =0;
  std::string final_str;
  for (int i=0; i<this->files; i++) {
    char* tmp = this->infiles[i];
    std::string str(tmp);
    int sz = str.size();
    c_Data_sz(i) = sz;
    final_str.append(str);
  }
  
  int ssz = final_str.size()+1;
  c_Data_sz(this->files) = ssz;

  Message c_Data(((char*) final_str.c_str()), ssz);
  
  if ( theChannel.sendID(dbTag, commitTag, c_Data_sz) < 0) {
    opserr << "DRMLoadPatternWrapper::sendSelf L.156 failed to sendID2 \n";
    return -1;
  }

  if (theChannel.sendMsg(dbTag, commitTag, c_Data) ) {
    opserr << "DRMLoadPatternWrapper::sendSelf L.161 failed to sendMsg \n";
    return -1;
  }

  return 0;

}

int DRMLoadPatternWrapper::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{

  int dbTag = this->getDbTag();
  
  static ID i_Data(22);
  if (theChannel.recvID(dbTag, commitTag, i_Data) < 0) {
    opserr << "DRMLoadPatternWrapper::recvSelf L.176 failed to recvID \n";
    return -1;
  } 
  

  this->setTag(i_Data(0));

  this->fileData_size = i_Data(1);
  
  this->file_data = new int[this->fileData_size];
  for (int i=0; i<this->fileData_size; i++)
    file_data[i] = i_Data(2+i);
  this->files = i_Data(17);
  this->num_steps = i_Data(18);
  this->nd1 = i_Data(19);
  this->nd2 = i_Data(20);
  this->steps_to_cache = i_Data(21);
  
  this->myPattern = 0;
  this->initialized = false;
  
  static Vector d_Data(11);
  
  if (theChannel.recvVector(dbTag, commitTag, d_Data) < 0) {
    opserr << "DRMLoadPatternWrapper::recvSelf L.200 failed to recvVector \n";
    return -1;
  } 


  this->dt = d_Data(0);
  this->drm_box_crds = new double[6];
  for (int i=0; i<6; i++)
    this->drm_box_crds[i] = d_Data(1+i);
  this->eleD = new double[3];
  for (int i=0; i<3; i++)
    this->eleD[i] = d_Data(7+i);

  this->factor = d_Data(10);

  static ID c_Data_sz(this->files+1);
  if (theChannel.recvID(dbTag, commitTag, c_Data_sz) < 0) {
    opserr << "DRMLoadPatternWrapper::recvSelf L.217 failed to recvID2 \n";
    return -1;
  } 


  int ssz = c_Data_sz(this->files);
  char *stor = new char[ssz];
  Message c_Data(stor, ssz);

  if (theChannel.recvMsg(dbTag, commitTag, c_Data) < 0) {
    opserr << "DRMLoadPatternWrapper::recvSelf L.227 failed to recvMsg \n";
    return -1;
  } 
  std::string final_str(c_Data.getData());
  
  int pos =0;

  this->infiles = new char*[this->files];

  for (int i=0; i<this->files; i++) {    
    int tmp_sz = c_Data_sz(i);
    std::string tmp = final_str.substr(pos, tmp_sz);
    char* tmpc = new char[tmp.size()+1];
    strcpy(tmpc,tmp.c_str());
    this->infiles[i] = tmpc;
    pos += tmp_sz;
  }
  
  this->cleanUpAfterMySelf = true;
  delete [] stor;
  return 0;
}

LoadPattern* DRMLoadPatternWrapper::getCopy() {
  return new DRMLoadPatternWrapper(this->getTag(), this->factor, this->infiles, this->files, this->dt,
				   this->num_steps, this->file_data, this->fileData_size,
				   this->nd1, this->nd2,
				   this->drm_box_crds, this->eleD,
				   this->steps_to_cache);
}


void DRMLoadPatternWrapper::Print(OPS_Stream &s, int flag) {
  opserr << "DRMLoadPattern::Print() - not yet implemented\n";    
  for (int i=0; i<this->files; i++)
    printf("infiles @i %d %s \n",i,this->infiles[i]);
  printf(" %d %f %d %f %d %d %d %d \n",fileData_size, factor,files,dt,num_steps,nd1,nd2,steps_to_cache);
  for (int i=0; i<6; i++)
    printf("drmbiox @i %d %f  \n",i,this->drm_box_crds[i]);
  for (int i=0; i<3; i++)
    printf("eleD @i %d %f \n",i,this->eleD[i]);
}
