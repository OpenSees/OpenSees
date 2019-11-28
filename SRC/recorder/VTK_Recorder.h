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
                                                                        
// $Revision: 1.0 $
// $Date: 2015-11-12 $
                                                                        
#ifndef VTK_Recorder_h
#define VTK_Recorder_h

// Written: fmk
//
// Description: This file contains the class definition for 
// VTK_Recorder. A VTK_Recorder is used to store all responses in pvd format.


#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <ID.h>
#include <Recorder.h>

class Node;
class Element;

class OutputData {

 public:
  OutputData();
  OutputData &operator=(const OutputData &other);

  bool disp,
    disp2,
    disp3,
    vel,
    vel2,
    vel3,
    accel,
    accel2,
    accel3,
    reaction,
    reaction2,
    reaction3,
    mass,
    unbalancedLoad;
  int modes[10];
};


class VTK_Recorder: public Recorder
{
public:
  OutputData outputData;
  typedef std::vector<std::string> EleData;
    
public:
  VTK_Recorder(const char *filename, const OutputData& ndata,
	       const std::vector<EleData>& edata, int ind=2, int pre=10, double dt=0);
  VTK_Recorder();
  ~VTK_Recorder();
  
  int record(int commitTag, double timeStamp);
  int restart();
  int domainChanged();    
  int setDomain(Domain &theDomain);
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
 protected:
  virtual int initialize();
  bool initDone;
  
  virtual int vtu();
  virtual void addEleData(const EleData& edata) {eledata.push_back(edata);}
  std::vector<EleData> eledata;
  
 private:
  virtual void indent();
  virtual void incrLevel() {indentlevel++;}
  virtual void decrLevel() {indentlevel--;}
  void getfilename(const char* name);
  
  
 private:
  int indentsize, precision, indentlevel;
  char quota;    
  Domain *theDomain;
  bool echoTimeFlag;   
  int dataFlag;        
  double nextTimeStampToRecord;
  double deltaT;
  
  char *name;
  int counter;
  
  int ndm; // max ndm of nodes, 3 min for 3d viewing
  int ndf; // max ndf of all nodes
  
  std::ofstream thePVDFile;
  std::ofstream theVTUFile;
  
  std::map<int,int>theNodeMapping; // output requires points indexed at 0
  std::map<int,int>theEleMapping; // output requires points indexed at 0
  
  std::vector<int>theNodeTags;
  std::vector<int>theEleTags;
  std::vector<int>theEleClassTags;
  std::vector<int>theEleVtkTags;
  std::vector<int>theEleVtkOffsets;
  
 public:
  enum VtkType {
    VTK_VERTEX=1,VTK_POLY_VERTEX=2,VTK_LINE=3,VTK_POLY_LINE=4,
    VTK_TRIANGLE=5,VTK_TRIANGLE_STRIP=6,VTK_POLYGON=7,
    VTK_PIXEL=8,VTK_QUAD=9,VTK_TETRA=10,VTK_VOXEL=11,
    VTK_HEXAHEDRON=12,VTK_WEDGE=12,VTK_PYRAMID=14,
    VTK_QUADRATIC_EDGE=21,VTK_QUADRATIC_TRIANGLE=22,
    VTK_QUADRATIC_QUAD=23,VTK_QUADRATIC_TETRA=24,
    VTK_QUADRATIC_HEXAHEDRON=25
  };
  static void setVTKType();
  
 private:
  static std::map<int,VtkType> vtktypes;
  int numNode;
  int numElement;
  int maxNDM;
  int maxNDF;
};

#endif
