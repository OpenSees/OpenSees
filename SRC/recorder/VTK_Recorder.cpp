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
// Save all data into paraview format

#include "VTK_Recorder.h"
#include <sstream>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Message.h>

#include <classTags.h>
#include <iostream>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif


OutputData::OutputData()
{
    disp = false;
    disp2 = false;
    disp3 = false;
    vel = false;
    vel2 = false;
    vel3 = false;
    accel = false;
    accel2 = false;
    accel3 = false;
    reaction = false;
    reaction2 = false;
    reaction3 = false;
    mass = false;
    unbalancedLoad = false;

    for (int i=0; i<10; i++) {
      modes[i] = 0;
    }
}
OutputData &
OutputData::operator=(const OutputData &other) 
{
  // first check we are not trying v = v
  if (this != &other) {
    disp = other.disp;
    disp2 = other.disp2;
    disp3 = other.disp3;
    vel = other.vel;
    vel2 = other.vel2;
    vel3 = other.vel3;
    accel = other.accel;
    accel2 = other.accel2;
    accel3 = other.accel3;
    reaction = other.reaction;
    reaction2 = other.reaction2;
    reaction3 = other.reaction3;
    for (int i=0; i<10; i++) {
      modes[i] = 0;
    }
  }
  return *this;
}


std::map<int,VTK_Recorder::VtkType> VTK_Recorder::vtktypes;


void* OPS_VTK_Recorder()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if(numdata < 1) {
	opserr<<"WARNING: insufficient number of arguments\n";
	return 0;
    }

    // filename
    const char* name = OPS_GetString();

    // plotting options
    numdata = OPS_GetNumRemainingInputArgs();
    int indent=2;
    int precision = 10;
    OutputData outputData;
    std::vector<VTK_Recorder::EleData> eledata;
    double dT = 0.0;
    double rTolDt = 0.00001;

    while(numdata > 0) {
	const char* type = OPS_GetString();

	if(strcmp(type, "disp") == 0) {
	    outputData.disp = true;
	} else if(strcmp(type, "disp2") == 0) {
	  outputData.disp2 = true;
	} else if(strcmp(type, "disp3") == 0) {
	  outputData.disp3 = true;

	} else if(strcmp(type, "vel") == 0) {
	    outputData.vel = true;

	    /*
	} else if(strcmp(type, "vel2") == 0) {
	    outputData.vel2 = true;
	} else if(strcmp(type, "vel3") == 0) {
	    outputData.vel3 = true;
	    */

	} else if(strcmp(type, "accel") == 0) {
	    outputData.accel = true;
	    /*
	} else if(strcmp(type, "accel2") == 0) {
	    outputData.accel2 = true;
	} else if(strcmp(type, "accel3") == 0) {
	    outputData.accel3 = true;
	    */

	} else if(strcmp(type, "reaction") == 0) {
	    outputData.reaction = true;
	} else if(strcmp(type, "reaction2") == 0) {
	    outputData.reaction2 = true;
	} else if(strcmp(type, "reaction3") == 0) {
	    outputData.reaction3 = true;


	} else if(strcmp(type, "mass") == 0) {
	    outputData.mass = true;
	} else if(strcmp(type, "unbalancedLoad") == 0) {
	    outputData.unbalancedLoad = true;
	} else if(strcmp(type, "eigen") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 1) {
		opserr<<"WARNING: eigen needs 'numEigenvector'\n";
		return 0;
	    }
	    numdata = 1;
	    int mode;
	    if(OPS_GetIntInput(&numdata,&mode) < 0) {
		opserr << "WARNING: failed to read numeigen\n";
		return 0;
	    }
	    for (int i=0; i<10; i++)
	      if (outputData.modes[i] == 0)
		outputData.modes[i] = mode;

	} else if(strcmp(type, "-precision") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 1) {
		opserr<<"WARNING: needs precision \n";
		return 0;
	    }
	    numdata = 1;
	    if(OPS_GetIntInput(&numdata,&precision) < 0) {
		opserr << "WARNING: failed to read precision\n";
		return 0;
	    }

	} else if(strcmp(type, "eleResponse") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 1) {
		opserr<<"WANRING: elementResponse needs 'argc','argv'\n";
		return 0;
	    }
	    VTK_Recorder::EleData edata;
	    numdata = OPS_GetNumRemainingInputArgs();
	    edata.resize(numdata);
	    for(int i=0; i<numdata; i++) {
		edata[i] = OPS_GetString();
	    }
	    eledata.push_back(edata);
	} else if(strcmp(type, "-dT") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 1) {
		opserr<<"WARNING: needs dT \n";
		return 0;
	    }
	    numdata = 1;
	    if(OPS_GetDoubleInput(&numdata,&dT) < 0) {
		opserr << "WARNING: failed to read dT\n";
		return 0;
	    }
	    if (dT < 0) dT = 0;
	} else if(strcmp(type, "-rTolDt") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 1) {
		opserr<<"WARNING: needs rTolDt \n";
		return 0;
	    }
	    numdata = 1;
	    if(OPS_GetDoubleInput(&numdata,&rTolDt) < 0) {
		opserr << "WARNING: failed to read rTolDt\n";
		return 0;
	    }
	    if (rTolDt < 0) rTolDt = 0;
	}
	numdata = OPS_GetNumRemainingInputArgs();
    }

    // create recorder
    return new VTK_Recorder(name,outputData,eledata,indent,precision,dT, rTolDt);
}

VTK_Recorder::VTK_Recorder(const char *inputName, 
			   const OutputData& outData,
			   const std::vector<EleData>& edata, 
			   int ind, int pre, double dt, double rTolDt)
    :Recorder(RECORDER_TAGS_VTK_Recorder), 
     indentsize(ind), 
     precision(pre),
     indentlevel(0), 
     quota('\"'), 
     theDomain(0),
     nextTimeStampToRecord(0.0),
     deltaT(dt),
     relDeltaTTol(rTolDt),
     counter(0),
     initializationDone(false),
     sendSelfCount(0)
{
  outputData = outData;

  name = new char[strlen(inputName+1)];
  strcpy(name, inputName);

  //
  // mkdir of same name to store the vtu files
  //

#if defined(_WIN32)
  _mkdir(name);
#else 
  mkdir(name, 0777); // notice that 777 is different than 0777
#endif

  //
  // set all the possible VTK types for use in init
  //

  VTK_Recorder::setVTKType();

  initDone = false;

  //
  // open pvd file
  //

  char *filename = new char[strlen(name) + 5];
  sprintf(filename, "%s.pvd",name);
  
  thePVDFile.close();
  thePVDFile.open(filename, std::ios::out);

  if(thePVDFile.fail()) {
    opserr<<"WARNING: Failed to open vtd file "<< filename<< "\n";
  }
  thePVDFile.precision(precision);
  thePVDFile << std::scientific;
  
  //
  // spit out header to file
  //
  thePVDFile << "<?xml version="<<quota<<"1.0"<<quota<<"?>\n";
  thePVDFile << "<VTKFile type=\"Collection\" version=\"1.0\" \n";
  thePVDFile << "byte_order=\"LittleEndian\" \n";
  thePVDFile << "compressor=\"vtkZLibDataCompressor\">\n";
  
  thePVDFile<<"<Collection>\n";
}

VTK_Recorder::VTK_Recorder()
  :Recorder(RECORDER_TAGS_VTK_Recorder),
   indentsize(0), 
   precision(0),
   indentlevel(0), 
   quota('\"'), 
   theDomain(0),
   nextTimeStampToRecord(0.0),
   deltaT(0.0),
   relDeltaTTol(0.00001),
   counter(0),
   initializationDone(false),
   sendSelfCount(0)   
{
  name = NULL;

  //
  // set all the possible VTK types for use in init
  //

  VTK_Recorder::setVTKType();

  initDone = false;
}


VTK_Recorder::~VTK_Recorder()
{
  //
  // write out last bits and close the vtd file
  //

  thePVDFile << "</Collection>\n </VTKFile>\n";
  thePVDFile.close();
}

int
VTK_Recorder::record(int ctag, double timeStamp)
{
  if (initializationDone == false) {
    this->initialize();
    initializationDone = true;
  }

  // where relDeltaTTol is the maximum reliable ratio between analysis time step and deltaT
  // and provides tolerance for floating point precision (see floating-point-tolerance-for-recorder-time-step.md)
    if (deltaT == 0.0 || timeStamp - nextTimeStampToRecord >= -deltaT * relDeltaTTol) {
    
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;
  
    //
    // add a line to pvd file
    //

    char *filename = new char[2*strlen(name)+26];

    // process p0 writes the pvd file, part for each process including itself 0
    if (sendSelfCount >= 0) {
      for (int i=0; i<= sendSelfCount; i++) {
	sprintf(filename, "%s/%s%d%020d.vtu",name, name, i, counter);    
	thePVDFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"" << i 
		   <<"\"" << " file=\"" << filename << "\"/>\n";
      }
    }


    //
    // write vtu file
    //

    return vtu();
  }

  return 0;
}

int
VTK_Recorder::restart()
{
  return 0;
}

int
VTK_Recorder::domainChanged()
{
  this->initialize();
  return 0;
}

int
VTK_Recorder::setDomain(Domain& domain)
{
  theDomain = &domain;
  return 0;
}


int
VTK_Recorder::vtu()
{
  if (theDomain == 0) {
    opserr << "WARNING: failed to get domain -- VTK_Recorder::vtu\n";
    return -1;
  }
  
  char *filename = new char[2*strlen(name)+26];
  if (sendSelfCount < 0) {
    sprintf(filename, "%s/%s%d%020d.vtu",name, name, -sendSelfCount, counter);    
  } else {
    sprintf(filename, "%s/%s%d%020d.vtu",name, name, 0, counter);    
  }
  
  counter ++;
  
  std::ofstream theFileVTU;
  theFileVTU.open(filename, std::ios::out);
  
  if(theFileVTU.fail()) {
    opserr<<"WARNING: Failed to open file "<<filename<<"\n";
    return -1;
  }
  
  theFileVTU.precision(precision);
  theFileVTU << std::scientific;
  
  // header
  theFileVTU<<"<?xml version="<<quota<<"1.0"<<quota<<"?>\n";
  theFileVTU<<"<VTKFile type="<<quota<<"UnstructuredGrid"<<quota;
  theFileVTU<<" version="<<quota<<"1.0"<<quota;
  theFileVTU<<" byte_order="<<quota<<"LittleEndian"<<quota;
  theFileVTU<<">\n";
  this->incrLevel();
  this->indent();
  theFileVTU<<"<UnstructuredGrid>\n";
  
  // Piece
  this->incrLevel();
  this->indent();
  theFileVTU<<"<Piece NumberOfPoints=\"" << numNode <<"\" NumberOfCells=\"" << numElement << "\">\n";

  //
  // POINT DATA
  // 

  theFileVTU<<"<PointData>\n";
  this->incrLevel();

  // node tags
  theFileVTU<<"<DataArray type=\"Int64\" Name=\"Node Tag\" format=\"ascii\">\n";
  this->incrLevel();
  for (auto i : theNodeTags)
    theFileVTU << i << " ";
  theFileVTU<<"\n</DataArray>\n";


  // node displacements
  if (outputData.disp == true) {
    theFileVTU<<"<DataArray type=\"Float64\" Name=\"Disp\" NumberOfComponents=\"" << maxNDF << "\" format=\"ascii\">\n";
    for (auto i : theNodeTags) {
      Node *theNode=theDomain->getNode(i);
      const Vector &output=theNode->getDisp();
      int numDOF = output.Size();
      for (int i=0; i<numDOF; i++) 
	theFileVTU << output(i) << " ";
      for (int i=numDOF; i<maxNDF; i++)
	theFileVTU << 0.0 << " ";
      theFileVTU << "\n";
    }
    theFileVTU<<"\n</DataArray>\n";
  }

  if (outputData.disp2 == true) {
    theFileVTU<<"<DataArray type=\"Float64\" Name=\"Disp2\" NumberOfComponents=\"" << 2 << "\" format=\"ascii\">\n";
    for (auto i : theNodeTags) {
      Node *theNode=theDomain->getNode(i);
      const Vector &output=theNode->getDisp();
      int numDOF = output.Size();
      for (int i=0; i<2; i++) 
	if (i < numDOF) 
	  theFileVTU << output(i) << " ";
	else
	  theFileVTU << 0.0 << " ";
      theFileVTU << "\n";
    }
    theFileVTU<<"\n</DataArray>\n";
  }

  if (outputData.disp3 == true) {
    theFileVTU<<"<DataArray type=\"Float64\" Name=\"Disp3\" NumberOfComponents=\"" << 3 << "\" format=\"ascii\">\n";
    for (auto i : theNodeTags) {
      Node *theNode=theDomain->getNode(i);
      const Vector &output=theNode->getDisp();
      int numDOF = output.Size();
      for (int i=0; i<3; i++) 
	if (i < numDOF) 
	  theFileVTU << output(i) << " ";
	else
	  theFileVTU << 0.0 << " ";
      theFileVTU << "\n";
    }
    theFileVTU<<"\n</DataArray>\n";
  }

  //
  // node vel
  //

  if (outputData.vel == true) {
    theFileVTU<<"<DataArray type=\"Float64\" Name=\"Vel\" NumberOfComponents=\"" << maxNDF << "\" format=\"ascii\">\n";
    for (auto i : theNodeTags) {
      Node *theNode=theDomain->getNode(i);
      const Vector &output=theNode->getVel();
      int numDOF = output.Size();
      for (int i=0; i<numDOF; i++) 
	theFileVTU << output(i) << " ";
      for (int i=numDOF; i<maxNDF; i++)
	theFileVTU << 0.0 << " ";
      theFileVTU << "\n";
    }
    theFileVTU<<"\n</DataArray>\n";
  }

  //
  // node accel
  //

  if (outputData.accel == true) {
    theFileVTU<<"<DataArray type=\"Float64\" Name=\"Accel\" NumberOfComponents=\"" << maxNDF << "\" format=\"ascii\">\n";
    for (auto i : theNodeTags) {
      Node *theNode=theDomain->getNode(i);
      const Vector &output=theNode->getAccel();
      int numDOF = output.Size();
      for (int i=0; i<numDOF; i++) 
	theFileVTU << output(i) << " ";
      for (int i=numDOF; i<maxNDF; i++)
	theFileVTU << 0.0 << " ";
      theFileVTU << "\n";
    }
    theFileVTU<<"\n</DataArray>\n";
  }

  // 
  // switch to element data
  //

  theFileVTU<<"</PointData>\n<CellData>\n";

  // ele tags
  theFileVTU<<"<DataArray type=\"Int64\" Name=\"Element Tag\" format=\"ascii\">\n";
  this->incrLevel();
  for (auto i : theEleTags)
    theFileVTU << i << " ";
  theFileVTU<<"\n</DataArray>\n";

  // ele class tags
  theFileVTU<<"<DataArray type=\"Int64\" Name=\"Element Class\" format=\"ascii\">\n";
  this->incrLevel();
  for (auto i : theEleClassTags)
    theFileVTU << i << " ";
  theFileVTU<<"\n</DataArray>\n";

  theFileVTU<<"</CellData>\n";

  //
  // points - output nodal coords
  //

  this->incrLevel();
  this->indent();
  theFileVTU<<"<Points>\n";
  theFileVTU<<"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (auto i : theNodeTags) {
    Node *theNode=theDomain->getNode(i);
    const Vector &crd=theNode->getCrds();
    int numCrd = crd.Size();
    for (int i=0; i<numCrd; i++) 
      theFileVTU << crd(i) << " ";
    for (int i=numCrd; i<3; i++)
      theFileVTU << 0.0 << " ";
    theFileVTU << "\n";
  }
  theFileVTU<<"</DataArray>\n";
  theFileVTU<<"</Points>\n";

  //
  // cells - output element connectivity, offsets and types
  //

  theFileVTU<<"<Cells>\n";

  // connectivity
  theFileVTU<<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (auto i : theEleTags) {
    Element *theEle=theDomain->getElement(i);
    if (theEle != 0) {
      const ID &theNodes=theEle->getExternalNodes();
      int numNode = theNodes.Size();
      for (int i=0; i<numNode; i++) {
	int nodeTag = theNodes(i);
	auto nodeID = theNodeMapping[nodeTag];
	theFileVTU << nodeID << " ";
      }
      theFileVTU << "\n";
    }
  }
  theFileVTU<<"</DataArray>\n";

  // offset
  theFileVTU<<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (auto i : theEleVtkOffsets)
    theFileVTU << i << " ";
  theFileVTU<<"\n</DataArray>\n";

  // types
  theFileVTU<<"<DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
  for (auto i : theEleVtkTags)
    theFileVTU << i << " ";
  theFileVTU<<"\n</DataArray>\n";

  theFileVTU<<"</Cells>\n";

  this->indent();
  /*
    theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
    theFileVTU<<" Name="<<quota<<"Points"<<quota;
    theFileVTU<<" NumberOfComponents="<<quota<<3<<quota;
    theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
    
    // points coordinates
    this->incrLevel();
    for(int i=0; i<(int)nodes.size(); i++) {
	const Vector& crds = nodes[i]->getCrds();
	this->indent();
	for(int j=0; j<3; j++) {
	    if(j < crds.Size()) {
		theFileVTU<<crds(j)<<' ';
	    } else {
		theFileVTU<<0.0<<' ';
	    }
	}
	theFileVTU<<std::endl;
    }
    */

    /*
    // node velocity
    if(nodedata.vel) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"Velocity"<<quota;
	theFileVTU<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getTrialVel();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFileVTU<<vel(j)<<' ';
		} else {
		    theFileVTU<<0.0<<' ';
		}
	    }
	    theFileVTU<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }

    // node displacement
    if(nodedata.disp) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"Displacement"<<quota;
	theFileVTU<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getTrialDisp();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFileVTU<<vel(j)<<' ';
		} else {
		    theFileVTU<<0.0<<' ';
		}
	    }
	    theFileVTU<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }

    // node incr displacement
    if(nodedata.incrdisp) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"IncrDisplacement"<<quota;
	theFileVTU<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getIncrDisp();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFileVTU<<vel(j)<<' ';
		} else {
		    theFileVTU<<0.0<<' ';
		}
	    }
	    theFileVTU<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }

    // node acceleration
    if(nodedata.accel) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"Acceleration"<<quota;
	theFileVTU<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getTrialAccel();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFileVTU<<vel(j)<<' ';
		} else {
		    theFileVTU<<0.0<<' ';
		}
	    }
	    theFileVTU<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }

    // node pressure
    if(nodedata.pressure) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"Pressure"<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    double pressure = 0.0;
	    Pressure_Constraint* thePC = theDomain->getPressure_Constraint(nodes[i]->getTag());
	    if(thePC != 0) {
		pressure = thePC->getPressure();
	    }
	    this->indent();
	    theFileVTU<<pressure<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }

    // node reaction
    if(nodedata.reaction) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"Reaction"<<quota;
	theFileVTU<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getReaction();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFileVTU<<vel(j)<<' ';
		} else {
		    theFileVTU<<0.0<<' ';
		}
	    }
	    theFileVTU<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }

    // node unbalanced load
    if(nodedata.unbalanced) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"UnbalancedLoad"<<quota;
	theFileVTU<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getUnbalancedLoad();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFileVTU<<vel(j)<<' ';
		} else {
		    theFileVTU<<0.0<<' ';
		}
	    }
	    theFileVTU<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }

    // node mass
    if(nodedata.mass) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"NodeMass"<<quota;
	theFileVTU<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Matrix& mat = nodes[i]->getMass();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < mat.noRows()) {
		    theFileVTU<<mat(j,j)<<' ';
		} else {
		    theFileVTU<<0.0<<' ';
		}
	    }
	    theFileVTU<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }

    // node eigen vector
    for(int k=0; k<nodedata.numeigen; k++) {
	this->indent();
	theFileVTU<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFileVTU<<" Name="<<quota<<"EigenVector"<<k+1<<quota;
	theFileVTU<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFileVTU<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Matrix& eigens = nodes[i]->getEigenvectors();
	    if(k >= eigens.noCols()) {
		opserr<<"WARNING: eigenvector "<<k+1<<" is too large\n";
		return -1;
	    }
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < eigens.noRows()) {
		    theFileVTU<<eigens(j,k)<<' ';
		} else {
		    theFileVTU<<0.0<<' ';
		}
	    }
	    theFileVTU<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFileVTU<<"</DataArray>\n";
    }




    */


    // footer
    this->decrLevel();
    this->indent();
    theFileVTU<<"</Piece>\n";

    this->decrLevel();
    this->indent();
    theFileVTU<<"</UnstructuredGrid>\n";

    this->decrLevel();
    this->indent();
    theFileVTU<<"</VTKFile>\n";

    theFileVTU.close();

    return 0;
}

void
VTK_Recorder::indent() {
    for(int i=0; i<indentlevel*indentsize; i++) {
	theVTUFile<<' ';
    }
}

int
VTK_Recorder::sendSelf(int commitTag, Channel &theChannel)
{
  sendSelfCount++;

  static ID idData(2+14+1);
  int fileNameLength = 0;
  if (name != 0)
    fileNameLength = strlen(name);

  idData(0) = fileNameLength;
  idData(1) = -sendSelfCount; // -sendSelfCount indicates a process other than P0
  idData(2) = outputData.disp;
  idData(3) = outputData.disp2;
  idData(4) = outputData.disp3;
  idData(5) = outputData.vel;
  idData(6) = outputData.vel2;
  idData(7) = outputData.vel3;
  idData(8) = outputData.accel;
  idData(9) = outputData.accel2;
  idData(10) = outputData.accel3;
  idData(11) = outputData.reaction;
  idData(12) = outputData.reaction2;
  idData(13) = outputData.reaction3;
  idData(14) = outputData.mass;
  idData(15) = outputData.unbalancedLoad;

  idData(16) = precision;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "FileStream::sendSelf() - failed to send id data\n";
    return -1;
  }

  if (fileNameLength != 0) {
    Message theMessage(name, fileNameLength);
    if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
      opserr << "FileStream::sendSelf() - failed to send message\n";
      return -1;
    }
  }

  return 0;
}

int
VTK_Recorder::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(2+14+1);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "FileStream::recvSelf() - failed to recv id data\n";
    return -1;
  }

  int fileNameLength = idData(0);
  sendSelfCount = idData(1);

  outputData.disp = idData(2);
  outputData.disp2 = idData(3);
  outputData.disp3 = idData(4);
  outputData.vel = idData(5);
  outputData.vel2 = idData(6);
  outputData.vel3 = idData(7);
  outputData.accel = idData(8);
  outputData.accel2 = idData(9);
  outputData.accel3 = idData(10);
  outputData.reaction = idData(11);
  outputData.reaction2 = idData(12);
  outputData.reaction3 = idData(13);
  outputData.mass = idData(14);
  outputData.unbalancedLoad = idData(15);

  precision = idData(16);

  if (fileNameLength != 0) {
    if (name != 0)
      delete [] name;

    name = new char[fileNameLength+1];

    if (name == 0) {
      opserr << "FileStream::recvSelf() - out of memory\n";
      return -1;
    }

    Message theMessage(name, fileNameLength);
    if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
      opserr << "FileStream::recvSelf() - failed to recv message\n";
      return -1;
    }
    
    name[fileNameLength]='\0';
  }

  return 0;
}

void
VTK_Recorder::setVTKType()
{
    if (vtktypes.empty() == false) {
	return;
    }
    //    vtktypes[ELE_TAG_Subdomain] = VTK_POLY_VERTEX;
    vtktypes[ELEMENT_TAGS_WrapperElement] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ElasticBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_ModElasticBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_Beam2d] = VTK_LINE;
    vtktypes[ELE_TAG_beam2d02] = VTK_LINE;
    vtktypes[ELE_TAG_beam2d03] = VTK_LINE;
    vtktypes[ELE_TAG_beam2d04] = VTK_LINE;
    vtktypes[ELE_TAG_beam3d01] = VTK_LINE;
    vtktypes[ELE_TAG_beam3d02] = VTK_LINE;
    vtktypes[ELE_TAG_Truss] = VTK_LINE;
    vtktypes[ELE_TAG_TrussSection] = VTK_LINE;
    vtktypes[ELE_TAG_CorotTruss] = VTK_LINE;
    vtktypes[ELE_TAG_CorotTrussSection] = VTK_LINE;
    vtktypes[ELE_TAG_fElmt05] = VTK_LINE;
    vtktypes[ELE_TAG_fElmt02] = VTK_LINE;
    vtktypes[ELE_TAG_MyTruss] = VTK_LINE;
    vtktypes[ELE_TAG_ZeroLength] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthSection] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthND] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthContact2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthContact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthContactASDimplex] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthContactNTS2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthInterface2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_CoupledZeroLength] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthRocking] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_NLBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_NLBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_LargeDispBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_FourNodeQuad] = VTK_QUAD;
    vtktypes[ELE_TAG_FourNodeQuad3d] = VTK_QUAD;
    vtktypes[ELE_TAG_Tri31] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_SixNodeTri] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_BeamWithHinges2d] = VTK_LINE;
    vtktypes[ELE_TAG_BeamWithHinges3d] = VTK_LINE;
    vtktypes[ELE_TAG_EightNodeBrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_TwentyNodeBrick] = VTK_QUADRATIC_HEXAHEDRON;
    vtktypes[ELE_TAG_EightNodeBrick_u_p_U] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_TwentyNodeBrick_u_p_U] = VTK_QUADRATIC_HEXAHEDRON;
    vtktypes[ELE_TAG_FourNodeQuadUP] = VTK_QUAD;
    vtktypes[ELE_TAG_TotalLagrangianFD20NodeBrick] = VTK_QUADRATIC_HEXAHEDRON;
    vtktypes[ELE_TAG_TotalLagrangianFD8NodeBrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_EightNode_LDBrick_u_p] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_EightNode_Brick_u_p] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_TwentySevenNodeBrick] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BrickUP] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_Nine_Four_Node_QuadUP] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Twenty_Eight_Node_BrickUP] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Twenty_Node_Brick] = VTK_QUADRATIC_HEXAHEDRON;
    vtktypes[ELE_TAG_BBarFourNodeQuadUP] = VTK_QUAD;
    vtktypes[ELE_TAG_BBarBrickUP] = VTK_QUAD;
    vtktypes[ELE_TAG_PlateMITC4] = VTK_QUAD;
    vtktypes[ELE_TAG_ShellMITC4] = VTK_QUAD;
    vtktypes[ELE_TAG_ShellMITC9] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ASDShellQ4] = VTK_QUAD;
    vtktypes[ELE_TAG_ASDShellT3] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_Plate1] = VTK_QUAD;
    vtktypes[ELE_TAG_Brick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_BbarBrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_FLBrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_EnhancedQuad] = VTK_QUAD;
    vtktypes[ELE_TAG_ConstantPressureVolumeQuad] = VTK_QUAD;
    vtktypes[ELE_TAG_NineNodeMixedQuad] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_NineNodeQuad] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_EightNodeQuad] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DispBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_TimoshenkoBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_DispBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_DispBeamColumnWarping3d] = VTK_LINE;
    vtktypes[ELE_TAG_HingedBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_HingedBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_TwoPointHingedBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_TwoPointHingedBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_OnePointHingedBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_OnePointHingedBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_BeamColumnJoint2d] = VTK_QUAD;
    vtktypes[ELE_TAG_BeamColumnJoint3d] = VTK_QUAD;
    vtktypes[ELE_TAG_ForceBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_ForceBeamColumnWarping2d] = VTK_LINE;
    vtktypes[ELE_TAG_ForceBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticForceBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticForceBeamColumnWarping2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticForceBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_ForceBeamColumnCBDI2d] = VTK_LINE;
    vtktypes[ELE_TAG_ForceBeamColumnCBDI3d] = VTK_LINE;
    vtktypes[ELE_TAG_DispBeamColumn2dInt] = VTK_LINE;
    vtktypes[ELE_TAG_InternalSpring] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SimpleJoint2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Joint2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Joint3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ElastomericBearingPlasticity3d] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericBearingPlasticity2d] = VTK_LINE;
    vtktypes[ELE_TAG_TwoNodeLink] = VTK_LINE;
    vtktypes[ELE_TAG_ActuatorCorot] = VTK_LINE;
    vtktypes[ELE_TAG_Actuator] = VTK_LINE;
    vtktypes[ELE_TAG_Adapter] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ElastomericBearingBoucWen2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericBearingBoucWen3d] = VTK_LINE;
    vtktypes[ELE_TAG_FlatSliderSimple2d] = VTK_LINE;
    vtktypes[ELE_TAG_FlatSliderSimple3d] = VTK_LINE;
    vtktypes[ELE_TAG_FlatSlider2d] = VTK_LINE;
    vtktypes[ELE_TAG_FlatSlider3d] = VTK_LINE;
    vtktypes[ELE_TAG_SingleFPSimple2d] = VTK_LINE;
    vtktypes[ELE_TAG_SingleFPSimple3d] = VTK_LINE;
    vtktypes[ELE_TAG_SingleFP2d] = VTK_LINE;
    vtktypes[ELE_TAG_SingleFP3d] = VTK_LINE;
    vtktypes[ELE_TAG_DoubleFPSimple2d] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DoubleFPSimple3d] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DoubleFP2d] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DoubleFP3d] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_TripleFPSimple2d] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFPSimple3d] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFP2d] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFP3d] = VTK_LINE;
    vtktypes[ELE_TAG_MultiFP2d] = VTK_LINE;
    vtktypes[ELE_TAG_MultiFP3d] = VTK_LINE;
    vtktypes[ELE_TAG_GenericClient] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_GenericCopy] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_PY_MACRO2D] = VTK_LINE;
    vtktypes[ELE_TAG_SimpleContact2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SimpleContact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BeamContact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SurfaceLoad] = VTK_QUAD;
    vtktypes[ELE_TAG_BeamContact2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BeamEndContact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SSPquad] = VTK_QUAD;
    vtktypes[ELE_TAG_SSPquadUP] = VTK_QUAD;
    vtktypes[ELE_TAG_SSPbrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_SSPbrickUP] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_BeamContact2Dp] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BeamContact3Dp] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BeamEndContact3Dp] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Quad4FiberOverlay] = VTK_QUAD;
    vtktypes[ELE_TAG_Brick8FiberOverlay] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_QuadBeamEmbedContact] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DispBeamColumn2dThermal] = VTK_LINE;
    vtktypes[ELE_TAG_TPB1D] = VTK_LINE;
    vtktypes[ELE_TAG_TFP_Bearing] = VTK_LINE;
    vtktypes[ELE_TAG_TFP_Bearing2d] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFrictionPendulum] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFrictionPendulumX] = VTK_LINE;
    vtktypes[ELE_TAG_PFEMElement2D] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_FourNodeQuad02] = VTK_QUAD;
    vtktypes[ELE_TAG_cont2d01] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_cont2d02] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_CST] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_Truss2] = VTK_LINE;
    vtktypes[ELE_TAG_CorotTruss2] = VTK_LINE;
    vtktypes[ELE_Tag_ZeroLengthImpact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_PFEMElement3D] = VTK_TETRA;
    vtktypes[ELE_TAG_PFEMElement2DCompressible] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_PFEMElement2DBubble] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_PFEMElement2Dmini] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_ElasticTimoshenkoBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticTimoshenkoBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericBearingUFRP2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericBearingUFRP3d] = VTK_LINE;
    vtktypes[ELE_TAG_RJWatsonEQS2d] = VTK_LINE;
    vtktypes[ELE_TAG_RJWatsonEQS3d] = VTK_LINE;
    vtktypes[ELE_TAG_HDR] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericX] = VTK_LINE;
    vtktypes[ELE_TAG_LeadRubberX] = VTK_LINE;
    vtktypes[ELE_TAG_PileToe3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_N4BiaxialTruss] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ShellDKGQ] = VTK_QUAD;
    vtktypes[ELE_TAG_ShellNLDKGQ] = VTK_QUAD;
    vtktypes[ELE_TAG_MultipleShearSpring] = VTK_LINE;
    vtktypes[ELE_TAG_MultipleNormalSpring] = VTK_LINE;
    vtktypes[ELE_TAG_KikuchiBearing] = VTK_LINE;
    vtktypes[ELE_TAG_ComponentElement2d] = VTK_LINE;
    vtktypes[ELE_TAG_YamamotoBiaxialHDR] = VTK_LINE;
    vtktypes[ELE_TAG_MVLEM] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SFI_MVLEM] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_MVLEM_3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SFI_MVLEM_3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_E_SFI_MVLEM_3D] = VTK_POLY_VERTEX;
	vtktypes[ELE_TAG_E_SFI] = VTK_POLY_VERTEX;
	vtktypes[ELE_TAG_MEFI] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_PFEMElement2DFIC] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_TaylorHood2D] = VTK_QUADRATIC_TRIANGLE;
    vtktypes[ELE_TAG_PFEMElement2DQuasi] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_MINI] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_CatenaryCable] = VTK_LINE;
    vtktypes[ELE_TAG_FourNodeTetrahedron] = VTK_TETRA;
    vtktypes[ELE_TAG_PFEMElement3DBubble] = VTK_TETRA;
    vtktypes[ELE_TAG_TriSurfaceLoad] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_ShellANDeS] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_ShellDKGT] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_ShellNLDKGT] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_InertiaTruss] = VTK_LINE;
    vtktypes[ELE_TAG_ASDAbsorbingBoundary2D] = VTK_QUAD;
    vtktypes[ELE_TAG_ASDAbsorbingBoundary3D] = VTK_HEXAHEDRON;
}


int
VTK_Recorder::initialize()
{
  theNodeMapping.clear();
  theEleMapping.clear();
  theNodeTags.clear();
  theEleTags.clear();
  theEleClassTags.clear();
  theEleVtkTags.clear();
  theEleVtkOffsets.clear();


  //
  // create a list of node tags and a mapping for node tags to vtk points
  //    while at it determine max spatial dimension of mesh and max number of nodal dof
  //

  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  
  numNode = 0;
  maxNDM = 0;
  maxNDF = 0;

  while ((theNode = theNodes()) != 0) {

    int nodeTag = theNode->getTag();
    const Vector &crd=theNode->getCrds();
    if (crd.Size() > maxNDM)
      maxNDM = crd.Size();
    const Vector &disp=theNode->getTrialDisp();
    if (disp.Size() > maxNDF)
      maxNDF = disp.Size();

    theNodeMapping[nodeTag]=numNode;
    theNodeTags.push_back(nodeTag);
    numNode++;
  }

  //
  // create a list of ele tags and a mapping form element to cells
  //

  ElementIter &theElements = theDomain->getElements();
  Element *theElement;
  
  numElement = 0;
  int offset = 0;
  std::map<int,VTK_Recorder::VtkType>::iterator it;
  while ((theElement = theElements()) != 0) {
    int eleTag = theElement->getTag();
    int classTag = theElement->getClassTag();
    it = vtktypes.find(classTag);
    if (it != vtktypes.end()) {
      int vtkType = vtktypes[classTag];    
      //      if (vtkType != 0) {
      theEleMapping[eleTag]=numElement;
      theEleTags.push_back(eleTag);
      theEleClassTags.push_back(classTag);
      theEleVtkTags.push_back(vtkType);
      const ID &theNodes=theElement->getExternalNodes();
      int numNode = theNodes.Size();
      offset += numNode;
      theEleVtkOffsets.push_back(offset);
      numElement++;
    } else {
      if (classTag != ELE_TAG_Subdomain)
	opserr << "VTK_Recorder::init Element: " << eleTag << " of unknown vtk type\n";
    }
  }

  //
  // get information on nodal constraints .. may need to be done at each record step!
  //

  /*
    int nodeID = nodPtr->getTag();
    std::multimap<int,SP_Constraint*>::iterator first = allSPs.lower_bound(nodeID);
    std::multimap<int,SP_Constraint*>::iterator last = allSPs.upper_bound(nodeID);
    for(std::multimap<int,SP_Constraint*>::iterator it=first; it!=last; it++) {
    spPtr = it->second;
    const ID &id = dofPtr->getID();
    int dof = spPtr->getDOF_Number();		
    if (id(dof) == -2) {
    dofPtr->setID(spPtr->getDOF_Number(),-1);
    countDOF--;	
    } else {
    opserr << "WARNING PlainHandler::handle() - ";
    opserr << " multiple single pointconstraints at DOF " << dof;
    opserr << " for node " << spPtr->getNodeTag() << endln;
    }
    }
  */

  initDone = true;

  return 0;
}

int VTK_Recorder::flush(void) {
  if (thePVDFile.is_open() && thePVDFile.good()) {
    thePVDFile.flush();
  }
  if (theVTUFile.is_open() && theVTUFile.good()) {
    theVTUFile.flush();
  }
  return 0;
}
