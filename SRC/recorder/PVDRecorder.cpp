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

#include "PVDRecorder.h"
#include <sstream>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <Pressure_Constraint.h>
#include <Pressure_ConstraintIter.h>
#include <Matrix.h>
#include <classTags.h>
#include <NodeIter.h>

#include "PFEMElement/BackgroundDef.h"
#include "PFEMElement/Particle.h"
#include "PFEMElement/ParticleGroup.h"

std::map<int,PVDRecorder::VtkType> PVDRecorder::vtktypes;

void* OPS_PVDRecorder()
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
    PVDRecorder::NodeData nodedata;
    std::vector<PVDRecorder::EleData> eledata;
    double dT = 0.0;
    double rTolDt = 0.00001;
    while(numdata > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type, "disp") == 0) {
	    nodedata.disp = true;
	} else if(strcmp(type, "vel") == 0) {
	    nodedata.vel = true;
	} else if(strcmp(type, "accel") == 0) {
	    nodedata.accel = true;
	} else if(strcmp(type, "incrDisp") == 0) {
	    nodedata.incrdisp = true;
	} else if(strcmp(type, "reaction") == 0) {
	    nodedata.reaction = true;
	} else if(strcmp(type, "pressure") == 0) {
	    nodedata.pressure = true;
	} else if(strcmp(type, "unbalancedLoad") == 0) {
	    nodedata.unbalanced = true;
	} else if(strcmp(type, "mass") == 0) {
	    nodedata.mass = true;
	} else if(strcmp(type, "eigen") == 0) {
	    numdata = OPS_GetNumRemainingInputArgs();
	    if(numdata < 1) {
		opserr<<"WARNING: eigen needs 'numEigenvector'\n";
		return 0;
	    }
	    numdata = 1;
	    if(OPS_GetIntInput(&numdata,&nodedata.numeigen) < 0) {
		opserr << "WARNING: failed to read numeigen\n";
		return 0;
	    }

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
	    PVDRecorder::EleData edata;
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
    return new PVDRecorder(name,nodedata,eledata,indent,precision,dT, rTolDt);
}

PVDRecorder::PVDRecorder(const char *name, const NodeData& ndata,
			 const std::vector<EleData>& edata, int ind, int pre,
			 double dt, double rTolDt)
    :Recorder(RECORDER_TAGS_PVDRecorder), indentsize(ind), precision(pre),
     indentlevel(0), pathname(), basename(),
     timestep(), timeparts(), theFile(), quota('\"'), parts(),
     nodedata(ndata), eledata(edata), theDomain(0), partnum(),
     dT(dt), relDeltaTTol(rTolDt), nextTime(0.0)
{
    PVDRecorder::setVTKType();
    getfilename(name);
}

PVDRecorder::PVDRecorder()
    :Recorder(RECORDER_TAGS_PVDRecorder)
{
}


PVDRecorder::~PVDRecorder()
{
}

// PVD
// part 0 - all nodes
// part 1 - all particles
// part n - element type n
int
PVDRecorder::record(int ctag, double timestamp)
{
    if (dT == 0.0 || timestamp - nextTime >= -dT * relDeltaTTol) {
      if (dT > 0.0) {
        nextTime = timestamp + dT;
      }


      if (dT > 0) {
        nextTime = timestamp+dT;
      }

      if(precision==0)
         return 0;

      // get current time
      timestep.push_back(timestamp);

      // save vtu file
      if(vtu() < 0) return -1;

      // save pvd file
      if(pvd() < 0) return -1;
    }
    return 0;
}

int
PVDRecorder::restart()
{
    timestep.clear();
    timeparts.clear();
    return 0;
}

int
PVDRecorder::domainChanged()
{
    return 0;
}

int
PVDRecorder::setDomain(Domain& domain)
{
    theDomain = &domain;
    return 0;
}

int
PVDRecorder::pvd()
{
    // open pvd file
    theFile.close();
    std::string pvdname = pathname+basename+".pvd";

    theFile.open(pvdname.c_str(), std::ios::trunc|std::ios::out);
    if(theFile.fail()) {
	opserr<<"WARNING: Failed to open file "<<pvdname.c_str()<<"\n";
	return -1;
    }
    theFile.precision(precision);
    theFile << std::scientific;

    // header
    theFile<<"<?xml version="<<quota<<"1.0"<<quota<<"?>\n";
    theFile<<"<VTKFile type="<<quota<<"Collection"<<quota;
    theFile<<" compressor="<<quota<<"vtkZLibDataCompressor"<<quota;
    theFile<<">\n";

    // collection
    this->incrLevel();
    this->indent();
    theFile<<"<Collection>\n";

    // all data files
    this->incrLevel();
    for(int i=0; i<(int)timestep.size(); i++) {
	double t = timestep[i];
	const ID& partno = timeparts[i];
	for(int j=0; j<partno.Size(); j++) {
	    this->indent();
	    theFile<<"<DataSet timestep="<<quota<<t<<quota;
	    theFile<<" group="<<quota<<quota;
	    theFile<<" part="<<quota<<partno(j)<<quota;
	    theFile<<" file="<<quota<<basename.c_str();
	    theFile<<"/"<<basename.c_str()<<"_T"<<t<<"_P";
	    theFile<<partno(j)<<".vtu"<<quota;
	    theFile<<"/>\n";
	}
    }

    // end colloection
    this->decrLevel();
    this->indent();
    theFile<<"</Collection>\n";

    // end VTKFile
    this->decrLevel();
    this->indent();
    theFile<<"</VTKFile>\n";

    theFile.close();

    return 0;
}

int
PVDRecorder::vtu()
{
    if (theDomain == 0) {
	opserr << "WARNING: failed to get domain -- PVDRecorder::vtu\n";
	return -1;
    }
    // get node ndf
    NodeIter& theNodes = theDomain->getNodes();
    Node* theNode = 0;
    int nodendf = 0;
    while ((theNode = theNodes()) != 0) {
	if(nodendf < theNode->getNumberDOF()) {
	    nodendf = theNode->getNumberDOF();
	}
    }
    if (nodendf < 3) {
	nodendf = 3;
    } else if (nodendf > 3) {
        nodendf = 3;
    }

    // get parts
    this->getParts();

    // get background mesh
    VInt gtags;
    TaggedObjectIter& meshes = OPS_getAllMesh();
    Mesh* mesh = 0;
    while((mesh = dynamic_cast<Mesh*>(meshes())) != 0) {
	ParticleGroup* group = dynamic_cast<ParticleGroup*>(mesh);
	if (group == 0) {
	    continue;
	}
        gtags.push_back(group->getTag());
    }


    // part 0: all nodes
    ID partno(0, (int)parts.size()+(int)gtags.size()+1);
    partno[0] = 0;
    if (this->savePart0(nodendf) < 0) {
        return -1;
    }

    // particle parts
    for (int i=0; i<(int)gtags.size(); ++i) {
        partno[1+i] = 1+i;
        if (this->savePartParticle(1+i, gtags[i],nodendf) < 0) {
            return -1;
        }
    }

    // save other parts
    // int index = 1;
    for(std::map<int,ID>::iterator it=parts.begin(); it!=parts.end(); it++) {
	// int& no = partnum[it->first];
	// if (no == 0) {
	//     no = (int)partnum.size();
	// }
	int no = partno.Size();
	partno[no] = no;
	if(this->savePart(no,it->first,nodendf) < 0) return -1;
    }


    timeparts.push_back(partno);

    // clear parts
    parts.clear();

    return 0;
}

void
PVDRecorder::getParts()
{
    if (theDomain == 0) {
	opserr<<"WARNING: setDomain has not been called -- PVDRecorder\n";
	return;
    }

    ElementIter* eiter = &(theDomain->getElements());
    Element* theEle = 0;
    while((theEle = (*eiter)()) != 0) {
	int ctag = theEle->getClassTag();
	int etag = theEle->getTag();
	parts[ctag].insert(etag);
    }
}

int
PVDRecorder::savePart0(int nodendf)
{
    if (theDomain == 0) {
	opserr<<"WARNING: setDomain has not been called -- PVDRecorder\n";
	return -1;
    }

    // get time and part
    std::stringstream ss;
    ss.precision(precision);
    ss << std::scientific;
    ss << 0 << ' ' << timestep.back();
    std::string stime, spart;
    ss >> spart >> stime;

    // open file
    theFile.close();
    std::string vtuname = pathname+basename+"/"+basename+"_T"+stime+"_P"+spart+".vtu";
    theFile.open(vtuname.c_str(), std::ios::trunc|std::ios::out);
    if(theFile.fail()) {
	opserr<<"WARNING: Failed to open file "<<vtuname.c_str()<<"\n";
	return -1;
    }
    theFile.precision(precision);
    theFile << std::scientific;

    // header
    theFile<<"<?xml version="<<quota<<"1.0"<<quota<<"?>\n";
    theFile<<"<VTKFile type="<<quota<<"UnstructuredGrid"<<quota;
    theFile<<" version="<<quota<<"1.0"<<quota;
    theFile<<" byte_order="<<quota<<"LittleEndian"<<quota;
    theFile<<" compressor="<<quota<<"vtkZLibDataCompressor"<<quota;
    theFile<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<"<UnstructuredGrid>\n";

    // get pressure nodes
    ID ptags(0,theDomain->getNumPCs());
    Pressure_ConstraintIter& thePCs = theDomain->getPCs();
    Pressure_Constraint* thePC = 0;
    while ((thePC = thePCs()) != 0) {
	Node* pnode = thePC->getPressureNode();
	if (pnode != 0) {
	    ptags.insert(pnode->getTag());
	}
    }

    // get all nodes except pressure nodes
    std::vector<Node*> nodes;
    NodeIter& theNodes = theDomain->getNodes();
    Node* theNode = 0;
    while ((theNode = theNodes()) != 0) {
	int nd = theNode->getTag();
	if (ptags.getLocationOrdered(nd) < 0) {
	    nodes.push_back(theNode);
	}
    }

    // Piece
    this->incrLevel();
    this->indent();
    theFile<<"<Piece NumberOfPoints="<<quota<<(int)nodes.size()<<quota;
    theFile<<" NumberOfCells="<<quota<<1<<quota<<">\n";

    // points
    this->incrLevel();
    this->indent();
    theFile<<"<Points>\n";

    // points header
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
    theFile<<" Name="<<quota<<"Points"<<quota;
    theFile<<" NumberOfComponents="<<quota<<3<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";

    // points coordinates
    this->incrLevel();
    for(int i=0; i<(int)nodes.size(); i++) {
	const Vector& crds = nodes[i]->getCrds();
	this->indent();
	for(int j=0; j<3; j++) {
	    if(j < crds.Size()) {
		theFile<<crds(j)<<' ';
	    } else {
		theFile<<0.0<<' ';
	    }
	}
	theFile<<std::endl;
    }

    // points footer
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";
    this->decrLevel();
    this->indent();
    theFile<<"</Points>\n";

    // cells
    this->indent();
    theFile<<"<Cells>\n";

    // connectivity
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"connectivity"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    for(int i=0; i<(int)nodes.size(); i++) {
	this->indent();
	theFile<<i<<std::endl;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // offsets
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"offsets"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<(int)nodes.size()<<std::endl;
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // types
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"types"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<VTK_POLY_VERTEX<<std::endl;
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // cells footer
    this->decrLevel();
    this->indent();
    theFile<<"</Cells>\n";

    // point data
    this->indent();
    theFile<<"<PointData>\n";

    // node tags
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"NodeTag"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    for(int i=0; i<(int)nodes.size(); i++) {
	this->indent();
	theFile<<nodes[i]->getTag()<<std::endl;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // node velocity
    if(nodedata.vel) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Velocity"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getTrialVel();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node displacement
    if(nodedata.disp) {
	// all displacement
    // this->indent();
	// theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	// theFile<<" Name="<<quota<<"AllDisplacement"<<quota;
	// theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	// theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	// this->incrLevel();
	// for(int i=0; i<(int)nodes.size(); i++) {
	//     const Vector& vel = nodes[i]->getTrialDisp();
	//     this->indent();
	//     for(int j=0; j<nodendf; j++) {
	// 	if(j < vel.Size()) {
	// 	    theFile<<vel(j)<<' ';
	// 	} else {
	// 	    theFile<<0.0<<' ';
	// 	}
	//     }
	//     theFile<<std::endl;
	// }
	// this->decrLevel();
	// this->indent();
	// theFile<<"</DataArray>\n";

    // displacement
    this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Displacement"<<quota;
	theFile<<" NumberOfComponents="<<quota<<3<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getTrialDisp();
	    this->indent();
	    for(int j=0; j<3; j++) {
		if(j < vel.Size() && j < nodes[i]->getCrds().Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node incr displacement
    if(nodedata.incrdisp) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"IncrDisplacement"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getIncrDisp();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node acceleration
    if(nodedata.accel) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Acceleration"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getTrialAccel();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node pressure
    if(nodedata.pressure) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Pressure"<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    double pressure = 0.0;
	    Pressure_Constraint* thePC = theDomain->getPressure_Constraint(nodes[i]->getTag());
	    if(thePC != 0) {
		pressure = thePC->getPressure();
	    }
	    this->indent();
	    theFile<<pressure<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node reaction
    if(nodedata.reaction) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Reaction"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getReaction();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node unbalanced load
    if(nodedata.unbalanced) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"UnbalancedLoad"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Vector& vel = nodes[i]->getUnbalancedLoad();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node mass
    if(nodedata.mass) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"NodeMass"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)nodes.size(); i++) {
	    const Matrix& mat = nodes[i]->getMass();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < mat.noRows()) {
		    theFile<<mat(j,j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node eigen vector
    for(int k=0; k<nodedata.numeigen; k++) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"EigenVector"<<k+1<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
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
		    theFile<<eigens(j,k)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // point data footer
    this->decrLevel();
    this->indent();
    theFile<<"</PointData>\n";

    // cell data
    this->indent();
    theFile<<"<CellData>\n";

    // element tags
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"ElementTag"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<0<<std::endl;
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // cell data footer
    this->decrLevel();
    this->indent();
    theFile<<"</CellData>\n";

    // footer
    this->decrLevel();
    this->indent();
    theFile<<"</Piece>\n";

    this->decrLevel();
    this->indent();
    theFile<<"</UnstructuredGrid>\n";

    this->decrLevel();
    this->indent();
    theFile<<"</VTKFile>\n";

    theFile.close();

    return 0;
}

int
PVDRecorder::savePartParticle(int pno, int bgtag, int nodendf)
{
    if (theDomain == 0) {
	opserr<<"WARNING: setDomain has not been called -- PVDRecorder\n";
	return -1;
    }

    // get time and part
    std::stringstream ss;
    ss.precision(precision);
    ss << std::scientific;
    ss << pno << ' ' << timestep.back();
    std::string stime, spart;
    ss >> spart >> stime;

    // open file
    theFile.close();
    std::string vtuname = pathname+basename+"/"+basename+"_T"+stime+"_P"+spart+".vtu";
    theFile.open(vtuname.c_str(), std::ios::trunc|std::ios::out);
    if(theFile.fail()) {
	opserr<<"WARNING: Failed to open file "<<vtuname.c_str()<<"\n";
	return -1;
    }
    theFile.precision(precision);
    theFile << std::scientific;

    // header
    theFile<<"<?xml version="<<quota<<"1.0"<<quota<<"?>\n";
    theFile<<"<VTKFile type="<<quota<<"UnstructuredGrid"<<quota;
    theFile<<" version="<<quota<<"1.0"<<quota;
    theFile<<" byte_order="<<quota<<"LittleEndian"<<quota;
    theFile<<" compressor="<<quota<<"vtkZLibDataCompressor"<<quota;
    theFile<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<"<UnstructuredGrid>\n";

    // get particles in group
    VParticle particles;
    ParticleGroup* group = dynamic_cast<ParticleGroup*>(OPS_getMesh(bgtag));
    if (group == 0) {
        opserr << "WARNING: particle group "<<bgtag<<"doesn't exist\n";
        return -1;
    }
    for(int j=0; j<group->numParticles(); j++) {
	Particle* p = group->getParticle(j);
	if(p == 0) continue;
	particles.push_back(p);
    }

    // Piece
    this->incrLevel();
    this->indent();
    theFile<<"<Piece NumberOfPoints="<<quota<<(int)particles.size()<<quota;
    theFile<<" NumberOfCells="<<quota<<1<<quota<<">\n";

    // points
    this->incrLevel();
    this->indent();
    theFile<<"<Points>\n";

    // points header
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
    theFile<<" Name="<<quota<<"Points"<<quota;
    theFile<<" NumberOfComponents="<<quota<<3<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";

    // points coordinates
    this->incrLevel();
    for(int i=0; i<(int)particles.size(); i++) {
	const VDouble& crds = particles[i]->getCrds();
	this->indent();
	for(int j=0; j<3; j++) {
	    if(j < (int)crds.size()) {
		theFile<<crds[j]<<' ';
	    } else {
		theFile<<0.0<<' ';
	    }
	}
	theFile<<std::endl;
    }

    // points footer
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";
    this->decrLevel();
    this->indent();
    theFile<<"</Points>\n";

    // cells
    this->indent();
    theFile<<"<Cells>\n";

    // connectivity
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"connectivity"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    for(int i=0; i<(int)particles.size(); i++) {
	this->indent();
	theFile<<i<<std::endl;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // offsets
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"offsets"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<(int)particles.size()<<std::endl;
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // types
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"types"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<VTK_POLY_VERTEX<<std::endl;
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // cells footer
    this->decrLevel();
    this->indent();
    theFile<<"</Cells>\n";

    // point data
    this->indent();
    theFile<<"<PointData>\n";

    // node tags
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"NodeTag"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    for(int i=0; i<(int)particles.size(); i++) {
	this->indent();
	theFile<<particles[i]->getTag()<<std::endl;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // node velocity
    if(nodedata.vel) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Velocity"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    const VDouble& vel = particles[i]->getVel();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < (int)vel.size()) {
		    theFile<<vel[j]<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node displacement
    if(nodedata.disp) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Displacement"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		theFile<<0.0<<' ';
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node incr displacement
    if(nodedata.incrdisp) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"IncrDisplacement"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		theFile<<0.0<<' ';
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node acceleration
    if(nodedata.accel) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Acceleration"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		theFile<<0.0<<' ';
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node pressure
    if(nodedata.pressure) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Pressure"<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    double pressure = particles[i]->getPressure();
	    this->indent();
	    theFile<<pressure<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node reaction
    if(nodedata.reaction) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Reaction"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		theFile<<0.0<<' ';
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node unbalanced load
    if(nodedata.unbalanced) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"UnbalancedLoad"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		theFile<<0.0<<' ';
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node mass
    if(nodedata.mass) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"NodeMass"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		theFile<<0.0<<' ';
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node eigen vector
    for(int k=0; k<nodedata.numeigen; k++) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"EigenVector"<<k+1<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<(int)particles.size(); i++) {
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		theFile<<0.0<<' ';
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // point data footer
    this->decrLevel();
    this->indent();
    theFile<<"</PointData>\n";

    // cell data
    this->indent();
    theFile<<"<CellData>\n";

    // element tags
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"ElementTag"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<0<<std::endl;
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // cell data footer
    this->decrLevel();
    this->indent();
    theFile<<"</CellData>\n";

    // footer
    this->decrLevel();
    this->indent();
    theFile<<"</Piece>\n";

    this->decrLevel();
    this->indent();
    theFile<<"</UnstructuredGrid>\n";

    this->decrLevel();
    this->indent();
    theFile<<"</VTKFile>\n";

    theFile.close();

    return 0;
}

int
PVDRecorder::savePart(int partno, int ctag, int nodendf)
{
    if (theDomain == 0) {
	opserr<<"WARNING: setDomain has not been called -- PVDRecorder\n";
	return -1;
    }

    // get time and part
    std::stringstream ss;
    ss.precision(precision);
    ss << std::scientific;
    ss << partno << ' ' << timestep.back();
    std::string stime, spart;
    ss >> spart >> stime;

    // open file
    theFile.close();
    std::string vtuname = pathname+basename+"/"+basename+"_T"+stime+"_P"+spart+".vtu";
    theFile.open(vtuname.c_str(), std::ios::trunc|std::ios::out);
    if(theFile.fail()) {
	opserr<<"WARNING: Failed to open file "<<vtuname.c_str()<<"\n";
	return -1;
    }
    theFile.precision(precision);
    theFile << std::scientific;

    // header
    theFile<<"<?xml version="<<quota<<"1.0"<<quota<<"?>\n";
    theFile<<"<VTKFile type="<<quota<<"UnstructuredGrid"<<quota;
    theFile<<" version="<<quota<<"1.0"<<quota;
    theFile<<" byte_order="<<quota<<"LittleEndian"<<quota;
    theFile<<" compressor="<<quota<<"vtkZLibDataCompressor"<<quota;
    theFile<<">\n";
    this->incrLevel();
    this->indent();
    theFile<<"<UnstructuredGrid>\n";

    // get nodes
    const ID& eletags = parts[ctag];
    ID ndtags(0,eletags.Size()*3);
    std::vector<Element*> eles(eletags.Size());
    int numelenodes = 0;
    int increlenodes = 1;
    for(int i=0; i<eletags.Size(); i++) {
	eles[i] = theDomain->getElement(eletags(i));
	if (eles[i] == 0) {
	    opserr<<"WARNING: element "<<eletags(i)<<" is not defined--pvdRecorder\n";
	    return -1;
	}
	const ID& elenodes = eles[i]->getExternalNodes();
	if(numelenodes == 0) {
	    numelenodes = elenodes.Size();
	    if(ctag==ELE_TAG_PFEMElement2D||
	       ctag==ELE_TAG_PFEMElement2DCompressible||
	       ctag==ELE_TAG_PFEMElement2DBubble||
	       ctag==ELE_TAG_PFEMElement2Dmini ||
	       ctag==ELE_TAG_MINI ||
	       ctag==ELE_TAG_PFEMElement2DQuasi) {
		numelenodes = 3;
		increlenodes = 2;
	    } else if (ctag==ELE_TAG_TaylorHood2D) {
		numelenodes = 6;
		increlenodes = 1;
	    } else if (ctag==ELE_TAG_PFEMElement3DBubble) {
		numelenodes = 4;
		increlenodes = 2;
	    }
	}
	for(int j=0; j<numelenodes; j++) {
	    ndtags.insert(elenodes(j*increlenodes));
	}
    }

    // Piece
    this->incrLevel();
    this->indent();
    theFile<<"<Piece NumberOfPoints="<<quota<<ndtags.Size()<<quota;
    theFile<<" NumberOfCells="<<quota<<eletags.Size()<<quota<<">\n";

    // points
    this->incrLevel();
    this->indent();
    theFile<<"<Points>\n";

    // points header
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
    theFile<<" Name="<<quota<<"Points"<<quota;
    theFile<<" NumberOfComponents="<<quota<<3<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";

    // points coordinates
    this->incrLevel();
    std::vector<Node*> nodes(ndtags.Size());
    for(int i=0; i<ndtags.Size(); i++) {
	nodes[i] = theDomain->getNode(ndtags(i));
	if(nodes[i] == 0) {
	    opserr<<"WARNING: Node "<<ndtags(i)<<" is not defined -- pvdRecorder\n";
	    return -1;
	}
	const Vector& crds = nodes[i]->getCrds();
	this->indent();
	for(int j=0; j<3; j++) {
	    if(j < crds.Size()) {
		theFile<<crds(j)<<' ';
	    } else {
		theFile<<0.0<<' ';
	    }
	}
	theFile<<std::endl;
    }

    // points footer
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";
    this->decrLevel();
    this->indent();
    theFile<<"</Points>\n";

    // cells
    this->indent();
    theFile<<"<Cells>\n";

    // connectivity
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"connectivity"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    for(int i=0; i<eletags.Size(); i++) {
	const ID& elenodes = eles[i]->getExternalNodes();
	this->indent();
	if (ctag==ELE_TAG_TaylorHood2D) {

	    // for 2nd order element, the order of mid nodes
	    // is different to VTK
	    int vtkOrder[] = {0,1,2,5,3,4};
	    for(int j=0; j<numelenodes; j++) {
		theFile<<ndtags.getLocationOrdered(elenodes(vtkOrder[j]*increlenodes))<<' ';
	    }

	} else {

	    for(int j=0; j<numelenodes; j++) {
		theFile<<ndtags.getLocationOrdered(elenodes(j*increlenodes))<<' ';
	    }
	}
	theFile<<std::endl;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // offsets
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"offsets"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    int offset = numelenodes;
    for(int i=0; i<eletags.Size(); i++) {
	this->indent();
	theFile<<offset<<std::endl;
	offset += numelenodes;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // types
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"types"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    int type = vtktypes[ctag];
    if (type == 0) {
	opserr<<"WARNING: the element type cannot be assigned a VTK type\n";
	return -1;
    }
    for(int i=0; i<eletags.Size(); i++) {
	this->indent();
	theFile<<type<<std::endl;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // cells footer
    this->decrLevel();
    this->indent();
    theFile<<"</Cells>\n";

    // point data
    this->indent();
    theFile<<"<PointData>\n";

    // node tags
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"NodeTag"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    for(int i=0; i<ndtags.Size(); i++) {
	this->indent();
	theFile<<ndtags(i)<<std::endl;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // node velocity
    if(nodedata.vel) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Velocity"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    const Vector& vel = nodes[i]->getTrialVel();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node displacement
    if(nodedata.disp) {
	// all displacement
    // this->indent();
	// theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	// theFile<<" Name="<<quota<<"AllDisplacement"<<quota;
	// theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	// theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	// this->incrLevel();
	// for(int i=0; i<ndtags.Size(); i++) {
	//     const Vector& vel = nodes[i]->getTrialDisp();
	//     this->indent();
	//     for(int j=0; j<nodendf; j++) {
	// 	if(j < vel.Size()) {
	// 	    theFile<<vel(j)<<' ';
	// 	} else {
	// 	    theFile<<0.0<<' ';
	// 	}
	//     }
	//     theFile<<std::endl;
	// }
	// this->decrLevel();
	// this->indent();
	// theFile<<"</DataArray>\n";

    // displacement
    this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Displacement"<<quota;
	theFile<<" NumberOfComponents="<<quota<<3<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    const Vector& vel = nodes[i]->getTrialDisp();
	    this->indent();
	    for(int j=0; j<3; j++) {
		if(j < vel.Size() && j < nodes[i]->getCrds().Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node incr displacement
    if(nodedata.incrdisp) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"IncrDisplacement"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    const Vector& vel = nodes[i]->getIncrDisp();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node acceleration
    if(nodedata.accel) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Acceleration"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    const Vector& vel = nodes[i]->getTrialAccel();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node pressure
    if(nodedata.pressure) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Pressure"<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    double pressure = 0.0;
	    Pressure_Constraint* thePC = theDomain->getPressure_Constraint(ndtags(i));
	    if(thePC != 0) {
		pressure = thePC->getPressure();
	    }
	    this->indent();
	    theFile<<pressure<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node reaction
    if(nodedata.reaction) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"Reaction"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    const Vector& vel = nodes[i]->getReaction();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node unbalanced load
    if(nodedata.unbalanced) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"UnbalancedLoad"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    const Vector& vel = nodes[i]->getUnbalancedLoad();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < vel.Size()) {
		    theFile<<vel(j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node mass
    if(nodedata.mass) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"NodeMass"<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    const Matrix& mat = nodes[i]->getMass();
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < mat.noRows()) {
		    theFile<<mat(j,j)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // node eigen vector
    for(int k=0; k<nodedata.numeigen; k++) {
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<"EigenVector"<<k+1<<quota;
	theFile<<" NumberOfComponents="<<quota<<nodendf<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int i=0; i<ndtags.Size(); i++) {
	    const Matrix& eigens = nodes[i]->getEigenvectors();
	    if(k >= eigens.noCols()) {
		opserr<<"WARNING: eigenvector "<<k+1<<" is too large\n";
		return -1;
	    }
	    this->indent();
	    for(int j=0; j<nodendf; j++) {
		if(j < eigens.noRows()) {
		    theFile<<eigens(j,k)<<' ';
		} else {
		    theFile<<0.0<<' ';
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // point data footer
    this->decrLevel();
    this->indent();
    theFile<<"</PointData>\n";

    // cell data
    this->indent();
    theFile<<"<CellData>\n";

    // element tags
    this->incrLevel();
    this->indent();
    theFile<<"<DataArray type="<<quota<<"Int64"<<quota;
    theFile<<" Name="<<quota<<"ElementTag"<<quota;
    theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
    this->incrLevel();
    for(int i=0; i<eletags.Size(); i++) {
	this->indent();
	theFile<<eletags(i)<<std::endl;
    }
    this->decrLevel();
    this->indent();
    theFile<<"</DataArray>\n";

    // element response
    for(int i=0; i<(int)eledata.size(); i++) {

	if(eletags.Size() == 0) break;

	// check data
	int argc = (int)eledata[i].size();
	if(argc == 0) continue;
	std::vector<const char*> argv(argc);
	for(int j=0; j<argc; j++) {
	    argv[j] = eledata[i][j].c_str();
	}
	const Vector* data =theDomain->getElementResponse(eletags(0),&(argv[0]),argc);
	if(data==0) continue;
	int eressize = data->Size();
	if(eressize == 0) continue;

	// save data
	this->indent();
	theFile<<"<DataArray type="<<quota<<"Float64"<<quota;
	theFile<<" Name="<<quota<<eles[0]->getClassType();
	for(int j=0; j<argc; j++) {
	    theFile<<argv[j];
	}
	theFile<<quota;
	theFile<<" NumberOfComponents="<<quota<<eressize<<quota;
	theFile<<" format="<<quota<<"ascii"<<quota<<">\n";
	this->incrLevel();
	for(int j=0; j<eletags.Size(); j++) {
	    data=theDomain->getElementResponse(eletags(j),&(argv[0]),argc);
	    if(data==0) {
		opserr<<"WARNING: can't get response for element "<<eletags(j)<<"\n";
		return -1;
	    }
	    this->indent();
	    for(int k=0; k<eressize; k++) {
		if (k>=data->Size()) {
		    theFile<<0.0<<" ";
		} else {
		    theFile<<(*data)(k)<<" ";
		}
	    }
	    theFile<<std::endl;
	}
	this->decrLevel();
	this->indent();
	theFile<<"</DataArray>\n";
    }

    // cell data footer
    this->decrLevel();
    this->indent();
    theFile<<"</CellData>\n";

    // footer
    this->decrLevel();
    this->indent();
    theFile<<"</Piece>\n";

    this->decrLevel();
    this->indent();
    theFile<<"</UnstructuredGrid>\n";

    this->decrLevel();
    this->indent();
    theFile<<"</VTKFile>\n";

    theFile.close();

    return 0;
}

void
PVDRecorder::indent() {
    for(int i=0; i<indentlevel*indentsize; i++) {
	theFile<<' ';
    }
}

int
PVDRecorder::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
PVDRecorder::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
PVDRecorder::setVTKType()
{
    if (vtktypes.empty() == false) {
	return;
    }
    vtktypes[ELE_TAG_Subdomain] = VTK_POLY_VERTEX;
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
    vtktypes[ELE_TAG_PFEMContact2D] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_PFEMContact3D] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_InertiaTruss] = VTK_LINE;
    vtktypes[ELE_TAG_ASDAbsorbingBoundary2D] = VTK_QUAD;
    vtktypes[ELE_TAG_ASDAbsorbingBoundary3D] = VTK_HEXAHEDRON;
}

void
PVDRecorder::getfilename(const char* name)
{
    // use string
    std::string fname(name);
    
    // no slash at all
    std::size_t found = fname.find_last_of("/\\");
    if (found == std::string::npos) {
	pathname = "./";
	basename = fname;
	return;
    }

    // remove trailing slash
    if (found == fname.length()-1) {
	fname = fname.substr(0,fname.length()-1);
	found = fname.find_last_of("/\\");
    }
    
    // only trailing slash
    if(found == std::string::npos) {
	pathname = "./";
	basename = fname;
	return;
    }

    // more slash
    pathname = fname.substr(0,found+1);
    basename = fname.substr(found+1);
}

int PVDRecorder::flush(void) {
  if (theFile.is_open() && theFile.good()) {
    theFile.flush();
  }
  return 0;
}
