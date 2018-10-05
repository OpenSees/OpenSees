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

// $Revision$
// $Date$
// $URL$

// Written: Minjie Zhu
//
// Description: This class defines the BackgroundMesh
//

#include "BackgroundMesh.h"
#include "ParticleGroup.h"
#include "PFEMElement2DBubble.h"
#include "PFEMElement3DBubble.h"
#include "PFEMElement2Dmini.h"
#include "PFEMElement2DCompressible.h"
#include <elementAPI.h>
#include <cmath>
#ifdef _LINUX
#include <Timer.h>
#endif
#include <omp.h>
#include <iostream>
#include <Node.h>
#include <NodeIter.h>
#include <Domain.h>
#include <algorithm>
#include <Pressure_Constraint.h>
#include <TriangleMeshGenerator.h>
#include <TetMeshGenerator.h>
#include <time.h>
#include <fstream>
#include <SP_Constraint.h>
#include <Matrix.h>

int BackgroundMesh::FLUID = 1;
int BackgroundMesh::STRUCTURE = 2;
int BackgroundMesh::FIXED = 3;
int BackgroundMesh::EMPTY = 4;
int BackgroundMesh::FSI = 5;

static BackgroundMesh bgmesh;

BackgroundMesh& OPS_getBgMesh()
{
    return bgmesh;
}

// OPS_BgMesh
int OPS_BgMesh()
{
    int ndm = OPS_GetNDM();

    // check input
    if(OPS_GetNumRemainingInputArgs() < 2*ndm+1) {
	std::cerr<<"WARNING: basicsize? lower? upper? <-tol tol? -meshtol tol? -wave wavefilename? numl? locs? -freesurface -numsub numsub? -fixed numnodes? fixNodes? -structure numnodes? structuralNodes?>\n";
	return -1;
    }

    // get basicsize
    double size;
    if (OPS_GetNumRemainingInputArgs() < 1) {
	std::cerr << "WARNING: need basic size\n";
	return -1;
    }
    int num = 1;
    if (OPS_GetDoubleInput(&num, &size) < 0) {
	std::cerr << "WARNING: failed to get basic size\n";
	return -1;
    }
    if (size <= 0) {
	std::cerr << "WARNING: basic size <= 0\n";
	return -1;
    }
    bgmesh.setBasicSize(size);

    // get range
    VDouble lower(ndm), upper(ndm);
    double* ptr = &lower[0];
    if (OPS_GetDoubleInput(&ndm, ptr) < 0) {
	std::cerr << "WARNING: failed to get min\n";
	return -1;
    }
    ptr = &upper[0];
    if (OPS_GetDoubleInput(&ndm, ptr) < 0) {
	std::cerr << "WARNING: failed to get max\n";
	return -1;
    }
    bgmesh.setRange(lower,upper);

    // get tolerance
    while (OPS_GetNumRemainingInputArgs()>0) {

	const char* opt = OPS_GetString();

	if (strcmp(opt, "-tol") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING: need tol\n";
		return -1;
	    }
	    double tol;
	    if (OPS_GetDoubleInput(&num, &tol) < 0) {
		std::cerr << "WARNING: failed to read tolerance\n";
		return -1;
	    }
	    bgmesh.setTol(tol);
	} else if (strcmp(opt, "-meshtol") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING: need tol\n";
		return -1;
	    }
	    double tol;
	    if (OPS_GetDoubleInput(&num, &tol) < 0) {
		std::cerr << "WARNING: failed to read mesh tolerance\n";
		return -1;
	    }
	    bgmesh.setMeshTol(tol);

	} else if (strcmp(opt, "-wave") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 2) {
		opserr << "WARNING: need wavefilename and numl\n";
		return -1;
	    }

	    const char* wavefilename = OPS_GetString();
	    if(bgmesh.setFile(wavefilename) < 0) {
		return -1;
	    }

	    int numl;
	    num = 1;
	    if (OPS_GetIntInput(&num, &numl)) {
		opserr << "WARNING: failed to read numl\n";
		return -1;
	    }

	    num = numl*ndm;
	    if (OPS_GetNumRemainingInputArgs() < num) {
		opserr << "WARNING: insufficient number of locations\n";
		return -1;
	    }

	    VDouble locs(num);
	    if (OPS_GetDoubleInput(&num, &locs[0]) < 0) {
		std::cerr << "WARNING: failed to read wave recording locations\n";
		return -1;
	    }
	    bgmesh.setLocs(locs);
	} else if (strcmp(opt, "-freesurface") == 0) {
	    bgmesh.setFreeSurf();
	} else if (strcmp(opt, "-numsub") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING: need numsub\n";
		return -1;
	    }
	    int numsub;
	    num = 1;
	    if (OPS_GetIntInput(&num, &numsub)) {
		opserr << "WARNING: failed to read numsub\n";
		return -1;
	    }
	    bgmesh.setNumSub(numsub);
	} else if (strcmp(opt, "-structure") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING: need numnodes\n";
		return -1;
	    }
	    int numnodes;
	    num = 1;
	    if (OPS_GetIntInput(&num, &numnodes)) {
		opserr << "WARNING: failed to read numnodes\n";
		return -1;
	    }
	    if (OPS_GetNumRemainingInputArgs() < numnodes) {
		opserr << "WARNING: insufficient number of structural nodes\n";
		return -1;
	    }
	    VInt snodes(numnodes);
	    if (OPS_GetIntInput(&numnodes, &snodes[0])) {
		opserr << "WARNING: failed to read structural nodes\n";
		return -1;
	    }
	    bgmesh.addStructuralNodes(snodes);
	} else if (strcmp(opt, "-fixed") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING: need numnodes\n";
		return -1;
	    }
	    int numnodes;
	    num = 1;
	    if (OPS_GetIntInput(&num, &numnodes)) {
		opserr << "WARNING: failed to read numnodes\n";
		return -1;
	    }
	    if (OPS_GetNumRemainingInputArgs() < numnodes) {
		opserr << "WARNING: insufficient number of fixed nodes\n";
		return -1;
	    }
	    VInt fixednodes(numnodes);
	    if (OPS_GetIntInput(&numnodes, &fixednodes[0])) {
		opserr << "WARNING: failed to read fixed nodes\n";
		return -1;
	    }
	    bgmesh.addFixed(fixednodes);
	}
    }

    // turn off disp in PFEM elements
    PFEMElement2DBubble::dispon = false;
    PFEMElement3DBubble::dispon = false;
    PFEMElement2DCompressible::dispon = false;
    PFEMElement2Dmini::dispon = false;

    // bg mesh
    if (bgmesh.remesh(true) < 0) {
	opserr << "WARNING: failed to create background mesh\n";
	return -1;	
    }

    return 0;
}

BackgroundMesh::BackgroundMesh()
    :lower(), upper(), bcells(),
     tol(1e-10), meshtol(0.1), bsize(-1.0),
     numave(2), numsub(4), recorders(),locs(),
     currentTime(0.0), theFile(),
     freesurf(false), atomsP(0.0),
     fixedNodes(), structuralNodes()
{
}

BackgroundMesh::~BackgroundMesh()
{
    for (unsigned int i=0; i<recorders.size(); ++i) {
	if (recorders[i] != 0) {
	    delete recorders[i];
	}
    }
    recorders.clear();
}

void
BackgroundMesh::addRecorder(Recorder* recorder)
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;

    recorder->setDomain(*domain);
    recorders.push_back(recorder);
}

void
BackgroundMesh::setRange(const VDouble& l, const VDouble& u)
{
    nearIndex(l,lower);
    nearIndex(u,upper);
}

int
BackgroundMesh::setFile(const char* name)
{
    theFile.open(name, std::ios::trunc|std::ios::out);
    if(theFile.fail()) {
	opserr<<"WARNING: Failed to open file "<<name<<"\n";
	return -1;
    }

    int precision = 12;
    theFile.precision(precision);
    theFile << std::scientific;

    return 0;
}

void
BackgroundMesh::addStructuralNodes(VInt& snodes)
{
    for (unsigned int i=0; i<snodes.size(); ++i) {
	structuralNodes.insert(snodes[i]);
    }
}

void
BackgroundMesh::getIndex(const VDouble& crds, double incr, VInt& index) const
{
    index.resize(crds.size());
    for (unsigned int i=0; i<crds.size(); ++i) {
	double crd = crds[i]/bsize + incr;
	index[i] = (int)floor(crd);
    }
}

void
BackgroundMesh::lowerIndex(const VDouble& crds, VInt& index) const
{
    getIndex(crds,0.0,index);
}

void
BackgroundMesh::upperIndex(const VDouble& crds, VInt& index) const
{
    getIndex(crds,1.0,index);
}

void
BackgroundMesh::nearIndex(const VDouble& crds, VInt& index) const
{
    getIndex(crds,0.5,index);
}

void
BackgroundMesh::getCrds(const VInt& index, VDouble& crds) const
{
    crds.resize(index.size(), 0.0);
    for (unsigned int i=0; i<crds.size(); ++i) {
	crds[i] = index[i]*bsize;
    }
}

// get corners from index
// count num+1 in each direction
void
BackgroundMesh::getCorners(const VInt& index, int num, VVInt& indices) const
{
    int ndm = OPS_GetNDM();
    int counter = 0;

    if (ndm == 2) {
	indices.resize((num+1)*(num+1));
	for (int j=index[1]; j<=index[1]+num; j+=1) {
	    for (int i=index[0]; i<=index[0]+num; i+=1) {
		indices[counter].resize(ndm);
		indices[counter][0] = i;
		indices[counter][1] = j;
		++counter;
	    }
	}
    } else if (ndm == 3) {
	indices.resize((num+1)*(num+1)*(num+1));
	for (int k=index[2]; k<=index[2]+num; k+=1) {
	    for (int j=index[1]; j<=index[1]+num; j+=1) {
		for (int i=index[0]; i<=index[0]+num; i+=1) {
		    indices[counter].resize(ndm);
		    indices[counter][0] = i;
		    indices[counter][1] = j;
		    indices[counter][2] = k;
		    ++counter;
		}
	    }
	}
    }
    
}

// gather particles from minind to maxind (not included)
// if checkfsi = true, skip fluid cells
void
BackgroundMesh::gatherParticles(const VInt& minind, const VInt& maxind,
				VParticle& pts, bool checkfsi)
{
    int ndm = OPS_GetNDM();
    pts.clear();
    VInt index(ndm);
    if (ndm == 2) {
	for (int i=minind[0]; i<maxind[0]; ++i) {
	    index[0] = i;
	    for (int j=minind[1]; j<maxind[1]; ++j) {
		index[1] = j;
		std::map<VInt,BCell>::iterator it = bcells.find(index);
		if (it != bcells.end()) {
		    BCell& cell = it->second;
		    if (checkfsi && cell.type==FLUID) continue;
		    pts.insert(pts.end(), cell.pts.begin(), cell.pts.end());
		}
	    }
	}
    } else if (ndm == 3) {
	for (int i=minind[0]; i<maxind[0]; ++i) {
	    index[0] = i;
	    for (int j=minind[1]; j<maxind[1]; ++j) {
		index[1] = j;
		for (int k=minind[2]; k<maxind[2]; ++k) {
		    index[2] = k;
		    std::map<VInt,BCell>::iterator it = bcells.find(index);
		    if (it != bcells.end()) {
			BCell& cell = it->second;
			if (checkfsi && cell.type==FLUID) continue;
			pts.insert(pts.end(), cell.pts.begin(), cell.pts.end());
		    }
		}
	    }
	}
    }
}

double
BackgroundMesh::QuinticKernel(double q, double h, int ndm)
{
    static double pi = 3.141592653589793;
    if (q<0 || q>2) return 0.0;
    double aD = 0.0;
    if (ndm == 2) {
	aD = 7.0/(4*pi*h*h);
    } else if (ndm == 3) {
	aD = 7.0/(8*pi*h*h*h);
    }
    double a = 1.0-q/2.0;

    return aD*a*a*a*a*(2*q+1);
}

int
BackgroundMesh::preNForTri(double x1, double y1, double x2, double y2,
			   double x3, double y3, VDouble& coeff)
{
    coeff.resize(9,0.0);

    coeff[0] = x2*y3-x3*y2;
    coeff[1] = x3*y1-x1*y3;
    coeff[2] = x1*y2-x2*y1;

    coeff[3] = y2-y3;
    coeff[4] = y3-y1;
    coeff[5] = y1-y2;

    coeff[6] = x3-x2;
    coeff[7] = x1-x3;
    coeff[8] = x2-x1;

    double A = coeff[0]+coeff[1]+coeff[2];

    if (A<0 || fabs(A) < 1e-15) {
	opserr << "A <= 0\n";
	return -1;
    }

    for (unsigned int i=0; i<coeff.size(); ++i) {
	coeff[i] /= A;
    }

    return 0;
}

int
BackgroundMesh::preNForTet(const VDouble& crds1, const VDouble& crds2,
			   const VDouble& crds3, const VDouble& crds4,
			   VVDouble& coeff)
{
    int ndm = OPS_GetNDM();
    if (ndm != 3) {
	return 0;
    }
    if ((int)crds1.size() < ndm) {
	return 0;
    }
    if ((int)crds2.size() < ndm) {
	return 0;
    }
    if ((int)crds3.size() < ndm) {
	return 0;
    }
    if ((int)crds4.size() < ndm) {
	return 0;
    }
    Matrix Jmat(4,4), Jfact(4,4);
    Jmat(0,0) = 1.0; Jmat(0,1) = 1.0; Jmat(0,2) = 1.0; Jmat(0,3) = 1.0;
    for (int j=0; j<ndm; ++j) {
	Jmat(j+1,0) = crds1[j];
	Jmat(j+1,1) = crds2[j];
	Jmat(j+1,2) = crds3[j];
	Jmat(j+1,3) = crds4[j];
    }

    PFEMElement3DBubble::cofactor(Jmat,Jfact);

    coeff.resize(Jfact.noRows());
    double vol = 0.0;
    for (int i=0; i<Jfact.noRows(); ++i) {
	coeff[i].resize(Jfact.noCols());
	for (int j=0; j<Jfact.noCols(); ++j) {
	    coeff[i][j] = Jfact(i,j);
	}
	vol += coeff[i][0];
    }

    if (vol < 0 || fabs(vol) < 1e-15) {
	opserr<<"vol "<<vol<<" <= 0\n";
	return -1;
    }

    for (int i=0; i<Jfact.noRows(); ++i) {
	for (int j=0; j<Jfact.noCols(); ++j) {
	    coeff[i][j] /= vol;
	}
    }

    return 0;
}

void
BackgroundMesh::getNForTri(const VDouble& coeff, double x, double y, VDouble& N)
{
    N.resize(3,0.0);

    for (unsigned int i=0; i<N.size(); ++i) {
	double val = coeff[i]+coeff[i+3]*x+coeff[i+6]*y;
	if (fabs(val) < tol) {
	    // make sure it's in the tri
	    N[i] = tol;
	}
	N[i] = val;
    }
}

void
BackgroundMesh::getNForTet(const VVDouble& coeff, const VDouble& crds,
			   VDouble& N)
{
    if (crds.size() != 3) {
	return;
    }
    if (coeff.size() != 4) {
	return;
    }
    
    N.resize(4,0.0);
    VDouble col(4);
    col[0] = 1.0;
    for (unsigned int i=0; i<crds.size(); ++i) {
	col[i+1] = crds[i];
    }

    for (unsigned int i=0; i<coeff.size(); ++i) {
	if (coeff[i].size() != 4) {
	    return;
	}
	N[i] = dotVDouble(coeff[i],col);
	if (fabs(N[i]) < tol) {
	    N[i] = tol;
	}
    }
}

void
BackgroundMesh::getNForRect(double x0, double y0, double hx, double hy,
			    double x, double y, VDouble& N)
{
    // compute local coordinate of the particle
    double xl = (x-x0)/hx;
    double yl = (y-y0)/hy;

    // map to [-1, 1]
    xl = xl*2-1;
    yl = yl*2-1;

    // shape function
    N.resize(4);
    N[0] = (1-xl)*(1-yl)/4.0;
    N[1] = (1+xl)*(1-yl)/4.0;
    N[2] = (1+xl)*(1+yl)/4.0;
    N[3] = (1-xl)*(1+yl)/4.0;
}

void
BackgroundMesh::getNForRect(double x0, double y0, double z0,
			    double hx, double hy, double hz,
			    double x, double y, double z,
			    VDouble& N)
{
    // compute local coordinate of the particle
    double xl = (x-x0)/hx;
    double yl = (y-y0)/hy;
    double zl = (z-z0)/hz;

    // map to [-1, 1]
    xl = xl*2-1;
    yl = yl*2-1;
    zl = zl*2-1;

    // shape function
    N.resize(8);
    N[0] = (1-xl)*(1-yl)*(1-zl)/8.0;
    N[1] = (1+xl)*(1-yl)*(1-zl)/8.0;
    N[2] = (1+xl)*(1+yl)*(1-zl)/8.0;
    N[3] = (1-xl)*(1+yl)*(1-zl)/8.0;
    N[4] = (1-xl)*(1-yl)*(1+zl)/8.0;
    N[5] = (1+xl)*(1-yl)*(1+zl)/8.0;
    N[6] = (1+xl)*(1+yl)*(1+zl)/8.0;
    N[7] = (1-xl)*(1+yl)*(1+zl)/8.0;
}

int
BackgroundMesh::clearBackground()
{

    // remove elements
    clearGridEles();

    // remove cells
    clearGrid();

    return 0;
}

void
BackgroundMesh::clearGridEles()
{
    // remove elements
    TaggedObjectIter& meshes = OPS_getAllMesh();
    Mesh* mesh = 0;
    while ((mesh = dynamic_cast<Mesh*>(meshes())) != 0) {
	ParticleGroup* group = dynamic_cast<ParticleGroup*>(mesh);
	if (group == 0) {
	    continue;
	}

	// remove elements
	group->clearEles();
    }
}

void
BackgroundMesh::clearGrid()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;

    // remove cells
    std::map<VInt,BNode> fixedbnodes;
    for (std::map<VInt,BNode>::iterator it=bnodes.begin(); it!=bnodes.end(); ++it) {
	BNode& bnode = it->second;
	const VInt& tags = bnode.tags;
	const VInt& types = bnode.type;

	BNode newbnode;
	for (unsigned int i=0; i<types.size(); ++i) {
	    if (types[i] == FLUID) {
		// remove node
		Node* nd = domain->removeNode(tags[i]);
		if (nd != 0) {
		    delete nd;
		}
		
		// remove pc
		Pressure_Constraint* pc = domain->removePressure_Constraint(tags[i]);
		if (pc != 0) {
		    delete pc;
		}
	    } else if (types[i] == FIXED) {
		newbnode.tags.push_back(tags[i]);
		newbnode.crdsn.push_back(bnode.crdsn[i]);
		newbnode.vn.push_back(bnode.vn[i]);
		newbnode.dvn.push_back(bnode.dvn[i]);
		newbnode.pn.push_back(bnode.pn[i]);
		newbnode.dpn.push_back(bnode.dpn[i]);
		newbnode.type.push_back(types[i]);
	    }
	}

	if (newbnode.tags.empty() == false) {
	    fixedbnodes[it->first] = newbnode;	    
	}
    }

    bnodes.clear();
    bcells.clear();
    bnodes = fixedbnodes;
}

bool
BackgroundMesh::inEle(const VDouble& N)
{
    // out
    for (unsigned int j=0; j<N.size(); ++j) {

	if (N[j]<0) {
	    // j+1
	    return false;
	}
    }

    return true;
}

int
BackgroundMesh::solveLine(const VDouble& p1, const VDouble& dir,
			  int dim, double crd, double& k)
{
    // check
    if (p1.size()!=dir.size()) {
	std::cerr << "WARNING: sizes are not compatible -- BgMesh::solveLine\n";
	return -1;
    }
    if (dim<0 || dim>=(int)dir.size()) {
	std::cerr << "WARNING: dim is out of range -- BgMesh::solveLine\n";
	return -1;
    }

    // solve k
    if (dir[dim] == 0.0) {
	k = -1.0;
    } else {
	k = (crd-p1[dim])/dir[dim];
    }

    return 0;
}

int
BackgroundMesh::remesh(bool init)
{
    // clear and check
    if (bsize <= 0.0) {
	std::cerr << "WARNING: basic mesh size has not been set -- BgMesh::addParticles\n";
	return -1;
    }

#ifdef _LINUX
    Timer timer;
    timer.start();
#endif

    // move particles
    if (moveParticles() < 0) {
	std::cerr << "WARNING: failed to move particles\n";
	return -1;
    }

#ifdef _LINUX
    timer.pause();
    std::cout<<"time for move particles = "<<timer.getReal()<<"\n";
    timer.start();
#endif

    // clear background
    clearBackground();

    // add structure
    if (addStructure() < 0) {
	std::cerr << "WARNING: failed to add structure\n";
	return -1;
    }

#ifdef _LINUX
    timer.pause();
    std::cout<<"time for add structure = "<<timer.getReal()<<"\n";
    timer.start();
#endif
    // add particles
    if (addParticles(init) < 0) {
	std::cerr << "WARNING: failed to add particles\n";
	return -1;
    }

#ifdef _LINUX
    timer.pause();
    std::cout<<"time for add particles = "<<timer.getReal()<<"\n";
    timer.start();
#endif

    // create grid nodes
    if (gridNodes() < 0) {
	std::cerr << "WARNING: failed to create grid nodes\n";
	return -1;
    }

#ifdef _LINUX
    timer.pause();
    std::cout<<"time for grid nodes = "<<timer.getReal()<<"\n";
    timer.start();
#endif

    // move particles in fixed cells
    if (moveFixedParticles()) {
    	opserr << "WARNING: failed to move particles in fixed cells";
    	return -1;
    }

#ifdef _LINUX
    timer.pause();
    std::cout<<"time for moving fixed particles = "<<timer.getReal()<<"\n";
    timer.start();
#endif

    // create grid elements
    if (gridFluid() < 0) {
	std::cerr << "WARNING: failed to create fluid elements\n";
	return -1;
    }

#ifdef _LINUX
    timer.pause();
    std::cout<<"time for fluid eles = "<<timer.getReal()<<"\n";
    timer.start();
#endif

    // create FSI elements
    if (gridFSI() < 0) {
    	std::cerr << "WARNING: failed to create FSI elements\n";
    	return -1;
    }

#ifdef _LINUX
    timer.pause();
    std::cout<<"time for fsi eles = "<<timer.getReal()<<"\n";
    timer.start();
#endif

    // if (freeSurface() < 0) {
    // 	opserr << "WARNING: failed to add pressures on free surface\n";
    // 	return -1;
    // }

    // timer.pause();
    // std::cout<<"time for free surface = "<<timer.getReal()<<"\n";
    // timer.start();

    if (record(init) < 0) {
	std::cerr << "WARNING: failed to record\n";
	return -1;
    }

#ifdef _LINUX
    timer.pause();
    std::cout<<"time for recording = "<<timer.getReal()<<"\n";
    timer.start();
#endif

    return 0;
}

// for each fixed node
// the list is unique
// get current states
// get nearest grid
// check if on that grid, if not skip it
// check if already a node, error
// add to bnode with type = FIXED
int
BackgroundMesh::addFixed(const VInt& ndtags)
{
    // get domain
    int ndm = OPS_GetNDM();
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    // add node
    int ndtag = Mesh::nextNodeTag();
    for (unsigned int i=0; i<ndtags.size(); ++i) {

	// if already there
	if (fixedNodes.find(ndtags[i]) != fixedNodes.end()) continue;

	// get node
	Node* nd = domain->getNode(ndtags[i]);
	if (nd == 0) continue;

	// nodal data
	const Vector& crds = nd->getCrds();
    	const Vector& disp = nd->getTrialDisp();
    	const Vector& vel = nd->getTrialVel();
	const Vector& accel = nd->getTrialAccel();

	if (crds.Size()!=ndm || disp.Size()<ndm) {
	    continue;
	}

	// nodal crds
	VDouble crdsn(ndm), vn(ndm), dvn(ndm);
	for (int j=0; j<ndm; ++j) {
	    crdsn[j] = crds(j)+disp(j);
	    vn[j] = vel(j);
	    dvn[j] = accel(j);
	}

	// near index
	VInt index;
	nearIndex(crdsn, index);

	// if not on a grid point
	VDouble gridcrds;
	getCrds(index,gridcrds);
	gridcrds -= crdsn;
	if (normVDouble(gridcrds) > meshtol*bsize) continue;

	// create pressure constraint
	Pressure_Constraint* pc = domain->getPressure_Constraint(nd->getTag());
	if(pc != 0) {
	    pc->setDomain(domain);
	} else {

	    // create pressure node
	    Node* pnode = 0;
	    if (ndm == 2) {
		pnode = new Node(ndtag++, 1, crds[0], crds[1]);
	    } else if (ndm == 3) {
		pnode = new Node(ndtag++, 1, crds[0], crds[1], crds[2]);
	    }
	    if (pnode == 0) {
		std::cerr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
		return -1;
	    }
	    if (domain->addNode(pnode) == false) {
		std::cerr << "WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
		delete pnode;
		return -1;
	    }

	    pc = new Pressure_Constraint(nd->getTag(), pnode->getTag());
	    if(pc == 0) {
		std::cerr<<"WARNING: no enough memory for Pressure_Constraint\n";
		return -1;
	    }
	    if (domain->addPressure_Constraint(pc) == false) {
		std::cerr << "WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
		delete pc;
		return -1;
	    }
	}

	// nodal pressures
	double pressure = pc->getPressure();
	double pdot = pc->getPdot();

	// create bnode
	BNode& bnode = bnodes[index];
	if (bnode.tags.empty() == false) {
	    opserr << "WARNING: two fixed nodes are at same location -- addFixed\n";
	    return -1;
	}

	bnode.addNode(nd->getTag(),crdsn,vn,dvn,pressure,pdot,FIXED);
	
	fixedNodes.insert(ndtags[i]);
    }

    return 0;
}

// for each structural node
// if it's in fixed, skip it
// get current states
// find nearest bnode
// if bnode is fixed, check if too close to it
// if bnode is empty, clear it
// add structure node to bnodes
// get min and max of FSI area
// set cells close to the structures as STRUCTURE
// set grids close to the structures as EMPTY unless it's FIXED or STRUCTURE
// set FSI area bnodes and cells
int
BackgroundMesh::addStructure()
{
    // get domain
    int ndm = OPS_GetNDM();
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    
    // add all structural nodes to the background
    int ndtag = Mesh::nextNodeTag();
    for (std::set<int>::iterator it=structuralNodes.begin();
	 it!=structuralNodes.end(); ++it) {

	// if fixed
	if (fixedNodes.find(*it) != fixedNodes.end()) continue;

	// get node
	Node* nd = domain->getNode(*it);
	if (nd == 0) continue;

	// nodal data
	const Vector& crds = nd->getCrds();
    	const Vector& disp = nd->getTrialDisp();
    	const Vector& vel = nd->getTrialVel();
	const Vector& accel = nd->getTrialAccel();

	if (crds.Size()!=ndm || disp.Size()<ndm) {
	    continue;
	}

	// create pressure constraint
	Pressure_Constraint* pc = domain->getPressure_Constraint(nd->getTag());
	if(pc != 0) {
	    pc->setDomain(domain);
	} else {

	    // create pressure node
	    Node* pnode = 0;
	    if (ndm == 2) {
		pnode = new Node(ndtag++, 1, crds[0], crds[1]);
	    } else if (ndm == 3) {
		pnode = new Node(ndtag++, 1, crds[0], crds[1], crds[2]);
	    }
	    if (pnode == 0) {
		std::cerr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
		return -1;
	    }
	    if (domain->addNode(pnode) == false) {
		std::cerr << "WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
		delete pnode;
		return -1;
	    }

	    pc = new Pressure_Constraint(nd->getTag(), pnode->getTag());
	    if(pc == 0) {
		std::cerr<<"WARNING: no enough memory for Pressure_Constraint\n";
		return -1;
	    }
	    if (domain->addPressure_Constraint(pc) == false) {
		std::cerr << "WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
		delete pc;
		return -1;
	    }
	}

	// nodal pressures
	double pressure = pc->getPressure();
	double pdot = pc->getPdot();

	// nodal crds
	VDouble crdsn(ndm), vn(ndm), dvn(ndm);
	for (int i=0; i<ndm; ++i) {
	    crdsn[i] = crds(i)+disp(i);
	    vn[i] = vel(i);
	    dvn[i] = accel(i);
	}

	// near index
	VInt index;
	nearIndex(crdsn, index);

	// add structural node to the bnode
	BNode& bnode = bnodes[index];

	// check with fixed bnodes and empty bnodes
	bool ignore = false;
	for (unsigned int i=0; i<bnode.type.size(); ++i) {
	    if (bnode.type[i] == FIXED) {
		// ignore the structural node
		ignore = true;
		break;
		
	    } else if (bnode.type[i] == EMPTY) {
		bnode.clear();
		break;
		
	    } else if (bnode.type[i] == FLUID) {
		bnode.clear();
		break;
	    }
	}

	// ignore structural node
	if (ignore) continue;

	bnode.addNode(nd->getTag(),crdsn,vn,dvn,pressure,pdot,STRUCTURE);

	// set empty bnodes
	VInt ind = index;
	ind -= 1;
	VVInt indices;
	getCorners(ind, 2, indices);
	for (unsigned int i=0; i<indices.size(); ++i) {
	    BNode& bnd = bnodes[indices[i]];
	    if (bnd.size() == 0) {
		bnd.addNode(EMPTY);
	    } else if (bnd.type[0] == FLUID) {
		bnd.clear();
		bnd.addNode(EMPTY);
	    }
	}

	// set STRUCTURE cells
	getCorners(ind, 1, indices);
	for (unsigned int i=0; i<indices.size(); ++i) {
	    BCell& bcell = bcells[indices[i]];
	    bcell.type = STRUCTURE;
	    
	    // set corners
	    if (bcell.bnodes.empty()) {

		VVInt corners;
		getCorners(indices[i], 1, corners);
		
		for (unsigned int j=0; j<corners.size(); ++j) {
		    BNode& bnode = bnodes[corners[j]];
		    bcell.bnodes.push_back(&bnode);
		    bcell.bindex.push_back(corners[j]);
		}
	    }
	}

	// set FSI bnodes
	ind -= 1;
	getCorners(ind, 4, indices);
	for (unsigned int i=0; i<indices.size(); ++i) {
	    BNode& bnd = bnodes[indices[i]];
	    if (bnd.size() == 0) {
		bnd.addNode(FLUID);
	    }
	}

	// set STRUCTURE cells
	getCorners(ind, 3, indices);
	for (unsigned int i=0; i<indices.size(); ++i) {
	    BCell& bcell = bcells[indices[i]];
	    if (bcell.type == FLUID) {
		bcell.type = FSI;
	    }
	    
	    // set corners
	    if (bcell.bnodes.empty()) {

		VVInt corners;
		getCorners(indices[i], 1, corners);
		
		for (unsigned int j=0; j<corners.size(); ++j) {
		    BNode& bnode = bnodes[corners[j]];
		    bcell.bnodes.push_back(&bnode);
		    bcell.bindex.push_back(corners[j]);
		}
	    }
	}
    }
    
    return 0;
}

// for each particle
// find lower grid
// check if out of range, if out, remove it
// get the containing cell
// set bnodes of the cell
// add the particle to the cell
// add bnodes to the cell
int
BackgroundMesh::addParticles(bool init)
{
    // for all particles
    TaggedObjectIter& meshes = OPS_getAllMesh();
    Mesh* mesh = 0;
    while((mesh = dynamic_cast<Mesh*>(meshes())) != 0) {
	ParticleGroup* group = dynamic_cast<ParticleGroup*>(mesh);
	if (group == 0) {
	    continue;
	}

	// remove particles
	VInt rm(group->numParticles(), 0);

	// for all particles
	for (int j=0; j<group->numParticles(); j++) {

	    // get particle
	    Particle* p = group->getParticle(j);
	    if (p == 0) continue;

	    // get particle coordinates
	    const VDouble& crds = p->getCrds();

	    // get index
	    VInt index;
	    lowerIndex(crds, index);

	    // if out of range
	    for (unsigned int i=0; i<index.size(); ++i) {
		if (index[i]<lower[i] || index[i]>=upper[i]) {
		    rm[j] = 1;
		    break;
		}
	    }
	    if (rm[j] == 1) continue;

	    // get bcell
	    BCell& bcell = bcells[index];

	    // if initial check structure cell
	    if (bcell.type == STRUCTURE) {
		rm[j] = 1;
		continue;
	    }

	    // add particles
	    bcell.add(p);

	    // add bnodes of the cell
	    if (bcell.bnodes.empty()) {

		// get corners
		VVInt indices;
		getCorners(index,1,indices);
		
		// set corners
  		for (unsigned int i=0; i<indices.size(); ++i) {
		    BNode& bnode = bnodes[indices[i]];
		    if (bnode.size() == 0) {
			bnode.addNode(FLUID);
		    }
		    bcell.bnodes.push_back(&bnode);
		    bcell.bindex.push_back(indices[i]);
		}
	    }

	}

	// remove out of range particles
	group->removeParticles(rm);
    }

    return 0;
}



int
BackgroundMesh::gridNodes()
{
    // get domain
    int ndm = OPS_GetNDM();
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    // vector of iterators
    std::vector<std::map<VInt,BNode>::iterator> iters;
    iters.reserve(bnodes.size());
    for (std::map<VInt,BNode>::iterator it=bnodes.begin(); it!=bnodes.end(); ++it) {
	iters.push_back(it);
    }

    // each cell
    int ndtag = Mesh::nextNodeTag();
    std::vector<Node*> newnodes(iters.size(),0), newpnodes(iters.size(), 0);
    std::vector<Pressure_Constraint*> newpcs(iters.size(), 0);

    int res = 0;

#pragma omp parallel for
    for (unsigned int j=0; j<iters.size(); ++j) {

	// get iterator
	std::map<VInt,BNode>::iterator it = iters[j];

	// get cell
	const VInt& index = it->first;
	BNode& bnode = it->second;
	if (bnode.size() == 0) {
	    opserr << "WARNING: bnode.size() = 0 -- gridNodes\n";
	    continue;
	}
	if (bnode.size() > 1) continue;
	if (bnode.type[0] != FLUID) continue;

	// coordinates
	VDouble crds;
	getCrds(index,crds);

	// get particles
	VParticle pts;
	VInt minind = index;
	VInt maxind = index;
	minind -= numave;
	maxind += numave;
	gatherParticles(minind,maxind,pts);

	// get information
	double wt = 0.0, pre = 0.0, pdot = 0.0;
	VDouble vel(ndm), accel(ndm);
	for (unsigned int i=0; i<pts.size(); ++i) {

	    // get particle
	    if (pts[i] == 0) continue;

	    // particle coordinates
	    const VDouble& pcrds = pts[i]->getCrds();

	    // distance from particle to current location
	    VDouble dist = pcrds;
	    dist -= crds;
	    double q = normVDouble(dist) / (bsize);

	    // weight for the particle
	    double w = QuinticKernel(q, bsize, ndm);

	    // add weight
	    wt += w;

	    // add pressure
	    pre +=  pts[i]->getPressure() * w;
	    pdot += pts[i]->getPdot() * w;

	    // add velocity
	    const VDouble& pvel = pts[i]->getVel();
	    for (int k=0; k<ndm; k++) {
		vel[k] += w*pvel[k];
	    }

	    // add acceleration
	    const VDouble& paccel = pts[i]->getAccel();
	    for (int k=0; k<ndm; k++) {
		accel[k] += w*paccel[k];
	    }
	}
	if (wt > 0) {
	    pre /= wt;
	    pdot /= wt;
	    vel /= wt;
	    accel /= wt;
	}

	// create node
	Node* node = 0;
	if (ndm == 2) {
	    node = new Node(ndtag+2*j, ndm, crds[0], crds[1]);
	} else if (ndm == 3) {
	    node = new Node(ndtag+2*j, ndm, crds[0], crds[1], crds[2]);
	}
	if (node == 0) {
	    std::cerr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
	    res = -1;
	    continue;
	}

	if (wt > 0) {
	    Vector vvel;
	    toVector(vel, vvel);
	    Vector vaccel;
	    toVector(accel, vaccel);
	    node->setTrialVel(vvel);
	    node->setTrialAccel(vaccel);
	    node->commitState();
	}

	// add to newnodes
	newnodes[j] = node;

	// set the bnode
	bnode.setNode(0,node->getTag(),crds,vel,accel,pre,pdot,FLUID);

	// set pressure
	Pressure_Constraint* thePC = domain->getPressure_Constraint(node->getTag());
	Node* pnode = 0;
	if(thePC != 0) {
	    thePC->setDomain(domain);
	    pnode = thePC->getPressureNode();
	    if (pnode == 0) {
		std::cerr << "WARNING: pressure does not exist -- BgMesh::gridNodes\n";
		res = -1;
		continue;
	    }
	} else {

	    // create pressure node
	    if (ndm == 2) {
		pnode = new Node(ndtag+2*j+1, 1, crds[0], crds[1]);
	    } else if (ndm == 3) {
		pnode = new Node(ndtag+2*j+1, 1, crds[0], crds[1], crds[2]);
	    }
	    if (pnode == 0) {
		std::cerr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
		res = -1;
		continue;
	    }
	    newpnodes[j] = pnode;

	    thePC = new Pressure_Constraint(node->getTag(), pnode->getTag());
	    if(thePC == 0) {
		std::cerr<<"WARNING: no enough memory for Pressure_Constraint\n";
		res = -1;
		continue;
	    }
	    newpcs[j] = thePC;

	}
	if (wt > 0) {
	    Vector newvel = pnode->getVel();
	    newvel.Zero();
	    newvel(0) = pre;
	    Vector newaccel = pnode->getAccel();
	    newaccel.Zero();
	    newaccel(0) = pdot;
	    pnode->setTrialVel(newvel);
	    pnode->setTrialAccel(newaccel);
	    pnode->commitState();
	}

	
    }

    if (res < 0) {
	return -1;
    }

    // add nodes and pcs to domain
    for (unsigned int i=0; i<newnodes.size(); ++i) {
	if (newnodes[i] == 0) continue;

	// add to domain
	if (domain->addNode(newnodes[i]) == false) {
	    std::cerr<<"WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
	    delete newnodes[i];
	    return -1;
	}
    }
    for (unsigned int i=0; i<newpnodes.size(); ++i) {
	if (newpnodes[i] == 0) continue;

	// add to domain
	if (domain->addNode(newpnodes[i]) == false) {
	    std::cerr<<"WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
	    delete newpnodes[i];
	    return -1;
	}
    }
    for (unsigned int i=0; i<newpcs.size(); ++i) {
	if (newpcs[i] == 0) continue;

	// add to domain
	if(domain->addPressure_Constraint(newpcs[i]) == false) {
	    std::cerr<<"WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
	    delete newpcs[i];
	    return -1;
	}
    }

    return 0;
}

int
BackgroundMesh::moveFixedParticles()
{
    // check each cell
    for (std::map<VInt, BCell>::iterator it = bcells.begin(); it != bcells.end(); ++it) {

	// get cell
	const VInt& index = it->first;
	BCell& cell = it->second;

	// empty cell
	if (cell.pts.empty()) {
	    continue;
	}

	// check if STRUCTURE cell
	if (cell.type != STRUCTURE) {
	    continue;
	}

	// check neighbor cells
	VInt ind = index;
	VVInt indices;
	ind -= 1;
	getCorners(ind, 2, indices);

	// give each cell a score
	VInt scores(indices.size());
	for (unsigned int i=0; i<indices.size(); ++i) {
	    std::map<VInt,BCell>::iterator cellit = bcells.find(indices[i]);
	    if (cellit == bcells.end()) {
		scores[i] = 1;
		continue;
	    }
	    if (cellit->second.type == STRUCTURE) {
		scores[i] = -1;
	    } else {
		scores[i] = 1;
	    }
	}
	int cellmap[9][2] = {{1,3},{0,2},{1,5},{0,6},
			     {4,4},{2,8},{3,7},{6,8},
			     {5,7}};
	for (int i=0; i<9; ++i) {
	    if (scores[i] > 0) {
		scores[i] += scores[cellmap[i][0]];
		scores[i] += scores[cellmap[i][1]];
	    }
	}

	// find the cell with highest score
	int high = -1;
	for (unsigned int i=0; i<scores.size(); ++i) {
	    if (high < scores[i]) {
		high = scores[i];
		ind = indices[i];
	    }
	}
	if (high < 0) continue;
	
	// move the particles
	VDouble crds;
	getCrds(index,crds);
	for (int i = 0; i < (int)cell.pts.size(); ++i) {
	    Particle* pt = cell.pts[i];
	    const VDouble& pcrds = pt->getCrds();

	    // move the particle
	    VDouble newcrds;
	    getCrds(ind,newcrds);
	    VDouble disp = pcrds;
	    disp -= crds[0];
	    newcrds += disp;
	    pt->moveTo(newcrds,0.0);

	    // add particles to the new cell
	    std::map<VInt,BCell>::iterator cellit = bcells.find(ind);
	    if (cellit != bcells.end()) {
		cellit->second.pts.push_back(pt);
	    }
	}

	cell.pts.clear();
    }
    return 0;
}

int
BackgroundMesh::gridFluid()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();

    // store cells in a vector
    std::vector<BCell*> cells;
    VVInt indices;
    cells.reserve(bcells.size());
    indices.reserve(bcells.size());
    for (std::map<VInt,BCell>::iterator it=bcells.begin(); it!=bcells.end(); ++it) {
	indices.push_back(it->first);
	cells.push_back(&(it->second));
    }

    // create elements in each cell
    int numele = 0;
    int numelenodes = 0;
    if (ndm == 2) {
	numele = 2;
	numelenodes = 3;
    } else if (ndm == 3) {
	numele = 6;
	numelenodes = 4;
    }
    VVInt elends(numele*cells.size());
    VInt gtags(numele*cells.size());
#pragma omp parallel for
    for (unsigned int j=0; j<cells.size(); ++j) {

	// structural cell
	if (cells[j]->type == STRUCTURE ||
	    cells[j]->type == FSI) continue;

	// find the group of this mesh
	std::map<int,int> numpts;
	for (unsigned int i=0; i<cells[j]->pts.size(); ++i) {
	    numpts[cells[j]->pts[i]->getGroupTag()] += 1;
	}
	int num=0;
	int gtag=0;
	for (std::map<int,int>::iterator it=numpts.begin(); it!=numpts.end(); ++it) {
	    if (num < it->second) {
		num = it->second;
		gtag = it->first;
	    }
	}
	for (int i=0; i<numele; ++i) {
	    gtags[numele*j+i] = gtag;
	}

	// add to elenodes
	VInt cnodes(cells[j]->bnodes.size());
	for (unsigned int i=0; i<cnodes.size(); ++i) {
	    if (cells[j]->bnodes[i]->tags.empty()) {
		opserr << "WARNING: failed to be fluid node -- gridFluid\n";
		continue;
	    }
	    cnodes[i] = cells[j]->bnodes[i]->tags[0];
	}
	for (int i=0; i<numele; ++i) {
	    elends[numele*j+i].resize(numelenodes);
	}
	if (ndm == 2) {
	    elends[numele*j][0] = cnodes[0];
	    elends[numele*j][1] = cnodes[3];
	    elends[numele*j][2] = cnodes[2];

	    elends[numele*j+1][0] = cnodes[0];
	    elends[numele*j+1][1] = cnodes[1];
	    elends[numele*j+1][2] = cnodes[3];
	    
	} else if (ndm == 3) {
	    elends[numele*j][0] = cnodes[0];
	    elends[numele*j][1] = cnodes[7];
	    elends[numele*j][2] = cnodes[4];
	    elends[numele*j][3] = cnodes[5];

	    elends[numele*j+1][0] = cnodes[0];
	    elends[numele*j+1][1] = cnodes[1];
	    elends[numele*j+1][2] = cnodes[7];
	    elends[numele*j+1][3] = cnodes[5];

	    elends[numele*j+2][0] = cnodes[0];
	    elends[numele*j+2][1] = cnodes[1];
	    elends[numele*j+2][2] = cnodes[3];
	    elends[numele*j+2][3] = cnodes[7];

	    elends[numele*j+3][0] = cnodes[0];
	    elends[numele*j+3][1] = cnodes[3];
	    elends[numele*j+3][2] = cnodes[2];
	    elends[numele*j+3][3] = cnodes[7];

	    elends[numele*j+4][0] = cnodes[6];
	    elends[numele*j+4][1] = cnodes[0];
	    elends[numele*j+4][2] = cnodes[7];
	    elends[numele*j+4][3] = cnodes[4];

	    elends[numele*j+5][0] = cnodes[6];
	    elends[numele*j+5][1] = cnodes[0];
	    elends[numele*j+5][2] = cnodes[2];
	    elends[numele*j+5][3] = cnodes[7];
	}
	
    }

    // get particle group tags
    std::map<int,ID> elenodes;
    for (unsigned int i=0; i<elends.size(); ++i) {

	// no elenodes, no element
	if (elends[i].empty()) continue;

	// if all nodes are fluid and fixed nodes
	ID& nds = elenodes[gtags[i]];
	for (unsigned int j=0; j<elends[i].size(); ++j) {
	    nds[nds.Size()] = elends[i][j];
	}
    }

    // create elements
    for (std::map<int,ID>::iterator it=elenodes.begin(); it!=elenodes.end(); ++it) {
	ParticleGroup* group = dynamic_cast<ParticleGroup*>(OPS_getMesh(it->first));
	if (group == 0) {
	    std::cerr << "WARNING: failed to get particle group -- BgMesh::gridFluid\n";
	    return -1;
	}
	group->setEleNodes(it->second);

	if (group->newElements(it->second) < 0) {
	    std::cerr << "WARNING: failed to create elements for mesh ";
	    std::cerr << group->getTag()<<" -- BgMesh::gridFluid\n";
	    return -1;
	}
    }
    
    return 0;
}

int
BackgroundMesh::gridFSI()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();

    TriangleMeshGenerator gen;
    TetMeshGenerator tetgen;

    // gather bnodes
    std::map<VInt,BNode*> fsibnodes;
    for (std::map<VInt,BCell>::iterator it=bcells.begin(); it!=bcells.end(); ++it) {
	
	// only for structural and FSI cells
	BCell& bcell = it->second;	
	if (bcell.type == FLUID) continue;

	// get bnode
	for (unsigned int j=0; j<bcell.bnodes.size(); ++j) {
	    BNode* bnode = bcell.bnodes[j];
	    VInt bindex = bcell.bindex[j];
	    if (bnode == 0) {
		opserr << "WARNING: failed to get bnode -- gridFSI\n";
		return -1;
	    }
	    fsibnodes[bindex] = bnode;
	}
    }

    // add points
    VInt ndtags, ndtypes;
    VVInt ndindex;
    VDouble min, max;
    for (std::map<VInt,BNode*>::iterator it=fsibnodes.begin(); it!=fsibnodes.end(); ++it) {

	VInt bindex = it->first;
	BNode* bnode = it->second;
	
	VInt& tags = bnode->tags;
	VVDouble& crdsn = bnode->crdsn;
	VInt& type = bnode->type;

	for (unsigned int i=0; i<tags.size(); ++i) {
	    if (type[i] == EMPTY) continue;
	    
	    ndtags.push_back(tags[i]);
	    ndtypes.push_back(type[i]);
	    ndindex.push_back(bindex);
	    
	    if (ndm == 2) {
		gen.addPoint(crdsn[i][0],crdsn[i][1]);
	    } else if (ndm == 3) {
		tetgen.addPoint(crdsn[i][0],crdsn[i][1],crdsn[i][2],0);
	    }

	    // max, min coordinates
	    if (min.empty() || max.empty()) {
		min = crdsn[i];
		max = crdsn[i];
	    } else {
		for (unsigned int j=0; j<crdsn.size(); ++j) {
		    if (min[j] > crdsn[i][j]) min[j] = crdsn[i][j];
		    if (max[j] < crdsn[i][j]) max[j] = crdsn[i][j];
		}
	    }
	}
    }

    // no mesh
    int numpoints = 0;
    if (ndm == 2) {
	numpoints = gen.getNumPoints();
    } else if (ndm == 3) {
	numpoints = tetgen.getNumPoints();
    }
    if (gen.getNumPoints() < 3 &&
	tetgen.getNumPoints() < 4 ) {
	return 0;
    }
    if (min.empty() || max.empty()) {
	return 0;
    }

    // add extra points to avoid small triangles on the edge
    min -= bsize;
    max += bsize;
    if (ndm == 2) {
	gen.addPoint(min[0],min[1]);
	gen.addPoint(max[0],min[1]);
	gen.addPoint(min[0],max[1]);
	gen.addPoint(max[0],max[1]);
    } else if (ndm == 3) {
	tetgen.addPoint(min[0],min[1],min[2],0);
	tetgen.addPoint(min[0],max[1],min[2],0);
	tetgen.addPoint(min[0],max[1],max[2],0);
	tetgen.addPoint(min[0],min[1],max[2],0);
	tetgen.addPoint(max[0],min[1],min[2],0);
	tetgen.addPoint(max[0],max[1],min[2],0);
	tetgen.addPoint(max[0],max[1],max[2],0);
	tetgen.addPoint(max[0],min[1],max[2],0);
    }

    // mesh
    if (ndm == 2) {
	gen.remesh(-1.0);
    } else if (ndm == 3) {
	tetgen.remesh(-1.0);
    }

    // get triangles or tetrahedrons
    int numele = 0;
    if (ndm == 2) {
	numele = gen.getNumTriangles();
    } else if (ndm == 3) {
	numele = tetgen.getNumTets();
    }
    
    VVInt elends(numele);
    VInt gtags(numele);
#pragma omp parallel for
    for (int i=0; i<numele; ++i) {

	// get points
	VInt tri;
	if (ndm == 2) {
	    tri.resize(3);
	    gen.getTriangle(i,tri[0],tri[1],tri[2]);
	} else if (ndm == 3) {
	    tri.resize(4);
	    tetgen.getTet(i,tri[0],tri[1],tri[2],tri[3]);
	}

	// check if connect to extra points
	bool extra = false;
	for (unsigned int j=0; j<tri.size(); ++j) {
	    if (tri[j] >= numpoints) {
		extra = true;
		break;
	    }
	}
	if (extra) continue;

	// check if connect to fluid only
	bool fluid = true;
	for (unsigned int j=0; j<tri.size(); ++j) {
	    if (ndtypes[tri[j]] == STRUCTURE) {
		fluid = false;
		break;
	    }
	}

	// get min and max ind
	VInt maxind = ndindex[tri[0]];
	VInt minind = ndindex[tri[0]];
	for (int k=0; k<ndm; ++k) {
	    for (unsigned int j=1; j<tri.size(); ++j) {
		if (ndindex[tri[j]][k] < ndindex[tri[0]][k]) {
		    minind[k] = ndindex[tri[j]][k];
		} else if(ndindex[tri[j]][k] > ndindex[tri[0]][k]) {
		    maxind[k] = ndindex[tri[j]][k];
		}
	    }
	}

	// check if triangle out side of FSI
	bool outside = false;
	if (fluid) {
	    VInt currind(ndm);
	    if (ndm == 2) {
		for (int j=minind[0]; j<maxind[0]; ++j) {
		    for (int k=minind[1]; k<maxind[1]; ++k) {
			currind[0] = j;
			currind[1] = k;
			std::map<VInt,BCell>::iterator it = bcells.find(currind);
			if (it != bcells.end()) {
			    if (it->second.type == FLUID) {
				outside = true;
				break;
			    }
			}
			if (outside) break;
		    }
 		}
	    } else if (ndm == 3) {
		for (int j=minind[0]; j<maxind[0]; ++j) {
		    for (int k=minind[1]; k<maxind[1]; ++k) {
			for (int l=minind[2]; l<maxind[2]; ++l) {
			    currind[0] = j;
			    currind[1] = k;
			    currind[2] = l;
			    std::map<VInt,BCell>::iterator it = bcells.find(currind);
			    if (it != bcells.end()) {
				if (it->second.type == FLUID) {
				    outside = true;
				    break;
				}
			    }
			    if (outside) break;
			}
		    }
 		}
 	    }
 	}
	if (outside) continue;

	// gather particles
	VParticle tripts;
	minind -= 1;
	maxind += 1;
	gatherParticles(minind,maxind,tripts,true);
	if (tripts.empty()) continue;

	// find the group of this mesh
	std::map<int,int> numpts;
	for (unsigned int j=0; j<tripts.size(); ++j) {
	    numpts[tripts[j]->getGroupTag()] += 1;
	}
	int num=0;
	for (std::map<int,int>::iterator it=numpts.begin(); it!=numpts.end(); ++it) {
	    if (num < it->second) {
		num = it->second;
		gtags[i] = it->first;
	    }
	}

	// add to elenodes
	elends[i].resize(tri.size());
	for (unsigned int j=0; j<tri.size(); ++j) {
	    elends[i][j] = ndtags[tri[j]];
	}
    }

    // get particle group tags
    std::map<int,ID> elenodes;
    for (unsigned int i=0; i<elends.size(); ++i) {

	// no elenodes, no element
	if (elends[i].empty()) continue;

	// if all nodes are fluid and fixed nodes
	ID& nds = elenodes[gtags[i]];
	for (unsigned int j=0; j<elends[i].size(); ++j) {
	    nds[nds.Size()] = elends[i][j];
	}
    }

    // create elements
    for (std::map<int,ID>::iterator it=elenodes.begin(); it!=elenodes.end(); ++it) {
	ParticleGroup* group = dynamic_cast<ParticleGroup*>(OPS_getMesh(it->first));
	if (group == 0) {
	    std::cerr << "WARNING: failed to get particle group -- BgMesh::gridFSI\n";
	    return -1;
	}
	group->addEleNodes(it->second);

	if (group->newElements(it->second) < 0) {
	    std::cerr << "WARNING: failed to create elements for mesh ";
	    std::cerr << group->getTag()<<" -- BgMesh::gridFSI\n";
	    return -1;
	}
    }

    return 0;
}

int
BackgroundMesh::gridEles()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();

    TriangleMeshGenerator gen;
    TetMeshGenerator tetgen;

    // add points
    VInt ndtags;
    ndtags.reserve(bnodes.size()*1.05);
    for (std::map<VInt,BNode>::iterator it=bnodes.begin(); it!=bnodes.end(); ++it) {
	BNode& bnode = it->second;
	for (unsigned int i=0; i<bnode.crdsn.size(); ++i) {
	    ndtags.push_back(bnode.tags[i]);
	    if (ndm == 2) {
		gen.addPoint(bnode.crdsn[i][0],bnode.crdsn[i][1]);
	    } else if (ndm == 3) {
		tetgen.addPoint(bnode.crdsn[i][0],bnode.crdsn[i][1],bnode.crdsn[i][2],0);
	    }
	}
    }

    // mesh
    if (ndm == 2) {
	gen.remesh(-1.0);
    } else if (ndm == 3) {
	tetgen.remesh(-1.0);
    }

    // get triangles or tetrahedrons
    int numele = 0;
    if (ndm == 2) {
	numele = gen.getNumTriangles();
    } else if (ndm == 3) {
	numele = tetgen.getNumTets();
    }
    
    VVInt elends(numele);
    VInt gtags(numele);
#pragma omp parallel for
    for (int i=0; i<numele; ++i) {

	// get points
	VInt tri;
	if (ndm == 2) {
	    tri.resize(3);
	    gen.getTriangle(i,tri[0],tri[1],tri[2]);
	} else if (ndm == 3) {
	    tri.resize(4);
	    tetgen.getTet(i,tri[0],tri[1],tri[2],tri[3]);
	}

	// get point crds 
	VVDouble ptcrds(tri.size());
	for (unsigned int j=0; j<tri.size(); ++j) {
	    ptcrds[j].resize(ndm);
	    int mark;
	    if (ndm == 2) {
		gen.getPoint(tri[j], ptcrds[j][0], ptcrds[j][1], mark);
	    } else if (ndm == 3) {
		tetgen.getPoint(tri[j], ptcrds[j][0], ptcrds[j][1], ptcrds[j][2], mark);
	    }
	}

	// get index for points
	// if the element is too large
	bool large = false;
	VVInt indices(tri.size());
	for (unsigned int j=0; j<indices.size(); ++j) {
	    nearIndex(ptcrds[j], indices[j]);
	}
	for (int ii=0; ii<(int)indices.size(); ++ii) {
	    for (int j=0; j<(int)indices.size()-1; ++j) {
		for (int k=0; k<ndm; ++k) {
		    if (fabs(indices[ii][k]-indices[j][k]) > 3) {
			large = true;
			break;
		    }
		}
		if (large) break;
	    }
	    if (large) break;
	}
	if (large) {
	    continue;
	}

	// precalculate shape functions
	VDouble coeff;
	VVDouble tetcoeff;
	bool zerovol = false;
	if (ndm == 2) {
	    if (preNForTri(ptcrds[0][0],ptcrds[0][1],
			   ptcrds[1][0],ptcrds[1][1],
			   ptcrds[2][0],ptcrds[2][1],
			   coeff) < 0) {
		zerovol = true;
	    }
	} else if (ndm == 3) {
	    if (preNForTet(ptcrds[0],ptcrds[1],ptcrds[2],ptcrds[3],
			   tetcoeff) < 0) {
		std::cout<<ptcrds[0];
		std::cout<<ptcrds[1];
		std::cout<<ptcrds[2];
		std::cout<<ptcrds[3];
		zerovol = true;
	    }
	}

	// zero volumne, no element
	if (zerovol) {
	    continue;
	}

	// get index range and point crds
	VInt minind, maxind;
	for (unsigned int j=0; j<tri.size(); ++j) {
	    VInt low, up;
	    lowerIndex(ptcrds[j], low);
	    upperIndex(ptcrds[j], up);
	    if (minind.empty()) minind = low;
	    if (maxind.empty()) maxind = up;
	    for (int k=0; k<ndm; ++k) {
		if (minind[k] > low[k]) {
		    minind[k] = low[k];
		}
		if (maxind[k] < up[k]) {
		    maxind[k] = up[k];
		}
	    }
	}

	// get particles
	VParticle pts;
	gatherParticles(minind,maxind,pts);

	// check which particles are in the triangle
	VParticle tripts;
	for (unsigned int j=0; j<pts.size(); ++j) {

	    // get shape function
	    const VDouble& pcrds = pts[j]->getCrds();
	    VDouble N;
	    if (ndm == 2) {
		getNForTri(coeff,pcrds[0],pcrds[1],N);
	    } else if (ndm == 3) {
		getNForTet(tetcoeff,pcrds,N);
	    }
	    if (inEle(N)) {
		// in tri
		tripts.push_back(pts[j]);
	    }
	}

	// no particles, no element
	if (tripts.empty()) continue;

	// find the group of this mesh
	std::map<int,int> numpts;
	for (unsigned int j=0; j<tripts.size(); ++j) {
	    numpts[tripts[j]->getGroupTag()] += 1;
	}
	int num=0;
	for (std::map<int,int>::iterator it=numpts.begin(); it!=numpts.end(); ++it) {
	    if (num < it->second) {
		num = it->second;
		gtags[i] = it->first;
	    }
	}

	// add to elenodes
	elends[i].resize(tri.size());
	for (unsigned int j=0; j<tri.size(); ++j) {
	    elends[i][j] = ndtags[tri[j]];
	}
    }

    // get particle group tags
    std::map<int,ID> elenodes;
    for (unsigned int i=0; i<elends.size(); ++i) {

	// no elenodes, no element
	if (elends[i].empty()) continue;

	// if all nodes are fluid and fixed nodes
	ID& nds = elenodes[gtags[i]];
	for (unsigned int j=0; j<elends[i].size(); ++j) {
	    nds[nds.Size()] = elends[i][j];
	}
    }

    // create elements
    for (std::map<int,ID>::iterator it=elenodes.begin(); it!=elenodes.end(); ++it) {
	ParticleGroup* group = dynamic_cast<ParticleGroup*>(OPS_getMesh(it->first));
	if (group == 0) {
	    std::cerr << "WARNING: failed to get particle group -- BgMesh::gridEles\n";
	    return -1;
	}
	group->setEleNodes(it->second);

	if (group->newElements(it->second) < 0) {
	    std::cerr << "WARNING: failed to create elements for mesh ";
	    std::cerr << group->getTag()<<" -- BgMesh::gridEles\n";
	    return -1;
	}
    }

    return 0;
}

// time, vx, vy, (vz), top, bottom, left, right, dt, numIter
// 

int
BackgroundMesh::record(bool init)
{

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();

    double timestamp = domain->getCurrentTime();
    currentTime = timestamp;

    // record
    for (unsigned int i=0; i<recorders.size(); ++i) {
	if (recorders[i] != 0) {
	    if (recorders[i]->record(domain->getCommitTag(),currentTime)) {
		std::cerr << "WARNING: failed to record -- BgMesh::gridEles\n";
		return -1;
	    }
	}
    }
    if (init) {
	return 0;
    }

    // record time
    if(theFile.good() == false) {
	return 0;
    }
    theFile << timestamp << " ";

    // record wave height and velocity
    for (unsigned int i=0; i<locs.size(); i+=2) {

	// lower index
	VDouble crds(2);
	crds[0] = locs[i];
	crds[1] = locs[i+1];

	VInt index;
	this->lowerIndex(crds,index);

	// shape function
	VDouble N;
	this->getNForRect(index[0]*bsize,index[1]*bsize,
			  bsize,bsize,
			  crds[0],crds[1],N);

	// velocity
	VDouble vel(ndm);

	// get corners
	VVInt indices;
	getCorners(index,1,indices);
	for (unsigned int i=0; i<indices.size(); ++i) {

	    // get crds
	    getCrds(indices[i], crds);

	    // check bnode
	    std::map<VInt,BNode>::iterator it = bnodes.find(indices[i]);
	    if (it == bnodes.end()) continue;

	    // get bnode
	    BNode& bnode = it->second;
	    if (bnode.vn.empty()) {
		continue;
	    }

	    // get vel
	    for (int j = 0; j < ndm; ++j) {
		vel[j] += N[i]*bnode.vn[0][j];
	    }
	}

	// get the lowest cell
	double bottom = 0.0;
	for (int iy=lower[1]; iy<=upper[1]; ++iy) {
	    VInt ind = index;
	    ind[1] = iy;
	    if (bcells.find(ind) != bcells.end()) {
		bottom = iy*bsize;
		break;
	    }
	}

	// get the most left cell
	double left = 0.0;
	for (int ix=lower[0]; ix<=upper[0]; ++ix) {
	    VInt ind = index;
	    ind[0] = ix;
	    if (bcells.find(ind) != bcells.end()) {
		left = ix*bsize;
		break;
	    }
	}
	// get the highest cell
	double top = 0.0;
	for (int iy=upper[1]; iy>=lower[1]; --iy) {
	    VInt ind = index;
	    ind[1] = iy;
	    if (bcells.find(ind) != bcells.end()) {
		top = (iy+1)*bsize;
		break;
	    }
	}

	// get the most right cell
	double right = 0.0;
	for (int ix=upper[0]; ix>=lower[0]; --ix) {
	    VInt ind = index;
	    ind[0] = ix;
	    if (bcells.find(ind) != bcells.end()) {
		right = (ix+1)*bsize;
		break;
	    }
	}

	// record velocity and wave
	for (int j=0; j<ndm; ++j) {
	    theFile << vel[j] << " ";
	}
	theFile << top << " " << bottom << " ";
	theFile << left << " " << right << " ";
    }

    // time step and iteration
    //theFile << ops_Dt << " " << OPS_numIter()<<"\n";
    theFile << ops_Dt << "\n";
    theFile.flush();

    return 0;
}

int
BackgroundMesh::moveParticles()
{
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();
    double dt = domain->getCurrentTime() - currentTime;

    // get current disp and velocity
    for (std::map<VInt,BNode>::iterator it=bnodes.begin(); it!=bnodes.end(); ++it) {
	BNode& bnode = it->second;
	VInt& tags = bnode.tags;

	for (int i=0; i<bnode.size(); ++i) {
	    Node* nd = domain->getNode(tags[i]);
	    Pressure_Constraint* pc = domain->getPressure_Constraint(tags[i]);

	    if (pc != 0) {
		bnode.pn[i] = pc->getPressure();
		bnode.dpn[i] = pc->getPdot();
	    }
	    if (nd != 0) {
		const Vector& vel = nd->getTrialVel();
		const Vector& accel = nd->getTrialAccel();
		for (int j=0; j<ndm; ++j) {
		    bnode.vn[i][j] = vel(j);
		    bnode.dvn[i][j] = accel(j);
		}
	    }
	}
    }

    // store cells in a vector
    std::vector<BCell*> cells;
    VVInt indices;
    cells.reserve(bcells.size());
    indices.reserve(bcells.size());
    for (std::map<VInt,BCell>::iterator it=bcells.begin(); it!=bcells.end(); ++it) {
	indices.push_back(it->first);
	cells.push_back(&(it->second));
    }

    // move particles in each cell
    int res = 0;
#pragma omp parallel for
    for (unsigned int j=0; j<cells.size(); ++j) {

	// get particles in cell
	const VParticle& pts = cells[j]->pts;

	// move the particle
	for (unsigned int i=0; i<pts.size(); ++i) {

	    // set update state
	    if (pts[i] == 0) continue;
	    pts[i]->needUpdate(dt);

	    // convect the particle
	    if (convectParticle(pts[i],indices[j],numsub) < 0) {
		std::cerr << "WARNING: failed to convect particle";
		std::cerr << " -- BgMesh::moveParticles\n";
		res = -1;
		continue;
	    }
	}
    }

    if (res < 0) return -1;

    return 0;
}

int
BackgroundMesh::convectParticle(Particle* pt, VInt index, int nums)
{

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    int ndm = OPS_GetNDM();

    // check dt
    double dt = pt->getDt();
    if (dt <= 0) {
	return 0;
    }

    // get corners
    VVInt indices;
    VVDouble crds;
    VInt fixed;
    VVDouble vels, dvns;
    VDouble pns, dpns;

    // convect in a cell
    double subdt = dt / nums;
    for (int n=0; n<nums; ++n) {

	// particle crds
	const VDouble& pcrds = pt->getCrds();

	// new index
	VInt newIndex;
	lowerIndex(pcrds,newIndex);

	// update corners
	if (n==0 || newIndex != index) {
	    index = newIndex;
	    
	    VVInt temp;
	    getCorners(index,1,temp);
	    indices = temp;
	    indices[2] = temp[3];
	    indices[3] = temp[2];
	    if (ndm == 3) {
		indices[6] = temp[7];
		indices[7] = temp[6];
	    }

	    // get corner coordinates, fixed, and velocities
	    crds.assign(indices.size(),VDouble());
	    fixed.assign(indices.size(),0);
	    vels.assign(indices.size(),VDouble());
	    dvns.assign(indices.size(),VDouble());
	    pns.assign(indices.size(), 0.0);
	    dpns.assign(indices.size(), 0.0);

	    for (unsigned int i=0; i<indices.size(); ++i) {

		// get crds
		getCrds(indices[i], crds[i]);

		// check bnode
		std::map<VInt,BNode>::iterator it = bnodes.find(indices[i]);
		if (it == bnodes.end()) continue;

		// get bnode
		BNode& bnode = it->second;
		
		// get fixed
		for (int j=0; j<bnode.size(); ++j) {
		    if (bnode.type[j] != FLUID) {
			fixed[i] = 1;
			break;
		    }
		}
		if (fixed[i] == 1) {
		    continue;
		}

		// get vn and dvn
		if (bnode.tags.size() != 1) {
		    std::cerr << "WARNING: fluid bnode tags.size() != 1 ";
		    std::cerr << "-- BgMesh::convectParticle\n";
		    return -1;
		}
		vels[i] = bnode.vn[0];
		dvns[i] = bnode.dvn[0];
		pns[i] = bnode.pn[0];
		dpns[i] = bnode.dpn[0];
	    }
	}

	// get particle velocity
	VDouble pvel;
	if (interpolate(pt,indices,vels,dvns,pns,dpns,crds,fixed,pvel) < 0) {
	    std::cerr << "WARNING: failed to interpolate particle velocity";
	    std::cerr << "-- BgMesh::convectParticle\n";
	    return -1;
	}

	// original->dest
	VDouble pdest = pvel;
	pdest *= subdt;
	pdest += pcrds;

	// check which boundary is fixed
	VInt boundfix;
	if (ndm == 2) {
	    boundfix.resize(4);
	    boundfix[0] = fixed[0] && fixed[3];
	    boundfix[1] = fixed[0] && fixed[1];
	    boundfix[2] = fixed[1] && fixed[2];
	    boundfix[3] = fixed[2] && fixed[3];
	    
	} else if (ndm == 3) {
	    boundfix.resize(6);
	    boundfix[0] = fixed[0] && fixed[3] && fixed[4] && fixed[7];
	    boundfix[1] = fixed[0] && fixed[1] && fixed[4] && fixed[5];
	    boundfix[2] = fixed[0] && fixed[1] && fixed[2] && fixed[3];
	    boundfix[3] = fixed[1] && fixed[2] && fixed[5] && fixed[6];
	    boundfix[4] = fixed[2] && fixed[3] && fixed[6] && fixed[7];
	    boundfix[5] = fixed[4] && fixed[5] && fixed[6] && fixed[7];
	}

	// xlow,ylow,zlow,xup,yup,zup
	VDouble bound(boundfix.size());
	if (ndm == 2) {
	    bound[0] = crds[0][0];
	    bound[1] = crds[0][1];
	    bound[2] = crds[2][0];
	    bound[3] = crds[2][1];
	    
	} else if (ndm == 3) {
	    bound[0] = crds[0][0];
	    bound[1] = crds[0][1];
	    bound[2] = crds[0][2];
	    bound[3] = crds[6][0];
	    bound[4] = crds[6][1];
	    bound[5] = crds[6][2];
	}

	// lower boundaries
	for (int i=0; i<ndm; ++i) {
	    if (boundfix[i]==1 && pdest[i]<bound[i]+bsize*meshtol) {
		double dist = bound[i] - pdest[i];
		if (dist < bsize*meshtol) {
		    dist = bsize*meshtol;
		} else if (dist > (1-meshtol)*bsize) {
		    dist = bsize*(1-meshtol);
		}
		pdest[i] = bound[i] + dist;
	    }
	}

	// upper boundaries
	for (int i=ndm; i<2*ndm; ++i) {
	    if (boundfix[i]==1 && pdest[i-ndm]>bound[i]-bsize*meshtol) {
		double dist = pdest[i-ndm] - bound[i];
		if (dist < bsize*meshtol) {
		    dist = bsize*meshtol;
		} else if (dist > (1-meshtol)*bsize) {
		    dist = bsize*(1-meshtol);
		}
		pdest[i-ndm] = bound[i] - dist;
	    }
	}

	// move the particles
	pdest -= pcrds;
	for (int i=0; i<ndm; ++i) {
	    if (pdest[i] > (1-meshtol)*bsize) {
		pdest[i] = (1-meshtol)*bsize;
	    } else if (pdest[i] < -(1-meshtol)*bsize) {
		pdest[i] = -(1-meshtol)*bsize;
	    }
	}
	pt->move(pdest,subdt);
    }

    return 0;
}

int
BackgroundMesh::interpolate(Particle* pt, const VVInt& index,
			    const VVDouble& vels, const VVDouble& dvns,
			    const VDouble& pns, const VDouble& dpns,
			    const VVDouble& crds,
			    const VInt& fixed, VDouble& pvel)
{
    int ndm = OPS_GetNDM();
    
    // check
    if (ndm == 2) {
	if (index.size() != 4) return 0;
	if (vels.size() != 4) return 0;
	if (pns.size() != 4) return 0;
	if (dpns.size() != 4) return 0;
	if (crds.size() != 4) return 0;
	if (fixed.size() != 4) return 0;
    } else if (ndm == 3) {
	if (index.size() != 8) return 0;
	if (vels.size() != 8) return 0;
	if (pns.size() != 8) return 0;
	if (dpns.size() != 8) return 0;
	if (crds.size() != 8) return 0;
	if (fixed.size() != 8) return 0;
    }
    
    const VDouble& pcrds = pt->getCrds();
    if ((int)pcrds.size() != ndm) {
	opserr << "WARNING: pcrds.size() != ndm -- BgMesh::interpolate\n";
	return -1;
    }

    // get shape functions for pt
    VDouble N;
    if (ndm == 2) {
	double hx = (crds[1][0]+crds[2][0])/2.0-crds[0][0];
	double hy = (crds[2][1]+crds[3][1])/2.0-crds[0][1];
	getNForRect(crds[0][0], crds[0][1], hx, hy, pcrds[0], pcrds[1], N);
    } else if (ndm == 3) {
	double hx = (crds[1][0]+crds[2][0])/2.0-crds[0][0];
	double hy = (crds[2][1]+crds[3][1])/2.0-crds[0][1];
	double hz = (crds[4][2]+crds[5][2]+crds[6][2]+crds[7][2])/4.0-crds[0][2];
	getNForRect(crds[0][0],crds[0][1],crds[0][2],
		    hx,hy,hz,pcrds[0],pcrds[1],pcrds[2],N);
    }

    // average velocity
    pvel.resize(ndm, 0.0);
    VDouble pdvn(ndm);
    double ppre=0.0, pdp=0.0;
    double Nsum = 0.0;
    for (unsigned int j=0; j<vels.size(); ++j) {

	// no vel at this corner
	if (vels[j].empty()) continue;

	// free point
	for (int k=0; k<ndm; ++k) {
	    pvel[k] += N[j]*vels[j][k];
	    if (pt->isUpdated() == false) {
		pdvn[k] += N[j]*dvns[j][k];
	    }
	}
	if (pt->isUpdated() == false) {
	    ppre += N[j]*pns[j];
	    pdp += N[j]*dpns[j];
	}
	Nsum += N[j];
    }

    if (Nsum > 0) {
	pvel /= Nsum;
	if (pt->isUpdated() == false) {
	    pdvn /= Nsum;
	    ppre /= Nsum;
	    pdp /= Nsum;
	    pt->setVel(pvel);
	    pt->setAccel(pdvn);
	    pt->setPressure(ppre);
	    pt->setPdot(pdp);
	}
    } else {
	pvel = pt->getVel();
    }

    return 0;
}

void
BackgroundMesh::getFreeSurface(VInt& fsnodes)
{
    for (std::map<VInt,BNode>::iterator it=bnodes.begin(); it!=bnodes.end(); ++it) {

	// check if fixed
	VInt index = it->first;
	BNode& bnode = it->second;
	// if (bnode.type != FLUID) {
	//     continue;
	// }
	if (bnode.tags.size() != 1) {
	    continue;
	}

	// get neighbors corners
	index -= 1;
	VVInt indices;
	getCorners(index,1,indices);
	bool free = false;
	for (unsigned int i=0; i<indices.size(); ++i) {
	    std::map<VInt,BCell>::iterator it = bcells.find(indices[i]);
	    if (it == bcells.end()) {
		free = true;
		break;
	    }
	    if (it->second.pts.empty()) {
		free = true;
		break;
	    }
	}

	// add free surface node
	if (free) {
	    fsnodes.push_back(bnode.tags[0]);
	}
    }
}

int
BackgroundMesh::freeSurface()
{
    if (!freesurf) {
	return 0;
    }

    // get domain
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return -1;

    // clear previous sps
    // for (int i=0; i<(int)sps.size(); i++) {
    // 	SP_Constraint* sp = domain->removeSP_Constraint(sps[i]);
    // 	if (sp != 0) {
    // 	    delete sp;
    // 	}
    // }
    // sps.clear();

    // get free surface
    VInt fsnodes;
    getFreeSurface(fsnodes);
    if (fsnodes.empty()) {
	return 0;
    }

    // add SP_Constraint for pressure nodes
    for (unsigned int i=0; i<fsnodes.size(); ++i) {
	Pressure_Constraint* pc = domain->getPressure_Constraint(fsnodes[i]);
	if (pc == 0) {
	    continue;
	}
	pc->setFreeSurf();

	Node* node = pc->getPressureNode();
	if (node == 0) {
	    continue;
	}

	Vector vel = node->getTrialVel();
	vel.Zero();
	node->setTrialVel(vel);
	node->commitState();

	// SP_Constraint* sp = new SP_Constraint(node->getTag(), 0, atomsP, true);
	// if (sp == 0) {
	//     opserr<<"WARING: run out of memory -- BackgroundMesh::freeSurface\n";
	//     return -1;
	// }
	// if (domain->addSP_Constraint(sp) == false) {
	//     opserr<<"WARNING: failed to add sp to domain -- BackgroundMesh::freeSurface\n";
	//     delete sp;
	//     return -1;
	// }
	// sps.push_back(sp->getTag());
    }

    return 0;
}
