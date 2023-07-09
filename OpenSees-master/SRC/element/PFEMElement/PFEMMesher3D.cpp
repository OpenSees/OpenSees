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
// $Date: 2013-2-25 11:25:05 $
// $Source: SRC/element/PFEMElement/PFEMMesher3D.h,v $
                                                                        
// Written: Minjie Zhu
// Created: Feb. 25
//
// Description: The class PFEMMesher3D encapsulates all necessary
// routines for PFEM analysis for fluid-structure interaction
// The triangulation uses the 
// tetgen program written by Hang Si http://tetgen.berlios.de
//                                                    

//
// What: "@(#) PFEMMesher3D.cpp, revA"

#include "PFEMMesher3D.h"
#include <Matrix.h>
#include <ID.h>
#include <Domain.h>
#include <Node.h>
#include <Element.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <Pressure_Constraint.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <Pressure_ConstraintIter.h>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <PFEMElement3D.h>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <iostream>

PFEMMesher3D::PFEMMesher3D()
    :bound(6), avesize(0.0)
{
}

PFEMMesher3D::~PFEMMesher3D()
{
}



// PLC
int 
PFEMMesher3D::discretize(int startnode, const PointVec& points, const FacetVec facets,
                         const PointVec& holes, double maxvol, int ndf, 
                         const std::vector<int>& fix, const std::vector<double>& mass,
                         Domain* theDomain, int& endnode)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher3D::discretize\n";
        return -1;
    }

    // structs
    tetgenio in, out;

    // input points
    in.numberofpoints = points.size();
    if(in.numberofpoints > 0) {
        in.pointlist = new double[in.numberofpoints*3];
    } else {
        in.pointlist = NULL;
    }
    for(int i=0; i<in.numberofpoints; i++) {
        const Point& point = points[i];
        in.pointlist[3*i] = point.x();
        in.pointlist[3*i+1] = point.y();
        in.pointlist[3*i+2] = point.z();
        //opserr<<"point "<<i<<": "<<point.x()<<" "<<point.y()<<" "<<point.z()<<"\n";
    }

    // input facets
    in.numberoffacets = facets.size();
    if(in.numberoffacets > 0) {
        in.facetlist = new tetgenio::facet[in.numberoffacets];
    } else {
        in.facetlist = NULL;
    }
    for(int i=0; i<in.numberoffacets; i++) {
        const Facet& facet = facets[i];
        const PolygonVec& polygons = facet.polygons();
        const PointVec& holes = facet.holes();

        // new facet
        tetgenio::facet& tfacet = in.facetlist[i];

        // polygons
        tfacet.numberofpolygons = polygons.size();
        if(tfacet.numberofpolygons > 0) {
            tfacet.polygonlist = new tetgenio::polygon[tfacet.numberofpolygons];
        } else {
            tfacet.polygonlist = NULL;
        }
        for(int j=0; j<tfacet.numberofpolygons; j++) {
            tetgenio::polygon& tpolygon = tfacet.polygonlist[j];
            const Polygon& polygon = polygons[j];
            tpolygon.numberofvertices = polygon.size();
            if(tpolygon.numberofvertices > 0) {
                tpolygon.vertexlist = new int[tpolygon.numberofvertices];
            } else {
                tpolygon.vertexlist = NULL;
            }
            for(int k=0; k<tpolygon.numberofvertices; k++) {
                tpolygon.vertexlist[k] = polygon[k];
            }
        }

        // holes
        tfacet.numberofholes = holes.size();
        if(tfacet.numberofholes > 0) {
            tfacet.holelist = new double[tfacet.numberofholes*3];
        } else {
            tfacet.holelist = NULL;
        }
        for(int j=0; j<tfacet.numberofholes; j++) {
            const Point& hole = holes[j];
            tfacet.holelist[3*j] = hole.x();
            tfacet.holelist[3*j+1] = hole.y();
            tfacet.holelist[3*j+2] = hole.z();
        }
    }

    // input holes
    in.numberofholes = holes.size();
    if(in.numberofholes > 0) {
        in.holelist = new double[in.numberofholes*3];
    } else {
        in.holelist = NULL;
    }
    for(int i=0; i<in.numberofholes; i++) {
        const Point& hole = holes[i];
        in.holelist[3*i] = hole.x();
        in.holelist[3*i+1] = hole.y();
        in.holelist[3*i+2] = hole.z();
    }

    // conforming Delaunay triangulation
    std::stringstream ss;
    ss << "Qzqpa"<< maxvol;
    char* s = (char*)ss.str().c_str();
    opserr<<"Start Delaunay Triangulation -- If got segmentation fault, please check inputs \n";
    tetrahedralize(s, &in, &out);
    opserr<<"Finish Delaunay Triangulation\n";

    // read outputs
    endnode = startnode-1;
    for(int i=0; i<out.numberofpoints; i++) {
        const double& x = out.pointlist[3*i];
        const double& y = out.pointlist[3*i+1];
        const double& z = out.pointlist[3*i+2];

        // create nodes
        Node* theNode = new Node(++endnode, ndf, x, y, z);
        if (theNode == 0) {
            opserr << "WARNING ran out of memory creating node\n";
            opserr << "node: " << endnode << "\n";
            opserr << "PFEMMesher3D::discretize\n";
            return -1;
        }
        if(theDomain->addNode(theNode) == false) {
            opserr << "WARNING failed to add node to the domain\n";
            opserr << "node: " << endnode << "\n";
            opserr << "PFEMMesher3D::discretize\n";
            delete theNode; // otherwise memory leak
            return -1;
        }
        
        // fix
        int numfix = fix.size();
        if(numfix > ndf) numfix = ndf;
        for(int j=0; j<numfix; j++) {
            if(fix[j] != 0) {
                SP_Constraint* theSP = new SP_Constraint(endnode, j, 0.0, true);
                if(theSP == 0) {
                    opserr << "WARNING ran out of memory creating SP_Constraint\n";
                    opserr << "node: " << endnode << "\n";
                    opserr << "PFEMMesher3D::discretize\n";
                    return -1;
                }
                if(theDomain->addSP_Constraint(theSP) == false) {
                    opserr << "WARNING failed to add SP_Constraint to the domain\n";
                    opserr << "node: " << endnode << "\n";
                    opserr << "PFEMMesher3D::discretize\n";
                    delete theSP; // otherwise memory leak
                    return -1;
                }
            }
        }

        // mass
        int nummass = mass.size();
        if(nummass > 0) {
            if(nummass > ndf) nummass = ndf;
            Matrix theMass(ndf, ndf);
            for(int j=0; j<nummass; j++) {
                theMass(j,j) = mass[j];
            }
            if(theNode->setMass(theMass) < 0) {
                opserr<<"WARNING: failed to set mass of node "<<endnode;
                opserr<<" -- PFEMMesher3D::discretize\n";
                return -1;
            }
        }
    }

    return 0;
}

// cube
int 
PFEMMesher3D::discretize(int startnode, const Point& pt, const dvector& hs, 
                         const ivector& ns, int ndf, const ivector& fix,
                         const dvector& mass, const dvector& vel,
                         Domain* theDomain, int& endnode)
{
    if(hs.size() < 3) return 0;
    if(ns.size() < 3) return 0;

    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher3D::discretize\n";
        return -1;
    }

    for(int i=0; i<(int)ns.size(); i++) {
        if(ns[i] < 0) {
            opserr<<"WARNING: negative number of divisions";
            opserr<<"-- PFEMMesher3D::discretize\n";
            return -1;
        }
    }

    double x0 = pt.x();
    double y0 = pt.y();
    double z0 = pt.z();

    // velocity
    int numvel = vel.size();
    if(numvel > ndf) numvel = ndf;
    Vector tvel(ndf);
    for(int i=0; i<numvel; i++) {
        tvel(i) = vel[i];
    }

    // fix
    int numfix = fix.size();
    if(numfix > ndf) numfix = ndf;

    // mass
    int nummass = mass.size();
    if(nummass > ndf) nummass = ndf;
    Matrix theMass(ndf, ndf);
    for(int i=0; i<nummass; i++) {
        theMass(i,i) = mass[i];
    }

    int newtag = startnode-1;
    for(int i=0; i<ns[0]; i++) {
        double x = x0+hs[0]*i;
        for(int j=0; j<ns[1]; j++) {
            double y = y0+hs[1]*j;
            for(int k=0; k<ns[2]; k++) {
                double z = z0+hs[2]*k;

                // define nodes
                Node* theNode = new Node(++newtag, ndf, x, y, z);
                if (theNode == 0) {
                    opserr << "WARNING ran out of memory creating node\n";
                    opserr << "node: " << newtag << "\n";
                    opserr << "PFEMMesher3D::discretize\n";
                    return -1;
                }
                if(theDomain->addNode(theNode) == false) {
                    opserr << "WARNING failed to add node to the domain\n";
                    opserr << "node: " << newtag << "\n";
                    opserr << "PFEMMesher3D::discretize\n";
                    delete theNode; // otherwise memory leak
                    return -1;
                }

                // initial velocity
                if(numvel > 0) {
                    theNode->setTrialVel(tvel);
                    theNode->commitState();
                }

                // fix
                for(int ifix=0; ifix<numfix; ifix++) {
                    if(fix[ifix] != 0) {
                        SP_Constraint* theSP = new SP_Constraint(newtag, ifix, 0.0, true);
                        if(theSP == 0) {
                            opserr << "WARNING ran out of memory creating SP_Constraint\n";
                            opserr << "node: " << newtag << "\n";
                            opserr << "PFEMMesher3D::discretize\n";
                            return -1;
                        }
                        if(theDomain->addSP_Constraint(theSP) == false) {
                            opserr << "WARNING failed to add SP_Constraint to the domain\n";
                            opserr << "node: " << newtag << "\n";
                            opserr << "PFEMMesher3D::discretize\n";
                            delete theSP; // otherwise memory leak
                            return -1;
                        }
                    }
                }

                // mass
                if(nummass > 0) {
                    if(theNode->setMass(theMass)  < 0) {
                        opserr<<"WARNING: failed to set mass of node "<<newtag;
                        opserr<<" -- PFEMMesher3D::discretize\n";
                        return -1;
                    }
                    
                }
            }
        }
    }

    endnode = newtag;
    return 0;
}


int 
PFEMMesher3D::doTriangulation(const std::vector<int>& nodes, double alpha, double volthresh,
                              const std::vector<int>& addnodes, Domain* theDomain, 
                              std::vector<int>& eles)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher3D::doTriangulation\n";
        return -1;
    }

    const int ndm = 3;
    const int nptet = 4;

    // structs
    tetgenio in, out;

    // number of nodes
    int numnodes=0, numadd=0;
    for(int i=0; i<(int)nodes.size()/2; i++) {
        numnodes += nodes[2*i+1] - nodes[2*i] + 1;
    }
    for(int i=0; i<(int)addnodes.size()/2; i++) {
        numadd += addnodes[2*i+1] - addnodes[2*i] + 1;
    }
    in.numberofpoints = numnodes+numadd;
    if(in.numberofpoints < nptet) {
        return 0;
    }

    // get space for pointlist
    in.pointlist = new double[in.numberofpoints*ndm];
    if(in.pointlist == 0) {
        opserr<<"WARNING: no enough memory -- PFEMMesher3D::doTriangulation\n";
        return -1;
    }

    // copy all nodes
    int loc = 0, numpoints = 0;
    int nodesloc = 0, addnodesloc = 0;
    ID p2nd(0, in.numberofpoints);
    for(int i=0; i<in.numberofpoints; i++) {

        // get node
        int ndtag = 0;
        if(i<numnodes) {
            ndtag = nodes[2*nodesloc]+loc;
            //opserr<<"ndtag = "<<ndtag<<"\n";
            if(nodes[2*nodesloc+1]-nodes[2*nodesloc]==loc) {
                nodesloc++;
                loc = 0;
            } else {
                loc++;
            }
        } else {
            ndtag = addnodes[2*addnodesloc]+loc;
            if(addnodes[2*addnodesloc+1]-addnodes[2*addnodesloc]==loc) {
                addnodesloc++;
                loc = 0;
            } else {
                loc++;
            }
        }
 
        Node* node = theDomain->getNode(ndtag);
        if(node == 0) continue;

        // get coordinates
        const Vector& coord = node->getCrds();
        if(coord.Size() < ndm) {
            opserr<<"WARNING: 3d is required -- PFEMMesher3D::doTriangulation\n";
            return -1;
        }
        const Vector& disp = node->getTrialDisp();
        if(disp.Size() < ndm) {
            opserr<<"WARNING: 3 ndf is required -- PFEMMesher3D::doTriangulation\n";
            return -1;
        }

        // set pointlist
        for(int j=0; j<ndm; j++) {
            in.pointlist[ndm*numpoints+j] = coord(j) + disp(j);
        }
        p2nd[numpoints++] = ndtag;

    }
    in.numberofpoints = numpoints;

    // Delaunay Triangulation
    char s[] = "Qz";
    opserr<<"Start Delaunay Triangulation -- If got segmentation fault, please check inputs \n";
    tetrahedralize(s, &in, &out);
    opserr<<"Finish Delaunay Triangulation\n";

    // no outputs
    if(out.numberoftetrahedra < 1) {
        return 0;
    }

    // do alpha shape test
    eles.clear();
    if(alpha > 0) {

        // need mesh functions
        tetgenmesh mesh;

        // radius and average size of tetrahedra
        avesize = 0.0;
        Vector radius(out.numberoftetrahedra);
        Vector volume(out.numberoftetrahedra);
        for(int i=0; i<out.numberoftetrahedra; i++) {

            // tetrahedra points
            ID pt(nptet);
            for(int j=0; j<nptet; j++) {
                pt(j) = out.tetrahedronlist[out.numberofcorners*i+j];
            }

            // nodal coordinates
            double* ppointer[4];
            Matrix pcoord(nptet,ndm);
            for(int j=0; j<nptet; j++) {
                for(int k=0; k<ndm; k++) {
                    pcoord(j,k) = out.pointlist[ndm*pt(j)+k];
                }
                ppointer[j] = &(out.pointlist[ndm*pt(j)]);
            }

            // size of tetrahedra
            double he = -1.0;
            for(int j=0; j<nptet; j++) {
                for(int k=j+1; k<nptet; k++) {
                    double h = 0.0;
                    for(int l=0; l<ndm; l++) {
                        h += (pcoord(j,l)-pcoord(k,l))*(pcoord(j,l)-pcoord(k,l));
                    }
                    if(h<he || he==-1.0) {
                        he = h;
                    }
                }
            }
            avesize += sqrt(he);

            // volume
            double A[4][4], D;
            int indx[4];
            for(int j=0; j<nptet-1; j++) {
                for(int k=0; k<ndm; k++) {
                    A[j][k] = pcoord(j,k) - pcoord(3,k);
                }
            }
            mesh.lu_decmp(A,3,indx,&D,0);
            volume(i) = fabs(A[indx[0]][0]*A[indx[1]][1]*A[indx[2]][2]) / 6.0;

            // radius
            double* pradius = &radius(i);
            mesh.circumsphere(ppointer[0],ppointer[1],ppointer[2],ppointer[3],NULL,pradius);
        }
        avesize /= out.numberoftetrahedra;

        // alpha test
        double totalvolume = 0.0;
        for(int i=0; i<out.numberoftetrahedra; i++) {
            // coplanar factor
            // double coplane = volume(i) / pow(avesize,3);
            // if(q < 1e-6) {
            //     // opserr<<"coplanar: volume = "<<volume(i)<<", q = "<<q<<"\n";
            //     continue;
            // }
            double vsphere = 4./3.*3.14*radius(i)*radius(i)*radius(i);
            if(radius(i) / avesize <= alpha && volume(i)/vsphere > volthresh) {
                // pass the test
                totalvolume += volume(i);

                // check if all nodes are additional nodes
                int add[4] = {-1,-1,-1,-1};
                for(int j=0; j<4; j++) {
                    int tag = p2nd(out.tetrahedronlist[out.numberofcorners*i+j]);
                    for(int k=0; k<(int)addnodes.size()/2; k++) {
                        if(tag>=addnodes[2*k] && tag<=addnodes[2*k+1]) {
                            add[j] = k;
                            break;
                        }
                    }
                }

                // add ele
                // if(add[0]==add[1] && add[1]==add[2] && add[2]==add[3] && add[3]!=-1) continue;
                if(add[0]!=-1 && add[1]!=-1 && add[2]!=-1 && add[3]!=-1) continue;
                for(int j=0; j<4; j++) {
                    int tag = p2nd(out.tetrahedronlist[out.numberofcorners*i+j]);
                    //opserr<<tag<<" ";
                    eles.push_back(tag);
                }
            }
        }

    } else if(alpha < 0) {

        // copy triangles
        for(int i=0; i<out.numberoftetrahedra; i++) {
            for(int j=0; j<4; j++) {
                int tag = p2nd(out.tetrahedronlist[out.numberofcorners*i+j]);
                eles.push_back(tag);
            }
        }
    }

    return 0;
}

int 
PFEMMesher3D::doTriangulation(int startele, const std::vector<int>& nodes, double alpha, 
                              double volthresh, const std::vector<int>& addnodes, Domain* theDomain,
                              double rho, double mu, double b1, double b2, double b3)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }
    //Timer timer;
    // do triangulation
    std::vector<int> eles;
    int res = doTriangulation(nodes, alpha,
                              volthresh, addnodes, 
                              theDomain, eles);
    if(res < 0) {
        opserr<<"WARNING: failed to do triangulation --";
        opserr<<"PFEMMesher3D::soTriangulation\n";
        return res;
    }

    // add PFEM elements
    //timer.start();
    int numeles = eles.size()/4;
    if(numeles == 0) return 0;
    int etag = startele-1;
    for(int i=0; i<numeles; i++) {
        //opserr<<eles(4*i)<<" "<<eles(4*i+1)<<" "<<eles(4*i+2)<<" "<<eles(4*i+3)<<"\n";
        PFEMElement3D* theEle = new PFEMElement3D(++etag, eles[4*i], eles[4*i+1], eles[4*i+2],
                                                  eles[4*i+3], rho, mu, b1, b2, b3);

        if(theEle == 0) {
            opserr<<"WARNING: no enough memory -- ";
            opserr<<" -- PFEMMesher3D::doTriangulation\n";
            return -1;
        }
        if(theDomain->addElement(theEle) == false) {
            opserr<<"WARNING: failed to add element to domain -- ";
            opserr<<" -- PFEMMesher3D::doTriangulation\n";
            delete theEle;
            return -1;
        }
    }
    //timer.pause();
    //opserr<<"create PFEM elements :"<<timer;
    return etag;
}

int 
PFEMMesher3D::save(const char* filename, const ID& snodes, int step, Domain* theDomain)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher3D::save\n";
        return -1;
    }

    // Pressure_Constraints
    std::map<int, std::pair<int,Node*> > nodes;
    std::map<int, Node*> pnodes;
    Pressure_ConstraintIter& thePCs = theDomain->getPCs();
    Pressure_Constraint* thePC = 0;
    int numisolated = 0;
    while((thePC = thePCs()) != 0) {
        int ntag = thePC->getTag();
        // int ptag = thePC->getPressureNode();
        Node* nnd = theDomain->getNode(ntag);
        // Node* pnd = theDomain->getNode(ptag);
        Node* pnd = thePC->getPressureNode();
        int ptag = pnd->getTag();
        if(nnd!=0 && pnd!=0) {
            std::pair<int,Node*>& nnode = nodes[ntag];
            nnode.second = nnd;
            pnodes[ptag] = pnd;
            if(thePC->isFluid()) {
                nnode.first = 0;
            } else if(thePC->isStructure()) {
                nnode.first = 3;
            } else if(thePC->isInterface()) {
                nnode.first = 2;
            } else if(thePC->isIsolated()) {
                nnode.first = 1;
                numisolated++;
            }
        }
    }

    // structure nodes
    NodeIter& theNodes = theDomain->getNodes();
    Node* theNode = 0;
    while((theNode = theNodes()) != 0) {
        int tag = theNode->getTag();
        bool haveit = false;
        int stype = 0;
        for(int i=0; i<snodes.Size()/2; i++) {
            if(tag>=snodes(2*i) && tag<=snodes(2*i+1)) {
                haveit = true;
                stype = i;
                break;
            }
        }
        if(haveit) {
            std::pair<int,Node*>& node = nodes[tag];
            node.second = theNode;
            node.first = stype+4;
        }
    }

    // elements
    std::map<int,Element*> elements;
    ElementIter& theEles = theDomain->getElements();
    Element* theEle = 0;
    while((theEle = theEles()) != 0) {
        const ID& ntags = theEle->getExternalNodes();
        bool haveit = true;
        for(int i=0; i<ntags.Size(); i++) {
            if(nodes.find(ntags(i))==nodes.end() && pnodes.find(ntags(i))==pnodes.end()) {
                haveit = false;
                break;
            }
        }
        if(haveit) {
            elements[theEle->getTag()] = theEle;
        }
    }

    // open file
    std::stringstream ss;
    ss << filename << "-" << step << ".msh";
    std::ofstream file(ss.str().c_str());

    // write msh file header
    file << "$MeshFormat\n";
    file << "2.2 0 8\n";
    file << "$EndMeshFormat\n";

    // write nodes
    int startnode = 0;
    file << "$Nodes\n";
    file << nodes.size() << "\n";
    for(std::map<int,std::pair<int,Node*> >::iterator it=nodes.begin(); it!=nodes.end(); it++) {
        int ntag = it->first;
        if(ntag<0 && startnode==0) startnode = -ntag;
        Node* theNode = it->second.second;
        const Vector& coord = theNode->getCrds();
        const Vector& disp = theNode->getTrialDisp();
        if(disp.Size()>2 && coord.Size()>2) {
            file << ntag+startnode << " ";
            file << coord(0)+disp(0) << " ";
            file << coord(1)+disp(1) << " ";
            file << coord(2)+disp(2) << "\n";
        } else if(disp.Size()>1 && coord.Size()>1) {
            file << ntag+startnode << " ";
            file << coord(0)+disp(0) << " ";
            file << coord(1)+disp(1) << " ";
            file << 0.0 << "\n";
        } else if(disp.Size()>0 && coord.Size()>0) {
            file << ntag+startnode << " ";
            file << coord(0)+disp(0) << " ";
            file << 0.0 << " ";
            file << 0.0 << "\n";
        } else {
            opserr<<"WARNING: size of disp or coord equal to zero ";
            opserr<<"-- PFEMMesher3D::save\n";
        }
    }
    file << "$EndNodes\n";

    // write elements
    file << "$Elements\n";
    file << elements.size()+numisolated << "\n";
    int startele = 0;
    int endele = 0;
    for(std::map<int,Element*>::iterator it=elements.begin(); it!=elements.end(); it++) {
        int etag = it->first;
        if(etag<0 && startele==0) startele = -etag;
        Element* theEle = it->second;
        const ID& ntags = theEle->getExternalNodes();
        const char* type = theEle->getClassType();
        int numnodes = ntags.Size();

        file << etag+startele << " ";
        if(strcmp(type, "PFEMElement3D") == 0) {
            file << " 4 0 ";
            for(int i=0; i<numnodes/2; i++) {
                file << ntags(2*i) << " ";
            }
            file << "\n";
        } else if(numnodes == 1) {
            file << " 15 0 ";
            file << ntags(0) << "\n";
        } else if(numnodes == 2) {
            file << " 1 0 ";
            for(int i=0; i<numnodes; i++) {
                file << ntags(i) << " ";
            }
            file << "\n";
        } else if(numnodes == 3) {
            file << " 2 0 ";
            for(int i=0; i<numnodes; i++) {
                file << ntags(i) << " ";
            }
            file << "\n";
        } else if(numnodes == 4) {
            file << " 3 0 ";
            for(int i=0; i<numnodes; i++) {
                file << ntags(i) << " ";
            }
            file << "\n";
        }
        endele = etag+startele;
    }

    for(std::map<int,std::pair<int,Node*> >::iterator it=nodes.begin(); it!=nodes.end(); it++) {
        int ntag = it->first;
        int type = it->second.first;
        if(type == 1) {
            file << ++endele << " 15 0 "<<ntag+startnode<<"\n";
        }
    }
    file << "$EndElements\n";

    // write model view
    file << "$NodeData\n";
    file << 1 << "\n";
    file << "\"model-" << step <<"\"\n";
    file << 1 << "\n";
    file << theDomain->getCurrentTime() << "\n";
    file << 3 << "\n";
    file << step << "\n";
    file << 1 << "\n";
    file << nodes.size() << "\n";
    for(std::map<int,std::pair<int,Node*> >::iterator it=nodes.begin(); it!=nodes.end(); it++) {
        int ntag = it->first;
        int type = it->second.first;
        file << ntag+startnode << " " << type << "\n";
    }
    file << "$EndNodeData\n";

    // write pressure view
    file << "$NodeData\n";
    file << 1 << "\n";
    file << "\"pressure-" << step << "\"\n";
    file << 1 << "\n";
    file << theDomain->getCurrentTime() << "\n";
    file << 3 << "\n";
    file << step << "\n";
    file << 1 << "\n";
    file << nodes.size() << "\n";
    for(std::map<int,std::pair<int,Node*> >::iterator it=nodes.begin(); it!=nodes.end(); it++) {
        int ntag = it->first;
        Pressure_Constraint* thePC = theDomain->getPressure_Constraint(ntag);
        if(thePC == 0) {
            file << ntag << " " << 0 << "\n";
        } else {
            // int ptag = thePC->getPressureNode();
            // Node* pnode = pnodes[ptag];
            Node* pnode = thePC->getPressureNode();
            if(pnode == 0) {
                file << ntag << " " << 0 << "\n";
            } else {
                const Vector& vel = pnode->getTrialVel();
                file << ntag << " " << vel(0) << "\n";
            }
            
        }
    }
    file << "$EndNodeData\n";

    // write displacement view
    file << "$NodeData\n";
    file << 1 << "\n";
    file << "\"displacement-" << step <<"\"\n";
    file << 1 << "\n";
    file << theDomain->getCurrentTime() << "\n";
    file << 3 << "\n";
    file << step << "\n";
    file << 3 << "\n";
    file << nodes.size() << "\n";
    for(std::map<int,std::pair<int,Node*> >::iterator it=nodes.begin(); it!=nodes.end(); it++) {
        int ntag = it->first;
        Node* theNode = it->second.second;
        const Vector& disp = theNode->getTrialDisp();
        if(disp.Size() > 2) {
            file << ntag << " " << disp(0) << " ";
            file << disp(1) << " " << disp(2) << "\n";
        } else if(disp.Size() > 1) {
            file << ntag << " " << disp(0) << " ";
            file << disp(1) << " " << 0.0 << "\n";
        } else if(disp.Size() > 0) {
            file << ntag << " " << disp(0) << " ";
            file << 0.0 << " " << 0.0 << "\n";
        } else {
            opserr<<"WARNING: size of disp equal to zero ";
            opserr<<"-- PFEMMesher3D::save\n";
        }
    }
    file << "$EndNodeData\n";

    // write velocity view
    file << "$NodeData\n";
    file << 1 << "\n";
    file << "\"velocity-" << step <<"\"\n";
    file << 1 << "\n";
    file << theDomain->getCurrentTime() << "\n";
    file << 3 << "\n";
    file << step << "\n";
    file << 3 << "\n";
    file << nodes.size() << "\n";
    for(std::map<int,std::pair<int,Node*> >::iterator it=nodes.begin(); it!=nodes.end(); it++) {
        int ntag = it->first;
        Node* theNode = it->second.second;
        const Vector& vel = theNode->getTrialVel();
        if(vel.Size() > 2) {
            file << ntag << " " << vel(0) << " ";
            file << vel(1) << " " << vel(2) << "\n";
        } else if(vel.Size() > 1) {
            file << ntag << " " << vel(0) << " ";
            file << vel(1) << " " << 0.0 << "\n";
        } else if(vel.Size() > 0) {
            file << ntag << " " << vel(0) << " ";
            file << 0.0 << " " << 0.0 << "\n";
        } else {
            opserr<<"WARNING: size of vel equal to zero ";
            opserr<<"-- PFEMMesher3D::save\n";
        }
    }
    file << "$EndNodeData\n";

    // write acceleration view
    file << "$NodeData\n";
    file << 1 << "\n";
    file << "\"acceleration-" << step << "\"\n";
    file << 1 << "\n";
    file << theDomain->getCurrentTime() << "\n";
    file << 3 << "\n";
    file << step << "\n";
    file << 3 << "\n";
    file << nodes.size() << "\n";
    for(std::map<int,std::pair<int,Node*> >::iterator it=nodes.begin(); it!=nodes.end(); it++) {
        int ntag = it->first;
        Node* theNode = it->second.second;
        const Vector& accel = theNode->getTrialAccel();
        if(accel.Size() > 2) {
            file << ntag << " " << accel(0) << " ";
            file << accel(1) << " " << accel(2) << "\n";
        } else if(accel.Size() > 1) {
            file << ntag << " " << accel(0) << " ";
            file << accel(1) << " " << 0.0 << "\n";
        } else if(accel.Size() > 0) {
            file << ntag << " " << accel(0) << " ";
            file << 0.0 << " " << 0.0 << "\n";
        } else {
            opserr<<"WARNING: size of accel equal to zero ";
            opserr<<"-- PFEMMesher3D::save\n";
        }
    }
    file << "$EndNodeData\n";
    file.close();

    return 0;
}

void
PFEMMesher3D::setBoundary(double x1, double y1, double z1, double x2, double y2, double z2)
{
    bound(0) = x1;
    bound(1) = y1;
    bound(2) = z2;
    bound(3) = x2;
    bound(4) = y2;
    bound(5) = z2;
}

void 
PFEMMesher3D::removeOutBoundNodes(const ID& nodes, Domain* theDomain)
{
    // check nodes
    std::vector<int> removelist;
    for(int i=0; i<nodes.Size()/2; i++) {
        for(int tag=nodes(2*i); tag<=nodes(2*i+1); tag++) {
            Node* theNode = theDomain->getNode(tag);
            if(theNode == 0) continue;
            const Vector& coord = theNode->getCrds();
            if(coord.Size() < 3) {
                continue;
            }
            const Vector& disp = theNode->getTrialDisp();
            if(disp.Size() < 3) {
                continue;
            }
            double x = coord(0);        // initial
            double y = coord(1);
            double z = coord(2);
            x += disp(0);               // current
            y += disp(1);
            z += disp(2);
        
            // if out of boundary
            if(x<bound(0) || x>bound(3) || y<bound(1) || y>bound(4) || z<bound(2) || z>bound(5)) {
                removelist.push_back(tag);
            }
        }
    }

    // remove nodes
    for(int i=0; i<(int)removelist.size(); i++) {
        Node* theNode = theDomain->removeNode(removelist[i]);
        if(theNode != 0) {
            delete theNode;
        }
        Pressure_Constraint* thePC = theDomain->removePressure_Constraint(removelist[i]);
        if(thePC != 0) {
            delete thePC;
        }
    }
}

