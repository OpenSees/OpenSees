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
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <PFEMElement3D.h>
#include <map>
#include <set>
#include <fstream>

PFEMMesher3D::PFEMMesher3D()
    :PI(3.1415926535897932384626433), bound(6)
{
}

PFEMMesher3D::~PFEMMesher3D()
{
}

// int 
// PFEMMesher3D::discretize (const Vector& points, const ID& segments, 
//                           double maxarea, int ndf, const ID& fix, const Vector& mass, 
//                           Domain* theDomain, ID& nodes)
// {
//     if(theDomain == 0) {
//         opserr<<"WARNING: null domain";
//         opserr<<" -- PFEMMesher3D::discretize\n";
//         return -1;
//     }

//     // structs
//     triangulateio in, out, vout;
//     initializeTri(in);
//     initializeTri(out);
//     initializeTri(vout);

//     // number of segments
//     in.numberofsegments = segments.Size()/2;

//     // no segments
//     if(in.numberofsegments < 1) {
//         return 0;
//     }

//     // numbering points
//     in.numberofpoints = points.Size()/2;

//     // no points
//     if(in.numberofpoints < 2) {
//         return 0;
//     }
    
//     // get space for segmentlist 
//     in.segmentlist = (int*) calloc(segments.Size(), sizeof(int));
//     if(in.segmentlist == NULL) {
//         opserr<<"WARNING: no enough memory -- PFEMMesher3D::discretize\n";
//         return -1;
//     }

//     // copy segments
//     int maxpt=0;
//     for(int i=0; i<segments.Size(); i++) {
//         if(segments(i) > maxpt) maxpt = segments(i);
//         in.segmentlist[i] = segments(i);
//     }
//     if(maxpt >= in.numberofpoints) {
//         opserr<<"WARNING: no point "<<maxpt<<" in the point list -- ";
//         opserr<<"PFEMMesher3D::discretize\n";
//         freeTri(in);
//         return -1;
//     }

//     // get space for pointlist
//     in.pointlist = (double*) calloc(points.Size(), sizeof(double));
//     if(in.pointlist == NULL) {
//         opserr<<"WARNING: no enough memory -- PFEMMesher3D::discretize\n";
//         return -1;
//     }

//     // copy points
//     for(int i=0; i<points.Size(); i++) {
//         in.pointlist[i] = points(i);
//     }

//     // conforming Delaunay triangulation
//     char s[100];
//     sprintf(s,"Qzqpa%.60f",maxarea);
//     triangulate(s, &in, &out, &vout);
//     freeTri(in);

//     // read outputs
//     nodes = ID();
//     if(out.numberofpoints > 0) {
//         nodes.resize(out.numberofpoints);
//     }
//     for(int i=0; i<out.numberofpoints; i++) {
//         const double& x = out.pointlist[2*i];
//         const double& y = out.pointlist[2*i+1];

//         // create nodes
//         int newtag = findNodeTag(theDomain);
//         Node* theNode = new Node(newtag, ndf, x, y);
//         if (theNode == 0) {
//             opserr << "WARNING ran out of memory creating node\n";
//             opserr << "node: " << newtag << "\n";
//             opserr << "PFEMMesher3D::discretize\n";
//             return -1;
//         }
//         if(theDomain->addNode(theNode) == false) {
//             opserr << "WARNING failed to add node to the domain\n";
//             opserr << "node: " << newtag << "\n";
//             opserr << "PFEMMesher3D::discretize\n";
//             delete theNode; // otherwise memory leak
//             return -1;
//         }
//         nodes(i) = newtag;
        
//         // fix
//         int numfix = fix.Size();
//         if(numfix > ndf) numfix = ndf;
//         for(int i=0; i<numfix; i++) {
//             if(fix(i) != 0) {
//                 SP_Constraint* theSP = new SP_Constraint(newtag, i, 0.0, true);
//                 if(theSP == 0) {
//                     opserr << "WARNING ran out of memory creating SP_Constraint\n";
//                     opserr << "node: " << newtag << "\n";
//                     opserr << "PFEMMesher3D::discretize\n";
//                     return -1;
//                 }
//                 if(theDomain->addSP_Constraint(theSP) == false) {
//                     opserr << "WARNING failed to add SP_Constraint to the domain\n";
//                     opserr << "node: " << newtag << "\n";
//                     opserr << "PFEMMesher3D::discretize\n";
//                     delete theSP; // otherwise memory leak
//                     return -1;
//                 }
//             }
//         }

//         // mass
//         int nummass = mass.Size();
//         if(nummass > 0) {
//             if(nummass > ndf) nummass = ndf;
//             Matrix theMass(ndf, ndf);
//             for(int i=0; i<nummass; i++) {
//                 theMass(i,i) = mass(i);
//             }
//             if(theNode->setMass(theMass) < 0) {
//                 opserr<<"WARNING: failed to set mass of node "<<newtag;
//                 opserr<<" -- PFEMMesher3D::discretize\n";
//                 return -1;
//             }
//         }
//     }

//     // free out and vout
//     freeTri(out);
//     freeTri(vout);

//     return 0;
// }

// int 
// PFEMMesher3D::discretize(double x1, double y1, double hx, double hy, double angle,
//                    int nx, int ny, int ndf, const ID& fix, const Vector& mass, 
//                    Domain* theDomain, ID& nodes)
// {
//     if(theDomain == 0) {
//         opserr<<"WARNING: null domain";
//         opserr<<" -- PFEMMesher3D::discretize\n";
//         return -1;
//     }

//     if(nx<0 || ny<0) return 0;
//     double sy = sin((angle+90)*PI/180.0);
//     double cy = cos((angle+90)*PI/180.0);
    
//     // foreach line
//     double lhx = cy*hy;
//     double lhy = sy*hy;

//     for(int i=0; i<=ny; i++) {
//         double lx = i*lhx+x1;
//         double ly = i*lhy+y1;
//         discretize(lx, ly, hx, angle, nx, ndf, fix, mass, theDomain, nodes);
//     }


//     return 0;
// }

// int 
// PFEMMesher3D::discretize(int startnode, Vector point, double h, Vector dir, 
//                          int num, int ndf, const ID& fix, const Vector& mass, 
//                          Domain* theDomain)
// {
//     if(theDomain == 0) {
//         opserr<<"WARNING: null domain";
//         opserr<<" -- PFEMMesher3D::discretize\n";
//         return -1;
//     }

//     if(num < 0) return 0;

//     double s = sin(angle*PI/180.0);
//     double c = cos(angle*PI/180.0);

//     double hx = c*h;
//     double hy = s*h;

//     for(int i=0; i<=num; i++) {
//         double x = i*hx+x1;
//         double y = i*hy+y1;

//         // define node
//         int newtag = findNodeTag(theDomain);
//         Node* theNode = new Node(newtag, ndf, x, y);
//         if (theNode == 0) {
//             opserr << "WARNING ran out of memory creating node\n";
//             opserr << "node: " << newtag << "\n";
//             opserr << "PFEMMesher3D::discretize\n";
//             return -1;
//         }
//         if(theDomain->addNode(theNode) == false) {
//             opserr << "WARNING failed to add node to the domain\n";
//             opserr << "node: " << newtag << "\n";
//             opserr << "PFEMMesher3D::discretize\n";
//             delete theNode; // otherwise memory leak
//             return -1;
//         }
//         nodes.insert(newtag);

//         // fix
//         int numfix = fix.Size();
//         if(numfix > ndf) numfix = ndf;
//         for(int i=0; i<numfix; i++) {
//             if(fix(i) != 0) {
//                 SP_Constraint* theSP = new SP_Constraint(newtag, i, 0.0, true);
//                 if(theSP == 0) {
//                     opserr << "WARNING ran out of memory creating SP_Constraint\n";
//                     opserr << "node: " << newtag << "\n";
//                     opserr << "PFEMMesher3D::discretize\n";
//                     return -1;
//                 }
//                 if(theDomain->addSP_Constraint(theSP) == false) {
//                     opserr << "WARNING failed to add SP_Constraint to the domain\n";
//                     opserr << "node: " << newtag << "\n";
//                     opserr << "PFEMMesher3D::discretize\n";
//                     delete theSP; // otherwise memory leak
//                     return -1;
//                 }
//             }
//         }

//         // mass
//         int nummass = mass.Size();
//         if(nummass > 0) {
//             if(nummass > ndf) nummass = ndf;
//             Matrix theMass(ndf, ndf);
//             for(int i=0; i<nummass; i++) {
//                 theMass(i,i) = mass(i);
//             }
//             if(theNode->setMass(theMass) < 0) {
//                 opserr<<"WARNING: failed to set mass of node "<<newtag;
//                 opserr<<" -- PFEMMesher3D::discretize\n";
//                 return -1;
//             }
//         }
//     }

//     return 0;
// }

// int
// PFEMMesher3D::doTriangulation(const ID& nodes, double alpha, 
//                         const ID& addNodes, Domain* theDomain, 
//                         ID& eles, double rho, double mu,
//                         double b1, double b2)
// {
//     if(theDomain == 0) {
//         opserr<<"WARNING: null domain";
//         opserr<<" -- PFEMMesher3D::doTriangulation\n";
//         return -1;
//     }

//     // do triangulation
//     int res = doTriangulation(nodes, alpha, addNodes, theDomain, eles);
//     if(res < 0) {
//         return res;
//     }

//     // add PFEM elements
//     int numeles = eles.Size()/3;
//     ID eletags(numeles);
//     for(int i=0; i<numeles; i++) {
//         int tag = PFEMMesher3D::findEleTag(theDomain);
//         PFEMElement2D* theEle = new PFEMElement2D(tag, eles(3*i), eles(3*i+1), eles(3*i+2),
//                                                   rho, mu, b1, b2);

//         if(theEle == 0) {
//             opserr<<"WARNING: no enough memory -- ";
//             opserr<<" -- delaunay2D::doTriangulation\n";
//             return -1;
//         }
//         if(theDomain->addElement(theEle) == false) {
//             opserr<<"WARNING: failed to add element to domain -- ";
//             opserr<<" -- delaunay2D::doTriangulation\n";
//             return -1;
//         }
//         eletags(i) = tag;
//     }
//     eles = eletags;

//     return 0;
    
// }

int
PFEMMesher3D::doTriangulation(int startnode, int endnode, double alpha, 
                              int startanode, int endanode, Domain* theDomain, 
                              ID& eles)
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
    int numnodes = endnode-startnode+1;
    int numadd = endanode-startanode+1;
    if(numnodes < 0) numnodes = 0;
    if(numadd < 0) numadd = 0;
    in.numberofpoints = numnodes + numadd;
    if(in.numberofpoints < nptet) {
        return 0;
    }

    // get space for pointlist
    in.pointlist = new double[in.numberofpoints*ndm];
    if(in.pointlist == 0) {
        opserr<<"WARNING: no enough memory -- PFEMMesher3D::doTriangulation\n";
        return -1;
    }

    // copy nodes to pointlist, use current coordinates
    int loc = 0;
    int numpoints = 0;
    ID point2node(0, in.numberofpoints);
    for(int i=0; i<in.numberofpoints; i++) {

        // pointer to node
        int ndtag = 0;
        if(i<numnodes) {
            ndtag = startnode + i;
        } else {
            ndtag = startanode + i - numnodes;
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
        double x = coord(0);        // initial
        double y = coord(1);
        double z = coord(2);
        x += disp(0);               // current
        y += disp(1);
        z += disp(2);

        in.pointlist[loc++] = x;
        in.pointlist[loc++] = y;
        in.pointlist[loc++] = z;

        point2node[numpoints++] = ndtag;
    }
    in.numberofpoints = numpoints;

    // Delaunay Triangulation
    char s[] = "Qvz";
    tetrahedralize(s, &in, &out);

    // no outputs
    if(out.numberoftetrahedra < 1) {
        return 0;
    }

    // do alpha shape test
    if(alpha > 0) {

        // radius and average size of tetrahedra
        double avesize = 0.0;
        Vector radius(out.numberoftetrahedra);
        for(int i=0; i<out.numberoftetrahedra; i++) {

            // circumcenter of tetrahedra
            double& xc = out.vpointlist[ndm*i];
            double& yc = out.vpointlist[ndm*i+1];
            double& zc = out.vpointlist[ndm*i+2];

            // tetrahedra points
            ID pt(nptet);
            for(int j=0; j<nptet; j++) {
                pt(j) = out.tetrahedronlist[out.numberofcorners*i+j];
            }

            // nodal coordinates
            Vector x(nptet), y(nptet), z(nptet);
            for(int j=0; j<nptet; j++) {
                x(j) = out.pointlist[ndm*pt(j)];
                y(j) = out.pointlist[ndm*pt(j)+1];
                z(j) = out.pointlist[ndm*pt(j)+2];
            }

            // size of tetrahedra
            double he = -1.0;
            for(int j=0; j<nptet; j++) {
                for(int k=j+1; k<nptet; k++) {
                    double h = (x[j]-x[k])*(x[j]-x[k])+(y[j]-y[k])*(y[j]-y[k])+(z[j]-z[k])*(z[j]-z[k]);
                    if(h<he || he==-1.0) {
                        he = h;
                    }
                }
            }
            avesize += sqrt(he);
            // double a1 =  det(x[1], y[1], z[1], x[2], y[2], z[2], x[3], y[3], z[3]);
            // double b1 = -det( 1.0, y[1], z[1],  1.0, y[2], z[2],  1.0, y[3], z[3]);
            // double c1 = -det(x[1],  1.0, z[1], x[2],  1.0, z[2], x[3],  1.0, z[3]);
            // double d1 = -det(x[1], y[1],  1.0, x[2], y[2],  1.0, x[3], y[3],  1.0);
            // double volume = (a1+b1*x[0]+c1*y[0]+d1*z[0])/6.0;
            // avesize += pow(volume, 1.0/3.0);

            // radius
            radius(i) = sqrt((xc-x[0])*(xc-x[0])+(yc-y[0])*(yc-y[0])+(zc-z[0])*(zc-z[0]));
        }
        avesize /= out.numberoftetrahedra;

        // alpha test
        int num = 0;
        eles.resize(out.numberoftetrahedra*nptet);
        for(int i=0; i<out.numberoftetrahedra; i++) {
            if(radius(i) / avesize <= alpha) {
                // pass the test

                // check if all nodes are additional nodes
                bool add = true;
                for(int j=0; j<nptet; j++) {
                    int tag = point2node(out.tetrahedronlist[out.numberofcorners*i+j]);
                    if(tag>=startnode && tag<=endnode) {
                        add = false;
                        break;
                    }
                }
                    
                // add ele
                if(!add) {
                    for(int j=0; j<nptet; j++) {
                        int tag = point2node(out.tetrahedronlist[out.numberofcorners*i+j]);
                        eles(num++) = tag;
                    }
                }
            }
        }
        if(num == 0) {
            eles = ID();
        } else {
            eles.resize(num);
        }
        
    } else if(alpha < 0) {

        // eles
        eles.resize(out.numberoftetrahedra*nptet);

        // copy triangles
        for(int i=0; i<out.numberoftetrahedra; i++) {
            for(int j=0; j<nptet; j++) {
                int tag = point2node(out.tetrahedronlist[out.numberofcorners*i+j]);
                eles(nptet*i+j) = tag;
            }
        }

    }
    
    return 0;
}

// int 
// PFEMMesher3D::save(const char* filename, const ID& nodes, Domain* theDomain)
// {
//     if(theDomain == 0) {
//         opserr<<"WARNING: null domain";
//         opserr<<" -- PFEMMesher3D::save\n";
//         return -1;
//     }

//     // write nodes
//     std::ofstream file2(std::string(filename).append(".node").c_str());
//     file2<<theDomain->getCurrentTime()<<" "<<0<<" "<<0<<" ";
//     file2<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<"\n";

//     std::map<int,int> nodeIndex;
//     int index = 0;

//     // pressure constraints
//     Pressure_ConstraintIter& thePCs = theDomain->getPCs();
//     Pressure_Constraint* thePC = 0;
//     while((thePC = thePCs()) != 0) {
//         int ntag = thePC->getTag();
//         Node* node = theDomain->getNode(ntag);
//         int ptag = thePC->getPressureNode();
//         Node* pnode = theDomain->getNode(ptag);
//         if(node!=0 || pnode!=0) {
//             nodeIndex[ntag] = index++;
//             const Vector& coord = node->getCrds();
//             const Vector& disp = node->getDisp();
//             const Vector& vel = node->getVel();
//             const Vector& accel = node->getAccel();
//             const Vector& pvel = pnode->getVel();
//             file2<<coord(0)+disp(0)<<" "<<coord(1)+disp(1)<<" "<<1<<" ";
//             file2<<vel(0)<<" "<<vel(1)<<" "<<accel(0)<<" "<<accel(1);
//             file2<<" "<<pvel(2)<<"\n";
//         }
//     }

//     // nodes
//     for(int i=0; i<nodes.Size(); i++) {
//         int tag = nodes(i);
//         Node* theNode = theDomain->getNode(tag);
//         if(theNode != 0) {
//             nodeIndex[tag] = index++;
//             const Vector& coord = theNode->getCrds();
//             const Vector& disp = theNode->getDisp();
//             const Vector& vel = theNode->getVel();
//             const Vector& accel = theNode->getAccel();
//             file2<<coord(0)+disp(0)<<" "<<coord(1)+disp(1)<<" "<<2<<" ";
//             file2<<vel(0)<<" "<<vel(1)<<" "<<accel(0)<<" "<<accel(1);
//             file2<<" "<<0<<"\n";

//         }
//     }
//     file2.close();

//     // write edges
//     std::ofstream file1(std::string(filename).append(".edge").c_str());
//     typedef std::set< std::pair<int,int> > EdgeSet;
//     typedef EdgeSet::iterator EdgeSetIter;
//     typedef std::pair<EdgeSetIter, bool> EdgeSetRet;
//     EdgeSet edges;
//     EdgeSetRet ret;
//     ElementIter& theEles = theDomain->getElements();
//     Element* theEle = 0;
//     while((theEle = theEles()) != 0) {
//         const ID& allntags = theEle->getExternalNodes();
//         ID ntags(0,allntags.Size());
//         for(int i=0; i<allntags.Size(); i++) {
//             if(nodeIndex.find(allntags(i)) != nodeIndex.end()) {
//                 ntags[ntags.Size()] = allntags(i);
//             }
//         }
//         if(ntags.Size()==1) continue;
//         if(ntags.Size()>2) {
//             ntags[ntags.Size()] = ntags(0);
//         }
//         for(int i=0; i<ntags.Size()-1; i++) {
//             ID nd(0,2);
//             for(int j=0; j<2; j++) {
//                 nd.insert(ntags(i+j));
//             }
//             ret = edges.insert(std::make_pair(nd(0),nd(1)));
//             if(ret.second == true) {
//                 file1<<nodeIndex[nd(0)]<<" "<<nodeIndex[nd(1)]<<"\n";
//             }
//         }
//     }
//     file1.close();

//     return 0;
// }

// void
// PFEMMesher3D::setBoundary(double x1, double y1, double x2, double y2)
// {
//     bound(0) = x1;
//     bound(1) = y1;
//     bound(2) = x2;
//     bound(3) = y2;
// }

// void 
// PFEMMesher3D::removeOutBoundNodes(const ID& nodes, ID& nodes2, Domain* theDomain)
// {
//     // check nodes
//     ID removelist(0, nodes.Size());
//     for(int i=0; i<nodes.Size(); i++) {
//         int tag = nodes(i);
//         Node* theNode = theDomain->getNode(tag);
//         const Vector& coord = theNode->getCrds();
//         if(coord.Size() < 2) {
//             continue;
//         }
//         const Vector& disp = theNode->getTrialDisp();
//         if(disp.Size() < 2) {
//             continue;
//         }
//         double x = coord(0);        // initial
//         double y = coord(1);
//         x += disp(0);               // current
//         y += disp(1);
        
//         // if out of boundary
//         if(x<bound(0) || x>bound(2) || y<bound(1) || y>bound(3)) {
//             removelist.insert(tag);
//         } else {
//             nodes2.insert(tag);
//         }
//     }

//     // remove nodes
//     for(int i=0; i<removelist.Size(); i++) {
//         Node* theNode = theDomain->removeNode(removelist(i));
//         if(theNode != 0) {
//             delete theNode;
//         }
//         Pressure_Constraint* thePC = theDomain->removePressure_Constraint(removelist(i));
//         if(thePC != 0) {
//             delete thePC;
//         }
//     }
// }

// int 
// PFEMMesher3D::findNodeTag(Domain* theDomain)
// {
//     NodeIter& theNodes = theDomain->getNodes();
//     Node* firstNode = theNodes();
//     if(firstNode == 0) {
//         return -1;
//     }
//     int tag = firstNode->getTag();
//     return tag-1;
// }

// int 
// PFEMMesher3D::findEleTag(Domain* theDomain)
// {
//     ElementIter& theEles = theDomain->getElements();
//     Element* firstEle = theEles();
//     if(firstEle == 0) {
//         return -1;
//     }
//     int tag = firstEle->getTag();
//     return tag-1;
// }

