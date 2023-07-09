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
// $Date: 2015-10-19 $

// Written: Minjie Zhu
// Created: Oct. 19
//
// Description: The class TriMesh is for meshing triangles
//

#include "TriMesh.h"
#include <Node.h>
#include <NodeIter.h>
#include <Domain.h>
#include <string.h>
#include "TriangleMeshGenerator.h"
#include <elementAPI.h>
#include <Element.h>
#include "LineMesh.h"
#include <set>
#include <map>
#include <vector>
#include <Pressure_Constraint.h>
#include "BackgroundDef.h"
#include <cmath>

int OPS_TriMesh() {
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "WARNING: want tag? numlines? ltags? id? ndf? size? eleType? eleArgs?\n";
        return -1;
    }

    // get tag and number lines
    int num = 2;
    int idata[2];
    if (OPS_GetIntInput(&num, idata) < 0) {
        opserr << "WARNING: failed to read mesh tag and number of lines\n";
        return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < idata[1] + 3) {
        opserr << "WARNING: want ltags? id? ndf? size? <eleType? eleArgs?>\n";
        return -1;
    }

    // create mesh
    TriMesh *mesh = new TriMesh(idata[0]);
    if (OPS_addMesh(mesh) == false) {
        opserr << "WARNING: failed to add mesh\n";
        return -1;
    }

    // get line tags
    num = idata[1];
    ID ltags(num);
    if (OPS_GetIntInput(&num, &ltags(0)) < 0) {
        opserr << "WARNING: failed to read line tags\n";
        return -1;
    }
    mesh->setLineTags(ltags);

    // get id and ndf
    num = 2;
    int data[2];
    if (OPS_GetIntInput(&num, data) < 0) {
        opserr << "WARNING: failed to read id and ndf\n";
        return -1;
    }
    mesh->setID(data[0]);
    mesh->setNdf(data[1]);

    // get size
    num = 1;
    double size;
    if (OPS_GetDoubleInput(&num, &size) < 0) {
        opserr << "WARNING: failed to read mesh size\n";
        return -1;
    }
    mesh->setMeshsize(size);

    // set eleArgs
    if (mesh->setEleArgs() < 0) {
        opserr << "WARNING: failed to set element arguments\n";
        return -1;
    }

    // mesh
    if (mesh->mesh() < 0) {
        opserr << "WARNING: failed to do triangular mesh\n";
        return -1;
    }

    return 0;
}

TriMesh::TriMesh(int tag)
        : Mesh(tag, 3), ltags() {
}

TriMesh::~TriMesh() {
}

int
TriMesh::mesh() {
    Domain *domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return -1;
    }

    // check
    double size = this->getMeshsize();
    if (size <= 0) {
        opserr << "WARNING: mesh size <= 0\n";
        return -1;
    }
    if (ltags.Size() == 0) return 0;

    // get nodes and elements from lines
    ID ndtags, lndtags;
    std::set<ID> elenodeset;
    for (int i = 0; i < ltags.Size(); ++i) {
        LineMesh *line = dynamic_cast<LineMesh *>(OPS_getMesh(ltags(i)));
        if (line == 0) {
            opserr << "WARNING: mesh " << ltags(i) << " is not a line\n";
            return -1;
        }
        const ID &tags = line->getNodeTags();
        for (int j = 0; j < tags.Size(); ++j) {
            ndtags.insert(tags(j));
            lndtags.insert(tags(j));
        }
        const ID &newtags = line->getNewNodeTags();
        for (int j = 0; j < newtags.Size(); ++j) {
            ndtags.insert(newtags(j));
        }
        const ID &elends = line->getEleNodes();
        for (int j = 0; j < elends.Size() / 2; ++j) {
            ID twonds(2);
            twonds(0) = elends(2 * j);
            twonds(1) = elends(2 * j + 1);
            elenodeset.insert(twonds);
        }
    }
    if (ndtags.Size() < 3) {
        opserr << "WARNING: input number of nodes < 3\n";
        return -1;
    }
    if (lndtags.Size() < 3) {
        opserr << "WARNING: input number of lnodes < 3\n";
        return -1;
    }
    if (elenodeset.size() < 3) {
        opserr << "WARNING: input number of segments < 3\n";
        return -1;
    }
    this->setNodeTags(ndtags);

    // direction of the plane if 3D
    int ndm = OPS_GetNDM();
    VVDouble dir(ndm);
    if (ndm == 3) {
        for (int j = 0; j < lndtags.Size() - 2; ++j) {

            for (int i = 0; i < ndm; ++i) {
                // get node
                Node *theNode = domain->getNode(lndtags[j + i]);
                if (theNode == 0) {
                    opserr << "WARNING: node " << lndtags(j + i) << " is not defined\n";
                    return -1;
                }
                const Vector &crds = theNode->getCrds();
                const Vector &disp = theNode->getTrialDisp();
                if (crds.Size() < ndm) {
                    opserr << "WARNING: ndm < 3\n";
                    return -1;
                }
                if (disp.Size() >= ndm) {
                    dir[i].resize(ndm, 0.0);
                    for (int k = 0; k < ndm; ++k) {
                        dir[i][k] = disp(k) + crds(k);
                    }
                }
            }

            // direction 2->0, 2->1
            dir[0] -= dir[2];
            dir[1] -= dir[2];

            // check if same direction
            if (fabs(dotVDouble(dir[0], dir[1]) - 1) < 1e-3) {
                continue;
            }

            // normal direction
            VDouble nor;
            crossVDouble(dir[0], dir[1], nor);

            // direction perpenticular to 2->0 and norm
            crossVDouble(dir[0], nor, dir[1]);

            // normalize
            dir[0] /= normVDouble(dir[0]);
            dir[1] /= normVDouble(dir[1]);

            break;
        }
    }

    // get segments from lines
    ID segs(elenodeset.size() * 2);
    int iseg = 0;
    for (std::set<ID>::iterator it = elenodeset.begin(); it != elenodeset.end(); ++it) {
        const ID &twonds = *it;
        segs(2 * iseg) = ndtags.getLocationOrdered(twonds(0));
        segs(2 * iseg + 1) = ndtags.getLocationOrdered(twonds(1));
        ++iseg;
    }

    // calling mesh generator
    TriangleMeshGenerator gen;
    int nodecounter = nextNodeTag();
    for (int i = 0; i < ndtags.Size(); i++) {

        // get node
        Node *theNode = domain->getNode(ndtags(i));
        if (theNode == 0) {
            opserr << "WARNING: node " << ndtags(i) << " is not defined\n";
            return -1;
        }
        const Vector &crds = theNode->getCrds();
        const Vector &disp = theNode->getTrialDisp();
        if (crds.Size() < ndm) {
            opserr << "WARNING: ndm < " << ndm << "\n";
            return -1;
        }

        VDouble vcrds(ndm);
        if (disp.Size() >= ndm) {
            for (int j = 0; j < ndm; ++j) {
                vcrds[j] = crds(j) + disp(j);
            }
        }

        // x and y
        double x = vcrds[0];
        double y = vcrds[1];
        if (ndm == 3) {
            VDouble temp = vcrds;
            temp -= dir[2];
            x = dotVDouble(temp, dir[0]);
            y = dotVDouble(temp, dir[1]);
        }

        // add point
        gen.addPoint(x, y);

        // add pc
        if (this->isFluid()) {
            // create pressure constraint
            Pressure_Constraint *thePC = domain->getPressure_Constraint(ndtags(i));
            if (thePC != 0) {
                thePC->setDomain(domain);
            } else {

                // create pressure node
                Node *pnode = 0;
                if (ndm == 2) {
                    pnode = new Node(nodecounter++, 1, vcrds[0], vcrds[1]);
                } else if (ndm == 3) {
                    pnode = new Node(nodecounter++, 1, vcrds[0], vcrds[1], vcrds[2]);
                }
                if (pnode == 0) {
                    opserr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
                    return -1;
                }
                if (domain->addNode(pnode) == false) {
                    opserr << "WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
                    delete pnode;
                    return -1;
                }

                thePC = new Pressure_Constraint(ndtags(i), pnode->getTag());
                if (thePC == 0) {
                    opserr << "WARNING: no enough memory for Pressure_Constraint\n";
                    return -1;
                }
                if (domain->addPressure_Constraint(thePC) == false) {
                    opserr << "WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
                    delete thePC;
                    return -1;
                }
            }
        }
    }
    for (int i = 0; i < segs.Size() / 2; ++i) {
        // add segment
        gen.addSegment(segs(2 * i), segs(2 * i + 1), 0);
    }

    // meshing
    gen.mesh(size * size * 0.5, false);

    // get points and create nodes
    int nump = gen.getNumPoints();
    if (nump == 0) {
        opserr << "WARNING: no nodes is meshed\n";
        return -1;
    }
    ID newndtags(nump - ndtags.Size());
    ID allndtags(nump);
    for (int i = 0; i < ndtags.Size(); ++i) {
        allndtags(i) = ndtags(i);
    }

    for (int i = ndtags.Size(); i < nump; ++i) {
        // get point
        Vector crds(ndm);
        int mark = 0;
        gen.getPoint(i, crds(0), crds(1), mark);

        // 3D transformatin
        if (ndm == 3) {
            VDouble vcrds = dir[2];
            VDouble incr = dir[0];
            incr *= crds(0);
            vcrds += incr;
            incr = dir[1];
            incr *= crds(1);
            vcrds += incr;
            toVector(vcrds, crds);
        }

        // create node
        Node *node = newNode(nodecounter++, crds);
        if (node == 0) {
            opserr << "WARNING: failed to create node\n";
            return -1;
        }
        if (domain->addNode(node) == false) {
            opserr << "WARNING: failed to add node to domain\n";
            delete node;
            return -1;
        }
        allndtags(i) = node->getTag();
        newndtags(i - ndtags.Size()) = node->getTag();

        // add pc
        if (this->isFluid()) {
            // create pressure constraint
            Pressure_Constraint *thePC = domain->getPressure_Constraint(node->getTag());
            if (thePC != 0) {
                thePC->setDomain(domain);
            } else {

                // create pressure node
                Node *pnode = 0;
                if (ndm == 2) {
                    pnode = new Node(nodecounter++, 1, crds(0), crds(1));
                } else {
                    pnode = new Node(nodecounter++, 1, crds(0), crds(1), crds(2));
                }
                if (pnode == 0) {
                    opserr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
                    return -1;
                }
                if (domain->addNode(pnode) == false) {
                    opserr << "WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
                    delete pnode;
                    return -1;
                }

                thePC = new Pressure_Constraint(node->getTag(), pnode->getTag());
                if (thePC == 0) {
                    opserr << "WARNING: no enough memory for Pressure_Constraint\n";
                    return -1;
                }
                if (domain->addPressure_Constraint(thePC) == false) {
                    opserr << "WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
                    delete thePC;
                    return -1;
                }
            }
        }
    }
    this->setNewNodeTags(newndtags);

    // get triangles
    int numtri = gen.getNumTriangles();
    if (numtri == 0) return 0;

    ID elenodes;
    if (this->getNumEleNodes() == 3) {
        elenodes.resize(numtri * 3);

        for (int i = 0; i < numtri; i++) {
            int p1, p2, p3;
            gen.getTriangle(i, p1, p2, p3);
            elenodes(3 * i) = allndtags(p1);
            elenodes(3 * i + 1) = allndtags(p2);
            elenodes(3 * i + 2) = allndtags(p3);
        }
    } else if (this->getNumEleNodes() == 2) {
        elenodes.resize(numtri * 6);
        for (int i = 0; i < numtri; i++) {
            int p1, p2, p3;
            gen.getTriangle(i, p1, p2, p3);
            elenodes(6 * i) = allndtags(p1);
            elenodes(6 * i + 1) = allndtags(p2);
            elenodes(6 * i + 2) = allndtags(p2);
            elenodes(6 * i + 3) = allndtags(p3);
            elenodes(6 * i + 4) = allndtags(p1);
            elenodes(6 * i + 5) = allndtags(p3);
        }
    }
    this->setEleNodes(elenodes);

    // create elements
    if (this->newElements(elenodes) < 0) {
        opserr << "WARNING: failed to create elements\n";
        return -1;
    }

    return 0;
}


// Every node may be in multiple mesh
// Do triangulation of all nodes
// For each triangle:
//     If 3 nodes share 1 mesh, 
//        use that mesh for the triangle
//     If 3 nodes share multiple mesh, 
//        use the mesh with eleArgs and lowest id
//     If 3 nodes are in different mesh, 
//        use the mesh with lowest id
//     If the selected mesh id >= 0, skip triangle
//     If the selected mesh has no eleArgs, skip triangle
int
TriMesh::remesh(double alpha) {
    if (OPS_GetNDM() != 2) {
        opserr << "WARNING: TriMesh::remesh is only for 2D problem\n";
        return -1;
    }
    Domain *domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return -1;
    }

    // get all nodes
    std::map<int, std::vector<int> > ndinfo;
    TaggedObjectIter &meshes = OPS_getAllMesh();
    ID nodetags;
    Mesh *msh = 0;
    while ((msh = dynamic_cast<Mesh *>(meshes())) != 0) {

        int id = msh->getID();
        int mtag = msh->getTag();

        if (id == 0) continue;

        // get bound nodes
        const ID &tags = msh->getNodeTags();
        for (int i = 0; i < tags.Size(); ++i) {
            std::vector<int> &info = ndinfo[tags(i)];
            info.push_back(mtag);
            info.push_back(id);
            info.push_back(msh->hasEleArgs());
            nodetags.insert(tags(i));
        }

        // get internal nodes
        const ID &newtags = msh->getNewNodeTags();
        for (int i = 0; i < newtags.Size(); ++i) {
            std::vector<int> &info = ndinfo[newtags(i)];
            info.push_back(mtag);
            info.push_back(id);
            info.push_back(msh->hasEleArgs());
            nodetags.insert(newtags(i));
        }
    }

    if (nodetags.Size() == 0) return 0;

    // calling mesh generator
    TriangleMeshGenerator gen;
    int nodecounter = nextNodeTag();
    for (int i = 0; i < nodetags.Size(); ++i) {

        // get node
        Node *theNode = domain->getNode(nodetags(i));
        if (theNode == 0) {
            opserr << "WARNING: node " << nodetags(i) << " is not defined\n";
            return -1;
        }
        const Vector &crds = theNode->getCrds();
        const Vector &disp = theNode->getTrialDisp();
        if (crds.Size() < 2 || disp.Size() < 2) {
            opserr << "WARNING: ndm < 2 or ndf < 2\n";
            return -1;
        }

        // create pc if not
        Pressure_Constraint *thePC = domain->getPressure_Constraint(nodetags(i));
        if (thePC != 0) {
            thePC->setDomain(domain);
        } else {

            // create pressure node
            Node *pnode = 0;
            pnode = new Node(nodecounter++, 1, crds[0], crds[1]);

            if (pnode == 0) {
                opserr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
                return -1;
            }
            if (domain->addNode(pnode) == false) {
                opserr << "WARNING: failed to add node to domain -- BgMesh::gridNodes\n";
                delete pnode;
                return -1;
            }

            thePC = new Pressure_Constraint(nodetags(i), pnode->getTag());
            if (thePC == 0) {
                opserr << "WARNING: no enough memory for Pressure_Constraint\n";
                return -1;
            }
            if (domain->addPressure_Constraint(thePC) == false) {
                opserr << "WARNING: failed to add PC to domain -- BgMesh::gridNodes\n";
                delete thePC;
                return -1;
            }
        }

        // add point
        gen.addPoint(crds(0) + disp(0), crds(1) + disp(1));
    }

    // meshing
    gen.remesh(alpha);

    // get elenodes
    std::map<int, ID> meshelenodes;
    for (int i = 0; i < gen.getNumTriangles(); i++) {

        // get points
        int p1, p2, p3;
        gen.getTriangle(i, p1, p2, p3);

        // get nodes
        int nds[3];
        nds[0] = nodetags(p1);
        nds[1] = nodetags(p2);
        nds[2] = nodetags(p3);

        // check if all nodes in same mesh
        std::vector<int> &info1 = ndinfo[nds[0]];
        int mtag = 0, id = 0, hasele = false;
        int same_mtag = 0, same_id = 0, same_hasele=false;
        bool same = false;
        for (int k = 0; k < (int) info1.size() / 3; ++k) {
            // check if any mesh of node 1 is same for another two nodes
            mtag = info1[3 * k];
            id = info1[3 * k + 1];
            hasele = info1[3 * k + 2];

            int num = 0;
            for (int j = 1; j < 3; ++j) {
                std::vector<int> &infoj = ndinfo[nds[j]];
                for (int kj = 0; kj < (int) infoj.size() / 3; ++kj) {
                    int mtagj = infoj[3 * kj];
                    if (mtag == mtagj) {
                        ++num;
                        break;
                    }
                }

            }
            if (num == 2) {
                // get first same
                if (!same) {
                    same_mtag = mtag;
                    same_id = id;
                    same_hasele = hasele;
                    same = true;
                } else {
                    // if already has same with ele, then next same has to have smaller id
                    if (same_hasele) {
                        if (same_id > id) {
                            same_mtag = mtag;
                            same_id = id;
                            same_hasele = hasele;
                        }
                    } else {
                        // if already has same without ele, then next same has ele 
                        if (hasele) {
                            same_mtag = mtag;
                            same_id = id;
                            same_hasele = hasele;
                        } else {
                            // if already has same without ele, then next same without ele must have smaller id
                            if (same_id > id) {
                                same_mtag = mtag;
                                same_id = id;
                                same_hasele = hasele;
                            }
                        }
                    }
                }
            }
        }

        // nodes in different mesh
        if (same) {
            mtag = same_mtag;
            id = same_id;
            hasele = same_hasele;
        } else {
            // find the mesh with lowest id
            mtag = 0;
            id = 0;
            hasele = false;
            for (int j = 0; j < 3; ++j) {
                std::vector<int> &info = ndinfo[nds[j]];
                for (int k = 0; k < (int) info.size() / 3; ++k) {
                    if (!info[3 * k + 2]) {
                        continue;
                    }
                    if (id == 0 || info[3 * k + 1] < id) {
                        mtag = info[3 * k];
                        id = info[3 * k + 1];
                        hasele = info[3 * k + 2];
                    }
                }
            }
        }

        // if all connected to structure
        if (id >= 0) continue;

        // if no ele is associated
        if (!hasele) continue;

        // add elenodes to its mesh
        ID &elenodes = meshelenodes[mtag];
        for (int j = 0; j < 3; ++j) {
            elenodes[elenodes.Size()] = nds[j];
        }
    }

    // creat elements
    for (std::map<int, ID>::iterator it = meshelenodes.begin();
         it != meshelenodes.end(); ++it) {

        int mtag = it->first;
        ID &elenodes = it->second;

        msh = OPS_getMesh(mtag);
        if (msh != 0) {

            int id = msh->getID();

            // remove mesh for id<0
            if (id < 0) {
                if (msh->clearEles() < 0) {
                    opserr << "WARNING: failed to clear element in mesh" << mtag << "\n";
                    return -1;
                }

                if (msh->newElements(elenodes) < 0) {
                    opserr << "WARNING: failed to create new elements in mesh" << mtag << "\n";
                    return -1;
                }
            }
        }
    }

    return 0;
}
