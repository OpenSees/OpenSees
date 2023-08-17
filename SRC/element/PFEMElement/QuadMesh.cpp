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
// Description: The class QuadMesh is for meshing quadangles
//

#include "QuadMesh.h"
#include <Node.h>
#include <NodeIter.h>
#include <Domain.h>
#include <string.h>
#include "QuadMeshGenerator.h"
#include <elementAPI.h>
#include <Element.h>
#include "LineMesh.h"
#include <set>
#include <map>
#include <vector>
#include <Pressure_Constraint.h>
#include "BackgroundDef.h"
#include <cmath>

int OPS_QuadMesh() {
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
    QuadMesh *mesh = new QuadMesh(idata[0]);
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
        opserr << "WARNING: failed to do quad mesh\n";
        return -1;
    }

    return 0;
}

QuadMesh::QuadMesh(int tag)
        : Mesh(tag, 4), ltags() {
}

QuadMesh::~QuadMesh() {
}

int
QuadMesh::mesh() {
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
    ID ndtags;
    std::vector<ID> lines;
    for (int i = 0; i < ltags.Size(); ++i) {
        LineMesh *line = dynamic_cast<LineMesh *>(OPS_getMesh(ltags(i)));
        if (line == 0) {
            opserr << "WARNING: mesh " << ltags(i) << " is not a line\n";
            return -1;
        }
        const ID &tags = line->getNodeTags();
        const ID &newtags = line->getNewNodeTags();

        ID linenodes(tags.Size() + newtags.Size());
        linenodes(0) = tags(0);
        linenodes(linenodes.Size() - 1) = tags(1);
        for (int j = 0; j < tags.Size(); ++j) {
            ndtags.insert(tags(j));
        }
        for (int j = 0; j < newtags.Size(); ++j) {
            ndtags.insert(newtags(j));
            linenodes(j + 1) = newtags(j);
        }
        lines.push_back(linenodes);
    }
    if (ndtags.Size() < 4) {
        opserr << "WARNING: input number of nodes < 4\n";
        return -1;
    }
    if (lines.size() != 4) {
        opserr << "WARNING: input number of lines != 4\n";
        return -1;
    }
    this->setNodeTags(ndtags);

    // get point index
    for (unsigned int i = 0; i < lines.size(); ++i) {
        for (int j = 0; j < lines[i].Size(); ++j) {
            int tag = lines[i](j);
            lines[i](j) = ndtags.getLocationOrdered(tag);
        }
    }

    // calling mesh generator
    QuadMeshGenerator gen;
    int nodecounter = nextNodeTag();
    int ndm = OPS_GetNDM();
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

        Vector vcrds(ndm);
        if (disp.Size() >= ndm) {
            for (int j = 0; j < ndm; ++j) {
                vcrds(j) = crds(j) + disp(j);
            }
        }

        // add point
        gen.addPoint(vcrds);

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
                    pnode = new Node(nodecounter++, 1, vcrds(0), vcrds(1));
                } else if (ndm == 3) {
                    pnode = new Node(nodecounter++, 1, vcrds(0), vcrds(1), vcrds(2));
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
    for (unsigned int i = 0; i < lines.size(); ++i) {
        // add line
        gen.addLine(lines[i]);
    }

    // meshing
    if (gen.mesh(size) < 0) return -1;

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
        gen.getPoint(i, crds);

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

    // get quads
    int numquad = gen.getNumQuads();
    if (numquad == 0) return 0;
    ID elenodes(numquad * 4);
    if (this->getNumEleNodes() == 4) {

        for (int i = 0; i < numquad; i++) {
            ID pts;
            gen.getQuad(i, pts);
            for (int j = 0; j < pts.Size(); ++j) {
                elenodes(4 * i + j) = allndtags(pts(j));
            }
        }
    } else if (this->getNumEleNodes() == 2) {

        elenodes.resize(numquad * 10);

        for (int i = 0; i < numquad; i++) {
            ID pts;
            gen.getQuad(i, pts);
            elenodes(10 * i) = allndtags(pts(0));
            elenodes(10 * i + 1) = allndtags(pts(1));
            elenodes(10 * i + 2) = allndtags(pts(1));
            elenodes(10 * i + 3) = allndtags(pts(2));
            elenodes(10 * i + 4) = allndtags(pts(2));
            elenodes(10 * i + 5) = allndtags(pts(3));
            elenodes(10 * i + 6) = allndtags(pts(3));
            elenodes(10 * i + 7) = allndtags(pts(0));
            elenodes(10 * i + 8) = allndtags(pts(0));
            elenodes(10 * i + 9) = allndtags(pts(2));
        }

    } else if (this->getNumEleNodes() == 3) {

        elenodes.resize(numquad * 6);

        for (int i = 0; i < numquad; i++) {
            ID pts;
            gen.getQuad(i, pts);
            elenodes(6 * i) = allndtags(pts(0));
            elenodes(6 * i + 1) = allndtags(pts(1));
            elenodes(6 * i + 2) = allndtags(pts(2));
            elenodes(6 * i + 3) = allndtags(pts(0));
            elenodes(6 * i + 4) = allndtags(pts(2));
            elenodes(6 * i + 5) = allndtags(pts(3));
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


