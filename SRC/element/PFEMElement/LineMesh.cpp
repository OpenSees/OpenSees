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
// Description: The class LineMesh is for meshing lines.
//

#include "LineMesh.h"
#include <Node.h>
#include <NodeIter.h>
#include <Domain.h>
#include <elementAPI.h>
#include <Element.h>
#include <LineMeshGenerator.h>

int OPS_LineMesh() {
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "WARNING: want tag? numnodes? ndtags? id? ndf? size? <eleType? eleArgs?>\n";
        return -1;
    }

    // get tag and number of nodes
    int num = 2;
    int idata[2];
    if (OPS_GetIntInput(&num, idata) < 0) {
        opserr << "WARNING: failed to read mesh tag, and number of nodes\n";
        return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < idata[1] + 3) {
        opserr << "WARNING: want ndtags? id? ndf? size? <eleType? eleArgs?>\n";
        return -1;
    }

    // create mesh
    LineMesh *mesh = new LineMesh(idata[0]);
    if (OPS_addMesh(mesh) == false) {
        opserr << "WARNING: failed to add mesh\n";
        return -1;
    }

    // get node tags
    num = idata[1];
    ID ndtags(num);
    if (OPS_GetIntInput(&num, &ndtags(0)) < 0) {
        opserr << "WARNING: failed to read node tags\n";
        return -1;
    }
    mesh->setNodeTags(ndtags);

    // get id and ndf
    num = 2;
    int id[2];
    if (OPS_GetIntInput(&num, id) < 0) {
        opserr << "WARNING: failed to read mesh id and node ndf\n";
        return -1;
    }
    mesh->setID(id[0]);
    mesh->setNdf(id[1]);

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
        opserr << "WARNING: failed to mesh line\n";
        return -1;
    }

    return 0;
}

LineMesh::LineMesh(int tag)
        : Mesh(tag, 2) {
}

LineMesh::~LineMesh() {
}


int
LineMesh::mesh() {
    Domain *domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "WARNING: domain is not created\n";
        return -1;
    }

    // get data
    const ID &tags = this->getNodeTags();
    double meshsize = this->getMeshsize();

    // line mesh generator
    LineMeshGenerator gen;

    // add points and lines
    for (int i = 0; i < tags.Size(); ++i) {
        Node *node = domain->getNode(tags(i));
        if (node == 0) {
            opserr << "WARNING: node " << tags(i) << " does not exist\n";
            return -1;
        }

        // current crds of the node
        Vector crds = node->getCrds();
        const Vector &disp = node->getTrialDisp();
        if (disp.Size() >= crds.Size()) {
            for (int j = 0; j < crds.Size(); ++j) {
                crds(j) += disp(j);
            }
        }

        // add the node
        gen.addPoint(crds);

        // add the line
        if (i > 0) {
            ID line(2);
            line(0) = i - 1;
            line(1) = i;
            gen.addLine(line);
        }
    }

    // mesh the lines
    if (gen.mesh(meshsize) < 0) {
        opserr << "WARNIGN: failed to mesh line\n";
        return -1;
    }

    // get points and create nodes
    int nump = gen.getNumPoints();
    if (nump == 0) {
        opserr << "WARNING: no nodes is meshed\n";
        return -1;
    }
    ID newndtags(nump - tags.Size());
    ID allndtags(nump);
    for (int i = 0; i < tags.Size(); ++i) {
        allndtags(i) = tags(i);
    }

    int nodecounter = nextNodeTag();
    for (int i = tags.Size(); i < nump; ++i) {
        Vector crds;
        gen.getPoint(i, crds);
        Node *node = newNode(nodecounter++, crds);
        if (node == 0) {
            opserr << "WARING: failed to create node\n";
            return -1;
        }
        if (domain->addNode(node) == false) {
            opserr << "WARNING: failed to add node to domain\n";
            delete node;
            return -1;
        }
        allndtags(i) = node->getTag();
        newndtags(i - tags.Size()) = node->getTag();
    }
    this->setNewNodeTags(newndtags);

    // get lines 
    int numlines = gen.getNumLines();
    if (numlines == 0) return 0;
    ID elenodes(numlines * 2);

    for (int i = 0; i < numlines; ++i) {
        ID pts;
        gen.getLine(i, pts);
        elenodes(2 * i) = allndtags(pts(0));
        elenodes(2 * i + 1) = allndtags(pts(1));
        if (elenodes(2 * i) > elenodes(2 * i + 1)) {
            elenodes(2 * i) = allndtags(pts(1));
            elenodes(2 * i + 1) = allndtags(pts(0));
        }

    }
    this->setEleNodes(elenodes);

    // create elemnts
    if (this->newElements(elenodes) < 0) {
        opserr << "WARNING: failed to create elements\n";
        return -1;
    }

    return 0;
}

