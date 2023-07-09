/* ******************************************************************
***
**    OpenSees - Open System for Earthquake Engineering Simulation **
**          Pacific Earthquake Engineering Research Center **
** **
** **
** (C) Copyright 1999, The Regents of the University of California **
** All Rights Reserved. **
** **
** Commercial use of this program without express permission of the **
** University of California, Berkeley, is strictly prohibited.  See **
** file 'COPYRIGHT'  in main directory for information on usage and **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES. **
** **
** Developed by: **
**   Frank McKenna (fmckenna@ce.berkeley.edu) **
**   Gregory L. Fenves (fenves@ce.berkeley.edu) **
**   Filip C. Filippou (filippou@ce.berkeley.edu) **
** **
** ******************************************************************
*/

#include "Flume.h"

#include <Domain.h>
#include <Node.h>
#include <elementAPI.h>

#include <cmath>

int OPS_Flume() {
    // get inputs
    int ndm = OPS_GetNDM();

    if (OPS_GetNumRemainingInputArgs() < 2 * ndm + 4) {
        opserr << "WARNING: want tag?, id?, ndf?, x0?, y0?, <z0?>, "
                  "L?, <B?>, H?, size?, <-top>\n";
        return -1;
    }

    int num = 3;
    std::vector<int> idata(num);
    if (OPS_GetIntInput(&num, &idata[0]) < 0) {
        opserr << "WARNING: failed to get tag, id, ndf \n";
        return -1;
    }
    if (idata[2] <= 0) {
        opserr << "WARING: ndf <= 0\n";
        return -1;
    }

    num = ndm;
    std::vector<double> crds(num);
    if (OPS_GetDoubleInput(&num, &crds[0]) < 0) {
        opserr << "WARNING: failed to get crds\n";
        return -1;
    }

    std::vector<double> data(num);
    if (OPS_GetDoubleInput(&num, &data[0]) < 0) {
        opserr << "WARNING: failed to get data\n";
        return -1;
    }
    for (int i = 0; i < ndm; ++i) {
        if (data[i] <= 0) {
            opserr << "WARNING: dimension cannot <= 0\n";
            return -1;
        }
    }

    num = 1;
    double size;
    if (OPS_GetDoubleInput(&num, &size) < 0) {
        opserr << "WARNING: failed to get size\n";
        return -1;
    }
    if (size <= 0) {
        opserr << "WARNING: size <= 0\n";
        return -1;
    }

    bool top = false;
    if (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();
        if (strcmp(opt, "-top") == 0) {
            top = true;
        }
    }

    // create mesh
    Flume* mesh = new Flume(idata[0], crds, data, top);
    if (OPS_addMesh(mesh) == false) {
        opserr << "WARNING: failed to add mesh\n";
        return -1;
    }

    mesh->setID(idata[1]);
    mesh->setNdf(idata[2]);
    mesh->setMeshsize(size);

    // mesh
    if (mesh->mesh() < 0) {
        opserr << "WARNING: failed to create flume\n";
        return -1;
    }

    return 0;
}

Flume::Flume(int tag, const std::vector<double>& c,
             const std::vector<double>& d, bool t)
    : Mesh(tag, 0), crds(c), dimensions(d), top(t) {}

Flume::~Flume() {}

Node* Flume::create_node(const std::vector<double>& crds,
                         int& nodeTag) {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "no domain - create_node\n";
        return 0;
    }
    Node* node = 0;
    int ndf = this->getNdf();
    if (crds.size() == 2) {
        node = new Node(nodeTag, ndf, crds[0], crds[1]);
    } else if (crds.size() == 3) {
        node = new Node(nodeTag, ndf, crds[0], crds[1], crds[2]);
    }
    if (node == 0) {
        opserr << "failed to create node - create_node\n";
        return 0;
    }
    if (domain->addNode(node) == false) {
        opserr << "WARNING: node " << nodeTag
               << "already exists - create_node\n";
        delete node;
        return 0;
    }
    ++nodeTag;
    return node;
}

int Flume::create_line(Node* nd1, Node* nd2, int& nodeTag, int dir) {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "no domain\n";
        return -1;
    }
    if (dir < 0 || dir > 2) {
        opserr << "dir not correct - create_line\n";
        return -1;
    }
    if (nd1 == 0 || nd2 == 0) {
        opserr << "nd1 or nd2 invalid - create_line\n";
        return -1;
    }

    double size = this->getMeshsize();

    const Vector& crds1 = nd1->getCrds();
    const Vector& crds2 = nd2->getCrds();
    std::vector<double> curr(crds1.Size());
    for (int i = 0; i < crds1.Size(); ++i) {
        curr[i] = crds1(i);
    }
    curr[dir] = crds1[dir] + size;
    while (curr[dir] < crds2[dir] - size / 2.0) {
        if (create_node(curr, nodeTag) == 0) {
            return -1;
        }
        curr[dir] += size;
    }

    return 0;
}

int Flume::create_face(Node* nd1, Node* nd2, int& nodeTag, int dir1,
                       int dir2) {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "no domain - create_face\n";
        return -1;
    }
    if (dir1 < 0 || dir1 > 2) {
        opserr << "dir1 invalid - create_face\n";
        return -1;
    }
    if (dir2 < 0 || dir2 > 2) {
        opserr << "dir2 invalid - create_face\n";
        return -1;
    }
    if (nd1 == 0 || nd2 == 0) {
        opserr << "nd1 or nd2 invalid - create_face\n";
        return -1;
    }

    double size = this->getMeshsize();

    const Vector& crds1 = nd1->getCrds();
    const Vector& crds2 = nd2->getCrds();
    std::vector<double> curr(crds1.Size());
    for (int i = 0; i < crds1.Size(); ++i) {
        curr[i] = crds1(i);
    }
    curr[dir1] = crds1[dir1] + size;
    while (curr[dir1] < crds2[dir1] - size * 0.5) {
        curr[dir2] = crds1[dir2] + size;
        while (curr[dir2] < crds2[dir2] - size * 0.5) {
            if (create_node(curr, nodeTag) == 0) {
                return -1;
            }
            curr[dir2] += size;
        }
        curr[dir1] += size;
    }

    return 0;
}

int Flume::mesh() {
    Domain* domain = OPS_GetDomain();
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

    // create corner nodes
    int nodeTag = nextNodeTag();
    int initNodeTag = nodeTag;
    int ndm = (int)crds.size();
    std::vector<double> curr(ndm);
    std::vector<Node*> nodes;
    if (ndm == 2) {
        for (int j = 0; j < 2; ++j) {
            curr[1] = crds[1] + j * dimensions[1];
            for (int i = 0; i < 2; ++i) {
                curr[0] = crds[0] + i * dimensions[0];
                Node* node = create_node(curr, nodeTag);
                if (node == 0) {
                    return -1;
                }
                nodes.push_back(node);
            }
        }
    } else if (ndm == 3) {
        for (int k = 0; k < 2; ++k) {
            curr[2] = crds[2] + k * dimensions[2];  // H
            for (int j = 0; j < 2; ++j) {
                curr[1] = crds[1] + j * dimensions[1];
                for (int i = 0; i < 2; ++i) {
                    curr[0] = crds[0] + i * dimensions[0];
                    Node* node = create_node(curr, nodeTag);
                    if (node == 0) {
                        return -1;
                    }
                    nodes.push_back(node);
                }
            }
        }
    }

    // create edge nodes
    if (ndm == 2) {
        if (create_line(nodes[0], nodes[1], nodeTag, 0) < 0) {
            return -1;
        }
        if (top && create_line(nodes[2], nodes[3], nodeTag, 0) < 0) {
            return -1;
        }
        if (create_line(nodes[0], nodes[2], nodeTag, 1) < 0) {
            return -1;
        }
        if (create_line(nodes[1], nodes[3], nodeTag, 1) < 0) {
            return -1;
        }
    } else if (ndm == 3) {
        if (create_line(nodes[0], nodes[1], nodeTag, 0) < 0) {
            return -1;
        }
        if (create_line(nodes[2], nodes[3], nodeTag, 0) < 0) {
            return -1;
        }
        if (create_line(nodes[0], nodes[2], nodeTag, 1) < 0) {
            return -1;
        }
        if (create_line(nodes[1], nodes[3], nodeTag, 1) < 0) {
            return -1;
        }

        if (create_line(nodes[4], nodes[5], nodeTag, 0) < 0) {
            return -1;
        }
        if (create_line(nodes[6], nodes[7], nodeTag, 0) < 0) {
            return -1;
        }
        if (create_line(nodes[4], nodes[6], nodeTag, 1) < 0) {
            return -1;
        }
        if (create_line(nodes[5], nodes[7], nodeTag, 1) < 0) {
            return -1;
        }

        if (create_line(nodes[0], nodes[4], nodeTag, 2) < 0) {
            return -1;
        }
        if (create_line(nodes[1], nodes[5], nodeTag, 2) < 0) {
            return -1;
        }
        if (create_line(nodes[2], nodes[6], nodeTag, 2) < 0) {
            return -1;
        }
        if (create_line(nodes[3], nodes[7], nodeTag, 2) < 0) {
            return -1;
        }
    }

    // create face nodes
    if (ndm == 3) {
        if (create_face(nodes[0], nodes[3], nodeTag, 0, 1) < 0) {
            return -1;
        }
        if (create_face(nodes[0], nodes[5], nodeTag, 0, 2) < 0) {
            return -1;
        }
        if (create_face(nodes[1], nodes[7], nodeTag, 1, 2) < 0) {
            return -1;
        }
        if (create_face(nodes[2], nodes[7], nodeTag, 0, 2) < 0) {
            return -1;
        }
        if (create_face(nodes[0], nodes[6], nodeTag, 1, 2) < 0) {
            return -1;
        }
        if (top &&
            create_face(nodes[4], nodes[7], nodeTag, 0, 1) < 0) {
            return -1;
        }
    }

    ID newndtags(nodeTag - initNodeTag);
    for (int i = 0; i < newndtags.Size(); ++i) {
        newndtags(i) = initNodeTag + i;
    }
    this->setNewNodeTags(newndtags);

    return 0;
}
