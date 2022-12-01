#include "BeamBrick.h"

#include <Domain.h>
#include <Node.h>
#include <elementAPI.h>

#include <cmath>

int OPS_BeamBrick() {
    // get inputs
    int ndm = OPS_GetNDM();

    if (OPS_GetNumRemainingInputArgs() < 2 * ndm + 4) {
        opserr << "WARNING: want tag?, id?, ndf?, x0?, y0?, <z0?>, "
                  "L?, <B?>, H?, size?\n";
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

    // create mesh
    BeamBrick* mesh = new BeamBrick(idata[0], crds, data);
    if (OPS_addMesh(mesh) == false) {
        opserr << "WARNING: failed to add mesh\n";
        return -1;
    }

    mesh->setID(idata[1]);
    mesh->setNdf(idata[2]);
    mesh->setMeshsize(size);

    // set eleArgs
    if (mesh->setEleArgs() < 0) {
        opserr << "WARNING: failed to set element arguments\n";
        return -1;
    }

    // check type
    int type = mesh->getEleType();
    if (type > 0 && type != ELE_TAG_ElasticBeam3d &&
        type != ELE_TAG_ForceBeamColumn3d &&
        type != ELE_TAG_DispBeamColumn3d) {
        opserr
            << "WARNING: element must be elasticBeamColumn, "
               "dispBeamColumn, or forceBeamcolumn for BeamBrick\n";
        return -1;
    }

    // mesh
    if (mesh->mesh() < 0) {
        opserr << "WARNING: failed to create beam brick\n";
        return -1;
    }

    return 0;
}

BeamBrick::BeamBrick(int tag, const std::vector<double>& c,
                     const std::vector<double>& d)
    : Flume(tag, c, d, true), elenodes(), linenodes() {}

BeamBrick::~BeamBrick() {}

int BeamBrick::create_line(Node* nd1, Node* nd2, int& nodeTag,
                           int dir) {
    // check
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        opserr << "no domain - create_line\n";
        return -1;
    }
    if (dir < 0 || dir > 2) {
        opserr << "dir invalid - create_line\n";
        return -1;
    }
    if (nd1 == 0 || nd2 == 0) {
        opserr << "nd1 or nd2 invalid - create_line\n";
        return -1;
    }

    // linenodes
    auto& lnodes = linenodes[std::make_pair(nd1->getTag(), dir)];

    // mesh size
    double size = this->getMeshsize();

    // coordinates
    const Vector& crds1 = nd1->getCrds();
    const Vector& crds2 = nd2->getCrds();
    std::vector<double> curr(crds1.Size());
    for (int i = 0; i < crds1.Size(); ++i) {
        curr[i] = crds1(i);
    }
    if (fabs(crds1[dir] - crds2[dir]) < 1.2 * size) {
        opserr << "nd1 or nd2 too close - create_line\n";
        return -1;
    }
    if (crds1[dir] < crds2[dir]) {
        curr[dir] = crds1[dir] + size;
    } else {
        curr[dir] = crds2[dir] + size;
    }

    // add ele nodes
    elenodes.push_back(nd1->getTag());

    // add line nodes
    lnodes.push_back(nd1->getTag());

    // loop from nd1 to nd2
    while ((crds1[dir] < crds2[dir] &&
            curr[dir] < crds2[dir] - size * 0.1) ||
           (crds2[dir] < crds1[dir] &&
            curr[dir] < crds1[dir] - size * 0.1)) {
        if (create_node(curr, nodeTag) == 0) {
            opserr << "failed to create node - create_line\n";
            return -1;
        }
        elenodes.push_back(nodeTag - 1);
        elenodes.push_back(nodeTag - 1);
        curr[dir] += size;

        lnodes.push_back(nodeTag - 1);
    }

    elenodes.push_back(nd2->getTag());
    lnodes.push_back(nd2->getTag());

    return 0;
}

int BeamBrick::create_face(Node* nd1, Node* nd2, int& nodeTag,
                           int dir1, int dir2) {
    // check
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

    // mesh size
    double size = this->getMeshsize();

    // coordinates
    const Vector& crds1 = nd1->getCrds();
    std::vector<double> curr(crds1.Size());
    for (int i = 0; i < crds1.Size(); ++i) {
        curr[i] = crds1(i);
    }
    curr[dir1] = crds1[dir1] + size;

    // line nodes
    auto& lnodes1 = linenodes[std::make_pair(nd1->getTag(), dir1)];
    auto& lnodes2 = linenodes[std::make_pair(nd1->getTag(), dir2)];
    auto& lnodes3 = linenodes[std::make_pair(lnodes2.back(), dir1)];
    auto& lnodes4 = linenodes[std::make_pair(lnodes1.back(), dir2)];

    if (lnodes1.size() != lnodes3.size()) {
        opserr << "lnodes1 and lnodes3 not match - create_face\n";
        return -1;
    }
    if (lnodes2.size() != lnodes4.size()) {
        opserr << "lnodes2 and lnodes4 not match - create_face\n";
        return -1;
    }

    // map of node tags
    std::vector<std::vector<int>> checkboard(lnodes1.size());
    for (int i = 0; i < (int)lnodes1.size(); ++i) {
        checkboard[i].resize(lnodes2.size());
        checkboard[i][0] = lnodes1[i];
        checkboard[i][lnodes2.size() - 1] = lnodes3[i];
    }
    for (int j = 0; j < (int)lnodes2.size(); ++j) {
        checkboard[0][j] = lnodes2[j];
        checkboard[lnodes1.size() - 1][j] = lnodes4[j];
    }

    // loop through the face
    for (int i = 1; i < (int)lnodes1.size() - 1; ++i) {
        curr[dir2] = crds1[dir2] + size;
        for (int j = 1; j < (int)lnodes2.size() - 1; ++j) {
            // create node at i,j
            if (create_node(curr, nodeTag) == 0) {
                opserr << "failed to create node - create_face\n";
                return -1;
            }
            curr[dir2] += size;
            checkboard[i][j] = nodeTag - 1;

            // add elenodes
            elenodes.push_back(checkboard[i - 1][j]);
            elenodes.push_back(nodeTag - 1);
            elenodes.push_back(checkboard[i][j - 1]);
            elenodes.push_back(nodeTag - 1);

            if (i == (int)lnodes1.size() - 2) {
                elenodes.push_back(nodeTag - 1);
                elenodes.push_back(checkboard[i + 1][j]);
            }
        }
        curr[dir1] += size;

        elenodes.push_back(checkboard[i][lnodes2.size() - 2]);
        elenodes.push_back(checkboard[i][lnodes2.size() - 1]);
    }

    return 0;
}

int BeamBrick::mesh() {
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

    // init elenodes
    elenodes.clear();
    linenodes.clear();

    // create corner nodes
    int nodeTag = nextNodeTag();
    int initNodeTag = nodeTag;
    const auto& crds = this->getCrds();
    const auto& dimensions = this->getDimensions();
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
        if (create_line(nodes[2], nodes[3], nodeTag, 0) < 0) {
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
        if (create_face(nodes[4], nodes[7], nodeTag, 0, 1) < 0) {
            return -1;
        }
    }

    ID newndtags(nodeTag - initNodeTag);
    for (int i = 0; i < newndtags.Size(); ++i) {
        newndtags(i) = initNodeTag + i;
    }
    this->setNewNodeTags(newndtags);

    ID newelenodes(elenodes.size());
    for (int i = 0; i < newelenodes.Size(); ++i) {
        newelenodes(i) = elenodes[i];
    }
    this->setEleNodes(newelenodes);

    // create elements
    if (this->newElements(newelenodes) < 0) {
        opserr << "WARNING: failed to create elements\n";
        return -1;
    }

    return 0;
}
