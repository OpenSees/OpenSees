#include "BeamDisk.h"

#include <Domain.h>
#include <Node.h>
#include <elementAPI.h>

#include <cmath>

double PI = atan(1.0) * 4.0;

int OPS_BeamDisk() {
    // get inputs
    int ndm = OPS_GetNDM();
    if (ndm != 3) {
        opserr << "WARNING: BeamDisk is only for 3D\n";
        return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < 10) {
        opserr << "WARNING: want tag?, id?, ndf?, xc?, yc?, zc?, "
                  "thk? radius? size? dir?\n";
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

    std::vector<double> crds(num);
    if (OPS_GetDoubleInput(&num, &crds[0]) < 0) {
        opserr << "WARNING: failed to get center crds\n";
        return -1;
    }

    num = 4;
    std::vector<double> data(num);
    if (OPS_GetDoubleInput(&num, &data[0]) < 0) {
        opserr << "WARNING: failed to get data\n";
        return -1;
    }
    if (data[0] < 0) {
        opserr << "WARNING: thickness cannot < 0\n";
        return -1;
    }
    if (data[1] < 0) {
        opserr << "WARNING: radius cannot < 0\n";
        return -1;
    }
    if (data[2] < 0) {
        opserr << "WARNING: mesh size cannot < 0\n";
        return -1;
    }
    if (data[3] <= 0 || data[3] > 3) {
        opserr << "WARNING: dir has to be 1, 2, or 3\n";
        return -1;
    }
    --data[3];

    // create mesh
    BeamDisk* mesh = new BeamDisk(idata[0], crds, data);
    if (OPS_addMesh(mesh) == false) {
        opserr << "WARNING: failed to add mesh\n";
        return -1;
    }

    mesh->setID(idata[1]);
    mesh->setNdf(idata[2]);
    mesh->setMeshsize(data[2]);

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
        opserr << "WARNING: element must be elasticBeamColumn, "
                  "dispBeamColumn, or forceBeamcolumn for BeamDisk\n";
        return -1;
    }

    // mesh
    if (mesh->mesh() < 0) {
        opserr << "WARNING: failed to create beam brick\n";
        return -1;
    }

    return 0;
}

BeamDisk::BeamDisk(int tag, const std::vector<double>& center,
                   const std::vector<double>& data)
    : Flume(tag, center, data, true) {}

BeamDisk::~BeamDisk() {}

int BeamDisk::mesh() {
    // dimensions
    const auto& dimensions = this->getDimensions();
    double thk = dimensions[0];
    double radius = dimensions[1];
    double size = dimensions[2];
    int dir = (int)dimensions[3];

    // nlayers
    int nlayers = int(ceil(thk / size)) + 1;
    std::vector<std::vector<Node*>> layerNodes(nlayers);
    std::vector<int> elenodes;

    // directions
    int ndm = OPS_GetNDM();
    int dir1 = dir + 1;
    int dir2 = dir + 2;
    if (dir1 >= ndm) dir1 -= ndm;
    if (dir2 >= ndm) dir2 -= ndm;

    // vertical thickness size
    double vsize = thk / (nlayers - 1);

    // angle size
    double angle = size / radius;
    int nangle = int(ceil(0.5 * PI / angle));
    double asize = 0.5 * PI / nangle;
    int na = 4 * nangle;

    // loop all angles
    int nodeTag = nextNodeTag();
    int initNodeTag = nodeTag;
    std::vector<double> currCrds(ndm);
    const auto& crds = getCrds();
    for (int j = 0; j < na; ++j) {
        Node* prev = 0;
        for (int i = 0; i < nlayers; ++i) {
            // nodetags for this layer
            auto& ndtags = layerNodes[i];

            // loop all angles
            double angle = j * asize;
            currCrds[dir] = crds[dir] + i * vsize;
            currCrds[dir1] = crds[dir1] + cos(angle) * radius;
            currCrds[dir2] = crds[dir2] + sin(angle) * radius;
            Node* node = create_node(currCrds, nodeTag);
            if (node == 0) {
                return -1;
            }
            ndtags.push_back(node);

            // add elements
            if (prev != 0) {
                elenodes.push_back(prev->getTag());
                elenodes.push_back(node->getTag());
            }

            // set prev
            prev = node;
        }
    }

    // connect circles
    for (int i = 0; i < nlayers; ++i) {
        // nodetags for this layer
        auto& ndtags = layerNodes[i];

        // elenodes
        elenodes.push_back(ndtags.front()->getTag());
        for (int j = 1; j < (int)ndtags.size(); ++j) {
            elenodes.push_back(ndtags[j]->getTag());
            elenodes.push_back(ndtags[j]->getTag());
        }
        elenodes.push_back(ndtags.front()->getTag());
    }

    // top and bottom face
    for (int i = 0; i < nlayers; i += nlayers - 1) {
        // nodetags for this layer
        auto& ndtags = layerNodes[i];

        // loop half circle
        std::vector<double> curr(ndm);
        curr[dir] = crds[dir] + i * vsize;
        std::map<int, int> prevList;
        for (int j = 1; j < 2 * nangle; ++j) {
            prevList[2 * nangle + j] =
                ndtags[2 * nangle + j]->getTag();
        }
        for (int j = 1; j < nangle; ++j) {
            // center nodes
            curr[dir1] = ndtags[3 * nangle]->getCrds()[dir1];
            curr[dir2] = ndtags[3 * nangle + j]->getCrds()[dir2];
            int prev1 = nodeTag;
            int prev2 = nodeTag;
            if (create_node(curr, nodeTag) == 0) {
                return -1;
            }
            elenodes.push_back(prevList[3 * nangle]);
            elenodes.push_back(nodeTag - 1);
            prevList[3 * nangle] = nodeTag - 1;

            // side nodes
            for (int k = 1; k < j; ++k) {
                curr[dir1] = ndtags[3 * nangle + k]->getCrds()[dir1];
                if (create_node(curr, nodeTag) == 0) {
                    return -1;
                }
                elenodes.push_back(prev1);
                elenodes.push_back(nodeTag - 1);
                prev1 = nodeTag - 1;

                elenodes.push_back(prevList[3 * nangle + k]);
                elenodes.push_back(nodeTag - 1);
                prevList[3 * nangle + k] = nodeTag - 1;

                curr[dir1] = ndtags[3 * nangle - k]->getCrds()[dir1];
                if (create_node(curr, nodeTag) == 0) {
                    return -1;
                }
                elenodes.push_back(prev2);
                elenodes.push_back(nodeTag - 1);
                prev2 = nodeTag - 1;

                elenodes.push_back(prevList[3 * nangle - k]);
                elenodes.push_back(nodeTag - 1);
                prevList[3 * nangle - k] = nodeTag - 1;
            }

            // ele for ends
            elenodes.push_back(prev1);
            elenodes.push_back(ndtags[3 * nangle + j]->getTag());
            elenodes.push_back(prev2);
            elenodes.push_back(ndtags[3 * nangle - j]->getTag());
        }

        for (int j = 1; j < 2 * nangle; ++j) {
            prevList[j] = ndtags[j]->getTag();
        }
        for (int j = 1; j < nangle + 1; ++j) {
            // center nodes
            curr[dir1] = ndtags[nangle]->getCrds()[dir1];
            curr[dir2] = ndtags[nangle + j]->getCrds()[dir2];
            int prev1 = nodeTag;
            int prev2 = nodeTag;
            if (create_node(curr, nodeTag) == 0) {
                return -1;
            }
            elenodes.push_back(prevList[nangle]);
            elenodes.push_back(nodeTag - 1);
            prevList[nangle] = nodeTag - 1;

            // side nodes
            for (int k = 1; k < j; ++k) {
                curr[dir1] = ndtags[nangle + k]->getCrds()[dir1];
                if (create_node(curr, nodeTag) == 0) {
                    return -1;
                }
                elenodes.push_back(prev1);
                elenodes.push_back(nodeTag - 1);
                prev1 = nodeTag - 1;

                elenodes.push_back(prevList[nangle + k]);
                elenodes.push_back(nodeTag - 1);
                prevList[nangle + k] = nodeTag - 1;

                curr[dir1] = ndtags[nangle - k]->getCrds()[dir1];
                if (create_node(curr, nodeTag) == 0) {
                    return -1;
                }
                elenodes.push_back(prev2);
                elenodes.push_back(nodeTag - 1);
                prev2 = nodeTag - 1;

                elenodes.push_back(prevList[nangle - k]);
                elenodes.push_back(nodeTag - 1);
                prevList[nangle - k] = nodeTag - 1;
            }

            // ele for ends
            elenodes.push_back(prev1);
            elenodes.push_back(ndtags[nangle + j]->getTag());
            elenodes.push_back(prev2);
            elenodes.push_back(ndtags[nangle - j]->getTag());
        }

        for (int j = 1; j < 2 * nangle; ++j) {
            elenodes.push_back(prevList[2 * nangle + j]);
            elenodes.push_back(prevList[2 * nangle - j]);
        }
    }

    // add new nodes
    ID newndtags(nodeTag - initNodeTag);
    for (int i = 0; i < newndtags.Size(); ++i) {
        newndtags(i) = initNodeTag + i;
    }
    this->setNewNodeTags(newndtags);

    // add new elements
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