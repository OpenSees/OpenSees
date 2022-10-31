#include "BeamDisk.h"

#include <Domain.h>
#include <Node.h>
#include <elementAPI.h>

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
