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

#include <elementAPI.h>

#include <cmath>

#include "PFEMElement2DBubble.h"
#include "PFEMElement2DCompressible.h"
#include "PFEMElement2Dmini.h"
#include "PFEMElement3DBubble.h"
#include "ParticleGroup.h"
#ifdef _LINUX
#include <Timer.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include <ConstantSeries.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <Matrix.h>
#include <Node.h>
#include <NodeIter.h>
#include <PFEMContact2D.h>
#include <Pressure_Constraint.h>
#include <SP_Constraint.h>
#include <TetMeshGenerator.h>
#include <TriangleMeshGenerator.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <iostream>

static BackgroundMesh bgmesh;

BackgroundMesh& OPS_getBgMesh() { return bgmesh; }

// OPS_BgMesh
int OPS_BgMesh() {
    int ndm = OPS_GetNDM();

    // check input
    if (OPS_GetNumRemainingInputArgs() < 2 * ndm + 1) {
        opserr << "WARNING: basicsize? lower? upper? <-tol tol? "
                  "-wave wavefilename? numl? locs? -numsub numsub? "
                  "-dispon? "
                  "-structure sid? ?numnodes? structuralNodes? "
                  "-alphaS sid? alphaS?"
                  "-contact kdoverAd? thk? mu? beta? Dc? alpha? E? rho?>";
        return -1;
    }

    // get basicsize
    double size;
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: need basic size\n";
        return -1;
    }
    int num = 1;
    if (OPS_GetDoubleInput(&num, &size) < 0) {
        opserr << "WARNING: failed to get basic size\n";
        return -1;
    }
    if (size <= 0) {
        opserr << "WARNING: basic size <= 0\n";
        return -1;
    }
    bgmesh.setBasicSize(size);

    // get range
    VDouble lower(ndm), upper(ndm);
    double* ptr = &lower[0];
    if (OPS_GetDoubleInput(&ndm, ptr) < 0) {
        opserr << "WARNING: failed to get min\n";
        return -1;
    }
    ptr = &upper[0];
    if (OPS_GetDoubleInput(&ndm, ptr) < 0) {
        opserr << "WARNING: failed to get max\n";
        return -1;
    }
    bgmesh.setRange(lower, upper);
    bgmesh.setDispOn(false);

    // get tolerance
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();

        if (strcmp(opt, "-tol") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: need tol\n";
                return -1;
            }
            double tol;
            if (OPS_GetDoubleInput(&num, &tol) < 0) {
                opserr << "WARNING: failed to read tolerance\n";
                return -1;
            }
            bgmesh.setTol(tol);

        } else if (strcmp(opt, "-wave") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 2) {
                opserr << "WARNING: need wavefilename and numl\n";
                return -1;
            }

            const char* wavefilename = OPS_GetString();
            if (bgmesh.setFile(wavefilename) < 0) {
                return -1;
            }

            int numl;
            num = 1;
            if (OPS_GetIntInput(&num, &numl) < 0) {
                opserr << "WARNING: failed to read numl\n";
                return -1;
            }
            if (numl > 0) {
                num = numl * ndm;
                if (OPS_GetNumRemainingInputArgs() < num) {
                    opserr
                        << "WARNING: insufficient number of locations\n";
                    return -1;
                }

                VDouble locs(num);
                if (OPS_GetDoubleInput(&num, &locs[0]) < 0) {
                    opserr << "WARNING: failed to read wave recording "
                              "locations\n";
                    return -1;
                }
                bgmesh.setLocs(locs);
            }

        } else if (strcmp(opt, "-numsub") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: need numsub\n";
                return -1;
            }
            int numsub;
            num = 1;
            if (OPS_GetIntInput(&num, &numsub) < 0) {
                opserr << "WARNING: failed to read numsub\n";
                return -1;
            }
            if (numsub > 0) {
                bgmesh.setNumSub(numsub);
            }
        } else if (strcmp(opt, "-structure") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 2) {
                opserr << "WARNING: need sid and numnodes\n";
                return -1;
            }
            int sid, numnodes;
            num = 1;
            if (OPS_GetIntInput(&num, &sid) < 0) {
                opserr << "WARNING: failed to read sid\n";
                return -1;
            }
            if (OPS_GetIntInput(&num, &numnodes) < 0) {
                opserr << "WARNING: failed to read numnodes\n";
                return -1;
            }
            if (OPS_GetNumRemainingInputArgs() < numnodes) {
                opserr << "WARNING: insufficient number of structural "
                          "nodes\n";
                return -1;
            }
            if (numnodes > 0) {
                VInt snodes(numnodes);
                if (OPS_GetIntInput(&numnodes, &snodes[0]) < 0) {
                    opserr << "WARNING: failed to read structural nodes\n";
                    return -1;
                }
                bgmesh.addStructuralNodes(snodes, sid);
            }

        } else if (strcmp(opt, "-contact") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 8) {
                opserr << "WARNING: need kdoverAd, thk, mu, beta, Dc, "
                          "alpha, E, rho\n";
                return -1;
            }
            num = 8;
            VDouble data(num);
            if (OPS_GetDoubleInput(&num, &data[0]) < 0) {
                opserr << "WARNING: failed to get kdoverAd, thk, mu, "
                          "beta, Dc, alpha, E, rho\n";
                return -1;
            }
            bgmesh.setContactData(data);

        } else if (strcmp(opt, "-alphaS") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 2) {
                opserr << "WARNING: need sid alphaS\n";
                return -1;
            }
            num = 1;
            int sid;
            if (OPS_GetIntInput(&num, &sid) < 0) {
                opserr << "WARNING: failed to get sid\n";
                return -1;
            }

            double alpha = 0.0;
            if (OPS_GetDoubleInput(&num, &alpha) < 0) {
                opserr << "WARNING: failed to get alphaS\n";
                return -1;
            }
            if (alpha <= 0) {
                alpha = 0.0;
            } else if (alpha > 1) {
                alpha = 1.0;
            }
            bgmesh.setAlphaS(sid, alpha);

        } else if (strcmp(opt, "-dispOn") == 0) {
            bgmesh.setDispOn(true);
        }
    }

    // bg mesh
    if (bgmesh.remesh(true) < 0) {
        opserr << "WARNING: failed to create background mesh\n";
        return -1;
    }

    return 0;
}

BackgroundMesh::BackgroundMesh()
    : lower(),
      upper(),
      bcells(),
      bnodes(),
      tol(1e-10),
      bsize(-1.0),
      numave(2),
      numsub(4),
      recorders(),
      locs(),
      currentTime(0.0),
      theFile(),
      structuralNodes(),
      contactData(8),
      contactEles(),
      dispon(true),
      alphaS() {}

BackgroundMesh::~BackgroundMesh() {
    for (int i = 0; i < (int)recorders.size(); ++i) {
        if (recorders[i] != 0) {
            delete recorders[i];
        }
    }
    recorders.clear();
}

void BackgroundMesh::addRecorder(Recorder* recorder) {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;

    recorder->setDomain(*domain);
    recorders.push_back(recorder);
}

void BackgroundMesh::setRange(const VDouble& l, const VDouble& u) {
    nearIndex(l, lower);
    nearIndex(u, upper);
}

int BackgroundMesh::setFile(const char* name) {
    theFile.open(name, std::ios::trunc | std::ios::out);
    if (theFile.fail()) {
        opserr << "WARNING: Failed to open file " << name << "\n";
        return -1;
    }

    int precision = 12;
    theFile.precision(precision);
    theFile << std::scientific;

    return 0;
}

void BackgroundMesh::setDispOn(bool on) {
    dispon = on;
    PFEMElement2DBubble::dispon = on;
    PFEMElement3DBubble::dispon = on;
    PFEMElement2DCompressible::dispon = on;
    PFEMElement2Dmini::dispon = on;
}

void BackgroundMesh::addStructuralNodes(VInt& snodes, int sid) {
    VInt& curr = structuralNodes[sid];
    for (int i = 0; i < (int)snodes.size(); ++i) {
        curr.push_back(snodes[i]);
    }
}

void BackgroundMesh::getIndex(const VDouble& crds, double incr,
                              VInt& index) const {
    index.resize(crds.size());
    for (int i = 0; i < (int)crds.size(); ++i) {
        double crd = crds[i] / bsize + incr;
        index[i] = (int)floor(crd);
    }
}

void BackgroundMesh::lowerIndex(const VDouble& crds, VInt& index) const {
    getIndex(crds, 0.0, index);
}

void BackgroundMesh::upperIndex(const VDouble& crds, VInt& index) const {
    getIndex(crds, 1.0, index);
}

void BackgroundMesh::nearIndex(const VDouble& crds, VInt& index) const {
    getIndex(crds, 0.5, index);
}

void BackgroundMesh::getCrds(const VInt& index, VDouble& crds) const {
    crds.resize(index.size(), 0.0);
    for (int i = 0; i < (int)crds.size(); ++i) {
        crds[i] = index[i] * bsize;
    }
}

// get corners from index
// count num+1 in each direction
// 2   3
// -----
// |   |
// |   |
// -----
// 0   1
void BackgroundMesh::getCorners(const VInt& index, int num,
                                VVInt& indices) const {
    int ndm = OPS_GetNDM();
    int counter = 0;

    if (ndm == 2) {
        indices.resize((num + 1) * (num + 1));
        for (int j = index[1]; j <= index[1] + num; j += 1) {
            for (int i = index[0]; i <= index[0] + num; i += 1) {
                indices[counter].resize(ndm);
                indices[counter][0] = i;
                indices[counter][1] = j;
                ++counter;
            }
        }
    } else if (ndm == 3) {
        indices.resize((num + 1) * (num + 1) * (num + 1));
        for (int k = index[2]; k <= index[2] + num; k += 1) {
            for (int j = index[1]; j <= index[1] + num; j += 1) {
                for (int i = index[0]; i <= index[0] + num; i += 1) {
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
void BackgroundMesh::gatherParticles(const VInt& minind,
                                     const VInt& maxind, VParticle& pts,
                                     bool checkfsi) {
    int ndm = OPS_GetNDM();
    pts.clear();
    VInt index(ndm);
    if (ndm == 2) {
        for (int i = minind[0]; i < maxind[0]; ++i) {
            index[0] = i;
            for (int j = minind[1]; j < maxind[1]; ++j) {
                index[1] = j;
                std::map<VInt, BCell>::iterator it = bcells.find(index);
                if (it != bcells.end()) {
                    BCell& cell = it->second;
                    if (checkfsi && cell.getType() == BACKGROUND_FLUID) {
                        continue;
                    }
                    pts.insert(pts.end(), cell.getPts().begin(),
                               cell.getPts().end());
                }
            }
        }
    } else if (ndm == 3) {
        for (int i = minind[0]; i < maxind[0]; ++i) {
            index[0] = i;
            for (int j = minind[1]; j < maxind[1]; ++j) {
                index[1] = j;
                for (int k = minind[2]; k < maxind[2]; ++k) {
                    index[2] = k;
                    std::map<VInt, BCell>::iterator it =
                        bcells.find(index);
                    if (it != bcells.end()) {
                        BCell& cell = it->second;
                        if (checkfsi &&
                            cell.getType() == BACKGROUND_FLUID) {
                            continue;
                        }
                        pts.insert(pts.end(), cell.getPts().begin(),
                                   cell.getPts().end());
                    }
                }
            }
        }
    }
}

double BackgroundMesh::QuinticKernel(double q, double h, int ndm) {
    static double pi = 3.141592653589793;
    if (q < 0 || q > 2) return 0.0;
    double aD = 0.0;
    if (ndm == 2) {
        aD = 7.0 / (4 * pi * h * h);
    } else if (ndm == 3) {
        aD = 7.0 / (8 * pi * h * h * h);
    }
    double a = 1.0 - q / 2.0;

    return aD * a * a * a * a * (2 * q + 1);
}

int BackgroundMesh::preNForTri(double x1, double y1, double x2, double y2,
                               double x3, double y3, VDouble& coeff) {
    coeff.resize(9, 0.0);

    coeff[0] = x2 * y3 - x3 * y2;
    coeff[1] = x3 * y1 - x1 * y3;
    coeff[2] = x1 * y2 - x2 * y1;

    coeff[3] = y2 - y3;
    coeff[4] = y3 - y1;
    coeff[5] = y1 - y2;

    coeff[6] = x3 - x2;
    coeff[7] = x1 - x3;
    coeff[8] = x2 - x1;

    double A = coeff[0] + coeff[1] + coeff[2];

    if (A < 0 || fabs(A) < 1e-15) {
        // opserr << "A <= 0\n";
        return -1;
    }

    for (int i = 0; i < (int)coeff.size(); ++i) {
        coeff[i] /= A;
    }

    return 0;
}

int BackgroundMesh::preNForTet(const VDouble& crds1, const VDouble& crds2,
                               const VDouble& crds3, const VDouble& crds4,
                               VVDouble& coeff) {
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
    Matrix Jmat(4, 4), Jfact(4, 4);
    Jmat(0, 0) = 1.0;
    Jmat(0, 1) = 1.0;
    Jmat(0, 2) = 1.0;
    Jmat(0, 3) = 1.0;
    for (int j = 0; j < ndm; ++j) {
        Jmat(j + 1, 0) = crds1[j];
        Jmat(j + 1, 1) = crds2[j];
        Jmat(j + 1, 2) = crds3[j];
        Jmat(j + 1, 3) = crds4[j];
    }

    PFEMElement3DBubble::cofactor(Jmat, Jfact);

    coeff.resize(Jfact.noRows());
    double vol = 0.0;
    for (int i = 0; i < Jfact.noRows(); ++i) {
        coeff[i].resize(Jfact.noCols());
        for (int j = 0; j < Jfact.noCols(); ++j) {
            coeff[i][j] = Jfact(i, j);
        }
        vol += coeff[i][0];
    }

    if (vol < 0 || fabs(vol) < 1e-14) {
        // opserr<<"vol "<<vol<<" <= 0\n";
        return -1;
    }

    for (int i = 0; i < Jfact.noRows(); ++i) {
        for (int j = 0; j < Jfact.noCols(); ++j) {
            coeff[i][j] /= vol;
        }
    }

    return 0;
}

void BackgroundMesh::getNForTri(const VDouble& coeff, double x, double y,
                                VDouble& N) {
    N.resize(3, 0.0);

    for (int i = 0; i < (int)N.size(); ++i) {
        double val = coeff[i] + coeff[i + 3] * x + coeff[i + 6] * y;
        if (fabs(val) < tol) {
            // make sure it's in the tri
            N[i] = tol;
        }
        N[i] = val;
    }
}

void BackgroundMesh::getNForTet(const VVDouble& coeff, const VDouble& crds,
                                VDouble& N) {
    if (crds.size() != 3) {
        return;
    }
    if (coeff.size() != 4) {
        return;
    }

    N.resize(4, 0.0);
    VDouble col(4);
    col[0] = 1.0;
    for (int i = 0; i < (int)crds.size(); ++i) {
        col[i + 1] = crds[i];
    }

    for (int i = 0; i < (int)coeff.size(); ++i) {
        if (coeff[i].size() != 4) {
            return;
        }
        N[i] = dotVDouble(coeff[i], col);
        if (fabs(N[i]) < tol) {
            N[i] = tol;
        }
    }
}

void BackgroundMesh::getNForRect(double x0, double y0, double hx,
                                 double hy, double x, double y,
                                 VDouble& N) {
    // compute local coordinate of the particle
    double xl = (x - x0) / hx;
    double yl = (y - y0) / hy;

    // map to [-1, 1]
    xl = xl * 2 - 1;
    yl = yl * 2 - 1;

    // shape function
    N.resize(4);
    N[0] = (1 - xl) * (1 - yl) / 4.0;
    N[1] = (1 + xl) * (1 - yl) / 4.0;
    N[2] = (1 + xl) * (1 + yl) / 4.0;
    N[3] = (1 - xl) * (1 + yl) / 4.0;
}

void BackgroundMesh::getNForRect(double x0, double y0, double z0,
                                 double hx, double hy, double hz, double x,
                                 double y, double z, VDouble& N) {
    // compute local coordinate of the particle
    double xl = (x - x0) / hx;
    double yl = (y - y0) / hy;
    double zl = (z - z0) / hz;

    // map to [-1, 1]
    xl = xl * 2 - 1;
    yl = yl * 2 - 1;
    zl = zl * 2 - 1;

    // shape function
    N.resize(8);
    N[0] = (1 - xl) * (1 - yl) * (1 - zl) / 8.0;
    N[1] = (1 + xl) * (1 - yl) * (1 - zl) / 8.0;
    N[2] = (1 + xl) * (1 + yl) * (1 - zl) / 8.0;
    N[3] = (1 - xl) * (1 + yl) * (1 - zl) / 8.0;
    N[4] = (1 - xl) * (1 - yl) * (1 + zl) / 8.0;
    N[5] = (1 + xl) * (1 - yl) * (1 + zl) / 8.0;
    N[6] = (1 + xl) * (1 + yl) * (1 + zl) / 8.0;
    N[7] = (1 - xl) * (1 + yl) * (1 + zl) / 8.0;
}

void BackgroundMesh::clearAll() {
    clearBackground();
    lower.clear();
    upper.clear();
    bcells.clear();
    bnodes.clear();
    tol = 1e-10;
    bsize = -1.0;
    numave = 2;
    numsub = 4;

    for (int i = 0; i < (int)recorders.size(); ++i) {
        if (recorders[i] != 0) {
            delete recorders[i];
        }
    }
    recorders.clear();
    locs.clear();
    currentTime = 0.0;
    theFile.close();
    structuralNodes.clear();
    for (int i = 0; i < (int)contactData.size(); ++i) {
        contactData[i] = 0.0;
    }
    contactEles.clear();
    dispon = true;
    alphaS.clear();
}

int BackgroundMesh::clearBackground() {
    // remove elements
    clearGridEles();

    // remove cells
    clearGrid();

    return 0;
}

void BackgroundMesh::clearGridEles() {
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

void BackgroundMesh::clearGrid() {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return;

    // remove cells
    for (std::map<VInt, BNode>::iterator it = bnodes.begin();
         it != bnodes.end(); ++it) {
        BNode& bnode = it->second;
        const VInt& tags = bnode.getTags();
        int type = bnode.getType();

        for (int i = 0; i < (int)tags.size(); ++i) {
            if (type == BACKGROUND_FLUID) {
                // remove node
                Node* nd = domain->removeNode(tags[i]);
                if (nd != 0) {
                    delete nd;
                }

                // remove pc
                Pressure_Constraint* pc =
                    domain->removePressure_Constraint(tags[i]);
                if (pc != 0) {
                    delete pc;
                }
            }
        }
    }

    bnodes.clear();
    bcells.clear();
}

bool BackgroundMesh::inEle(const VDouble& N) {
    // out
    for (int j = 0; j < (int)N.size(); ++j) {
        if (N[j] < 0) {
            // j+1
            return false;
        }
    }

    return true;
}

int BackgroundMesh::solveLine(const VDouble& p1, const VDouble& dir,
                              int dim, double crd, double& k) {
    // check
    if (p1.size() != dir.size()) {
        opserr
            << "WARNING: sizes are not compatible -- BgMesh::solveLine\n";
        return -1;
    }
    if (dim < 0 || dim >= (int)dir.size()) {
        opserr << "WARNING: dim is out of range -- BgMesh::solveLine\n";
        return -1;
    }

    // solve k
    if (dir[dim] == 0.0) {
        k = -1.0;
    } else {
        k = (crd - p1[dim]) / dir[dim];
    }

    return 0;
}

int BackgroundMesh::remesh(bool init) {
    // clear and check
    if (bsize <= 0.0) {
        opserr << "WARNING: basic mesh size has not been set -- "
                  "BgMesh::addParticles\n";
        return -1;
    }

#ifdef _LINUX
    Timer timer;
    timer.start();
#endif

    // move particles
    if (moveParticles() < 0) {
        opserr << "WARNING: failed to move particles\n";
        return -1;
    }

#ifdef _LINUX
    timer.pause();
    opserr << "time for move particles = " << timer.getReal() << "\n";
    timer.start();
#endif

    // clear background
    clearBackground();

    // add structure
    if (addStructure() < 0) {
        opserr << "WARNING: failed to add structure\n";
        return -1;
    }

#ifdef _LINUX
    timer.pause();
    opserr << "time for add structure = " << timer.getReal() << "\n";
    timer.start();
#endif
    // add particles
    if (addParticles() < 0) {
        opserr << "WARNING: failed to add particles\n";
        return -1;
    }

#ifdef _LINUX
    timer.pause();
    opserr << "time for add particles = " << timer.getReal() << "\n";
    timer.start();
#endif

    // move particles in fixed cells
    if (moveFixedParticles()) {
        opserr << "WARNING: failed to move particles in fixed cells";
        return -1;
    }

#ifdef _LINUX
    timer.pause();
    opserr << "time for moving fixed particles = " << timer.getReal()
           << "\n";
    timer.start();
#endif

    // create grid nodes
    if (gridNodes() < 0) {
        opserr << "WARNING: failed to create grid nodes\n";
        return -1;
    }

#ifdef _LINUX
    timer.pause();
    opserr << "time for grid nodes = " << timer.getReal() << "\n";
    timer.start();
#endif

    // create grid elements
    if (gridFluid() < 0) {
        opserr << "WARNING: failed to create fluid elements\n";
        return -1;
    }

#ifdef _LINUX
    timer.pause();
    opserr << "time for fluid eles = " << timer.getReal() << "\n";
    timer.start();
#endif

    // create FSI elements
    if (gridFSInoDT() < 0) {
        opserr << "WARNING: failed to create FSI elements\n";
        return -1;
    }

#ifdef _LINUX
    timer.pause();
    opserr << "time for fsi eles = " << timer.getReal() << "\n";
    timer.start();
#endif

    if (record(init) < 0) {
        opserr << "WARNING: failed to record\n";
        return -1;
    }

#ifdef _LINUX
    timer.pause();
    opserr << "time for recording = " << timer.getReal() << "\n";
    timer.start();
#endif

    return 0;
}

// for each structural node
// get current states
// find nearest bnode
// if bnode is empty, clear it
// add structure node to bnodes
// get min and max of FSI area
// set cells close to the structures as STRUCTURE
// set grids close to the structures as EMPTY unless it's STRUCTURE
// set FSI area bnodes and cells
int BackgroundMesh::addStructure() {
    // get domain
    int ndm = OPS_GetNDM();
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    // add all structural nodes to the background
    int ndtag = Mesh::nextNodeTag();
    std::set<int> allnodes;
    for (std::map<int, VInt>::iterator it = structuralNodes.begin();
         it != structuralNodes.end(); ++it) {
        int sid = it->first;
        const VInt& snodes = it->second;

        for (int k = 0; k < (int)snodes.size(); ++k) {
            // check if already there
            std::pair<std::set<int>::iterator, bool> res =
                allnodes.insert(snodes[k]);
            if (res.second == false) {
                continue;
            }

            // get node
            Node* nd = domain->getNode(snodes[k]);
            if (nd == 0) continue;

            // nodal data
            const Vector& crds = nd->getCrds();
            const Vector& disp = nd->getTrialDisp();
            const Vector& vel = nd->getTrialVel();
            const Vector& accel = nd->getTrialAccel();

            if (crds.Size() != ndm || disp.Size() < ndm) {
                continue;
            }

            // create pressure constraint
            Pressure_Constraint* pc =
                domain->getPressure_Constraint(nd->getTag());
            if (pc != 0) {
                pc->setDomain(domain);
            } else {
                // create pressure node
                Node* pnode = 0;
                if (ndm == 2) {
                    pnode = new Node(ndtag++, 1, crds[0], crds[1]);
                } else if (ndm == 3) {
                    pnode =
                        new Node(ndtag++, 1, crds[0], crds[1], crds[2]);
                }
                if (pnode == 0) {
                    opserr << "WARNING: run out of memory -- "
                              "BgMesh::gridNodes\n";
                    return -1;
                }
                if (domain->addNode(pnode) == false) {
                    opserr << "WARNING: failed to add node to domain -- "
                              "BgMesh::gridNodes\n";
                    delete pnode;
                    return -1;
                }

                pc =
                    new Pressure_Constraint(nd->getTag(), pnode->getTag());
                if (pc == 0) {
                    opserr << "WARNING: no enough memory for "
                              "Pressure_Constraint\n";
                    return -1;
                }
                if (domain->addPressure_Constraint(pc) == false) {
                    opserr << "WARNING: failed to add PC to domain -- "
                              "BgMesh::gridNodes\n";
                    delete pc;
                    return -1;
                }
            }

            // nodal pressures
            double pressure = pc->getPressure();
            double pdot = pc->getPdot();

            // nodal crds
            VDouble crdsn(ndm), vn(ndm), dvn(ndm);
            for (int i = 0; i < ndm; ++i) {
                crdsn[i] = crds(i) + disp(i);
                vn[i] = vel(i);
                dvn[i] = accel(i);
            }

            // near index
            VInt index;
            nearIndex(crdsn, index);

            // add structural node to the bnode
            BNode& bnode = bnodes[index];
            bnode.addNode(nd->getTag(), crdsn, vn, dvn, pressure, pdot,
                          BACKGROUND_STRUCTURE, sid);

            // set fixed  bnodes
            VInt ind = index;
            ind -= 1;
            VVInt indices;
            getCorners(ind, 2, indices);
            for (int i = 0; i < (int)indices.size(); ++i) {
                BNode& bnd = bnodes[indices[i]];
                if (bnd.size() == 0) {
                    bnd.setType(BACKGROUND_FIXED);
                }
            }

            // set STRUCTURE cells
            getCorners(ind, 1, indices);
            for (int i = 0; i < (int)indices.size(); ++i) {
                BCell& bcell = bcells[indices[i]];
                bcell.setType(BACKGROUND_STRUCTURE);

                // set corners
                if (bcell.getNodes().empty()) {
                    VVInt corners;
                    getCorners(indices[i], 1, corners);

                    for (int j = 0; j < (int)corners.size(); ++j) {
                        BNode& bnd = bnodes[corners[j]];
                        bcell.addNode(&bnd, corners[j]);
                    }
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
int BackgroundMesh::addParticles() {
    // for all particles
    TaggedObjectIter& meshes = OPS_getAllMesh();
    Mesh* mesh = 0;
    while ((mesh = dynamic_cast<Mesh*>(meshes())) != 0) {
        ParticleGroup* group = dynamic_cast<ParticleGroup*>(mesh);
        if (group == 0) {
            continue;
        }

        // check group tag
        if (group->getTag() == contact_tag) {
            opserr << "WARNING: the particle group tag " << contact_tag;
            opserr << " is reserved for internal use. Please select a "
                      "different one\n";
            return -1;
        }

        // remove particles
        VInt rm(group->numParticles(), 0);

        // for all particles
        for (int j = 0; j < group->numParticles(); j++) {
            // get particle
            Particle* p = group->getParticle(j);
            if (p == 0) continue;

            // get particle coordinates
            const VDouble& crds = p->getCrds();

            // get index
            VInt index;
            lowerIndex(crds, index);

            // if out of range
            for (int i = 0; i < (int)index.size(); ++i) {
                if (index[i] < lower[i] || index[i] >= upper[i]) {
                    rm[j] = 1;
                    break;
                }
            }
            if (rm[j] == 1) continue;

            // get bcell
            BCell& bcell = bcells[index];

            // add particles
            bcell.add(p);

            // if initial check structure cell
            if (bcell.getType() == BACKGROUND_STRUCTURE) {
                continue;
            }

            // add bnodes of the cell
            if (bcell.getNodes().empty()) {
                // get corners
                VVInt indices;
                getCorners(index, 1, indices);

                // set corners
                for (int i = 0; i < (int)indices.size(); ++i) {
                    BNode& bnode = bnodes[indices[i]];
                    if (bnode.size() == 0) {
                        bnode.setType(BACKGROUND_FLUID);
                    }
                    bcell.addNode(&bnode, indices[i]);
                }
            }
        }

        // remove out of range particles
        group->removeParticles(rm);
    }

    return 0;
}

int BackgroundMesh::gridNodes() {
    // get domain
    int ndm = OPS_GetNDM();
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    // vector of iterators
    std::vector<std::map<VInt, BNode>::iterator> iters;
    iters.reserve(bnodes.size());
    for (std::map<VInt, BNode>::iterator it = bnodes.begin();
         it != bnodes.end(); ++it) {
        iters.push_back(it);
    }

    // vector of new objects
    int ndtag = Mesh::nextNodeTag();
    std::vector<Node*> newnodes(iters.size(), 0),
        newpnodes(iters.size(), 0);
    std::vector<Pressure_Constraint*> newpcs(iters.size(), 0);

    int res = 0;

#pragma omp parallel for
    for (int j = 0; j < (int)iters.size(); ++j) {
        // get iterator
        std::map<VInt, BNode>::iterator it = iters[j];

        // get cell
        const VInt& index = it->first;
        BNode& bnode = it->second;
        if (bnode.getType() == BACKGROUND_FIXED) {
            continue;
        }

        // coordinates
        VDouble crds;
        getCrds(index, crds);

        // get particles
        VParticle pts;
        VInt minind = index;
        VInt maxind = index;
        minind -= numave;
        maxind += numave;
        gatherParticles(minind, maxind, pts);

        // get information
        double wt = 0.0, pre = 0.0, pdot = 0.0;
        VDouble vel(ndm), accel(ndm);
        for (int i = 0; i < (int)pts.size(); ++i) {
            // get particle
            if (pts[i] == 0) {
                continue;
            }

            // particle coordinates
            const VDouble& pcrds = pts[i]->getCrds();

            // distance from particle to current location
            VDouble dist = pcrds;
            dist -= crds;
            double q = normVDouble(dist) / (bsize * numave);

            // weight for the particle
            double w = QuinticKernel(q, bsize, ndm);

            // check velocity
            const VDouble& pvel = pts[i]->getVel();

            // add pressure
            pre += pts[i]->getPressure() * w;
            pdot += pts[i]->getPdot() * w;

            // add velocity
            for (int k = 0; k < ndm; k++) {
                vel[k] += w * pvel[k];
            }

            // add acceleration
            const VDouble& paccel = pts[i]->getAccel();
            for (int k = 0; k < ndm; k++) {
                accel[k] += w * paccel[k];
            }

            wt += w;
        }

        // get nodal states
        if (wt > 0) {
            pre /= wt;
            pdot /= wt;
            vel /= wt;
            accel /= wt;
        }

        // update pressure for structural nodes
        if (bnode.getType() == BACKGROUND_STRUCTURE) {
            auto& tags = bnode.getTags();
            for (int i = 0; i < (int)bnode.size(); ++i) {
                Pressure_Constraint* pc =
                    domain->getPressure_Constraint(tags[i]);
                if (pc == 0) {
                    opserr << "WARNING: structural node " << tags[i];
                    opserr << " has not pc associated\n";
                    continue;
                }
                pc->setPressure(pre);
                pc->setPdot(pdot);
            }

            continue;
        }

        // create node
        Node* node = 0;
        if (ndm == 2) {
            node = new Node(ndtag + 2 * j, ndm, crds[0], crds[1]);
        } else if (ndm == 3) {
            node = new Node(ndtag + 2 * j, ndm, crds[0], crds[1], crds[2]);
        }
        if (node == 0) {
            opserr << "WARNING: run out of memory -- BgMesh::gridNodes\n";
            res = -1;
            continue;
        }

        if (wt > 0) {
            Vector vec;
            toVector(vel, vec);
            node->setTrialVel(vec);
            toVector(accel, vec);
            node->setTrialAccel(vec);
            node->commitState();
        }

        // add to newnodes
        newnodes[j] = node;

        // set the bnode
        bnode.addNode(node->getTag(), crds, vel, accel, pre, pdot,
                      BACKGROUND_FLUID);

        // set pressure
        Pressure_Constraint* thePC =
            domain->getPressure_Constraint(node->getTag());
        Node* pnode = 0;
        if (thePC != 0) {
            thePC->setDomain(domain);
            pnode = thePC->getPressureNode();
            if (pnode == 0) {
                opserr << "WARNING: pressure does not exist -- "
                          "BgMesh::gridNodes\n";
                res = -1;
                continue;
            }
        } else {
            // create pressure node
            if (ndm == 2) {
                pnode = new Node(ndtag + 2 * j + 1, 1, crds[0], crds[1]);
            } else if (ndm == 3) {
                pnode = new Node(ndtag + 2 * j + 1, 1, crds[0], crds[1],
                                 crds[2]);
            }
            if (pnode == 0) {
                opserr
                    << "WARNING: run out of memory -- BgMesh::gridNodes\n";
                res = -1;
                continue;
            }
            newpnodes[j] = pnode;

            thePC =
                new Pressure_Constraint(node->getTag(), pnode->getTag());
            if (thePC == 0) {
                opserr << "WARNING: no enough memory for "
                          "Pressure_Constraint\n";
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
    for (int i = 0; i < (int)newnodes.size(); ++i) {
        if (newnodes[i] == 0) continue;

        // add to domain
        if (domain->addNode(newnodes[i]) == false) {
            opserr << "WARNING: failed to add node to domain -- "
                      "BgMesh::gridNodes\n";
            delete newnodes[i];
            return -1;
        }
    }
    for (int i = 0; i < (int)newpnodes.size(); ++i) {
        if (newpnodes[i] == 0) continue;

        // add to domain
        if (domain->addNode(newpnodes[i]) == false) {
            opserr << "WARNING: failed to add node to domain -- "
                      "BgMesh::gridNodes\n";
            delete newpnodes[i];
            return -1;
        }
    }
    for (int i = 0; i < (int)newpcs.size(); ++i) {
        if (newpcs[i] == 0) continue;

        // add to domain
        if (domain->addPressure_Constraint(newpcs[i]) == false) {
            opserr << "WARNING: failed to add PC to domain -- "
                      "BgMesh::gridNodes\n";
            delete newpcs[i];
            return -1;
        }
    }

    return 0;
}

int BackgroundMesh::moveFixedParticles() {
    int ndm = OPS_GetNDM();

    // check each cell
    for (std::map<VInt, BCell>::iterator it = bcells.begin();
         it != bcells.end(); ++it) {
        // get cell
        const VInt& index = it->first;
        BCell& cell = it->second;

        // empty cell
        if (cell.getPts().empty()) {
            continue;
        }

        // check if BACKGROUND_STRUCTURE cell
        if (cell.getType() != BACKGROUND_STRUCTURE) {
            continue;
        }

        // check neighbor cells
        VInt ind = index;
        VVInt indices;
        ind -= 1;
        getCorners(ind, 2, indices);

        // give each cell a score
        VInt scores(indices.size());
        for (int i = 0; i < (int)indices.size(); ++i) {
            std::map<VInt, BCell>::iterator cellit =
                bcells.find(indices[i]);
            if (cellit == bcells.end()) {
                // empty cell
                scores[i] = 2;
                continue;
            }
            if (cellit->second.getType() == BACKGROUND_STRUCTURE) {
                scores[i] = -1;
            } else {
                scores[i] = 1;
            }
        }
        VVInt cellmap;
        if (ndm == 2) {
            cellmap.resize(9);
            for (int i = 0; i < (int)cellmap.size(); ++i) {
                cellmap[i].resize(3);
                cellmap[i][0] = i;
            }
            cellmap[0][1] = 1;
            cellmap[0][2] = 3;
            cellmap[1][1] = 0;
            cellmap[1][2] = 2;
            cellmap[2][1] = 1;
            cellmap[2][2] = 5;
            cellmap[3][1] = 0;
            cellmap[3][2] = 6;
            cellmap[4][1] = 4;
            cellmap[4][2] = 4;
            cellmap[5][1] = 2;
            cellmap[5][2] = 8;
            cellmap[6][1] = 3;
            cellmap[6][2] = 7;
            cellmap[7][1] = 6;
            cellmap[7][2] = 8;
            cellmap[8][1] = 5;
            cellmap[8][2] = 7;

        } else if (ndm == 3) {
            cellmap.resize(27);

            int cm0[] = {0, 1, 3, 9};
            cellmap[0].assign(cm0, cm0 + 4);

            int cm1[] = {0, 1, 2, 10};
            cellmap[1].assign(cm1, cm1 + 4);

            int cm2[] = {1, 2, 5, 11};
            cellmap[2].assign(cm2, cm2 + 4);

            int cm3[] = {0, 3, 6, 12};
            cellmap[3].assign(cm3, cm3 + 4);

            int cm4[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
            cellmap[4].assign(cm4, cm4 + 9);

            int cm5[] = {2, 5, 8, 14};
            cellmap[5].assign(cm5, cm5 + 4);

            int cm6[] = {3, 6, 7, 15};
            cellmap[6].assign(cm6, cm6 + 4);

            int cm7[] = {6, 7, 8, 16};
            cellmap[7].assign(cm7, cm7 + 4);

            int cm8[] = {5, 8, 7, 17};
            cellmap[8].assign(cm8, cm8 + 4);

            int cm9[] = {9, 10, 12, 0, 18};
            cellmap[9].assign(cm9, cm9 + 5);

            int cm10[] = {9, 10, 11, 0, 1, 2, 18, 19, 20};
            cellmap[10].assign(cm10, cm10 + 9);

            int cm11[] = {10, 11, 14, 2, 20};
            cellmap[11].assign(cm11, cm11 + 5);

            int cm12[] = {9, 12, 15, 0, 3, 6, 18, 21, 24};
            cellmap[12].assign(cm12, cm12 + 9);

            int cm13[] = {13, 13, 13, 13};
            cellmap[13].assign(cm13, cm13 + 4);

            int cm14[] = {11, 14, 17, 2, 5, 8, 20, 23, 26};
            cellmap[14].assign(cm14, cm14 + 9);

            int cm15[] = {12, 15, 16, 6, 24};
            cellmap[15].assign(cm15, cm15 + 5);

            int cm16[] = {15, 16, 17, 6, 7, 8, 24, 25, 26};
            cellmap[16].assign(cm16, cm16 + 9);

            int cm17[] = {14, 17, 16, 8, 26};
            cellmap[17].assign(cm17, cm17 + 5);

            int cm18[] = {18, 19, 21, 9};
            cellmap[18].assign(cm18, cm18 + 4);

            int cm19[] = {18, 19, 20, 10};
            cellmap[19].assign(cm19, cm19 + 4);

            int cm20[] = {19, 20, 23, 11};
            cellmap[20].assign(cm20, cm20 + 4);

            int cm21[] = {18, 21, 24, 12};
            cellmap[21].assign(cm21, cm21 + 4);

            int cm22[] = {18, 19, 20, 21, 22, 23, 24, 25, 26};
            cellmap[22].assign(cm22, cm22 + 9);

            int cm23[] = {20, 23, 26, 14};
            cellmap[23].assign(cm23, cm23 + 4);

            int cm24[] = {21, 24, 25, 15};
            cellmap[24].assign(cm24, cm24 + 4);

            int cm25[] = {24, 25, 26, 16};
            cellmap[25].assign(cm25, cm25 + 4);

            int cm26[] = {23, 25, 26, 17};
            cellmap[26].assign(cm26, cm26 + 4);
        }

        // get final scores
        VInt finalscores(scores.size(), 0);
        for (int j = 0; j < (int)finalscores.size(); ++j) {
            // calculate score
            for (int k = 0; k < (int)cellmap[j].size(); ++k) {
                finalscores[j] += scores[cellmap[j][k]];
            }

            // not move to a structural cell
            if (scores[j] < 0) {
                // self is structure, score always < 0
                finalscores[j] = -1;
            } else {
                // self is not a structure, but score always >= 0
                if (finalscores[j] < 0) {
                    finalscores[j] = 0;
                }
            }
        }

        // find the cell with highest score
        int high = -100;
        ind = index;
        for (int i = 0; i < (int)finalscores.size(); ++i) {
            if (high < finalscores[i]) {
                high = finalscores[i];
                ind = indices[i];
            }
        }
        if (ind == index) continue;

        // get new  cell
        auto& new_cell = bcells[ind];
        if (new_cell.getType() == BACKGROUND_STRUCTURE) {
            continue;
        }

        // if an new empty cell
        if (new_cell.getNodes().empty()) {
            // get corners
            indices.clear();
            getCorners(ind, 1, indices);

            // set corners
            for (int k = 0; k < (int)indices.size(); ++k) {
                BNode& bnode = bnodes[indices[k]];
                if (bnode.size() == 0) {
                    bnode.setType(BACKGROUND_FLUID);
                }
                new_cell.addNode(&bnode, indices[k]);
            }
        }

        // get averate structural velocity
        VDouble svel(ndm);
        int numv = 0;
        for (auto* bnode : cell.getNodes()) {
            if (bnode->getType() != BACKGROUND_STRUCTURE) {
                continue;
            }
            auto& vel = bnode->getVel();
            for (int i = 0; i < bnode->size(); ++i) {
                svel += vel[i];
                ++numv;
            }
        }
        svel /= numv;
        double mag_svel = normVDouble(svel);

        // move the particles
        VDouble crds;
        getCrds(index, crds);
        auto& pts = cell.getPts();
        for (int i = 0; i < (int)pts.size(); ++i) {
            Particle* pt = pts[i];
            auto pcrds = pt->getCrds();

            // move the particle
            VDouble newcrds;
            getCrds(ind, newcrds);
            pcrds -= crds;
            pcrds += newcrds;
            pt->moveTo(pcrds, 0.0);
            auto& pvel = pt->getVel();
            double mag_pvel = normVDouble(pvel);
            if (mag_svel > mag_pvel) {
                pt->setVelOnly(svel);
            }

            new_cell.add(pt);
        }

        cell.clearParticles();
    }
    return 0;
}

int BackgroundMesh::gridFluid() {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();

    // store cells in a vector
    std::vector<BCell*> cells;
    VVInt indices;
    cells.reserve(bcells.size());
    indices.reserve(bcells.size());
    for (std::map<VInt, BCell>::iterator it = bcells.begin();
         it != bcells.end(); ++it) {
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
    VVInt elends(numele * cells.size());
    VInt gtags(numele * cells.size());
#pragma omp parallel for
    for (int j = 0; j < (int)cells.size(); ++j) {
        // structural cell
        if (cells[j]->getType() == BACKGROUND_STRUCTURE) continue;

        // find the group of this mesh
        std::map<int, int> numpts;
        auto& pts = cells[j]->getPts();
        for (int i = 0; i < (int)pts.size(); ++i) {
            numpts[pts[i]->getGroupTag()] += 1;
        }
        int num = 0;
        int gtag = 0;
        for (std::map<int, int>::iterator it = numpts.begin();
             it != numpts.end(); ++it) {
            if (num < it->second) {
                num = it->second;
                gtag = it->first;
            }
        }
        for (int i = 0; i < numele; ++i) {
            gtags[numele * j + i] = gtag;
        }

        // add to elenodes
        auto& bnodes = cells[j]->getNodes();
        VInt cnodes(bnodes.size());
        for (int i = 0; i < (int)cnodes.size(); ++i) {
            if (bnodes[i]->getTags().empty()) {
                opserr
                    << "WARNING: failed to be fluid node -- gridFluid\n";
                continue;
            }
            auto& tags = bnodes[i]->getTags();
            cnodes[i] = tags[0];
        }
        for (int i = 0; i < numele; ++i) {
            elends[numele * j + i].resize(numelenodes);
        }
        if (ndm == 2) {
            elends[numele * j][0] = cnodes[0];
            elends[numele * j][1] = cnodes[3];
            elends[numele * j][2] = cnodes[2];

            elends[numele * j + 1][0] = cnodes[0];
            elends[numele * j + 1][1] = cnodes[1];
            elends[numele * j + 1][2] = cnodes[3];

        } else if (ndm == 3) {
            elends[numele * j][0] = cnodes[0];
            elends[numele * j][1] = cnodes[7];
            elends[numele * j][2] = cnodes[4];
            elends[numele * j][3] = cnodes[5];

            elends[numele * j + 1][0] = cnodes[0];
            elends[numele * j + 1][1] = cnodes[1];
            elends[numele * j + 1][2] = cnodes[7];
            elends[numele * j + 1][3] = cnodes[5];

            elends[numele * j + 2][0] = cnodes[0];
            elends[numele * j + 2][1] = cnodes[1];
            elends[numele * j + 2][2] = cnodes[3];
            elends[numele * j + 2][3] = cnodes[7];

            elends[numele * j + 3][0] = cnodes[0];
            elends[numele * j + 3][1] = cnodes[3];
            elends[numele * j + 3][2] = cnodes[2];
            elends[numele * j + 3][3] = cnodes[7];

            elends[numele * j + 4][0] = cnodes[6];
            elends[numele * j + 4][1] = cnodes[0];
            elends[numele * j + 4][2] = cnodes[7];
            elends[numele * j + 4][3] = cnodes[4];

            elends[numele * j + 5][0] = cnodes[6];
            elends[numele * j + 5][1] = cnodes[0];
            elends[numele * j + 5][2] = cnodes[2];
            elends[numele * j + 5][3] = cnodes[7];
        }
    }

    // get particle group tags
    std::map<int, ID> elenodes;
    for (int i = 0; i < (int)elends.size(); ++i) {
        // no elenodes, no element
        if (elends[i].empty()) continue;

        // if all nodes are fluid and fixed nodes
        ID& nds = elenodes[gtags[i]];
        for (int j = 0; j < (int)elends[i].size(); ++j) {
            nds[nds.Size()] = elends[i][j];
        }
    }

    // create elements
    for (std::map<int, ID>::iterator it = elenodes.begin();
         it != elenodes.end(); ++it) {
        ParticleGroup* group =
            dynamic_cast<ParticleGroup*>(OPS_getMesh(it->first));
        if (group == 0) {
            opserr << "WARNING: failed to get particle group -- "
                      "BgMesh::gridFluid\n";
            return -1;
        }
        group->setEleNodes(it->second);

        if (group->newElements(it->second) < 0) {
            opserr << "WARNING: failed to create elements for mesh ";
            opserr << group->getTag() << " -- BgMesh::gridFluid\n";
            return -1;
        }
    }

    return 0;
}

int BackgroundMesh::gridFSInoDT() {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();

    // store cells in a vector
    std::vector<BCell*> cells;
    for (auto& item : bcells) {
        auto& cell = item.second;
        if (cell.getType() == BACKGROUND_STRUCTURE) {
            cells.push_back(&cell);
        }
    }

    // create elements in each cell
    int numele = 0;
    if (ndm == 2) {
        numele = 2;
    } else if (ndm == 3) {
        numele = 6;
    }
    VVInt elends(numele * cells.size());
    VInt gtags(numele * cells.size());
#pragma omp parallel for
    for (int j = 0; j < (int)cells.size(); ++j) {
        // get indices
        auto& cindices = cells[j]->getIndices();

        // get bnodes
        auto& cbnodes = cells[j]->getNodes();

        // node tags and sids
        VInt tags(cbnodes.size(), -1), sids(cbnodes.size()),
            types(cbnodes.size());
        VVDouble ndcrds(cbnodes.size());
        for (int i = 0; i < (int)tags.size(); ++i) {
            types[i] = cbnodes[i]->getType();
            if (types[i] != BACKGROUND_FIXED) {
                auto& ctags = cbnodes[i]->getTags();
                tags[i] = ctags[0];
                auto& ids = cbnodes[i]->getSid();
                sids[i] = ids[0];
                auto& crdsn = cbnodes[i]->getCrds();
                ndcrds[i] = crdsn[0];
            }
        }

        // get fluid and structural local indices
        VInt flocal, slocal, fixlocal;
        for (int i = 0; i < (int)types.size(); ++i) {
            int type = types[i];
            if (type == BACKGROUND_STRUCTURE) {
                slocal.push_back(i);
            } else if (type == BACKGROUND_FLUID) {
                flocal.push_back(i);
            } else {
                fixlocal.push_back(i);
            }
        }

        // contact cell: must fully engaged
        if (ndm == 2 && slocal.size() == 4) {
            VInt csnds(3), csids(3);
            if (slocal.size() == 4) {
                // two contact elements
                csnds[0] = tags[0];
                csnds[1] = tags[3];
                csnds[2] = tags[2];
                csids[0] = sids[0];
                csids[1] = sids[3];
                csids[2] = sids[2];
                if (createContact(csnds, csids, elends[numele * j])) {
                    gtags[numele * j] = contact_tag;
                }

                csnds[0] = tags[0];
                csnds[1] = tags[1];
                csnds[2] = tags[3];
                csids[0] = sids[0];
                csids[1] = sids[1];
                csids[2] = sids[3];
                if (createContact(csnds, csids, elends[numele * j + 1])) {
                    gtags[numele * j + 1] = contact_tag;
                }
                continue;
            }
        } else if (ndm == 3 && slocal.size() == 8) {
            // 3D contact element: TODO
            continue;
        }

        // get surrounding cells
        VInt min = cindices[0];
        min -= 1;
        VInt max = cindices[0];
        max += 2;

        // gather particles
        VParticle pts;
        gatherParticles(min, max, pts);
        if (pts.empty()) {
            continue;
        }

        // find the group of this mesh
        std::map<int, int> numpts;
        for (int i = 0; i < (int)pts.size(); ++i) {
            numpts[pts[i]->getGroupTag()] += 1;
        }
        int num = 0;
        int gtag = 0;
        for (std::map<int, int>::iterator it = numpts.begin();
             it != numpts.end(); ++it) {
            if (num < it->second) {
                num = it->second;
                gtag = it->first;
            }
        }
        for (int i = 0; i < numele; ++i) {
            gtags[numele * j + i] = gtag;
        }

        // add to elenodes
        if (ndm == 2) {
            if (types[1] == BACKGROUND_STRUCTURE &&
                types[2] == BACKGROUND_STRUCTURE) {
                if (types[0] != BACKGROUND_FIXED) {
                    if (!check_area(ndcrds[0], ndcrds[1], ndcrds[2])) {
                        elends[numele * j].push_back(tags[0]);
                        elends[numele * j].push_back(tags[1]);
                        elends[numele * j].push_back(tags[2]);
                    }
                }
                if (types[3] != BACKGROUND_FIXED) {
                    if (!check_area(ndcrds[1], ndcrds[3], ndcrds[2])) {
                        elends[numele * j + 1].push_back(tags[1]);
                        elends[numele * j + 1].push_back(tags[3]);
                        elends[numele * j + 1].push_back(tags[2]);
                    }
                }

            } else if (types[0] == BACKGROUND_STRUCTURE &&
                       types[3] == BACKGROUND_STRUCTURE) {
                if (types[1] != BACKGROUND_FIXED) {
                    if (!check_area(ndcrds[0], ndcrds[1], ndcrds[3])) {
                        elends[numele * j].push_back(tags[0]);
                        elends[numele * j].push_back(tags[1]);
                        elends[numele * j].push_back(tags[3]);
                    }
                }
                if (types[2] != BACKGROUND_FIXED) {
                    if (!check_area(ndcrds[0], ndcrds[3], ndcrds[2])) {
                        elends[numele * j + 1].push_back(tags[0]);
                        elends[numele * j + 1].push_back(tags[3]);
                        elends[numele * j + 1].push_back(tags[2]);
                    }
                }

            } else if (types[1] != BACKGROUND_FIXED &&
                       types[2] != BACKGROUND_FIXED) {
                if (types[0] != BACKGROUND_FIXED) {
                    if (!check_area(ndcrds[0], ndcrds[1], ndcrds[2])) {
                        elends[numele * j].push_back(tags[0]);
                        elends[numele * j].push_back(tags[1]);
                        elends[numele * j].push_back(tags[2]);
                    }
                }
                if (types[3] != BACKGROUND_FIXED) {
                    if (!check_area(ndcrds[1], ndcrds[3], ndcrds[2])) {
                        elends[numele * j + 1].push_back(tags[1]);
                        elends[numele * j + 1].push_back(tags[3]);
                        elends[numele * j + 1].push_back(tags[2]);
                    }
                }

            } else if (types[0] != BACKGROUND_FIXED &&
                       types[3] != BACKGROUND_FIXED) {
                if (types[1] != BACKGROUND_FIXED) {
                    if (!check_area(ndcrds[0], ndcrds[1], ndcrds[3])) {
                        elends[numele * j].push_back(tags[0]);
                        elends[numele * j].push_back(tags[1]);
                        elends[numele * j].push_back(tags[3]);
                    }
                }
                if (types[2] != BACKGROUND_FIXED) {
                    if (!check_area(ndcrds[0], ndcrds[3], ndcrds[2])) {
                        elends[numele * j + 1].push_back(tags[0]);
                        elends[numele * j + 1].push_back(tags[3]);
                        elends[numele * j + 1].push_back(tags[2]);
                    }
                }
            }

        } else if (ndm == 3) {
            VInt prism(6);
            VVInt tets;
            bool incl1, incl4;

            int face[6][4] = {{1, 2, 5, 6}, {0, 3, 4, 7}, {1, 3, 4, 6},
                              {0, 2, 5, 7}, {2, 3, 4, 5}, {0, 1, 6, 7}};
            int prisms[12][6] = {{5, 4, 6, 1, 0, 2}, {6, 7, 5, 2, 3, 1},
                                 {4, 6, 7, 0, 2, 3}, {7, 5, 4, 3, 1, 0},
                                 {4, 5, 1, 6, 7, 3}, {1, 0, 4, 3, 2, 6},
                                 {0, 4, 5, 2, 6, 7}, {5, 1, 0, 7, 3, 2},
                                 {3, 1, 5, 2, 0, 4}, {5, 7, 3, 4, 6, 2},
                                 {1, 5, 7, 0, 4, 6}, {7, 3, 1, 6, 2, 0}};

            // if find a splitting
            bool find = false;

            // check which face is structure
            for (int i = 0; i < 6; ++i) {
                find = true;
                for (int k = 0; k < 4; ++k) {
                    if (types[face[i][k]] != BACKGROUND_STRUCTURE) {
                        find = false;
                    }
                }

                // if it's structure
                if (find) {
                    // set prism
                    for (int k = 0; k < 2; ++k) {
                        for (int m = 0; m < 6; ++m) {
                            prism[m] = prisms[2 * i + k][m];
                        }

                        // check if include corner 1
                        incl1 = false;
                        if (types[prism[1]] != BACKGROUND_FIXED) {
                            incl1 = true;
                        }

                        // check if include corner 4
                        incl4 = false;
                        if (types[prism[4]] != BACKGROUND_FIXED) {
                            incl4 = true;
                        }

                        // split the prism
                        splitPrism(prism, tets, incl1, incl4);

                        // add ele nodes
                        for (int m = 0; m < (int)tets.size(); ++m) {
                            if (!check_vol(ndcrds[tets[m][0]],
                                           ndcrds[tets[m][1]],
                                           ndcrds[tets[m][2]],
                                           ndcrds[tets[m][3]])) {
                                for (int p : tets[m]) {
                                    elends[numele * j + m].push_back(
                                        tags[p]);
                                }
                            }
                        }
                    }
                    break;
                }
            }

            // check which face is not fixed
            if (!find) {
                for (int i = 0; i < 6; ++i) {
                    find = true;
                    for (int k = 0; k < 4; ++k) {
                        if (types[face[i][k]] == BACKGROUND_FIXED) {
                            find = false;
                        }
                    }

                    // if it's not fixed
                    if (find) {
                        // set prism
                        for (int k = 0; k < 2; ++k) {
                            for (int m = 0; m < 6; ++m) {
                                prism[m] = prisms[2 * i + k][m];
                            }

                            // check if include corner 1
                            incl1 = false;
                            if (types[prism[1]] != BACKGROUND_FIXED) {
                                incl1 = true;
                            }

                            // check if include corner 4
                            incl4 = false;
                            if (types[prism[4]] != BACKGROUND_FIXED) {
                                incl4 = true;
                            }

                            // split the prism
                            splitPrism(prism, tets, incl1, incl4);

                            // add ele nodes
                            for (int m = 0; m < (int)tets.size(); ++m) {
                                if (!check_vol(ndcrds[tets[m][0]],
                                               ndcrds[tets[m][1]],
                                               ndcrds[tets[m][2]],
                                               ndcrds[tets[m][3]])) {
                                    for (int p : tets[m]) {
                                        elends[numele * j + m].push_back(
                                            tags[p]);
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }
    }

    // get particle group tags
    std::map<int, ID> elenodes;
    int nextEletag = Mesh::nextEleTag();
    VInt oldContactEles = contactEles;
    VInt removedEles(oldContactEles.size(), 1);
    contactEles.clear();
    for (int i = 0; i < (int)elends.size(); ++i) {
        // no elenodes, no element
        if (elends[i].empty()) continue;

        // contact element
        if (gtags[i] == contact_tag) {
            if (ndm == 2) {
                if (elends[i].size() != 3) {
                    opserr << "WARNING: 2D contact should have 3 "
                              "nodes\n";
                    return -1;
                }

                // check if exists
                bool created = false;
                for (int j = 0; j < (int)oldContactEles.size(); ++j) {
                    if (removedEles[j] == 0) continue;
                    Element* ele = domain->getElement(oldContactEles[j]);
                    if (ele == 0) continue;
                    const ID& contactNodes = ele->getExternalNodes();
                    if (contactNodes(0) == elends[i][0] &&
                        contactNodes(1) == elends[i][1] &&
                        contactNodes(2) == elends[i][2]) {
                        created = true;
                        removedEles[j] = 0;
                        contactEles.push_back(oldContactEles[j]);
                        break;
                    }
                }
                if (created) continue;

                if (contactData[0] <= 0 || contactData[1] <= 0 ||
                    contactData[2] < 0 || contactData[3] < 0 ||
                    contactData[4] <= 0 || contactData[6] <= 0 ||
                    contactData[7] <= 0) {
                    opserr << "WARNING: contact data is not "
                              "correctly "
                              "set\n";
                    return -1;
                }
                Element* ele = new PFEMContact2D(
                    nextEletag, elends[i][0], elends[i][1], elends[i][2],
                    contactData[0], contactData[1], contactData[2],
                    contactData[3], contactData[4], contactData[5],
                    contactData[6], contactData[7]);
                if (ele == 0) {
                    opserr << "WARNING: failed to create contact "
                              "element\n";
                    return -1;
                }
                if (domain->addElement(ele) == false) {
                    opserr << "WARNING: failed to add element "
                           << nextEletag << "\n";
                    delete ele;
                    return -1;
                }
                contactEles.push_back(nextEletag);
                nextEletag += 1;

            } else if (ndm == 3) {
                if (elends[i].size() != 4) {
                    opserr << "WARNING: 3D contact should have 4 "
                              "nodes\n";
                    return -1;
                }
                opserr << "WARNING: 3D contact element hasn't been "
                          "developed\n";
            }
            continue;
        }

        // if all nodes are fluid nodes
        ID& nds = elenodes[gtags[i]];
        for (int j = 0; j < (int)elends[i].size(); ++j) {
            nds[nds.Size()] = elends[i][j];
        }
    }

    // remove old contact elements
    for (int i = 0; i < (int)oldContactEles.size(); ++i) {
        if (removedEles[i] == 1) {
            opserr << "old contact element " << oldContactEles[i] << ": ";
            Element* ele = domain->removeElement(oldContactEles[i]);
            opserr << ele->getExternalNodes() << " is removed\n";
            if (ele != 0) {
                delete ele;
            }
        }
    }

    // create elements
    for (auto& item : elenodes) {
        ParticleGroup* group =
            dynamic_cast<ParticleGroup*>(OPS_getMesh(item.first));
        if (group == 0) {
            opserr << "WARNING: failed to get particle group -- "
                      "BgMesh::gridFluid\n";
            return -1;
        }
        group->setEleNodes(item.second);

        if (group->newElements(item.second) < 0) {
            opserr << "WARNING: failed to create elements for mesh ";
            opserr << group->getTag() << " -- BgMesh::gridFluid\n";
            return -1;
        }
    }

    return 0;
}

int BackgroundMesh::gridFSI() {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();

    TriangleMeshGenerator gen;
    TetMeshGenerator tetgen;

    // gather bnodes
    std::map<VInt, BNode*> fsibnodes;
    for (std::map<VInt, BCell>::iterator it = bcells.begin();
         it != bcells.end(); ++it) {
        // only for structural cells
        BCell& bcell = it->second;
        if (bcell.getType() == BACKGROUND_FLUID) continue;

        // get bnode
        auto& bnodes = bcell.getNodes();
        auto& bindices = bcell.getIndices();
        for (int j = 0; j < (int)bnodes.size(); ++j) {
            BNode* bnode = bnodes[j];
            VInt bindex = bindices[j];
            if (bnode == 0) {
                opserr << "WARNING: failed to get bnode -- gridFSI\n";
                return -1;
            }
            fsibnodes[bindex] = bnode;
        }
    }

    // add points
    VInt ndtags, ndtypes, ndsids;
    VVInt ndindex;
    VDouble min, max;
    for (std::map<VInt, BNode*>::iterator it = fsibnodes.begin();
         it != fsibnodes.end(); ++it) {
        VInt bindex = it->first;
        BNode* bnode = it->second;

        const VInt& tags = bnode->getTags();
        const VVDouble& crdsn = bnode->getCrds();
        int type = bnode->getType();
        const VInt& sid = bnode->getSid();

        for (int i = -1; i < (int)tags.size(); ++i) {
            // -1 is only for fixed bnode
            if (i == -1 && type != BACKGROUND_FIXED) {
                continue;
            }
            if (i == -1) {
                ndtags.push_back(0);
                ndsids.push_back(0);
            } else {
                ndtags.push_back(tags[i]);
                ndsids.push_back(sid[i]);
            }
            ndtypes.push_back(type);
            ndindex.push_back(bindex);

            VDouble crds;
            if (i == -1) {
                getCrds(bindex, crds);
            } else {
                crds = crdsn[i];
            }

            if (ndm == 2) {
                gen.addPoint(crds[0], crds[1]);
            } else if (ndm == 3) {
                tetgen.addPoint(crds[0], crds[1], crds[2], 0);
            }

            // max, min coordinates
            if (min.empty() || max.empty()) {
                min = crds;
                max = crds;
            } else {
                for (int j = 0; j < ndm; ++j) {
                    if (min[j] > crds[j]) min[j] = crds[j];
                    if (max[j] < crds[j]) max[j] = crds[j];
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
    if (gen.getNumPoints() < 3 && tetgen.getNumPoints() < 4) {
        return 0;
    }
    if (min.empty() || max.empty()) {
        return 0;
    }

    // add extra points to avoid small triangles on the edge
    min -= 10 * bsize;
    max += 10 * bsize;
    if (ndm == 2) {
        gen.addPoint(min[0], min[1]);
        gen.addPoint(max[0], min[1]);
        gen.addPoint(min[0], max[1]);
        gen.addPoint(max[0], max[1]);
    } else if (ndm == 3) {
        tetgen.addPoint(min[0], min[1], min[2], 0);
        tetgen.addPoint(min[0], max[1], min[2], 0);
        tetgen.addPoint(min[0], max[1], max[2], 0);
        tetgen.addPoint(min[0], min[1], max[2], 0);
        tetgen.addPoint(max[0], min[1], min[2], 0);
        tetgen.addPoint(max[0], max[1], min[2], 0);
        tetgen.addPoint(max[0], max[1], max[2], 0);
        tetgen.addPoint(max[0], min[1], max[2], 0);
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
    for (int i = 0; i < numele; ++i) {
        // get points
        VInt tri;
        if (ndm == 2) {
            tri.resize(3);
            gen.getTriangle(i, tri[0], tri[1], tri[2]);
        } else if (ndm == 3) {
            tri.resize(4);
            tetgen.getTet(i, tri[0], tri[1], tri[2], tri[3]);
        }

        // check if connect to extra points
        bool extra = false;
        for (int j = 0; j < (int)tri.size(); ++j) {
            if (tri[j] >= numpoints) {
                extra = true;
                break;
            }
        }
        if (extra) continue;

        // get fluid and structural indices
        VVInt findex, sindex;
        bool fixed = false;
        for (int j = 0; j < (int)tri.size(); ++j) {
            if (ndtypes[tri[j]] == BACKGROUND_STRUCTURE) {
                sindex.push_back(ndindex[tri[j]]);
            } else if (ndtypes[tri[j]] == BACKGROUND_FIXED) {
                fixed = true;
            } else {
                findex.push_back(ndindex[tri[j]]);
            }
        }
        if (fixed) continue;

        // if all structure, create contact elements
        if (sindex.size() == tri.size()) {
            VInt snds(tri.size()), sids(tri.size());
            for (int j = 0; j < (int)tri.size(); ++j) {
                snds[j] = ndtags[tri[j]];
                sids[j] = ndsids[tri[j]];
            }
            if (createContact(snds, sids, elends[i]) == 0) {
                gtags[i] = contact_tag;
            }
            continue;
        }

        // get min and max ind
        VInt maxind = ndindex[tri[0]];
        VInt minind = ndindex[tri[0]];
        for (int k = 0; k < ndm; ++k) {
            for (int j = 1; j < (int)tri.size(); ++j) {
                if (minind[k] > ndindex[tri[j]][k]) {
                    minind[k] = ndindex[tri[j]][k];
                } else if (maxind[k] < ndindex[tri[j]][k]) {
                    maxind[k] = ndindex[tri[j]][k];
                }
            }
        }

        // check if triangle out side of FSI
        bool outside = false;
        if (findex.size() == tri.size()) {
            VInt currind(ndm);
            if (ndm == 2) {
                for (int j = minind[0]; j < maxind[0]; ++j) {
                    for (int k = minind[1]; k < maxind[1]; ++k) {
                        currind[0] = j;
                        currind[1] = k;
                        std::map<VInt, BCell>::iterator it =
                            bcells.find(currind);
                        if (it == bcells.end()) {
                            outside = true;
                            break;
                        } else {
                            if (it->second.getType() == BACKGROUND_FLUID) {
                                outside = true;
                                break;
                            }
                        }
                        if (outside) break;
                    }
                }
            } else if (ndm == 3) {
                for (int j = minind[0]; j < maxind[0]; ++j) {
                    for (int k = minind[1]; k < maxind[1]; ++k) {
                        for (int l = minind[2]; l < maxind[2]; ++l) {
                            currind[0] = j;
                            currind[1] = k;
                            currind[2] = l;
                            std::map<VInt, BCell>::iterator it =
                                bcells.find(currind);
                            if (it == bcells.end()) {
                                outside = true;
                                break;
                            } else {
                                if (it->second.getType() ==
                                    BACKGROUND_FLUID) {
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

        // point coordinates
        VVDouble ptcrds(tri.size());
        for (int j = 0; j < (int)tri.size(); ++j) {
            ptcrds[j].resize(ndm);
            int mark;
            if (ndm == 2) {
                gen.getPoint(tri[j], ptcrds[j][0], ptcrds[j][1], mark);
            } else if (ndm == 3) {
                tetgen.getPoint(tri[j], ptcrds[j][0], ptcrds[j][1],
                                ptcrds[j][2], mark);
            }
        }

        // gather particles
        VParticle tripts;
        for (int j = 0; j < (int)findex.size(); ++j) {
            VInt ind = findex[j];
            ind -= 1;
            VVInt indices;
            getCorners(ind, 1, indices);
            for (int k = 0; k < (int)indices.size(); ++k) {
                std::map<VInt, BCell>::iterator cellit =
                    bcells.find(indices[k]);
                if (cellit == bcells.end()) continue;
                if (cellit->second.getType() == BACKGROUND_STRUCTURE)
                    continue;
                if (cellit->second.getPts().empty()) continue;
                const VParticle& pts = cellit->second.getPts();
                tripts.insert(tripts.end(), pts.begin(), pts.end());
            }
        }

        if (tripts.empty()) {
            continue;
        }

        // in-element check
        bool zerovol = false;
        if (ndm == 2) {
            check_area(ptcrds[0], ptcrds[1], ptcrds[2]);
        } else {
            check_vol(ptcrds[0], ptcrds[1], ptcrds[2], ptcrds[3]);
        }
        if (zerovol) {
            continue;
        }

        VParticle newtripts;
        newtripts = tripts;

        if (newtripts.empty()) continue;

        // find the group of this mesh
        std::map<int, int> numpts;
        for (int j = 0; j < (int)newtripts.size(); ++j) {
            numpts[newtripts[j]->getGroupTag()] += 1;
        }
        int num = 0;
        for (std::map<int, int>::iterator it = numpts.begin();
             it != numpts.end(); ++it) {
            if (num < it->second) {
                num = it->second;
                gtags[i] = it->first;
            }
        }

        // add to elenodes
        elends[i].resize(tri.size());
        for (int j = 0; j < (int)tri.size(); ++j) {
            elends[i][j] = ndtags[tri[j]];
        }
    }

    // get particle group tags
    std::map<int, ID> elenodes;
    int nextEletag = Mesh::nextEleTag();
    VInt oldContactEles = contactEles;
    VInt removedEles(oldContactEles.size(), 1);
    contactEles.clear();
    for (int i = 0; i < (int)elends.size(); ++i) {
        // no elenodes, no element
        if (elends[i].empty()) continue;

        // contact element
        if (gtags[i] == contact_tag) {
            if (ndm == 2) {
                if (elends[i].size() != 3) {
                    opserr << "WARNING: 2D contact should have 3 "
                              "nodes\n";
                    return -1;
                }

                // check if exists
                bool created = false;
                for (int j = 0; j < (int)oldContactEles.size(); ++j) {
                    if (removedEles[j] == 0) continue;
                    Element* ele = domain->getElement(oldContactEles[j]);
                    if (ele == 0) continue;
                    const ID& contactNodes = ele->getExternalNodes();
                    if (contactNodes(0) == elends[i][0] &&
                        contactNodes(1) == elends[i][1] &&
                        contactNodes(2) == elends[i][2]) {
                        created = true;
                        removedEles[j] = 0;
                        contactEles.push_back(oldContactEles[j]);
                        break;
                    }
                }
                if (created) continue;

                if (contactData[0] <= 0 || contactData[1] <= 0 ||
                    contactData[2] < 0 || contactData[3] < 0 ||
                    contactData[4] <= 0 || contactData[6] <= 0 ||
                    contactData[7] <= 0) {
                    opserr << "WARNING: contact data is not "
                              "correctly "
                              "set\n";
                    return -1;
                }
                Element* ele = new PFEMContact2D(
                    nextEletag, elends[i][0], elends[i][1], elends[i][2],
                    contactData[0], contactData[1], contactData[2],
                    contactData[3], contactData[4], contactData[5],
                    contactData[6], contactData[7]);
                if (ele == 0) {
                    opserr << "WARNING: failed to create contact "
                              "element\n";
                    return -1;
                }
                if (domain->addElement(ele) == false) {
                    opserr << "WARNING: failed to add element "
                           << nextEletag << "\n";
                    delete ele;
                    return -1;
                }
                contactEles.push_back(nextEletag);
                nextEletag += 1;

            } else if (ndm == 3) {
                if (elends[i].size() != 4) {
                    opserr << "WARNING: 3D contact should have 4 "
                              "nodes\n";
                    return -1;
                }
                opserr << "WARNING: 3D contact element hasn't been "
                          "developed\n";
            }
            continue;
        }

        // if all nodes are fluid and fixed nodes
        ID& nds = elenodes[gtags[i]];
        for (int j = 0; j < (int)elends[i].size(); ++j) {
            nds[nds.Size()] = elends[i][j];
        }
    }

    // remove old contact elements
    for (int i = 0; i < (int)oldContactEles.size(); ++i) {
        if (removedEles[i] == 1) {
            opserr << "old contact element " << oldContactEles[i] << ": ";
            Element* ele = domain->removeElement(oldContactEles[i]);
            opserr << ele->getExternalNodes() << " is removed\n";
            if (ele != 0) {
                delete ele;
            }
        }
    }

    // create elements
    for (std::map<int, ID>::iterator it = elenodes.begin();
         it != elenodes.end(); ++it) {
        ParticleGroup* group =
            dynamic_cast<ParticleGroup*>(OPS_getMesh(it->first));
        if (group == 0) {
            opserr << "WARNING: failed to get particle group -- "
                      "BgMesh::gridFSI\n";
            return -1;
        }
        group->addEleNodes(it->second);

        if (group->newElements(it->second) < 0) {
            opserr << "WARNING: failed to create elements for mesh ";
            opserr << group->getTag() << " -- BgMesh::gridFSI\n";
            return -1;
        }
    }

    return 0;
}

// time, vx, vy, (vz), top, bottom, left, right, dt, numIter
//

int BackgroundMesh::record(bool init) {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();

    double timestamp = domain->getCurrentTime();
    currentTime = timestamp;

    // record
    for (int i = 0; i < (int)recorders.size(); ++i) {
        if (recorders[i] != 0) {
            if (recorders[i]->record(domain->getCommitTag(),
                                     currentTime)) {
                opserr << "WARNING: failed to record -- "
                          "BgMesh::gridEles\n";
                return -1;
            }
        }
    }
    if (init) {
        return 0;
    }

    // record time
    if (theFile.good() == false) {
        return 0;
    }
    theFile << timestamp << " ";

    // record wave height and velocity
    for (int i = 0; i < (int)locs.size(); i += ndm) {
        // lower index
        VDouble crds(ndm);
        for (int j = 0; j < ndm; ++j) {
            crds[j] = locs[i + j];
        }

        VInt index;
        this->lowerIndex(crds, index);

        // shape function
        VDouble N;
        double hh = bsize;
        if (ndm == 2) {
            getNForRect(index[0] * bsize, index[1] * bsize, hh, hh,
                        crds[0], crds[1], N);
        } else if (ndm == 3) {
            getNForRect(index[0] * bsize, index[1] * bsize,
                        index[2] * bsize, hh, hh, hh, crds[0], crds[1],
                        crds[2], N);
        }

        // velocity
        VDouble vel(ndm);

        // get corners
        VVInt indices;
        getCorners(index, 1, indices);
        for (int k = 0; k < (int)indices.size(); ++k) {
            // get crds
            getCrds(indices[k], crds);

            // check bnode
            std::map<VInt, BNode>::iterator bit = bnodes.find(indices[k]);
            if (bit == bnodes.end()) continue;

            // get bnode
            BNode& bnode = bit->second;
            auto& vn = bnode.getVel();
            if (vn.empty()) {
                continue;
            }

            // get vel
            for (int j = 0; j < ndm; ++j) {
                vel[j] += N[k] * vn[0][j];
            }
        }

        // get the lowest cell
        double bottom = 0.0;
        for (int iy = lower[1]; iy <= upper[1]; ++iy) {
            VInt ind = index;
            ind[1] = iy;
            auto it = bcells.find(ind);
            if (it != bcells.end() &&
                it->second.getPts().empty() == false) {
                bottom = iy * bsize + 0.5 * bsize;
                break;
            }
        }

        // get the most left cell
        double left = 0.0;
        for (int ix = lower[0]; ix <= upper[0]; ++ix) {
            VInt ind = index;
            ind[0] = ix;
            auto it = bcells.find(ind);
            if (it != bcells.end() &&
                it->second.getPts().empty() == false) {
                left = ix * bsize + 0.5 * bsize;
                break;
            }
        }
        // get the highest cell
        double top = 0.0;
        for (int iy = upper[1]; iy >= lower[1]; --iy) {
            VInt ind = index;
            ind[1] = iy;
            auto it = bcells.find(ind);
            if (it != bcells.end() &&
                it->second.getPts().empty() == false) {
                top = (iy + 0.5) * bsize;
                break;
            }
        }

        // get the most right cell
        double right = 0.0;
        for (int ix = upper[0]; ix >= lower[0]; --ix) {
            VInt ind = index;
            ind[0] = ix;
            auto it = bcells.find(ind);
            if (it != bcells.end() &&
                it->second.getPts().empty() == false) {
                right = (ix + 0.5) * bsize;
                break;
            }
        }

        // record velocity and wave
        for (int j = 0; j < ndm; ++j) {
            theFile << vel[j] << " ";
        }
        theFile << top << " " << bottom << " ";
        theFile << left << " " << right << " ";
    }

    // time step and iteration
    // theFile << ops_Dt << " " << OPS_numIter()<<"\n";
    theFile << ops_Dt << "\n";
    theFile.flush();

    return 0;
}

int BackgroundMesh::moveParticles() {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;
    int ndm = OPS_GetNDM();
    double dt = domain->getCurrentTime() - currentTime;

    // get current disp and velocity
    for (std::map<VInt, BNode>::iterator it = bnodes.begin();
         it != bnodes.end(); ++it) {
        BNode& bnode = it->second;
        VInt& tags = bnode.getTags();

        for (int i = 0; i < (int)bnode.size(); ++i) {
            Node* nd = domain->getNode(tags[i]);
            Pressure_Constraint* pc =
                domain->getPressure_Constraint(tags[i]);

            if (pc != 0) {
                auto& pn = bnode.getPressure();
                auto& dpn = bnode.getPdot();
                pn[i] = pc->getPressure();
                dpn[i] = pc->getPdot();
            }
            if (nd != 0) {
                const Vector& vel = nd->getTrialVel();
                const Vector& accel = nd->getTrialAccel();
                auto& vn = bnode.getVel();
                auto& dvn = bnode.getAccel();
                for (int j = 0; j < ndm; ++j) {
                    vn[i][j] = vel(j);
                    dvn[i][j] = accel(j);
                }
            }
        }
    }

    // store cells in a vector
    std::vector<BCell*> cells;
    VVInt indices;
    cells.reserve(bcells.size());
    indices.reserve(bcells.size());
    for (std::map<VInt, BCell>::iterator it = bcells.begin();
         it != bcells.end(); ++it) {
        indices.push_back(it->first);
        cells.push_back(&(it->second));
    }

    // move particles in each cell
    int res = 0;
#pragma omp parallel for
    for (int j = 0; j < (int)cells.size(); ++j) {
        // get particles in cell
        const VParticle& pts = cells[j]->getPts();

        // move the particle
        for (int i = 0; i < (int)pts.size(); ++i) {
            // set update state
            if (pts[i] == 0) continue;
            pts[i]->needUpdate(dt);

            // convect the particle
            if (convectParticle(pts[i], indices[j], numsub) < 0) {
                opserr << "WARNING: failed to convect particle";
                opserr << " -- BgMesh::moveParticles\n";
                res = -1;
                continue;
            }
        }
    }

    if (res < 0) return -1;

    return 0;
}

int BackgroundMesh::convectParticle(Particle* pt, VInt index, int nums) {
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    int ndm = OPS_GetNDM();

    // check dt
    double dt = pt->getDt();
    if (dt <= 0) {
        return 0;
    }

    // get corners
    VVInt indices, ndtags;
    VVDouble crds;
    std::vector<BackgroundType> types;
    VVDouble vels, dvns;
    VDouble pns, dpns;

    // convect in a cell
    while (pt->getDt() > 0) {
        // get subdt
        double subdt = dt / nums;
        if (pt->getDt() < subdt) {
            subdt = pt->getDt();
        }

        // particle crds
        const VDouble& pcrds = pt->getCrds();

        // new index
        VInt newIndex;
        lowerIndex(pcrds, newIndex);

        // update corners
        if (indices.empty() || newIndex != index) {
            index = newIndex;

            VVInt temp;
            getCorners(index, 1, temp);
            indices = temp;
            indices[2] = temp[3];
            indices[3] = temp[2];
            if (ndm == 3) {
                indices[6] = temp[7];
                indices[7] = temp[6];
            }

            // get corner coordinates, types, and velocities
            // pressures, alphas for structural damping
            ndtags.assign(indices.size(), VInt());
            crds.assign(indices.size(), VDouble());
            types.assign(indices.size(), BACKGROUND_FIXED);
            vels.assign(indices.size(), VDouble());
            dvns.assign(indices.size(), VDouble());
            pns.assign(indices.size(), 0.0);
            dpns.assign(indices.size(), 0.0);

            for (int i = 0; i < (int)indices.size(); ++i) {
                // get crds
                getCrds(indices[i], crds[i]);

                // check bnode
                auto it = bnodes.find(indices[i]);
                if (it == bnodes.end()) continue;

                // get bnode
                BNode& bnode = it->second;
                types[i] = bnode.getType();
                if (types[i] == BACKGROUND_FIXED) continue;

                // get node tags
                ndtags[i] = bnode.getTags();

                // get vn and dvn
                if (bnode.getTags().size() < 1) {
                    opserr << "WARNING: fluid bnode tags.size() < 1 ";
                    opserr << "-- BgMesh::convectParticle\n";
                    return -1;
                }
                auto& vn = bnode.getVel();
                auto& dvn = bnode.getAccel();
                auto& pn = bnode.getPressure();
                auto& dpn = bnode.getPdot();
                vels[i] = vn[0];
                dvns[i] = dvn[0];
                pns[i] = pn[0];
                dpns[i] = dpn[0];
            }
        }

        // get particle velocity and move
        if (interpolate(pt, indices, vels, dvns, pns, dpns, crds, types,
                        ndtags, subdt) < 0) {
            opserr << "WARNING: failed to interpolate particle "
                      "velocity";
            opserr << "-- BgMesh::convectParticle\n";
            return -1;
        }
    }

    return 0;
}

int BackgroundMesh::interpolate(Particle* pt, const VVInt& index,
                                const VVDouble& vels, const VVDouble& dvns,
                                const VDouble& pns, const VDouble& dpns,
                                const VVDouble& crds,
                                const std::vector<BackgroundType>& types,
                                const VVInt& ndtags, double dt) {
    int ndm = OPS_GetNDM();
    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    // check
    if (ndm == 2) {
        if (index.size() != 4) return 0;
        if (vels.size() != 4) return 0;
        if (pns.size() != 4) return 0;
        if (dpns.size() != 4) return 0;
        if (crds.size() != 4) return 0;
        if (types.size() != 4) return 0;
    } else if (ndm == 3) {
        if (index.size() != 8) return 0;
        if (vels.size() != 8) return 0;
        if (pns.size() != 8) return 0;
        if (dpns.size() != 8) return 0;
        if (crds.size() != 8) return 0;
        if (types.size() != 8) return 0;
    }

    const VDouble& pcrds = pt->getCrds();
    if ((int)pcrds.size() != ndm) {
        opserr << "WARNING: pcrds.size() != ndm -- "
                  "BgMesh::interpolate\n";
        return -1;
    }

    // get shape functions for pt
    VDouble N;
    if (ndm == 2) {
        double hx = (crds[1][0] + crds[2][0]) / 2.0 - crds[0][0];
        double hy = (crds[2][1] + crds[3][1]) / 2.0 - crds[0][1];
        getNForRect(crds[0][0], crds[0][1], hx, hy, pcrds[0], pcrds[1], N);
    } else if (ndm == 3) {
        double hx = (crds[1][0] + crds[2][0]) / 2.0 - crds[0][0];
        double hy = (crds[2][1] + crds[3][1]) / 2.0 - crds[0][1];
        double hz =
            (crds[4][2] + crds[5][2] + crds[6][2] + crds[7][2]) / 4.0 -
            crds[0][2];
        getNForRect(crds[0][0], crds[0][1], crds[0][2], hx, hy, hz,
                    pcrds[0], pcrds[1], pcrds[2], N);
    }

    // particle velocity
    VDouble pvel(ndm);
    VDouble pdvn(ndm);
    double ppre = 0.0, pdp = 0.0;
    double Nvsum = 0.0, Npsum = 0.0;
    for (int j = 0; j < (int)vels.size(); ++j) {
        // no vel at this corner
        if (vels[j].empty()) continue;

        // interpolate
        if (types[j] == BACKGROUND_FLUID) {
            // check surrounding nodes
            VInt ind = index[j];
            VVInt indices;
            ind -= 1;
            getCorners(ind, 2, indices);

            VDouble svel(ndm);
            double alphas = 0.0;
            int num_svel = 0, num_sid = 0;
            for (auto& indi : indices) {
                auto it = bnodes.find(indi);
                int size = 0;
                if (it != bnodes.end()) {
                    size = it->second.size();
                }
                for (int i = 0; i < size; ++i) {
                    auto& v = it->second.getVel();
                    svel += v[i];
                    ++num_svel;

                    auto& sid = it->second.getSid();
                    auto it_alpha = alphaS.find(sid[i]);
                    if (it_alpha != alphaS.end()) {
                        alphas += it_alpha->second;
                        ++num_sid;
                    }
                }
            }
            if (num_svel > 0) {
                svel /= num_svel;
            }
            if (num_sid > 0) {
                alphas /= num_sid;
            }

            // check surrounding cells
            getCorners(ind, 1, indices);
            VVInt sindices;
            for (int k = 0; k < (int)indices.size(); ++k) {
                auto it = bcells.find(indices[k]);
                if (it != bcells.end() &&
                    it->second.getType() == BACKGROUND_STRUCTURE) {
                    sindices.push_back(indices[k]);
                }
            }
            VInt num_equal(ndm);
            for (int k = 0; k < ndm; ++k) {
                for (int m = 0; m < (int)sindices.size(); ++m) {
                    int num = 0;
                    for (int n = 0; n < (int)sindices.size(); ++n) {
                        if (sindices[m][k] == sindices[n][k]) {
                            num += 1;
                        }
                    }
                    if (num_equal[k] < num) {
                        num_equal[k] = num;
                    }
                }
            }

            // interpolate
            for (int k = 0; k < ndm; ++k) {
                bool fixk = (ndm == 2 && num_equal[k] >= 2) ||
                            (ndm == 3 && num_equal[k] >= 4);
                if (num_svel > 0 && fixk && fabs(svel[k]) < tol) {
                    // k is fixed structure
                    pvel[k] += N[j] * svel[k];
                } else if (num_svel > 0 &&
                           fabs(svel[k]) > fabs(vels[j][k])) {
                    // k is faster structure
                    pvel[k] += N[j] * svel[k];
                } else if (num_svel > 0 &&
                           fabs(svel[k]) < fabs(vels[j][k])) {
                    // k is slower structure
                    pvel[k] += N[j] * (alphas * svel[k] +
                                       (1.0 - alphas) * vels[j][k]);
                } else {
                    // all others
                    pvel[k] += N[j] * vels[j][k];
                }
                if (pt->isUpdated() == false) {
                    pdvn[k] += N[j] * dvns[j][k];
                }
            }
            Nvsum += N[j];
        }
        if (pt->isUpdated() == false) {
            ppre += N[j] * pns[j];
            pdp += N[j] * dpns[j];
        }
        Npsum += N[j];
    }

    if (Nvsum > 0) {
        pvel /= Nvsum;
        pdvn /= Nvsum;
    } else {
        pvel = pt->getVel();
        pdvn = pt->getAccel();
    }

    if (Npsum > 0) {
        ppre /= Npsum;
        pdp /= Npsum;
    } else {
        ppre = pt->getPressure();
        pdp = pt->getPdot();
    }

    // particle displacement
    // cannot travel more than one cell
    VDouble newpcrds;
    bool travel_cell = false;
    while (newpcrds.empty() || travel_cell) {
        newpcrds = pvel;
        newpcrds *= dt;
        travel_cell = false;
        for (int i = 0; i < ndm; ++i) {
            if (newpcrds[i] > bsize) {
                travel_cell = true;
                dt *= 0.5;
                break;
            }
        }
    }

    // new particle coordinates
    newpcrds += pcrds;

    // // find new cell
    // VInt newindex;
    // lowerIndex(newpcrds, newindex);
    // auto it = bnodes.find(newindex);

    // // if new cell is structure
    // if (it != bnodes.end()) {
    //     if (it->second.getType() == BACKGROUND_STRUCTURE) {
    //         // check each direction
    //         for (int i = 0; i < ndm; ++i) {
    //             int diff = newindex[i] - index[0][i];
    //             if (diff == 0) continue;
    //             double out_disp = 0.0;
    //             if (diff > 0) {
    //                 out_disp = newpcrds[i] - crds[0][i] - bsize;
    //                 newpcrds[i] = crds[0][i] + bsize - out_disp;
    //                 pvel[i] = -pvel[i];
    //             } else {
    //                 out_disp = crds[0][i] - newpcrds[i];
    //                 newpcrds[i] = crds[0][i] + out_disp;
    //                 pvel[i] = -pvel[i];
    //             }
    //         }
    //     }
    // }

    // update particle
    if (pt->isUpdated() == false) {
        pt->setVel(pvel);
        pt->setAccel(pdvn);
        pt->setPressure(ppre);
        pt->setPdot(pdp);
    }

    // original->dest
    pt->moveTo(newpcrds, dt);

    return 0;
}

int BackgroundMesh::interpolate(const VVDouble& values, const VDouble& N,
                                VDouble& newvalue) {
    if (N.size() != values.size()) {
        opserr << "WARNING: sizes of shape function and nodal values "
                  "don't match\n";
        return -1;
    }
    if (N.empty()) {
        opserr << "WARNING: no shape functions\n";
        return -1;
    }
    if (values[0].empty()) {
        opserr << "WARNING: no nodal values\n";
        return -1;
    }

    VDouble temp(values[0].size());
    newvalue.assign(values[0].size(), 0.0);
    for (int i = 0; i < (int)N.size(); ++i) {
        if (values[i].size() != values[0].size()) {
            opserr << "WARNING: dimensions of nodal values are "
                      "different\n";
            newvalue.clear();
            return -1;
        }
        temp = values[i];
        temp *= N[i];
        newvalue += temp;
    }
    return 0;
}

int BackgroundMesh::interpolate(const VDouble& values, const VDouble& N,
                                double& newvalue) {
    if (N.size() != values.size()) {
        opserr << "WARNING: sizes of shape function and nodal values "
                  "don't match\n";
        return -1;
    }
    if (N.empty()) {
        opserr << "WARNING: no shape functions\n";
        return -1;
    }

    newvalue = 0.0;
    for (int i = 0; i < (int)N.size(); ++i) {
        newvalue += values[i] * N[i];
    }
    return 0;
}

int BackgroundMesh::createContact(const VInt& ndtags, const VInt& sids,
                                  VInt& elends) {
    // check inputs
    int ndm = OPS_GetNDM();
    if (ndtags.size() != sids.size()) {
        return 1;
    }
    if (ndm == 2) {
        if (ndtags.size() != 3) {
            opserr << "WARNING: 2D contact needs 3 nodes\n";
            return -1;
        }
    } else if (ndm == 3) {
        if (ndtags.size() != 4) {
            opserr << "WARNING: 3D contact needs 4 nodes\n";
            return -1;
        }
    }

    // get groups
    std::map<int, VInt> grp;
    for (int i = 0; i < (int)sids.size(); ++i) {
        grp[sids[i]].push_back(ndtags[i]);
    }

    if (grp.size() == 1) {
        // from same structure
        return 1;
    }

    // get secondary node
    int secondary = 0;
    int id = 0;
    bool find = false;
    for (std::map<int, VInt>::iterator it = grp.begin(); it != grp.end();
         ++it) {
        VInt& nds = it->second;
        if (nds.size() == 1) {
            // secondary node with largest sid
            if (!find || (id < it->first)) {
                id = it->first;
                secondary = nds[0];
                find = true;
            }
        } else if (find && id < it->first) {
            // if primary nodes have larger sid
            find = false;
        }
    }
    if (!find) return 1;

    // index for secondary node
    int index = 0;
    for (int i = 0; i < (int)ndtags.size(); ++i) {
        if (ndtags[i] == secondary) {
            index = i + 1;
            if (index >= (int)ndtags.size()) {
                index -= ndtags.size();
            }
            break;
        }
    }

    // get primary nodes
    elends.clear();
    for (int i = 0; i < (int)ndtags.size() - 1; ++i) {
        elends.push_back(ndtags[index]);
        index += 1;
        if (index >= (int)ndtags.size()) {
            index -= ndtags.size();
        }
    }
    elends.push_back(secondary);

    return 0;
}

void BackgroundMesh::setContactData(const VDouble& data) {
    contactData = data;
}

void BackgroundMesh::splitPrism(const VInt& prism, VVInt& tets, bool incl1,
                                bool incl4) {
    if (prism.size() != 6) return;
    tets.clear();
    VInt tet(4);

    // first tetehedron
    if (incl4) {
        tet[0] = prism[3];
        tet[1] = prism[4];
        tet[2] = prism[2];
        tet[3] = prism[5];
        tets.push_back(tet);
    }

    // second tetehedron
    if (incl1 && incl4) {
        tet[0] = prism[3];
        tet[1] = prism[1];
        tet[2] = prism[2];
        tet[3] = prism[4];
        tets.push_back(tet);
    }

    // third tetehedron
    if (incl1) {
        tet[0] = prism[3];
        tet[1] = prism[0];
        tet[2] = prism[2];
        tet[3] = prism[1];
        tets.push_back(tet);
    }
}

bool BackgroundMesh::check_area(const VDouble& ndcrds1,
                                const VDouble& ndcrds2,
                                const VDouble& ndcrds3) {
    VDouble coeff;
    bool zerovol = false;
    if (preNForTri(ndcrds1[0], ndcrds1[1], ndcrds2[0], ndcrds2[1],
                   ndcrds3[0], ndcrds3[1], coeff) < 0) {
        zerovol = true;
    }
    return zerovol;
}

bool BackgroundMesh::check_vol(const VDouble& ndcrds1,
                               const VDouble& ndcrds2,
                               const VDouble& ndcrds3,
                               const VDouble& ndcrds4) {
    VVDouble coeff;
    bool zerovol = false;
    if (preNForTet(ndcrds1, ndcrds2, ndcrds3, ndcrds4, coeff) < 0) {
        zerovol = true;
    }

    return zerovol;
}
