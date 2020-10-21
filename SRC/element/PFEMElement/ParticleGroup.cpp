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
// $Date: 2016-1-27  $

// Written: Minjie Zhu
//
// Description: This class defines the ParticleGroup class
//

#include "ParticleGroup.h"
#include <elementAPI.h>
#include <fstream>
#include <string>
#include <sstream>

int OPS_ParticleGroup() {
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING: tag? type? pts?\n";
        return -1;
    }

    // get mesh tag
    int tag;
    int num = 1;
    if (OPS_GetIntInput(&num, &tag) < 0) {
        opserr << "WARNING: failed to get mesh tag\n";
        return -1;
    }

    // get particle group
    ParticleGroup *group = dynamic_cast<ParticleGroup *>(OPS_getMesh(tag));
    if (group == 0) {
        group = new ParticleGroup(tag);
        if (group == 0) {
            opserr << "WARNING: failed to create/get new particle group\n";
            return -1;
        }
        if (OPS_addMesh(group) == false) {
            opserr << "WARNING: failed to add particle group\n";
            return -1;
        }
    }

    // get geometry
    int ndm = OPS_GetNDM();
    VDouble p1(ndm), p2(ndm), p3(ndm), p4(ndm);
    VDouble p5(ndm), p6(ndm), p7(ndm), p8(ndm);
    VDouble pointdata;
    VInt nump(3);
    const char *geotype = OPS_GetString();
    if (strcmp(geotype, "quad") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 4 * ndm + 2) {
            opserr << "WARNING: insufficient args\n";
            return -1;
        }

        // node coord
        if (OPS_GetDoubleInput(&ndm, &p1[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for first point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p2[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for second point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p3[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for third point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p4[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for fouth point\n";
            return -1;
        }

        // num of particles
        int numdata = 2;
        if (OPS_GetIntInput(&numdata, &nump[0]) < 0) {
            opserr << "WARNING: failed to get particle mesh size\n";
            return -1;
        }
    } else if (strcmp(geotype, "cube") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 8 * ndm + 3) {
            opserr << "WARNING: insufficient args -- ";
            opserr << "bottom p1, p2, p3, p4 and top p5, p6, p7, p8\n";
            return -1;
        }

        // node coord
        if (OPS_GetDoubleInput(&ndm, &p1[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for first bottom point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p2[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for second bottom point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p3[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for third bottom point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p4[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for fouth bottom point\n";
            return -1;
        }

        // node coord
        if (OPS_GetDoubleInput(&ndm, &p5[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for first top point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p6[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for second top point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p7[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for third top point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p8[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for fouth top point\n";
            return -1;
        }

        // num of particles
        int numdata = 3;
        nump.resize(numdata);
        if (OPS_GetIntInput(&numdata, &nump[0]) < 0) {
            opserr << "WARNING: failed to get particle mesh size\n";
            return -1;
        }

    } else if (strcmp(geotype, "tri") == 0) {

        if (OPS_GetNumRemainingInputArgs() < 3 * ndm + 2) {
            opserr << "WARNING: insufficient args\n";
            return -1;
        }

        // node coord
        if (OPS_GetDoubleInput(&ndm, &p1[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for first point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p2[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for second point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p3[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for third point\n";
            return -1;
        }

        // num of particles
        int numdata = 2;
        if (OPS_GetIntInput(&numdata, &nump[0]) < 0) {
            opserr << "WARNING: failed to get particle mesh size\n";
            return -1;
        }
    } else if (strcmp(geotype, "line") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 2 * ndm + 1) {
            opserr << "WARNING: insufficient args\n";
            return -1;
        }

        // node coord
        if (OPS_GetDoubleInput(&ndm, &p1[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for first point\n";
            return -1;
        }
        if (OPS_GetDoubleInput(&ndm, &p2[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for second point\n";
            return -1;
        }

        // num of particles
        int numdata = 1;
        if (OPS_GetIntInput(&numdata, &nump[0]) < 0) {
            opserr << "WARNING: failed to get number of particles\n";
            return -1;
        }
    } else if (strcmp(geotype, "pointlist") == 0) {
        int numdata = OPS_GetNumRemainingInputArgs();
        if (numdata < 1) {
            group->pointlist(pointdata);
            numdata = (int) pointdata.size(); 
            if (OPS_SetDoubleOutput(&numdata, &pointdata[0], false) < 0) {
                opserr << "WARNING: failed to set output\n";
                return -1;
            }
            return 0;
        }

        // number of points
        int num_point = 0;
        numdata = 1;
        if (OPS_GetIntInput(&numdata, &num_point) < 0) {
            opserr << "WARNING: failed to get number of points\n";
            return -1;
        }

        // check input
        numdata = num_point * (4 * ndm + 1);
        if (OPS_GetNumRemainingInputArgs() < numdata) {
            opserr << "WARNING: insufficient input for " 
            << num_point << " points: [x1n, y1, <z1n>, x1, y1, <z1> "
            << "vx1, vy1, <vz1>, ax1, ay1, <az1>, p1, x2n, ...]\n";
            return -1;
        }

        // node coord
        pointdata.resize(numdata);
        if (OPS_GetDoubleInput(&numdata, &pointdata[0]) < 0) {
            opserr << "WARNING: failed to get cooridnates for points\n";
            return -1;
        }

    } else {
        opserr << "WARNING: unknown geometry type\n";
        return -1;
    }

    // set ele args
    if (group->setEleArgs() < 0) {
        opserr << "WARNING: failed to set element arguments\n";
        return -1;
    }

    // intial velocity and pressure
    VDouble vel0(ndm);
    double p0 = 0;
    const char *eletype = 0;

    // type
    if (OPS_GetNumRemainingInputArgs() > 0) {
        eletype = OPS_GetString();

        if (strcmp(eletype, "-vel") == 0) {
            // initial velocity
            if (OPS_GetNumRemainingInputArgs() < ndm) {
                opserr << "WARNING: -vel insufficient args\n";
                return -1;
            }

            // node vel
            if (OPS_GetDoubleInput(&ndm, &vel0[0]) < 0) {
                opserr << "WARNING: failed to get initial velocity\n";
                return -1;
            }

        } else if (strcmp(eletype, "-pressure") == 0) {
            // initial pressure
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: -pressure insufficient args\n";
                return -1;
            }

            // node pressure
            int numdata = 1;
            if (OPS_GetDoubleInput(&numdata, &p0) < 0) {
                opserr << "WARNING: failed to get initial pressure\n";
                return -1;
            }
        }
    }

    // generate particles
    if (strcmp(geotype, "quad") == 0) {
        group->qua_d(p1, p2, p3, p4, nump[0], nump[1], vel0, p0);
    } else if (strcmp(geotype, "cube") == 0) {
        VVDouble pts(8);
        pts[0] = p1;
        pts[1] = p2;
        pts[2] = p3;
        pts[3] = p4;
        pts[4] = p5;
        pts[5] = p6;
        pts[6] = p7;
        pts[7] = p8;
        group->cube(pts, nump, vel0, p0);
    } else if (strcmp(geotype, "tri") == 0) {
        group->tri(p1, p2, p3, nump[0], nump[1], vel0, p0);
    } else if (strcmp(geotype, "line") == 0) {
        group->line(p1, p2, nump[0], vel0, p0);
    } else if (strcmp(geotype, "pointlist") == 0) {
        group->pointlist(pointdata, ndm);
    }

    return 0;
}

ParticleGroup::ParticleGroup(int tag)
        : Mesh(tag, 1), particles() {
}

ParticleGroup::~ParticleGroup() {
    for (int i = 0; i < (int) particles.size(); i++) {
        Particle *p = particles[i];
        if (p != 0)
            delete p;
    }
    particles.clear();
}

void ParticleGroup::addParticle(const VDouble &coord, const VDouble &vel, double p) {
    Particle *particle = new Particle;
    particles.push_back(particle);

    particle->moveTo(coord, 0.0);
    particle->setVel(vel);
    particle->setPressure(p);
    VDouble accel = vel;
    accel *= 0.0;
    particle->setAccel(accel);

    particle->setGroupTag(this->getTag());
}

void ParticleGroup::addParticle(const VDouble &coordn,
                                const VDouble &coord,
                                const VDouble &vel,
                                const VDouble &accel,
                                double p) {
    Particle *particle = new Particle;
    particles.push_back(particle);

    particle->moveTo(coordn, 0.0);
    particle->setVel(vel);
    particle->moveTo(coord, 0.0);
    particle->setPressure(p);
    particle->setAccel(accel);
    particle->setGroupTag(this->getTag());
}

void ParticleGroup::removeParticles(const VInt &rm) {
    if (rm.size() != particles.size()) return;
    VParticle newp;
    for (unsigned int i = 0; i < particles.size(); ++i) {
        Particle *p = particles[i];
        if (p != 0 && rm[i] == 0) {
            newp.push_back(p);
        } else if (p != 0) {
            delete p;
        }
    }
    particles = newp;
}

int ParticleGroup::line(const VDouble &p1, const VDouble &p2, int num,
                        const VDouble &vel0, double p0) {
    if (num <= 0)
        return 0;

    if (p1.size() != p2.size())
        return -1;

    VDouble p1p2 = p2;
    p1p2 -= p1;
    p1p2 /= num;

    VDouble crds(p1);
    VDouble vel(crds.size());
    for (int i = 0; i < (int) vel.size(); i++) {
        if (i < (int) vel0.size()) {
            vel[i] = vel0[i];
        }
    }

    for (int i = 0; i <= num; i++) {
        this->addParticle(crds, vel, p0);
        crds += p1p2;
    }

    return 0;
}

int ParticleGroup::pointlist(VDouble &pointdata) {
    int ndm = OPS_GetNDM();
    pointdata.clear();
    pointdata.reserve(particles.size() * (4 * ndm + 1));
    for (auto particle : particles) {
        auto tag = particle->getTag();
        const auto& crdsn = particle->getCrdsn(); 
        const auto& crds = particle->getCrds(); 
        const auto& vel = particle->getVel(); 
        const auto& accel = particle->getAccel(); 
        double p = particle->getPressure(); 
        pointdata.push_back(tag);
        for (int j = 0; j < ndm; ++j) {
            pointdata.push_back(crdsn[j]);
        }
        for (int j = 0; j < ndm; ++j) {
            pointdata.push_back(crds[j]);
        }
        for (int j = 0; j < ndm; ++j) {
            pointdata.push_back(vel[j]);
        }
        for (int j = 0; j < ndm; ++j) {
            pointdata.push_back(accel[j]);
        }
        pointdata.push_back(p);
    }

    return 0;
}

int ParticleGroup::pointlist(const VDouble &pointdata, int ndm) {
    VDouble crdsn(ndm);
    VDouble crds(ndm);
    VDouble vel(ndm);
    VDouble accel(ndm);
    double p0 = 0.0;
    for (int i = 0; i < (int)pointdata.size(); i += 4 * ndm + 1) {
        for (int j = 0; j < ndm; ++j) {
            crdsn[j] = pointdata[i + j];
        }
        for (int j = 0; j < ndm; ++j) {
            crds[j] = pointdata[i + ndm + j];
        }
        for (int j = 0; j < ndm; ++j) {
            vel[j] = pointdata[i + 2 * ndm + j];
        }
        for (int j = 0; j < ndm; ++j) {
            accel[j] = pointdata[i + 3 * ndm + j];
        }
        p0 = pointdata[i + 4 * ndm];
        this->addParticle(crdsn, crds, vel, accel, p0);
    }

    return 0;
}

int ParticleGroup::qua_d(const VDouble &p1, const VDouble &p2, const VDouble &p3,
                         const VDouble &p4, int m, int n, const VDouble &vel0, double p0) {

    if (m <= 0 || n <= 0)
        return 0;
    if (p1.size() != p2.size())
        return -1;
    if (p3.size() != p4.size())
        return -1;
    if (p1.size() != p4.size())
        return -1;

    // line 12
    VDouble p1p2 = p2;
    p1p2 -= p1;
    p1p2 /= m;

    // line 43
    VDouble p4p3 = p3;
    p4p3 -= p4;
    p4p3 /= m;

    // each line  between 12 and 43
    VDouble crds12 = p1p2;
    crds12 /= 2.0;
    crds12 += p1;

    VDouble crds43 = p4p3;
    crds43 /= 2.0;
    crds43 += p4;

    for (int i = 1; i <= m; i++) {

        // line 12 to 43
        VDouble p1243 = crds43;
        p1243 -= crds12;
        p1243 /= 2 * n;

        VDouble plow = crds12;
        plow += p1243;
        VDouble phigh = crds43;
        phigh -= p1243;
        if (this->line(plow, phigh, n - 1, vel0, p0) < 0) {
            return -1;
        }

        // incr
        crds12 += p1p2;
        crds43 += p4p3;
    }

    return 0;
}

int ParticleGroup::tri(const VDouble &p1, const VDouble &p2, const VDouble &p3,
                       int m, int n, const VDouble &vel0, double p0) {

    if (m <= 0 || n <= 0)
        return 0;
    if (p1.size() != p2.size())
        return -1;
    if (p3.size() != p1.size())
        return -1;

    // the mesh size along edge 1-2 and 1-3
    double h1 = 1.0 / m;
    double h2 = 1.0 / n;

    // initial vel and pressure
    VDouble vel(p1.size());
    for (int i = 0; i < (int) vel.size(); i++) {
        if (i < (int) vel0.size()) {
            vel[i] = vel0[i];
        }
    }

    // using area coordinates to generate particles' coordinates
    VDouble crds;
    VDouble temp;
    for (int i = 0; i < m; i++) {
        double L1 = (i + 0.5) * h1;
        for (int j = 0; j < n; j++) {
            double L2 = (j + 0.5) * h2;
            double L3 = 1.0 - L1 - L2;
            if (L3 < -1e-6)
                continue;

            crds = p1;
            crds *= L1;

            temp = p2;
            temp *= L2;
            crds += temp;

            temp = p3;
            temp *= L3;
            crds += temp;

            this->addParticle(crds, vel, p0);
        }
    }

    return 0;
}

int
ParticleGroup::cube(const VVDouble &pts, const VInt &num,
                    const VDouble &vel0, double p0) {
    if (pts.size() != 8) {
        opserr << "WARNING: pts.size() != 8 -- ParticleGroup::cube\n";
        return -1;
    }
    if (num.size() != 3) {
        opserr << "WARNING: num.size() != 3 -- ParticleGroup::cube\n";
        return -1;
    }

    // line 15,26,37,48
    VVDouble dirs(4), crds(4);
    for (unsigned int i = 0; i < dirs.size(); ++i) {
        dirs[i] = pts[i + 4];
        dirs[i] -= pts[i];
        dirs[i] /= num[2];
        crds[i] = dirs[i];
        crds[i] /= 2.0;
        crds[i] += pts[i];
    }

    // create each plane of particles
    for (int i = 0; i < num[2]; ++i) {
        if (qua_d(crds[0], crds[1], crds[2], crds[3], num[0], num[1], vel0, p0) < 0) {
            opserr << "WARNING: failed to create particles -- ParticleGroup::cube\n";
            return -1;
        }

        for (unsigned int j = 0; j < crds.size(); ++j) {
            crds[j] += dirs[j];
        }
    }

    return 0;
}
