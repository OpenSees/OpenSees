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

ParticleGroup::ParticleGroup()
    :particles(), type(0), prop(), ndf(0)
{

}

ParticleGroup::~ParticleGroup()
{
    for (int i=0; i<(int)particles.size(); i++) {
	Particle* p = particles[i];
	if (p != 0) delete p;
    }
    particles.clear();

    if (type != 0) delete [] type;
}

void
ParticleGroup::addParticle(const Vector& coord, const Vector& vel, double p)
{
    Particle* particle = new Particle;
    particles.push_back(particle);

    particle->moveTo(coord);
    particle->setVel(vel);
    particle->setPressure(p);
    particle->setGroup(this);
}

void
ParticleGroup::setType(const char* eletype)
{
    if (type != 0) {
	delete [] type;
    }

    type = new char[strlen(eletype)+1];
    strcpy(type, eletype);
    
    // set velocity ndf
    if (strcmp(type,"PFEMElement2D") == 0
	||strcmp(type,"PFEMElement2DBubble") == 0
	||strcmp(type,"PFEMElement2DCompressible") == 0) {
	ndf = 2;
    } else {
	opserr<<"WARNING: particle type "<<type<<" is unknown\n";
    }
}

void
ParticleGroup::print() {
    for(int i=0; i<this->numParticles(); i++) {
	Particle* p = particles[i];
	if (p != 0) {
	    p->print();
	}
    }
}

int
ParticleGroup::line(const Vector& p1, const Vector& p2, int num,
		    const Vector& vel0, double p0)
{
    if(num <= 0) return 0;
    
    if(p1.Size() != p2.Size()) return -1;
    
    if(ndf == 0) {
	opserr<<"WARNING: type of particle has not been set or unkown\n";
	return -1;
    }
    
    Vector p1p2 = p2-p1;
    p1p2 /= num;

    Vector crds(p1);
    Vector vel(ndf);
    for(int i=0; i<vel.Size(); i++) {
	if (i < vel0.Size()) {
	    vel(i) = vel0(i);
	}
    }
    for(int i=0; i<=num; i++) {
	this->addParticle(crds, vel, p0);
	crds += p1p2;
    }
    
    return 0;
}

int
ParticleGroup::point(const Vector& p1, const Vector& vel0, double p0)
{
    if(ndf == 0) {
	opserr<<"WARNING: type of particle has not been set\n";
	return -1;
    }

    Vector crds(p1);
    Vector vel(ndf);
    for(int i=0; i<vel.Size(); i++) {
	if (i < vel0.Size()) {
	    vel(i) = vel0(i);
	}
    }
	
    this->addParticle(crds, vel, p0);
    
    return 0;
}

int
ParticleGroup::qua_d(const Vector& p1, const Vector& p2, const Vector& p3,
		    const Vector& p4, int m, int n, const Vector& vel0, double p0)
{

    if (m<=0 || n<=0) return 0;
    if(p1.Size() != p2.Size()) return -1;
    if(p3.Size() != p4.Size()) return -1;
    if(p1.Size() != p4.Size()) return -1;

    // line 12
    Vector p1p2 = p2-p1;
    p1p2 /= m;

    // line 43
    Vector p4p3 = p3-p4;
    p4p3 /= m;

    // each line  between 12 and 43
    Vector crds12 = p1+p1p2/2.0;
    Vector crds43 = p4+p4p3/2.0;
    for(int i=1; i<=m; i++) {

	// line 12 to 43
	Vector p1243 = crds43-crds12;
	p1243 /= 2*n;
	if (this->line(crds12+p1243, crds43-p1243, n-1, vel0, p0) < 0) {
	    return -1;
	}

	// incr
	crds12 += p1p2;
	crds43 += p4p3;
    }
	
    return 0;
}

int
ParticleGroup::tri(const Vector& p1, const Vector& p2, const Vector& p3,
		   int m, int n, const Vector& vel0, double p0)
{

    if (m<=0 || n<=0) return 0;
    if(p1.Size() != p2.Size()) return -1;
    if(p3.Size() != p1.Size()) return -1;

    // the mesh size along edge 1-2 and 1-3
    double h1 = 1.0/m;
    double h2 = 1.0/n;

    // initial vel and pressure
    Vector crds(p1);
    Vector vel(ndf);
    for(int i=0; i<vel.Size(); i++) {
	if (i < vel0.Size()) {
	    vel(i) = vel0(i);
	}
    }

    // using area coordinates to generate particles' coordinates
    for(int i=0; i<m; i++) {
	double L1 = (i+0.5)*h1;
	for (int j=0; j<n; j++) {
	    double L2 = (j+0.5)*h2;
	    double L3 = 1.0-L1-L2;
	    if (L3 < -1e-6) continue;

	    crds.Zero();
	    crds.addVector(0.0, p1, L1);
	    crds.addVector(1.0, p2, L2);
	    crds.addVector(1.0, p3, L3);
	    
	    this->addParticle(crds, vel, p0);
	}
    }

    return 0;
}
