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

#ifndef ParticleGroup_h
#define ParticleGroup_h


#include "BackgroundDef.h"
#include "Particle.h"
#include "Mesh.h"

class ParticleGroup : public Mesh
{
public:
    explicit ParticleGroup(int tag);
    ~ParticleGroup();

    // particles
    void addParticle(const VDouble& coord, const VDouble& vel, double p);
    void addParticle(const VDouble& coordn, const VDouble& coord,
                     const VDouble& vel,
                     const VDouble& accel, double p);
    void removeParticles(const VInt& rm);
    int numParticles() const {return (int)particles.size();}
    Particle* getParticle(int i) {
        return (i>=0&&i<numParticles())? particles[i]:0;
    }

    // dummy mesh
    int mesh(){return 0;}

    // return particles
    int pointlist(VDouble& pointdata);

    // create particles
    int pointlist(const VDouble& pointdata, int ndm);
    int line(const VDouble& p1, const VDouble& p2, int num,
	     const VDouble& vel0, double p0);
    int qua_d(const VDouble& p1, const VDouble& p2, const VDouble& p3,
	     const VDouble& p4, int m, int n, const VDouble& vel0, double p0);
    int tri(const VDouble& p1, const VDouble& p2,
	    const VDouble& p4, int m, int n, const VDouble& vel0, double p0);
    int cube(const VVDouble& pts, const VInt& num, const VDouble& vel0, double p0);

private:

    VParticle particles;
};

#endif
