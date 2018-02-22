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


#include <vector>
#include "Particle.h"
#include <string.h>

class ParticleGroup
{
public:
    ParticleGroup();
    ~ParticleGroup();

    void addParticle(const Vector& coord, const Vector& vel, double p);
    int numParticles() const {return (int)particles.size();}
    Particle* getParticle(int i) {return (i>=0&&i<this->numParticles())? particles[i]:0;}


    void setType(const char* eletype);
    const char* getType() const {return type;}

    void setProp(const Vector& p) {prop = p;}
    const Vector& getProp() const {return prop;}

    int getNDF() const {return ndf;}

    // printing
    void print();

    // a line of particles
    int point(const Vector& p1, const Vector& vel0, double p0);
    int line(const Vector& p1, const Vector& p2, int num,
	     const Vector& vel0, double p0);
    int qua_d(const Vector& p1, const Vector& p2, const Vector& p3,
	     const Vector& p4, int m, int n, const Vector& vel0, double p0);
    int tri(const Vector& p1, const Vector& p2,
	    const Vector& p4, int m, int n, const Vector& vel0, double p0);
    
private:
    
    std::vector<Particle*> particles;
    char* type;
    Vector prop;
    int ndf;
};

#endif
