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
// Description: This class defines the Particle class
//

#ifndef Particle_h
#define Particle_h

#include <Vector.h>

class ParticleGroup;

class Particle
{
public:
    Particle():coord(),velocity(),pressure(0),group(0) {}
    ~Particle() {}

    void moveTo(const Vector& coord) {
	this->coord = coord;
    }

    void move(const Vector& disp) {
	this->coord += disp;
    }

    void setVel(const Vector& vel) {
	this->velocity = vel;
    }

    void setPressure(double p) {
	pressure = p;
    }

    void setGroup(ParticleGroup* g) {group = g;}

    void print();

    const Vector& getCrds() const {return coord;}
    const Vector& getVel() const {return velocity;}
    double getPressure() const {return pressure;}
    ParticleGroup* getGroup() {return group;}
    
private:
    Vector coord;
    Vector velocity;
    double pressure;
    ParticleGroup* group;
};


#endif
