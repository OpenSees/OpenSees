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

#include "BackgroundDef.h"

class Particle {

public:
    Particle();

    ~Particle();

    void moveTo(const VDouble &coord, double subdt) {
        this->coord = coord;
        dt -= subdt;
        if (dt < 0) dt = 0.0;
    }

    void move(const VDouble &disp, double subdt) {
        this->coord += disp;
        dt -= subdt;
        if (dt < 0) dt = 0.0;
    }

    void setVel(const VDouble &vel) {
        if (!updated) {
            this->velocity = vel;
            this->coordn = this->coord;
            updated = true;
        }
    }

    void incrVel(const VDouble &dv) {
        if (!updated) {
            this->velocity += dv;
            this->coordn = this->coord;
            updated = true;
        }
    }

    void setPressure(double p) {
        pressure = p;
    }

    void setAccel(const VDouble &accel) {
        this->accel = accel;
    }

    void setPdot(double pdot) {
        this->pdot = pdot;
    }

    void setGroupTag(int tag) {
        this->gtag = tag;
    }

    void needUpdate(double dt) {
        updated = false;
        this->dt = dt;
    }

    void setFixed() {
        fixed = true;
    }

    const VDouble &getCrds() const { return coord; }

    const VDouble &getCrdsn() const { return coordn; }

    const VDouble &getVel() const { return velocity; }

    const VDouble &getAccel() const { return accel; }

    double getPressure() const { return pressure; }

    double getPdot() const { return pdot; }

    int getGroupTag() const { return gtag; }

    bool isUpdated() const { return updated; }

    bool isFixed() const { return fixed; }

    double getDt() const { return dt; }

    size_t getTag() const {return tag;}

private:
    VDouble coord, coordn;
    VDouble velocity;
    VDouble accel;
    double pressure, pdot;
    int gtag;
    bool updated;
    double dt;
    bool fixed;
    size_t tag;

    static size_t curr_tag;
};


#endif
