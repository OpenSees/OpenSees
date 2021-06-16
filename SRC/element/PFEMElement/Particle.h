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
            updated = true;
        }
    }

    void setVelOnly(const VDouble &vel) { this->velocity = vel; }

    void setPressure(double p) { pressure = p; }

    void setAccel(const VDouble &accel) { this->accel = accel; }

    void setPdot(double pdot) { this->pdot = pdot; }

    void setGroupTag(int tag) { this->gtag = tag; }

    void needUpdate(double dt) {
        updated = false;
        this->dt = dt;
    }

    const VDouble &getCrds() const { return coord; }

    const VDouble &getVel() const { return velocity; }

    const VDouble &getAccel() const { return accel; }

    double getPressure() const { return pressure; }

    double getPdot() const { return pdot; }

    int getGroupTag() const { return gtag; }

    bool isUpdated() const { return updated; }

    double getDt() const { return dt; }

    size_t getTag() const { return tag; }

   private:
    VDouble coord;
    VDouble velocity;
    VDouble accel;
    double pressure, pdot;
    int gtag;
    bool updated;
    double dt;
    size_t tag;

    static size_t curr_tag;
};

#endif
