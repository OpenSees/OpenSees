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
// Description: type definitions for background mesh method
//

#ifndef BackgroundDef_h
#define BackgroundDef_h

#include <vector>
#include <set>
#include <map>
#include <Vector.h>
#include <math.h>
#include <iostream>

class Particle;
class ParticleGroup;
class Node;

// type defs
typedef std::vector<int> VInt;
typedef std::vector<bool> VBool;
typedef std::vector<VInt> VVInt;
typedef std::vector<double> VDouble;
typedef std::vector<VDouble> VVDouble;
typedef std::vector<Particle*> VParticle;
typedef std::vector<VParticle> VVParticle;

// functions for VDouble
const VDouble& operator+=(VDouble& v1, const VDouble& v2);
const VDouble& operator-=(VDouble& v1, const VDouble& v2);
const VDouble& operator+=(VDouble& v1, double val);
const VDouble& operator-=(VDouble& v1, double val);
const VDouble& operator*=(VDouble& v1, double val);
const VDouble& operator/=(VDouble& v1, double val);
const VInt& operator+=(VInt& v1, int val);
const VInt& operator-=(VInt& v1, int val);
std::ostream& operator<<(std::ostream& os, const VDouble& v);
std::ostream& operator<<(std::ostream& os, const VInt& v);
void toVDouble(const Vector& vec, VDouble& res);
void toVector(const VDouble& v, Vector& res);
double normVDouble(const VDouble& v);
double dotVDouble(const VDouble& v1, const VDouble& v2);
void crossVDouble(const VDouble& v1, const VDouble& v2, VDouble& res);

#endif
