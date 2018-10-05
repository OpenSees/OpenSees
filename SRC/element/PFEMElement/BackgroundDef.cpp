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

#include "BackgroundDef.h"
#include <cmath>
#include <elementAPI.h>

////// Overload functions ////////////////
const VDouble& operator+=(VDouble& v1, const VDouble& v2)
{
    if (v1.size() > v2.size()) return v1;
    for (unsigned int i=0; i<v1.size(); i++) {
	v1[i] += v2[i];
    }
    return v1;
}

const VDouble& operator-=(VDouble& v1, const VDouble& v2)
{
    if (v1.size() > v2.size()) return v1;
    for (unsigned int i=0; i<v1.size(); i++) {
	v1[i] -= v2[i];
    }
    return v1;
}

const VDouble& operator+=(VDouble& v1, double val)
{
    for (unsigned int i=0; i<v1.size(); i++) {
	v1[i] += val;
    }
    return v1;
}

const VDouble& operator-=(VDouble& v1, double val)
{
    for (unsigned int i=0; i<v1.size(); i++) {
	v1[i] -= val;
    }
    return v1;
}

const VDouble& operator*=(VDouble& v1, double val)
{
    for (unsigned int i=0; i<v1.size(); i++) {
	v1[i] *= val;
    }
    return v1;
}

const VDouble& operator/=(VDouble& v1, double val)
{
    for (unsigned int i=0; i<v1.size(); i++) {
	v1[i] /= val;
    }
    return v1;
}

const VInt& operator+=(VInt& v1, int val)
{
    for (unsigned int i=0; i<v1.size(); i++) {
	v1[i] += val;
    }
    return v1;
}

const VInt& operator-=(VInt& v1, int val)
{
    for (unsigned int i=0; i<v1.size(); i++) {
	v1[i] -= val;
    }
    return v1;
}

void toVDouble(const Vector& vec, VDouble& res)
{
    res.resize(vec.Size());

    for (int i=0; i<vec.Size(); i++) {
	res[i] = vec(i);
    }
}

void toVector(const VDouble& v, Vector& res)
{
    res.resize(v.size());

    for (unsigned int i=0; i<v.size(); i++) {
	res(i) = v[i];
    }

}

double normVDouble(const VDouble& v)
{
    double res = 0.0;
    for (unsigned int i=0; i<v.size(); i++) {
	res += v[i]*v[i];
    }
    return sqrt(res);
}

std::ostream& operator<<(std::ostream& os, const VDouble& v)
{
    for (unsigned int i=0; i<v.size(); i++) {
	os << v[i] << " ";
    }
    os << "\n";

    return os;
}

std::ostream& operator<<(std::ostream& os, const VInt& v)
{
    for (VInt::size_type i=0; i<v.size(); i++) {
	os << v[i] << " ";
    }
    os << "\n";

    return os;
}

double dotVDouble(const VDouble& v1, const VDouble& v2)
{
    if (v1.size() != v2.size()) {
	return 0.0;
    }
    double res = 0.0;
    for (unsigned int i=0; i<v1.size(); ++i) {
	res += v1[i]*v2[i];
    }
    return res;
}

void crossVDouble(const VDouble& v1, const VDouble& v2, VDouble& res)
{
    res.resize(3, 0.0);
    if (v1.size()==2 && v2.size()==2) {
	res[2] = v1[0]*v2[1]-v1[1]*v2[0];
    } else if (v1.size()==3 && v2.size()==3) {
	res[0] = v1[1]*v2[2]-v1[2]*v2[1];
	res[1] = v1[2]*v2[0]-v1[0]*v2[2];
	res[2] = v1[0]*v2[1]-v1[1]*v2[0];
    }
}
