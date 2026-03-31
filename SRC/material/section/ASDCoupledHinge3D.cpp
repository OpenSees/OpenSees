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
                                                                        
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ASDCoupledHinge3D.cpp,v $
                                                                        
                                                                        
// File: ~/section/ASDCoupledHinge3D.C
//
// Written: Diego Talledo
// Created: April 2021
// Revision: A
//
// Purpose: This file contains the implementation for the ASDCoupledHinge3D class. 

//#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <classTags.h>
#include <ASDCoupledHinge3D.h>
#include <MaterialResponse.h>
#include <ID.h>
#include <classTags.h>
#include <elementAPI.h>
#include <Pinching4Material.h>
#include <ElasticMaterial.h>
#include <Message.h>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <vector>
#include <array>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//#define _DBG_COUPLEDSEC3D

ASDCoupledHinge3DDomainData::ASDCoupledHinge3DDomainData(int nN, int nTheta, int nData) :
    size(nN * nTheta * nData), numberAxial(nN), numberTheta(nTheta), numberData(nData), theVector(size)
{
}

double ASDCoupledHinge3DDomainData::getValue(int i, int j, int k)
{
    int idx = i * (numberTheta * numberData) + k * numberTheta + j;
    return theVector(idx);
}

void ASDCoupledHinge3DDomainData::setValue(int i, int j, int k, double val)
{
    int idx = i * (numberTheta * numberData) + k * numberTheta + j;
    theVector(idx) = val;
}

void ASDCoupledHinge3DDomainData::getRangeN(double& Nmin, double& Nmax)
{
    Nmin = this->getValue(0, 0, 0);
    Nmax = this->getValue(numberAxial - 1, 0, 0);
}

void ASDCoupledHinge3DDomainData::print(void)
{
    opserr << "DomainData object with a total of " << size << " elements" << endln;
    opserr << "  nAxial = " << numberAxial << " - nTheta = " << numberTheta << " - nData = " << numberData << endln;
    for (int i = 0; i < numberAxial; i++) 
    {
        opserr << "i = " << i << " (N = " << this->getValue(i, 0, 0) << ")\n";
        for (int j = 0; j < numberTheta; j++) 
        {
            opserr << "(" << this->getValue(i,j,0);
            for (int k = 1; k < numberData; k++) 
            {
                opserr << ", " << this->getValue(i, j, k);
            }
            opserr << "), ";
        }
        opserr << endln;
    }
}

int ASDCoupledHinge3DDomainData::getMyMzForNAndDirection(double N, double theta, double& My, double& Mz) {
    // transform theta in the range [0, 2pi[
    while (theta > 2 * M_PI)
        theta -= 2 * M_PI;
    Vector direction(2);
    direction(0) = cos(theta);
    direction(1) = sin(theta);

    // Find N first
    int i = 0;
    bool found = false;
#ifdef _DBG_COUPLEDSEC3D
    opserr << "\n\ngetMyMzForNAndDirection: \n"
        "N = " << N << "\n"
        "theta = " << theta << "\n"
        "direction = [" << direction(0) << ", " << direction(1) << "]\n";
    opserr << "First check within the limits: " << this->getValue(0, 0, 0) << " and " << this->getValue(numberAxial - 1, 0, 0) << "\n";
#endif

    if ((N <= this->getValue(0, 0, 0)) || (N >= this->getValue(numberAxial - 1, 0, 0))) 
    {
        // N is greater than Nmax or lower than Nmin -> return 0
#ifdef _DBG_COUPLEDSEC3D
        opserr << "N = " << N << " is greater than " << this->getValue(numberAxial - 1, 0, 0) << " or lower than " << this->getValue(0, 0, 0) << "\n";
#endif
        My = 0;
        Mz = 0;
        return 0;
    }
    while ((!found) && (i < numberAxial)) 
    {
        if (this->getValue(i, 0, 0) > N) 
            found = true;
        i += 1;
    }

    if (found) 
    {
#ifdef _DBG_COUPLEDSEC3D
        opserr << "N = " << N << " found checking the domain\n";
#endif
    }
    else 
    {
        opserr << "Not possible to find N during interpolation. An errror occurred\n";
        return -1;
    }

    int i1 = i - 2;
    int i2 = i - 1;

    double N1 = this->getValue(i1, 0, 0);
    double N2 = this->getValue(i2, 0, 0);

#ifdef _DBG_COUPLEDSEC3D
    opserr << "N = " << N << " between " << i1 << " and " << i2 << " with N1 = " << N1 << " and N2 = " << N2 << "\n";
#endif

    //opserr << "found N between " << N1 << "(" << i1 << ") and " << N2 << "(" << i2 << ")\n";
    //opserr << "My-Mz domain for N1:\n";
    //for (int k = 0; k < numberTheta; k++) {
    //    opserr << "(" << this->getValue(i1, k, 1) << ", " << this->getValue(i1, k, 2) << ") ";
    //}
    //opserr << "\n";

    // Find My and Mz for given "theta"
    double My1_j, My2_j, My_j, Mz1_j, Mz2_j, Mz_j;
    double My1_jp1, My2_jp1, My_jp1, Mz1_jp1, Mz2_jp1, Mz_jp1;
    found = false;
    double theta_j = 0;
    double theta_jp1 = 0;
    Vector dir_j(2);
    Vector dir_jp1(2);
    int j = -1;
    int jp1 = 0;
    int inverse = 0;
    double jXjp1, jXdirection, directionXjp1;
    while ((!found) && (j < numberTheta)) 
    { 
        j += 1;
        jp1 += 1;
        if (jp1 >= numberTheta) jp1 = 0;
        //opserr << "j = " << j << " - jp1 = " << jp1 << "\n";
        My1_j = this->getValue(i1, j, 1);
        My2_j = this->getValue(i2, j, 1);
        My_j = (My2_j - My1_j) / (N2 - N1) * (N - N1) + My1_j;
        Mz1_j = this->getValue(i1, j, 2);
        Mz2_j = this->getValue(i2, j, 2);
        Mz_j = (Mz2_j - Mz1_j) / (N2 - N1) * (N - N1) + Mz1_j;
        My1_jp1 = this->getValue(i1, jp1, 1);
        My2_jp1 = this->getValue(i2, jp1, 1);
        My_jp1 = (My2_jp1 - My1_jp1) / (N2 - N1) * (N - N1) + My1_jp1;
        Mz1_jp1 = this->getValue(i1, jp1, 2);
        Mz2_jp1 = this->getValue(i2, jp1, 2);
        Mz_jp1 = (Mz2_jp1 - Mz1_jp1) / (N2 - N1) * (N - N1) + Mz1_jp1;
        dir_j(0) = My_j;
        dir_j(1) = Mz_j;
        dir_j.Normalize();
        theta_j = atan2(Mz_j, My_j);
        // transform theta in the range [0, 2pi[
        if (theta_j < 0)
            theta_j += 2 * M_PI;
        dir_jp1(0) = My_jp1;
        dir_jp1(1) = Mz_jp1;
        dir_jp1.Normalize();
        theta_jp1 = atan2(Mz_jp1, My_jp1);
        if (theta_jp1 < 0)
            theta_jp1 += 2 * M_PI;
#ifdef _DBG_COUPLEDSEC3D
        opserr << "j = " << j << " - jp1 = " << jp1 << "\n";
        opserr << "theta_j = " << theta_j <<  " - theta_jp1 = " << theta_jp1 << " - theta = " << theta << "\n";
        opserr << "direction j = [" << dir_j(0) << ", " << dir_j(1) << "]\n";
        opserr << "direction j+1 = [" << dir_jp1(0) << ", " << dir_jp1(1) << "]\n";

#endif
        jXjp1 = dir_j(0) * dir_jp1(1) - dir_j(1) * dir_jp1(0);
        jXdirection = dir_j(0) * direction(1) - dir_j(1) * direction(0);
        directionXjp1 = direction(0) * dir_jp1(1) - direction(1) * dir_jp1(0);
        if ((jXjp1 * jXdirection >= 0) && (directionXjp1 * jXjp1 >= 0))
        {
#ifdef _DBG_COUPLEDSEC3D
            opserr << "wanted direction found!";
            opserr << " j = " << j << " - jp1 = " << jp1 << endln;
#endif
            found = true;
            if (jXjp1 > 0)
                inverse = 1;
        }
        //if ((theta >= theta_j) && (theta < theta_jp1)) 
        //{
        //    //opserr << "theta_j = " << theta_j << " - theta = " << theta << " - theta_jp1 = " << theta_jp1 << "\n";
        //    found = true;
        //}
        //else if ((theta <= theta_j) && (theta > theta_jp1)) 
        //{
        //    //opserr << "theta_jp1 = " << theta_jp1 << " - theta = " << theta << " - theta_j = " << theta_j << "\n";
        //    found = true;
        //    inverse = 1;
        //}
    }

    if (found) 
    {
        //opserr << "found theta\n";
        //opserr << "between: " << theta_j << "(" << j << ") and " << theta_jp1 << "(" << jp1 << ") - inverse: " << inverse << "\n";
    }
    else 
    {
#ifdef _DBG_COUPLEDSEC3D
        opserr << "\n\n";
        opserr << "This is the whole domain (i1 = " << i1 << " - i2 = " << i2 << "): \n";
        this->print();
        opserr << "\n\n";
        opserr << "j = " << j << " - jp1 = " << jp1 << endln;
        opserr << "My1_j = " << My1_j << " - My2_j = " << My2_j << " - My_j (interp) = " << My_j << "\n";
        opserr << "Mz1_j = " << Mz1_j << " - Mz2_j = " << Mz2_j << " - Mz_j (interp) = " << Mz_j << "\n";
        opserr << "My1_jp1 = " << My1_jp1 << " - My2_jp1 = " << My2_jp1 << " - My_jp1 (interp) = " << My_jp1 << "\n";
        opserr << "Mz1_jp1 = " << Mz1_jp1 << " - Mz2_jp1 = " << Mz2_jp1 << " - Mz_jp1 (interp) = " << Mz_jp1 << "\n";
        opserr << "theta_j = " << theta_j << " - theta_jp1 = " << theta_jp1 << "\n";
        opserr << "Looked theta = " << theta << "\n";
#endif
        return -1;
    }
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Values from which interpolate:\n";
    opserr << "My : " << My_j << " and " << My_jp1 << "\n";
    opserr << "Mz : " << Mz_j << " and " << Mz_jp1 << "\n";
#endif
    double alfa = jXdirection / jXjp1;
    double beta = directionXjp1 / jXjp1;
    double sum = alfa + beta;
    alfa /= sum;
    beta /= sum;
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Check alfa, beta\n";
    opserr << alfa << " + " << beta << " = " << alfa+beta << "\n";
#endif
    My = My_j * beta + My_jp1 * alfa;
    Mz = Mz_j * beta + Mz_jp1 * alfa;
   /* if (inverse == 0) 
    {
        My = (My_jp1 - My_j) / (theta_jp1 - theta_j) * (theta - theta_j) + My_j;
        Mz = (Mz_jp1 - Mz_j) / (theta_jp1 - theta_j) * (theta - theta_j) + Mz_j;
    } else 
    {
        My = (My_j - My_jp1) / (theta_j - theta_jp1) * (theta - theta_jp1) + My_jp1;
        Mz = (Mz_j - Mz_jp1) / (theta_j - theta_jp1) * (theta - theta_jp1) + Mz_jp1;
    }*/

#ifdef _DBG_COUPLEDSEC3D
    opserr << "Interpolated values:\n";
    opserr << "My = " << My << " Mz = " << Mz << "\n";
#endif

    return 0;
}


// utilities
namespace {

    int replacePlaceholderWithValue(std::string& theString, const char* placeholder, double value)
    {
        static char buffer[100];
        //static const char* g1 = "%.6g";
        static const char* g1 = "%.6e";
        //static const char* g2 = "%.6g.0";

        std::size_t size = strlen(placeholder);
        std::size_t index = 0;
        while (true) {
            /* Locate the substring to replace. */
            index = theString.find(placeholder, index);
            if (index == std::string::npos) break;

            /* Make the replacement. */
            //const char* gg = (value - int(value) == 0) ? g2 : g1;
            const char* gg = g1;
            sprintf(buffer, gg, value);
            theString.replace(index, size, buffer);

            /* Advance index forward so the next iteration doesn't pick it up as well. */
            index += size;
        }
        return 0;
    }

    // Serialization
#define SerializerFormatDouble std::setprecision(std::numeric_limits<double>::digits10 + 1)
    class Serializer
    {
    private:
        std::stringstream ss;

    public:
        Serializer() = default;
        Serializer(const char* c)
            : ss(c)
        {}
        inline std::string str() const {
            return ss.str();
        }
        explicit operator bool() const {
            return !ss.fail();
        }
        bool operator!() const {
            return ss.fail();
        }

    public:
        inline Serializer& operator << (std::size_t x) {
            ss << x << '\n';
            return *this;
        }
        inline Serializer& operator << (int x) {
            ss << x << '\n';
            return *this;
        }
        inline Serializer& operator << (double x) {
            ss << SerializerFormatDouble << x << '\n';
            return *this;
        }
        inline Serializer& operator << (const std::string& x) {
            ss << x.length() << ' ' << x << '\n';
            return *this;
        }
        inline Serializer& operator << (const std::vector<int>& x) {
            ss << x.size() << '\n';
            for (auto i : x)
                ss << i << '\n';
            return *this;
        }
        inline Serializer& operator << (const Vector& x) {
            ss << x.Size() << '\n';
            for (int i = 0; i < x.Size(); ++i)
                ss << x(i) << '\n';
            return *this;
        }

    public:
        inline Serializer& operator >> (std::size_t& x) {
            ss >> x;
            return *this;
        }
        inline Serializer& operator >> (int& x) {
            ss >> x;
            return *this;
        }
        inline Serializer& operator >> (double& x) {
            ss >> x;
            return *this;
        }
        inline Serializer& operator >> (std::string& x) {
            std::size_t n;
            ss >> n; // needed to make it work even when string is not the first entry
            char dummy;
            ss.read(&dummy, 1); // 1 white space
            x.resize(n);
            ss.read(&x[0], n);
            return *this;
        }
        inline Serializer& operator >> (std::vector<int>& x) {
            std::size_t n;
            if (!(ss >> n))
                return *this;
            x.resize(n);
            for (std::size_t i = 0; i < n; ++i)
                ss >> x[i];
            return *this;
        }
        inline Serializer& operator >> (Vector& x) {
            int n;
            if (!(ss >> n))
                return *this;
            x.resize(n);
            for (int i = 0; i < n; ++i)
                ss >> x(i);
            return *this;
        }
    };

}


void* OPS_ASDCoupledHinge3D() {
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDCoupledHinge3D - Developed by: Diego Talledo, Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARNING insufficient arguments\n"
            "Want: section ASDCoupledHinge3D tag? Ktor? Kvy? Kvz? Kax? -initialFlexuralStiffness \"Ky?\" \"Kz?\""
            "<-simpleStrengthDomain Nmin? Nmax? MyMax? MzMax? NMMax?> <-strengthDomainByPoints nN? nTheta? {listN} {listMy} {listMz}> <-hardening as?>"
            "-thetaP \"thetaPy?\" \"thetaPz?\" -thetaPC \"thetaPCy?\" \"thetaPCz?\""
            "<-hystereticParams rDispP? rForceP? uForceP? rDispN? rForceN? uForceN?>"
            "<-damageUnloadStiffness gk1? gk2? gk3? gk4? gklim?>"
            "<-damageReloadStiffness gd1? gd2? gd3? gd4? gdlim?>"
            "<-damageForce gf1? gf2? gf3? gf4? gflim?>"
            "<-damageType type? ge?>\n";
        return 0;
    }
	    
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid ASDCoupledHinge3D tag" << endln;
        return 0;
    }

    numdata = 4;
    double data[6];
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
        opserr << "WARNING invalid double inputs\n"
            "section ASDCoupledHinge3D: " << tag << endln;
        return 0;
    }
    double Ktor = data[0];
    double Kv_y = data[1];
    double Kv_z = data[2];
    double Kax = data[3];

    //opserr << "Parsed mandatory data\n";

    // start parsing all optional parameters:
    // string tcl expression for initial stiffness: supported parameters for expressions:
    // __N__ : committed axial load
    // __Vy__ : commited shear load in local y direction
    // __Vz__ : commited shear load in local z direction
    // __My__ : commited moment in local y direciton
    // __Mz__ : commited moment in local z direction
    // __M__ : infers by itself the direction y or z
    // __V__ : infers by itself the direction y or z

    // an expression for effective stiffness in y and z local directions
    bool yieldStiffnessDefined = false;
    std::string theRawInitialStiffnessExpressionY;
    std::string theRawInitialStiffnessExpressionZ;
    std::string theInitialStiffnessExpressionY;
    std::string theInitialStiffnessExpressionZ;

    // modality of strength domain providing (one and only one of "simplified" or "by points" need to be provided)
    int strengthDomainMode = STRENGTH_DOMAIN_UNDEFINED;
    // simplified:
    double Nmin = 0; // compression is negative
    double Nmax = 0;
    double MyMax = 0;
    double MzMax = 0;
    double NMMax = 0;
    // by points
    // number of discretizations in N axis
    // number of discretization about angle in My-Mz plane
    int nN = 0;
    int nTheta = 0;
    int n = 0;
    Vector dataN(0);
    Vector dataMy(0);
    Vector dataMz(0);
    ASDCoupledHinge3DDomainData strengthDomain;
    // Domain Values for My and Mz positive and negative when N = 0
    double My_u_p = 0.0;
    double My_u_n = 0.0;
    double Mz_u_p = 0.0;
    double Mz_u_n = 0.0;

    // hardening by default
    double a_s = 1.001;

    // Hysteretic behavior parameters by default
    double rDispP = 0.2, rForceP = 0.5, uForceP = 0.0;
    double rDispN = 0.2, rForceN = 0.5, uForceN = 0.0;
    double gk1 = 0.0, gk2 = 0.0, gk3 = 0.0, gk4 = 0.0, gklim = 0.0;
    double gd1 = 0.0, gd2 = 0.0, gd3 = 0.0, gd4 = 0.0, gdlim = 0.0;
    double gf1 = 0.0, gf2 = 0.0, gf3 = 0.0, gf4 = 0.0, gflim = 0.0;
    double ge = 100.0;
    double dmgType = 1; // cycle by default

    // an expression for pre-cap plastic deformation in y and z local directions
    bool preCapDeformationDefined = false;
    std::string theRawThetaPExpressionY;
    std::string theRawThetaPExpressionZ;
    std::string theThetaPExpressionY;
    std::string theThetaPExpressionZ;

    // an expression for post-cap plastic deformation in y and z local directions
    bool postCapDeformationDefined = false;
    std::string theRawThetaPCExpressionY;
    std::string theRawThetaPCExpressionZ;
    std::string theThetaPCExpressionY;
    std::string theThetaPCExpressionZ;

    //opserr << "Parsing other data\n";

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        if (strcmp(type, "-initialFlexuralStiffness") == 0) {
            //opserr << "Parsing initialFlexuralStiffness\n";
            numdata = 2;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-initialFlexuralStiffness \"Ky?\" \"Kz?\""
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            theRawInitialStiffnessExpressionY = OPS_GetString();
            theRawInitialStiffnessExpressionZ = OPS_GetString();
            theInitialStiffnessExpressionY = theRawInitialStiffnessExpressionY;
            theInitialStiffnessExpressionZ = theRawInitialStiffnessExpressionZ;
            //opserr << theInitialStiffnessExpressionY.c_str() << "\n";
            //opserr << theInitialStiffnessExpressionZ.c_str() << "\n";

            replacePlaceholderWithValue(theInitialStiffnessExpressionY, "__N__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionZ, "__N__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionY, "__M__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionZ, "__M__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionY, "__V__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionZ, "__V__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionY, "__My__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionZ, "__My__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionY, "__Vz__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionZ, "__Vz__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionY, "__Mz__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionZ, "__Mz__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionY, "__Vy__", 0.0);
            replacePlaceholderWithValue(theInitialStiffnessExpressionZ, "__Vy__", 0.0);
#ifdef _DBG_COUPLEDSEC3D
            opserr << theInitialStiffnessExpressionY.c_str() << "\n";
            opserr << theInitialStiffnessExpressionZ.c_str() << "\n";
#endif

            yieldStiffnessDefined = true;
        } else if (strcmp(type, "-thetaP") == 0) {
            //opserr << "Parsing thetaP\n";
            numdata = 2;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-thetaP \"thetaPy?\" \"thetaPz?\""
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            theRawThetaPExpressionY = OPS_GetString();
            theRawThetaPExpressionZ = OPS_GetString();
            theThetaPExpressionY = theRawThetaPExpressionY;
            theThetaPExpressionZ = theRawThetaPExpressionZ;
            //opserr << theThetaPExpressionY.c_str() << "\n";
            //opserr << theThetaPExpressionZ.c_str() << "\n";

            replacePlaceholderWithValue(theThetaPExpressionY, "__N__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionZ, "__N__", 1.1);
            replacePlaceholderWithValue(theThetaPExpressionY, "__M__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionZ, "__M__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionY, "__V__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionZ, "__V__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionY, "__My__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionZ, "__My__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionY, "__Vz__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionZ, "__Vz__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionY, "__Mz__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionZ, "__Mz__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionY, "__Vy__", 0.0);
            replacePlaceholderWithValue(theThetaPExpressionZ, "__Vy__", 0.0);
#ifdef _DBG_COUPLEDSEC3D
            opserr << theThetaPExpressionY.c_str() << "\n";
            opserr << theThetaPExpressionZ.c_str() << "\n";
#endif
            preCapDeformationDefined = true;
        } else if (strcmp(type, "-thetaPC") == 0) {
            //opserr << "Parsing thetaPC\n";
            numdata = 2;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-thetaPC \"thetaPCy?\" \"thetaPCz?\""
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            theRawThetaPCExpressionY = OPS_GetString();
            theRawThetaPCExpressionZ = OPS_GetString();
            theThetaPCExpressionY = theRawThetaPCExpressionY;
            theThetaPCExpressionZ = theRawThetaPCExpressionZ;
            //opserr << theThetaPCExpressionY.c_str() << "\n";
            //opserr << theThetaPCExpressionZ.c_str() << "\n";

            replacePlaceholderWithValue(theThetaPCExpressionY, "__N__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionZ, "__N__", 1.1);
            replacePlaceholderWithValue(theThetaPCExpressionY, "__M__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionZ, "__M__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionY, "__V__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionZ, "__V__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionY, "__My__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionZ, "__My__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionY, "__Vz__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionZ, "__Vz__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionY, "__Mz__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionZ, "__Mz__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionY, "__Vy__", 0.0);
            replacePlaceholderWithValue(theThetaPCExpressionZ, "__Vy__", 0.0);

            //opserr << theThetaPCExpressionY.c_str() << "\n";
            //opserr << theThetaPCExpressionZ.c_str() << "\n";

            postCapDeformationDefined = true;
        } else  if (strcmp(type, "-hardening") == 0) {
            numdata = 1;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-hardening a_s?"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, data) < 0) {
                opserr << "WARNING invalid double inputs -hardening\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            a_s = data[0];
        } else  if (strcmp(type, "-hystereticParams") == 0) {
            numdata = 6;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-hystereticParams rDispP? rForceP? uForceP? rDispN? rForceN? uForceN?"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, data) < 0) {
                opserr << "WARNING invalid double inputs -hystereticParams\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            rDispP = data[0]; rForceP = data[1]; uForceP = data[2];
            rDispN = data[3]; rForceN = data[4]; uForceN = data[5];
        } else  if (strcmp(type, "-damageUnloadStiffness") == 0) {
            numdata = 5;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-damageUnloadStiffness gk1? gk2? gk3? gk4? gklim?"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, data) < 0) {
                opserr << "WARNING invalid double inputs -damageUnloadStiffness\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            gk1 = data[0]; gk2 = data[1]; gk3 = data[2]; gk4 = data[3]; gklim = data[4];
        } else  if (strcmp(type, "-damageReloadStiffness") == 0) {
            numdata = 5;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-damageUnloadStiffness gk1? gk2? gk3? gk4? gklim?"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, data) < 0) {
                opserr << "WARNING invalid double inputs -damageReloadStiffness\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            gd1 = data[0]; gd2 = data[1]; gd3 = data[2]; gd4 = data[3]; gdlim = data[4];
        } else  if (strcmp(type, "-damageForce") == 0) {
            numdata = 5;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-damageForce gf1? gf2? gf3? gf4? gflim?"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, data) < 0) {
                opserr << "WARNING invalid double inputs -damageForce\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            gf1 = data[0]; gf2 = data[1]; gf3 = data[2]; gf4 = data[3]; gflim = data[4];
        } else  if (strcmp(type, "-damageType") == 0) {
            const char* type = OPS_GetString();
            if (strcmp(type, "cycle") == 0 || strcmp(type, "Cycle") == 0 || strcmp(type, "DamageCycle") == 0 || strcmp(type, "damageCycle") == 0) {
                dmgType = 1;
            }
            else if (strcmp(type, "energy") == 0 || strcmp(type, "Energy") == 0 || strcmp(type, "DamageEnergy") == 0 || strcmp(type, "damageEnergy") == 0) {
                dmgType = 0;
            }
            else {
                opserr << "WARNING invalid type of damage calculation specified\n";
                opserr << "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            numdata = 1;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-damageType type? ge?"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, data) < 0) {
                opserr << "WARNING invalid double inputs -damageType\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            ge = data[0];
        } else if (strcmp(type, "-simpleStrengthDomain") == 0) {
            if (strengthDomainMode != STRENGTH_DOMAIN_UNDEFINED) {
                opserr << "only one strength domain need to be defined!\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            // simplified strength domain, giving few data.
            // we assume a double symetric - rect - section so that MinMoment = - MaxMoment
            // we assume also that NMMax is corresponding to both My and Mz
            numdata = 5;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-simpleStrengthDomain Nmin? Nmax? MyMax? MzMax? NMMax?"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            if (OPS_GetDoubleInput(&numdata, data) < 0) {
                opserr << "WARNING invalid double inputs -simpleStrengthDomain\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            Nmin = data[0]; 
            Nmax = data[1];
            MyMax = data[2];
            MzMax = data[3];
            NMMax = data[4];

            strengthDomainMode = STRENGTH_DOMAIN_SIMPLIFIED;

            // create a simplified domain. In future with some interpolation
            nN = 3;
            nTheta = 4;

            // Create the domain as a DomainData object
            strengthDomain = ASDCoupledHinge3DDomainData(nN, nTheta, 3);
            for (int j = 0; j < nTheta; j++) {
                strengthDomain.setValue(0, j, 0, Nmin);
                strengthDomain.setValue(2, j, 0, Nmax);
                for (int k = 1; k < 3; k++) {
                    strengthDomain.setValue(0, j, k, 0.0);
                    strengthDomain.setValue(2, j, k, 0.0);
                }
            }
            // write points of maximum moment
            for (int j = 0; j < nTheta; j++) {
                strengthDomain.setValue(1, j, 0, NMMax);
            }
            strengthDomain.setValue(1, 0, 1, MyMax);
            strengthDomain.setValue(1, 0, 2, 0.0);
            strengthDomain.setValue(1, 1, 1, 0.0);
            strengthDomain.setValue(1, 1, 2, MzMax);
            strengthDomain.setValue(1, 2, 1, -MyMax);
            strengthDomain.setValue(1, 2, 2, 0.0);
            strengthDomain.setValue(1, 3, 1, 0.0);
            strengthDomain.setValue(1, 3, 2, -MzMax);

#ifdef _DBG_COUPLEDSEC3D
            opserr << "strengthDomain created: " << endln;
            strengthDomain.print();
#endif
        } else if (strcmp(type, "-strengthDomainByPoints") == 0) {
            //opserr << "strengthDomainByPoints -> reading data\n";
            // strength domain given by points
            // it need to be correctly structured
            if (strengthDomainMode != STRENGTH_DOMAIN_UNDEFINED) {
                opserr << "only one strength domain need to be defined!\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            if (OPS_GetNumRemainingInputArgs() < 5) {
                opserr << "WARNING insufficient arguments\n"
                    "<-strengthDomainByPoints nN? nTheta? {listN}? {listMy}? {listMz}?>\n";
                return 0;
            }

            numdata = 1;
            // Read the data the provides the number of discretizations
            if (OPS_GetIntInput(&numdata, &nN) < 0) {
                opserr << "WARNING invalid ASDCoupledHinge3D nN" << endln;
                return 0;
            }
            if (OPS_GetIntInput(&numdata, &nTheta) < 0) {
                opserr << "WARNING invalid ASDCoupledHinge3D nTheta" << endln;
                return 0;
            }

            // total number of data I need in each list
            n = nN * nTheta;

            // read N list
            int sizeN;
            if (OPS_GetDoubleListInput(&sizeN, &dataN)) {
                opserr << "WARNING An error occurred reading the list of dataN" << endln;
                return 0;
            }
            
            // read My list
            int sizeMy;
            if (OPS_GetDoubleListInput(&sizeMy, &dataMy)) {
                opserr << "WARNING An error occurred reading the list of dataMy" << endln;
                return 0;
            }

            // read Mz list
            int sizeMz;
            if (OPS_GetDoubleListInput(&sizeMz, &dataMz)) {
                opserr << "WARNING An error occurred reading the list of dataMz" << endln;
                return 0;
            }

            if (sizeN != sizeMy || sizeN != sizeMz) {
                opserr << "WARNING ASDCoupledHinge3D - vector containing data ";
                opserr << "points for N, My and Mz are not of the same size\n";
                return 0;
            } 

            if (sizeN != n) {
                opserr << "WARNING ASDCoupledHinge3D - vector containing data ";
                opserr << "Number of elements passed should be equal to " << n << "(" << nN << " x " << nTheta << ")f\n";
                return 0;
            }

            // Create the domain as a DomainData object
            strengthDomain = ASDCoupledHinge3DDomainData(nN, nTheta, 3);

            // write all elements of the domain data
            int idx;
            for (int i = 0; i < nN; i++) {
                for (int j = 0; j < nTheta; j++) {
                    idx = i * nTheta + j;
                    strengthDomain.setValue(i, j, 0, dataN(idx));
                    strengthDomain.setValue(i, j, 1, dataMy(idx));
                    strengthDomain.setValue(i, j, 2, dataMz(idx));
                }
            }
#ifdef _DBG_COUPLEDSEC3D
            opserr << "strengthDomain created: " << endln;
            strengthDomain.print();
#endif
            strengthDomainMode = STRENGTH_DOMAIN_BY_POINTS;
        }
    }

    opserr << "Parsed the command. Creating the materials.\n";

    // Check that a domain is created
    if (strengthDomainMode == STRENGTH_DOMAIN_UNDEFINED) {
        opserr << "A domain needs to be created providing simpleStrengthDomain or strengthDomainByPoints\n";
        return 0;
    }

    // get Domain Values for My and Mz positive and negative when N = 0
    double tmp;
    strengthDomain.getMyMzForNAndDirection(0.0, 0.0, My_u_p, tmp);
    strengthDomain.getMyMzForNAndDirection(0.0, M_PI, My_u_n, tmp);
    strengthDomain.getMyMzForNAndDirection(0.0, M_PI / 2, tmp, Mz_u_p);
    strengthDomain.getMyMzForNAndDirection(0.0, M_PI * 3 / 2.0, tmp, Mz_u_n);

#ifdef _DBG_COUPLEDSEC3D
    opserr << "Domain intersections after creation" << endln;
    opserr << "For N = 0 \nMy+ = " << My_u_p << endln;
    opserr << "My- = " << My_u_n << endln;
    opserr << "Mz+ = " << Mz_u_p << endln;
    opserr << "Mz- = " << Mz_u_n << endln;

    opserr << "Domain is this one: \n";
    strengthDomain.print();
    opserr << endln;
#endif

    // Create an elastic material for axial
    UniaxialMaterial* matAxial = new ElasticMaterial(0, Kax, 0.0, Kax);
    if (matAxial == 0) {
        opserr << "ASDCoupledHinge3D could not create the new axial material" << endln;
        return 0;
    }
    // Create an elastic material for shear_y
    UniaxialMaterial* matShearY = new ElasticMaterial(0, Kv_y, 0.0, Kv_y);
    if (matShearY == 0) {
        opserr << "ASDCoupledHinge3D could not create the new shearY material" << endln;
        return 0;
    }
    // Create an elastic material for shear_z
    UniaxialMaterial* matShearZ = new ElasticMaterial(0, Kv_z, 0.0, Kv_z);
    if (matShearZ == 0) {
        opserr << "ASDCoupledHinge3D could not create the new shearZ material" << endln;
        return 0;
    }
    // Create an elastic material for shear_z
    UniaxialMaterial* matTorsion = new ElasticMaterial(0, Ktor, 0.0, Ktor);
    if (matTorsion == 0) {
        opserr << "ASDCoupledHinge3D could not create the new torsion material" << endln;
        return 0;
    }
    opserr << "Materials for Torsion, Axial and Shear (Y and Z) defined\n";

    // Create pinching4 material for My
    double EI;
    if (OPS_EvalDoubleStringExpression(theInitialStiffnessExpressionY.c_str(), EI) < 0) {
        opserr << "Error evaluating expression:\n" << theInitialStiffnessExpressionY.c_str() << "\n";
    }
    double th_cap;
    if (OPS_EvalDoubleStringExpression(theThetaPExpressionY.c_str(), th_cap) < 0) {
        opserr << "Error evaluating expression:\n" << theThetaPExpressionY.c_str() << "\n";
    }
    double th_pcap;
    if (OPS_EvalDoubleStringExpression(theThetaPCExpressionY.c_str(), th_pcap) < 0) {
        opserr << "Error evaluating expression:\n" << theThetaPCExpressionY.c_str() << "\n";
    }
    
#ifdef _DBG_COUPLEDSEC3D
    opserr << theInitialStiffnessExpressionY.c_str() << " = " << EI << endln;
    opserr << theThetaPExpressionY.c_str() << " = " << th_cap << endln;
    opserr << theThetaPCExpressionY.c_str() << " = " << th_pcap << endln;
#endif
    // Yield strength from Domain
    double f2p = My_u_p;
    double d2p = My_u_p / EI;
    // First point in the same line (arbitraly choosen 1/5)
    double f1p = f2p / 5.0;
    double d1p = d2p / 5.0;
    // Hardening point - cap
    double f3p = a_s * f2p;
    double d3p = th_cap; // alternativelly: d2p + th_cap; (but user need to know that th_cap_p need to be input)
    // Ultimate point
    double f4p = 0.1 * My_u_p; // To be choosen by user later ???
    double d4p = th_pcap; // alternativelly: d3p + th_pcap; (but user need to know that th_pcap_p need to be input)
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Initial spring My + : (" << f1p << ", " << d1p << "), (" << f2p << ", " << d2p << "), (" << f3p << ", " << d3p << "), (" << f4p << ", " << d4p << ")\n";
#endif
    // Yield strength
    double f2n = My_u_n;
    double d2n = My_u_n / EI;
    // First point in the same line (arbitraly choosen 1/5)
    double f1n = f2n / 5.0;
    double d1n = d2n / 5.0;
    // Hardening point - cap
    double f3n = a_s * f2n;
    double d3n = -th_cap; // alternativelly: d2n - th_cap; (but user need to know that th_cap_p need to be input)
    // Ultimate point
    double f4n = 0.1 * My_u_n; // To be choosen by user in future versions; For now 10% of My
    double d4n = -th_pcap; // alternatively: d3n - th_pcap; (but user need to know that th_pcap_p need to be input)
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Initial spring My - : (" << f1n << ", " << d1n << "), (" << f2n << ", " << d2n << "), (" << f3n << ", " << d3n << "), (" << f4n << ", " << d4n << ")\n";
#endif

    // Create a Pinching4 material for My
    double mdp = rDispP;
    double mfp = rForceP;
    double msp = uForceP;
    double mdn = rDispN;
    double mfn = rForceN;
    double msn = uForceN;
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Hysteretic parameters used for pinching4\n";
    opserr << "mdp = " << mdp << " - mfp = " << mfp << " - msp = " << msp << endln;
    opserr << "mdn = " << mdn << " - mfn = " << mfn << " - msn = " << msn << endln;
    opserr << "gk1 = " << gk1 << " - gk2 = " << gk2 << " - gk3 = " << gk3 << " - gk4 = " << gk4 << " - gklim = " << gklim << endln;
    opserr << "gd1 = " << gd1 << " - gd2 = " << gd2 << " - gd3 = " << gd3 << " - gd4 = " << gd4 << " - gdlim = " << gdlim << endln;
    opserr << "gf1 = " << gf1 << " - gf2 = " << gf2 << " - gf3 = " << gf3 << " - gf4 = " << gf4 << " - gflim = " << gflim << endln;
    opserr << "ge = " << ge << " - dmgType = " << dmgType << endln;
#endif
    UniaxialMaterial* matMy = new Pinching4Material(0,f1p,d1p,f2p,d2p,f3p,d3p,f4p,d4p,f1n,d1n,f2n,d2n,f3n,d3n,f4n,d4n,mdp,mfp,msp,mdn,mfn,msn,gk1,gk2,gk3,gk4,gklim,gd1,gd2,gd3,gd4,gdlim,gf1,gf2,gf3,gf4,gflim,ge,dmgType);
#ifdef _DBG_COUPLEDSEC3D
    matMy->Print(opserr, 2);
    opserr << endln;
#endif
    
    // Create pinching4 material for Mz
    if (OPS_EvalDoubleStringExpression(theInitialStiffnessExpressionZ.c_str(), EI) < 0) {
        opserr << "Error evaluating expression:\n" << theInitialStiffnessExpressionZ.c_str() << "\n";
    }
    if (OPS_EvalDoubleStringExpression(theThetaPExpressionZ.c_str(), th_cap) < 0) {
        opserr << "Error evaluating expression:\n" << theThetaPExpressionZ.c_str() << "\n";
    }
    if (OPS_EvalDoubleStringExpression(theThetaPCExpressionZ.c_str(), th_pcap) < 0) {
        opserr << "Error evaluating expression:\n" << theThetaPCExpressionZ.c_str() << "\n";
    }
#ifdef _DBG_COUPLEDSEC3D
    opserr << theInitialStiffnessExpressionZ.c_str() << " = " << EI << endln;
    opserr << theThetaPExpressionZ.c_str() << " = " << th_cap << endln;
    opserr << theThetaPCExpressionZ.c_str() << " = " << th_pcap << endln;
#endif
    // Yield strength
    f2p = Mz_u_p;
    d2p = Mz_u_p / EI;
    // First point in the same line (arbitraly choosen 1/5)
    f1p = f2p / 5.0;
    d1p = d2p / 5.0;
    // Hardening point - cap
    f3p = a_s * f2p;
    d3p = th_cap; // alternativelly: d2p + th_cap; (but user need to know that th_cap_p need to be input)
    // Ultimate point
    f4p = 0.1 * Mz_u_p; // To be choosen by user later ???
    d4p = th_pcap; // alternativelly: d3p + th_pcap; (but user need to know that th_pcap_p need to be input)
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Initial spring Mz + : (" << f1p << ", " << d1p << "), (" << f2p << ", " << d2p << "), (" << f3p << ", " << d3p << "), (" << f4p << ", " << d4p << ")\n";
#endif
    // Yield strength
    f2n = Mz_u_n;
    d2n = Mz_u_n / EI;
    // First point in the same line (arbitraly choosen 1/5)
    f1n = f2n / 5.0;
    d1n = d2n / 5.0;
    // Hardening point - cap
    f3n = a_s * f2n;
    d3n = -th_cap; // alternativelly: d2n - th_cap; (but user need to know that th_cap_p need to be input)
    // Ultimate point
    f4n = 0.1 * Mz_u_n; // To be choosen by user later ???
    d4n = -th_pcap; // alternativelly d3n - th_pcap; (but user need to know that th_pcap_p need to be input)
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Initial spring Mz - : (" << f1n << ", " << d1n << "), (" << f2n << ", " << d2n << "), (" << f3n << ", " << d3n << "), (" << f4n << ", " << d4n << ")\n";
#endif

    // Create a Pinching4 material for Mz - For now same as Y direction for histeretic properties
    UniaxialMaterial* matMz = new Pinching4Material(0, f1p, d1p, f2p, d2p, f3p, d3p, f4p, d4p, f1n, d1n, f2n, d2n, f3n, d3n, f4n, d4n, mdp, mfp, msp, mdn, mfn, msn, gk1, gk2, gk3, gk4, gklim, gd1, gd2, gd3, gd4, gdlim, gf1, gf2, gf3, gf4, gflim, ge, dmgType);
#ifdef _DBG_COUPLEDSEC3D
    matMz->Print(opserr, 2);
    opserr << endln;
#endif
    // Create the new section object

    // Dati da passare in qualche modo
    /*std::string theRawInitialStiffnessExpressionY;
    std::string theRawInitialStiffnessExpressionZ;
    std::string theRawThetaPExpressionY;
    std::string theRawThetaPExpressionZ;
    std::string theRawThetaPCExpressionY;
    std::string theRawThetaPCExpressionZ;
    double a_s;*/
    ASDCoupledHinge3D *theSection = new ASDCoupledHinge3D(tag, matTorsion, matAxial, matShearY, matShearZ, matMy, matMz, strengthDomain, theRawInitialStiffnessExpressionY, theRawInitialStiffnessExpressionZ, theRawThetaPExpressionY, theRawThetaPExpressionZ, theRawThetaPCExpressionY, theRawThetaPCExpressionZ, a_s);

    // Now I delete them because I've already copied them inside the secion object.
    delete matAxial;
    delete matMy;
    delete matMz;
    delete matShearY;
    delete matShearZ;
    delete matTorsion;

#ifdef _DBG_COUPLEDSEC3D
    opserr << "ASDCoupledHinge3D Section created\n\n";
#endif

    return theSection;
}

// constructors:
ASDCoupledHinge3D::ASDCoupledHinge3D(
    int tag, 
    UniaxialMaterial* theTorsionMaterial, 
    UniaxialMaterial* theAxialMaterial, 
    UniaxialMaterial* theShearYMaterial, 
    UniaxialMaterial* theShearZMaterial, 
    UniaxialMaterial* theMomentYMaterial,
    UniaxialMaterial* theMomentZMaterial,
    const ASDCoupledHinge3DDomainData &theStrengthDomain, 
    const std::string &theRawInitialStiffnessExpressionY,
    const std::string &theRawInitialStiffnessExpressionZ, 
    const std::string &theRawThetaPExpressionY, 
    const std::string &theRawThetaPExpressionZ, 
    const std::string &theRawThetaPCExpressionY, 
    const std::string &theRawThetaPCExpressionZ, 
    double a_s_i)
    : SectionForceDeformation(tag, SEC_TAG_ASDCoupledHinge3D)
    , strengthDomain(theStrengthDomain)
    , rawInitialStiffnessExpressionY(theRawInitialStiffnessExpressionY)
    , rawInitialStiffnessExpressionZ(theRawInitialStiffnessExpressionZ)
    , rawThetaPExpressionY(theRawThetaPExpressionY)
    , rawThetaPExpressionZ(theRawThetaPExpressionZ)
    , rawThetaPCExpressionY(theRawThetaPCExpressionY)
    , rawThetaPCExpressionZ(theRawThetaPCExpressionZ)
    , a_s(a_s_i)
{
    // Create a copy of the materials passed to the constructor
    if ((!theAxialMaterial) || (!theMomentYMaterial) || (!theMomentZMaterial) || (!theTorsionMaterial) || (!theShearYMaterial) || (!theShearZMaterial)) {
        opserr << "ASDCoupledHinge3D::ASDCoupledHinge3D " << tag << " -- null uniaxial material passed\n";
        exit(-1);
    }
    axialMaterial = theAxialMaterial->getCopy();
    MyMaterial = theMomentYMaterial->getCopy();
    MzMaterial = theMomentZMaterial->getCopy();
    torsionMaterial = theTorsionMaterial->getCopy();
    VyMaterial = theShearYMaterial->getCopy();
    VzMaterial = theShearZMaterial->getCopy();

    // setup parameters
    setupParameters();

    // Compute the needed tolerances that are relative to the domain:
    // tolN = (Nmax - Nmin)*0.001; 1 per thousand of deltaN
    double Nmin, Nmax;
    strengthDomain.getRangeN(Nmin, Nmax);
    tolN = (Nmax - Nmin) * 0.001;
    double MminY, MmaxY, MminZ, MmaxZ;
    setUncoupledStrengthDomainforAxial(0.0, MmaxY, MminY, MmaxZ, MminZ);
    MmaxAbs = std::max(std::max(MmaxY, std::abs(MminY)), std::max(MmaxZ, std::abs(MminZ)));
    tolM = 0.0001 * MmaxAbs;

#ifdef _DBG_COUPLEDSEC3D
    opserr << "The copied ultimate domain: \n";
    strengthDomain.print();
#endif
}

// constructor for blank object that recvSelf needs to be invoked upon
ASDCoupledHinge3D::ASDCoupledHinge3D():
    SectionForceDeformation(0, SEC_TAG_ASDCoupledHinge3D),
    otherDbTag(0), axialMaterial(0), MyMaterial(0), MzMaterial(0), VyMaterial(0), VzMaterial(0), torsionMaterial(0)
{

}
//
// destructor:
ASDCoupledHinge3D::~ASDCoupledHinge3D()
{
    if (axialMaterial) {
        delete axialMaterial;
    }
    if (MyMaterial) {
        delete MyMaterial;
    }
    if (MzMaterial) {
        delete MzMaterial;
    }
    if (VyMaterial) {
        delete VyMaterial;
    }
    if (VzMaterial) {
        delete VzMaterial;
    }
    if (torsionMaterial) {
        delete torsionMaterial;
    }
}

int ASDCoupledHinge3D::setTrialSectionDeformation (const Vector &def)
{
  int ret = 0;

#ifdef  ASD_HINGE_NUM_TANG

  static constexpr double pert = 1.0e-9;
  auto lam_pert = [this](UniaxialMaterial* mat, double strain, int order) {
      mat->setTrialStrain(strain + pert);
      double sh = mat->getStress();
      int ret = mat->setTrialStrain(strain);
      double ki = (sh - mat->getStress()) / pert;
      m_num_tang(order) = ki;
      return ret;
};
  ret += lam_pert(axialMaterial, def(0), 0);
  ret += lam_pert(MyMaterial, def(1), 1);
  ret += lam_pert(MzMaterial, def(2), 2);
  ret += lam_pert(VyMaterial, def(3), 3);
  ret += lam_pert(VzMaterial, def(4), 4);
  ret += lam_pert(torsionMaterial, def(5), 5);
#ifdef _DBG_COUPLEDSEC3D
  opserr << "ASDCoupledHinge3D::setTrialSectionDeformation\n";
  opserr << "My tangent: " << m_num_tang(1);
  opserr << " Mz tangent: " << m_num_tang(2) << endln;
  opserr << "set trial strain to My law: " << def(1) << endln;
  opserr << "-- My law: " << endln;
  MyMaterial->Print(opserr, 3);
  opserr << "\n";
  opserr << "set trial strain to Mz law: " << def(2) << endln;
  opserr << "-- Mz law: " << endln;
  MzMaterial->Print(opserr, 3);
  opserr << "\n";
#endif

#else

  ret += axialMaterial->setTrialStrain(def(0));
  ret += MyMaterial->setTrialStrain(def(1));
  ret += MzMaterial->setTrialStrain(def(2));
  ret += VyMaterial->setTrialStrain(def(3));
  ret += VzMaterial->setTrialStrain(def(4));
  ret += torsionMaterial->setTrialStrain(def(5));

#endif //  ASD_HINGE_NUM_TANG

  return ret;
}

const Vector &
ASDCoupledHinge3D::getSectionDeformation(void)
{
    static Vector e(6);
    e(0) = axialMaterial->getStrain();
    e(1) = MyMaterial->getStrain();
    e(2) = MzMaterial->getStrain();
    e(3) = VyMaterial->getStrain();
    e(4) = VzMaterial->getStrain();
    e(5) = torsionMaterial->getStrain();
    return e;
}

const Matrix &
ASDCoupledHinge3D::getSectionTangent(void)
{
    static Matrix k(6, 6);
#ifdef ASD_HINGE_NUM_TANG
    for (int i = 0; i < 6; ++i)
        k(i, i) = m_num_tang(i);
#else
    k(0, 0) = axialMaterial->getTangent();
    k(1, 1) = MyMaterial->getTangent();
    k(2, 2) = MzMaterial->getTangent();
    k(3, 3) = VyMaterial->getTangent();
    k(4, 4) = VzMaterial->getTangent();
    k(5, 5) = torsionMaterial->getTangent();
#endif // ASD_HINGE_NUM_TANG

    return k;
}

const Matrix &
ASDCoupledHinge3D::getInitialTangent(void)
{
    static Matrix k(6, 6);
    k(0, 0) = axialMaterial->getInitialTangent();
    k(1, 1) = MyMaterial->getInitialTangent();
    k(2, 2) = MzMaterial->getInitialTangent();
    k(3, 3) = VyMaterial->getInitialTangent();
    k(4, 4) = VzMaterial->getInitialTangent();
    k(5, 5) = torsionMaterial->getInitialTangent();
    return k;
}

const Matrix &
ASDCoupledHinge3D::getSectionFlexibility(void)
{
    static Matrix f(6, 6);
    
#ifdef ASD_HINGE_NUM_TANG

    for (int i = 0; i < 6; ++i) {
        double ki = m_num_tang(i);
        if (ki == 0.0)
            ki = 1.0e-14;
        double fi = 1.0 / ki;
        f(i, i) = fi;
    }

#else

    double k;
    k = axialMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(0, 0) = 1.e14;
    }
    else {
        f(0, 0) = 1 / k;
    }
    k = MyMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(1, 1) = 1.e14;
    }
    else {
        f(1, 1) = 1 / k;
    }
    k = MzMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(2, 2) = 1.e14;
    }
    else {
        f(2, 2) = 1 / k;
    }
    k = VyMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(3, 3) = 1.e14;
    }
    else {
        f(3, 3) = 1 / k;
    }
    k = VzMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(4, 4) = 1.e14;
    }
    else {
        f(4, 4) = 1 / k;
    }
    k = torsionMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(5, 5) = 1.e14;
    }
    else {
        f(5, 5) = 1 / k;
    }

#endif
  
    return f;
}

const Matrix &
ASDCoupledHinge3D::getInitialFlexibility(void)
{
    static Matrix f(6, 6);

    double k;
    k = axialMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(0, 0) = 1.e14;
    }
    else {
        f(0, 0) = 1 / k;
    }
    k = MyMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(1, 1) = 1.e14;
    }
    else {
        f(1, 1) = 1 / k;
    }
    k = MzMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(2, 2) = 1.e14;
    }
    else {
        f(2, 2) = 1 / k;
    }
    k = VyMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(3, 3) = 1.e14;
    }
    else {
        f(3, 3) = 1 / k;
    }
    k = VzMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(4, 4) = 1.e14;
    }
    else {
        f(4, 4) = 1 / k;
    }
    k = torsionMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        f(5, 5) = 1.e14;
    }
    else {
        f(5, 5) = 1 / k;
    }

    return f;
}

const Vector &
ASDCoupledHinge3D::getStressResultant(void)
{
    static Vector s(6);
    s(0) = axialMaterial->getStress();
    s(1) = MyMaterial->getStress();
    s(2) = MzMaterial->getStress();
    s(3) = VyMaterial->getStress();
    s(4) = VzMaterial->getStress();
    s(5) = torsionMaterial->getStress();
    return s;
}

SectionForceDeformation *
ASDCoupledHinge3D::getCopy(void)
{
  ASDCoupledHinge3D *theCopy = 0;

  theCopy = new ASDCoupledHinge3D(this->getTag(), torsionMaterial, axialMaterial, VyMaterial, VzMaterial, MyMaterial, MzMaterial, strengthDomain, rawInitialStiffnessExpressionY, rawInitialStiffnessExpressionZ, rawThetaPExpressionY, rawThetaPExpressionZ, rawThetaPCExpressionY, rawThetaPCExpressionZ, a_s);
  
  if (theCopy == 0) {
    opserr << "ASDCoupledHinge3D::getCopy -- failed to allocate copy\n";
    exit(-1);
  }
		
  return theCopy;
}

const ID&
ASDCoupledHinge3D::getType ()
{
    auto lam = []() {
        ID x(6);
        x(0) = SECTION_RESPONSE_P;
        x(1) = SECTION_RESPONSE_MY;
        x(2) = SECTION_RESPONSE_MZ;
        x(3) = SECTION_RESPONSE_VY;
        x(4) = SECTION_RESPONSE_VZ;
        x(5) = SECTION_RESPONSE_T;
        return x;
    };
    static ID theCode = lam();
    return theCode;
}

int
ASDCoupledHinge3D::getOrder () const
{
    return 6;
}

void
ASDCoupledHinge3D::updateLaws(void)
{
    int err = 0;

    // Get forces from materials
    double N = axialMaterial->getStress();
    double My = MyMaterial->getStress();
    double Mz = MzMaterial->getStress();
    double Vy = VyMaterial->getStress();
    double Vz = VzMaterial->getStress();
#ifdef _DBG_COUPLEDSEC3D
    opserr << "ASDCoupledHinge3D::updateLaws\n";
    opserr << "My getted from getStress: " << My << "\n";
    opserr << "Mz getted from getStress: " << Mz << "\n";
#endif

    // Idea: get ky and kz and multiply by elastic stiffness in order to see direction? I am not sure because elastic stiffness can change during analysis...?
    double ky = MyMaterial->getStrain();
    double kz = MzMaterial->getStrain();
    if (ky >= 0)
        My = par_f1p_Y.getValue() / par_d1p_Y.getValue() * ky;
    else
        My = par_f1n_Y.getValue() / par_d1n_Y.getValue() * ky;
    if (kz >= 0)
        Mz = par_f1p_Z.getValue() / par_d1p_Z.getValue() * kz;
    else
        Mz = par_f1n_Z.getValue() / par_d1n_Z.getValue() * kz;

    // use a relvative measure for tolerance: lets use (Nmax - Nmin)*0.01 (1% of deltaN)
    double My_u_p, My_u_n, Mz_u_p, Mz_u_n;

    //I need to do this each time because scaling does not work, so we need to compute the two laws at each load increment

    // get the uncoupled laws in Y and Z
    //setUncoupledStrengthDomainforAxial(N, My_u_p, My_u_n, Mz_u_p, Mz_u_n);
//#ifdef _DBG_COUPLEDSEC3D
//    opserr << "For N = " << N << " the domain is intercepted at: My = " << My_u_p << " My- = " << My_u_n << " Mz+ = " << Mz_u_p << " Mz- = " << Mz_u_n << "\n";
//    opserr << "New values for strength:\n";
//    opserr << "For N = " << N << "\nMy+ = " << My_u_p << endln;
//    opserr << "My- = " << My_u_n << endln;
//    opserr << "Mz+ = " << Mz_u_p << endln;
//    opserr << "Mz- = " << Mz_u_n << endln;
//#endif

    // Check moments are never exactly zero to avoid numerical problems
    if (std::abs(My) < tolM) {
        if (My >= 0)
            My = tolM;
        else
            My = -tolM;
    }
    if (std::abs(Mz) < tolM) {
        if (Mz >= 0)
            Mz = tolM;
        else
            Mz = -tolM;
    }
    // compute the direction of moments
    double theta = atan2(Mz, My);
    if (theta < 0)
        theta += 2 * M_PI;
    // Now compute the points on the domain with the same direction
#ifdef _DBG_COUPLEDSEC3D
    opserr << "point 1 Y " << par_f1p_Y.getValue() << " " << par_d1p_Y.getValue() << endln;
    opserr << "point 1 Z " << par_f1p_Z.getValue() << " " << par_d1p_Z.getValue() << endln;
    opserr << "ky = " << ky << " kz = " << kz << endln;
    opserr << "My = " << My << " (tolM = " << tolM << ")" << endln;
    opserr << "Mz = " << Mz << " (tolM = " << tolM << ")" << endln;
#endif
    strengthDomain.getMyMzForNAndDirection(N, theta, My_u_p, Mz_u_p);
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Both moments above tolerance -> biaxial bending\n";
    opserr << "My = " << My << " - Mz = " << Mz << "\n";
    opserr << "angle theta = " << theta << "\n";
    opserr << "getMyMz for N = " << N << " : My = " << My_u_p << " - Mz = " << Mz_u_p << "\n";
#endif
    theta += M_PI;
    if (theta > 2 * M_PI)
        theta -= 2 * M_PI;
    strengthDomain.getMyMzForNAndDirection(N, theta, My_u_n, Mz_u_n);
#ifdef _DBG_COUPLEDSEC3D
    opserr << "getMyMz for N = " << N << " : My = " << My_u_n << " - Mz = " << Mz_u_n << "\n";
    opserr << "opposite angle theta = " << theta << "\n";
#endif
    double tmp;
    if (My_u_p < My_u_n) {
        tmp = My_u_p;
        My_u_p = My_u_n;
        My_u_n = tmp;
    }
    if (Mz_u_p < Mz_u_n) {
        tmp = Mz_u_p;
        Mz_u_p = Mz_u_n;
        Mz_u_n = tmp;
    }
#ifdef _DBG_COUPLEDSEC3D
    opserr << "The intersection points with the domain define the following strengths: \n";
    opserr << "New values for strength:\n";
    opserr << "For N = " << N << "\nMy+ = " << My_u_p << endln;
    opserr << "My- = " << My_u_n << endln;
    opserr << "Mz+ = " << Mz_u_p << endln;
    opserr << "Mz- = " << Mz_u_n << endln;
#endif

    // evalaute the new strings
    std::string initialStiffnessExpressionY = rawInitialStiffnessExpressionY;
    std::string initialStiffnessExpressionZ = rawInitialStiffnessExpressionZ;
    std::string thetaPExpressionY = rawThetaPExpressionY;
    std::string thetaPExpressionZ = rawThetaPExpressionZ;
    std::string thetaPCExpressionY = rawThetaPCExpressionY;
    std::string thetaPCExpressionZ = rawThetaPCExpressionZ;

    // Read the actual values
    N = axialMaterial->getStress();
    My = MyMaterial->getStress();
    Mz = MzMaterial->getStress();
    Vy = VyMaterial->getStress();
    Vz = VzMaterial->getStress();

    // substitute values
    replacePlaceholderWithValue(initialStiffnessExpressionY, "__N__", N);
    replacePlaceholderWithValue(initialStiffnessExpressionZ, "__N__", N);
    replacePlaceholderWithValue(initialStiffnessExpressionY, "__M__", My);
    replacePlaceholderWithValue(initialStiffnessExpressionZ, "__M__", Mz);
    replacePlaceholderWithValue(initialStiffnessExpressionY, "__V__", Vz);
    replacePlaceholderWithValue(initialStiffnessExpressionZ, "__V__", Vy);
    replacePlaceholderWithValue(initialStiffnessExpressionY, "__My__", My);
    replacePlaceholderWithValue(initialStiffnessExpressionZ, "__My__", My);
    replacePlaceholderWithValue(initialStiffnessExpressionY, "__Vz__", Vz);
    replacePlaceholderWithValue(initialStiffnessExpressionZ, "__Vz__", Vz);
    replacePlaceholderWithValue(initialStiffnessExpressionY, "__Mz__", Mz);
    replacePlaceholderWithValue(initialStiffnessExpressionZ, "__Mz__", Mz);
    replacePlaceholderWithValue(initialStiffnessExpressionY, "__Vy__", Vy);
    replacePlaceholderWithValue(initialStiffnessExpressionZ, "__Vy__", Vy);
    replacePlaceholderWithValue(thetaPExpressionY, "__N__", N);
    replacePlaceholderWithValue(thetaPExpressionZ, "__N__", N);
    replacePlaceholderWithValue(thetaPExpressionY, "__M__", My);
    replacePlaceholderWithValue(thetaPExpressionZ, "__M__", Mz);
    replacePlaceholderWithValue(thetaPExpressionY, "__V__", Vz);
    replacePlaceholderWithValue(thetaPExpressionZ, "__V__", Vy);
    replacePlaceholderWithValue(thetaPExpressionY, "__My__", My);
    replacePlaceholderWithValue(thetaPExpressionZ, "__My__", My);
    replacePlaceholderWithValue(thetaPExpressionY, "__Vz__", Vz);
    replacePlaceholderWithValue(thetaPExpressionZ, "__Vz__", Vz);
    replacePlaceholderWithValue(thetaPExpressionY, "__Mz__", Mz);
    replacePlaceholderWithValue(thetaPExpressionZ, "__Mz__", Mz);
    replacePlaceholderWithValue(thetaPExpressionY, "__Vy__", Vy);
    replacePlaceholderWithValue(thetaPExpressionZ, "__Vy__", Vy);
    replacePlaceholderWithValue(thetaPCExpressionY, "__N__", N);
    replacePlaceholderWithValue(thetaPCExpressionZ, "__N__", N);
    replacePlaceholderWithValue(thetaPCExpressionY, "__M__", My);
    replacePlaceholderWithValue(thetaPCExpressionZ, "__M__", Mz);
    replacePlaceholderWithValue(thetaPCExpressionY, "__V__", Vz);
    replacePlaceholderWithValue(thetaPCExpressionZ, "__V__", Vy);
    replacePlaceholderWithValue(thetaPCExpressionY, "__My__", My);
    replacePlaceholderWithValue(thetaPCExpressionZ, "__My__", My);
    replacePlaceholderWithValue(thetaPCExpressionY, "__Vz__", Vz);
    replacePlaceholderWithValue(thetaPCExpressionZ, "__Vz__", Vz);
    replacePlaceholderWithValue(thetaPCExpressionY, "__Mz__", Mz);
    replacePlaceholderWithValue(thetaPCExpressionZ, "__Mz__", Mz);
    replacePlaceholderWithValue(thetaPCExpressionY, "__Vy__", Vy);
    replacePlaceholderWithValue(thetaPCExpressionZ, "__Vy__", Vy);

    // update My material
    double EI;
    if (OPS_EvalDoubleStringExpression(initialStiffnessExpressionY.c_str(), EI) < 0) {
        opserr << "\n\nError evaluating expression:\n" << initialStiffnessExpressionY.c_str() << "\n";
        opserr << "Raw string: " << rawInitialStiffnessExpressionY.c_str() << "\n";
        opserr << "N Value = " << N << "\n";
        opserr << "Replaced String: " << initialStiffnessExpressionY.c_str() << "\n";
        opserr << "Replacing again the string: \n";
        replacePlaceholderWithValue(initialStiffnessExpressionY, "__N__", N);
    }
    double th_cap_p;
    if (OPS_EvalDoubleStringExpression(thetaPExpressionY.c_str(), th_cap_p) < 0) {
        opserr << "Error evaluating expression:\n" << thetaPExpressionY.c_str() << "\n";
    }
    double th_pcap_p;
    if (OPS_EvalDoubleStringExpression(thetaPCExpressionY.c_str(), th_pcap_p) < 0) {
        opserr << "Error evaluating expression:\n" << thetaPCExpressionY.c_str() << "\n";
    }

#ifdef _DBG_COUPLEDSEC3D
    opserr << "Initial stiffness evaluation:\n" << rawInitialStiffnessExpressionY.c_str() << endln;
    opserr << initialStiffnessExpressionY.c_str() << " = " << EI << endln;
    opserr << thetaPExpressionY.c_str() << " = " << th_cap_p << endln;
    opserr << thetaPCExpressionY.c_str() << " = " << th_pcap_p << endln;
#endif

    // My material
    par_f1p_Y.update(My_u_p / 5.0);
    par_f1n_Y.update(My_u_n / 5.0);
    par_d1p_Y.update((My_u_p / EI) / 5.0);
    par_d1n_Y.update((My_u_n / EI) / 5.0);
    par_f2p_Y.update(My_u_p);
    par_f2n_Y.update(My_u_n);
    par_d2p_Y.update(My_u_p / EI);
    par_d2n_Y.update(My_u_n / EI);
    par_f3p_Y.update(My_u_p * a_s);
    par_f3n_Y.update(My_u_n * a_s);
    par_d3p_Y.update(par_d2p_Y.getValue() + th_cap_p);
    par_d3n_Y.update(par_d2n_Y.getValue() - th_cap_p);
    par_f4p_Y.update(0.1 * My_u_p);
    par_f4n_Y.update(0.1 * My_u_n);
    // Be sure that point 4 is after point 3 otherwise Pinching4 works anyway but with unexpected behavior
    double d4;
    d4 = par_d2p_Y.getValue() + th_pcap_p;
    d4 = d4 > par_d3p_Y.getValue() ? d4 : par_d3p_Y.getValue() * 1.1;
    par_d4p_Y.update(d4);
    d4 = par_d2n_Y.getValue() - th_pcap_p;
    d4 = d4 < par_d3n_Y.getValue() ? d4 : par_d3n_Y.getValue() * 1.1;
    par_d4n_Y.update(d4);

#ifdef _DBG_COUPLEDSEC3D
    opserr << "\nLaw MyMaterial updated ***:\n";
    MyMaterial->Print(opserr, 2);
#endif

    // update Mz material
    if (OPS_EvalDoubleStringExpression(initialStiffnessExpressionZ.c_str(), EI) < 0) {
        opserr << "Error evaluating expression:\n" << initialStiffnessExpressionZ.c_str() << "\n";
    }
    if (OPS_EvalDoubleStringExpression(thetaPExpressionZ.c_str(), th_cap_p) < 0) {
        opserr << "Error evaluating expression:\n" << thetaPExpressionZ.c_str() << "\n";
    }
    if (OPS_EvalDoubleStringExpression(thetaPCExpressionZ.c_str(), th_pcap_p) < 0) {
        opserr << "Error evaluating expression:\n" << thetaPCExpressionZ.c_str() << "\n";
    }
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Initial stiffness evaluation:\n" << rawInitialStiffnessExpressionZ.c_str() << endln;
    opserr << initialStiffnessExpressionZ.c_str() << " = " << EI << endln;
    opserr << thetaPExpressionZ.c_str() << " = " << th_cap_p << endln;
    opserr << thetaPCExpressionZ.c_str() << " = " << th_pcap_p << endln;
#endif

    // Mz material
    par_f1p_Z.update(Mz_u_p / 5.0);
    par_f1n_Z.update(Mz_u_n / 5.0);
    par_d1p_Z.update((Mz_u_p / EI) / 5.0);
    par_d1n_Z.update((Mz_u_n / EI) / 5.0);
    par_f2p_Z.update(Mz_u_p);
    par_f2n_Z.update(Mz_u_n);
    par_d2p_Z.update(Mz_u_p / EI);
    par_d2n_Z.update(Mz_u_n / EI);
    par_f3p_Z.update(Mz_u_p * a_s);
    par_f3n_Z.update(Mz_u_n * a_s);
    par_d3p_Z.update(par_d2p_Z.getValue() + th_cap_p);
    par_d3n_Z.update(par_d2n_Z.getValue() - th_cap_p);
    par_f4p_Z.update(0.1 * Mz_u_p);
    par_f4n_Z.update(0.1 * Mz_u_n);
    // Be sure that point 4 is after point 3 otherwise Pinching4 works anyway but with unexpected behavior
    d4 = par_d2p_Z.getValue() + th_pcap_p;
    d4 = d4 > par_d3p_Z.getValue() ? d4 : par_d3p_Z.getValue() * 1.1;
    par_d4p_Z.update(d4);
    d4 = par_d2n_Z.getValue() - th_pcap_p;
    d4 = d4 < par_d3n_Z.getValue() ? d4 : par_d3n_Z.getValue() * 1.1;
    par_d4n_Z.update(d4);

#ifdef _DBG_COUPLEDSEC3D
    opserr << "\nLaw MzMaterial updated ***:\n";
    MzMaterial->Print(opserr, 2);
#endif
}

int
ASDCoupledHinge3D::commitState(void)
{
    int err = 0;

    err += axialMaterial->commitState();
    err += MyMaterial->commitState();
    err += MzMaterial->commitState();
    err += VyMaterial->commitState();
    err += VzMaterial->commitState();
    err += torsionMaterial->commitState();

#ifdef _DBG_COUPLEDSEC3D
    opserr << "ASDCoupleHinge3D::commitState\n\t" << "N = " << axialMaterial->getStress() << endln;
#endif
    if (err == 0) {
        updateLaws();

#ifdef _DBG_COUPLEDSEC3D
        opserr << "\nLaws updated after commit: " << endln;
        opserr << "My\n";
        MyMaterial->Print(opserr, 2);
        opserr << "Mz\n";
        MzMaterial->Print(opserr, 2);
        opserr << "Finished update " << endln;
#endif

        //double My, Mz;
        //My = MyMaterial->getStress();
        //Mz = MzMaterial->getStress();

        //MyMaterial->setTrialStrain(MyMaterial->getStrain());
        //MzMaterial->setTrialStrain(MzMaterial->getStrain());

        //// compare My with MyMaterial->getStress();
        //double errY = abs(MyMaterial->getStress() - My);
        //double errZ = abs(MzMaterial->getStress() - Mz);

        //errExplicit = std::max(errY, errZ) / MmaxAbs;

        //// revert to last commit
        //err += MyMaterial->revertToLastCommit();
        //err += MzMaterial->revertToLastCommit();
#ifdef _DBG_COUPLEDSEC3D
        opserr << "\nReverted to Last commit\n";
        opserr << "My\n";
        MyMaterial->Print(opserr, 2);
        opserr << "Mz\n";
        MzMaterial->Print(opserr, 2);
#endif
    } 
#ifdef _DBG_COUPLEDSEC3D
    opserr << "ASDCoupledHinge3D::commitState (updated laws)\n";
    opserr << "My tangent: " << m_num_tang(1);
    opserr << " Mz tangent: " << m_num_tang(2) << endln;
    opserr << "-- My law: " << endln;
    MyMaterial->Print(opserr, 3);
    opserr << "\n";
    opserr << "-- Mz law: " << endln;
    MzMaterial->Print(opserr, 3);
    opserr << "\n";
#endif
    return err;
}

int
ASDCoupledHinge3D::revertToLastCommit(void)
{
    int err = 0;

    err += axialMaterial->revertToLastCommit();
    err += MyMaterial->revertToLastCommit();
    err += MzMaterial->revertToLastCommit();
    err += VyMaterial->revertToLastCommit();
    err += VzMaterial->revertToLastCommit();
    err += torsionMaterial->revertToLastCommit();

#ifdef _DBG_COUPLEDSEC3D
    opserr << "ASDCoupleHinge3D::revertToLastCommit\n\t" << "N = " << axialMaterial->getStress() << endln;
#endif
    if (err == 0) {
        updateLaws();
    }

    

    return err;
}	

int
ASDCoupledHinge3D::revertToStart(void)
{
    int err = 0;

    err += axialMaterial->revertToStart();
    err += MyMaterial->revertToStart();
    err += MzMaterial->revertToStart();
    err += VyMaterial->revertToStart();
    err += VzMaterial->revertToStart();
    err += torsionMaterial->revertToStart();

#ifdef _DBG_COUPLEDSEC3D
    opserr << "ASDCoupleHinge3D::revertToStart\n\t" << "N = " << axialMaterial->getStress() << endln;
#endif
    if (err == 0) {
        updateLaws();
    }

    return err;
}
//
//
int
ASDCoupledHinge3D::sendSelf(int cTag, Channel &theChannel)
{
    // my db tag
    if (otherDbTag == 0)
        otherDbTag = theChannel.getDbTag();

    // materials' class and db tags
    std::vector<int> sub_classes(6);
    std::vector<int> sub_dbtags(6);
    std::array<UniaxialMaterial*, 6> subs = { axialMaterial, MyMaterial, MzMaterial, torsionMaterial, VyMaterial, VzMaterial };
    for (int i = 0; i < subs.size(); ++i) {
        sub_classes[i] = subs[i]->getClassTag();
        sub_dbtags[i] = subs[i]->getDbTag();
        if (sub_dbtags[i] == 0) {
            sub_dbtags[i] = theChannel.getDbTag();
            if (sub_dbtags[i] != 0)
                subs[i]->setDbTag(sub_dbtags[i]);
        }
    }

    // serializer
    Serializer ser;

    // serialize everything
    if (!(ser
        // misc info
        << getTag()
        << otherDbTag
        << sub_classes
        << sub_dbtags
        // settings
        << tolN
        << tolM
        << MmaxAbs
        << errExplicit
        // strength Domain
        << strengthDomain.size
        << strengthDomain.numberAxial
        << strengthDomain.numberTheta
        << strengthDomain.numberData
        << strengthDomain.theVector
        // strings
        << rawInitialStiffnessExpressionY
        << rawInitialStiffnessExpressionZ
        << rawThetaPExpressionY
        << rawThetaPExpressionZ
        << rawThetaPCExpressionY
        << rawThetaPCExpressionZ
        << a_s
        ))
    {
        opserr << "ASDCoupledHinge3D::sendSelf() - failed to serialize data\n";
        return -1;
    }

    // get message string and size
    std::string msg_string = ser.str();
    std::vector<char> msg_data(msg_string.size() + 1);
    std::copy(msg_string.begin(), msg_string.end(), msg_data.begin());
    msg_data.back() = '\0';
    int msg_data_size = static_cast<int>(msg_string.size());

    // send message size
    ID idata(1);
    idata(0) = msg_data_size;
    if (theChannel.sendID(otherDbTag, cTag, idata) < 0) {
        opserr << "ASDCoupledHinge3D::sendSelf() - failed to send message size\n";
        return -1;
    }

    // send message
    Message msg(msg_data.data(), msg_data_size);
    if (theChannel.sendMsg(otherDbTag, cTag, msg) < 0) {
        opserr << "ASDCoupledHinge3D::sendSelf() - failed to send message\n";
        return -1;
    }

    // send materials
    for (int i = 0; i < subs.size(); ++i) {
        if (subs[i]->sendSelf(cTag, theChannel) < 0) {
            opserr << "ASDCoupledHinge3D::sendSelf() - failed to send material " << i << "\n";
            return -1;
        }
    }

    // done
    return 0;
}


int
ASDCoupledHinge3D::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // recv message size
    ID idata(1);
    if (theChannel.recvID(0, cTag, idata) < 0) {
        opserr << "ASDCoupledHinge3D::recvSelf() - failed to recv message size\n";
        return -1;
    }
    int msg_data_size = idata(0);

    // recv message
    std::vector<char> msg_data(static_cast<size_t>(msg_data_size) + 1);
    Message msg(msg_data.data(), msg_data_size);
    if (theChannel.recvMsg(0, cTag, msg) < 0) {
        opserr << "ASDCoupledHinge3D::recvSelf() - failed to recv message\n";
        return -1;
    }
    msg_data.back() = '\0';

    // serializer
    Serializer ser(msg_data.data());

    // aux data for de-serialziation
    int my_tag;
    // materials' class and db tags
    std::vector<int> sub_classes(6);
    std::vector<int> sub_dbtags(6);

    // de-serialize everything
    if (!(ser
        // misc info
        >> my_tag
        >> otherDbTag
        >> sub_classes
        >> sub_dbtags
        // settings
        >> tolN
        >> tolM
        >> MmaxAbs
        >> errExplicit
        // strength Domain
        >> strengthDomain.size
        >> strengthDomain.numberAxial
        >> strengthDomain.numberTheta
        >> strengthDomain.numberData
        >> strengthDomain.theVector
        // strings
        >> rawInitialStiffnessExpressionY
        >> rawInitialStiffnessExpressionZ
        >> rawThetaPExpressionY
        >> rawThetaPExpressionZ
        >> rawThetaPCExpressionY
        >> rawThetaPCExpressionZ
        >> a_s
        ))
    {
        opserr << "ASDCoupledHinge3D::recvSelf() - failed to de-serialize data\n";
        return -1;
    }

    // set tag
    setTag(my_tag);

    // recv materials
    std::array<UniaxialMaterial**, 6> subs = { &axialMaterial, &MyMaterial, &MzMaterial, &torsionMaterial, &VyMaterial, &VzMaterial };
    for (int i = 0; i < subs.size(); ++i) {
        UniaxialMaterial *m = theBroker.getNewUniaxialMaterial(sub_classes[i]);
        m->setDbTag(sub_dbtags[i]);
        if (m->recvSelf(cTag, theChannel, theBroker) < 0) {
            opserr << "ASDCoupledHinge3D::recvSelf() - failed to recv material " << i << "\n";
            return -1;
        }
        *(subs[i]) = m;
    }

    // setup parameters
    setupParameters();

    // done
    return 0;
}

void
ASDCoupledHinge3D::Print(OPS_Stream &s, int flag)
{
    s << "ASDCoupledHinge3D" << endln;
}

int
ASDCoupledHinge3D::getVariable(const char *argv, Information &info)
{
    info.theDouble = 0.0;
    const Vector& e = this->getSectionDeformation();
    const ID& code = this->getType();
    if (strcmp(argv, "axialStrain") == 0) {
        for (int i = 0; i < 6; i++)
            if (code(i) == SECTION_RESPONSE_P)
                info.theDouble += e(i);
    }
    else if (strcmp(argv, "curvatureZ") == 0) {
        for (int i = 0; i < 6; i++)
            if (code(i) == SECTION_RESPONSE_MZ)
                info.theDouble += e(i);
    }
    else if (strcmp(argv, "curvatureY") == 0) {
        for (int i = 0; i < 6; i++)
            if (code(i) == SECTION_RESPONSE_MY)
                info.theDouble += e(i);
    }
    else
        return -1;
    return 0;
}

void ASDCoupledHinge3D::resetStrengthDomain(double& My_u_p, double& My_u_n, double& Mz_u_p, double& Mz_u_n)
{
    setUncoupledStrengthDomainforAxial(0.0, My_u_p, My_u_n, Mz_u_p, Mz_u_n);

    return;
}

void ASDCoupledHinge3D::setUncoupledStrengthDomainforAxial(const double N, double& My_u_p, double& My_u_n, double& Mz_u_p, double& Mz_u_n)
{
    double tmp;
    strengthDomain.getMyMzForNAndDirection(N, 0.0, My_u_p, tmp);
    strengthDomain.getMyMzForNAndDirection(N, M_PI, My_u_n, tmp);
    strengthDomain.getMyMzForNAndDirection(N, M_PI / 2, tmp, Mz_u_p);
    strengthDomain.getMyMzForNAndDirection(N, M_PI * 3 / 2.0, tmp, Mz_u_n);
}

void ASDCoupledHinge3D::setupParameters()
{
    // create the parameters for My and Mz material: 
    const char* argv[] = { "par" };
    // My material
    // Point 1 positive
    argv[0] = "f1p";
    MyMaterial->setParameter(argv, 1, par_f1p_Y);
    argv[0] = "d1p";
    MyMaterial->setParameter(argv, 1, par_d1p_Y);
    // Point 2 postive
    argv[0] = "f2p";
    MyMaterial->setParameter(argv, 1, par_f2p_Y);
    argv[0] = "d2p";
    MyMaterial->setParameter(argv, 1, par_d2p_Y);
    // Point 3 positive
    argv[0] = "f3p";
    MyMaterial->setParameter(argv, 1, par_f3p_Y);
    argv[0] = "d3p";
    MyMaterial->setParameter(argv, 1, par_d3p_Y);
    // Point 4 postive
    argv[0] = "f4p";
    MyMaterial->setParameter(argv, 1, par_f4p_Y);
    argv[0] = "d4p";
    MyMaterial->setParameter(argv, 1, par_d4p_Y);
    // Point 1 negative
    argv[0] = "f1n";
    MyMaterial->setParameter(argv, 1, par_f1n_Y);
    argv[0] = "d1n";
    MyMaterial->setParameter(argv, 1, par_d1n_Y);
    // Point 2 negative
    argv[0] = "f2n";
    MyMaterial->setParameter(argv, 1, par_f2n_Y);
    argv[0] = "d2n";
    MyMaterial->setParameter(argv, 1, par_d2n_Y);
    // Point 3 negative
    argv[0] = "f3n";
    MyMaterial->setParameter(argv, 1, par_f3n_Y);
    argv[0] = "d3n";
    MyMaterial->setParameter(argv, 1, par_d3n_Y);
    // Point 4 negative
    argv[0] = "f4n";
    MyMaterial->setParameter(argv, 1, par_f4n_Y);
    argv[0] = "d4n";
    MyMaterial->setParameter(argv, 1, par_d4n_Y);
    // Mz material
    // Point 1 positive
    argv[0] = "f1p";
    MzMaterial->setParameter(argv, 1, par_f1p_Z);
    argv[0] = "d1p";
    MzMaterial->setParameter(argv, 1, par_d1p_Z);
    // Point 2 postive
    argv[0] = "f2p";
    MzMaterial->setParameter(argv, 1, par_f2p_Z);
    argv[0] = "d2p";
    MzMaterial->setParameter(argv, 1, par_d2p_Z);
    // Point 3 positive
    argv[0] = "f3p";
    MzMaterial->setParameter(argv, 1, par_f3p_Z);
    argv[0] = "d3p";
    MzMaterial->setParameter(argv, 1, par_d3p_Z);
    // Point 4 postive
    argv[0] = "f4p";
    MzMaterial->setParameter(argv, 1, par_f4p_Z);
    argv[0] = "d4p";
    MzMaterial->setParameter(argv, 1, par_d4p_Z);
    // Point 1 negative
    argv[0] = "f1n";
    MzMaterial->setParameter(argv, 1, par_f1n_Z);
    argv[0] = "d1n";
    MzMaterial->setParameter(argv, 1, par_d1n_Z);
    // Point 2 negative
    argv[0] = "f2n";
    MzMaterial->setParameter(argv, 1, par_f2n_Z);
    argv[0] = "d2n";
    MzMaterial->setParameter(argv, 1, par_d2n_Z);
    // Point 3 negative
    argv[0] = "f3n";
    MzMaterial->setParameter(argv, 1, par_f3n_Z);
    argv[0] = "d3n";
    MzMaterial->setParameter(argv, 1, par_d3n_Z);
    // Point 4 negative
    argv[0] = "f4n";
    MzMaterial->setParameter(argv, 1, par_f4n_Z);
    argv[0] = "d4n";
    MzMaterial->setParameter(argv, 1, par_d4n_Z);
}

Response* ASDCoupledHinge3D::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    // check responses for this class
    if (argc > 0) {
        if (strcmp(argv[0], "explicitError") == 0) {
            output.tag("SectionOutput");
            output.attr("secType", this->getClassType());
            output.attr("secTag", this->getTag());
            output.tag("ResponseType", "error");
            double err = 0.0;
            Response *theResponse = new MaterialResponse(this, 100, err);
            output.endTag();
            return theResponse;
        }
    }
    // otherwise fallbackt to baseclass
    return SectionForceDeformation::setResponse(argv, argc, output);
}

int ASDCoupledHinge3D::getResponse(int responseID, Information& sectInfo)
{
    if (responseID == 100) {
        return sectInfo.setDouble(errExplicit);
    }

    return SectionForceDeformation::getResponse(responseID, sectInfo);
}
