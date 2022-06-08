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

#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <classTags.h>
#include <ASDCoupledHinge3D.h>
#include <MaterialResponse.h>
#include <ID.h>
#include <math.h>

#include <classTags.h>
#include <elementAPI.h>
#include <Pinching4Material.h>
#include <ElasticMaterial.h>

//#define _DBG_COUPLEDSEC3D

const double DomainData::pi = acos(-1.0);

DomainData::DomainData() :
    numberAxial(0), numberTheta(0), numberData(0), size(0)
{
    theVector = new Vector(size);
}

DomainData::DomainData(int nN, int nTheta, int nData) :
    numberAxial(nN), numberTheta(nTheta), numberData(nData)
{
    size = numberAxial * numberTheta * numberData;
    theVector = new Vector(size);
    theVector->Zero();
}

DomainData::~DomainData()
{
    delete theVector;
}

double DomainData::getValue(int i, int j, int k)
{
    int idx = i * (numberTheta * numberData) + k * numberTheta + j;
    return (*theVector)(idx);
}

void DomainData::setValue(int i, int j, int k, double val)
{
    int idx = i * (numberTheta * numberData) + k * numberTheta + j;
    (*theVector)(idx) = val;
}

void DomainData::getRangeN(double& Nmin, double& Nmax)
{
    Nmin = this->getValue(0, 0, 0);
    Nmax = this->getValue(numberAxial - 1, 0, 0);
}

void DomainData::print(void)
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

DomainData* DomainData::getCopy(void)
{
    DomainData* theCopy = new DomainData(numberAxial, numberTheta, numberData);
    // Deep copy of the data
    for (int i = 0; i < numberAxial; i++) 
    {
        for (int j = 0; j < numberTheta; j++) 
        {
            for (int k = 0; k < numberData; k++) 
            {
//#ifdef _DBG_COUPLEDSEC3D
//                opserr << "copying value (i=" << i << ", j=" << j << ", k=" << k << ") val = " << this->getValue(i, j, k) << " into new domain\n";
//#endif
                theCopy->setValue(i, j, k, this->getValue(i, j, k));
            }
        }
    }
//#ifdef _DBG_COUPLEDSEC3D
//    opserr << "Original domain\n";
//    this->print();
//    opserr << "Copied domain\n";
//    theCopy->print();
//#endif
    return theCopy;
}

int DomainData::getMyMzForNAndDirection(double N, double theta, double& My, double& Mz) {
    // transform theta in the range [0, 2pi[
    while (theta > 2 * pi)
        theta -= 2 * pi;
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
    while ((!found) && (j <= numberTheta)) 
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
            theta_j += 2 * pi;
        dir_jp1(0) = My_jp1;
        dir_jp1(1) = Mz_jp1;
        dir_jp1.Normalize();
        theta_jp1 = atan2(Mz_jp1, My_jp1);
        if (theta_jp1 < 0)
            theta_jp1 += 2 * pi;
#ifdef _DBG_COUPLEDSEC3D
        opserr << "j = " << j << " - jp1 = " << jp1 << "\n";
        opserr << "theta_j = " << theta_j <<  " - theta_jp1 = " << theta_jp1 << " - theta = " << theta << "\n";
        opserr << "direction j = [" << dir_j(0) << ", " << dir_j(1) << "]\n";
        opserr << "direction j+1 = [" << dir_jp1(0) << ", " << dir_jp1(1) << "]\n";

#endif
        double jXjp1 = dir_j(0) * dir_jp1(1) - dir_j(1) * dir_jp1(0);
        double jXdirection = dir_j(0) * direction(1) - dir_j(1) * direction(0);
        double directionXjp1 = direction(0) * dir_jp1(1) - direction(1) * dir_jp1(0);
        if ((jXjp1 * jXdirection > 0) && (directionXjp1 * jXjp1 > 0))
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
    if (inverse == 0) 
    {
        My = (My_jp1 - My_j) / (theta_jp1 - theta_j) * (theta - theta_j) + My_j;
        Mz = (Mz_jp1 - Mz_j) / (theta_jp1 - theta_j) * (theta - theta_j) + Mz_j;
    } else 
    {
        My = (My_j - My_jp1) / (theta_j - theta_jp1) * (theta - theta_jp1) + My_jp1;
        Mz = (Mz_j - Mz_jp1) / (theta_j - theta_jp1) * (theta - theta_jp1) + Mz_jp1;
    }

#ifdef _DBG_COUPLEDSEC3D
    opserr << "Interpolated values:\n";
    opserr << "My = " << My << " Mz = " << Mz << "\n";
#endif

    return 0;
}

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

void* OPS_ASDCoupledHinge3D()
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDCoupledHinge3D - Developed by: Diego Talledo, Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARNING insufficient arguments\n"
            "Want: section ASDCoupledHinge3D tag? Ktor? Kvy? Kvz? Kax? -initialFlexuralStiffness \"Ky?\" \"Kz?\""
            "<-simpleStrengthDomain Nmin? Nmax? MyMax? NMyMax? MzMax? NMzMax?> <-strengthDomainByPoints nN? nTheta? {listN} {listMy} {listMz}> <-hardening as?>"
            "-thetaP \"thetaPy?\" \"thetaPz?\" -thetaPC \"thetaPCy?\" \"thetaPCz?\"\n";
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
    double NMyMax = 0;
    double NMzMax = 0;
    // by points
    // number of discretizations in N axis
    // number of discretization about angle in My-Mz plane
    int nN = 0;
    int nTheta = 0;
    int n = 0;
    Vector dataN(0);
    Vector dataMy(0);
    Vector dataMz(0);
    DomainData* ultDomain = nullptr;
    // Domain Values for My and Mz positive and negative when N = 0
    double My_u_p = 0.0;
    double My_u_n = 0.0;
    double Mz_u_p = 0.0;
    double Mz_u_n = 0.0;

    // hardening
    double a_s = 1.001;

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
        }
        else if (strcmp(type, "-thetaPC") == 0) {
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
            //opserr << "Parsing hardening\n";
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
        } else if (strcmp(type, "-simpleStrengthDomain") == 0) {
            if (strengthDomainMode != STRENGTH_DOMAIN_UNDEFINED) {
                opserr << "only one strength domain need to be defined!\n"
                    "section ASDCoupledHinge3D: " << tag << endln;
                return 0;
            }
            // simplified strength domain, giving few data.
            // we assume a double symetric - rect - section so that MinMoment = - MaxMoment
            numdata = 6;
            if (OPS_GetNumRemainingInputArgs() < numdata) {
                opserr << "WARNING insufficient arguments\n"
                    "-simpleStrengthDomain Nmin? Nmax? MyMax? NMyMax? MzMax? NMzMax?"
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
            NMyMax = data[4];
            NMzMax = data[5];

            strengthDomainMode = STRENGTH_DOMAIN_SIMPLIFIED;

            // TO DO: check correctness of data and then create a simplified domain with some sort of interpolation

            opserr << "To be implemented !!! ERROR\n";
            return 0;
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

            //for (int i = 0; i < nN; i++) {
            //    int idx;
            //    for (int j = 0; j < nTheta; j++) {
            //        idx = i * nTheta + j;
            //        opserr << "dataN(idx): idx = " << idx << " - dataN(idx) = " << dataN(idx) << "\n";
            //        opserr << "dataMy(idx): idx = " << idx << " - dataMy(idx) = " << dataMy(idx) << "\n";
            //        opserr << "dataMz(idx): idx = " << idx << " - dataMz(idx) = " << dataMz(idx) << "\n";
            //        opserr << "\n";
            //    }
            //    opserr << "\n\n";
            //}

            // Create the domain as a DomainData object
            ultDomain = new DomainData(nN, nTheta, 3);

            // write all elements of the domain data
            int idx;
            for (int i = 0; i < nN; i++) {
                for (int j = 0; j < nTheta; j++) {
                    idx = i * nTheta + j;
                    ultDomain->setValue(i, j, 0, dataN(idx));
                    ultDomain->setValue(i, j, 1, dataMy(idx));
                    ultDomain->setValue(i, j, 2, dataMz(idx));
                }
            }
            //opserr << "ultimateDomain created: " << endln;
            //ultDomain->print();

            // get Domain Values for My and Mz positive and negative when N = 0
            double tmp;
            //opserr << "\n\nFinding for theta = 0 (should be positive My)\n";
            ultDomain->getMyMzForNAndDirection(0.0,0.0,My_u_p,tmp);
            //opserr << "\n\n\nFinding for theta = pi (should be negative My)\n";
            ultDomain->getMyMzForNAndDirection(0.0, DomainData::pi, My_u_n, tmp);
            //opserr << "\n\n\nFinding for theta = pi/2 (should be positive Mz)\n";
            ultDomain->getMyMzForNAndDirection(0.0, DomainData::pi/2, tmp, Mz_u_p);
            //opserr << "\n\n\nFinding for theta = 3pi/2 (should be negative Mz)\n";
            ultDomain->getMyMzForNAndDirection(0.0, DomainData::pi*3/2.0, tmp, Mz_u_n);
            
            /*My_u_n = 0.0;
            Mz_u_p = 0.0;
            Mz_u_n = 0.0;*/

#ifdef _DBG_COUPLEDSEC3D
            opserr << "For N = 0 \nMy+ = " << My_u_p << endln;
            opserr << "My- = " << My_u_n << endln;
            opserr << "Mz+ = " << Mz_u_p << endln;
            opserr << "Mz- = " << Mz_u_n << endln;

            opserr << "Domain is this one: \n";
            ultDomain->print();
            opserr << endln;
#endif

            strengthDomainMode = STRENGTH_DOMAIN_BY_POINTS;
        }
        
    }

    opserr << "Parsed the command. Creating the materials.\n";

    // Create an elastic material for axial
    UniaxialMaterial* matAxial = new ElasticMaterial(0, Kax, 0.0, Kax);
    if (matAxial == 0) {
        opserr << "Error TO DO" << endln;
        return 0;
    }
    // Create an elastic material for shear_y
    UniaxialMaterial* matShearY = new ElasticMaterial(0, Kv_y, 0.0, Kv_y);
    if (matShearY == 0) {
        opserr << "Error TO DO" << endln;
        return 0;
    }
    // Create an elastic material for shear_z
    UniaxialMaterial* matShearZ = new ElasticMaterial(0, Kv_z, 0.0, Kv_z);
    if (matShearZ == 0) {
        opserr << "Error TO DO" << endln;
        return 0;
    }
    // Create an elastic material for shear_z
    UniaxialMaterial* matTorsion = new ElasticMaterial(0, Ktor, 0.0, Ktor);
    if (matTorsion == 0) {
        opserr << "Error TO DO" << endln;
        return 0;
    }
    opserr << "Materials for Torsion, Axial and Shear (Y and Z) defined\n";

    // Create pinching4 material for My
    double EI;
    if (OPS_EvalDoubleTclStringExpression(theInitialStiffnessExpressionY.c_str(), EI) < 0) {
        opserr << "Error evaluating expression:\n" << theInitialStiffnessExpressionY.c_str() << "\n";
    }
    opserr << theInitialStiffnessExpressionY.c_str() << " = " << EI << endln;
    double th_cap_p;
    if (OPS_EvalDoubleTclStringExpression(theThetaPExpressionY.c_str(), th_cap_p) < 0) {
        opserr << "Error evaluating expression:\n" << theThetaPExpressionY.c_str() << "\n";
    }
    opserr << theThetaPExpressionY.c_str() << " = " << th_cap_p << endln;
    double th_pcap_p;
    if (OPS_EvalDoubleTclStringExpression(theThetaPCExpressionY.c_str(), th_pcap_p) < 0) {
        opserr << "Error evaluating expression:\n" << theThetaPCExpressionY.c_str() << "\n";
    }
    opserr << theThetaPCExpressionY.c_str() << " = " << th_pcap_p << endln;
#ifdef _DBG_COUPLEDSEC3D
#endif
    // Yield strength
    double f2p = My_u_p;
    double d2p = My_u_p / EI;
    // First point in the same line (arbitraly choosen 1/5)
    double f1p = f2p / 5.0;
    double d1p = d2p / 5.0;
    // Hardening point - cap
    double f3p = a_s * f2p;
    double d3p = d2p + th_cap_p;
    // Ultimate point
    double f4p = 0.1 * My_u_p; // To be choosen by user later ???
    double d4p = d2p + th_pcap_p;
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
    double d3n = d2n - th_cap_p;
    // Ultimate point
    double f4n = 0.1 * My_u_n; // To be choosen by user later ???
    double d4n = d2n - th_pcap_p;
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Initial spring My - : (" << f1n << ", " << d1n << "), (" << f2n << ", " << d2n << "), (" << f3n << ", " << d3n << "), (" << f4n << ", " << d4n << ")\n";
#endif

    // Create a Pinching4 material for My

    double mdp = 0.5;
    double mfp = 0.3;
    double msp = 0.0;
    double mdn = 0.5;
    double mfn = 0.3;
    double msn = 0.0;
    double gk1 = 0.0;
    double gk2 = 0.0;
    double gk3 = 0.0;
    double gk4 = 0.0;
    double gklim = 0.0;
    double gd1 = 0.0;
    double gd2 = 0.0;
    double gd3 = 0.0;
    double gd4 = 0.0;
    double gdlim = 0.0;
    double gf1 = 0.0;
    double gf2 = 0.0;
    double gf3 = 0.0;
    double gf4 = 0.0;
    double gflim = 0.0;
    double ge = 0.0;
    int dc = 1;
    UniaxialMaterial* matMy = new Pinching4Material(0,f1p,d1p,f2p,d2p,f3p,d3p,f4p,d4p,f1n,d1n,f2n,d2n,f3n,d3n,f4n,d4n,mdp,mfp,msp,mdn,mfn,msn,gk1,gk2,gk3,gk4,gklim,gd1,gd2,gd3,gd4,gdlim,gf1,gf2,gf3,gf4,gflim,ge,dc);
#ifdef _DBG_COUPLEDSEC3D
    matMy->Print(opserr, 2);
    opserr << endln;
#endif
    
    // Create pinching4 material for Mz
    if (OPS_EvalDoubleTclStringExpression(theInitialStiffnessExpressionZ.c_str(), EI) < 0) {
        opserr << "Error evaluating expression:\n" << theInitialStiffnessExpressionZ.c_str() << "\n";
    }
    opserr << theInitialStiffnessExpressionZ.c_str() << " = " << EI << endln;
    if (OPS_EvalDoubleTclStringExpression(theThetaPExpressionZ.c_str(), th_cap_p) < 0) {
        opserr << "Error evaluating expression:\n" << theThetaPExpressionZ.c_str() << "\n";
    }
    opserr << theThetaPExpressionZ.c_str() << " = " << th_cap_p << endln;
    if (OPS_EvalDoubleTclStringExpression(theThetaPCExpressionZ.c_str(), th_pcap_p) < 0) {
        opserr << "Error evaluating expression:\n" << theThetaPCExpressionZ.c_str() << "\n";
    }
    opserr << theThetaPCExpressionZ.c_str() << " = " << th_pcap_p << endln;
    // Yield strength
    f2p = Mz_u_p;
    d2p = Mz_u_p / EI;
    // First point in the same line (arbitraly choosen 1/5)
    f1p = f2p / 5.0;
    d1p = d2p / 5.0;
    // Hardening point - cap
    f3p = a_s * f2p;
    d3p = d2p + th_cap_p;
    // Ultimate point
    f4p = 0.1 * Mz_u_p; // To be choosen by user later ???
    d4p = d2p + th_pcap_p;
    opserr << "Initial spring Mz + : (" << f1p << ", " << d1p << "), (" << f2p << ", " << d2p << "), (" << f3p << ", " << d3p << "), (" << f4p << ", " << d4p << ")\n";
    // Yield strength
    f2n = Mz_u_n;
    d2n = Mz_u_n / EI;
    // First point in the same line (arbitraly choosen 1/5)
    f1n = f2n / 5.0;
    d1n = d2n / 5.0;
    // Hardening point - cap
    f3n = a_s * f2n;
    d3n = d2n - th_cap_p;
    // Ultimate point
    f4n = 0.1 * Mz_u_n; // To be choosen by user later ???
    d4n = d2n - th_pcap_p;
    opserr << "Initial spring Mz - : (" << f1n << ", " << d1n << "), (" << f2n << ", " << d2n << "), (" << f3n << ", " << d3n << "), (" << f4n << ", " << d4n << ")\n";

    // Create a Pinching4 material for Mz
    mdp = 0.5;
    mfp = 0.3;
    msp = 0.0;
    mdn = 0.5;
    mfn = 0.3;
    msn = 0.0;
    gk1 = 0.0;
    gk2 = 0.0;
    gk3 = 0.0;
    gk4 = 0.0;
    gklim = 0.0;
    gd1 = 0.0;
    gd2 = 0.0;
    gd3 = 0.0;
    gd4 = 0.0;
    gdlim = 0.0;
    gf1 = 0.0;
    gf2 = 0.0;
    gf3 = 0.0;
    gf4 = 0.0;
    gflim = 0.0;
    ge = 0.0;
    dc = 1;
    UniaxialMaterial* matMz = new Pinching4Material(0, f1p, d1p, f2p, d2p, f3p, d3p, f4p, d4p, f1n, d1n, f2n, d2n, f3n, d3n, f4n, d4n, mdp, mfp, msp, mdn, mfn, msn, gk1, gk2, gk3, gk4, gklim, gd1, gd2, gd3, gd4, gdlim, gf1, gf2, gf3, gf4, gflim, ge, dc);
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
    ASDCoupledHinge3D *theSection = new ASDCoupledHinge3D(tag, matTorsion, matAxial, matShearY, matShearZ, matMy, matMz, ultDomain, theRawInitialStiffnessExpressionY, theRawInitialStiffnessExpressionZ, theRawThetaPExpressionY, theRawThetaPExpressionZ, theRawThetaPCExpressionY, theRawThetaPCExpressionZ, a_s);

    // Now I delete them because I've already copied them inside the secion object.
    delete matAxial;
    delete matMy;
    delete matMz;
    delete matShearY;
    delete matShearZ;
    delete matTorsion;
    delete ultDomain;

    opserr << "Section created\n\n";

    return theSection;
}


#define maxOrder 10

// Assumes section order is less than or equal to maxOrder.
// Can increase if needed!!!
double ASDCoupledHinge3D::workArea[2*maxOrder*(maxOrder+1)];
int    ASDCoupledHinge3D::codeArea[maxOrder];

// constructors:
// (tag, matAxial, matMy, matMz, ultDomain, theRawInitialStiffnessExpressionY, theRawInitialStiffnessExpressionZ, theRawThetaPExpressionY, theRawThetaPExpressionZ, theRawThetaPCExpressionY, theRawThetaPCExpressionZ, a_s);
ASDCoupledHinge3D::ASDCoupledHinge3D(int tag, UniaxialMaterial* theTorsionMaterial, UniaxialMaterial* theAxialMaterial, UniaxialMaterial* theShearYMaterial, UniaxialMaterial* theShearZMaterial, UniaxialMaterial* theMomentYMaterial, UniaxialMaterial* theMomentZMaterial,
    DomainData* ultDomain, std::string theRawInitialStiffnessExpressionY, std::string theRawInitialStiffnessExpressionZ, std::string theRawThetaPExpressionY, std::string theRawThetaPExpressionZ, 
    std::string theRawThetaPCExpressionY, std::string theRawThetaPCExpressionZ, double a_s_i):
    SectionForceDeformation(tag, SEC_TAG_ASDCoupledHinge3D),
    matCodes(0), numMats(6), theCode(0), e(0), s(0), ks(0), fs(0), otherDbTag(0), axialMaterial(0), MyMaterial(0), MzMaterial(0),ultimateDomain(0)
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

    // copy the raw strings and the needed parameters for updating
    rawInitialStiffnessExpressionY = theRawInitialStiffnessExpressionY;
    rawInitialStiffnessExpressionZ = theRawInitialStiffnessExpressionZ;
    rawThetaPExpressionY = theRawThetaPExpressionY;
    rawThetaPExpressionZ = theRawThetaPExpressionZ;
    rawThetaPCExpressionY = theRawThetaPCExpressionY;
    rawThetaPCExpressionZ = theRawThetaPCExpressionZ;
    // hardening coefficient
    a_s = a_s_i;

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

    //// MovableObject::setParameter(const char **argv, int argc, Parameter &param)
    //opserr << "Before updating: \n";
    //MyMaterial->Print(opserr, 1);
    //par_d4p_Y.update(10.0);
    //opserr << "After updating: \n";
    //MyMaterial->Print(opserr, 1);
    //// It works! keep this!

    int order = 6;
    theCode = new ID(order);
    e = new Vector(workArea, order); // il dato è sharato
    //opserr << "Size of e " << e->Size();
    //opserr << e->operator[](0) << "\n"; // nota e[0] è il vettore intero (perché e è un array di vettori)
    //opserr << e[0][0] << "\n";
    s = new Vector(&workArea[maxOrder], order);
    ks = new Matrix(&workArea[2 * maxOrder], order, order);
    fs = new Matrix(&workArea[maxOrder * (maxOrder + 2)], order, order);
    matCodes = new ID(order);

    if (theCode == 0 || e == 0 || s == 0 || ks == 0 || fs == 0 || matCodes == 0) {
        opserr << "ASDCoupledHinge3D::ASDCoupledHinge3D " << tag << " -- out of memory\n";
        exit(-1);
    }

    // Create the codes for the section responses
    (*matCodes)(0) = SECTION_RESPONSE_P;
    (*matCodes)(1) = SECTION_RESPONSE_MY;
    (*matCodes)(2) = SECTION_RESPONSE_MZ;
    (*matCodes)(3) = SECTION_RESPONSE_VY;
    (*matCodes)(4) = SECTION_RESPONSE_VZ;
    (*matCodes)(5) = SECTION_RESPONSE_T;

    // Create and copy the new vectors for N, My, Mz of the ultimate domain
    ultimateDomain = ultDomain->getCopy();
    // Compute the needed tolerances that are relative to the domain:
    // tolN = (Nmax - Nmin)*0.001; 1 per thousand of deltaN
    double Nmin, Nmax;
    ultimateDomain->getRangeN(Nmin, Nmax);
    tolN = (Nmax - Nmin) * 0.001;
    // tolM0 = Mmax * 0.01
    double MminY, MmaxY, MminZ, MmaxZ;
    setUncoupledStrengthDomainforAxial(0.0, MmaxY, MminY, MmaxZ, MminZ);
    MmaxAbs = std::max(std::max(MmaxY, abs(MminY)), std::max(MmaxZ, abs(MminZ)));
    tolM = 0.001 * MmaxAbs;

#ifdef _DBG_COUPLEDSEC3D
    opserr << "The copied ultimate domain: \n";
    ultimateDomain->print();
#endif
}

// constructor for blank object that recvSelf needs to be invoked upon
ASDCoupledHinge3D::ASDCoupledHinge3D():
    SectionForceDeformation(0, SEC_TAG_ASDCoupledHinge3D),
    matCodes(0), numMats(3), theCode(0), e(0), s(0), ks(0), fs(0), otherDbTag(0), axialMaterial(0), MyMaterial(0), MzMaterial(0), VyMaterial(0), VzMaterial(0), torsionMaterial(0), ultimateDomain(0)
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
    if (ultimateDomain) {
        delete ultimateDomain;
    }

    if (e != 0) {
        delete e;
    }

    if (s != 0) {
        delete s;
    }

    if (ks != 0) {
        delete ks;
    }

    if (fs != 0) {
        delete fs;
    }

    if (theCode != 0) {
        delete theCode;
    }

    if (matCodes != 0) {
        delete matCodes;
    }
}

int ASDCoupledHinge3D::setTrialSectionDeformation (const Vector &def)
{
  int ret = 0;
  
  ret += axialMaterial->setTrialStrain(def(0));
  ret += MyMaterial->setTrialStrain(def(1));
  ret += MzMaterial->setTrialStrain(def(2));
  ret += VyMaterial->setTrialStrain(def(3));
  ret += VzMaterial->setTrialStrain(def(4));
  ret += torsionMaterial->setTrialStrain(def(5));
  
  return ret;
}

const Vector &
ASDCoupledHinge3D::getSectionDeformation(void)
{
    (*e)(0) = axialMaterial->getStrain();
    (*e)(1) = MyMaterial->getStrain();
    (*e)(2) = MzMaterial->getStrain();
    (*e)(3) = VyMaterial->getStrain();
    (*e)(4) = VzMaterial->getStrain();
    (*e)(5) = torsionMaterial->getStrain();

    return *e;
}

const Matrix &
ASDCoupledHinge3D::getSectionTangent(void)
{
    // Zero before assembly
    ks->Zero();

    (*ks)(0, 0) = axialMaterial->getTangent();
    (*ks)(1, 1) = MyMaterial->getTangent();
    (*ks)(2, 2) = MzMaterial->getTangent();
    (*ks)(3, 3) = VyMaterial->getTangent();
    (*ks)(4, 4) = VzMaterial->getTangent();
    (*ks)(5, 5) = torsionMaterial->getTangent();

    return *ks;
}

const Matrix &
ASDCoupledHinge3D::getInitialTangent(void)
{
    // Zero before assembly
    ks->Zero();

    (*ks)(0, 0) = axialMaterial->getInitialTangent();
    (*ks)(1, 1) = MyMaterial->getInitialTangent();
    (*ks)(2, 2) = MzMaterial->getInitialTangent();
    (*ks)(3, 3) = VyMaterial->getInitialTangent();
    (*ks)(4, 4) = VzMaterial->getInitialTangent();
    (*ks)(5, 5) = torsionMaterial->getInitialTangent();

    return *ks;
}

const Matrix &
ASDCoupledHinge3D::getSectionFlexibility(void)
{
    // Zero before assembly
    fs->Zero();
    
    double k;
    k = axialMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(0, 0) = 1.e14;
    }
    else {
        (*fs)(0, 0) = 1 / k;
    }
    k = MyMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(1, 1) = 1.e14;
    }
    else {
        (*fs)(1, 1) = 1 / k;
    }
    k = MzMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(2, 2) = 1.e14;
    }
    else {
        (*fs)(2, 2) = 1 / k;
    }
    k = VyMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(3, 3) = 1.e14;
    }
    else {
        (*fs)(3, 3) = 1 / k;
    }
    k = VzMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(4, 4) = 1.e14;
    }
    else {
        (*fs)(4, 4) = 1 / k;
    }
    k = torsionMaterial->getTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(5, 5) = 1.e14;
    }
    else {
        (*fs)(5, 5) = 1 / k;
    }
  
    return *fs;
}

const Matrix &
ASDCoupledHinge3D::getInitialFlexibility(void)
{
    // Zero before assembly
    fs->Zero();

    double k;
    k = axialMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(0, 0) = 1.e14;
    }
    else {
        (*fs)(0, 0) = 1 / k;
    }
    k = MyMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(1, 1) = 1.e14;
    }
    else {
        (*fs)(1, 1) = 1 / k;
    }
    k = MzMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(2, 2) = 1.e14;
    }
    else {
        (*fs)(2, 2) = 1 / k;
    }
    k = VyMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(3, 3) = 1.e14;
    }
    else {
        (*fs)(3, 3) = 1 / k;
    }
    k = VzMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(4, 4) = 1.e14;
    }
    else {
        (*fs)(4, 4) = 1 / k;
    }
    k = torsionMaterial->getInitialTangent();
    if (k == 0.0) {
        opserr << "ASDCoupledHinge3D::getSectionFlexibility -- singular section stiffness\n";
        (*fs)(5, 5) = 1.e14;
    }
    else {
        (*fs)(5, 5) = 1 / k;
    }

    return *fs;
}

const Vector &
ASDCoupledHinge3D::getStressResultant(void)
{
    (*s)(0) = axialMaterial->getStress();
    (*s)(1) = MyMaterial->getStress();
    (*s)(2) = MzMaterial->getStress();
    (*s)(3) = VyMaterial->getStress();
    (*s)(4) = VzMaterial->getStress();
    (*s)(5) = torsionMaterial->getStress();

    return *s;
}

SectionForceDeformation *
ASDCoupledHinge3D::getCopy(void)
{
  ASDCoupledHinge3D *theCopy = 0;

  theCopy = new ASDCoupledHinge3D(this->getTag(), torsionMaterial, axialMaterial, VyMaterial, VzMaterial, MyMaterial, MzMaterial, ultimateDomain, rawInitialStiffnessExpressionY, rawInitialStiffnessExpressionZ, rawThetaPExpressionY, rawThetaPExpressionZ, rawThetaPCExpressionY, rawThetaPCExpressionZ, a_s);
  
  if (theCopy == 0) {
    opserr << "ASDCoupledHinge3D::getCopy -- failed to allocate copy\n";
    exit(-1);
  }
		
  return theCopy;
}

const ID&
ASDCoupledHinge3D::getType ()
{
    (*theCode)(0) = (*matCodes)(0);
    (*theCode)(1) = (*matCodes)(1);
    (*theCode)(2) = (*matCodes)(2);
    (*theCode)(3) = (*matCodes)(3);
    (*theCode)(4) = (*matCodes)(4);
    (*theCode)(5) = (*matCodes)(5);

    return *theCode;
}

int
ASDCoupledHinge3D::getOrder () const
{
    int order = numMats;
    
    return order;
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

    // Idea: get ky and kz and multiply by elastic stiffness in order to see direction? I am not sure because elastic stiffness can change during analysis...?
    double ky = MyMaterial->getStrain();
    double kz = MzMaterial->getStrain();
    if (ky > 0)
        My = par_f1p_Y.getValue() / par_d1p_Y.getValue() * ky;
    else
        My = par_f1n_Y.getValue() / par_d1n_Y.getValue() * ky;
    if (kz > 0)
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
    if (abs(My) < tolM) {
        if (My >= 0)
            My = tolM;
        else
            My = -tolM;
    }
    if (abs(Mz) < tolM) {
        if (Mz >= 0)
            Mz = tolM;
        else
            Mz = -tolM;
    }
    // compute the direction of moments
    double theta = atan2(Mz, My);
    if (theta < 0)
        theta += 2 * DomainData::pi;
    // Now compute the points on the domain with the same direction
    ultimateDomain->getMyMzForNAndDirection(N, theta, My_u_p, Mz_u_p);
#ifdef _DBG_COUPLEDSEC3D
    opserr << "Both moments above tolerance -> biaxial bending\n";
    opserr << "My = " << My << " - Mz = " << Mz << "\n";
    opserr << "angle theta = " << theta << "\n";
    opserr << "getMyMz for N = " << N << " : My = " << My_u_p << " - Mz = " << Mz_u_p << "\n";
#endif
    theta += DomainData::pi;
    if (theta > 2 * DomainData::pi)
        theta -= 2 * DomainData::pi;
    ultimateDomain->getMyMzForNAndDirection(N, theta, My_u_n, Mz_u_n);
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
    if (OPS_EvalDoubleTclStringExpression(initialStiffnessExpressionY.c_str(), EI) < 0) {
        opserr << "\n\nError evaluating expression:\n" << initialStiffnessExpressionY.c_str() << "\n";
        opserr << "Raw string: " << rawInitialStiffnessExpressionY.c_str() << "\n";
        opserr << "N Value = " << N << "\n";
        opserr << "Replaced String: " << initialStiffnessExpressionY.c_str() << "\n";
        opserr << "Replacing again the string: \n";
        replacePlaceholderWithValue(initialStiffnessExpressionY, "__N__", N);
    }
    double th_cap_p;
    if (OPS_EvalDoubleTclStringExpression(thetaPExpressionY.c_str(), th_cap_p) < 0) {
        opserr << "Error evaluating expression:\n" << thetaPExpressionY.c_str() << "\n";
    }
    double th_pcap_p;
    if (OPS_EvalDoubleTclStringExpression(thetaPCExpressionY.c_str(), th_pcap_p) < 0) {
        opserr << "Error evaluating expression:\n" << thetaPCExpressionY.c_str() << "\n";
    }

#ifdef _DBG_COUPLEDSEC3D
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
    opserr << "\nLaw MyMaterial updated:\n";
    MyMaterial->Print(opserr, 2);
#endif

    // update Mz material
    if (OPS_EvalDoubleTclStringExpression(initialStiffnessExpressionZ.c_str(), EI) < 0) {
        opserr << "Error evaluating expression:\n" << initialStiffnessExpressionZ.c_str() << "\n";
    }
    if (OPS_EvalDoubleTclStringExpression(thetaPExpressionZ.c_str(), th_cap_p) < 0) {
        opserr << "Error evaluating expression:\n" << thetaPExpressionZ.c_str() << "\n";
    }
    if (OPS_EvalDoubleTclStringExpression(thetaPCExpressionZ.c_str(), th_pcap_p) < 0) {
        opserr << "Error evaluating expression:\n" << thetaPCExpressionZ.c_str() << "\n";
    }
#ifdef _DBG_COUPLEDSEC3D
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
    opserr << "\nLaw MzMaterial updated:\n";
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

    if (err == 0) {
        updateLaws();

        double My, Mz;
        My = MyMaterial->getStress();
        Mz = MzMaterial->getStress();

        MyMaterial->setTrialStrain(MyMaterial->getStrain());
        MzMaterial->setTrialStrain(MzMaterial->getStrain());

        // compare My with MyMaterial->getStress();
        double errY = abs(MyMaterial->getStress() - My);
        double errZ = abs(MzMaterial->getStress() - Mz);

        errExplicit = std::max(errY, errZ) / MmaxAbs;

        // revert to last commit
        err += MyMaterial->revertToLastCommit();
        err += MzMaterial->revertToLastCommit();
    } 
  
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
    // Massimo
    int res = 0;

    // Need otherDbTag since classTags ID and data ID may be the same size
    if (otherDbTag == 0) 
    otherDbTag = theChannel.getDbTag();
  
 //   // Create ID for tag and section order data
 //   static ID data(4);
 // 
 //   int order = this->getOrder();
 // 
 //   data(0) = this->getTag();
 //   data(1) = otherDbTag;
 //   data(2) = order;
 //   data(3) = numMats;

 //   // Send the tag and section order data
 //   res += theChannel.sendID(this->getDbTag(), cTag, data);
 //   if (res < 0) {
 //       opserr << "ASDCoupledHinge3D::sendSelf -- could not send data ID\n";
	//		    
 //       return res;
 //   }
 // 
 //   // Determine how many classTags there are and allocate ID vector
 //   // for the tags and section code
 //   int numTags = numMats;
 //   ID classTags(2*numTags + numMats);
 //// 
 //   // Loop over the UniaxialMaterials filling in class and db tags
 //   int i, dbTag;
 //   for (i = 0; i < numMats; i++) {
 //       classTags(i) = theAdditions[i]->getClassTag();
 //   
 //       dbTag = theAdditions[i]->getDbTag();
 //   
 //       if (dbTag == 0) {
 //           dbTag = theChannel.getDbTag();
 //           if (dbTag != 0)
 //       theAdditions[i]->setDbTag(dbTag);
 //       }
 //   
 //       classTags(i+numTags) = dbTag;
 //   }
 // 
 //   // Put the Section class and db tags into the ID vector
 //   if (theSection != 0) {
 //   classTags(numTags-1) = theSection->getClassTag();
 //   
 //   dbTag = theSection->getDbTag();
 //   
 //   if (dbTag == 0) {
 //       dbTag = theChannel.getDbTag();
 //       if (dbTag != 0)
 //   theSection->setDbTag(dbTag);
 //   }
 //   
 //   classTags(2*numTags-1) = dbTag;
 //   }
 // 
 //   // Put the UniaxialMaterial codes into the ID vector
 //   int j = 2*numTags;
 //   for (i = 0; i < numMats; i++, j++)
 //   classTags(j) = (*matCodes)(i);
 // 
 //   // Send the material class and db tags and section code
 //   res += theChannel.sendID(otherDbTag, cTag, classTags);
 //   if (res < 0) {
 //   opserr << "ASDCoupledHinge3D::sendSelf -- could not send classTags ID\n";
 //   return res;
 //   }

 //   // Ask the UniaxialMaterials to send themselves
 //   for (i = 0; i < numMats; i++) {
 //   res += theAdditions[i]->sendSelf(cTag, theChannel);
 //   if (res < 0) {
 //       opserr << "ASDCoupledHinge3D::sendSelf -- could not send UniaxialMaterial, i = " << i << endln;
 //       return res;
 //   }
 //   }
 // 
 //   // Ask the Section to send itself
 //   if (theSection != 0) {
 //   res += theSection->sendSelf(cTag, theChannel);
 //   if (res < 0) {
 //       opserr << "ASDCoupledHinge3D::sendSelf -- could not send SectionForceDeformation\n";
 //       return res;
 //   }
 //   }
 // 
    return res;
}


int
ASDCoupledHinge3D::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

 // // Create an ID and receive tag and section order
 // static ID data(5);
 // res += theChannel.recvID(this->getDbTag(), cTag, data);
 // if (res < 0) {
 //   opserr << "ASDCoupledHinge3D::recvSelf -- could not receive data ID\n";
 //   return res;
 // }
 // 
 // this->setTag(data(0));
 // otherDbTag = data(1);
 // int order = data(2);
 // int theSectionOrder = data(3);
 // numMats = data(4);

 // if (order > 0) {
 //   if (e == 0 || e->Size() != order) {
 //     if (e != 0) {
	//delete e;
	//delete s;
	//delete ks;
	//delete fs;
	//delete theCode;
 //     }
 //     e = new Vector(workArea, order);
 //     s = new Vector(&workArea[maxOrder], order);
 //     ks = new Matrix(&workArea[2*maxOrder], order, order);
 //     fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
 //     theCode = new ID(codeArea, order);
 //   }
 // }

 // if (numMats > 0) {
 //   if (matCodes == 0 || matCodes->Size() != numMats) {
 //     if (matCodes != 0)
	//delete matCodes;

 //     matCodes = new ID(numMats);
 //   }
 // }

 // // Determine how many classTags there are and allocate ID vector
 // int numTags = (theSectionOrder == 0) ? numMats : numMats + 1;
 // ID classTags(numTags*2 + numMats);
 // 
 // // Receive the material class and db tags
 // res += theChannel.recvID(otherDbTag, cTag, classTags);
 // if (res < 0) {
 //   opserr << "ASDCoupledHinge3D::recvSelf -- could not receive classTags ID\n";
 //   return res;
 // }

 // // Check if null pointer, allocate if so
 // if (theAdditions == 0) {
 //   theAdditions = new UniaxialMaterial *[numMats];
 //   if (theAdditions == 0) {
 //     opserr << "ASDCoupledHinge3D::recvSelf -- could not allocate UniaxialMaterial array\n";
 //     return -1;
 //   }
 //   // Set pointers to null ... will get allocated by theBroker
 //   for (int j = 0; j < numMats; j++)
 //     theAdditions[j] = 0;
 // }
 // 
 // // Loop over the UniaxialMaterials
 // int i, classTag;
 // for (i = 0; i < numMats; i++) {
 //   classTag = classTags(i);
 //   
 //   // Check if the UniaxialMaterial is null; if so, get a new one
 //   if (theAdditions[i] == 0)
 //     theAdditions[i] = theBroker.getNewUniaxialMaterial(classTag);
 //   
 //   // Check that the UniaxialMaterial is of the right type; if not, delete
 //   // the current one and get a new one of the right type
 //   else if (theAdditions[i]->getClassTag() != classTag) {
 //     delete theAdditions[i];
 //     theAdditions[i] = theBroker.getNewUniaxialMaterial(classTag);
 //   }
 //   
 //   // Check if either allocation failed
 //   if (theAdditions[i] == 0) {
 //     opserr << "ASDCoupledHinge3D::recvSelf -- could not get UniaxialMaterial, i = " << i << endln;
 //     return -1;
 //   }
 //   
 //   // Now, receive the UniaxialMaterial
 //   theAdditions[i]->setDbTag(classTags(i+numTags));
 //   res += theAdditions[i]->recvSelf(cTag, theChannel, theBroker);
 //   if (res < 0) {
 //     opserr << "ASDCoupledHinge3D::recvSelf -- could not receive UniaxialMaterial, i = " << i << endln;
 //     return res;
 //   }
 // }

 // // If there is no Section to receive, return
 // if (theSectionOrder != 0) {
	//
 //
	//classTag = classTags(numTags-1);
 // 
	//// Check if the Section is null; if so, get a new one
	//if (theSection == 0)
	//	theSection = theBroker.getNewSection(classTag);
 // 
	//// Check that the Section is of the right type; if not, delete
	//// the current one and get a new one of the right type
	//else if (theSection->getClassTag() != classTag) {
	//	delete theSection;
	//	theSection = theBroker.getNewSection(classTag);
	//}
 // 
	//// Check if either allocation failed
	//if (theSection == 0) {
	//	opserr << "ASDCoupledHinge3D::recvSelf -- could not get a SectionForceDeformation\n";
	//	return -1;
	//}

	//// Now, receive the Section
	//theSection->setDbTag(classTags(2*numTags-1));
	//res += theSection->recvSelf(cTag, theChannel, theBroker);
	//if (res < 0) {
	//	opserr << "ASDCoupledHinge3D::recvSelf -- could not receive SectionForceDeformation\n";
	//	return res;
	//}
 // }

 // // Fill in the section code
 // int j = 2*numTags;
 // for (i = 0; i < numMats; i++, j++)
 //   (*matCodes)(i) = classTags(j);

  return res;
}
//
//Response*
//ASDCoupledHinge3D::setResponse(const char **argv, int argc, OPS_Stream &output)
//{
//  
//    return SectionForceDeformation::setResponse(responseID, sectInfo);
//
//    return 0;
//}
//
////by SAJalali
//int
//ASDCoupledHinge3D::getResponse(int responseID, Information &sectInfo)
//{
//
//	return SectionForceDeformation::getResponse(responseID, sectInfo);
//
//}


void
ASDCoupledHinge3D::Print(OPS_Stream &s, int flag)
{
    s << "ASDCoupledHinge3D" << endln;
}


int
ASDCoupledHinge3D::getVariable(const char *argv, Information &info)
{

  info.theDouble = 0.0;
  int i;
  int order = numMats;

  const Vector &e = this->getSectionDeformation();
  const ID &code  = this->getType();

  if (strcmp(argv,"axialStrain") == 0) {
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_P)
	info.theDouble  += e(i);
  }  else if (strcmp(argv,"curvatureZ") == 0) {
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_MZ)
	info.theDouble += e(i);
  } else if (strcmp(argv,"curvatureY") == 0) {
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_MY)
	info.theDouble += e(i);
  } else 
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
    // get Domain Values for My and Mz positive and negative (uncoupled) when N = 0
    double tmp;
    //opserr << "\n\nFinding for theta = 0 (should be positive My)\n";
    ultimateDomain->getMyMzForNAndDirection(N, 0.0, My_u_p, tmp);
    //opserr << "\n\n\nFinding for theta = pi (should be negative My)\n";
    ultimateDomain->getMyMzForNAndDirection(N, DomainData::pi, My_u_n, tmp);
    //opserr << "\n\n\nFinding for theta = pi/2 (should be positive Mz)\n";
    ultimateDomain->getMyMzForNAndDirection(N, DomainData::pi / 2, tmp, Mz_u_p);
    //opserr << "\n\n\nFinding for theta = 3pi/2 (should be negative Mz)\n";
    ultimateDomain->getMyMzForNAndDirection(N, DomainData::pi * 3 / 2.0, tmp, Mz_u_n);

    return;
}

Response* ASDCoupledHinge3D::setResponse(const char** argv, int argc, OPS_Stream& output)
{

    Response* theResponse = 0;

    if (argc > 0) {


        if (strcmp(argv[0], "explicitError") == 0) {
            output.tag("SectionOutput");
            output.attr("secType", this->getClassType());
            output.attr("secTag", this->getTag());

            output.tag("ResponseType", "error");
            output.endTag();
            double err = 0.0;
            return theResponse = new MaterialResponse(this, 100, err);

        }

    }

    // If not a fiber response, call the base class method
    return SectionForceDeformation::setResponse(argv, argc, output);
}


int ASDCoupledHinge3D::getResponse(int responseID, Information& sectInfo)
{
    if (responseID == 100) {
        return sectInfo.setDouble(errExplicit);
    }

    return SectionForceDeformation::getResponse(responseID, sectInfo);
}
