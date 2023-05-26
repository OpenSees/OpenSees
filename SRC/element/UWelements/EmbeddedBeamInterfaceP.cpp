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

// Written: Alborz Ghofrani, Diego Turello, Pedro Arduino, U.Washington 
// Created: May 2017
// Description: This file contains the class definition for EmbeddedBeamInterfaceP.

#include <EmbeddedBeamInterfaceP.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CrdTransf.h>
#include <elementAPI.h>
#include <cmath>
#include <NodeIter.h>
#include <FileStream.h>

#ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
#endif

static int num_EmbeddedBeamInterfaceP = 0;
static const double m_Pi = 3.141592653589793;

void *
OPS_EmbeddedBeamInterfaceP(void)
{
    if (num_EmbeddedBeamInterfaceP == 0) {
        num_EmbeddedBeamInterfaceP++;
        opserr << "EmbeddedBeamInterfaceP element - Written: A.Ghofrani, D.Turello, P.Arduino, U.Washington\n";
    }

    Element *theElement = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "Want: EmbeddedBeamInterfaceP tag? \n";
        return 0;
    }

    int iData[1];
    int eleTag = 0;
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element EmbeddedBeamInterfaceP" << endln;
        return 0;
    }

    eleTag = iData[0];

    theElement = new EmbeddedBeamInterfaceP(eleTag);

    if (theElement == 0) {
        opserr << "WARNING could not create element of type EmbeddedBeamInterfaceP\n";
        return 0;
    }

    return theElement;
}


EmbeddedBeamInterfaceP::EmbeddedBeamInterfaceP(int tag) :
  Element(tag, ELE_TAG_EmbeddedBeamInterfaceP),
  theSolidTags(0), solidNodeTags(0), theBeamTags(0), beamNodeTags(0), theNodes(0),
  crdTransf(0)  
{

}

EmbeddedBeamInterfaceP::EmbeddedBeamInterfaceP(int tag, std::vector <int> beamTag, std::vector <int> solidTag, int crdTransfTag,
    std::vector <double>  beamRho, std::vector <double>  beamTheta, std::vector <double>  solidXi, std::vector <double>  solidEta,
    std::vector <double>  solidZeta, double radius, std::vector <double> area, std::vector <double> length, double penaltyParam, 
    bool writeConnectivity, const char * connectivityFN):
  Element(tag, ELE_TAG_EmbeddedBeamInterfaceP), m_beam_radius(radius),
  theSolidTags(0), solidNodeTags(0), theBeamTags(0), beamNodeTags(0), theNodes(0),
  m_ep(penaltyParam), mQa(3, 3), mQb(3, 3), mQc(3, 3),
  mBphi(3, 12), mBu(3, 12), mHf(3, 12), m_Ns(8), crdTransf(0)
{
    // get domain to access element tags and their nodes
    Domain &theDomain = *(OPS_GetDomain());

    m_numEmbeddedPoints = solidTag.size();
    theSolidTags = new int[m_numEmbeddedPoints];
    solidNodeTags = new int[8 * m_numEmbeddedPoints];
    theBeamTags = new int[m_numEmbeddedPoints];
    beamNodeTags = new int[2 * m_numEmbeddedPoints];
    m_beam_rho = m_beam_theta = m_solid_xi = m_solid_eta = m_solid_zeta = m_area = m_beamLength = Vector(m_numEmbeddedPoints);

    std::set <int> uniqueSolidNodeTags;
    std::set <int> uniqueBeamNodeTags;
    std::set <int> uniqueBeamTags;
    Element *theElement;
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        theSolidTags[ii] = solidTag[ii];
        theBeamTags[ii] = beamTag[ii];
        m_solid_xi(ii) = solidXi[ii];
        m_solid_eta(ii) = solidEta[ii];
        m_solid_zeta(ii) = solidZeta[ii];
        m_beam_rho(ii) = beamRho[ii];
        m_beam_theta(ii) = beamTheta[ii];
        m_area(ii) = area[ii];
        m_beamLength(ii) = length[ii];

        theElement = theDomain.getElement(solidTag[ii]);
        if (ii == 0)
            m_numSolidDOF = theElement->getNodePtrs()[0]->getNumberDOF();
        // opserr << "Point " << ii +1 << " : element " << solidTag[ii] << " at (" << solidXi[ii] << "," << solidEta[ii] << "," << solidZeta[ii] << ") , beam: " << beamTag << " at (" << beamRho[ii] << "," << beamTheta[ii] << ")" << endln;
        for (int jj = 0; jj < 8; jj++)
        {
            uniqueSolidNodeTags.insert(theElement->getNodePtrs()[jj]->getTag());
            solidNodeTags[ii * 8 + jj] = theElement->getNodePtrs()[jj]->getTag();
        }
        uniqueBeamTags.insert(beamTag[ii]);
        theElement = theDomain.getElement(beamTag[ii]);
        // opserr << "Point " << ii +1 << " : element " << solidTag[ii] << " at (" << solidXi[ii] << "," << solidEta[ii] << "," << solidZeta[ii] << ") , beam: " << beamTag << " at (" << beamRho[ii] << "," << beamTheta[ii] << ")" << endln;
        for (int jj = 0; jj < 2; jj++)
        {
            uniqueBeamNodeTags.insert(theElement->getNodePtrs()[jj]->getTag());
            beamNodeTags[ii * 2 + jj] = theElement->getNodePtrs()[jj]->getTag();
        }
    }


    m_numSolidNodes = (int)uniqueSolidNodeTags.size();
    m_numBeamNodes  = (int)uniqueBeamNodeTags.size();
    EBIP_numNodes   = m_numSolidNodes + m_numBeamNodes;
    EBIP_numDOF     = m_numSolidNodes * 3 + m_numBeamNodes * 6;

    m_solidInitDisp = Vector(m_numSolidNodes * 3);
    m_beamInitDisp  = Vector(m_numBeamNodes * 6);

    // opserr << m_numSolidNodes << endln;

    externalNodes = ID(EBIP_numNodes);
    theNodes = new Node*[EBIP_numNodes];

    int count = 0;
    for (std::set <int>::iterator it = uniqueSolidNodeTags.begin(); it != uniqueSolidNodeTags.end(); ++it)
    {
        m_solidNodeMap[*it] = count;
        externalNodes(count) = *it;

        theNodes[count] = theDomain.getNode(*it);

        Vector tempDisp = theNodes[count]->getDisp();
        m_solidInitDisp(count * 3 + 0) = tempDisp(0);
        m_solidInitDisp(count * 3 + 1) = tempDisp(1);
        m_solidInitDisp(count * 3 + 2) = tempDisp(2);

        count++;
    }

    int curCount = count;
    for (std::set <int>::iterator it = uniqueBeamNodeTags.begin(); it != uniqueBeamNodeTags.end(); ++it)
    {
        m_beamNodeMap[*it] = count - curCount;
        externalNodes(count) = *it;

        theNodes[count] = theDomain.getNode(*it);

        Vector tempDisp = theNodes[count]->getDisp();
        m_beamInitDisp((count - curCount) * 6 + 0) = tempDisp(0);
        m_beamInitDisp((count - curCount) * 6 + 1) = tempDisp(1);
        m_beamInitDisp((count - curCount) * 6 + 2) = tempDisp(2);
        m_beamInitDisp((count - curCount) * 6 + 3) = tempDisp(3);
        m_beamInitDisp((count - curCount) * 6 + 4) = tempDisp(4);
        m_beamInitDisp((count - curCount) * 6 + 5) = tempDisp(5);

        count++;
    }

    m_Lambda             = Vector(6 * m_numBeamNodes);
    m_InterfaceForces    = Vector(EBIP_numDOF);
    m_InterfaceStiffness = Matrix(EBIP_numDOF, EBIP_numDOF);

    mA   = Matrix(3 * m_numSolidNodes, 6 * m_numBeamNodes);
    mB   = Matrix(6 * m_numBeamNodes, 6 * m_numBeamNodes);
    // mAt  = Matrix(6 * m_numBeamNodes, 3 * m_numSolidNodes);
    // mBt  = Matrix(6 * m_numBeamNodes, 6 * m_numBeamNodes);
    // mAAt = Matrix(3 * m_numSolidNodes, 3 * m_numSolidNodes);
    // mBBt = Matrix(6 * m_numBeamNodes, 6 * m_numBeamNodes);
    // mABt = Matrix(3 * m_numSolidNodes, 6 * m_numBeamNodes);
    


    // get the coordinate transformation object
    crdTransf = OPS_getCrdTransf(crdTransfTag)->getCopy3d();


    if (writeConnectivity)
    {
        FileStream connFile(connectivityFN, APPEND);
        connFile << m_numSolidNodes << " " << (int)uniqueBeamTags.size();
        for (int ii = 0; ii < m_numSolidNodes; ii++)
            connFile << " " << externalNodes(ii);
        int curBeamTag = theBeamTags[0];
        connFile << " " << beamNodeTags[0] << " " << beamNodeTags[1];
        for (int ii = 1; ii < m_numEmbeddedPoints; ii++)
            if (theBeamTags[ii] == curBeamTag)
                continue;
            else
            {
                curBeamTag = theBeamTags[ii];
                connFile << " " << beamNodeTags[2 * ii] << " " << beamNodeTags[2 * ii + 1];
            }
        connFile << endln;
        connFile.close();
    }

}

EmbeddedBeamInterfaceP::EmbeddedBeamInterfaceP()
  : Element(0, ELE_TAG_EmbeddedBeamInterfaceP),
  theSolidTags(0), solidNodeTags(0), theBeamTags(0), beamNodeTags(0), theNodes(0),
    crdTransf(0)
{
  if (theSolidTags != 0)
    delete [] theSolidTags;
  if (solidNodeTags != 0)
    delete [] solidNodeTags;
  if (theBeamTags != 0)
    delete [] theBeamTags;
  if (beamNodeTags != 0)
    delete [] beamNodeTags;

  if (crdTransf != 0)
    delete crdTransf;
  if (theNodes != 0)
    delete [] theNodes;
}

EmbeddedBeamInterfaceP::~EmbeddedBeamInterfaceP()
{

}

int
EmbeddedBeamInterfaceP::getNumExternalNodes(void) const
{
    return EBIP_numNodes;
}

const ID&
EmbeddedBeamInterfaceP::getExternalNodes(void)
{
    return externalNodes;
}

Node **
EmbeddedBeamInterfaceP::getNodePtrs(void)
{
    return theNodes;
}

int
EmbeddedBeamInterfaceP::getNumDOF(void)
{
    return EBIP_numDOF;
}

int
EmbeddedBeamInterfaceP::revertToLastCommit(void)
{
    return 0;
}

int
EmbeddedBeamInterfaceP::revertToStart(void)
{
    return 0;
}


const Matrix&
EmbeddedBeamInterfaceP::getTangentStiff(void)
{
    return m_InterfaceStiffness;
}

const Matrix&
EmbeddedBeamInterfaceP::getInitialStiff(void)
{
    return this->getTangentStiff();
}

const Vector&
EmbeddedBeamInterfaceP::getResistingForce(void)
{
    m_InterfaceForces.Zero();
    Vector temp2(6 * m_numBeamNodes), temp(3 * m_numSolidNodes);

    temp  = mA * m_Lambda;
    temp2 = -1.0 * mB * m_Lambda;

    int II;
    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
    {
        II = m_numSolidDOF * (ii / 3) + (ii % 3);
        m_InterfaceForces(II) = temp(ii);
    }

    for (int ii = 0; ii < 6 * m_numBeamNodes; ii++)
        m_InterfaceForces(m_numSolidDOF * m_numSolidNodes + ii) = temp2(ii);

    return m_InterfaceForces;
}

int
EmbeddedBeamInterfaceP::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
EmbeddedBeamInterfaceP::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
    &theBroker)
{
    return 0;
}

int
EmbeddedBeamInterfaceP::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    return 0;
}

void
EmbeddedBeamInterfaceP::Print(OPS_Stream &s, int flag)
{
    return;
}

Response*
EmbeddedBeamInterfaceP::setResponse(const char **argv, int argc,
    OPS_Stream &s)
{
    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "globalForce") == 0)
    {
        return new ElementResponse(this, 1, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "displacement") == 0 || strcmp(argv[0], "disp") == 0)
    {
        return new ElementResponse(this, 2, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "beamForce") == 0 || strcmp(argv[0], "beamInteractionForce") == 0)
    {
        return new ElementResponse(this, 3, Vector(12 * (m_numBeamNodes - 1)));
    }
    else if (strcmp(argv[0], "solidForce") == 0 || strcmp(argv[0], "solidInteractionForce") == 0)
    {
        return new ElementResponse(this, 4, Vector(3 * m_numSolidNodes));
    }
    else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "locForce") == 0)
    {
        return new ElementResponse(this, 5, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "axialForce") == 0 || strcmp(argv[0], "locForceAxial") == 0)
    {
        return new ElementResponse(this, 6, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "radialForce") == 0 || strcmp(argv[0], "locForceNormal") == 0)
    {
        return new ElementResponse(this, 7, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "tangentialForce") == 0 || strcmp(argv[0], "locForceTangent") == 0)
    {
        return new ElementResponse(this, 8, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "localDisplacement") == 0 || strcmp(argv[0], "locDisp") == 0)
    {
        return new ElementResponse(this, 9, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "axialDisp") == 0 || strcmp(argv[0], "locDispAxial") == 0)
    {
        return new ElementResponse(this, 10, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "radialDisp") == 0 || strcmp(argv[0], "locDispNormal") == 0)
    {
        return new ElementResponse(this, 11, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "tangentialDisp") == 0 || strcmp(argv[0], "locDispTangent") == 0)
    {
        return new ElementResponse(this, 12, Vector(3 * m_numEmbeddedPoints));
    }
    else if (strcmp(argv[0], "beamLocalForce") == 0 || strcmp(argv[0], "beamInteractionLocalForce") == 0)
    {
        return new ElementResponse(this, 13, Vector(12 * (m_numBeamNodes - 1)));
    }
    else {
        opserr << "EmbeddedBeamInterfaceP Recorder, " << argv[0] << "is an unknown recorder request"
            << "  Element tag : " << this->getTag() << endln;
        return 0;
    }
}

int
EmbeddedBeamInterfaceP::getResponse(int responseID, Information &eleInformation)
{
    if (responseID == 1) // force
        return eleInformation.setVector(GetInteractionPtForce(0));
    else if (responseID == 2) // displacement
        return eleInformation.setVector(GetInteractionPtDisp(0));
    else if (responseID == 3) // beamForce
        return eleInformation.setVector(GetElementalForce(1));
    else if (responseID == 4) // solidForce
        return eleInformation.setVector(GetElementalForce(2));
    else if (responseID == 5) // localForce
        return eleInformation.setVector(GetInteractionPtForce(1));
    else if (responseID == 6) // axialForce
        return eleInformation.setVector(GetInteractionPtForce(2));
    else if (responseID == 7) // radialForce
        return eleInformation.setVector(GetInteractionPtForce(3));
    else if (responseID == 8) // tangentialForce
        return eleInformation.setVector(GetInteractionPtForce(4));
    else if (responseID == 9) // localDisplacement
        return eleInformation.setVector(GetInteractionPtDisp(1));
    else if (responseID == 10) // axialDisp
        return eleInformation.setVector(GetInteractionPtDisp(2));
    else if (responseID == 11) // radialDisp
        return eleInformation.setVector(GetInteractionPtDisp(3));
    else if (responseID == 12) // tangentialDisp
        return eleInformation.setVector(GetInteractionPtDisp(4));
    else if (responseID == 13) // beamLocalForce
        return eleInformation.setVector(GetElementalForce(3));
    else {
        opserr << "EmbeddedBeamInterfaceP, tag = " << this->getTag()
            << " -- unknown request" << endln;
        return -1;
    }
}

int
EmbeddedBeamInterfaceP::setParameter(const char **argv, int argc, Parameter &param)
{
    return 0;
}

int
EmbeddedBeamInterfaceP::updateParameter(int parameterID, Information &info)
{
    return 0;
}

void
EmbeddedBeamInterfaceP::setDomain(Domain *theDomain)
{
    for (int ii = 0; ii < m_numSolidNodes + m_numBeamNodes; ii++)
    {
        if (theNodes[ii] == 0)
        {
            opserr << "Could not find node " << externalNodes(ii) << "." << endln;
            return;
        }
        if (!((theNodes[ii]->getNumberDOF() == 3) || (theNodes[ii]->getNumberDOF() == 4)) && (ii < m_numSolidNodes))
        {
            opserr << "Solid node " << externalNodes(ii) << " has to have 3 or 4 degrees of freedom." << endln;
            return;
        }
        if ((theNodes[ii]->getNumberDOF() != 6) && (ii > m_numSolidNodes - 1))
        {
            opserr << "Beam node " << externalNodes(ii) << " has to have 6 degrees of freedom." << endln;
            return;
        }
    }
    
        


    // initialize the transformation
    if (crdTransf->initialize(theDomain->getNode(beamNodeTags[0]), theDomain->getNode(beamNodeTags[1])))
    {
        opserr << "EmbeddedBeamInterfaceP::setDomain(): Error initializing coordinate transformation";
        return;
    }

    m_beam_length = crdTransf->getInitialLength();
    if (m_beam_length < 1.0e-12)
    {
        opserr << "FATAL ERROR EmbeddedBeamInterfaceP (tag: " << this->getTag() << ") : "
            << "Beam element has zero length." << endln;
        return;
    }

    Vector initXAxis(3);
    Vector initYAxis(3);
    Vector initZAxis(3);
    crdTransf->getLocalAxes(initXAxis, initYAxis, initZAxis);
    // fill mQa
    for (int i = 0; i < 3; i++)
    {
        mQa(i, 0) = initXAxis(i);
        mQa(i, 1) = initYAxis(i);
        mQa(i, 2) = initZAxis(i);
    }
    // set mQb = mQa : beam column element requires zero initial twist
    // if mQa = mQb -> mchi = 0
    mQc = mQb = mQa;


    // calculate A and B
    mA.Zero();
    mB.Zero();
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));
        ComputeHf(mHf, m_beam_theta(ii));
        ComputeBphiAndBu(mBphi, mBu, m_beamLength(ii));

        Vector c2(3), c3(3);
        for (int jj = 0; jj < 3; jj++)
        {
            c2(jj) = mQc(jj, 1);
            c3(jj) = mQc(jj, 2);
        }
        Matrix Hb(3, 12);
        Hb = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;



        //Matrix mM(3, 3);
        //mM(0, 1) = -m_beam_radius * sin(m_beam_theta(ii));
        //mM(0, 2) = m_beam_radius * cos(m_beam_theta(ii));
        //mM(1, 0) = m_beam_radius * sin(m_beam_theta(ii));
        //mM(2, 0) = -m_beam_radius * cos(m_beam_theta(ii));
        //Matrix Hb(3, 12);
        ////Hb  = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
        //Hb = mBu - mQc * mM * mBphi;

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * ii + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * ii + 1]];

        for (int jj = 0; jj < 8; jj++)
        {
            int solidNodeInA = m_solidNodeMap[solidNodeTags[8 * ii + jj]];

            double curNs = m_Ns(jj);

            for (int kk = 0; kk < 6; kk++)
            {
                mA(3 * solidNodeInA + 0, 6 * beamNodeInA1 + kk) += curNs * mHf(0, kk) * m_area(ii);
                mA(3 * solidNodeInA + 1, 6 * beamNodeInA1 + kk) += curNs * mHf(1, kk) * m_area(ii);
                mA(3 * solidNodeInA + 2, 6 * beamNodeInA1 + kk) += curNs * mHf(2, kk) * m_area(ii);

                mA(3 * solidNodeInA + 0, 6 * beamNodeInA2 + kk) += curNs * mHf(0, kk + 6) * m_area(ii);
                mA(3 * solidNodeInA + 1, 6 * beamNodeInA2 + kk) += curNs * mHf(1, kk + 6) * m_area(ii);
                mA(3 * solidNodeInA + 2, 6 * beamNodeInA2 + kk) += curNs * mHf(2, kk + 6) * m_area(ii);
            }
        }

        Matrix HbT = Transpose(3, 12, Hb);
        Matrix temp = HbT * mHf * m_area(ii);

        for (int jj = 0; jj < 6; jj++)
        {
            for (int kk = 0; kk < 6; kk++)
            {
                mB(6 * beamNodeInA1 + jj, 6 * beamNodeInA1 + kk) += temp(jj, kk);
                mB(6 * beamNodeInA1 + jj, 6 * beamNodeInA2 + kk) += temp(jj, kk + 6);
                mB(6 * beamNodeInA2 + jj, 6 * beamNodeInA1 + kk) += temp(jj + 6, kk);
                mB(6 * beamNodeInA2 + jj, 6 * beamNodeInA2 + kk) += temp(jj + 6, kk + 6);
            }
        }

        /*FileStream tempFile("tempA.dat",APPEND);
        tempFile.setPrecision(15);
        tempFile << mA;
        FileStream tempFile2("tempB.dat",APPEND);
        tempFile2.setPrecision(15);
        tempFile2  << mB;*/
    }
    Matrix mAt (6 * m_numBeamNodes, 3 * m_numSolidNodes);
    Matrix mBt (6 * m_numBeamNodes, 6 * m_numBeamNodes);
    Matrix mAAt(3 * m_numSolidNodes, 3 * m_numSolidNodes);
    Matrix mBBt(6 * m_numBeamNodes, 6 * m_numBeamNodes);
    Matrix mABt(3 * m_numSolidNodes, 6 * m_numBeamNodes);

    mAt = Transpose(3 * m_numSolidNodes, 6 * m_numBeamNodes, mA);
    mBt = Transpose(6 * m_numBeamNodes, 6 * m_numBeamNodes, mB);
    mAAt = mA * mAt;
    mBBt = mB * mBt;
    mABt = mA * mBt;

    m_InterfaceStiffness.Zero();

    int II, JJ;
    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
        for (int jj = 0; jj < 3 * m_numSolidNodes; jj++)
        {
            II = m_numSolidDOF * (ii / 3) + (ii % 3);
            JJ = m_numSolidDOF * (jj / 3) + (jj % 3);
            m_InterfaceStiffness(II, JJ) = mAAt(ii, jj);
        }

    for (int ii = 0; ii < 3 * m_numSolidNodes; ii++)
        for (int jj = 0; jj < 6 * m_numBeamNodes; jj++)
        {
            II = m_numSolidDOF * (ii / 3) + (ii % 3);
            m_InterfaceStiffness(II, m_numSolidDOF * m_numSolidNodes + jj) = -1.0 * mABt(ii, jj);
            m_InterfaceStiffness(m_numSolidDOF * m_numSolidNodes + jj, II) = -1.0 * mABt(ii, jj);
        }

    for (int ii = 0; ii < 6 * m_numBeamNodes; ii++)
        for (int jj = 0; jj < 6 * m_numBeamNodes; jj++)
            m_InterfaceStiffness(m_numSolidDOF * m_numSolidNodes + ii, m_numSolidDOF * m_numSolidNodes + jj) = mBBt(ii, jj);

    m_InterfaceStiffness *= m_ep;

    this->DomainComponent::setDomain(theDomain);
    return;
}

int
EmbeddedBeamInterfaceP::update(void)
{
    Vector sDisp(3 * m_numSolidNodes), bDisp(6 * m_numBeamNodes);
    for (int ii = 0; ii < m_numSolidNodes; ii++)
    {
        Vector sDispCur   = theNodes[ii]->getTrialDisp();
        sDisp(3 * ii    ) = sDispCur(0) - m_solidInitDisp(ii * 3 + 0);
        sDisp(3 * ii + 1) = sDispCur(1) - m_solidInitDisp(ii * 3 + 1);
        sDisp(3 * ii + 2) = sDispCur(2) - m_solidInitDisp(ii * 3 + 2);
    }

    for (int ii = 0; ii < m_numBeamNodes; ii++)
    {
        Vector bDispCur = theNodes[m_numSolidNodes + ii]->getTrialDisp();
        for (int jj = 0; jj < 6; jj++)
            bDisp(6 * ii + jj) = bDispCur(jj);
        //  bDisp(6 * ii + jj) = bDispCur(jj) - m_beamInitDisp(ii * 6 + jj);
    }

    m_Lambda.Zero();
    m_Lambda.addMatrixTransposeVector(1.0, mA, sDisp, 1.0);
    m_Lambda.addMatrixTransposeVector(1.0, mB, bDisp, -1.0);
    m_Lambda *= m_ep;

    return 0;
}

int
EmbeddedBeamInterfaceP::commitState(void)
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "EmbeddedBeamInterfaceP::commitState() - failed in base class";
    }

    return retVal;
}

int EmbeddedBeamInterfaceP::updateShapeFuncs(double xi, double eta, double zeta, double rho, double L)
{
    if ((xi < -1.0) || (xi > 1.0) || (eta < -1.0) || (eta > 1.0) || (zeta < -1.0) || (zeta > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    if ((rho < -1.0) || (rho > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    double rho2 = rho * rho;
    double rho3 = rho * rho2;

    m_Ns(0) = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
    m_Ns(1) = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
    m_Ns(2) = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
    m_Ns(3) = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
    m_Ns(4) = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
    m_Ns(5) = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
    m_Ns(6) = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
    m_Ns(7) = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);

    m_Hb1 = 0.125 * (4.0 - 6.0 * rho + 2.0 * rho3);
    m_Hb3 = 0.125 * (4.0 + 6.0 * rho - 2.0 * rho3);
    m_Hb2 = 0.125 * L * (1.0 - rho - rho2 + rho3);
    m_Hb4 = 0.125 * L * (-1.0 - rho + rho2 + rho3);

    m_Nb1 = 0.5 * (1 - rho);
    m_Nb2 = 0.5 * (1 + rho);

    m_dH1 = 1.5 * (-1.0 + rho2);
    m_dH3 = 1.5 * (1.0 - rho2);
    m_dH2 = 0.25 * L * (-1.0 - 2.0 * rho + 3.0 * rho2);
    m_dH4 = 0.25 * L * (-1.0 + 2.0 * rho + 3.0 * rho2);

    return 0;
}

Matrix
EmbeddedBeamInterfaceP::Transpose(int dim1, int dim2, const Matrix &M)
{
    // copied from transpose function in Brick.cpp

    Matrix Mtran(dim2, dim1);

    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++)
            Mtran(j, i) = M(i, j);

    return Mtran;
}


Matrix
EmbeddedBeamInterfaceP::ComputeSkew(Vector th)
{
    Matrix skew_th(3, 3);

    skew_th(0, 0) = 0.0;
    skew_th(0, 1) = -th(2);
    skew_th(0, 2) = th(1);
    skew_th(1, 0) = th(2);
    skew_th(1, 1) = 0.0;
    skew_th(1, 2) = -th(0);
    skew_th(2, 0) = -th(1);
    skew_th(2, 1) = th(0);
    skew_th(2, 2) = 0.0;

    return skew_th;
}

void
EmbeddedBeamInterfaceP::ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu, double L)
{
    Matrix dummy1(3, 3);
    Matrix dummy2(3, 3);
    Matrix dummy3(3, 3);
    Matrix dummy4(3, 3);

    Bphi.Zero();
    Bu.Zero();

    // Compute Bphi(0:2, 3:5)
    dummy1.Zero();
    dummy2.Zero();
    dummy3.Zero();
    dummy4.Zero();

    // dummy1 = N1 * Qc*(E1 dyadic E1)
    dummy1(0, 0) = m_Nb1*mQc(0, 0);
    dummy1(1, 0) = m_Nb1*mQc(1, 0);
    dummy1(2, 0) = m_Nb1*mQc(2, 0);
    // dummy1 += dH2 * Qc*P1
    dummy1(0, 1) = m_dH2*mQc(0, 1) / L;  // dH2 * mQc(0:2,1:2)
    dummy1(1, 1) = m_dH2*mQc(1, 1) / L;
    dummy1(2, 1) = m_dH2*mQc(2, 1) / L;
    dummy1(0, 2) = m_dH2*mQc(0, 2) / L;
    dummy1(1, 2) = m_dH2*mQc(1, 2) / L;
    dummy1(2, 2) = m_dH2*mQc(2, 2) / L;
    // dummy2 = Qa^T
    dummy2 = Transpose(3, 3, mQa);
    // dummy3 = Qc * (N1*(E1 dyadic E1)+ dH2 * P1) * Qa^T
    dummy3 = dummy1*dummy2;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bphi(i, 3 + j) = dummy3(i, j);

    // Reuse parts of dummy1 and dummy2 to calculate Bu(0:2,0:2)
    // dummy1 += H1 * Qc*P1
    dummy1(0, 1) = m_Hb1*mQc(0, 1);  // H1 * mQc(0:2,1:2)
    dummy1(1, 1) = m_Hb1*mQc(1, 1);
    dummy1(2, 1) = m_Hb1*mQc(2, 1);
    dummy1(0, 2) = m_Hb1*mQc(0, 2);
    dummy1(1, 2) = m_Hb1*mQc(1, 2);
    dummy1(2, 2) = m_Hb1*mQc(2, 2);
    // dummy3 = Qc * (N1*(E1 dyadic E1)+ H1 * P1) * Qa^T
    dummy3 = dummy1*dummy2;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bu(i, j) = dummy3(i, j);

    // Reuse dummy2 and Compute Bphi(0:2, 0:2) and Bu(0:2, 3:5)
    dummy1.Zero();
    dummy3.Zero();
    // dummy1 = Qc*E^R*P1 (E^R is the skew symmetric meatrix for E1 cross product => E1 x a = [E^R].a)
    dummy1(0, 1) = mQc(0, 2);
    dummy1(0, 2) = -mQc(0, 1);
    dummy1(1, 1) = mQc(1, 2);
    dummy1(1, 2) = -mQc(1, 1);
    dummy1(2, 1) = mQc(2, 2);
    dummy1(2, 2) = -mQc(2, 1);
    // dummy3 = Qc*E^R*P1*Qa^T
    dummy3 = dummy1*dummy2;
    // Compute Bphi(0:2, 0:2)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bphi(i, j) = m_dH1 / L * dummy3(i, j);
    // Compute Bu(0:2, 3:5)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bu(i, 3 + j) = -m_Hb2 * dummy3(i, j);

    // Reuse dummy1 and Compute Bphi(0:2, 6:8) and Bu(0:2, 9:11)
    dummy2.Zero();
    dummy3.Zero();

    // dummy2 = Qb^T
    dummy2 = Transpose(3, 3, mQb);
    // dummy3 = Qc*E^R*P1*Qb^T
    dummy3 = dummy1*dummy2;
    // Compute Bphi(0:2, 6 : 8)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bphi(i, 6 + j) = m_dH3 / L * dummy3(i, j);
    // Compute Bu(0:2, 9:11)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bu(i, 9 + j) = -m_Hb4 * dummy3(i, j);

    // Reuse dummy2 and Compute Bphi(0:2, 9:11)
    dummy1.Zero();
    dummy3.Zero();
    // dummy1 = N2 * Qc*(E1 dyadic E1)
    dummy1(0, 0) = m_Nb2*mQc(0, 0);  // N2 * mQc(0:2,0)
    dummy1(1, 0) = m_Nb2*mQc(1, 0);
    dummy1(2, 0) = m_Nb2*mQc(2, 0);
    // dummy1 += dH4 * Qc*P1
    dummy1(0, 1) = m_dH4*mQc(0, 1) / L;     // dH4 * mQc(0:2,1:2)
    dummy1(1, 1) = m_dH4*mQc(1, 1) / L;
    dummy1(2, 1) = m_dH4*mQc(2, 1) / L;
    dummy1(0, 2) = m_dH4*mQc(0, 2) / L;
    dummy1(1, 2) = m_dH4*mQc(1, 2) / L;
    dummy1(2, 2) = m_dH4*mQc(2, 2) / L;
    // dummy3 = Qc * (N2*(E1 dyadic E1)+ dH4 * P1) * Qb^T
    dummy3 = dummy1*dummy2;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bphi(i, 9 + j) = dummy3(i, j);

    // Reuse parts of dummy1 and dummy2 to calculate Bu(0:2,6:8)
    // dummy1 += dH4 * Qc*P1
    dummy1(0, 1) = m_Hb3*mQc(0, 1);     // H3 * mQc(0:2,1:2)
    dummy1(1, 1) = m_Hb3*mQc(1, 1);
    dummy1(2, 1) = m_Hb3*mQc(2, 1);
    dummy1(0, 2) = m_Hb3*mQc(0, 2);
    dummy1(1, 2) = m_Hb3*mQc(1, 2);
    dummy1(2, 2) = m_Hb3*mQc(2, 2);
    // dummy3 = Qc * (N2*(E1 dyadic E1)+ H3 * P1) * Qb^T
    dummy3 = dummy1*dummy2;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Bu(i, 6 + j) = dummy3(i, j);

    //// new way to calculate Bu
    //Bu.Zero();
    //dummy1.Zero();
    //dummy2.Zero();
    //dummy3.Zero();
    //dummy4.Zero();

    //dummy2 = Transpose(3, 3, mQa);
    //dummy3 = Transpose(3, 3, mQb);

    //dummy1(0, 0) = m_Nb1;
    //dummy1(1, 1) = m_Hb1;
    //dummy1(2, 2) = m_Hb1;

    //dummy4 = mQc * dummy1 * dummy2;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bu(i, j) = dummy4(i, j);

    //dummy1.Zero();
    //dummy1(1, 2) = m_Hb2;
    //dummy1(2, 1) = -m_Hb2;

    //dummy4 = mQc * dummy1 * dummy2;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bu(i, j + 3) = dummy4(i, j);


    //dummy1.Zero();
    //dummy1(0, 0) = m_Nb2;
    //dummy1(1, 1) = m_Hb3;
    //dummy1(2, 2) = m_Hb3;

    //dummy4 = mQc * dummy1 * dummy3;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bu(i, j + 6) = dummy4(i, j);

    //dummy1.Zero();
    //dummy1(1, 2) = m_Hb4;
    //dummy1(2, 1) = -m_Hb4;

    //dummy4 = mQc * dummy1 * dummy3;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bu(i, j + 9) = dummy4(i, j);

    //// new way to calculate Bphi
    //Bphi.Zero();
    //dummy1.Zero();
    //dummy2.Zero();
    //dummy3.Zero();
    //dummy4.Zero();

    //dummy2 = Transpose(3, 3, mQa);
    //dummy3 = Transpose(3, 3, mQb);

    //dummy1(1, 2) = -m_dH1 / L;
    //dummy1(2, 1) =  m_dH1 / L;

    //dummy4 = dummy1 * dummy2;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bphi(i, j) = dummy4(i, j);

    //dummy1.Zero();
    //dummy1(0, 0) = m_Nb1;
    //dummy1(1, 1) = m_dH2 / L;
    //dummy1(2, 2) = m_dH2 / L;

    //dummy4 = dummy1 * dummy2;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bphi(i, j + 3) = dummy4(i, j);


    //dummy1.Zero();
    //dummy1(1, 2) = -m_dH3 / L;
    //dummy1(2, 1) =  m_dH3 / L;

    //dummy4 = dummy1 * dummy3;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bphi(i, j + 6) = dummy4(i, j);

    //dummy1.Zero();
    //dummy1(0, 0) = m_Nb2;
    //dummy1(1, 1) = m_dH4 / L;
    //dummy1(2, 2) = m_dH4 / L;

    //dummy4 = dummy1 * dummy3;

    //for (int i = 0; i < 3; i++)
    //    for (int j = 0; j < 3; j++)
    //        Bphi(i, j + 9) = dummy4(i, j);

    return;
}

void EmbeddedBeamInterfaceP::ComputeHf(Matrix & Hf, double theta)
{
    Hf.Zero();

    double oneOver2PiR = 0.5 / m_Pi / m_beam_radius;
    double oneOver2PiR2 = oneOver2PiR / m_beam_radius;
    for (int ii = 0; ii < 3; ii++)
    {
        for (int jj = 0; jj < 3; jj++)
        {
            Hf(ii, jj) = oneOver2PiR * m_Nb1 * mQa(jj, ii);
            Hf(ii, jj + 6) = oneOver2PiR * m_Nb2 * mQb(jj, ii);
        }
        Hf(0, ii + 3) = 2.0 * oneOver2PiR2 * m_Nb1 * (mQa(ii, 1) * sin(theta) - mQa(ii, 2) * cos(theta));
        Hf(1, ii + 3) = -oneOver2PiR2 * mQa(ii, 0) * m_Nb1 * sin(theta);
        Hf(2, ii + 3) = oneOver2PiR2 * mQa(ii, 0) * m_Nb1 * cos(theta);
        Hf(0, ii + 9) = 2.0 * oneOver2PiR2 * m_Nb2 * (mQb(ii, 1) * sin(theta) - mQb(ii, 2) * cos(theta));
        Hf(1, ii + 9) = -oneOver2PiR2 * mQb(ii, 0) * m_Nb2 * sin(theta);
        Hf(2, ii + 9) = oneOver2PiR2 * mQb(ii, 0) * m_Nb2 * cos(theta);
    }

    Hf = mQc * Hf;
    return;
}

Vector EmbeddedBeamInterfaceP::GetInteractionPtDisp(int flag)
{
    Vector res(3 * m_numEmbeddedPoints);
    Vector c1(3), c2(3), c3(3);
    Vector bDisp(6 * m_numBeamNodes);
    Matrix Hb(3, 12);
    Vector ptDisp(3);

    // update local coordinate system
    for (int ii = 0; ii < 3; ii++)
    {
        c1(ii) = mQc(ii, 0);
        c2(ii) = mQc(ii, 1);
        c3(ii) = mQc(ii, 2);
    }

    for (int ii = 0; ii < m_numBeamNodes; ii++)
    {
        Vector bDispCur = theNodes[m_numSolidNodes + ii]->getTrialDisp();
        for (int jj = 0; jj < 6; jj++)
            bDisp(6 * ii + jj) = bDispCur(jj);
    }


    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));

        // update the matrices that define kinematics of the points on the beam surface
        ComputeBphiAndBu(mBphi, mBu, m_beamLength(ii));

        Hb = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * ii + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * ii + 1]];

        Vector thisBeamDisp(12);
        for (int jj = 0; jj < 6; jj++)
        {
            thisBeamDisp(jj) = bDisp(6 * beamNodeInA1 + jj);
            thisBeamDisp(jj + 6) = bDisp(6 * beamNodeInA2 + jj);
        }

        ptDisp = Hb * thisBeamDisp;

        Vector tempVec(3); double temp;
        switch (flag)
        {
        case 0:
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = ptDisp(jj);
            break;
        case 1:
            ptDisp = Transpose(3, 3, mQc) * ptDisp;
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = ptDisp(jj);
            break;
        case 2:
            temp = ptDisp(0)*c1(0) + ptDisp(1)*c1(1) + ptDisp(2)*c1(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*c1(jj);
            break;
        case 3:
            tempVec = cos(m_beam_theta(ii)) * c2 + sin(m_beam_theta(ii)) * c3;
            temp = ptDisp(0)*tempVec(0) + ptDisp(1)*tempVec(1) + ptDisp(2)*tempVec(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*tempVec(jj);
            break;
        case 4:
            tempVec = cos(m_beam_theta(ii)) * c3 - sin(m_beam_theta(ii)) * c2;
            temp = ptDisp(0)*tempVec(0) + ptDisp(1)*tempVec(1) + ptDisp(2)*tempVec(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*tempVec(jj);
            break;
        }
    }

    return res;
}

Vector EmbeddedBeamInterfaceP::GetInteractionPtForce(int flag)
{
    Vector res(3 * m_numEmbeddedPoints);
    Vector ptForces(3);
    Matrix Hf(3, 12);
    Vector c1(3), c2(3), c3(3);

    // update local coordinate system
    for (int ii = 0; ii < 3; ii++)
    {
        c1(ii) = mQc(ii, 0);
        c2(ii) = mQc(ii, 1);
        c3(ii) = mQc(ii, 2);
    }

    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));

        // update the matrices that define kinematics of the points on the beam surface
        ComputeHf(Hf, m_beam_theta(ii));

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * ii + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * ii + 1]];

        Vector thisBeamForce(12);
        for (int jj = 0; jj < 6; jj++)
        {
            thisBeamForce(jj) = m_Lambda(6 * beamNodeInA1 + jj);
            thisBeamForce(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
        }

        ptForces = Hf * thisBeamForce * m_area(ii);

        Vector tempVec(3); double temp;
        switch (flag)
        {
        case 0:
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = ptForces(jj);
            break;
        case 1:
            ptForces = Transpose(3, 3, mQc) * ptForces;
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = ptForces(jj);
            break;
        case 2:
            temp = ptForces(0)*c1(0) + ptForces(1)*c1(1) + ptForces(2)*c1(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*c1(jj);
            break;
        case 3:
            tempVec = cos(m_beam_theta(ii)) * c2 + sin(m_beam_theta(ii)) * c3;
            temp = ptForces(0)*tempVec(0) + ptForces(1)*tempVec(1) + ptForces(2)*tempVec(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*tempVec(jj);
            break;
        case 4:
            tempVec = cos(m_beam_theta(ii)) * c3 - sin(m_beam_theta(ii)) * c2;
            temp = ptForces(0)*tempVec(0) + ptForces(1)*tempVec(1) + ptForces(2)*tempVec(2);
            for (int jj = 0; jj < 3; jj++)
                res(3 * ii + jj) = temp*tempVec(jj);
            break;
        }
    }

    return res;
}

Vector EmbeddedBeamInterfaceP::GetElementalForce(int flag)
{
    Vector res;

    if (flag == 1)
    {
        res.resize(12 * (m_numBeamNodes - 1));
        Matrix eleB(12, 12);
        Matrix Hf(3, 12), Hu(3, 12);
        Vector c2(3), c3(3);

        // update local coordinate system
        for (int ii = 0; ii < 3; ii++)
        {
            c2(ii) = mQc(ii, 1);
            c3(ii) = mQc(ii, 2);
        }

        int curBeamTag = theBeamTags[0];
        int curBeamCount = 0;
        for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
        {
            if (theBeamTags[ii] != curBeamTag)
            {
                curBeamTag = theBeamTags[ii];
                Vector thisBeamLambda(12);

                int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * (ii - 1) + 0]];
                int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * (ii - 1) + 1]];

                for (int jj = 0; jj < 6; jj++)
                {
                    thisBeamLambda(jj) = m_Lambda(6 * beamNodeInA1 + jj);
                    thisBeamLambda(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
                }

                Vector thisBeamForce = -1.0  * eleB * thisBeamLambda;

                for (int jj = 0; jj < 6; jj++)
                {
                    res(12 * curBeamCount + jj + 0) = thisBeamForce(jj);
                    res(12 * curBeamCount + jj + 6) = thisBeamForce(jj + 6);
                }

                eleB.Zero();
                curBeamCount++;
            }

            updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));
            ComputeHf(Hf, m_beam_theta(ii));
            ComputeBphiAndBu(mBphi, mBu, m_beamLength(ii));

            Vector c2(3), c3(3);
            for (int jj = 0; jj < 3; jj++)
            {
                c2(jj) = mQc(jj, 1);
                c3(jj) = mQc(jj, 2);
            }

            Hu = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
            eleB += Transpose(3, 12, Hu) * Hf * m_area(ii);
        }
        Vector thisBeamLambda(12);

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * (m_numEmbeddedPoints - 1) + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * (m_numEmbeddedPoints - 1) + 1]];

        for (int jj = 0; jj < 6; jj++)
        {
            thisBeamLambda(jj) = m_Lambda(6 * beamNodeInA1 + jj);
            thisBeamLambda(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
        }

        Vector thisBeamForce = -1.0  * eleB * thisBeamLambda;

        for (int jj = 0; jj < 6; jj++)
        {
            res(12 * curBeamCount + jj + 0) = thisBeamForce(jj);
            res(12 * curBeamCount + jj + 6) = thisBeamForce(jj + 6);
        }
    }
    else if (flag == 2)
    {
        res = mA * m_Lambda;
    }
    else if (flag == 3)
    {
        res.resize(12 * (m_numBeamNodes - 1));
        Matrix eleB(12, 12);
        Matrix Hf(3, 12), Hu(3, 12);
        Vector c2(3), c3(3);

        // update local coordinate system
        for (int ii = 0; ii < 3; ii++)
        {
            c2(ii) = mQc(ii, 1);
            c3(ii) = mQc(ii, 2);
        }

        int curBeamTag = theBeamTags[0];
        int curBeamCount = 0;
        for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
        {
            if (theBeamTags[ii] != curBeamTag)
            {
                curBeamTag = theBeamTags[ii];
                Vector thisBeamLambda(12);

                int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * (ii - 1) + 0]];
                int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * (ii - 1) + 1]];

                for (int jj = 0; jj < 6; jj++)
                {
                    thisBeamLambda(jj) = m_Lambda(6 * beamNodeInA1 + jj);
                    thisBeamLambda(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
                }

                Vector thisBeamForce = -1.0  * eleB * thisBeamLambda;

                for (int jj = 0; jj < 4; jj++)
                {
                    Vector temp1(3), temp2(3); Matrix Q(3, 3);
                    for (int kk = 0; kk < 3; kk++)
                        temp1(kk) = thisBeamForce(3 * jj + kk);
                    if (jj < 2)
                        Q = Transpose(3, 3, mQa);
                    else
                        Q = Transpose(3, 3, mQb);
                    temp2 = Q * temp1;
                    for (int kk = 0; kk < 3; kk++)
                        thisBeamForce(3 * jj + kk) = temp2(kk);
                }

                for (int jj = 0; jj < 6; jj++)
                {
                    res(12 * curBeamCount + jj + 0) = thisBeamForce(jj);
                    res(12 * curBeamCount + jj + 6) = thisBeamForce(jj + 6);
                }

                eleB.Zero();
                curBeamCount++;
            }

            updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii), m_beamLength(ii));
            ComputeHf(Hf, m_beam_theta(ii));
            ComputeBphiAndBu(mBphi, mBu, m_beamLength(ii));

            Vector c2(3), c3(3);
            for (int jj = 0; jj < 3; jj++)
            {
                c2(jj) = mQc(jj, 1);
                c3(jj) = mQc(jj, 2);
            }

            Hu = mBu - (m_beam_radius*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
            eleB += Transpose(3, 12, Hu) * Hf * m_area(ii);
        }
        Vector thisBeamLambda(12);

        int beamNodeInA1 = m_beamNodeMap[beamNodeTags[2 * (m_numEmbeddedPoints - 1) + 0]];
        int beamNodeInA2 = m_beamNodeMap[beamNodeTags[2 * (m_numEmbeddedPoints - 1) + 1]];

        for (int jj = 0; jj < 6; jj++)
        {
            thisBeamLambda(jj) = m_Lambda(6 * beamNodeInA1 + jj);
            thisBeamLambda(jj + 6) = m_Lambda(6 * beamNodeInA2 + jj);
        }

        Vector thisBeamForce = -1.0  * eleB * thisBeamLambda;

        for (int jj = 0; jj < 4; jj++)
        {
            Vector temp1(3), temp2(3); Matrix Q(3, 3);
            for (int kk = 0; kk < 3; kk++)
                temp1(kk) = thisBeamForce(3 * jj + kk);
            if (jj < 2)
                Q = Transpose(3, 3, mQa);
            else
                Q = Transpose(3, 3, mQb);
            temp2 = Q * temp1;
            for (int kk = 0; kk < 3; kk++)
                thisBeamForce(3 * jj + kk) = temp2(kk);
        }

        for (int jj = 0; jj < 6; jj++)
        {
            res(12 * curBeamCount + jj + 0) = thisBeamForce(jj);
            res(12 * curBeamCount + jj + 6) = thisBeamForce(jj + 6);
        }
    }
    else
    {
        opserr << "Unknown request ..." << endln;
        return Vector();
    }
    return res;
}
