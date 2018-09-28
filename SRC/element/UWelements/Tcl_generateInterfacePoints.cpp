#include <stdlib.h>
#include <string.h>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>

#include <Element.h>
#include <OPS_Stream.h>
#include <Domain.h>
#include <Node.h>
#include <ElementIter.h>
#include <ElasticBeam3d.h>
#include <CrdTransf.h>
#include <TclModelBuilder.h>
#include <elementAPI.h>
#include <FileStream.h>

#include "Tcl_generateInterfacePoints.h"

// Need the domain to get access to the elements
#ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
extern PartitionedDomain theDomain;
#else
extern Domain theDomain;
#endif

int
TclCommand_GenerateInterfacePoints(ClientData clientData, Tcl_Interp *interp, int argc,
    TCL_Char **argv)
{
    static int tagOffset = 0;
    // argument parsing variables
    int numArgsRemaining = argc - 1;
    int curArgPos        = 1;
    int res              = 0;

    int beamTag;
    double radius   = 1.0;
    int nP          = 20;
    int nL          = 20;
    int shape       = 1;
    int matTag      = 0;
    int baseNodeTag = 0;
    double penaltyParam = 1.0e12;

    bool radiusSet            = false;
    bool nPset                = false;
    bool nLset                = false;
    bool shapeSet             = false;
    bool solidEleSetDefined   = false;
    bool beamEleSetDefined    = false;
    bool crdTransfDefined     = false;
    bool debugFlag            = false;
    bool writeCoords          = false;
    bool writeConnectivity    = false;
    bool writeQuadInfo        = false;
    bool beamConnected        = false;
    bool baseNodeSet          = false;
    bool readingBeamEleTag    = false;
    bool readingSolidEleTag   = false;
    bool readingConBeamEleTag = false;

    int  lagrangeTag   = 8; // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty

    const char *         crdsFN;
    const char *         connectivityFN;
    const char *         quadNodeFN;
    const char *         quadElemFN;
    std::vector <int>    solidEleTags;
    std::vector <int>    beamEleTags;
    std::vector <int>    connectedBeamsTags;
    std::vector <int>    solidEleTagsInContact;
    std::vector <int>    beamEleTagsInContact;
    std::vector <double> contactPt_xi;
    std::vector <double> contactPt_eta;
    std::vector <double> contactPt_zeta;
    std::vector <double> contactPt_rho;
    std::vector <double> contactPt_theta;
    std::vector <double> contactPt_area;
    std::vector <double> contactPt_length;
    std::set    <int>    eleTagsInContact_unique;

    int transfTag = 1;

    if (numArgsRemaining < 2)
    {
        opserr << "Minimum requirements not met." << endln;
        return -1;
    }

    // Parse the arguments
    while (numArgsRemaining > 0)
    {
        if ((strcmp(argv[curArgPos], "-shape") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            numArgsRemaining--; curArgPos++;
            if (strcmp(argv[curArgPos], "circle") == 0)
            {
                shapeSet = true;
                shape = 1;
                numArgsRemaining--; curArgPos++;
            }
            else if (strcmp(argv[curArgPos], "square") == 0)
            {
                shapeSet = true;
                shape = 2;
                numArgsRemaining--; curArgPos++;
                opserr << "The procedure for square shapes is not implemented yet." << endln;
            }
            else if (strcmp(argv[curArgPos], "rectangle") == 0)
            {
                shapeSet = true;
                shape = 4;
                numArgsRemaining--; curArgPos++;
                opserr << "The procedure for rectangle shapes is not implemented yet." << endln;
            }
            else if (strcmp(argv[curArgPos], "hexagon") == 0)
            {
                shapeSet = true;
                shape = 6;
                numArgsRemaining--; curArgPos++;
                opserr << "The procedure for hexagon shapes is not implemented yet." << endln;
            }
            else {
                numArgsRemaining--; curArgPos++;
                opserr << "Unknown shape." << endln;
            }

        }
        else if ((strcmp(argv[curArgPos], "-radius") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            radiusSet = true;
            numArgsRemaining--; curArgPos++;
            radius = strtod(argv[curArgPos], NULL);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-crdTransf") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            crdTransfDefined = true;
            numArgsRemaining--; curArgPos++;
            transfTag = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-nP") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            nPset = true;
            numArgsRemaining--; curArgPos++;
            nP = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-nL") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            nLset = true;
            numArgsRemaining--; curArgPos++;
            nL = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-penaltyParam") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            numArgsRemaining--; curArgPos++;
            penaltyParam = strtod(argv[curArgPos], NULL);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-solidEle") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = true;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            solidEleSetDefined = true;
            numArgsRemaining--; curArgPos++;
            solidEleTags.push_back(strtol(argv[curArgPos], NULL, 10));
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-solidEleRange") == 0) && numArgsRemaining >= 3)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            solidEleSetDefined = true;
            numArgsRemaining--; curArgPos++;
            int startEle = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
            int endEle = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;

            for (int ii = startEle; ii <= endEle; ii++)
                solidEleTags.push_back(ii);
        }
        else if ((strcmp(argv[curArgPos], "-beamEle") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = true;
            beamEleSetDefined = true;
            numArgsRemaining--; curArgPos++;
            beamEleTags.push_back(strtol(argv[curArgPos], NULL, 10));
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-beamEleRange") == 0) && numArgsRemaining >= 3)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            beamEleSetDefined = true;
            numArgsRemaining--; curArgPos++;
            int startEle = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
            int endEle = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;

            for (int ii = startEle; ii <= endEle; ii++)
                beamEleTags.push_back(ii);
        }
        else if ((strcmp(argv[curArgPos], "-connectedBeams") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag     = false;
            readingConBeamEleTag = true;
            readingBeamEleTag = false;
            beamConnected     = true;
            numArgsRemaining--; curArgPos++;
            connectedBeamsTags.push_back(strtol(argv[curArgPos], NULL, 10));
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-file") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            writeCoords = true;
            numArgsRemaining--; curArgPos++;
            crdsFN = argv[curArgPos];
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-connectivity") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            writeConnectivity = true;
            numArgsRemaining--; curArgPos++;
            connectivityFN = argv[curArgPos];
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-penalty") == 0) || (strcmp(argv[curArgPos], "-lPenalty") == 0))
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 0;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-AL") == 0)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 1;
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-lagrange") == 0) || (strcmp(argv[curArgPos], "-gLagrange") == 0))
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 2;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-embed") == 0)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 3;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-gPenalty") == 0)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 4;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-gAL") == 0)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 5;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-contact") == 0)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = -1;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-ElastoPlastic") == 0)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 6;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-gContact") == 0)
        {
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 7;
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-mat") == 0) && (numArgsRemaining >= 2))
        {
            // This one is for the elastoplastic interface
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            numArgsRemaining--; curArgPos++;
            matTag = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-quadFileInfo") == 0) && (numArgsRemaining >= 3))
        {
            // This one is for the elastoplastic interface
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            writeQuadInfo = true;

            numArgsRemaining--; curArgPos++;
            quadNodeFN = argv[curArgPos];
            numArgsRemaining--; curArgPos++;
            quadElemFN = argv[curArgPos];
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-baseNode") == 0) && (numArgsRemaining >= 2))
        {
            // This one is for the elastoplastic interface
            readingSolidEleTag = false;
            readingConBeamEleTag = false;
            readingBeamEleTag = false;
            baseNodeSet = true;

            numArgsRemaining--; curArgPos++;
            baseNodeTag = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
        }
        else {
            if (readingSolidEleTag)
            {
                solidEleTags.push_back(strtol(argv[curArgPos], NULL, 10));
                numArgsRemaining--; curArgPos++;
                continue;
            }
            else if (readingConBeamEleTag)
            {
                connectedBeamsTags.push_back(strtol(argv[curArgPos], NULL, 10));
                numArgsRemaining--; curArgPos++;
                continue;
            }
            else if (readingBeamEleTag)
            {
                beamEleTags.push_back(strtol(argv[curArgPos], NULL, 10));
                numArgsRemaining--; curArgPos++;
                continue;
            }
            else {
                opserr << "Unknown argument : " << argv[curArgPos] << endln;
                res = -1;
                break;
            }
        }
    }

    // Get the coordinate transformation object
    CrdTransf *theTransf;
    if (crdTransfDefined)
         theTransf = OPS_GetCrdTransf(transfTag)->getCopy3d();
    else
    {
        opserr << "Coordinate transformation object tag not defined." << endln;
        return -1;
    }

    // Loop over all the elements and check if they are supported. If yes add them to the element pool
    ElementIter& theElements = theDomain.getElements();
    Element *theElement;
    // Record maximum element tag for generating new interface tags
    int maxTag = 0;
    while ((theElement = theElements()) != 0)
    {
        if (!solidEleSetDefined)
        {
            const char* eleClassType = theElement->getClassType();
            if ((strcmp(eleClassType, "Brick") == 0) ||
                (strcmp(eleClassType, "BbarBrick") == 0) ||
                (strcmp(eleClassType, "SSPbrick") == 0) ||
                (strcmp(eleClassType, "SSPbrickUP") == 0) ||
                (strcmp(eleClassType, "BbarBrick") == 0))
                solidEleTags.push_back(theElement->getTag());
        }
        if (theElement->getTag() > maxTag)
            maxTag = theElement->getTag();
    }
    int startTag = maxTag + 1;

    FileStream crdsFile, quadNodes, quadElems;
    if (writeCoords)
        crdsFile.setFile(crdsFN, APPEND);
    if (writeQuadInfo)
    {
        quadNodes.setFile(quadNodeFN, APPEND);
        quadElems.setFile(quadElemFN, APPEND);
    }

    // Define variables needed
    int contactElemCount = 0, contactPointCount = 0;
    double area;

    // Loop over beam elements
    for (int beamCount = 0; beamCount < beamEleTags.size(); beamCount++)
    {
        // Get the beam element and check if supported
        beamTag = beamEleTags[beamCount];
        theElement = theDomain.getElement(beamTag);
        const char* beamEleClassType = theElement->getClassType();
        if (!((strcmp(beamEleClassType, "ElasticBeam3d") == 0) ||
            (strcmp(beamEleClassType, "DispBeamColumn3d") == 0) ||
            (strcmp(beamEleClassType, "DispBeamColumn3dThermal") == 0) ||
            (strcmp(beamEleClassType, "ElasticTimoshenkoBeam3d") == 0) ||
            (strcmp(beamEleClassType, "ElasticForceBeamColumn3d") == 0) ||
            (strcmp(beamEleClassType, "ForceBeamColumn3d") == 0) ||
            (strcmp(beamEleClassType, "ForceBeamColumn3d") == 0)))
        {
            opserr << "Beam element of type " << theElement->getClassType() << "not supported." << endln;
            return -1;
        }
        
        Node** solidNodesPtr;
        Node** beamNodePtr;

        // get beam node locations (considering initial displacements)
        beamNodePtr = theDomain.getElement(beamTag)->getNodePtrs();
        double L1x = beamNodePtr[0]->getCrds()(0) + beamNodePtr[0]->getDisp()(0);
        double L1y = beamNodePtr[0]->getCrds()(1) + beamNodePtr[0]->getDisp()(1);
        double L1z = beamNodePtr[0]->getCrds()(2) + beamNodePtr[0]->getDisp()(2);
        double L2x = beamNodePtr[1]->getCrds()(0) + beamNodePtr[1]->getDisp()(0);
        double L2y = beamNodePtr[1]->getCrds()(1) + beamNodePtr[1]->getDisp()(1);
        double L2z = beamNodePtr[1]->getCrds()(2) + beamNodePtr[1]->getDisp()(2);

        // Use either the transformation object (next block of code) or calculate rotations and everything (below)

        // // calculate the rotation information(the centerline is calculated using a
        // // rotation operation)
        // double ux = L2x - L1x;
        // double uy = L2y - L1y;
        // double uz = L2z - L1z;
        // double L = sqrt(pow(ux, 2.0) + pow(uy, 2) + pow(uz, 2));
        // ux /= L;
        // uy /= L;
        // uz /= L;
        // 
        // double sintheta, costheta, tantheta;
        // if (ux == 0)
        // {
        // 	sintheta = 0;
        // 	costheta = 1;
        // }
        // else
        // {
        // 	tantheta = uy / ux;
        // 	costheta = 1 / sqrt(1.0 + pow(tantheta, 2.0));
        // 	sintheta = tantheta * costheta;
        // }
        // 
        // double singamma = uz;
        // double cosgamma = sqrt(1 - pow(singamma, 2));
        // 
        // // create the pile points
        // Vector t(nP);
        // for (int ii = 0; ii < nP; ii++)
        // {
        // 	t(ii) = 0 + (double)ii / nP * 2.0 * PI;
        // }
        // 
        // Vector cX(nL*nP), cY(nL*nP), cZ(nL*nP);
        // for (int ii = 0; ii < nL; ii++)
        // 	for (int jj = 0; jj < nP; jj++)
        // 	{
        // 		cX(ii * nP + jj) = radius * cos(t(jj));
        // 		cY(ii * nP + jj) = radius * sin(t(jj));
        // 		cZ(ii * nP + jj) = 0 + (double)ii / (nL - 1) * L;
        // 	}
        // 
        // Vector cXr = cX*costheta*singamma - cY*sintheta + cZ*costheta*cosgamma + L1x;
        // Vector cYr = cX*sintheta*singamma + cY*costheta + cZ*sintheta*cosgamma + L1y;
        // Vector cZr = cZ*singamma - cX*cosgamma + L1z;


        // use the geometric transformation object to generate beam surface points
        if (theTransf->initialize(beamNodePtr[0], beamNodePtr[1]))
        {
            opserr << "generateInterfacePoints: Error initializing coordinate transformation";
            return -1;
        }
        Vector loc_x(3), loc_y(3), loc_z(3);
        theTransf->getLocalAxes(loc_x, loc_y, loc_z);

        // create the pile points
        Vector t(nP);
        for (int ii = 0; ii < nP; ii++)
        {
            t(ii) = (0.5 + (double)ii) / nP * 2.0 * PI;
        }

        double L = sqrt(pow(L2x - L1x, 2.0) + pow(L2y - L1y, 2.0) + pow(L2z - L1z, 2.0));

        Vector cX(nL*nP), cY(nL*nP), cZ(nL*nP);
        Vector loc_rho(nL*nP), loc_theta(nL*nP);
        for (int ii = 0; ii < nL; ii++)
            for (int jj = 0; jj < nP; jj++)
            {
                cX(ii * nP + jj) = (1.0 + 2.0 * ii) / (2.0 * nL) * L;
                cY(ii * nP + jj) = radius * cos(t(jj));
                cZ(ii * nP + jj) = radius * sin(t(jj));

                loc_rho(ii * nP + jj) = 2.0 *  (1.0 + 2.0 * ii) / (2.0 * nL) - 1.0;
                loc_theta(ii * nP + jj) = t(jj);
            }

        Vector cXr = cX*loc_x(0) + cY*loc_y(0) + cZ*loc_z(0) + L1x;
        Vector cYr = cX*loc_x(1) + cY*loc_y(1) + cZ*loc_z(1) + L1y;
        Vector cZr = cX*loc_x(2) + cY*loc_y(2) + cZ*loc_z(2) + L1z;

        area = 2.0 * PI * radius * L / nP / nL;
        
        if (writeQuadInfo)
        {            
            double qt;

            for (int ii = 0; ii < nL + 1; ii++)
                for (int jj = 0; jj < nP + 1; jj++)
                {
                    qt = (double)jj / nP * 2.0 * PI;
                    double cXq = (double)ii / nL * L;
                    double cYq = radius * cos(qt);
                    double cZq = radius * sin(qt);

                    double qcXr = cXq*loc_x(0) + cYq*loc_y(0) + cZq*loc_z(0) + L1x;
                    double qcYr = cXq*loc_x(1) + cYq*loc_y(1) + cZq*loc_z(1) + L1y;
                    double qcZr = cXq*loc_x(2) + cYq*loc_y(2) + cZq*loc_z(2) + L1z;

                    quadNodes << tagOffset + beamCount * (nL + 1)*(nP + 1) + ii * (nP + 1) + jj << " " << qcXr << " " << qcYr << " " << qcZr << endln;
                }
        }

        if (debugFlag)
            for (int ii = 0; ii < nL; ii++)
                for (int jj = 0; jj < nP; jj++)
                    opserr << "point " << ii * nP + jj + 1 << " : " << cXr(ii * nP + jj) << " " << cYr(ii * nP + jj) << " " << cZr(ii * nP + jj) << endln;


        Vector tempX(8), tempY(8), tempZ(8);
        bool contactElemFlag = false;
        double xi, eta, zeta;
        bool inBounds = false;
        // find intersections and create interface elements.
        for (int jj = 0; jj < nP*nL; jj++)
        {
            for (int ii = 0; ii < solidEleTags.size(); ii++)
            {
                theElement = theDomain.getElement(solidEleTags[ii]);
                solidNodesPtr = theElement->getNodePtrs();
                for (int kk = 0; kk < 8; kk++)
                {
                    tempX(kk) = solidNodesPtr[kk]->getCrds()(0) + solidNodesPtr[kk]->getDisp()(0);
                    tempY(kk) = solidNodesPtr[kk]->getCrds()(1) + solidNodesPtr[kk]->getDisp()(1);
                    tempZ(kk) = solidNodesPtr[kk]->getCrds()(2) + solidNodesPtr[kk]->getDisp()(2);
                }
                invIsoMapping(tempX, tempY, tempZ, cXr(jj), cYr(jj), cZr(jj), xi, eta, zeta, inBounds);
                if (inBounds)
                {
                    contactElemFlag = true;
                    if (debugFlag)
                        opserr << "Beam tag : " << beamTag << ", Solid tag : " << solidEleTags[ii] << ", Real Coordinates = (" << cXr(jj) << "," << cYr(jj) << "," << cZr(jj) << "), Iso Coordinates = (" << xi << "," << eta << "," << zeta << ")" << endln;

                    eleTagsInContact_unique.insert(solidEleTags[ii]);
                    maxTag++;

                    if (lagrangeTag > 1)
                    {
                        // This is for a global element
                        solidEleTagsInContact.push_back(solidEleTags[ii]);
                        beamEleTagsInContact.push_back(beamEleTags[beamCount]);
                        contactPt_xi.push_back(xi);
                        contactPt_eta.push_back(eta);
                        contactPt_zeta.push_back(zeta);
                        contactPt_rho.push_back(loc_rho(jj));
                        contactPt_theta.push_back(loc_theta(jj));
                        contactPt_area.push_back(area);
                        contactPt_length.push_back(L);
                    }
                    else
                    {
                        // This is for a local element

                        // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange
                        if (lagrangeTag == 0) {
                            // theElement = new EmbeddedBeamInterface(maxTag, beamTag, solidEleTags[ii], transfTag, loc_rho(jj), loc_theta(jj), xi, eta, zeta, radius, area, writeConnectivity, connectivityFN);
			}
                        else if (lagrangeTag == 1)
 			{
                            // theElement = new EmbeddedBeamInterfaceAL(maxTag, beamTag, solidEleTags[ii], transfTag, loc_rho(jj), loc_theta(jj), xi, eta, zeta, radius, area);
			}
                        if (lagrangeTag == -1)
			{
                            //theElement = new EmbeddedBeamContact(maxTag, beamTag, solidEleTags[ii], transfTag, loc_rho(jj), loc_theta(jj), xi, eta, zeta, radius, area);
			}

                        // Create the return of the tcl command
                        char buffer[40];
                        sprintf(buffer, "%10i", maxTag);
                        Tcl_AppendResult(interp, buffer, NULL);

                        // Add the created element to the domain
                        if (theElement != 0)
                        {
                            if (theDomain.addElement(theElement) == false)
                            {
                                opserr << "WARNING could not add element with tag: " << theElement->getTag() << " and of type: "
                                    << theElement->getClassType() << " to the Domain\n";
                                delete theElement;
                                return -1;
                            }
                        }

                    }

                    if (writeCoords)
                        crdsFile << maxTag << "\t" << cXr(jj) << "\t" << cYr(jj) << "\t" << cZr(jj) << "\n";
                    
                    if (writeQuadInfo)
                    {
                        int setNum = (int)(jj / nP);
                        int ptNum = (int)(jj % nP);
                        int pt1 = tagOffset + beamCount * (nL + 1) * (nP + 1) + setNum * (nP + 1) + ptNum;
                        int pt3 = tagOffset + beamCount * (nL + 1) * (nP + 1) + (setNum + 1) * (nP + 1) + ptNum;

                        quadElems << contactPointCount << " " << pt1 << " " << pt1 + 1 << " " << pt3 + 1 << " " << pt3 << endln;
                    }

                    contactPointCount++;
                    break;
                }
            }
        }
    }

    tagOffset += contactPointCount;

    if (contactPointCount > 0)
    {
        // Create the global element

        // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange
        if (lagrangeTag > 1)
        {
            if (lagrangeTag == 2)
                theElement = new EmbeddedBeamInterfaceL(maxTag, beamEleTagsInContact, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, contactPt_area, contactPt_length, writeConnectivity, connectivityFN);
            else if (lagrangeTag == 3)
	    {
            //    if (beamConnected)
            //        theElement = new EmbeddedBeamInterfaceAL2(maxTag, beamEleTagsInContact, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, contactPt_area, writeConnectivity, connectivityFN);
            //    else
            //        theElement = new EmbeddedBeamInterfaceAL2(maxTag, beamEleTagsInContact, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, contactPt_area, writeConnectivity, connectivityFN);
            }
            else if (lagrangeTag == 4)
                theElement = new EmbeddedBeamInterfaceP(maxTag, beamEleTagsInContact, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, contactPt_area, contactPt_length, penaltyParam, writeConnectivity, connectivityFN);
            else if (lagrangeTag == 5)
            {
                // theElement = new EmbeddedBeamInterfaceAL2(maxTag, beamTag, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, area);
            //    theElement = new EmbeddedBeamInterfaceP2(maxTag, beamTag, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, area);
            }
            else if (lagrangeTag == 6)
            {
                theElement = new EmbeddedEPBeamInterface(maxTag, beamEleTagsInContact, solidEleTagsInContact, transfTag, matTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, contactPt_area, contactPt_length, writeConnectivity, connectivityFN);
            }
            else if (lagrangeTag == 7)
            {
            //    theElement = new EmbeddedBeamContactInterfaceP(maxTag, beamTag, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, area);
            // theElement = new EmbedBeamContactG(maxTag, beamTag, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, radius, area);
            }

            // Create the return of the tcl command
            char buffer[40];
            sprintf(buffer, "%10i", maxTag);
            Tcl_AppendResult(interp, buffer, NULL);

            // Add the created element to the domain
            if (theElement != 0)
            {
                if (theDomain.addElement(theElement) == false)
                {
                    opserr << "WARNING could not add element with tag: " << theElement->getTag() << " and of type: "
                        << theElement->getClassType() << " to the Domain\n";
                    delete theElement;
                    return -1;
                }
            }
        }
    }

    if (writeCoords)
        crdsFile.close();


    opserr << "Number of elements in contact: " << (int)eleTagsInContact_unique.size() << ", number of contact points: " << contactPointCount << endln;
    opserr << "EmbeddedBeamInterface Elements " << startTag << " to " << maxTag << " were created." << endln;

	if (debugFlag)
	{
		if (radiusSet)
			opserr << "Radius is set to " << radius << endln;
		else
			opserr << "No radius is set!!! Assuming 1.0" << endln;
		if (nPset)
			opserr << "nP is set to " << nP << endln;
		else
			opserr << "nP not set!!! Assuming " << nP << endln;
		if (nLset)
			opserr << "nL is set to " << nL << endln;
		else
			opserr << "nL not set!!! Assuming " << nL << endln;
		if (shapeSet)
			opserr << "Shape is set to " << shape << endln;
		else
			opserr << "Default shape (circle) is assumed." << endln;
		if (solidEleTags.size() > 0)
			opserr << "Here are the elements: \n     { ";
		for (int ii = 0; ii < solidEleTags.size(); ii++)
			opserr << solidEleTags[ii] << " ";
		if (solidEleTags.size() > 0)
			opserr << "}" << endln;
	}


	return 0;

}




























int
TclCommand_GenerateToeInterfacePoints(ClientData clientData, Tcl_Interp *interp, int argc,
    TCL_Char **argv)
{
    int numArgsRemaining = argc - 1;
    int curArgPos        = 1;
    int res              = 0;

    int beamTag, toeNodeTag;
    double radius = 1.0;
    int nP        = 100;
    int nR        = 20;
    int shape     = 1;

    bool radiusSet     = false;
    bool nPset         = false;
    bool nRset         = false;
    bool readingSolidEleTag = false;
    bool shapeSet      = false;
    bool solidEleSetDefined = false;
    bool debugFlag     = false;
    bool writeCoords   = false;
    int  lagrangeTag   = 3; // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange

    const char * crdsFN;
    std::vector <int> solidEleTags;
    std::vector <int> solidEleTagsInContact;
    std::vector <double> contactPt_xi;
    std::vector <double> contactPt_eta;
    std::vector <double> contactPt_zeta;
    std::vector <double> contactPt_rho;
    std::vector <double> contactPt_theta;
    std::vector <double> contactPt_radius;
    std::set    <int> eleTagsInContact_unique;

    int transfTag = 1;

    if (numArgsRemaining < 3)
    {
        opserr << "Need beam tag, geometric transformation tag and the toe node tag." << endln;
        return -1;
    }
    else
    {
        beamTag = strtol(argv[curArgPos], NULL, 10);
        numArgsRemaining--; curArgPos++;
        transfTag = strtol(argv[curArgPos], NULL, 10);
        numArgsRemaining--; curArgPos++;
        toeNodeTag = strtol(argv[curArgPos], NULL, 10);
        numArgsRemaining--; curArgPos++;
    }


    while (numArgsRemaining > 0)
    {
        if ((strcmp(argv[curArgPos], "-shape") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            numArgsRemaining--; curArgPos++;
            if (strcmp(argv[curArgPos], "circle") == 0)
            {
                shapeSet = true;
                shape = 1;
                numArgsRemaining--; curArgPos++;
            }
            else if (strcmp(argv[curArgPos], "square") == 0)
            {
                shapeSet = true;
                shape = 2;
                numArgsRemaining--; curArgPos++;
                opserr << "The procedure for square shapes not implemented yet." << endln;
            }
            else if (strcmp(argv[curArgPos], "rectangle") == 0)
            {
                shapeSet = true;
                shape = 4;
                numArgsRemaining--; curArgPos++;
                opserr << "The procedure for rectangle shapes not implemented yet." << endln;
            }
            else if (strcmp(argv[curArgPos], "hexagon") == 0)
            {
                shapeSet = true;
                shape = 6;
                numArgsRemaining--; curArgPos++;
                opserr << "The procedure for hexagon shapes not implemented yet." << endln;
            }
            else
            {
                numArgsRemaining--; curArgPos++;
                opserr << "Unknown shape." << endln;
            }

        }
        else if ((strcmp(argv[curArgPos], "-radius") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            radiusSet = true;
            numArgsRemaining--; curArgPos++;
            radius = strtod(argv[curArgPos], NULL);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-nP") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            nPset = true;
            numArgsRemaining--; curArgPos++;
            nP = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-nR") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            nRset = true;
            numArgsRemaining--; curArgPos++;
            nR = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-ele") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = true;
            solidEleSetDefined = true;
            numArgsRemaining--; curArgPos++;
            solidEleTags.push_back(strtol(argv[curArgPos], NULL, 10));
            numArgsRemaining--; curArgPos++;
        }
        else if ((strcmp(argv[curArgPos], "-eleRange") == 0) && numArgsRemaining >= 3)
        {
            readingSolidEleTag = false;
            solidEleSetDefined = true;
            numArgsRemaining--; curArgPos++;
            int startEle = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;
            int endEle = strtol(argv[curArgPos], NULL, 10);
            numArgsRemaining--; curArgPos++;

            for (int ii = startEle; ii <= endEle; ii++)
                solidEleTags.push_back(ii);
        }
        else if ((strcmp(argv[curArgPos], "-file") == 0) && numArgsRemaining >= 2)
        {
            readingSolidEleTag = false;
            writeCoords = true;
            numArgsRemaining--; curArgPos++;
            crdsFN = argv[curArgPos];
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-penalty") == 0)
        {
            readingSolidEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange
            lagrangeTag = 0;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-AL") == 0)
        {
            readingSolidEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange
            lagrangeTag = 1;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-lagrange") == 0)
        {
            readingSolidEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange
            lagrangeTag = 2;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-embed") == 0)
        {
            readingSolidEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange
            lagrangeTag = 3;
            numArgsRemaining--; curArgPos++;
        }
        else if (strcmp(argv[curArgPos], "-gPenalty") == 0)
        {
            readingSolidEleTag = false;
            // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange, 4 = global penalty
            lagrangeTag = 4;
            numArgsRemaining--; curArgPos++;
        }
        else
        {
            if (readingSolidEleTag)
            {
                solidEleTags.push_back(strtol(argv[curArgPos], NULL, 10));
                numArgsRemaining--; curArgPos++;
                continue;
            }
            else
            {
                opserr << "Unknown argument (" << argv[curArgPos] <<") ..." << endln;
                res = -1;
                break;
            }
        }
    }

    CrdTransf *theTransf = OPS_GetCrdTransf(transfTag)->getCopy3d();
    Element *theElement;
    theElement = theDomain.getElement(beamTag);
    const char* beamEleClassType = theElement->getClassType();
    if (!((strcmp(beamEleClassType, "ElasticBeam3d") == 0) ||
        (strcmp(beamEleClassType, "DispBeamColumn3d") == 0) ||
        (strcmp(beamEleClassType, "DispBeamColumn3dThermal") == 0) ||
        (strcmp(beamEleClassType, "ElasticTimoshenkoBeam3d") == 0) ||
        (strcmp(beamEleClassType, "ElasticForceBeamColumn3d") == 0) ||
        (strcmp(beamEleClassType, "ForceBeamColumn3d") == 0) ||
        (strcmp(beamEleClassType, "ForceBeamColumn3d") == 0)))
    {
        opserr << "Beam element of type " << theElement->getClassType() << "not supported." << endln;
        return -1;
    }


    ElementIter& theElements = theDomain.getElements();
    int maxTag = 0;
    while ((theElement = theElements()) != 0)
    {
        if (!solidEleSetDefined)
        {
            const char* eleClassType = theElement->getClassType();
            if ((strcmp(eleClassType, "Brick")     == 0) ||
                (strcmp(eleClassType, "BbarBrick") == 0) ||
                (strcmp(eleClassType, "SSPbrick") == 0) ||
                (strcmp(eleClassType, "SSPbrickUP") == 0) ||
                (strcmp(eleClassType, "BbarBrick") == 0))
                solidEleTags.push_back(theElement->getTag());
        }
        if (theElement->getTag() > maxTag)
            maxTag = theElement->getTag();
    }
    int startTag = maxTag + 1;




    Node** solidNodesPtr;
    Node** beamNodePtr;
    Node*  toeNode;

    beamNodePtr = theDomain.getElement(beamTag)->getNodePtrs();
    toeNode     = theDomain.getNode(toeNodeTag);

    double Tx = toeNode->getCrds()(0);
    double Ty = toeNode->getCrds()(1);
    double Tz = toeNode->getCrds()(2);

    
    // use the geometric transformation object to generate beam surface points
    if (theTransf->initialize(beamNodePtr[0], beamNodePtr[1]))
    {
        opserr << "generateInterfacePoints: Error initializing coordinate transformation";
        return -1;
    }
    Vector loc_x(3), loc_y(3), loc_z(3);
    theTransf->getLocalAxes(loc_x, loc_y, loc_z);
    

    // create the points on the toe footprint
    int n      = (int)ceil(nP / nR);
    int numPts = n*nR*(nR + 1) / 2;
    nP         = n * nR;

    Vector cX(numPts), cY(numPts), cZ(numPts);
    Vector loc_rho(numPts), loc_theta(numPts);
    Vector pt_radius(numPts);
    
    int pts_Count = 0;
    for (int ii = 0; ii < nR; ii++)
    {
        int this_nP = n * (ii + 1);
        Vector t(this_nP);
        double this_radius = (double)(ii+1) / (double)nR * radius;
        for (int jj = 0; jj < this_nP; jj++)
        {
            t(jj)              = 0 + (double)jj / this_nP * 2.0 * PI;
            cX(pts_Count + jj) = 0.0;
            cY(pts_Count + jj) = this_radius * cos(t(jj));
            cZ(pts_Count + jj) = this_radius * sin(t(jj));

            loc_theta(pts_Count + jj) = t(jj);
            loc_rho(pts_Count + jj)   = -1.0;
            pt_radius(pts_Count + jj) = this_radius;
        }
        pts_Count += this_nP;
    }

    Vector cXr = cX*loc_x(0) + cY*loc_y(0) + cZ*loc_z(0) + Tx;
    Vector cYr = cX*loc_x(1) + cY*loc_y(1) + cZ*loc_z(1) + Ty;
    Vector cZr = cX*loc_x(2) + cY*loc_y(2) + cZ*loc_z(2) + Tz;

    double area = PI * radius * radius / numPts;

    if (debugFlag)
        for (int ii = 0; ii < numPts; ii++)
                opserr << "point " << ii + 1 << " : " << cXr(ii) << " " << cYr(ii) << " " << cZr(ii) << endln;


    Vector tempX(8), tempY(8), tempZ(8);
    int contactElemCount = 0, contactPointCount = 0;
    bool contactElemFlag = false;
    double xi, eta, zeta;
    bool inBounds = false;




    std::fstream crdsFile;
    if (writeCoords)
        crdsFile.open(crdsFN, std::fstream::app);

    for (int jj = 0; jj < numPts; jj++)
    {
        for (int ii = 0; ii < solidEleTags.size(); ii++)
        {
            theElement = theDomain.getElement(solidEleTags[ii]);
            solidNodesPtr = theElement->getNodePtrs();

            for (int kk = 0; kk < 8; kk++)
            {
                if (debugFlag)
                    opserr << "    Node " << kk + 1 << " coordinates = " << solidNodesPtr[kk]->getCrds();
                tempX(kk) = solidNodesPtr[kk]->getCrds()(0);
                tempY(kk) = solidNodesPtr[kk]->getCrds()(1);
                tempZ(kk) = solidNodesPtr[kk]->getCrds()(2);
            }
            invIsoMapping(tempX, tempY, tempZ, cXr(jj), cYr(jj), cZr(jj), xi, eta, zeta, inBounds);
            if (inBounds)
            {
                contactElemFlag = true;
                if (debugFlag)
                    opserr << "Beam tag : " << beamTag << ", Solid tag : " << solidEleTags[ii] << ", Real Coordinates = (" << cXr(jj) << "," << cYr(jj) << "," << cZr(jj) << "), Iso Coordinates = (" << xi << "," << eta << "," << zeta << ")" << endln;


                maxTag++;
                if (lagrangeTag > 1)
                {
                    solidEleTagsInContact.push_back(solidEleTags[ii]);
                    contactPt_xi.push_back(xi);
                    contactPt_eta.push_back(eta);
                    contactPt_zeta.push_back(zeta);
                    contactPt_rho.push_back(loc_rho(jj));
                    contactPt_theta.push_back(loc_theta(jj));
                    contactPt_radius.push_back(radius);
                }
                else
                {
                    // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange
                    if (lagrangeTag == 0)
                    {
                        // theElement = new EmbeddedBeamInterface(maxTag, beamTag, solidEleTags[ii], transfTag, loc_rho(jj), loc_theta(jj), xi, eta, zeta, radius, area);
                    }
                    else if (lagrangeTag == 1)
                    {
                        //theElement = new EmbeddedBeamInterfaceAL(maxTag, beamTag, solidEleTags[ii], transfTag, loc_rho(jj), loc_theta(jj), xi, eta, zeta, radius, area);
                    }

                    char buffer[40];
                    sprintf(buffer, "%10i", maxTag);
                    Tcl_AppendResult(interp, buffer, NULL);

                    if (theElement != 0)
                    {
                        if (theDomain.addElement(theElement) == false)
                        {
                            opserr << "WARNING could not add element with tag: " << theElement->getTag() << " and of type: "
                                << theElement->getClassType() << " to the Domain\n";
                            delete theElement;
                            return -1;
                        }
                    }

                }

                if (writeCoords)
                    crdsFile << maxTag << "\t" << cXr(jj) << "\t" << cYr(jj) << "\t" << cZr(jj) << "\n";


                eleTagsInContact_unique.insert(solidEleTags[ii]);

                contactPointCount++;
                break;
            }
        }
    }

    // 0 = local penalty, 1 = local Augmented Lagrangian, 2 = global Lagrange, 3 = global embedded lagrange
    if (lagrangeTag > 1)
    {
        if (lagrangeTag == 2)
        {
            opserr << "This element is not implemented yet." << endln;
            return -1;
        }
        else if (lagrangeTag == 3)
{
            //theElement = new EmbeddedBeamToeP(maxTag, beamTag, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, contactPt_radius, radius, area);
} 
       else if (lagrangeTag == 4)
{
            //theElement = new EmbeddedBeamToeP(maxTag, beamTag, solidEleTagsInContact, transfTag, contactPt_rho, contactPt_theta, contactPt_xi, contactPt_eta, contactPt_zeta, contactPt_radius, radius, area);
}
        else
        {
            opserr << "This element is not implemented yet." << endln;
            return -1;
        }
        char buffer[40];
        sprintf(buffer, "%10i", maxTag);
        Tcl_AppendResult(interp, buffer, NULL);

        if (theElement != 0)
        {
            if (theDomain.addElement(theElement) == false)
            {
                opserr << "WARNING could not add element with tag: " << theElement->getTag() << " and of type: "
                    << theElement->getClassType() << " to the Domain\n";
                delete theElement;
                return -1;
            }
        }
    }

    if (writeCoords)
        crdsFile.close();

    opserr << "Number of elements in contact: " << (int)eleTagsInContact_unique.size() << ", number of contact points: " << contactPointCount << endln;
    opserr << "EmbeddedBeamInterface Elements " << startTag << " to " << maxTag << " were created." << endln;


    if (debugFlag)
    {
        if (radiusSet)
            opserr << "Radius is set to " << radius << endln;
        else
            opserr << "No radius is set!!! Assuming 1.0" << endln;
        if (nPset)
            opserr << "nP is set to " << nP << endln;
        else
            opserr << "nP not set!!! Assuming " << nP << endln;
        if (nRset)
            opserr << "nL is set to " << nR << endln;
        else
            opserr << "nL not set!!! Assuming " << nR << endln;
        if (shapeSet)
            opserr << "Shape is set to " << shape << endln;
        else
            opserr << "Default shape (circle) is assumed." << endln;
        if (solidEleTags.size() > 0)
            opserr << "Here are the elements: \n     { ";
        for (int ii = 0; ii < solidEleTags.size(); ii++)
            opserr << solidEleTags[ii] << " ";
        if (solidEleTags.size() > 0)
            opserr << "}" << endln;
    }


    return 0;

}

int 
invIsoMapping(Vector nodesX, Vector nodesY, Vector nodesZ, double Px, double Py, double Pz, double & xi, double & eta, double & zeta, bool & inBounds)
{


	// This function returns the parameters of the isoparametric element using
	// trilinear shape functions.nodesX, nodesY and nodesZ contain the
	// actual coordinates of the element and p_x, p_y and p_z are the actual
	// coordinates of the point for which the parameters of the isoparametric
	// shape functions are needed.Obviously any assumptions regarding the
	// trilinear shape functions apply here as well.

	inBounds = true;

	double x1 = nodesX(0); double y1 = nodesY(0); double z1 = nodesZ(0);
	double x2 = nodesX(1); double y2 = nodesY(1); double z2 = nodesZ(1);
	double x3 = nodesX(2); double y3 = nodesY(2); double z3 = nodesZ(2);
	double x4 = nodesX(3); double y4 = nodesY(3); double z4 = nodesZ(3);
	double x5 = nodesX(4); double y5 = nodesY(4); double z5 = nodesZ(4);
	double x6 = nodesX(5); double y6 = nodesY(5); double z6 = nodesZ(5);
	double x7 = nodesX(6); double y7 = nodesY(6); double z7 = nodesZ(6);
	double x8 = nodesX(7); double y8 = nodesY(7); double z8 = nodesZ(7);

	double a = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 - 8 * Px;
	double a_xi = -x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8;
	double a_eta = -x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8;
	double a_zeta = -x1 - x2 - x3 - x4 + x5 + x6 + x7 + x8;
	double a_xe = x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8;
	double a_xz = x1 - x2 - x3 + x4 - x5 + x6 + x7 - x8;
	double a_ez = x1 + x2 - x3 - x4 - x5 - x6 + x7 + x8;
	double a_xez = -x1 + x2 - x3 + x4 + x5 - x6 + x7 - x8;

	double b = y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 - 8 * Py;
	double b_xi = -y1 + y2 + y3 - y4 - y5 + y6 + y7 - y8;
	double b_eta = -y1 - y2 + y3 + y4 - y5 - y6 + y7 + y8;
	double b_zeta = -y1 - y2 - y3 - y4 + y5 + y6 + y7 + y8;
	double b_xe = y1 - y2 + y3 - y4 + y5 - y6 + y7 - y8;
	double b_xz = y1 - y2 - y3 + y4 - y5 + y6 + y7 - y8;
	double b_ez = y1 + y2 - y3 - y4 - y5 - y6 + y7 + y8;
	double b_xez = -y1 + y2 - y3 + y4 + y5 - y6 + y7 - y8;

	double c = z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 - 8 * Pz;
	double c_xi = -z1 + z2 + z3 - z4 - z5 + z6 + z7 - z8;
	double c_eta = -z1 - z2 + z3 + z4 - z5 - z6 + z7 + z8;
	double c_zeta = -z1 - z2 - z3 - z4 + z5 + z6 + z7 + z8;
	double c_xe = z1 - z2 + z3 - z4 + z5 - z6 + z7 - z8;
	double c_xz = z1 - z2 - z3 + z4 - z5 + z6 + z7 - z8;
	double c_ez = z1 + z2 - z3 - z4 - z5 - z6 + z7 + z8;
	double c_xez = -z1 + z2 - z3 + z4 + z5 - z6 + z7 - z8;

	Matrix R(3, 3);
	R(0, 0) = a_xi; R(0, 1) = a_eta; R(0, 2) = a_zeta;
	R(1, 0) = b_xi; R(1, 1) = b_eta; R(1, 2) = b_zeta;
	R(2, 0) = c_xi; R(2, 1) = c_eta; R(2, 2) = c_zeta;

	Vector res(3);
	res(0) = -a; res(1) = -b; res(2) = -c;

	Vector A(3);
	R.Solve(res, A);


	xi = A(0);
	eta = A(1);
	zeta = A(2);

	for (int ii = 1; ii <= 20; ii++)
	{
		res(0) = a + a_xi*xi + a_eta*eta + a_zeta*zeta + a_xe*xi*eta + a_xz*xi*zeta + a_ez*eta*zeta + a_xez*xi*eta*zeta;
		res(1) = b + b_xi*xi + b_eta*eta + b_zeta*zeta + b_xe*xi*eta + b_xz*xi*zeta + b_ez*eta*zeta + b_xez*xi*eta*zeta;
		res(2) = c + c_xi*xi + c_eta*eta + c_zeta*zeta + c_xe*xi*eta + c_xz*xi*zeta + c_ez*eta*zeta + c_xez*xi*eta*zeta;

		if (res.Norm() < 1.0e-10)
		{
			if ((xi < -1) || (xi > 1) || (eta < -1) || (eta > 1) || (zeta < -1) || (zeta > 1))
				inBounds = false;

			// opserr << "Calculated in " << ii - 1 << " steps." << endln;
			return 0;
		}

		R(0, 0) = a_xi + a_xe * eta + a_xz * zeta + a_xez * eta * zeta;
		R(0, 1) = a_eta + a_xe * xi + a_ez * zeta + a_xez * xi  * zeta;
		R(0, 2) = a_zeta + a_xz * xi + a_ez * eta + a_xez * xi  * eta;
		R(1, 0) = b_xi + b_xe * eta + b_xz * zeta + b_xez * eta * zeta;
		R(1, 1) = b_eta + b_xe * xi + b_ez * zeta + b_xez * xi  * zeta;
		R(1, 2) = b_zeta + b_xz * xi + b_ez * eta + b_xez * xi  * eta;
		R(2, 0) = c_xi + c_xe * eta + c_xz * zeta + c_xez * eta * zeta;
		R(2, 1) = c_eta + c_xe * xi + c_ez * zeta + c_xez * xi  * zeta;
		R(2, 2) = c_zeta + c_xz * xi + c_ez * eta + c_xez * xi  * eta;

		R.Solve(res, A);

		xi -= A(0);
		eta -= A(1);
		zeta -= A(2);

	}

	// opserr << "Did not converge in 20 steps." << endln;
	return -1;
}

