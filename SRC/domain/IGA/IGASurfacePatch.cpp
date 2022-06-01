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


#include <IGASurfacePatch.h>

//OpenSees stuff
#include <elementAPI.h>
#include "Node.h"
#include "ElementResponse.h"


//IGA stuff
#include "NurbsDers.h"
#include "IGAKLShell.h"
#include "IGAKLShell_BendingStrip.h"
#include "gaussQuadrature.h"


//STD library
#include <float.h>
#include <math.h>
#include <vector>

#include <set>
using namespace std;
template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}
#include <algorithm>


// IGASurfacePatch(int tag, int P_, int Q_, const Vector& uKnot_, const Vector& vKnot_, const Matrix& controlPts_);


void PrintSyntax()
{
    opserr << "IGA Patch $tag $nodeStartTag $P $Q $noPtsX $noPtsY \n\
        -type KLShell \n\
        -nonLinearGeometry [0 or 1] \n\
        -planeStressMatTags [list of tags] \n\
        -gFact $gx $gy $gz \n\
        -theta [list of thetas] \n\
        -thickness [list of layer thicknesses] \n\
        -uKnot [list of uKnots] \n\
        -vKnot [list of vKnots] \n\
        -controlPts [list coordinates of control points, u-direction first] \n\
        " << endln;
}


void* OPS_IGASurfacePatch()
{

    // opserr << "OPS_IGASurfacePatch !!" << endln;
    // opserr << "OPS_IGASurfacePatch OPS_GetNumRemainingInputArgs() = " << OPS_GetNumRemainingInputArgs() << endln;

    int tag = 0;
    int P = 0;
    int Q = 0;

    int noPtsX = 0;
    int noPtsY = 0;

    int numdata = 1;

    ShellType shtype = ShellType::KLShell;

    int sectionTag = 0;
    int nodeStartTag = 0;
    int elementStartTag = 1;

    int nonLinearGeometry = 1;
    Vector gFact(3);



    opserr << "Creating IGA Patch:" << endln;
    if (OPS_GetIntInput(&numdata, &tag) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "tag = " << tag << endln;
    if (OPS_GetIntInput(&numdata, &nodeStartTag) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "nodeStartTag = " << nodeStartTag << endln;
    if (OPS_GetIntInput(&numdata, &P) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "P = " << P << endln;
    if (OPS_GetIntInput(&numdata, &Q) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "Q = " << Q << endln;
    if (OPS_GetIntInput(&numdata, &noPtsX) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "noPtsX = " << noPtsX << endln;
    if (OPS_GetIntInput(&numdata, &noPtsY) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "noPtsY = " << noPtsY << endln;


    
    
    // opserr << "OPS_IGASurfacePatch OPS_GetNumRemainingInputArgs() = " << OPS_GetNumRemainingInputArgs() << endln;




    // get other inputs
    std::vector<double> theta_stdVector, thickness_stdVector , uKnot_stdVector, vKnot_stdVector, controlPts_stdVector;
    std::vector<int> matTags_stdVector;
    // int loc = 5; // Vector (knotVector) start position
    // int loc = 7;
    int loc = 8;
    while (OPS_GetNumRemainingInputArgs() > 0) {

        // opserr << "loc = " << loc << endln;

        // next arg
        const char* arg = OPS_GetString();
        loc++;

        // opserr << "arg = " << arg << endln;

        // check arg
        if (strcmp(arg, "-type") == 0) { //matTags
            // const char* type = OPS_GetString();
            // loc++;
            const char* typestring = OPS_GetString();
            opserr << "type = " << typestring << endln;
            if (strcmp(typestring, "KLShell") == 0)
            {
                shtype = ShellType::KLShell;
            }
            else if (strcmp(typestring, "KLShell_BendingStrip") == 0)
            {
                shtype = ShellType::KLShell_BendingStrip;
            }
            else
            {
                opserr << "IGASurfacePatch - Unknown shell of type " << typestring << endln;
            }
            loc++;
        }
        else if (strcmp(arg, "-planeStressMatTags") == 0) { //matTags
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int val;
                if (OPS_GetIntInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                matTags_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-theta") == 0) { //theta
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                theta_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-thickness") == 0) { //thickness
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                thickness_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-uKnot") == 0) { //uKnot
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                uKnot_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-vKnot") == 0) { //vKnot
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                vKnot_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-controlPts") == 0) { //controlPts
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                controlPts_stdVector.push_back(val);
                loc++;
            }
        }
        // else if (strcmp(arg, "-type") == 0) {
        //     const char* typestring = OPS_GetString();
        //     if (strcmp(typestring, "KLShell") == 0)
        //     {
        //         shtype = ShellType::KLShell;
        //     }
        //     else if (strcmp(typestring, "KLShell_BendingStrip") == 0)
        //     {
        //         shtype = ShellType::KLShell_BendingStrip;
        //     }
        //     else
        //     {
        //         opserr << "IGASurfacePatch - Unknown shell of type " << typestring << endln;
        //     }
        //     loc++;
        // }
        else if (strcmp(arg, "-sectionTag") == 0) {
            numdata = 1;

            if (OPS_GetIntInput(&numdata, &sectionTag) < 0) return 0;
            loc++;
        }
        else if (strcmp(arg, "-nodeStartTag") == 0) {
            numdata = 1;

            if (OPS_GetIntInput(&numdata, &nodeStartTag) < 0) return 0;
            loc++;
        }
        else if (strcmp(arg, "-elementStartTag") == 0) {
            numdata = 1;

            if (OPS_GetIntInput(&numdata, &elementStartTag) < 0) return 0;
            loc++;
        }
        else if (strcmp(arg, "-nonLinearGeometry") == 0) {
            numdata = 1;

            if (OPS_GetIntInput(&numdata, &nonLinearGeometry) < 0) return 0;
            loc++;
        }
        else if (strcmp(arg, "-gFact") == 0) {
            int i = 0;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                gFact(i) = val;
                i++;
                loc++;
            }
        }


    }

    ID matTags(&matTags_stdVector[0], (int)(matTags_stdVector.size()));
    Vector theta(&theta_stdVector[0], (int)(theta_stdVector.size()));
    Vector thickness(&thickness_stdVector[0], (int)(thickness_stdVector.size()));
    Vector uKnot(&uKnot_stdVector[0], (int)uKnot_stdVector.size());
    Vector vKnot(&vKnot_stdVector[0], (int)vKnot_stdVector.size());


    int controlPts_size = controlPts_stdVector.size();
    int M = controlPts_size / 4;
    int N = 4;

    // Matrix controlPts(&controlPts_stdVector[0], M, N);
    Matrix controlPts(&controlPts_stdVector[0], N, M);
    // opserr << "controlPts = " << controlPts ;
    IGASurfacePatch* patch = new IGASurfacePatch(tag, nodeStartTag, P, Q, noPtsX, noPtsY, nonLinearGeometry, gFact, matTags, theta, thickness, uKnot, vKnot, controlPts, shtype); //
    // opserr << "tag = " << tag << endln;
    // opserr << "\n\nprint command from OPS_IGASurfacePatch" << endln;
    // patch->Print(opserr);

    Domain* theDomain = OPS_GetDomain();
    theDomain->addElement(patch);

    // return 0;
    return patch;
}

// IGASurfacePatch::IGASurfacePatch(int tag, int P_, int Q_)
IGASurfacePatch::IGASurfacePatch(int tag, int P_, const Vector & uKnot_, const Matrix & controlPts_)
    :
    Subdomain(tag),
    P(P_),
    Q(0),
    R(0),
    uKnot(uKnot_),
    vKnot(0),
    wKnot(0),
    controlPts(controlPts_)
{
    opserr << "IGASurfacePatch::IGASurfacePatch() - 1D constructor" << endln;
    Print(opserr);

}


// IGASurfacePatch::IGASurfacePatch(int tag, int P_, int Q_, const Vector & uKnot_, const Vector & vKnot_, const Matrix & controlPts_)
IGASurfacePatch::IGASurfacePatch(int tag, int nodeStartTag_, int P_, int Q_, int noPtsX_, int noPtsY_, int nonLinearGeometry_, const Vector gFact_, const ID& matTags_, const Vector& theta_, const Vector& thickness_, const Vector& uKnot_, const Vector& vKnot_, const Matrix& controlPts_, ShellType shtype_)
    :
    Subdomain(tag),
    nodeStartTag(nodeStartTag_),
    P(P_),
    Q(Q_),
    R(0),
    noPtsX(noPtsX_), // testing
    noPtsY(noPtsY_), //testing
    gFact(gFact_),
    nonLinearGeometry(nonLinearGeometry_),
    noPts(noPtsX * noPtsY),
    matTags(matTags_),
    theta(theta_),
    thickness(thickness_),
    Zk(0),
    uKnot(uKnot_),
    vKnot(vKnot_),
    wKnot(0),
    noFuncs(0),
    controlPts(controlPts_),
    weights(0),
    shtype(shtype_),
    index(0),
    elRangeU(0),
    elRangeV(0),
    element(0),
    elConnU(0),
    elConnV(0),
    quadPoint(0),
    quadWeight(0)
{
    // opserr << "IGASurfacePatch::IGASurfacePatch() - 2D constructor" << endln;
    Print(opserr);

    generateIGA2DMesh(noElems, noElemsU, noElemsV);
    noFuncs = (P + 1) * (Q + 1);


    // Creating z vector with mid center coordinates of laminate plies

    // Vector Zk(thickness.Size()); // Vector with mid center coordinates of laminate plies measured from the top
    Zk = new Vector(thickness.Size()); // Pointer to Zk vector
    (*Zk)(0) = thickness(0) / 2; //First ply coordinate
    double rho = (OPS_getNDMaterial(getMatTag(0)))->getRho(); // Density of each laminate
    double Zbar = (*Zk)(0) * thickness(0) * rho; // Mid center z coordinate measured from the top, calculated on the next for loop
    double sumThicknessRho = thickness(0) * rho; // thickness*rho to get center of mass

    for (int i = 1; i < thickness.Size(); ++i)
    {
        rho = (OPS_getNDMaterial(getMatTag(i)))->getRho(); // Density of material
        (*Zk)(i) = (*Zk)(i - 1) + thickness(i - 1) / 2 + thickness(i) / 2;
        Zbar += (*Zk)(i) * thickness(i) * rho;
        sumThicknessRho += thickness(i) * rho;
    }

    Zbar /= sumThicknessRho;
    for (int i = 0; i < (*Zk).Size(); ++i)
    {
        (*Zk)(i) = Zbar - (*Zk)(i); //Distance from midcenter
    }

    // End creating Zk vector




    // // Defining Gauss Quadrature points
    // int quadorder = 4;
    // quadPoint = new Matrix(quadorder * quadorder, 2);
    // quadWeight = new Vector(quadorder * quadorder);
    // quadrature(quadorder, quadPoint, quadWeight);

    // int idU;
    // int idV;
    // Vector xiE(2);
    // Vector etaE(2);

    // double xi;
    // double eta;

    // double ptU;
    // double ptV;



    // opserr << endln << "Looping through the elements and showing the number, node span and gauss point coordinates" << endln << endln;
    // // Y aqui iteron sobre los elementos
    // for (int e = 0; e < noElems; ++e)
    // {
    //     opserr << "Element n: " << e << endln;

    //     idU = (*index)(e, 0);
    //     idV = (*index)(e, 1);

    //     xiE(0) = (*elRangeU)(idU, 0);
    //     xiE(1) = (*elRangeU)(idU, 1);

    //     etaE(0) = (*elRangeV)(idV, 0);
    //     etaE(1) = (*elRangeV)(idV, 1);

    //     // opserr << "Parametric global coordinates span of element in each direction:" << endln;

    //     // opserr << "xiE = " << xiE << endln;
    //     // opserr << "etaE = " << etaE << endln;



    //     opserr << "Element spans these nodes :";
    //     for (int i = 0; i < (*element).noCols(); ++i)
    //     {
    //         opserr << (*element)(e, i) << ", ";
    //     }
    //     opserr << endln;




    //     for (int gp = 0; gp < (*quadWeight).Size(); ++gp)
    //     {
    //         ptU = (*quadPoint)(gp, 0);
    //         ptV = (*quadPoint)(gp, 1);
    //         opserr << "Gauss point n: " << gp << endln;
    //         opserr << "ptU, ptV : " << ptU << ", " << ptV << endln;

    //         //Global parametric coordinates of gauss point [0-1], dont know if needed
    //         xi = parent2ParametricSpace(xiE, ptU);
    //         eta = parent2ParametricSpace(etaE, ptV);
    //         opserr << "xi, eta: " << xi << ", " << eta << endln;


    //     }
    //     opserr << endln;


    // }








    // xi = 0.46528357789851327; //Just for testing the derivatives, these will come from
    // eta = 0.46528357789851327; // the mapping from the parent to parametric space of the gauss points

    // noFuncs = (P + 1) * (Q + 1);     //Number of basis functions for the elements, these should be given to the element
    // Vector R(noFuncs);               //Value of each basis function R
    // Vector dRdxi(noFuncs);           //First derivative wrt to Xi
    // Vector dRdeta(noFuncs);          //First derivative wrt to Eta
    // Vector dR2dxi(noFuncs);          //Second derivative wrt to Xi
    // Vector dR2deta(noFuncs);         //Second derivative wrt to Eta
    // Vector dR2dxideta(noFuncs);      //Mixed second derivative wrt to Xi,Eta

    // Second derivatives of the basis functions
    // Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

    // generateIGA2DMesh();

    // opserr << "Finished Patch" << endln;

}

// IGASurfacePatch::IGASurfacePatch(int tag, int P_, int Q_)
IGASurfacePatch::IGASurfacePatch(int tag, int P_, int Q_, int R_, const Vector & uKnot_, const Vector & vKnot_, const Vector & wKnot_, const Matrix & controlPts_)
    :
    Subdomain(tag),
    P(P_),
    Q(Q_),
    R(0),
    uKnot(uKnot_),
    vKnot(vKnot_),
    wKnot(0),
    controlPts(controlPts_)
{

    opserr << "IGASurfacePatch::IGASurfacePatch() - 3D constructor" << endln;
    Print(opserr);



}

void
IGASurfacePatch::clearAll(void)
{
    opserr << "IGASurfacePatch::clearAll\n";
}

void IGASurfacePatch::setDomain(Domain *theDomain)
{
    opsDomain = theDomain;

    opserr << "IGASurfacePatch::setDomain ->" <<  endln;

    opserr << "nodeStartTag = " << nodeStartTag << endln;

    double xCoord;
    double yCoord;
    double zCoord;

    Vector weights(noPts);

    // Iterando sobre nodos, imprimo numero nodo global, numero local, coordenadas
    for (int i = 0; i < noPts; ++i)
    {
        // opserr << "controlPoint n: " << i << endln;
        xCoord = controlPts(0, i);
        yCoord = controlPts(1, i);
        zCoord = controlPts(2, i);

        weights(i) = controlPts(3, i);


        // opserr << "xCoord = " << xCoord << endln;
        // opserr << "yCoord = " << yCoord << endln;
        // opserr << "zCoord = " << zCoord << endln << endln;

        int tag = i + nodeStartTag; // Hacemos partir los nodos desde nodeStartTag, entregado desde python (para multipatches)
        int ndof = 3;
        Node * node = 0;

        if ((node = new Node( tag,  ndof,  xCoord,  yCoord,  zCoord)) != 0)
        {
            opsDomain->addNode(node);

            // Pipe: talvez aqui esta el error de transformation, comente esto 
            // this->addNode(node);
        }
        else
        {
            opserr << "IGASurfacePatch::setDomain - Out of memory creating new node. " << endln;
            return ;
        }
    }




    // Defining Gauss Quadrature points
    int quadorder;

    quadorder = (P + 1) * (Q + 1);
    // quadorder = (P+1)*(P+1);
    // quadorder = 6;

    // quadPoint = new Matrix(quadorder, 2);
    // quadWeight = new Vector(quadorder);

    // gaussQuad2dNurbs(P+1,Q+1, quadPoint, quadWeight);
    // gaussQuad2dNurbs(P+1,P+1, quadPoint, quadWeight);

    // opserr << "*quadWeight = " << *quadWeight << endln;
    // opserr << "*quadPoint = " << *quadPoint << endln;

    int idU;
    int idV;
    Vector xiE(2);
    Vector etaE(2);

    // double xi;
    // double eta;

    // double ptU;
    // double ptV;



    // opserr << endln << "Looping through the elements and showing the number, node span and gauss point coordinates" << endln << endln;
    // Y aqui itero sobre los elementos
    for (int e = 0; e < noElems; ++e)
    {
        // opserr << "Element n: " << e << endln;

        idU = (*index)(e, 0);
        idV = (*index)(e, 1);

        xiE(0) = (*elRangeU)(idU, 0);
        xiE(1) = (*elRangeU)(idU, 1);

        etaE(0) = (*elRangeV)(idV, 0);
        etaE(1) = (*elRangeV)(idV, 1);

        int tag = this->getTag() + e + 1;  // Hacemos partir los elementos desde el numero de patch en adelante
        int nnodes = (*element).noCols();
        ID nodes(nnodes);

        // opserr << "Element spans these nodes :";
        for (int i = 0; i < (*element).noCols(); ++i)
        {

            nodes(i) = (*element)(e, i) + nodeStartTag;    /// Ojo !!! con el nodeStartTag
            // opserr << "nodes(i) = " << nodes(i) << endln;
            // nodes(i) = (*element)(e, i);    /// Ojo !!! + 1
        }
        opserr << endln;


        Element* ele = 0 ;

        if (shtype == KLShell)
        {
            // Normal KLShell
            if ((ele = new IGAKLShell(tag, this, nodes, quadorder, xiE, etaE, matTags)) != 0)
            {
                opsDomain->addElement(ele);
            } else
            {
                opserr << "IGASurfacePatch::setDomain - Out of memory creating new element. " << endln;
                return ;
            }
        }
        else
        {
            // Bending Strip KLShell
            if ((ele = new IGAKLShell_BendingStrip(tag, this, nodes, quadorder, xiE, etaE, matTags)) != 0)
            {
                opsDomain->addElement(ele);
            } else
            {
                opserr << "IGASurfacePatch::setDomain - Out of memory creating new element. " << endln;
                return ;
            }
        }

        // for (int gp = 0; gp < (*quadWeight).Size(); ++gp)
        // {
        //     ptU = (*quadPoint)(gp, 0);
        //     ptV = (*quadPoint)(gp, 1);
        //     // opserr << "Gauss point n: " << gp << endln;
        //     // opserr << "ptU, ptV : " << ptU << ", " << ptV << endln;

        //     //Global parametric coordinates of gauss point [0-1], dont know if needed
        //     xi = parent2ParametricSpace(xiE, ptU);
        //     eta = parent2ParametricSpace(etaE, ptV);
        //     // opserr << "xi, eta: " << xi << ", " << eta << endln;


        // }


    }






}



bool IGASurfacePatch::buildConnectivity(int p, const Vector & knotVec, int nE, Matrix * elRange, Matrix * elConn)
{
    bool result = false;


    Matrix elKnotIndices(nE, 2);
    elKnotIndices.Zero(); //Se inicializan en 0?

    int element = 0;
    double previousKnotVal = 0;
    double currentKnotVal = 0;

    // opserr << "Starting first loop" << endln;
    for (int i = 0; i < knotVec.Size(); ++i)
    {
        currentKnotVal = knotVec(i);
        if (knotVec(i) != previousKnotVal)
        {
            (*elRange)(element, 0) = previousKnotVal;
            (*elRange)(element, 1) = currentKnotVal; // Hay alguna manera de settear una fila?
            elKnotIndices(element, 0) = i - 1;
            elKnotIndices(element, 1) = i;
            element += 1;
        }
        previousKnotVal = currentKnotVal;
    }
    // opserr << "Ending first loop" << endln;

    // opserr << "elKnotIndices = " << elKnotIndices << endln;


    int numRepeatedKnots = 0;

    // opserr << "Starting second loop" << endln;
    for (int e = 0; e < nE; ++e)
    {
        // opserr << "botInd = " << (elKnotIndices(e, 0) - p) << endln;
        // opserr << "topInd = " << (elKnotIndices(e, 0)) << endln;

        int botInd = elKnotIndices(e, 0) - p;
        int topInd = elKnotIndices(e, 0);

        // auto indices = arange<int>(elKnotIndices(e, 0) - p + 1, elKnotIndices(e, 0) + 1);

        Vector indices(abs(topInd - botInd) + 1);
        indices.Zero();
        for (int i = 0; i < indices.Size(); ++i)
        {
            indices(i) = botInd + i;
        }

        // opserr << "indices = " << indices << endln;

        // opserr << "Debug conec" << endln;
        // opserr << "indices[0] = " << indices[0] << endln;
        // opserr << "indices[1] = " << indices[1] << endln;
        Vector previousKnotVals(2);
        previousKnotVals(0) = knotVec(indices[0]);
        previousKnotVals(1) = knotVec(indices[1]);
        // opserr << "Debug conec end" << endln;

        Vector ones(p);
        ones += 1;
        Vector currentKnotVals = ones * knotVec(elKnotIndices(e, 0));

        int nonzero = 0;
        for (int i = 0; i < previousKnotVals.Size(); ++i)
        {
            if (previousKnotVals(i) != 0)
            {
                nonzero += 1;
            }
            if (previousKnotVals(i) != currentKnotVals(i))
            {
                nonzero += 1;
            }
        }

        if (previousKnotVals == currentKnotVals && nonzero > 1)
        {
            numRepeatedKnots += 1;
        }

        auto elConn_arange = arange<int>(elKnotIndices(e, 0) - p, elKnotIndices(e, 0) + 1);
        for (unsigned int i = 0; i < elConn_arange.size(); ++i)
        {
            (*elConn)(e, i) = elConn_arange[i];
        }
    }
    // opserr << "Ending second loop" << endln;



    //* Implementar *//
    // opserr << "Finished buildConnectivity" << endln;

    return result;
}

// Destructor
IGASurfacePatch::~IGASurfacePatch()
{
    // opserr << "IGASurfacePatch::~IGASurfacePatch" << endln;
    // opserr << "IGASurfacePatch::IGASurfacePatch" << endln;

    if (index!= 0) {
        delete index;
        index = 0;
    };
    if (elRangeU!= 0) {
        delete elRangeU;
        elRangeU = 0;
    };
    if (elRangeV!= 0) {
        delete elRangeV;
        elRangeV = 0;
    };
    if (element!= 0) {
        delete element;
        element = 0;
    };
    if (elConnU!= 0) {
        delete elConnU;
        elConnU = 0;
    };
    if (elConnV!= 0) {
        delete elConnV;
        elConnV = 0;
    };
    if (quadPoint!= 0) {
        delete quadPoint;
        quadPoint = 0;
    };
    if (quadWeight!= 0) {
        delete quadWeight;
        quadWeight = 0;
    };
    if (Zk!= 0) {
        delete Zk;
        Zk = 0;
    };
}


bool IGASurfacePatch::generateIGA2DMesh(int & noElems, int & noElemsU, int & noElemsV)//, Matrix & index, Matrix & elRangeU, Matrix & elRangeV, Matrix & elConnU, Matrix & elConnV)
{

    // opserr << "generateIGA2DMesh" << endln;
    bool result = false;

    // Using std::unique
    int uKnot_size = uKnot.Size();
    int vKnot_size = vKnot.Size();

    // Get unique elements of knot vectors to determine number of elements
    std::set<double> uniqueUKnots_set( &uKnot(0), &uKnot(0) + uKnot_size );
    std::set<double> uniqueVKnots_set( &vKnot(0), &vKnot(0) + vKnot_size );

    noElemsU = uniqueUKnots_set.size() - 1; //Now these are given by reference
    noElemsV = uniqueVKnots_set.size() - 1;

    opserr << "noElemsU = " << noElemsU << endln;
    opserr << "uKnot = " << uKnot << endln;

    opserr << "noElemsV = " << noElemsV << endln;
    opserr << "vKnot = " << vKnot << endln;




    // combine info from two directions to build the elements
    // element is numbered as follows
    //  5 | 6 | 7 | 8
    // ---------------
    //  1 | 2 | 3 | 4
    // for a 4x2 mesh

    // opserr << "Generating chan" << endln;

    Matrix chan(noPtsY, noPtsX);
    // chan.Zero();
    int count = 0;
    for (int i = 0; i < noPtsY; ++i)
    {
        for (int j = 0; j < noPtsX; ++j)
        {
            chan(i, j) = count;
            count++;
        }
    }



    // determine our element ranges and the corresponding
    // knot indices along each direction

    elRangeU = new Matrix(noElemsU, 2);
    elRangeV = new Matrix(noElemsV, 2);


    elConnU = new Matrix(noElemsU, P + 1);
    elConnV = new Matrix(noElemsV, Q + 1);

    // opserr << "Start Debugging " << endln;
    buildConnectivity(P, uKnot, noElemsU, elRangeU, elConnU);
    // opserr << "End Debugging " << endln;
    buildConnectivity(Q, vKnot, noElemsV, elRangeV, elConnV);




    noElems = noElemsU * noElemsV;
    element = new Matrix(noElems, (P + 1) * (Q + 1));

    int e = 0;
    for (int v = 0; v < noElemsV; ++v)
    {
        Vector vConn(Q + 1);
        for (int i = 0; i < vConn.Size(); ++i)
        {
            vConn(i) = (*elConnV)(v, i);
        }
        for (int u = 0; u < noElemsU; ++u)
        {
            int c = 0;

            Vector uConn(P + 1);
            for (int i = 0; i < uConn.Size(); ++i)
            {
                uConn(i) = (*elConnU)(u, i);
            }

            for (int i = 0; i < vConn.Size(); ++i)
            {
                for (int j = 0; j < uConn.Size(); ++j)
                {
                    (*element)(e, c) = chan((int)vConn(i), (int)uConn(j));
                    c++;
                }
            }
            e++;
        }
    }

    index = new Matrix(noElems, 2);


    count = 0;

    for (int j = 0; j < elRangeV->noRows(); ++j)
    {
        for (int i = 0; i < elRangeU->noRows(); ++i)
        {
            (*index)(count, 0) = i;
            (*index)(count, 1) = j;
            count++;
        }
    }







    //* Implementar *//
    result = true;

    return result;
}

int IGASurfacePatch::Nurbs2DBasis2ndDers(double xi, double eta, Vector & R, Vector & dRdxi, Vector & dRdeta, Vector & dR2dxi, Vector & dR2deta, Vector & dR2dxideta)
{
    bool result = false;


    // opserr << "weights = " << weights << endln;

    /* Return the 2D NURBS basis functions and first/second derivatives
       All non-zero basis functions and derivatives at point [xi,eta] are computed.

      We expect the function to be called as
      [R dRdxi dRdeta dR2dxi dR2deta dR2dxideta] = NURBS2DBasis2ndDers(...
                                xi, p, q, knotU,knotV, weights)

        xi           = point, [xi eta], where we want to interpolate
        knotU, knotV = knot vectors
        weights      = vector of weights

    Written originally for Matlab by Vinh Phu Nguyen, nvinhphu@gmail.com
    */

    //* Implementar *//

    // opserr << "knotU = " << this->uKnot << endln;



    int numKnotU = uKnot.Size();
    int numKnotV = vKnot.Size();

    int nU = numKnotU - 1 - P  - 1;
    int nV = numKnotV - 1 - Q  - 1;

    double tol = 100 * DBL_EPSILON;

    if (fabs(xi - uKnot[numKnotU - 1]) < tol)
        xi = uKnot[numKnotU - 1] - tol;

    if (fabs(eta - vKnot[numKnotV - 1]) < tol)
        eta = vKnot[numKnotV - 1] - tol;



    /* and evaluate the non-zero univariate B-spline basis functions
     * and first derivatives
     */

    // double *N      = (double *)malloc(sizeof(double) * (P + 1));
    // double *M      = (double *)malloc(sizeof(double) * (Q + 1));

    static Vector N(P+1);
    static Vector M(Q+1);
    N.resize(P+1);
    M.resize(Q+1);

    // double **dersN = init2DArray(nU + 1, P + 1);
    // double **dersM = init2DArray(nV + 1, Q + 1);

    static Matrix dersN(nU + 1, P + 1);
    static Matrix dersM(nV + 1, Q + 1);
    dersN.resize(nU + 1, P + 1);
    dersM.resize(nV + 1, Q + 1);



    int spanU      = FindSpan(nU, P, xi, uKnot);
    int spanV      = FindSpan(nV, Q, eta, vKnot);

    BasisFuns     (spanU, xi, P, uKnot, N);
    BasisFuns     (spanV, eta, Q, vKnot, M);

    dersBasisFuns (spanU, xi, P, nU, uKnot, dersN);
    dersBasisFuns (spanV, eta, Q, nV, vKnot, dersM);



    /* and create NURBS approximation */


    int i, j, k, c;

    int uind = spanU - P;
    int vind;

    double w      = 0.0; /* w = N_I w_I*/
    double dwdxi  = 0.0; /* first derivative of w w.r.t xi*/
    double d2wdxi = 0.0; /* second derivative of w w.r.t xi*/
    double dwdet  = 0.0; /* first derivative of w w.r.t eta*/
    double d2wdet = 0.0; /* second derivative of w w.r.t eta*/
    double d2wdxe = 0.0; /* second derivative of w w.r.t xi-eta*/
    double wi;



    for (j = 0; j <= Q; j++)
    {
        vind = spanV - Q + j;

        for (i = 0; i <= P; i++)
        {
            c   = uind + i + vind * (nU + 1);
            wi  = controlPts(3, c); //Weight of control point
            // wi  = weights(c); //Weight of control point

            w      += N[i]        * M[j] * wi;
            dwdxi  += dersN(1,i) * M[j] * wi;
            if (P >= 2)
            {
                d2wdxi += dersN(2,i) * M[j] * wi;
            }
            dwdet  += dersM(1,j) * N[i] * wi;
            if (Q >= 2)
            {
                d2wdet += dersM(2,j) * N[i] * wi;
            }
            d2wdxe += dersN(1,i) * dersM(1,j) * wi;
        }

    }


    uind = spanU - P;
    k    = 0;

    double fac;
    double NiMj, w3Inv, w2Inv;

    for (j = 0; j <= Q; j++)
    {
        vind = spanV - Q + j;

        for (i = 0; i <= P; i++)
        {
            c        = uind + i + vind * (nU + 1);
            wi = controlPts(3, c); //Weight of control point
            // wi = weights(c); //Weight of control point

            NiMj     = N[i] * M[j];
            w3Inv    = 1 / w / w / w;
            w2Inv    = w3Inv * w;
            fac      = wi * w2Inv;

            R(k)     = NiMj * fac * w; // NiMj * wi / w

            dRdxi(k) = (dersN(1,i) * M[j] * w - NiMj * dwdxi) * fac;
            dRdeta(k) = (dersM(1,j) * N[i] * w - NiMj * dwdet) * fac;

            if (P >= 2)
            {
                // dR2dxi(k) = wi * (dersN[2][i] * M[j] / w - 2 * dersN[1][i] * M[j] * dwdxi / w / w - N[i] * M[j] * d2wdxi / w / w + 2 * N[i] * M[j] * dwdxi * dwdxi / w / w / w);
                dR2dxi(k) = wi * (dersN(2,i) * M[j] / w - 2 * dersN(1,i) * M[j] * dwdxi * w2Inv - NiMj * d2wdxi * w2Inv + 2 * NiMj * dwdxi * dwdxi * w3Inv);
            }
            if (Q >= 2)
            {
                // dR2deta(k) = wi * (dersM[2][j] * N[i] / w - 2 * dersM[1][j] * N[i] * dwdet / w / w - N[i] * M[j] * d2wdet / w / w + 2 * N[i] * M[j] * dwdet * dwdet / w / w / w);
                dR2deta(k) = wi * (dersM(2,j) * N[i] / w - 2 * dersM(1,j) * N[i] * dwdet * w2Inv - NiMj * d2wdet * w2Inv + 2 * NiMj * dwdet * dwdet * w3Inv);
            }
            // dR2dxideta(k) = wi * (dersN[1][i] * dersM[1][j] / w - dersN[1][i] * M[j] * dwdet / w / w - N[i] * dersM[1][j] * dwdxi / w / w - N[i] * M[j] * d2wdxe / w / w + 2 * N[i] * M[j] * dwdxi * dwdet / w / w / w);
            dR2dxideta(k) = wi * (dersN(1,i) * dersM(1,j) / w - dersN(1,i) * M[j] * dwdet * w2Inv - N[i] * dersM(1,j) * dwdxi * w2Inv - NiMj * d2wdxe * w2Inv + 2 * NiMj * dwdxi * dwdet * w3Inv);

            k += 1;
        }
    }

    // free(N);
    // free(M);
    // free2Darray(dersN, (nU + 1));
    // free2Darray(dersM, (nV + 1));



    return result;
}

double IGASurfacePatch::parent2ParametricSpace(Vector range, double xibar)
{
    double xi = 0.5 * ((range(1) - range(0)) * xibar + range(1) + range(0));
    return xi;
}

int IGASurfacePatch::getNoFuncs()
{
    int noFuncs = (P + 1) * (Q + 1);
    return noFuncs;
}

ID IGASurfacePatch::getOrders()
{
    ID PQ(2);
    PQ(0) = P;
    PQ(1) = Q;

    return PQ;
}

int IGASurfacePatch::getNLayers()
{
    return thickness.Size();
}

int IGASurfacePatch::getMatTag(int iLayer)
{
    return matTags(iLayer);
}

double IGASurfacePatch::getThickness(int iLayer)
{
    return thickness(iLayer);
}

double IGASurfacePatch::getAngle(int iLayer)
{
    return theta(iLayer);
}

double IGASurfacePatch::getZk(int iLayer)
{
    return (*Zk)(iLayer);
}

bool IGASurfacePatch::getAnalysisType()
{
    if (nonLinearGeometry == 1)
    {
        return true;
    }
    return false;
}

Vector IGASurfacePatch::getGravityFactors()
{
    return gFact;
}







void IGASurfacePatch::Print(OPS_Stream & s, int flag)
{

    opserr << " P = " << P << endln;
    opserr << " Q = " << Q << endln;
    opserr << " R = " << R << endln << endln;
    opserr << " noPts = " << noPts << endln;
    opserr << " noPtsX = " << noPtsX << endln;
    opserr << " noPtsY = " << noPtsY << endln << endln;
    opserr << " uKnot = " << uKnot ;
    opserr << " vKnot = " << vKnot ;
    opserr << " wKnot = " << wKnot << endln;;
    opserr << " controlPts = " << controlPts << endln;
    opserr << "Finished printing" << endln;
    return;
}


int IGASurfacePatch::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    // opserr << "IGASurfacePatch::addLoad - Element tag = " << this->getTag() << endln;
    Element* theElement = 0;
    for (int ele = 0; ele < noElems; ++ele)
    {
        int tag = this->getTag() + 1 + ele;
        theElement = opsDomain->getElement(tag);
        if (theElement != 0)
        {
            theElement->addLoad(theLoad, loadFactor);
        } else
        {
            opserr << "IGASurfacePatch::addLoad - Element tag " << tag << " not found in main domain" << endln;
            return -1;
        }
        /* code */
    }
    return 0;
}

#define PARAM_ADVANCESTATE 1
#define PARAM_RESETSTRESSES 2

int IGASurfacePatch :: setParameter(const char **argv, int argc, Parameter &param)
{

    // opserr << "IGASurfacePatch :: setParameter called " << endln;
    // opserr << 'argv[0] = ' << argv[0] << endln; 
    // opserr << 'argv[1] = ' << argv[1] << endln; 
    if (argc < 1) // was 2
        return -1;

    int res = -1;

    Element* theElement = 0;
    for (int ele = 0; ele < noElems; ++ele)
    {
        int matRes = res;

        int tag = this->getTag() + 1 + ele;
        theElement = opsDomain->getElement(tag);
        if (theElement != 0)
        {
            // opserr << "Debugging" << endln;
            // opserr << "Element calling setParameter " << endln;
            matRes =  theElement->setParameter(argv, argc, param);
            if (matRes != -1) res = matRes;
        } else {
            opserr << "IGASurfacePatch :: setParameter - Element tag " << tag << " not found in main domain" << endln;
            res = -1;
        }
        /* code */
    }
    return res;
}

int IGASurfacePatch :: updateParameter(int parameterID, Information &info)
{
    int res = -1;
    int matRes = res;

    // opserr << "IGASurfacePatch :: updateParameter called"<< endln;

    if (parameterID == res)
    {
        return -1;
    } else {
        Element* theElement = 0;
        for (int ele = 0; ele < noElems; ++ele)
        {
            int matRes = res;

            int tag = this->getTag() + 1 + ele;
            theElement = opsDomain->getElement(tag);
            if (theElement != 0)
            {
                // opserr << " IGASurfacePatch :: updateParameter called" << endln;
                matRes =  theElement->updateParameter(parameterID, info);
                if (matRes != -1) res = matRes;
            } else {
                opserr << "IGASurfacePatch :: updateParameter - Element tag " << tag << " not found in main domain" << endln;
                res = -1;
            }
            /* code */
        }
    }
    return res;
}

#define RESPONSETYPE_IGASurfacePatch_NODELUMPED 1

// For MPCO recorder
// "IGAOrder" u-v orders [Vector 1-line 2-surface....] INT32
// "IGAKnot1P" xi coordinates gauss [Vectror od doubles) all
// "IGAKnot2P" eta coord gauss [Vectror od doubles) volumes and surface
// "IGAKnot3P" mu coord gauss [Vectror od doubles) only pr volumes
// "IGAWeight" (vector) for all
#define RESPONSETYPE_IGAOrder 2
#define RESPONSETYPE_IGAKnot1P 3
#define RESPONSETYPE_IGAKnot2P 4
#define RESPONSETYPE_IGAKnot3P 5
#define RESPONSETYPE_IGAWeight 6

Response*
IGASurfacePatch::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response* theResponse = 0;

  opserr << "IGASurfacePatch::setResponse - start argv list - argc = "  << argc << endln;
  for (int i = 0; i < argc; ++i)
  {
    opserr << "argv[" << i << "] = " << argv[i] << endln;
  }
  opserr << "IGASurfacePatch::setResponse - end argv list - argc = "  << argc << endln;

  // 1 - OPEN THE ELEMENT TAG
  output.tag("ElementOutput");
  output.attr("eleType", "IGASurfacePatch");
  output.attr("eleTag", this->getTag());
  int numNodes = this->getNumExternalNodes();
  const ID& nodes = this->getExternalNodes();
  static char nodeData[32];

  for (int i = 0; i < numNodes; i++) {
      sprintf(nodeData, "node%d", i + 1);
      output.attr(nodeData, nodes(i));
  }

  if (strcmp(argv[0], "material") == 0 ){
      
        // const Vector& data(numNodes);
        for (int i = 0; i < numNodes; i++) {
            sprintf(nodeData, "P%d", i + 1);
            output.tag("ResponseType", nodeData);
        }
        theResponse =  new ElementResponse(this, RESPONSETYPE_IGASurfacePatch_NODELUMPED, Vector(numNodes));
  } else if (strcmp(argv[0], "IGAOrder" ) == 0) {
    theResponse =  new ElementResponse(this, RESPONSETYPE_IGAOrder, ID(2));  // IGASurfacePatch is 2d
  } else if (strcmp(argv[0], "IGAKnot1P" ) == 0) {
    theResponse =  new ElementResponse(this, RESPONSETYPE_IGAKnot1P, Vector(noPtsX));
  } else if (strcmp(argv[0], "IGAKnot2P" ) == 0) {
    theResponse =  new ElementResponse(this, RESPONSETYPE_IGAKnot2P, Vector(noPtsY));
  } else if (strcmp(argv[0], "IGAKnot3P" ) == 0) {
    theResponse =  new ElementResponse(this, RESPONSETYPE_IGAKnot3P, Vector(0)); // No 3rd dim in this patch
  } else if (strcmp(argv[0], "IGAWeight" ) == 0) {
    theResponse =  new ElementResponse(this, RESPONSETYPE_IGAWeight, Vector(noPts));
  } 
  output.endTag();

  return theResponse;
}


int
IGASurfacePatch::getResponse(int responseID, Information &eleInfo)
{
    if(responseID == RESPONSETYPE_IGAOrder)
    {
        ID orders(2);
        orders(0) = P;
        orders(1) = Q;
        eleInfo.setID(orders);
    }
    else if(responseID == RESPONSETYPE_IGAKnot1P)
    {
        eleInfo.setVector(uKnot);
    }
    else if(responseID == RESPONSETYPE_IGAKnot2P)
    {
        eleInfo.setVector(vKnot);
    }
    else if(responseID == RESPONSETYPE_IGAKnot3P)
    {
        eleInfo.setVector(Vector(0));  // Has no 3rd dimension...
    }
    else if(responseID == RESPONSETYPE_IGAWeight)
    {
        eleInfo.setVector(weights);
    }
    return 0;
}
