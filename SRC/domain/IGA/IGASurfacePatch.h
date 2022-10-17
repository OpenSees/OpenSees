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

// Written: Felipe Elgueta and José A. Abell (UANDES, Chile) www.joseabell.com

// Description: This file contains the class definition for IGASurfacePatch.
//
//


#ifndef IGASurfacePatch_h
#define IGASurfacePatch_h


#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <ElementalLoad.h>
#include <Subdomain.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>

// Added by Pipe
#include <ErrorHandler.h>

// For parameters (damage model Paepegem)
#include <Information.h>
#include <Parameter.h>



typedef enum 
{
    KLShell, 
    KLShellLargeDef,
    KLShell_BendingStrip
} ShellType;

class IGASurfacePatch: public Subdomain
{
public:
    IGASurfacePatch(int tag, int P_, const Vector& uKnot_, const Matrix& controlPts_);
    IGASurfacePatch(int tag, int P_, int Q_, const Vector& uKnot_, const Vector& vKnot_, const Matrix& controlPts_); // I want to delete this, but execution fails if not declared
    IGASurfacePatch(int tag, int P_, int Q_, int noPtsX_, int noPtsY_, const Vector& uKnot_, const Vector& vKnot_, const Matrix& controlPts_, ShellType shtype, NDMaterial** ref_planestres_materials); // testing
    IGASurfacePatch(int tag, int nodeStartTag_, int P_, int Q_, int noPtsX_, int noPtsY_, int nonlinearGeometry_, const Vector gFact_, const ID& matTags_, const Vector&  theta_, const Vector& thickness_, const Vector& uKnot_, const Vector& vKnot_, const Matrix& controlPts_, ShellType shtype);
    IGASurfacePatch(int tag, int P_, int Q_, int R_, const Vector& uKnot_, const Vector& vKnot_, const Vector& wKnot_, const Matrix& controlPts_);
    
    virtual ~IGASurfacePatch();

    const char *getClassType(void) const {return "IGASurfacePatch";};

    // NURBS member functions
    int Nurbs2DBasis2ndDers(double xi, double eta, Vector& R, Vector& dRdxi, Vector& dRdeta, Vector& dR2dxi, Vector& dR2deta, Vector& dR2dxideta);
    double parent2ParametricSpace(Vector range, double xibar);
    int getNoFuncs();
    ID getOrders();

    // // Composite layer functions
    int getNLayers();
    int getMatTag(int iLayer);
    double getThickness(int iLayer); 
    double getAngle(int iLayer);
    double getZk(int iLayer);

    // Utility functions
    bool getAnalysisType(); // To know if perfmorming a linear o non-linear analysis, non-linear default
    Vector getGravityFactors(); // Function to obtain ther gravity in xyz



    // Subdomain/domain member function overloads.
    void Print(OPS_Stream &s, int flag = 0);
    void setDomain(Domain *theDomain);  

    void clearAll(void) ;

    virtual int addLoad(ElementalLoad *theLoad, double loadFactor);

    Response* setResponse( const char **argv, int argc, OPS_Stream &output );
    int getResponse( int responseID, Information &eleInfo );

private:
    bool generateIGA2DMesh(int & noElems, int & noElemsU, int & noElemsV);//, Matrix & index, Matrix & elRangeU, Matrix & elRangeV, Matrix & elConnU, Matrix & elConnV);
    bool buildConnectivity(int p, const Vector & knotVec, int nE, Matrix * elRange, Matrix * elConn);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

private:

    int nodeStartTag;       // First node tag (for multipatches)
    int P;                  // Order in the U direction
    int Q;                  // Order in the V direction
    int R;                  // Order in the W direction
    int noPtsX;             // Number of points in U direction
    int noPtsY;             // Number of points in V direction
    int nonLinearGeometry;
    int noPts;              // Number of control points (noPtsX*noPtsY)
    Vector gFact;           // Factor xyz for gravity
    Vector uKnot;           // Knot vector in the U direction
    Vector vKnot;           // Knot vector in the V direction
    Vector wKnot;       // Knot vector in the V direction
    Matrix controlPts;  // Control-points
    Vector weights;
    int noElemsU;        // Number of Elements in the U direction
    int noElemsV;        // Number of Elements in the V direction
    int noElems;        // Number of Elements
    int noFuncs;        //Number of shape functions


    ShellType shtype;
    ID matTags;

    Matrix* index;
    Matrix* elRangeU;
    Matrix* elRangeV;
    Matrix* element;
    Matrix* elConnU;
    Matrix* elConnV;

    Matrix* quadPoint;
    Vector* quadWeight;

    Domain * opsDomain;

    // Guardar vector de espesores
    Vector thickness; 
    // Guardar vector de orientaciones
    Vector theta; 
    // Guardar puntero a vector de posiciones
    Vector* Zk;



};

#endif
