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
// $Date: 2012-05-21 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/LayeredShellFiberSection.h,v $

// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Layered Shell Section
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

#include <SectionForceDeformation.h>


class LayeredShellFiberSection : public SectionForceDeformation{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    LayeredShellFiberSection( ) ;

    //full constructor
    LayeredShellFiberSection(   int tag, 
                                int iLayers, 
                                double *thickness, 
                                NDMaterial **fibers );

    const char *getClassType(void) const {return "LayeredShellFiberSection";};

    //destructor
    virtual ~LayeredShellFiberSection( ) ;

    //make a clone of this material
    SectionForceDeformation *getCopy( ) ;

    //mass per unit area
    double getRho() ;

    //send back order of strain in vector form
    int getOrder( ) const ;

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    //send back order of strain in vector form
    const ID& getType( ) ;

    //swap history variables
    int commitState( ) ; 

    //revert to last saved state
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    //get the strain and integrate plasticity equations
    int setTrialSectionDeformation( const Vector &strain_from_element ) ;

    //send back the strain
    const Vector& getSectionDeformation( ) ;

    //send back the stress 
    const Vector& getStressResultant( ) ;

    //send back the tangent 
    const Matrix& getSectionTangent( ) ;

    //send back the initial tangent 
    const Matrix& getInitialTangent( ) {return this->getSectionTangent();}

    //print out data
    void Print( OPS_Stream &s, int flag ) ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);


  private :

    int nLayers;
    //quadrature data
    double *sg;
    double *wg;

    double h ; //plate thickness

    NDMaterial **theFibers;  //pointers to the materials (fibers)

    Vector strainResultant ;

    static Vector stressResultant ;

    static Matrix tangent ;

    static ID array ;  

} ; //end of LayeredShellFiberSection declarations





