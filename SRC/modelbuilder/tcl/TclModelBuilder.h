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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclModelBuilder.h,v $
                                                                        
                                                                        
// File: ~/modelbuilder/tcl/TclModelBuilder.h
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the class definition for TclModelBuilder.
// A TclModelBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3 
// framework. currently these elements include:
//	1) linear-elastic 2 and 3d beam-column elements
//	2) non-linear material truss
//	3) non-linear 2 and 3d fiber-beam-column elements
//
// What: "@(#) TclModelBuilder.h, revA"

#ifndef TclModelBuilder_h
#define TclModelBuilder_h

#include <ModelBuilder.h>

class SectionForceDeformation;
class SectionRepres;
class UniaxialMaterial;
class NDMaterial;
class TaggedObjectStorage;

class CrdTransf2d;
class CrdTransf3d;

#include <tcl.h>
#include <tk.h>

class TclModelBuilder : public ModelBuilder
{
  public:
    TclModelBuilder(Domain &theDomain,Tcl_Interp *interp, int ndm, int ndf);
    ~TclModelBuilder();    

    int buildFE_Model(void);
    int getNDM(void) const;	
    int getNDF(void) const;

    // methods needed for the truss and fiber-beam elements for
    // adding/getting uniaxial material objects
    int addUniaxialMaterial(UniaxialMaterial &theMaterial);
    UniaxialMaterial *getUniaxialMaterial(int tag);

	// methods needed for the continuum elements and generic section
	// models to add/get ND material models
	int addNDMaterial(NDMaterial &theMaterial);
    NDMaterial *getNDMaterial(int tag);
    
    // methods needed for the nonlinear beam column elements to
    // add/get section objects
    int addSection(SectionForceDeformation &theSection);
    SectionForceDeformation *getSection(int tag);    
    int addSectionRepres(SectionRepres &theSectionRepres);
    SectionRepres *getSectionRepres(int tag);

    // methods needed for the frame elements to
    // add/get cordinate transformation objects
    int addCrdTransf2d(CrdTransf2d &theCrdTransf);
    int addCrdTransf3d(CrdTransf3d &theCrdTransf);
    CrdTransf2d *getCrdTransf2d(int tag);    
    CrdTransf3d *getCrdTransf3d(int tag);    

  protected:

  private:
    int ndm;	// space dimension of the mesh
    int ndf;	// number of degrees of freedom per node

    TaggedObjectStorage *theUniaxialMaterials;
    TaggedObjectStorage *theNDMaterials;
    TaggedObjectStorage *theSections;   
    TaggedObjectStorage *theSectionRepresents;
    TaggedObjectStorage *the2dGeomTransfs;
    TaggedObjectStorage *the3dGeomTransfs;
};

#endif







