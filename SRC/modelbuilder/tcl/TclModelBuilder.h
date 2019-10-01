/* ****************************************************************** **
**    Opensees - Open System for Earthquake Engineering Simulation    **
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
                                                                        
// $Revision: 1.9 $
// $Date: 2009-04-17 22:58:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclModelBuilder.h,v $
                                                                        
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
class NDMaterial;
class TaggedObjectStorage;
class YieldSurface_BC;
class YS_Evolution;
class PlasticHardeningMaterial;
class CyclicModel;	 //!!
class DamageModel;
class FrictionModel;

#include <tcl.h>

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
    // REMOVED    int addUniaxialMaterial(UniaxialMaterial &theMaterial);
    //            UniaxialMaterial *getUniaxialMaterial(int tag);

    // methods needed for the continuum elements and generic section
    // models to add/get ND material models

    //    int addNDMaterial(NDMaterial &theMaterial);
    //    NDMaterial *getNDMaterial(int tag);
    
    // methods needed for the nonlinear beam column elements to
    // add/get section objects
    int addSection(SectionForceDeformation &theSection);
    SectionForceDeformation *getSection(int tag);    
    int addSectionRepres(SectionRepres &theSectionRepres);
    SectionRepres *getSectionRepres(int tag);

    // methods needed for the yield surfaces
    int addYieldSurface_BC(YieldSurface_BC &theYS);
    YieldSurface_BC *getYieldSurface_BC(int tag);
    int addYS_EvolutionModel(YS_Evolution &theModel);
    YS_Evolution *getYS_EvolutionModel(int tag);
    int addPlasticMaterial(PlasticHardeningMaterial &theMaterial);
    PlasticHardeningMaterial *getPlasticMaterial(int tag);
    int addCyclicModel(CyclicModel &theModel); //!!
    CyclicModel *getCyclicModel(int tag); //!!
    int addDamageModel(DamageModel &theModel); //!!
    DamageModel *getDamageModel(int tag); //!!


    // methods needed for the friction models
    // int addFrictionModel(FrictionModel &theFrnMdl);
    // FrictionModel *getFrictionModel(int tag);

  private:
    int ndm;	// space dimension of the mesh
    int ndf;	// number of degrees of freedom per node

    //    TaggedObjectStorage *theUniaxialMaterials;
    TaggedObjectStorage *theNDMaterials;
    TaggedObjectStorage *theSections;   
    TaggedObjectStorage *theSectionRepresents;
    TaggedObjectStorage *theYieldSurface_BCs;
    TaggedObjectStorage *thePlasticMaterials;
    TaggedObjectStorage *theYS_EvolutionModels;
    TaggedObjectStorage *theCycModels; //!!
    //    TaggedObjectStorage *theDamageModels; //!!
    //    TaggedObjectStorage *theFrictionModels;
    TaggedObjectStorage *theLimitCurves;	//MRL

 protected:
    Tcl_Interp *theInterp;
};

#endif







