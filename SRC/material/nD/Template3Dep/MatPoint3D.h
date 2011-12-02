//=============================================================================
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              MatPoint3D.h
// CLASS:             MatPoint3D
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Gauss Point
// RETURN:
// VERSION:
// LANGUAGE:          C++.ver >= 3.0 (Borland.C++.ver=3.1||SUN.C++.ver=3.0.1)
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic
// PROGRAMMER:        Boris Jeremic
// DATE:              17 October 1994.
// UPDATE HISTORY:
//
//
//#                    Aug 2000 porting to OpenSees                            #
//=============================================================================
//
//
#ifndef MATPOINT3D_H
#define MATPOINT3D_H

#include <stresst.h>
#include <straint.h>

//#include <mmodel.h>
#include <EPState.h>
#include <NDMaterial.h>

class MatPoint3D
{
  private:
    short int  r_direction_point_number;
    short int  s_direction_point_number;
    short int  t_direction_point_number;
    double r;
    double s;
    double t;
    double rw;
    double sw;
    double tw;

  public: // dobro ovo su samo ustvari pointeri pa nema ekstra kopiranja !

    // *gpEPS is used to hold all state parameters and internal vars, instead of the 
    // following stresses, strain, and internal vars __Zhaohui 09-30-2000
    
    //Now no EPState needed. Each MatPoint has an NDMaterial
    //EPState        *gpEPS;
    
    //stresstensor * p_stress;
    //stresstensor * p_iterative_stress; // to be used for iterative nodal forces
    //double         q_ast_iterative;
    //straintensor * p_strain;
    //tensor  * TangentE;
    
    //Might be ElasticIsotropic3D or Template3Dep
    NDMaterial * matmodel;

  public:  
    //default constructor
    MatPoint3D(  short int INr_direction_point_number = 0,
               short int INs_direction_point_number = 0,
               short int INt_direction_point_number = 0,
               double r_coord = 0,
               double s_coord = 0,
               double t_coord = 0,
               double r_weight = 0,
               double s_weight = 0,
               double t_weight = 0,
               //EPState *eps    = 0,
               NDMaterial * p_mmodel = 0   
	       //stresstensor * p_INstress = 0,
               //stresstensor * p_INiterative_stress = 0,
               //double         IN_q_ast_iterative = 0.0,
               //straintensor * p_INstrain = 0,
               //tensor * p_Tangent_E_tensor = 0,
               );
        
    // Constructor 1
    ~MatPoint3D();

    void Initialize(short int INr_direction_point_number,
                    short int INs_direction_point_number,
                    short int INt_direction_point_number,
                    double r_coord,
                    double s_coord,
                    double t_coord,
                    double r_weight,
                    double s_weight,
                    double t_weight,
                    //stresstensor * p_INstress,
                    //stresstensor * p_INiterative_stress,
                    //double         IN_q_ast_iterative,
                    //straintensor * p_INstrain,
                    //tensor * p_Tangent_E_tensor,
                    //EPState  * EPS,
		    NDMaterial * p_mmodel
                   );                             


  public:
    short int GP_number_r(void) const;
    short int GP_number_s(void) const;
    short int GP_number_t(void) const;
   
    MatPoint3D * GP(void);

    double r_coordinate() const;
    double s_coordinate() const;
    double t_coordinate() const;
    
    double r_weight() const;
    double s_weight() const;
    double t_weight() const;
    
    //void setEPS(EPState *eps);
    //EPState *getEPS() const;
    NDMaterial* getNDMat() const;
    const char* getType (void) const;
    int getTag (void) const;
    double getrho() const;
    const stresstensor getStressTensor() const;
    const straintensor getStrainTensor() const;
    //Added Aug. 13, 2001 Joey
    const straintensor getPlasticStrainTensor() const;
    //Added 02-18-03 Joey
    double getpsi() const; //state parameter
    
    int commitState(void) ;
    int revertToLastCommit(void) ;
    int revertToStart(void) ;

    void report(char * msg) const;
    void reportpqtheta(char * msg) const;
    void reportTensor(char *msg) const;
};

#endif 
//
