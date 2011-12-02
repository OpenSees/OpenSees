//=============================================================================
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              MatPoint3D.cpp
// CLASS:             MatPoint3D
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Material Point
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
//
//=============================================================================
//
#ifndef MATPOINT3D_CPP
#define MATPOINT3D_CPP

#include "MatPoint3D.h"

//=============================================================================
// Constructor
MatPoint3D::MatPoint3D(short int INr_direction_point_number,
                       short int INs_direction_point_number,
                       short int INt_direction_point_number,
                       double r_coord,
                       double s_coord,
                       double t_coord,
                       double r_weight,
                       double s_weight,
                       double t_weight,
                       NDMaterial * p_INmatmodel
                       //stresstensor * p_INstress,
                       //stresstensor * p_INiterative_stress,
                       //double         IN_q_ast_iterative,
                       //straintensor * p_INstrain,
                       //tensor * p_Tangent_E_tensor,
                      )
  {
    r_direction_point_number = INr_direction_point_number;
    s_direction_point_number = INs_direction_point_number;
    t_direction_point_number = INt_direction_point_number;
    r = r_coord;
    s = s_coord;
    t = t_coord;
    rw = r_weight;
    sw = s_weight;
    tw = t_weight;

    //if ( eps ) {
    //    gpEPS = eps->newObj();  // Elasto-plastic 3D
    //}
    //else 
    //    gpEPS = 0;
    
    if (p_INmatmodel)
        matmodel = p_INmatmodel->getCopy();
    else 
        matmodel = 0;		

    //p_stress   = p_INstress; 
    //p_iterative_stress   = p_INiterative_stress;
    //q_ast_iterative      = IN_q_ast_iterative;
    //p_strain   = p_INstrain;
    //TangentE   = p_Tangent_E_tensor;

  }

////=============================================================================
MatPoint3D::~MatPoint3D() 
{
    if ( matmodel )
      //delete [] matmodel; //bug found: Inconsistent delete
      delete matmodel;

    //if ( GP )
    //  delete [] GP;

}

////=============================================================================
//MatPoint3D::MatPoint3D(EPState *eps) {
//
//    this->r_direction_point_number = 0;
//    this->s_direction_point_number = 0;
//    this->t_direction_point_number = 0;
//    this->r = 0;
//    this->s = 0;
//    this->t = 0;
//    this->rw = 0;
//    this->sw = 0;
//    this->tw = 0;
//	     
//    if ( eps ) 
//        gpEPS = eps->newObj();      
//    else
//        //gpEPS = new EPState();  //otherwise use default EPState __Zhaohui 09-30-2000
//        gpEPS = 0;  //otherwise use default EPState __Zhaohui 09-30-2000
//    
//}

//=============================================================================
void MatPoint3D::Initialize( short int INr_direction_point_number,
                             short int INs_direction_point_number,
                             short int INt_direction_point_number,
                             double r_coord,
                             double s_coord,
                             double t_coord,
                             double r_weight,
                             double s_weight,
                             double t_weight,
                             //EPState *EPS,
			     //stresstensor * p_INstress,
                             //stresstensor * p_INiterative_stress,
                             //double         IN_q_ast_iterative,
                             //straintensor * p_INstrain,
                             //tensor *       p_Tangent_E_tensor,
                             NDMaterial * p_INmmodel
                           )
  {
    this->r_direction_point_number = INr_direction_point_number;
    this->s_direction_point_number = INs_direction_point_number;
    this->t_direction_point_number = INt_direction_point_number;

    this->r = r_coord;
    this->s = s_coord;
    this->t = t_coord;
    this->rw = r_weight;
    this->sw = s_weight;
    this->tw = t_weight;

    //if (EPS)
    //   this->gpEPS = EPS->newObj();
    //else 
    //   this->gpEPS = 0;

    if ( p_INmmodel )
       this->matmodel = p_INmmodel;
    else 
       this->matmodel = 0;

    //this->p_stress   = p_INstress;
    //this->p_iterative_stress   = p_INiterative_stress;
    //q_ast_iterative       = IN_q_ast_iterative;
    //this->p_strain   = p_INstrain;
    //this->TangentE   = p_Tangent_E_tensor;
  }


//=============================================================================
//return MatPoint
MatPoint3D * MatPoint3D::GP(void) 
     {
       return this;
     }

//=============================================================================
short int MatPoint3D::GP_number_r(void) const {
    return r_direction_point_number; 
}


//=============================================================================
short int MatPoint3D::GP_number_s(void) const { 
    return s_direction_point_number; 
}

//=============================================================================
short int MatPoint3D::GP_number_t(void) const {
    return t_direction_point_number; 
}

//=============================================================================
double MatPoint3D::r_coordinate() const {
    return r; 
}

//=============================================================================
double MatPoint3D::s_coordinate() const {
    return s; 
}

//=============================================================================
double MatPoint3D::t_coordinate() const {
    return t; 
}

//=============================================================================
double MatPoint3D::r_weight() const {
    return rw; 
}

//=============================================================================
double MatPoint3D::s_weight() const {
    return sw; 
}

//=============================================================================
double MatPoint3D::t_weight() const {
    return tw; 
}

////=============================================================================
//void MatPoint3D::setEPS(EPState *eps) {
//
//    if ( eps ) 
//        gpEPS = eps->newObj();  
//    else 
//        gpEPS = 0;
//	//g3ErrorHandler->warning("MatPoint3D::setEPS  No initial values for EPState, using default");
//
//}


////=============================================================================
//EPState* MatPoint3D::getEPS() const {
//
//    return gpEPS;
//}
//


//=============================================================================
NDMaterial* MatPoint3D::getNDMat() const {

    return matmodel;
}

//=============================================================================
double MatPoint3D::getrho() const 
{
    return matmodel->getRho(); 
}

const char* MatPoint3D::getType (void) const
{
	return matmodel->getType();
}

int MatPoint3D::getTag (void) const
{
	return matmodel->getTag();
}

//=============================================================================
const stresstensor MatPoint3D::getStressTensor() const {

    return matmodel->getStressTensor();
}


//=============================================================================
const straintensor MatPoint3D::getStrainTensor() const {

    return matmodel->getStrainTensor();
}

//=============================================================================
const straintensor MatPoint3D::getPlasticStrainTensor() const {

    return matmodel->getPlasticStrainTensor();
}

//=============================================================================
double MatPoint3D::getpsi() const 
{
    return matmodel->getpsi();
}


//================================================================================
int MatPoint3D::commitState(void)
{
	int err;
	err = matmodel->commitState();
	return err;
}

//================================================================================
int MatPoint3D::revertToLastCommit(void)
{
	int err;
	err = matmodel->revertToLastCommit();
	return err;
}

//================================================================================
int MatPoint3D::revertToStart(void)
{
	int err;
	err = matmodel->revertToStart();
	return err;
}

//=============================================================================
void MatPoint3D::report(char *msg) const
  {
    //if ( msg ) ::printf("%s",msg);
    if ( msg )  opserr << msg;
    
    ::printf("\n\n\n---------------------------------------------------- \n");
    ::printf("Gauss point #r %d #s %d  #t %d \n",
                                 GP_number_r(), GP_number_s(), GP_number_t());

    ::printf("\tr->%.8e   s->%.8e   t->%.8e  \n",
     r_coordinate(),s_coordinate(),t_coordinate());

    ::printf("\tr_weight->%.8e   s_weight->%.8e   t_weight->%.8e  \n",
     r_weight(),s_weight(),t_weight());


    if ( matmodel ) 
      {
        opserr << (*matmodel);

        stresstensor tmpstress = matmodel->getStressTensor();
        tmpstress.report("stress at this Gauss-Legendre point\n");

        straintensor tmpstrain = matmodel->getStrainTensor();
        tmpstrain.report("strain at this Gauss-Legendre point\n");
      }
    else
      opserr << "Empty Material Model\n"; 

    //p_stress->report("stress at this Gauss-Legendre point\n");
    //p_iterative_stress->reportshortpqtheta("ITERATIVE stress at this Gauss-Legendre point\n");
    //::printf("ITERATIVE q_ast_iterative = %.8e  \n",q_ast_iterative);
    //p_strain->report("strain at this Gauss-Legendre point\n");
    //matmodel->report("material model at this Gauss-Legendre point\n");
    
    //if (gpEPS) 
    //  opserr << (*gpEPS);
    //else
    //  opserr << "Empty EPState\n"; 
    
  }

//=============================================================================
void MatPoint3D::reportpqtheta(char *msg) const
  {
    //if ( msg ) ::printf("%s",msg);
    //p_stress->reportshortpqtheta("");
    //p_stress->reportshortpqtheta("");
    if ( msg )  opserr << msg;
    
    //if ( gpEPS ) {
    //   ( gpEPS->getStress() ).reportshortpqtheta("");
    //   ( gpEPS->getStress() ).reportSHORTs1s2s3("");
    //}
   
    if ( matmodel ) {
       stresstensor tmp = matmodel->getStressTensor();
       tmp.reportshortpqtheta("");
       tmp.reportSHORTs1s2s3("");
    }

  }

//=============================================================================
void MatPoint3D::reportTensor(char *msg) const
  {
    //if ( msg ) ::printf("%s",msg);
    //p_stress->reportTensor("");

    if ( msg )  opserr << msg;
    
    //if (gpEPS )
    //   ( gpEPS->getStress() ).reportTensor("");

    if ( matmodel ) {
       stresstensor tmp = matmodel->getStressTensor();
       tmp.reportTensor("");
    }

  }


#endif




