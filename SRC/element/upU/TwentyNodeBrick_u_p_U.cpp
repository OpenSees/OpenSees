///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              TwentyNodeBrick_u_p_U.cpp
// CLASS:             TwentyNodeBrick_u_p_U
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Xiaoyan Wu
// DATE:              Sep. 2001
// UPDATE HISTORY:    Modified from EightNodeBrick_u_p_U.cpp Sep. 2001
//
//          01/16/2002    Xiaoyan
//          Add the permeability tensor and ks, kf  to the constructor
//
//
//  "Coupled system": Solid and fluid coexist.
//                    u-- Solid displacement
//                    p-- Pore pressure
//                    U-- Absolute fluid displacement
//
//                    31Oct2003. Qing fixed small inconsistencies in basic theory
//                               related to permeability...
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef TWENTYNODEBRICK_U_P_U_CPP
#define TWENTYNODEBRICK_U_P_U_CPP

#include <TwentyNodeBrick_u_p_U.h>
#define FixedOrder 3

// Changed to static data members on 01/16/2002

Matrix TwentyNodeBrick_u_p_U::K(140, 140);
Matrix TwentyNodeBrick_u_p_U::C(140, 140);
Matrix TwentyNodeBrick_u_p_U::M(140, 140);
Vector TwentyNodeBrick_u_p_U::p(140);
tensor TwentyNodeBrick_u_p_U::k(2,def_dim_2,0.0);
Node * TwentyNodeBrick_u_p_U::theNodes[20];

//=========================================================================
// Constructor. The dimension of K, C and M are(140,140)     Wxy 08/28/2001
//=========================================================================

TwentyNodeBrick_u_p_U::TwentyNodeBrick_u_p_U(int element_number,
                               int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
                               int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8,
                               int node_numb_9,  int node_numb_10, int node_numb_11, int node_numb_12,
                               int node_numb_13, int node_numb_14, int node_numb_15, int node_numb_16,
                               int node_numb_17, int node_numb_18, int node_numb_19, int node_numb_20,
                               NDMaterial * Globalmmodel, double b1, double b2,double b3,
             double nn, double alf, double rs, double rf,
             double permb_x,double permb_y,double permb_z,
             double kks, double kkf, double pp)
             // wxy added rs and rf for the solid and fluid density    08/28/2001

             //, EPState *InitEPS)  const char * type,
                               // Set it to 3 //int r_int_order, //int s_int_order, //int t_int_order,
             //tensor * IN_tangent_E,  //stresstensor * INstress, //stresstensor * INiterative_stress, //double * IN_q_ast_iterative, //straintensor * INstrain):  __ZHaohui 09-29-2000

  :Element(element_number, ELE_TAG_TwentyNodeBrick_u_p_U ),
   connectedExternalNodes(20), Q(140), bf(3),
   n(nn), alpha(alf), rho_s(rs), rho_f(rf), ks(kks), kf(kkf), pressure(pp), Ki(0)
  {
    //elem_numb = element_number;
    rho=(1-n)*rho_s+n*rho_f;
    bf(0) = b1;
    bf(1) = b2;
    bf(2) = b3;

    // permeability
    k.val(1,1)=permb_x;
    k.val(2,2)=permb_y;
    k.val(3,3)=permb_z;
    //k.print("k","\n test the permeability tensor k \n");    // out 06/26/2002

    determinant_of_Jacobian = 0.0;

    //r_integration_order = r_int_order;
    //s_integration_order = s_int_order;
    //t_integration_order = t_int_order;
    r_integration_order = FixedOrder; // Gauss-Legendre integration order in r direction
    s_integration_order = FixedOrder; // Gauss-Legendre integration order in s direction
    t_integration_order = FixedOrder; // Gauss-Legendre integration order in t direction

    //Not needed. Right now we have one NDMaterial for each material point
    //mmodel = Globalmmodel->getCopy( type ); // One global mat model

    int total_number_of_Gauss_points = r_integration_order*s_integration_order*t_integration_order;

    // according to ARM pp.61 default constructor will be called!!
    //MatPoint3D * matpoint = new MatPoint3D[total_number_of_Gauss_points];
    //prebaci sve u jednodimenzioni niz jer samo prvi stepen pointera moze da se pokriva
    //sa onim stosom derived * ->> base * !!

    if ( total_number_of_Gauss_points != 0 )
      {

  matpoint  = new MatPoint3D * [total_number_of_Gauss_points];

      }
    else
      {
  //GPstress = 0;//GPiterative_stress = 0;//GPq_ast_iterative  = 0; //GPstrain = 0;//GPtangent_E = 0;
        matpoint  = 0;
    }

    short where = 0;

    short GP_c_r;
        for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
        {
        double r = get_Gauss_p_c( r_integration_order, GP_c_r );
        double rw = get_Gauss_p_w( r_integration_order, GP_c_r );

    short GP_c_s;
        for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
        {
         double s = get_Gauss_p_c( s_integration_order, GP_c_s );
         double sw = get_Gauss_p_w( s_integration_order, GP_c_s );

     short GP_c_t;
        for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
         {
          double t = get_Gauss_p_c( t_integration_order, GP_c_t );
          double tw = get_Gauss_p_w( t_integration_order, GP_c_t );

                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;


                //cout << "where = " << where << endln;
                matpoint[where] = new MatPoint3D(GP_c_r,
                                                 GP_c_s,
                                                 GP_c_t,
                                                 r, s, t,
                                                 rw, sw, tw,
                                               //InitEPS,
                                                 Globalmmodel);
           //NMD);
           //&( GPstress[where] ), //&( GPiterative_stress[where] ), //IN_q_ast_iterative[where] ,//&( GPstrain[where] ),  //&( GPtangent_E[where] ),
                                         //&( (matpoint)->operator[](where) )
                                         // ugly syntax but it works! Still don't know what's wrong   // with the old style matpoint[where]
              }
          }
      }

      // Set connected external node IDs
      connectedExternalNodes(0) = node_numb_1;
      connectedExternalNodes(1) = node_numb_2;
      connectedExternalNodes(2) = node_numb_3;
      connectedExternalNodes(3) = node_numb_4;

      connectedExternalNodes(4) = node_numb_5;
      connectedExternalNodes(5) = node_numb_6;
      connectedExternalNodes(6) = node_numb_7;
      connectedExternalNodes(7) = node_numb_8;

      connectedExternalNodes( 8) = node_numb_9;
      connectedExternalNodes( 9) = node_numb_10;
      connectedExternalNodes(10) = node_numb_11;
      connectedExternalNodes(11) = node_numb_12;

      connectedExternalNodes(12) = node_numb_13;
      connectedExternalNodes(13) = node_numb_14;
      connectedExternalNodes(14) = node_numb_15;
      connectedExternalNodes(15) = node_numb_16;

      connectedExternalNodes(16) = node_numb_17;
      connectedExternalNodes(17) = node_numb_18;
      connectedExternalNodes(18) = node_numb_19;
      connectedExternalNodes(19) = node_numb_20;

      nodes_in_brick = 20;


}

////#############################################################################
//=================================================================================
// Default Constructor. The dimension of K, C and M are(140,140)     Wxy 08/28/2001
//=================================================================================

TwentyNodeBrick_u_p_U::TwentyNodeBrick_u_p_U ()
  :Element(0, ELE_TAG_TwentyNodeBrick_u_p_U ),
   connectedExternalNodes(20), Q(140), bf(3),
   n(0), alpha(1), rho_s(0.0),rho_f(0.0), ks(0.0), kf(0.0), pressure(0.0), mmodel(0), Ki(0)
{
     matpoint = 0;
}

////#############################################################################
//=========================================================================
// Destructor.
//=========================================================================

TwentyNodeBrick_u_p_U::~TwentyNodeBrick_u_p_U ()
{

    int total_number_of_Gauss_points = r_integration_order*s_integration_order*t_integration_order;

    int i;
    for (i = 0; i < total_number_of_Gauss_points; i++)
    {
  // Delete the NDMaterials at each integration point
  if (matpoint[i])
      delete matpoint[i];
    }

    // Delete the array of pointers to NDMaterial pointer arrays
    if (matpoint)
      delete [] matpoint;

    if (Ki != 0)
      delete Ki;
}

//=========================================================================
// Shape functions in "element global level". dimension are(60,3)
// Since we define the mass Matrix Mf or Ms as four order tensor this
// function will not be used.  Just keep here. Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::H_3D(double r1, double r2, double r3)
  {

    int dimension[] = {60,3};

    tensor H(2, dimension, 0.0);

    // influence of the node number 20
        H.val(58,1)=(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
        H.val(59,2)=(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
        H.val(60,3)=(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 19
        H.val(55,1)=(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
        H.val(56,2)=(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
        H.val(57,3)=(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 18
        H.val(52,1)=(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
        H.val(53,2)=(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
        H.val(54,3)=(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 17
        H.val(49,1)=(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
        H.val(50,2)=(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
        H.val(51,3)=(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;

    // influence of the node number 16
        H.val(46,1)=(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
        H.val(47,2)=(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
        H.val(48,3)=(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 15
        H.val(43,1)=(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
        H.val(44,2)=(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
        H.val(45,3)=(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    // influence of the node number 14
        H.val(40,1)=(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
        H.val(41,2)=(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
        H.val(42,3)=(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 13
        H.val(37,1)=(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;
        H.val(38,2)=(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;
        H.val(39,3)=(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;

    // influence of the node number 12
        H.val(34,1)=(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
        H.val(35,2)=(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
        H.val(36,3)=(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 11
        H.val(31,1)=(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
        H.val(32,2)=(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
        H.val(33,3)=(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    // influence of the node number 10
        H.val(28,1)=(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
        H.val(29,2)=(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
        H.val(30,3)=(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 9
        H.val(25,1)=(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;
        H.val(26,2)=(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;
        H.val(27,3)=(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;

    // influence of the node number 8
    H.val(22,1)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(43,1)+H.val(48,3)+H.val(60,3))/2.0;
    H.val(23,2)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(43,1)+H.val(48,3)+H.val(60,3))/2.0;
    H.val(24,3)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(43,1)+H.val(48,3)+H.val(60,3))/2.0;
    // influence of the node number 7
    H.val(19,1)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(42,3)+H.val(43,1)+H.val(57,3))/2.0;
    H.val(20,2)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(42,3)+H.val(43,1)+H.val(57,3))/2.0;
    H.val(21,3)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(42,3)+H.val(43,1)+H.val(57,3))/2.0;
    // influence of the node number 6
    H.val(16,1)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0 - (H.val(39,3)+H.val(42,3)+H.val(54,3))/2.0;
    H.val(17,2)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0 - (H.val(39,3)+H.val(42,3)+H.val(54,3))/2.0;
    H.val(18,3)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0 - (H.val(39,3)+H.val(42,3)+H.val(54,3))/2.0;
    // influence of the node number 5
    H.val(13,1)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0 - (H.val(39,3)+H.val(48,3)+H.val(51,3))/2.0;
    H.val(14,2)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0 - (H.val(39,3)+H.val(48,3)+H.val(51,3))/2.0;
    H.val(15,3)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0 - (H.val(39,3)+H.val(48,3)+H.val(51,3))/2.0;

    // influence of the node number 4
    H.val(10,1)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0 - (H.val(33,3)+H.val(36,3)+H.val(60,3))/2.0;
    H.val(11,2)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0 - (H.val(33,3)+H.val(36,3)+H.val(60,3))/2.0;
    H.val(12,3)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0 - (H.val(33,3)+H.val(36,3)+H.val(60,3))/2.0;
    // influence of the node number 3
    H.val(7,1) =(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0 - (H.val(30,3)+H.val(33,3)+H.val(57,3))/2.0;
    H.val(8,2) =(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0 - (H.val(30,3)+H.val(33,3)+H.val(57,3))/2.0;
    H.val(9,3) =(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0 - (H.val(30,3)+H.val(33,3)+H.val(57,3))/2.0;
    // influence of the node number 2
    H.val(4,1) =(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0 - (H.val(30,3)+H.val(54,3)+H.val(27,3))/2.0;
    H.val(5,2) =(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0 - (H.val(30,3)+H.val(54,3)+H.val(27,3))/2.0;
    H.val(6,3) =(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0 - (H.val(30,3)+H.val(54,3)+H.val(27,3))/2.0;
    // influence of the node number 1
    H.val(1,1) =(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0 - (H.val(36,3)+H.val(51,3)+H.val(27,3))/2.0;
    H.val(2,2) =(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0 - (H.val(36,3)+H.val(51,3)+H.val(27,3))/2.0;
    H.val(3,3) =(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0 - (H.val(36,3)+H.val(51,3)+H.val(27,3))/2.0;

    //         double sum = 0;
    //
    //   for (int i=1; i<=60 ; i++)
    //           {
    // //        sum+=H.cval(i,1);
    //       for (int j=1; j<= 1; j++)
    //          {
    //                    sum+=H.cval(i,1);
    //             ::printf( "  %+9.2e", H.cval(i,j) );
    //           }
    //            // ::printf( "  %d \n", i);
    //      }
    //       ::printf( " \n sum= %+6.2e\n", sum );


    //    printf("r1 = %lf, r2 = %lf, r3 = %lf\n", r1, r2, r3);
    //    H.print("h");

    return H;
  }
//==============================================================================
// Derivative of Shape functions in "element golbal level". dimension are(60,3)
// Xiaoyan added this function  08/28/2001
// Since we define the mass Matrix G1 or G2 as four order tensor this
// function will not be used.  Just keep here. Wxy 09/26/2001
//==============================================================================
tensor TwentyNodeBrick_u_p_U::dH_drst_at(double r1, double r2, double r3)
  {

    int dimension[] = {60,3}; // Xiaoyan changed from {60,3} to {24,3} for 8 nodes
                              // 3*8=24  07/12/00
    tensor dH(2, dimension, 0.0);

    // influence of the node number 20
    //    dH.val(58,1)=node_existance[20-1-8]*(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    dH.val(59,2)=node_existance[20-1-8]*(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    dH.val(60,3)=node_existance[20-1-8]*(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 19
    //    dH.val(55,1)=node_existance[19-1-8]*(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    dH.val(56,2)=node_existance[19-1-8]*(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    dH.val(57,3)=node_existance[19-1-8]*(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 18
    //    dH.val(52,1)=node_existance[18-1-8]*(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    dH.val(53,2)=node_existance[18-1-8]*(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    dH.val(54,3)=node_existance[18-1-8]*(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 17
    //    dH.val(49,1)=node_existance[17-1-8]*(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    dH.val(50,2)=node_existance[17-1-8]*(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    dH.val(51,3)=node_existance[17-1-8]*(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;

    // influence of the node number 16
    //    dH.val(46,1)=node_existance[16-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    dH.val(47,2)=node_existance[16-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    dH.val(48,3)=node_existance[16-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 15
    //    dH.val(43,1)=node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    //    dH.val(44,2)=node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    //    dH.val(45,3)=node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    // influence of the node number 14
    //    dH.val(40,1)=node_existance[14-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    dH.val(41,2)=node_existance[14-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    dH.val(42,3)=node_existance[14-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 13
    //    dH.val(37,1)=node_existance[13-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;
    //    dH.val(38,2)=node_existance[13-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;
    //    dH.val(39,3)=node_existance[13-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;

    // influence of the node number 12
    //    dH.val(34,1)=node_existance[12-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    dH.val(35,2)=node_existance[12-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    dH.val(36,3)=node_existance[12-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 11
    //    dH.val(31,1)=node_existance[11-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    //    dH.val(32,2)=node_existance[11-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    //    dH.val(33,3)=node_existance[11-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    // influence of the node number 10
    //    dH.val(28,1)=node_existance[10-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    dH.val(29,2)=node_existance[10-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    dH.val(30,3)=node_existance[10-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 9
    //    dH.val(25,1)=node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;
    //    dH.val(26,2)=node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;
    //    dH.val(27,3)=node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;
    //

    // 9-20 nodes commented by Xiaoyan  07/12/00

    // influence of the node number 8
    //    dH.val(22,1)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (dH.val(15)+dH.val(16)+dH.val(20))/2.0;
    //    dH.val(23,2)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (dH.val(15)+dH.val(16)+dH.val(20))/2.0;
    //    dH.val(24,3)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (dH.val(15)+dH.val(16)+dH.val(20))/2.0;
    dH.val(22,1)= (1.0-r2)*(1.0-r3)/8.0;// - (dH.val(43,1)+dH.val(48,3)+dH.val(60,3))/2.0;
    dH.val(23,2)=-(1.0+r1)*(1.0-r3)/8.0;// - (dH.val(43,1)+dH.val(48,3)+dH.val(60,3))/2.0;
    dH.val(24,3)=-(1.0+r1)*(1.0-r2)/8.0;// - (dH.val(43,1)+dH.val(48,3)+dH.val(60,3))/2.0;
    // influence of the node number 7
    dH.val(19,1)=-(1.0-r2)*(1.0-r3)/8.0;// - (dH.val(42,3)+dH.val(43,1)+dH.val(57,3))/2.0;
    dH.val(20,2)=-(1.0-r1)*(1.0-r3)/8.0;// - (dH.val(42,3)+dH.val(43,1)+dH.val(57,3))/2.0;
    dH.val(21,3)=-(1.0-r1)*(1.0-r2)/8.0;// - (dH.val(42,3)+dH.val(43,1)+dH.val(57,3))/2.0;
    // influence of the node number 6
    dH.val(16,1)=-(1.0+r2)*(1.0-r3)/8.0 ;//- (dH.val(39,3)+dH.val(42,3)+dH.val(54,3))/2.0;
    dH.val(17,2)= (1.0-r1)*(1.0-r3)/8.0 ;//- (dH.val(39,3)+dH.val(42,3)+dH.val(54,3))/2.0;
    dH.val(18,3)=-(1.0-r1)*(1.0+r2)/8.0 ;//- (dH.val(39,3)+dH.val(42,3)+dH.val(54,3))/2.0;
    // influence of the node number 5
    dH.val(13,1)= (1.0+r2)*(1.0-r3)/8.0 ;//- (dH.val(39,3)+dH.val(48,3)+dH.val(51,3))/2.0;
    dH.val(14,2)= (1.0+r1)*(1.0-r3)/8.0 ;//- (dH.val(39,3)+dH.val(48,3)+dH.val(51,3))/2.0;
    dH.val(15,3)=-(1.0+r1)*(1.0+r2)/8.0 ;//- (dH.val(39,3)+dH.val(48,3)+dH.val(51,3))/2.0;

    // influence of the node number 4
    dH.val(10,1)= (1.0-r2)*(1.0+r3)/8.0 ;//- (dH.val(33,3)+dH.val(36,3)+dH.val(60,3))/2.0;
    dH.val(11,2)=-(1.0+r1)*(1.0+r3)/8.0 ;//- (dH.val(33,3)+dH.val(36,3)+dH.val(60,3))/2.0;
    dH.val(12,3)= (1.0+r1)*(1.0-r2)/8.0 ;//- (dH.val(33,3)+dH.val(36,3)+dH.val(60,3))/2.0;
    // influence of the node number 3
    dH.val(7,1)=-(1.0-r2)*(1.0+r3)/8.0 ;//- (dH.val(30,3)+dH.val(33,3)+dH.val(57,3))/2.0;
    dH.val(8,2)=-(1.0-r1)*(1.0+r3)/8.0 ;//- (dH.val(30,3)+dH.val(33,3)+dH.val(57,3))/2.0;
    dH.val(9,3)= (1.0-r1)*(1.0-r2)/8.0 ;//- (dH.val(30,3)+dH.val(33,3)+dH.val(57,3))/2.0;
    // influence of the node number 2
    dH.val(4,1)=-(1.0+r2)*(1.0+r3)/8.0 ;//- (dH.val(30,3)+dH.val(54,3)+dH.val(27,3))/2.0;
    dH.val(5,2)= (1.0-r1)*(1.0+r3)/8.0 ;//- (dH.val(30,3)+dH.val(54,3)+dH.val(27,3))/2.0;
    dH.val(6,3)= (1.0-r1)*(1.0+r2)/8.0 ;//- (dH.val(30,3)+dH.val(54,3)+dH.val(27,3))/2.0;
    // influence of the node number 1
    dH.val(1,1)= (1.0+r2)*(1.0+r3)/8.0 ;//- (dH.val(36,3)+dH.val(51,3)+dH.val(27,3))/2.0;
    dH.val(2,2)= (1.0+r1)*(1.0+r3)/8.0 ;//- (dH.val(36,3)+dH.val(51,3)+dH.val(27,3))/2.0;
    dH.val(3,3)= (1.0+r1)*(1.0+r2)/8.0 ;//- (dH.val(36,3)+dH.val(51,3)+dH.val(27,3))/2.0;

                 // The second part were commented by Xiaoyan
    //         double sum = 0;
    //
    //   for (int i=1; i<=60 ; i++)
    //           {
    // //        sum+=dH.cval(i,1);
    //       for (int j=1; j<= 1; j++)
    //          {
    //                    sum+=dH.cval(i,1);
    //             ::printf( "  %+9.2e", dH.cval(i,j) );
    //           }
    //            // ::printf( "  %d \n", i);
    //      }
    //       ::printf( " \n sum= %+6.2e\n", sum );


    //    printf("r1 = %lf, r2 = %lf, r3 = %lf\n", r1, r2, r3);
    //    dH.print("h");

    return dH;
  }


//=========================================================================
// Shape functions in "element local level". dimension are(20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::interp_poli_at(double r1, double r2, double r3)
  {

    int dimension[] = {20};  // Xiaoyan changed from {20} to {8} for 8 nodes 07/12
    tensor h(1, dimension, 0.0);

      // influence of the node number 20
        h.val(20)=(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 19
        h.val(19)=(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 18
        h.val(18)=(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 17
        h.val(17)=(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;

    // influence of the node number 16
        h.val(16)=(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 15
        h.val(15)=(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    // influence of the node number 14
        h.val(14)=(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 13
        h.val(13)=(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;

    // influence of the node number 12
        h.val(12)=(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 11
        h.val(11)=(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    // influence of the node number 10
        h.val(10)=(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 9
        h.val( 9)=(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;

      // influence of the node number 8
    h.val(8)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (h.val(15)+h.val(16)+h.val(20))/2.0;
      // influence of the node number 7
    h.val(7)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0 - (h.val(14)+h.val(15)+h.val(19))/2.0;
      // influence of the node number 6
    h.val(6)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0 - (h.val(13)+h.val(14)+h.val(18))/2.0;
      // influence of the node number 5
    h.val(5)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0 - (h.val(13)+h.val(16)+h.val(17))/2.0;

      // influence of the node number 4
    h.val(4)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0 - (h.val(11)+h.val(12)+h.val(20))/2.0;
      // influence of the node number 3
    h.val(3)=(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0 - (h.val(10)+h.val(11)+h.val(19))/2.0;
      // influence of the node number 2
    h.val(2)=(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0 - (h.val(10)+h.val(18)+h.val(9))/2.0;
      // influence of the node number 1
    h.val(1)=(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0 - (h.val(12)+h.val(17)+h.val(9))/2.0;
    //    printf("r1 = %lf, r2 = %lf, r3 = %lf\n", r1, r2, r3);
    //    h.print("h");

    return h;
  }

//#############################################################################
//=========================================================================
// Derivarite of Shape functions  dimension are(20,3)     Wxy 09/26/2001
//=========================================================================

tensor TwentyNodeBrick_u_p_U::dh_drst_at(double r1, double r2, double r3)
  {

    int dimensions[] = {20,3};  // Changed from{20,3} to {8,3} Xiaoyan 07/12
    tensor dh(2, dimensions, 0.0);

    // influence of the node number 20
        dh.val(20,1) =   (1.0-r2)*(1.0-r3*r3)/4.0;
        dh.val(20,2) = - (1.0+r1)*(1.0-r3*r3)/4.0;
        dh.val(20,3) = - (1.0+r1)*(1.0-r2)*r3/2.0;
    // influence of the node number 19
        dh.val(19,1) = - (1.0-r2)*(1.0-r3*r3)/4.0;
        dh.val(19,2) = - (1.0-r1)*(1.0-r3*r3)/4.0;
        dh.val(19,3) = - (1.0-r1)*(1.0-r2)*r3/2.0;
    // influence of the node number 18
        dh.val(18,1) = - (1.0+r2)*(1.0-r3*r3)/4.0;
        dh.val(18,2) =   (1.0-r1)*(1.0-r3*r3)/4.0;
        dh.val(18,3) = - (1.0-r1)*(1.0+r2)*r3/2.0;
    // influence of the node number 17
        dh.val(17,1) =   (1.0+r2)*(1.0-r3*r3)/4.0;
        dh.val(17,2) =   (1.0+r1)*(1.0-r3*r3)/4.0;
        dh.val(17,3) = - (1.0+r1)*(1.0+r2)*r3/2.0;

    // influence of the node number 16
        dh.val(16,1) =   (1.0-r2*r2)*(1.0-r3)/4.0;
        dh.val(16,2) = - (1.0+r1)*r2*(1.0-r3)/2.0;
        dh.val(16,3) = - (1.0+r1)*(1.0-r2*r2)/4.0;
    // influnce of the node number 15
        dh.val(15,1) = - r1*(1.0-r2)*(1.0-r3)/2.0;
        dh.val(15,2) = - (1.0-r1*r1)*(1.0-r3)/4.0;
        dh.val(15,3) = - (1.0-r1*r1)*(1.0-r2)/4.0;
    // influence of the node number 14
        dh.val(14,1) = - (1.0-r2*r2)*(1.0-r3)/4.0;
        dh.val(14,2) = - (1.0-r1)*r2*(1.0-r3)/2.0;
        dh.val(14,3) = - (1.0-r1)*(1.0-r2*r2)/4.0;
    // influence of the node number 13
        dh.val(13,1) = - r1*(1.0+r2)*(1.0-r3)/2.0;
        dh.val(13,2) =   (1.0-r1*r1)*(1.0-r3)/4.0;
        dh.val(13,3) = - (1.0-r1*r1)*(1.0+r2)/4.0;

    // influence of the node number 12
        dh.val(12,1) =   (1.0-r2*r2)*(1.0+r3)/4.0;
        dh.val(12,2) = - (1.0+r1)*r2*(1.0+r3)/2.0;
        dh.val(12,3) =   (1.0+r1)*(1.0-r2*r2)/4.0;
    // influence of the node number 11
        dh.val(11,1) = - r1*(1.0-r2)*(1.0+r3)/2.0;
        dh.val(11,2) = - (1.0-r1*r1)*(1.0+r3)/4.0; // bug discovered 01 aug '95 2.0 -> 4.0
        dh.val(11,3) =   (1.0-r1*r1)*(1.0-r2)/4.0;
    // influence of the node number 10
        dh.val(10,1) = - (1.0-r2*r2)*(1.0+r3)/4.0;
        dh.val(10,2) = - (1.0-r1)*r2*(1.0+r3)/2.0;
        dh.val(10,3) =   (1.0-r1)*(1.0-r2*r2)/4.0;
    // influence of the node number 9
        dh.val(9,1)  = - r1*(1.0+r2)*(1.0+r3)/2.0;
        dh.val(9,2)  =   (1.0-r1*r1)*(1.0+r3)/4.0;
        dh.val(9,3)  =   (1.0-r1*r1)*(1.0+r2)/4.0;

      // influence of the node number 8
    dh.val(8,1)= (1.0-r2)*(1.0-r3)/8.0 - (dh.val(15,1)+dh.val(16,1)+dh.val(20,1))/2.0;
    dh.val(8,2)=-(1.0+r1)*(1.0-r3)/8.0 - (dh.val(15,2)+dh.val(16,2)+dh.val(20,2))/2.0;
    dh.val(8,3)=-(1.0+r1)*(1.0-r2)/8.0 - (dh.val(15,3)+dh.val(16,3)+dh.val(20,3))/2.0;
      // influence of the node number 7
    dh.val(7,1)=-(1.0-r2)*(1.0-r3)/8.0 - (dh.val(14,1)+dh.val(15,1)+dh.val(19,1))/2.0;
    dh.val(7,2)=-(1.0-r1)*(1.0-r3)/8.0 - (dh.val(14,2)+dh.val(15,2)+dh.val(19,2))/2.0;
    dh.val(7,3)=-(1.0-r1)*(1.0-r2)/8.0 - (dh.val(14,3)+dh.val(15,3)+dh.val(19,3))/2.0;
      // influence of the node number 6
    dh.val(6,1)=-(1.0+r2)*(1.0-r3)/8.0 - (dh.val(13,1)+dh.val(14,1)+dh.val(18,1))/2.0;
    dh.val(6,2)= (1.0-r1)*(1.0-r3)/8.0 - (dh.val(13,2)+dh.val(14,2)+dh.val(18,2))/2.0;
    dh.val(6,3)=-(1.0-r1)*(1.0+r2)/8.0 - (dh.val(13,3)+dh.val(14,3)+dh.val(18,3))/2.0;
      // influence of the node number 5
    dh.val(5,1)= (1.0+r2)*(1.0-r3)/8.0 - (dh.val(13,1)+dh.val(16,1)+dh.val(17,1))/2.0;
    dh.val(5,2)= (1.0+r1)*(1.0-r3)/8.0 - (dh.val(13,2)+dh.val(16,2)+dh.val(17,2))/2.0;
    dh.val(5,3)=-(1.0+r1)*(1.0+r2)/8.0 - (dh.val(13,3)+dh.val(16,3)+dh.val(17,3))/2.0;

      // influence of the node number 4
    dh.val(4,1)= (1.0-r2)*(1.0+r3)/8.0 - (dh.val(11,1)+dh.val(12,1)+dh.val(20,1))/2.0;
    dh.val(4,2)=-(1.0+r1)*(1.0+r3)/8.0 - (dh.val(11,2)+dh.val(12,2)+dh.val(20,2))/2.0;
    dh.val(4,3)= (1.0+r1)*(1.0-r2)/8.0 - (dh.val(11,3)+dh.val(12,3)+dh.val(20,3))/2.0;
      // influence of the node number 3
    dh.val(3,1)=-(1.0-r2)*(1.0+r3)/8.0 - (dh.val(10,1)+dh.val(11,1)+dh.val(19,1))/2.0;
    dh.val(3,2)=-(1.0-r1)*(1.0+r3)/8.0 - (dh.val(10,2)+dh.val(11,2)+dh.val(19,2))/2.0;
    dh.val(3,3)= (1.0-r1)*(1.0-r2)/8.0 - (dh.val(10,3)+dh.val(11,3)+dh.val(19,3))/2.0;
      // influence of the node number 2
    dh.val(2,1)=-(1.0+r2)*(1.0+r3)/8.0 - (dh.val(10,1)+dh.val(18,1)+dh.val(9,1))/2.0;
    dh.val(2,2)= (1.0-r1)*(1.0+r3)/8.0 - (dh.val(10,2)+dh.val(18,2)+dh.val(9,2))/2.0;
    dh.val(2,3)= (1.0-r1)*(1.0+r2)/8.0 - (dh.val(10,3)+dh.val(18,3)+dh.val(9,3))/2.0;
      // influence of the node number 1
    dh.val(1,1)= (1.0+r2)*(1.0+r3)/8.0 - (dh.val(12,1)+dh.val(17,1)+dh.val(9,1))/2.0;
    dh.val(1,2)= (1.0+r1)*(1.0+r3)/8.0 - (dh.val(12,2)+dh.val(17,2)+dh.val(9,2))/2.0;
    dh.val(1,3)= (1.0+r1)*(1.0+r2)/8.0 - (dh.val(12,3)+dh.val(17,3)+dh.val(9,3))/2.0;

    return dh;
  }

//CE Dynamic Allocation for brick3d


//=========================================================================
// Permeability tensor  dimension are(3,3)     Wxy 09/26/2001
//=========================================================================

// The permeability tensor has been moved in to the constructor 01/16/2002

// Xiaoyan added this function just want to test program. the values of k are not correct. 08/28/2001
/*tensor TwentyNodeBrick_u_p_U::k_at(double r1, double r2, double r3)
  {

    int k_dim[] = {3,3};
    tensor k(2, k_dim, 0.0);
    k.val(1,1)=r1;
    k.val(2,2)=r2;
    k.val(3,3)=r3;

    return k;
  }
*/
////#############################################################################
//Finite_Element * TwentyNodeBrick_u_p_U::new_el(int total)
//  {
//    TwentyNodeBrick_u_p_U *el_p;
//    el_p = new TwentyNodeBrick_u_p_U[total];
//    //DB//-------------------------------------------
//    //DB    for ( int i=0 ; i<total ; i++ )
//    //DB      {
//    //DB        el_p[i].report("derived TwentyNodeBrick_u_p_U\n");
//    //DB      }
//    //DB//-------------------------------------------
//    return el_p;
//  }

//=========================================================================
// Definition of Stiffness tensor Kep(20,3,3,8)     Wxy 09/26/2001
//=========================================================================

tensor TwentyNodeBrick_u_p_U::getStiffnessTensorKep()
  {
    int K_dim[] = {20,3,3,20};

    tensor Kep(4,K_dim,0.0);
    tensor Kkt(4,K_dim,0.0);


    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};
    tensor dh(2, dh_dim, 0.0);

    //    tensor Constitutive( 4, def_dim_4, 0.0);
    tensor Constitutive;

    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {20,3};
    tensor incremental_displacements(2,disp_dim,0.0); // \Delta u

    straintensor incremental_strain;
//    straintensor total_strain_at_GP;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
    tensor dhGlobal;


    short GP_c_r;
      for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        short GP_c_s;
          for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

          short GP_c_t;
          for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point fr04-Feb-2002om 3D array of short's
                where =
                   ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
    //dh.print("dh");
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                JacobianINVtemp = Jacobian.inverse();
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");  // changed form jk to kj 02/14/2002
                //        ::fprintf(stdout," # %d \n\n\n\n\n\n\n\n", El_count);
    //dhGlobal.print("dhGlobal");
                //weight
                weight = rw * sw * tw * det_of_Jacobian;
                // incremental straines at this Gauss point
                //GPstress[where].reportshortpqtheta("\n stress at GAUSS point in stiffness_tensor1\n");

                incremental_strain =
                     (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();

                Constitutive = (matpoint[where]->matmodel)->getTangentTensor();

                Kkt = dhGlobal("ib")*Constitutive("abcd");
                Kep = Kep + Kkt("aicd")*dhGlobal("jd")*weight;

    //Kk = Kk + dhGlobal("ib")*Constitutive("abcd")*dhGlobal("jd")*weight;
                //....K.print("K","\n\n K tensor \n");

    //Kmat = this->stiffness_matrix(Kk);
                //printf("K tensor max= %10.3e\n", Kmat.mmax());

                //convert constitutive and K to matrix and find min and max and print!



              }
          }
      }
    //Kk.print("K","\n\n K tensor \n");
    //K = Kk;
    return Kep;
  }

//=========================================================================
// Definition of Stiffness tensor G1(20,3,20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::getStiffnessTensorG1()  //(double rho_s, double n,)
  {
    //int M_dim[] = {20,3,3,20};
    int G_dim[] = {20,3,20};
    tensor G1(3,G_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dh(2, dh_dim, 0.0);

    int hp_dim[] = {20};
    tensor hp(1, hp_dim,0.0);     // Shape function.

    //int h_dim[] = {8,3};  // Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
//    tensor dH(2, h_dim, 0.0);
//    tensor Hp(1, Hp_dim,0.0);    // Xiaoyan added 08/27/2001

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
    tensor dhGlobal;

//    double RHO;
//    RHO= rho;    //global
//    double RHO_F=rho_f;
    double N=n;
    double ALPHA=alpha;

    short GP_c_r;
    for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

       short GP_c_s;
       for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

               short GP_c_t;
             for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
                hp= interp_poli_at(r,s,t);  // Assume the shape function of p (pressure) is the same as u.
                                // Xiaoyan 09/20/01
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //     Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //     printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates


                //weight
                weight = rw * sw * tw * det_of_Jacobian;

          //  tensor temp = H("ib")*H("kb");
    //temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

                G1 = G1 + dhGlobal("Ki")*hp("L")*weight * (ALPHA-N);
         //  printf("\n +++++++++++++++++++++++++ \n\n");
          //Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return G1;
  }

////#############################################################################
//=========================================================================
// Definition of Stiffness tensor G2(20,3,20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::getStiffnessTensorG2()  //(double rho_s, double n,)
  {
    //int M_dim[] = {20,3,3,20};
    int G_dim[] = {20,3,20};
    tensor G2(3,G_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dhU(2, dh_dim, 0.0);

    int hp_dim[] = {20};     // Xiaoyan changed from {60,3} to {24,3}
   tensor hp(1, hp_dim,0.0);     // Shape function.

    //int h_dim[] = {20,3};  // Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
//    tensor dH(2, h_dim, 0.0);
//    tensor Hp(1, Hp_dim,0.0);    // Xiaoyan added 08/27/2001

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
    tensor dhGlobal;

//    double RHO;
//    RHO= rho;    //global
//    double RHO_F=rho_f;
    double N=n;
//    double ALPHA=alpha;

    short GP_c_r;
    for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

          short GP_c_s;
          for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            short GP_c_t;
            for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dhU = dh_drst_at(r,s,t);  // Assume  the shape function of u and U are same. Xiaoyan 09/20/01
                hp= interp_poli_at(r,s,t); // Assume  the shape function of p is same as the u's. Xiaoyan 09/20/01
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dhU);
                //     Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dhU);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //     printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dhU("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates

                //weight
                weight = rw * sw * tw * det_of_Jacobian;
            //  printf("weight = %6.2e \n",weight);

    //M.print("M","BEFORE");

          //  tensor temp = H("ib")*H("kb");
    //temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

    G2 = G2 +  dhGlobal("Ki")*hp("L")*weight * N;
         //  printf("\n +++++++++++++++++++++++++ \n\n");
          //Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return G2;
  }
////#############################################################################
//=========================================================================
// Definition of Stiffness tensor P(20,20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::getStiffnessTensorP()
  {
    //int M_dim[] = {20,3,3,20};
    int G_dim[] = {20,20};
    tensor P(2,G_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dhU(2, dh_dim, 0.0);

    int hp_dim[] = {20};
    tensor hp(1, hp_dim,0.0);     // Shape function of pore presssure.

    //int h_dim[] = {8,3};  // Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
//    tensor dH(2, h_dim, 0.0);
//    tensor Hp(1, Hp_dim,0.0);    // Xiaoyan added 08/27/2001


    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
    tensor dhGlobal;

//    double RHO;
//    RHO= rho;    //global
//    double RHO_F=rho_f;
    double N=n;
    double ALPHA=alpha;
    double KS=ks;
    double KF=kf;
    double QQ= N/KF+(ALPHA-N)/KS;
//    double e;
//    double nu;
    short GP_c_r;
    for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

    short GP_c_s;
    for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

     short GP_c_t;
     for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dhU = dh_drst_at(r,s,t);  // Assume  the shape function of u and U are same. Xiaoyan 09/20/01
                hp= interp_poli_at(r,s,t); // Assume  the shape function of p is same as the u's. Xiaoyan 09/20/01
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dhU);
                //     Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dhU);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //     printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dhU("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates

                //weight
                weight = rw * sw * tw * det_of_Jacobian;
            //  printf("weight = %6.2e \n",weight);

    //M.print("M","BEFORE");

          //  tensor temp = H("ib")*H("kb");
    //temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

//    e = (matpoint[where]->matmodel)->getE();
//    cout<<"e=  "<<e<<endln;
//    nu = (matpoint[where]->matmodel)->getnu();
//    cout<<"nu=  "<<nu<<endln;
//    KS=e/(3*(1-2*nu));
//    cout<<"KS="<<KS<<endln;
//    QQ= N/KF+(ALPHA-N)/KS;
    P = P + hp("K") * hp("L")*QQ*weight;
         //  printf("\n +++++++++++++++++++++++++ \n\n");
          //Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return P;
  }

////#############################################################################
//=========================================================================
// Definition of Mass tensor Ms(20,3,3,20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::getMassTensorMs()  //(double rho_s, double n,)
  {
    int M_dim[] = {20,3,3,20};
    tensor Ms(4,M_dim,0.0);

    tensor I2("I", 2, def_dim_2);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {20};
    tensor H(1, h_dim, 0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

//    rho=(1-n)*rho_s+n*rho_f;
//    double RHO;
//    RHO= rho;    //global
    double RHO_S=rho_s;
    double N=n;

    short GP_c_r;
    for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

    short GP_c_s;
    for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

     short GP_c_t;
     for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //     Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //     printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                //                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates


                // printf("\n\nIN THE MASS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                // printf("  Mass_Tensor \n");
                // printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //                             GP_c_r,GP_c_s,GP_c_t);
                //
                H = interp_poli_at(r,s,t);

                //weight
                weight = rw * sw * tw * det_of_Jacobian;
            //  printf("weight = %6.2e \n",weight);

    //M.print("M","BEFORE");

          //  tensor temp = H("ib")*H("kb");
    //temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

            tensor temp = H("K") * I2("ij");
                  Ms = Ms + temp("Kij") * H("L") * ((1-N)*RHO_S *weight);
          //Ms.printshort("M");
              }
          }
      }
    //Ms.printshort("M");
    return Ms;
  }

////#############################################################################
//=========================================================================
// Definition of Mass tensor Mf(20,3,3,20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::getMassTensorMf()  //(double rho_s, double n,)
  {
    //int M_dim[] = {20,3,3,20};
    int M_dim[] = {20,3,3,20};
    tensor Mf(4,M_dim,0.0);

    tensor I2("I", 2, def_dim_2);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {20};
    tensor H(1, h_dim, 0.0);   // Shape fun. of u
    tensor HU(1, h_dim, 0.0);   // Shape fun. of U

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

//    rho=(1-n)*rho_s+n*rho_f;
//    double RHO;
//    RHO= rho;    //global
//    double RHO_S=rho_s;
    double RHO_F=rho_f;
    double N=n;


       short GP_c_r;
    for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

    short GP_c_s;
    for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

       short GP_c_t;
       for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //     Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //     printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                //                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates


                HU = interp_poli_at(r,s,t);   // now assume HU_3D(r,s,t) is the same as H_3D(r,s,t),
                                  // wxy 08/28/2001

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

            //  printf("weight = %6.2e \n",weight);

    //M.print("M","BEFORE");

          //  tensor temp = H("ib")*H("kb");
    //temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

    tensor temp = HU("K")* I2("ij");
    Mf = Mf + temp("Kij") * HU("L")* N * RHO_F * weight;
         //  printf("\n +++++++++++++++++++++++++ \n\n");
          //Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return Mf;
  }




//=========================================================================
// Definition of Damping tensor C1(20,3,3,20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::getDampTensorC1()  //(double rho_s, double n,)
  {
    int C_dim[] = {20,3,3,20};
    tensor C1(4,C_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {20};
    //int h_dim[] = {20,3};  // Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
    tensor H(1, h_dim, 0.0);
//    int k_dim[]={3,3};         // Xiaoyan added for permeability tensor 08/27/2001
//    tensor k(2,k_dim,0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

    double RHO_F=rho_f;
    double N=n;
    double G=9.8;

     short GP_c_r;
     for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        short GP_c_s;
        for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

        short GP_c_t;
        for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
//    k=k_at(r,s,t);    // wxy added for permeability.
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //     Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //     printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                //                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates


                // printf("\n\nIN THE MASS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                // printf("  Mass_Tensor \n");
                // printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //                             GP_c_r,GP_c_s,GP_c_t);
                //
                H = interp_poli_at(r,s,t);

                //weight

                weight = rw * sw * tw * det_of_Jacobian;

            //  printf("weight = %6.2e \n",weight);


          //  tensor temp = H("ib")*H("kb");
    //temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

    tensor k_inverse=k("mn").inverse();
    tensor temp=H("K")* k_inverse("ij");

    C1 = C1 + temp("Kij")*H("L")*weight * N*N*RHO_F*G;
         //  printf("\n +++++++++++++++++++++++++ \n\n");
          //Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return C1;
  }


//=========================================================================
// Definition of Damping tensor C2(20,3,3,20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::getDampTensorC2()  //(double rho_s, double n,)
  {
    //int M_dim[] = {20,3,3,20};
    int C_dim[] = {20,3,3,20};
    tensor C2(4,C_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {20};
    //int h_dim[] = {20,3};  // Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
    tensor H(1, h_dim, 0.0);
    tensor HU(1, h_dim, 0.0);
//   int k_dim[]={3,3};         // Xiaoyan added for permeability tensor 08/27/2001
//    tensor k(2,k_dim,0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

    double RHO_F=rho_f;
    double N=n;
    double G=9.8;


     short GP_c_r;
     for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

       short GP_c_s;
       for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

           short GP_c_t;
           for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
//    k=k_at(r,s,t);    // wxy added for permeability.
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //     Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //     printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                //                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates


                // printf("\n\nIN THE MASS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                // printf("  Mass_Tensor \n");
                // printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //                             GP_c_r,GP_c_s,GP_c_t);
                //
                H = interp_poli_at(r,s,t);
                HU = interp_poli_at(r,s,t);  // assume HU=H now . 08/28/2001


                //weight
                weight = rw * sw * tw * det_of_Jacobian;
            //  printf("weight = %6.2e \n",weight);

    //M.print("M","BEFORE");

          //  tensor temp = H("ib")*H("kb");
    //temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

    tensor k_inverse=k("mn").inverse();
    tensor temp=H("L")* k_inverse("ij");
    C2 = C2 + temp("Lij")*HU("K")*weight * N*N*RHO_F*G;
         //  printf("\n +++++++++++++++++++++++++ \n\n");
          //Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return C2;
  }


//=========================================================================
// Definition of Damping tensor C3(20,3,3,20)     Wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::getDampTensorC3()  //(double rho_s, double n,)
  {
    //int M_dim[] = {20,3,3,20};
    int C_dim[] = {20,3,3,20};
    tensor C3(4,C_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {20};
    //int h_dim[] = {20,3};  // Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
    tensor HU(1, h_dim, 0.0);
//    int k_dim[]={3,3};         // Xiaoyan added for permeability tensor 08/27/2001
//    tensor k(2,k_dim,0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

    double RHO_F=rho_f;
    double N=n;
    double G=9.8;

     short GP_c_r;
     for( GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        short GP_c_s;
        for( GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            short GP_c_t;
             for( GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
//    k=k_at(r,s,t);    // wxy added for permeability.
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //     Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //     printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                //                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates


                // printf("\n\nIN THE MASS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                // printf("  Mass_Tensor \n");
                // printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //                             GP_c_r,GP_c_s,GP_c_t);
                //
                // H = interp_poli_at(r,s,t);
                HU = interp_poli_at(r,s,t);  // assume HU=H now . 08/28/2001


                //weight
                weight = rw * sw * tw * det_of_Jacobian;
            //  printf("weight = %6.2e \n",weight);

    //M.print("M","BEFORE");

          //  tensor temp = H("ib")*H("kb");
    //temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

    tensor k_inverse=k("mn").inverse();
    tensor temp= HU("L")* k_inverse("ij");
    C3 = C3 + temp("Lij") *HU("K")*weight * N*N*RHO_F*G;
         //  printf("\n +++++++++++++++++++++++++ \n\n");
          //Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return C3;
  }

////#######################################################################
//=========================================================================
// Converting stiffness tensor to stiffness matrix K^ep 60*60
//=========================================================================

matrix TwentyNodeBrick_u_p_U::stiffness_matrixKep(const tensor  Kep)

  {
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    static matrix    Kepmatrix(60,60,0.0);

    int Ki=0;
    int Kj=0;

    int i;
    for ( i=1 ; i<=20 ; i++ )  // i<=20 for 20 nodes
      {
       int j;
       for ( j=1 ; j<=20 ; j++ )  // i<=20 for 20 nodes
          {
            int k;
            for ( k=1 ; k<=3 ; k++ )
              {
                int l;
        for ( l=1 ; l<=3 ; l++ )
                  {
                    Ki = k+3*(i-1);
                    Kj = l+3*(j-1);
                    //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                    Kepmatrix.val( Ki , Kj ) = Kep.cval(i,k,l,j);
                  }
              }
          }
      }
    return Kepmatrix;
  }

//=========================================================================
// Converting stiffness tensor to stiffness matrix G1 60*20
//=========================================================================

matrix TwentyNodeBrick_u_p_U::stiffness_matrixG1(const tensor  G1)

  {
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    static matrix    G1matrix(60,20,0.0);

    int Ki=0;
    int Kj=0;

    int i;
    for ( i=1 ; i<=20 ; i++ )
      {
        int k;
       for ( k=1 ; k<=3 ; k++ )
          {
            int l;
              for ( l=1 ; l<=20 ; l++ )
              {
                Ki = k+3*(i-1);
                Kj = l;
                //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                G1matrix.val( Ki , Kj ) = G1.cval(i,k,l);
              }
          }
      }
    return G1matrix;
  }

//=========================================================================
// Converting stiffness tensor to stiffness matrix G2 60*20
//=========================================================================


matrix TwentyNodeBrick_u_p_U::stiffness_matrixG2(const tensor  G2)

  {
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    static matrix    G2matrix(60,20,0.0);

    int Ki=0;
    int Kj=0;

    int i;
     for ( i=1 ; i<=20 ; i++ )  // i<=20 for 20 nodes
      {
        int k;
        for ( k=1 ; k<=3 ; k++ )
          {
            int l;
             for ( l=1 ; l<=20 ; l++ )
              {
                Ki = k+3*(i-1);
                Kj = l;
                //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                G2matrix.val( Ki , Kj ) = G2.cval(i,k,l);
              }
          }
      }
    return G2matrix;
  }


//=========================================================================
// Converting stiffness tensor to stiffness matrix P 20*20
//=========================================================================

matrix TwentyNodeBrick_u_p_U::stiffness_matrixP(const tensor P)

  {
    static matrix    Pmatrix(20,20,0.0);

    int i;
     for ( i=1 ; i<=20 ; i++ )  // i<=20 for 20 nodes
      {
        int j;
       for ( j=1 ; j<=20 ; j++ )
          {
             Pmatrix.val( i , j ) = P.cval(i,j);

          }
      }
    return Pmatrix;
  }

//=========================================================================
// Constructing the whole system K(140*140) (including Kep G1 and G2 56*56)   Wxy 09/26/2001
//=========================================================================

const Matrix &TwentyNodeBrick_u_p_U::getTangentStiff()

  {
    tensor tKep = getStiffnessTensorKep();
    tensor tG1  = getStiffnessTensorG1();
    tensor tG2  = getStiffnessTensorG2();
    tensor tP   = getStiffnessTensorP();
    matrix Kep = stiffness_matrixKep(tKep);
    matrix G1  = stiffness_matrixG1(tG1);
    matrix G2  = stiffness_matrixG2(tG2);
    matrix P   = stiffness_matrixP(tP);

//    Kep.write_standard("Kep_20upu.dat", "stiffness part of KupU");
//    G1.write_standard("G1_20upu.dat", "stiffness part of G1upU");
//    G2.write_standard("G2_20upu.dat", "stiffness part of G2upU");
//    P.write_standard("P_20upu.dat", "stiffness part of PupU");
    matrix G1t=G1.transpose();
    matrix G2t=G2.transpose();
//    G1t.write_standard("G1t.dat", "stiffness part of G1upU");
//    G2t.write_standard("G2t.dat", "stiffness part of G1upU");

    int i;
    for ( i=0 ; i<20 ; i++ )
     {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
            int n;
            for( n=0; n<3; n++)
          {
             int m;
             for(m=0; m<3; m++)
                    {
                       K(i*7+n ,  j*7+m)=Kep.val(3*i+n+1, 3*j+1+m);     // Add Kep to K
                     }
                }
           }
      }

 for ( i=0 ; i<20 ; i++ )
      {
        int n;
  for( n=0; n<3; n++)
          {
            int j;
            for ( j=0 ; j<20 ; j++ )
               {
            K(i*7+n ,  j*7+3)=-G1.val(3*i+n+1, j+1);     // Add -G1 to K
               }
          }
      }

 for ( i=0 ; i<20 ; i++ )
      {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
           int n;
           for( n=0; n<3; n++)
           {
                K(i*7+3 ,  j*7+n)=-G1t.val(i+1, 3*j+1+n);     // Add -G1T to K
        }
       }
      }


 for ( i=0 ; i<20 ; i++ )
      {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
               K(i*7+3 ,  j*7+3)=-P.val(i+1, j+1);     // Add -P to K
           }
      }


 for ( i=0 ; i<20 ; i++ )
      {
       int j;
       for ( j=0 ; j<20 ; j++ )
          {
             int n;
             for ( n=0; n<3; n++)
             {
                 K(i*7+3 ,  j*7+n+4)=-G2t.val(i+1, 3*j+1+n);     // Add -G2T to K
                }
           }
      }

 for ( i=0 ; i<20 ; i++ )
      {
        int n;
        for ( n=0; n<3; n++)
          {
            int j;
            for ( j=0 ; j<20 ; j++ )
        {
                 K(i*7+n+4 ,  j*7+3)=-G2.val(3*i+n+1, j+1);     // Add -G2 to K
        }
           }
      }

//     ofstream outK("K20.dat");  // want to check the whole K
//     K.Output(outK);

     return K;
  }


//=========================================================================
// Converting damping tensor to damping matrix C1(60*60)
//=========================================================================

matrix TwentyNodeBrick_u_p_U::damping_matrixC1(const tensor  C1)

  {
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    static matrix    C1matrix(60,60,0.0);

    int Ki=0;
    int Kj=0;

    int i;
    for ( i=1 ; i<=20 ; i++ )  //  i<=20 for 20 nodes
      {
        int j;
        for ( j=1 ; j<=20 ; j++ )  //  i<=20 for 20 nodes
          {
            int k;
            for ( k=1 ; k<=3 ; k++ )
              {
                int l;
          for ( l=1 ; l<=3 ; l++ )
                  {
                    Ki = k+3*(i-1);
                    Kj = l+3*(j-1);
                    //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                    C1matrix.val( Ki , Kj ) = C1.cval(i,k,l,j);
                  }
              }
          }
      }
    return C1matrix;
  }

////#######################################################################
//=========================================================================
// Converting damping tensor to damping matrix  C2(60*60)
//=========================================================================

matrix TwentyNodeBrick_u_p_U::damping_matrixC2(const tensor  C2)
  {
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    static matrix    C2matrix(60,60,0.0);

    int Ki=0;
    int Kj=0;

    int i;
    for ( i=1 ; i<=20 ; i++ )  //  i<=20 for 20 nodes
      {
        int j;
        for ( j=1 ; j<=20 ; j++ )  //  i<=20 for 20 nodes
          {
            int k;
            for ( k=1 ; k<=3 ; k++ )
              {
                int l;
          for ( l=1 ; l<=3 ; l++ )
                  {
                    Ki = k+3*(i-1);
                    Kj = l+3*(j-1);
                    //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                    C2matrix.val( Ki , Kj ) = C2.cval(i,k,l,j);
                  }
              }
          }
      }
    return C2matrix;
  }

//=========================================================================
// Converting damping tensor to damping matrix C3(60*60)
//=========================================================================

matrix TwentyNodeBrick_u_p_U::damping_matrixC3(const tensor  C3)
  {
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    static matrix    C3matrix(60,60,0.0);

    int Ki=0;
    int Kj=0;

    int i;
    for ( i=1 ; i<=20 ; i++ )  //  i<=20 for 20 nodes
      {
        int j;
        for ( j=1 ; j<=20 ; j++ )  //  i<=20 for 20 nodes
          {
            int k;
            for ( k=1 ; k<=3 ; k++ )
              {
                int l;
    for ( l=1 ; l<=3 ; l++ )
                  {
                    Ki = k+3*(i-1);
                    Kj = l+3*(j-1);
                    //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                    C3matrix.val( Ki , Kj ) = C3.cval(i,k,l,j);
                  }
              }
          }
      }
    return C3matrix;
  }

//=========================================================================
// Constructing Damping matrix of the whole system C(140*140) (including C1 C2 and C3)
//=========================================================================

const Matrix &TwentyNodeBrick_u_p_U::getDamp()

  {

    tensor tC1  = getDampTensorC1();
    tensor tC2  = getDampTensorC2();
    tensor tC3  = getDampTensorC3();
    matrix C1 = damping_matrixC1(tC1);
    matrix C2 = damping_matrixC2(tC2);
    matrix C3 = damping_matrixC3(tC3);
//    C1.write_standard("C1.dat", "Damping matrix");
//    C2.write_standard("C2.dat", "Damping matrix");
//    C3.write_standard("C3.dat", "Damping matrix");

    matrix C2t=C2.transpose();
//    C2t.write_standard("C2t.dat", "Damping matrix");

    int i;
    for ( i=0 ; i<20 ; i++ )
      {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
             int n;
             for ( n=0; n<3; n++)
              {
                 int m;
                 for ( m=0; m<3; m++)
                     {
                       C(i*7+n ,  j*7+m)=C1.val(3*i+n+1, 3*j+1+m);     // Add C1 to C
                           }
                     }
       }
      }

 for ( i=0 ; i<20 ; i++ )
      {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
              int n;
              for ( n=0; n<3; n++)
           {
                   int m;
             for ( m=0; m<3; m++)
                {
                        C(i*7+n ,  j*7+m+4)=-C2.val(3*i+n+1, 3*j+1+m);     // Add -C2 to C
                       }
                  }
           }
      }

 for ( i=0 ; i<20 ; i++ )
      {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
      int n;
      for ( n=0; n<3; n++)
               {
                 int m;
          for ( m=0; m<3; m++)
                   {
                     C(i*7+n+4 ,  j*7+m)=-C2t.val(3*i+n+1, 3*j+1+m);     // Add -C2t to C
                   }
              }
           }
      }

 for ( i=0 ; i<20 ; i++ )
      {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
            int n;
            for ( n=0; n<3; n++)
              {
                 int m;
          for ( m=0; m<3; m++)
                  {
                     C(i*7+n+4,  j*7+m+4)=C3.val(3*i+n+1, 3*j+1+m);     // Add C3 to C
                  }
               }
           }
      }

    return C;
  }

//=========================================================================
// Converting mass tensor to mass matrix Ms(60,60) ___Xiaoyan 08/27/2001
//=========================================================================
matrix TwentyNodeBrick_u_p_U::mass_matrixMs(const tensor  Ms)
  {
    static matrix    Msmatrix(60,60,0.0);

    int Ki=0;
    int Kj=0;

    int i;
    for ( i=1 ; i<=20 ; i++ )  //  i<=20 for 20 nodes
      {
        int j;
        for ( j=1 ; j<=20 ; j++ )  // i<=20 for 20 nodes
          {
            int k;
            for ( k=1 ; k<=3 ; k++ )
              {
                int l;
          for ( l=1 ; l<=3 ; l++ )
                  {
                    Ki = k+3*(i-1);
                    Kj = l+3*(j-1);
                    //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                    Msmatrix.val( Ki , Kj ) = Ms.cval(i,k,l,j);
                  }
              }
          }
      }
    return Msmatrix;
}

//=========================================================================
// Converting mass tensor to mass matrix Mf(60,60)___Xiaoyan 08/27/2001
//=========================================================================

matrix TwentyNodeBrick_u_p_U::mass_matrixMf(const tensor  Mf)
  {
    static matrix    Mfmatrix(60,60,0.0);

    int Ki=0;
    int Kj=0;

    int i;
    for ( i=1 ; i<=20 ; i++ )  // i<=20 for 20 nodes
      {
        int j;
        for ( j=1 ; j<=20 ; j++ )  //  i<=20 for 20 nodes
          {
            int k;
            for ( k=1 ; k<=3 ; k++ )
              {
                int l;
          for ( l=1 ; l<=3 ; l++ )
                  {
                    Ki = k+3*(i-1);
                    Kj = l+3*(j-1);
                    //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                    Mfmatrix.val( Ki , Kj ) = Mf.cval(i,k,l,j);
                  }
              }
          }
      }
    return Mfmatrix;
}
////#############################################################################

//#############################################################################
//=========================================================================
// Constructing Mass matrix of the whole system M(140,140) (including Ms and Mf)
//=========================================================================

const Matrix &TwentyNodeBrick_u_p_U::getMass()
 {
    tensor tMs  = getMassTensorMs();
    tensor tMf  = getMassTensorMf();
    matrix  Ms =  mass_matrixMs(tMs);
    matrix  Mf =  mass_matrixMf(tMf);
//    Ms.write_standard("Ms.dat", "Mass matrix");
//    Mf.write_standard("Mf.dat", "Mass matrix");

    int i;
    for ( i=0 ; i<20 ; i++ )
       {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
             int n;
             for ( n=0; n<3; n++)
               {
                 int m;
          for ( m=0; m<3; m++)
                  {
                     M(i*7+n ,  j*7+m)=Ms.val(3*i+n+1, 3*j+1+m);     // Add Ms to M
                  }
               }
           }
      }

 for ( i=0 ; i<20 ; i++ )
      {
        int j;
        for ( j=0 ; j<20 ; j++ )
          {
             int n;
             for ( n=0; n<3; n++)
               {
                 int m;
          for ( m=0; m<3; m++)
                   {
                     M(i*7+n+4 ,  j*7+m+4)=Mf.val(3*i+n+1, 3*j+1+m);     // Add Mf to M
                  }
               }
           }
      }

//     ofstream outM("M20.dat");  // want to check the whole M
//     M.Output(outM);

     return M;
  }

//#############################################################################


//=========================================================================
// Jacobian tensor J = dx/dr=dN/dr*x_i=dh*coordinate  wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::Jacobian_3D(tensor dh)
  {
                   // dh(20*3) Xiaoyan
     tensor N_C = Nodal_Coordinates();
     tensor Jacobian_3D = dh("ij") * N_C("ik");
     //         (3*3)      (20*3)      (20*3)
     return Jacobian_3D;
  }

//#############################################################################
//=========================================================================
// Jacobian inverse J^(-1)  wxy 09/26/2001
//=========================================================================
tensor TwentyNodeBrick_u_p_U::Jacobian_3Dinv(tensor dh)
  {
     //       dh ( 20*3)    // dh(20*3) Xiaoyan
     tensor N_C = Nodal_Coordinates();         // 20*3 Xiaoyan
     tensor Jacobian_3Dinv = (dh("ij") * N_C("ik")).inverse();
     return Jacobian_3Dinv;
  }


////#############################################################################
//=========================================================================
tensor TwentyNodeBrick_u_p_U::Nodal_Coordinates()
  {
    const int dimensions[] = {20,3};
    tensor N_coord(2, dimensions, 0.0);

    const Vector &nd1Crds = nd1Ptr->getCrds();
    const Vector &nd2Crds = nd2Ptr->getCrds();
    const Vector &nd3Crds = nd3Ptr->getCrds();
    const Vector &nd4Crds = nd4Ptr->getCrds();

    const Vector &nd5Crds = nd5Ptr->getCrds();
    const Vector &nd6Crds = nd6Ptr->getCrds();
    const Vector &nd7Crds = nd7Ptr->getCrds();
    const Vector &nd8Crds = nd8Ptr->getCrds();

    const Vector &nd9Crds  =  nd9Ptr->getCrds();
    const Vector &nd10Crds = nd10Ptr->getCrds();
    const Vector &nd11Crds = nd11Ptr->getCrds();
    const Vector &nd12Crds = nd12Ptr->getCrds();

    const Vector &nd13Crds = nd13Ptr->getCrds();
    const Vector &nd14Crds = nd14Ptr->getCrds();
    const Vector &nd15Crds = nd15Ptr->getCrds();
    const Vector &nd16Crds = nd16Ptr->getCrds();

    const Vector &nd17Crds = nd17Ptr->getCrds();
    const Vector &nd18Crds = nd18Ptr->getCrds();
    const Vector &nd19Crds = nd19Ptr->getCrds();
    const Vector &nd20Crds = nd20Ptr->getCrds();

    N_coord.val(1,1)=nd1Crds(0); N_coord.val(1,2)=nd1Crds(1); N_coord.val(1,3)=nd1Crds(2);
    N_coord.val(2,1)=nd2Crds(0); N_coord.val(2,2)=nd2Crds(1); N_coord.val(2,3)=nd2Crds(2);
    N_coord.val(3,1)=nd3Crds(0); N_coord.val(3,2)=nd3Crds(1); N_coord.val(3,3)=nd3Crds(2);
    N_coord.val(4,1)=nd4Crds(0); N_coord.val(4,2)=nd4Crds(1); N_coord.val(4,3)=nd4Crds(2);

    N_coord.val(5,1)=nd5Crds(0); N_coord.val(5,2)=nd5Crds(1); N_coord.val(5,3)=nd5Crds(2);
    N_coord.val(6,1)=nd6Crds(0); N_coord.val(6,2)=nd6Crds(1); N_coord.val(6,3)=nd6Crds(2);
    N_coord.val(7,1)=nd7Crds(0); N_coord.val(7,2)=nd7Crds(1); N_coord.val(7,3)=nd7Crds(2);
    N_coord.val(8,1)=nd8Crds(0); N_coord.val(8,2)=nd8Crds(1); N_coord.val(8,3)=nd8Crds(2);

    N_coord.val(9 ,1)=nd9Crds(0);  N_coord.val(9 ,2)=nd9Crds(1);  N_coord.val(9 ,3)=nd9Crds(2);
    N_coord.val(10,1)=nd10Crds(0); N_coord.val(10,2)=nd10Crds(1); N_coord.val(10,3)=nd10Crds(2);
    N_coord.val(11,1)=nd11Crds(0); N_coord.val(11,2)=nd11Crds(1); N_coord.val(11,3)=nd11Crds(2);
    N_coord.val(12,1)=nd12Crds(0); N_coord.val(12,2)=nd12Crds(1); N_coord.val(12,3)=nd12Crds(2);

    N_coord.val(13,1)=nd13Crds(0); N_coord.val(13,2)=nd13Crds(1); N_coord.val(13,3)=nd13Crds(2);
    N_coord.val(14,1)=nd14Crds(0); N_coord.val(14,2)=nd14Crds(1); N_coord.val(14,3)=nd14Crds(2);
    N_coord.val(15,1)=nd15Crds(0); N_coord.val(15,2)=nd15Crds(1); N_coord.val(15,3)=nd15Crds(2);
    N_coord.val(16,1)=nd16Crds(0); N_coord.val(16,2)=nd16Crds(1); N_coord.val(16,3)=nd16Crds(2);

    N_coord.val(17,1)=nd17Crds(0); N_coord.val(17,2)=nd17Crds(1); N_coord.val(17,3)=nd17Crds(2);
    N_coord.val(18,1)=nd18Crds(0); N_coord.val(18,2)=nd18Crds(1); N_coord.val(18,3)=nd18Crds(2);
    N_coord.val(19,1)=nd19Crds(0); N_coord.val(19,2)=nd19Crds(1); N_coord.val(19,3)=nd19Crds(2);
    N_coord.val(20,1)=nd20Crds(0); N_coord.val(20,2)=nd20Crds(1); N_coord.val(20,3)=nd20Crds(2);

    return N_coord;
  }

////#############################################################################
//=========================================================================
// Incremental displacement of solid part--u     02/10/2002
//=========================================================================
tensor TwentyNodeBrick_u_p_U::incr_dispDu()
  {
    const int dimensions[] = {20,3};
    tensor increment_dispDu(2, dimensions, 0.0);

    //for ( int i=0 ; i<8 ; i++ )   // Xiaoyan changed from 20 to 8 for 8 nodes
    //
    //  {
    //    // increment_disp.val(i+1,1) = nodes[ G_N_numbs[i] ].incremental_translation_x();
    //    // increment_disp.val(i+1,2) = nodes[ G_N_numbs[i] ].incremental_translation_y();
    //    // increment_disp.val(i+1,3) = nodes[ G_N_numbs[i] ].incremental_translation_z();
    //    // Xiaoyan changed to the following 09/27/00
    //    Vector IncremenDis = nodes[ G_N_numbs[i] ].getIncrDisp();
    //
    //    increment_disp.val(i+1,1) = IncremenDis(0);
    //    increment_disp.val(i+1,2) = IncremenDis(1);
    //    increment_disp.val(i+1,3) = IncremenDis(2);
    //
    //  }

    //Zhaohui using node pointers, which come from the Domain
    //const Vector &TotDis1 = nd1Ptr->getTrialDisp();
    //const Vector &incrdelDis1 = nd1Ptr->getIncrDisp();
    //Have to get IncrDeltaDisp, not IncrDisp for cumulation of incr_disp
    const Vector &IncrDis1 = nd1Ptr->getIncrDeltaDisp();
    const Vector &IncrDis2 = nd2Ptr->getIncrDeltaDisp();
    const Vector &IncrDis3 = nd3Ptr->getIncrDeltaDisp();
    const Vector &IncrDis4 = nd4Ptr->getIncrDeltaDisp();

    const Vector &IncrDis5 = nd5Ptr->getIncrDeltaDisp();
    const Vector &IncrDis6 = nd6Ptr->getIncrDeltaDisp();
    const Vector &IncrDis7 = nd7Ptr->getIncrDeltaDisp();
    const Vector &IncrDis8 = nd8Ptr->getIncrDeltaDisp();

    const Vector &IncrDis9  = nd9Ptr->getIncrDeltaDisp();
    const Vector &IncrDis10 = nd10Ptr->getIncrDeltaDisp();
    const Vector &IncrDis11 = nd11Ptr->getIncrDeltaDisp();
    const Vector &IncrDis12 = nd12Ptr->getIncrDeltaDisp();

    const Vector &IncrDis13 = nd13Ptr->getIncrDeltaDisp();
    const Vector &IncrDis14 = nd14Ptr->getIncrDeltaDisp();
    const Vector &IncrDis15 = nd15Ptr->getIncrDeltaDisp();
    const Vector &IncrDis16 = nd16Ptr->getIncrDeltaDisp();

    const Vector &IncrDis17 = nd17Ptr->getIncrDeltaDisp();
    const Vector &IncrDis18 = nd18Ptr->getIncrDeltaDisp();
    const Vector &IncrDis19 = nd19Ptr->getIncrDeltaDisp();
    const Vector &IncrDis20 = nd20Ptr->getIncrDeltaDisp();

// Get the first three incremental displacement for solid part. Xiaoyan 02/08/2002

    increment_dispDu.val(1,1)=IncrDis1(0);
    increment_dispDu.val(1,2)=IncrDis1(1);
    increment_dispDu.val(1,3)=IncrDis1(2);
//    increment_disp.val(1,4)=IncrDis1(3); increment_disp.val(1,5)=IncrDis1(4);increment_disp.val(1,6)=IncrDis1(5);
//    increment_disp.val(1,7)=IncrDis1(6);

    increment_dispDu.val(2,1)=IncrDis2(0);
    increment_dispDu.val(2,2)=IncrDis2(1);
    increment_dispDu.val(2,3)=IncrDis2(2);
//    increment_disp.val(2,4)=IncrDis2(3); increment_disp.val(2,5)=IncrDis2(4);increment_disp.val(2,6)=IncrDis2(5);
//    increment_disp.val(2,7)=IncrDis2(6);

    increment_dispDu.val(3,1)=IncrDis3(0);
    increment_dispDu.val(3,2)=IncrDis3(1);
    increment_dispDu.val(3,3)=IncrDis3(2);
//    increment_disp.val(3,4)=IncrDis3(3); increment_disp.val(3,5)=IncrDis3(4);increment_disp.val(3,6)=IncrDis3(5);
//    increment_disp.val(3,7)=IncrDis3(6);

    increment_dispDu.val(4,1)=IncrDis4(0);
    increment_dispDu.val(4,2)=IncrDis4(1);
    increment_dispDu.val(4,3)=IncrDis4(2);
//    increment_disp.val(4,4)=IncrDis4(3); increment_disp.val(4,5)=IncrDis4(4);increment_disp.val(4,6)=IncrDis4(5);
//    increment_disp.val(4,7)=IncrDis4(6);

    increment_dispDu.val(5,1)=IncrDis5(0);
    increment_dispDu.val(5,2)=IncrDis5(1);
    increment_dispDu.val(5,3)=IncrDis5(2);
//    increment_disp.val(5,4)=IncrDis5(3); increment_disp.val(5,5)=IncrDis5(4);increment_disp.val(5,6)=IncrDis5(5);
//    increment_disp.val(5,7)=IncrDis5(6);

    increment_dispDu.val(6,1)=IncrDis6(0);
    increment_dispDu.val(6,2)=IncrDis6(1);
    increment_dispDu.val(6,3)=IncrDis6(2);
//    increment_disp.val(6,4)=IncrDis6(3); increment_disp.val(6,5)=IncrDis6(4);increment_disp.val(4,6)=IncrDis6(5);
//    increment_disp.val(6,7)=IncrDis6(6);

    increment_dispDu.val(7,1)=IncrDis7(0);
    increment_dispDu.val(7,2)=IncrDis7(1);
    increment_dispDu.val(7,3)=IncrDis7(2);
//    increment_disp.val(7,4)=IncrDis7(3); increment_disp.val(7,5)=IncrDis7(4);increment_disp.val(7,6)=IncrDis7(5);
//    increment_disp.val(7,7)=IncrDis7(6);

    increment_dispDu.val(8,1)=IncrDis8(0);
    increment_dispDu.val(8,2)=IncrDis8(1);
    increment_dispDu.val(8,3)=IncrDis8(2);
//    increment_disp.val(8,4)=IncrDis8(3); increment_disp.val(8,5)=IncrDis8(4);increment_disp.val(8,6)=IncrDis8(5);
//    increment_disp.val(8,7)=IncrDis8(6);

    increment_dispDu.val(9,1)=IncrDis9(0);
    increment_dispDu.val(9,2)=IncrDis9(1);
    increment_dispDu.val(9,3)=IncrDis9(2);
//    increment_disp.val(9,4)=IncrDis9(3); increment_disp.val(9,5)=IncrDis9(4);increment_disp.val(9,6)=IncrDis9(5);
//    increment_disp.val(9,7)=IncrDis9(6);

    increment_dispDu.val(10,1)=IncrDis10(0);
    increment_dispDu.val(10,2)=IncrDis10(1);
    increment_dispDu.val(10,3)=IncrDis10(2);
//    increment_disp.val(10,4)=IncrDis10(3); increment_disp.val(10,5)=IncrDis10(4);increment_disp.val(10,6)=IncrDis10(5);
//    increment_disp.val(10,7)=IncrDis10(6);

    increment_dispDu.val(11,1)=IncrDis11(0);
    increment_dispDu.val(11,2)=IncrDis11(1);
    increment_dispDu.val(11,3)=IncrDis11(2);
//    increment_disp.val(11,4)=IncrDis11(3); increment_disp.val(11,5)=IncrDis11(4);increment_disp.val(11,6)=IncrDis11(5);
//    increment_disp.val(11,7)=IncrDis11(6);

    increment_dispDu.val(12,1)=IncrDis12(0);
    increment_dispDu.val(12,2)=IncrDis12(1);
    increment_dispDu.val(12,3)=IncrDis12(2);
//    increment_disp.val(12,4)=IncrDis12(3); increment_disp.val(12,5)=IncrDis12(4);increment_disp.val(12,6)=IncrDis12(5);
//    increment_disp.val(12,7)=IncrDis12(6);

    increment_dispDu.val(13,1)=IncrDis13(0);
    increment_dispDu.val(13,2)=IncrDis13(1);
    increment_dispDu.val(13,3)=IncrDis13(2);
//    increment_disp.val(13,4)=IncrDis13(3); increment_disp.val(13,5)=IncrDis13(4);increment_disp.val(13,6)=IncrDis13(5);
//    increment_disp.val(13,7)=IncrDis13(6);

    increment_dispDu.val(14,1)=IncrDis14(0);
    increment_dispDu.val(14,2)=IncrDis14(1);
    increment_dispDu.val(14,3)=IncrDis14(2);
//    increment_disp.val(14,4)=IncrDis14(3); increment_disp.val(14,5)=IncrDis14(4);increment_disp.val(14,6)=IncrDis14(5);
//    increment_disp.val(14,7)=IncrDis14(6);

    increment_dispDu.val(15,1)=IncrDis15(0);
    increment_dispDu.val(15,2)=IncrDis15(1);
    increment_dispDu.val(15,3)=IncrDis15(2);
//    increment_disp.val(15,4)=IncrDis15(3); increment_disp.val(15,5)=IncrDis15(4);increment_disp.val(15,6)=IncrDis15(5);
//    increment_disp.val(15,7)=IncrDis15(6);

    increment_dispDu.val(16,1)=IncrDis16(0);
    increment_dispDu.val(16,2)=IncrDis16(1);
    increment_dispDu.val(16,3)=IncrDis16(2);
//    increment_disp.val(16,4)=IncrDis16(3); increment_disp.val(16,5)=IncrDis16(4);increment_disp.val(16,6)=IncrDis16(5);
//    increment_disp.val(16,7)=IncrDis16(6);

    increment_dispDu.val(17,1)=IncrDis17(0);
    increment_dispDu.val(17,2)=IncrDis17(1);
    increment_dispDu.val(17,3)=IncrDis17(2);
//    increment_disp.val(17,4)=IncrDis17(3); increment_disp.val(17,5)=IncrDis17(4);increment_disp.val(17,6)=IncrDis17(5);
//    increment_disp.val(17,7)=IncrDis17(6);

    increment_dispDu.val(18,1)=IncrDis18(0);
    increment_dispDu.val(18,2)=IncrDis18(1);
    increment_dispDu.val(18,3)=IncrDis18(2);
//    increment_disp.val(18,4)=IncrDis18(3); increment_disp.val(18,5)=IncrDis18(4);increment_disp.val(18,6)=IncrDis18(5);
//    increment_disp.val(18,7)=IncrDis18(6);

    increment_dispDu.val(19,1)=IncrDis19(0);
    increment_dispDu.val(19,2)=IncrDis19(1);
    increment_dispDu.val(19,3)=IncrDis19(2);
//    increment_disp.val(19,4)=IncrDis19(3); increment_disp.val(19,5)=IncrDis19(4);increment_disp.val(19,6)=IncrDis19(5);
//    increment_disp.val(19,7)=IncrDis19(6);

    increment_dispDu.val(20,1)=IncrDis20(0);
    increment_dispDu.val(20,2)=IncrDis20(1);
    increment_dispDu.val(20,3)=IncrDis20(2);
//    increment_disp.val(20,4)=IncrDis20(3); increment_disp.val(20,5)=IncrDis20(4);increment_disp.val(20,6)=IncrDis20(5);
//    increment_disp.val(20,7)=IncrDis20(6);


    return increment_dispDu;
  }
////#############################################################################
//=========================================================================
// Incremental displacement of fluid part-- U     02/10/2002
//=========================================================================

//=========================================================================
// Incremental displacement of pressure part-- p         05/20/2002
//=========================================================================

tensor TwentyNodeBrick_u_p_U::incr_dispDp()
  {
    const int dimensions[] = {20};
    tensor increment_dispDp(1, dimensions, 0.0);

    //for ( int i=0 ; i<8 ; i++ )   // Xiaoyan changed from 20 to 8 for 8 nodes
    //
    //  {
    //    // increment_disp.val(i+1,1) = nodes[ G_N_numbs[i] ].incremental_translation_x();
    //    // increment_disp.val(i+1,2) = nodes[ G_N_numbs[i] ].incremental_translation_y();
    //    // increment_disp.val(i+1,3) = nodes[ G_N_numbs[i] ].incremental_translation_z();
    //    // Xiaoyan changed to the following 09/27/00
    //    Vector IncremenDis = nodes[ G_N_numbs[i] ].getIncrDisp();
    //
    //    increment_disp.val(i+1,1) = IncremenDis(0);
    //    increment_disp.val(i+1,2) = IncremenDis(1);
    //    increment_disp.val(i+1,3) = IncremenDis(2);
    //
    //  }

    //Zhaohui using node pointers, which come from the Domain
    //const Vector &TotDis1 = nd1Ptr->getTrialDisp();
    //const Vector &incrdelDis1 = nd1Ptr->getIncrDisp();
    //Have to get IncrDeltaDisp, not IncrDisp for cumulation of incr_disp
    const Vector &IncrDis1 = nd1Ptr->getIncrDeltaDisp();
    const Vector &IncrDis2 = nd2Ptr->getIncrDeltaDisp();
    const Vector &IncrDis3 = nd3Ptr->getIncrDeltaDisp();
    const Vector &IncrDis4 = nd4Ptr->getIncrDeltaDisp();

    const Vector &IncrDis5 = nd5Ptr->getIncrDeltaDisp();
    const Vector &IncrDis6 = nd6Ptr->getIncrDeltaDisp();
    const Vector &IncrDis7 = nd7Ptr->getIncrDeltaDisp();
    const Vector &IncrDis8 = nd8Ptr->getIncrDeltaDisp();

    const Vector &IncrDis9  = nd9Ptr->getIncrDeltaDisp();
    const Vector &IncrDis10 = nd10Ptr->getIncrDeltaDisp();
    const Vector &IncrDis11 = nd11Ptr->getIncrDeltaDisp();
    const Vector &IncrDis12 = nd12Ptr->getIncrDeltaDisp();

    const Vector &IncrDis13 = nd13Ptr->getIncrDeltaDisp();
    const Vector &IncrDis14 = nd14Ptr->getIncrDeltaDisp();
    const Vector &IncrDis15 = nd15Ptr->getIncrDeltaDisp();
    const Vector &IncrDis16 = nd16Ptr->getIncrDeltaDisp();

    const Vector &IncrDis17 = nd17Ptr->getIncrDeltaDisp();
    const Vector &IncrDis18 = nd18Ptr->getIncrDeltaDisp();
    const Vector &IncrDis19 = nd19Ptr->getIncrDeltaDisp();
    const Vector &IncrDis20 = nd20Ptr->getIncrDeltaDisp();

// Get the fourth incremental displacement for pressure part. Xiaoyan 05/20/2002

    increment_dispDp.val(1)=IncrDis1(3);
    increment_dispDp.val(2)=IncrDis2(3);
    increment_dispDp.val(3)=IncrDis3(3);
    increment_dispDp.val(4)=IncrDis4(3);
    increment_dispDp.val(5)=IncrDis5(3);
    increment_dispDp.val(6)=IncrDis6(3);
    increment_dispDp.val(7)=IncrDis7(3);
    increment_dispDp.val(8)=IncrDis8(3);
    increment_dispDp.val(9)=IncrDis9(3);
    increment_dispDp.val(10)=IncrDis10(3);
    increment_dispDp.val(11)=IncrDis11(3);
    increment_dispDp.val(12)=IncrDis12(3);
    increment_dispDp.val(13)=IncrDis13(3);
    increment_dispDp.val(14)=IncrDis14(3);
    increment_dispDp.val(15)=IncrDis15(3);
    increment_dispDp.val(16)=IncrDis16(3);
    increment_dispDp.val(17)=IncrDis17(3);
    increment_dispDp.val(18)=IncrDis18(3);
    increment_dispDp.val(19)=IncrDis19(3);
    increment_dispDp.val(20)=IncrDis20(3);

    return increment_dispDp;
  }

////#############################################################################

//=========================================================================
// Incremental displacement of fluid part-- U     02/10/2002
//=========================================================================
tensor TwentyNodeBrick_u_p_U::incr_dispDU()
  {
    const int dimensions[] = {20,3};
    tensor increment_dispDU(2, dimensions, 0.0);

    //for ( int i=0 ; i<8 ; i++ )   // Xiaoyan changed from 20 to 8 for 8 nodes
    //
    //  {
    //    // increment_disp.val(i+1,1) = nodes[ G_N_numbs[i] ].incremental_translation_x();
    //    // increment_disp.val(i+1,2) = nodes[ G_N_numbs[i] ].incremental_translation_y();
    //    // increment_disp.val(i+1,3) = nodes[ G_N_numbs[i] ].incremental_translation_z();
    //    // Xiaoyan changed to the following 09/27/00
    //    Vector IncremenDis = nodes[ G_N_numbs[i] ].getIncrDisp();
    //
    //    increment_disp.val(i+1,1) = IncremenDis(0);
    //    increment_disp.val(i+1,2) = IncremenDis(1);
    //    increment_disp.val(i+1,3) = IncremenDis(2);
    //
    //  }

    //Zhaohui using node pointers, which come from the Domain
    //const Vector &TotDis1 = nd1Ptr->getTrialDisp();
    //const Vector &incrdelDis1 = nd1Ptr->getIncrDisp();
    //Have to get IncrDeltaDisp, not IncrDisp for cumulation of incr_disp
    const Vector &IncrDis1 = nd1Ptr->getIncrDeltaDisp();
    const Vector &IncrDis2 = nd2Ptr->getIncrDeltaDisp();
    const Vector &IncrDis3 = nd3Ptr->getIncrDeltaDisp();
    const Vector &IncrDis4 = nd4Ptr->getIncrDeltaDisp();

    const Vector &IncrDis5 = nd5Ptr->getIncrDeltaDisp();
    const Vector &IncrDis6 = nd6Ptr->getIncrDeltaDisp();
    const Vector &IncrDis7 = nd7Ptr->getIncrDeltaDisp();
    const Vector &IncrDis8 = nd8Ptr->getIncrDeltaDisp();

    const Vector &IncrDis9  = nd9Ptr->getIncrDeltaDisp();
    const Vector &IncrDis10 = nd10Ptr->getIncrDeltaDisp();
    const Vector &IncrDis11 = nd11Ptr->getIncrDeltaDisp();
    const Vector &IncrDis12 = nd12Ptr->getIncrDeltaDisp();

    const Vector &IncrDis13 = nd13Ptr->getIncrDeltaDisp();
    const Vector &IncrDis14 = nd14Ptr->getIncrDeltaDisp();
    const Vector &IncrDis15 = nd15Ptr->getIncrDeltaDisp();
    const Vector &IncrDis16 = nd16Ptr->getIncrDeltaDisp();

    const Vector &IncrDis17 = nd17Ptr->getIncrDeltaDisp();
    const Vector &IncrDis18 = nd18Ptr->getIncrDeltaDisp();
    const Vector &IncrDis19 = nd19Ptr->getIncrDeltaDisp();
    const Vector &IncrDis20 = nd20Ptr->getIncrDeltaDisp();

// Get the Last three incremental displacement for fluid part. Xiaoyan 02/08/2002

//    increment_disp.val(1,1)=IncrDis1(0); increment_disp.val(1,2)=IncrDis1(1);increment_disp.val(1,3)=IncrDis1(2);
//    increment_disp.val(1,4)=IncrDis1(3);
    increment_dispDU.val(1,1)=IncrDis1(4);
    increment_dispDU.val(1,2)=IncrDis1(5);
    increment_dispDU.val(1,3)=IncrDis1(6);

//    increment_disp.val(2,1)=IncrDis2(0); increment_disp.val(2,2)=IncrDis2(1);increment_disp.val(2,3)=IncrDis2(2);
//    increment_disp.val(2,4)=IncrDis2(3);
    increment_dispDU.val(2,1)=IncrDis2(4);
    increment_dispDU.val(2,2)=IncrDis2(5);
    increment_dispDU.val(2,3)=IncrDis2(6);

//    increment_disp.val(3,1)=IncrDis3(0); increment_disp.val(3,2)=IncrDis3(1);increment_disp.val(3,3)=IncrDis3(2);
//    increment_disp.val(3,4)=IncrDis3(3);
    increment_dispDU.val(3,1)=IncrDis3(4);
    increment_dispDU.val(3,2)=IncrDis3(5);
    increment_dispDU.val(3,3)=IncrDis3(6);

//    increment_disp.val(4,1)=IncrDis4(0); increment_disp.val(4,2)=IncrDis4(1);increment_disp.val(4,3)=IncrDis4(2);
//    increment_disp.val(4,4)=IncrDis4(3);
    increment_dispDU.val(4,1)=IncrDis4(4);
    increment_dispDU.val(4,2)=IncrDis4(5);
    increment_dispDU.val(4,3)=IncrDis4(6);

//    increment_disp.val(5,1)=IncrDis5(0); increment_disp.val(5,2)=IncrDis5(1);increment_disp.val(5,3)=IncrDis5(2);
//    increment_disp.val(5,4)=IncrDis5(3);
    increment_dispDU.val(5,1)=IncrDis5(4);
    increment_dispDU.val(5,2)=IncrDis5(5);
    increment_dispDU.val(5,3)=IncrDis5(6);

//    increment_disp.val(6,1)=IncrDis6(0); increment_disp.val(6,2)=IncrDis6(1);increment_disp.val(6,3)=IncrDis6(2);
//    increment_disp.val(6,4)=IncrDis6(3);
    increment_dispDU.val(6,1)=IncrDis6(4);
    increment_dispDU.val(4,2)=IncrDis6(5);
    increment_dispDU.val(6,3)=IncrDis6(6);

//    increment_disp.val(7,1)=IncrDis7(0); increment_disp.val(7,2)=IncrDis7(1);increment_disp.val(7,3)=IncrDis7(2);
//    increment_disp.val(7,4)=IncrDis7(3);
    increment_dispDU.val(7,1)=IncrDis7(4);
    increment_dispDU.val(7,2)=IncrDis7(5);
    increment_dispDU.val(7,3)=IncrDis7(6);

//    increment_disp.val(8,1)=IncrDis8(0); increment_disp.val(8,2)=IncrDis8(1);increment_disp.val(8,3)=IncrDis8(2);
//    increment_disp.val(8,4)=IncrDis8(3);
    increment_dispDU.val(8,1)=IncrDis8(4);
    increment_dispDU.val(8,2)=IncrDis8(5);
    increment_dispDU.val(8,3)=IncrDis8(6);

//    increment_disp.val(9,1)=IncrDis9(0); increment_disp.val(9,2)=IncrDis9(1);increment_disp.val(9,3)=IncrDis9(2);
//    increment_disp.val(9,4)=IncrDis9(3);
    increment_dispDU.val(9,1)=IncrDis9(4);
    increment_dispDU.val(9,2)=IncrDis9(5);
    increment_dispDU.val(9,3)=IncrDis9(6);

//    increment_disp.val(10,1)=IncrDis10(0); increment_disp.val(10,2)=IncrDis10(1);increment_disp.val(10,3)=IncrDis10(2);
//    increment_disp.val(10,4)=IncrDis10(3);
    increment_dispDU.val(10,1)=IncrDis10(4);
    increment_dispDU.val(10,2)=IncrDis10(5);
    increment_dispDU.val(10,3)=IncrDis10(6);

//    increment_disp.val(11,1)=IncrDis11(0); increment_disp.val(11,2)=IncrDis11(1);increment_disp.val(11,3)=IncrDis11(2);
//    increment_disp.val(11,4)=IncrDis11(3);
    increment_dispDU.val(11,1)=IncrDis11(4);
    increment_dispDU.val(11,2)=IncrDis11(5);
    increment_dispDU.val(11,3)=IncrDis11(6);

//    increment_disp.val(12,1)=IncrDis12(0); increment_disp.val(12,2)=IncrDis12(1);increment_disp.val(12,3)=IncrDis12(2);
//    increment_disp.val(12,4)=IncrDis12(3);
    increment_dispDU.val(12,1)=IncrDis12(4);
    increment_dispDU.val(12,2)=IncrDis12(5);
    increment_dispDU.val(12,3)=IncrDis12(6);

//    increment_disp.val(13,1)=IncrDis13(0); increment_disp.val(13,2)=IncrDis13(1);increment_disp.val(13,3)=IncrDis13(2);
//    increment_disp.val(13,4)=IncrDis13(3);
    increment_dispDU.val(13,1)=IncrDis13(4);
    increment_dispDU.val(13,2)=IncrDis13(5);
    increment_dispDU.val(13,3)=IncrDis13(6);

//    increment_disp.val(14,1)=IncrDis14(0); increment_disp.val(14,2)=IncrDis14(1);increment_disp.val(14,3)=IncrDis14(2);
//    increment_disp.val(14,4)=IncrDis14(3);
    increment_dispDU.val(14,1)=IncrDis14(4);
    increment_dispDU.val(14,2)=IncrDis14(5);
    increment_dispDU.val(14,3)=IncrDis14(6);

//    increment_disp.val(15,1)=IncrDis15(0); increment_disp.val(15,2)=IncrDis15(1);increment_disp.val(15,3)=IncrDis15(2);
//    increment_disp.val(15,4)=IncrDis15(3);
    increment_dispDU.val(15,1)=IncrDis15(4);
    increment_dispDU.val(15,2)=IncrDis15(5);
    increment_dispDU.val(15,3)=IncrDis15(6);

//    increment_disp.val(16,1)=IncrDis16(0); increment_disp.val(16,2)=IncrDis16(1);increment_disp.val(16,3)=IncrDis16(2);
//    increment_disp.val(16,4)=IncrDis16(3);
    increment_dispDU.val(16,1)=IncrDis16(4);
    increment_dispDU.val(16,2)=IncrDis16(5);
    increment_dispDU.val(16,3)=IncrDis16(6);

//    increment_disp.val(17,1)=IncrDis17(0); increment_disp.val(17,2)=IncrDis17(1);increment_disp.val(17,3)=IncrDis17(2);
//    increment_disp.val(17,4)=IncrDis17(3);
    increment_dispDU.val(17,1)=IncrDis17(4);
    increment_dispDU.val(17,2)=IncrDis17(5);
    increment_dispDU.val(17,3)=IncrDis17(6);

//    increment_disp.val(18,1)=IncrDis18(0); increment_disp.val(18,2)=IncrDis18(1);increment_disp.val(18,3)=IncrDis18(2);
//    increment_disp.val(18,4)=IncrDis18(3);
    increment_dispDU.val(18,1)=IncrDis18(4);
    increment_dispDU.val(18,2)=IncrDis18(5);
    increment_dispDU.val(18,3)=IncrDis18(6);

//    increment_disp.val(19,1)=IncrDis19(0); increment_disp.val(19,2)=IncrDis19(1);increment_disp.val(19,3)=IncrDis19(2);
//    increment_disp.val(19,4)=IncrDis19(3);
    increment_dispDU.val(19,1)=IncrDis19(4);
    increment_dispDU.val(19,2)=IncrDis19(5);
    increment_dispDU.val(19,3)=IncrDis19(6);

//    increment_disp.val(20,1)=IncrDis20(0); increment_disp.val(20,2)=IncrDis20(1);increment_disp.val(20,3)=IncrDis20(2);
//    increment_disp.val(20,4)=IncrDis20(3);
    increment_dispDU.val(20,1)=IncrDis20(4);
    increment_dispDU.val(20,2)=IncrDis20(5);
    increment_dispDU.val(20,3)=IncrDis20(6);


    return increment_dispDU;
  }


////#############################################################################
//=========================================================================
// Total displacement of solid part --u    02/10/2002
//=========================================================================
tensor TwentyNodeBrick_u_p_U::total_dispDu()
  {
    const int dimensions[] = {20,3};
    tensor total_dispDu(2, dimensions, 0.0);

    // Get first three total displacement for solid part. Xiaoyan 02/08/2002
    const Vector &TotDis1 = nd1Ptr->getTrialDisp();
    opserr<<"\ntot node " << nd1Ptr->getTag() <<" ux "<< TotDis1(0) <<" uy "<< TotDis1(1) << " uz "<< TotDis1(2)
                        << endln;

    const Vector &TotDis2 = nd2Ptr->getTrialDisp();
    opserr << "tot node " << nd2Ptr->getTag() << " ux " << TotDis2(0) <<" uy "<< TotDis2(1) << " uz "<< TotDis2(2)
                        << endln;

    const Vector &TotDis3 = nd3Ptr->getTrialDisp();
    opserr << "tot node " << nd3Ptr->getTag() << " ux " << TotDis3(0) <<" uy "<< TotDis3(1) << " uz "<< TotDis3(2)
                        << endln;

    const Vector &TotDis4 = nd4Ptr->getTrialDisp();
    opserr << "tot node " << nd4Ptr->getTag() << " ux " << TotDis4(0) <<"uy "<< TotDis4(1) << " uz "<< TotDis4(2)
                        << endln;

    const Vector &TotDis5 = nd5Ptr->getTrialDisp();
    opserr << "tot node " << nd5Ptr->getTag() << " ux " << TotDis5(0) <<" uy "<< TotDis5(1) << " uz "<< TotDis5(2)
                        << endln;

    const Vector &TotDis6 = nd6Ptr->getTrialDisp();
    opserr << "tot node " << nd6Ptr->getTag() << " ux " << TotDis6(0) <<" uy "<< TotDis6(1) << " uz "<< TotDis6(2)
                        << endln;

    const Vector &TotDis7 = nd7Ptr->getTrialDisp();
    opserr << "tot node " << nd7Ptr->getTag() << " ux " << TotDis7(0) <<" uy "<< TotDis7(1) << " uz "<< TotDis7(2)
                        << endln;

    const Vector &TotDis8 = nd8Ptr->getTrialDisp();
    opserr << "tot node " << nd8Ptr->getTag() << " ux " << TotDis8(0) <<" uy "<< TotDis8(1) << " uz "<< TotDis8(2)
                        << endln;

    const Vector &TotDis9 = nd9Ptr->getTrialDisp();
    opserr << "tot node " << nd9Ptr->getTag() << " ux " << TotDis9(0) <<" uy "<< TotDis9(1) << " uz "<< TotDis9(2)
                        << endln;

    const Vector &TotDis10 = nd10Ptr->getTrialDisp();
    opserr << "tot node " << nd10Ptr->getTag() << " ux " << TotDis10(0) <<" uy "<< TotDis10(1) << " uz "<< TotDis10(2)
                        << endln;

    const Vector &TotDis11 = nd11Ptr->getTrialDisp();
    opserr << "tot node " << nd11Ptr->getTag() << " ux " << TotDis11(0) <<" uy "<< TotDis11(1) << " uz "<< TotDis11(2)
                        << endln;

    const Vector &TotDis12 = nd12Ptr->getTrialDisp();
    opserr << "tot node " << nd12Ptr->getTag() << " ux " << TotDis12(0) <<" uy "<< TotDis12(1) << " uz "<< TotDis12(2)
                        << endln;


    const Vector &TotDis13 = nd13Ptr->getTrialDisp();
    opserr << "tot node " << nd13Ptr->getTag() << " ux " << TotDis13(0) <<" uy "<< TotDis13(1) << " uz "<< TotDis13(2)
                        << endln;

    const Vector &TotDis14 = nd14Ptr->getTrialDisp();
    opserr << "tot node " << nd14Ptr->getTag() << " ux " << TotDis14(0) <<" uy "<< TotDis14(1) << " uz "<< TotDis14(2)
                        << endln;

    const Vector &TotDis15 = nd15Ptr->getTrialDisp();
    opserr << "tot node " << nd15Ptr->getTag() << " ux " << TotDis15(0) <<" uy "<< TotDis15(1) << " uz "<< TotDis15(2)
                        << endln;

    const Vector &TotDis16 = nd16Ptr->getTrialDisp();
    opserr << "tot node " << nd16Ptr->getTag() << " ux " << TotDis16(0) <<" uy "<< TotDis16(1) << " uz "<< TotDis16(2)
                        << endln;

    const Vector &TotDis17 = nd17Ptr->getTrialDisp();
    opserr << "tot node " << nd17Ptr->getTag() << " ux " << TotDis17(0) <<" uy "<< TotDis17(1) << " uz "<< TotDis17(2)
                        << endln;

    const Vector &TotDis18 = nd18Ptr->getTrialDisp();
    opserr << "tot node " << nd18Ptr->getTag() << " ux " << TotDis18(0) <<" uy "<< TotDis18(1) << " uz "<< TotDis18(2)
                        << endln;

    const Vector &TotDis19 = nd19Ptr->getTrialDisp();
    opserr << "tot node " << nd19Ptr->getTag() << " ux " << TotDis19(0) <<" uy "<< TotDis19(1) << " uz "<< TotDis19(2)
                        << endln;

    const Vector &TotDis20 = nd20Ptr->getTrialDisp();
    opserr << "tot node " << nd20Ptr->getTag() << " ux " << TotDis20(0) <<" uy "<< TotDis20(1) << " uz "<< TotDis20(2)
                        << endln;



    total_dispDu.val(1,1)=TotDis1(0); total_dispDu.val(1,2)=TotDis1(1);total_dispDu.val(1,3)=TotDis1(2);
//    total_disp.val(1,4)=TotDis1(3); total_disp.val(1,5)=TotDis1(4);total_disp.val(1,6)=TotDis1(5);
//    total_disp.val(1,7)=TotDis1(6);

    total_dispDu.val(2,1)=TotDis2(0); total_dispDu.val(2,2)=TotDis2(1);total_dispDu.val(2,3)=TotDis2(2);
//    total_disp.val(2,4)=TotDis2(3); total_disp.val(2,5)=TotDis2(4);total_disp.val(2,6)=TotDis2(5);
//    total_disp.val(2,7)=TotDis2(6);

    total_dispDu.val(3,1)=TotDis3(0); total_dispDu.val(3,2)=TotDis3(1);total_dispDu.val(3,3)=TotDis3(2);
//    total_disp.val(3,4)=TotDis3(3); total_disp.val(3,5)=TotDis3(4);total_disp.val(3,6)=TotDis3(5);
//    total_disp.val(3,7)=TotDis3(6);

    total_dispDu.val(4,1)=TotDis4(0); total_dispDu.val(4,2)=TotDis4(1);total_dispDu.val(4,3)=TotDis4(2);
//    total_disp.val(4,4)=TotDis4(3); total_disp.val(4,5)=TotDis4(4);total_disp.val(4,6)=TotDis4(5);
//    total_disp.val(4,7)=TotDis4(6);

    total_dispDu.val(5,1)=TotDis5(0); total_dispDu.val(5,2)=TotDis5(1);total_dispDu.val(5,3)=TotDis5(2);
//    total_disp.val(5,4)=TotDis5(3); total_disp.val(5,5)=TotDis5(4);total_disp.val(5,6)=TotDis5(5);
//    total_disp.val(5,7)=TotDis5(6);

    total_dispDu.val(6,1)=TotDis6(0); total_dispDu.val(6,2)=TotDis6(1);total_dispDu.val(6,3)=TotDis6(2);
//    total_disp.val(6,4)=TotDis6(3); total_disp.val(6,5)=TotDis6(4);total_disp.val(6,6)=TotDis6(5);
//    total_disp.val(6,7)=TotDis6(6);

    total_dispDu.val(7,1)=TotDis7(0); total_dispDu.val(7,2)=TotDis7(1);total_dispDu.val(7,3)=TotDis7(2);
//    total_disp.val(7,4)=TotDis7(3); total_disp.val(7,5)=TotDis7(4);total_disp.val(7,6)=TotDis7(5);
//    total_disp.val(7,7)=TotDis7(6);

    total_dispDu.val(8,1)=TotDis8(0); total_dispDu.val(8,2)=TotDis8(1);total_dispDu.val(8,3)=TotDis8(2);
//    total_disp.val(8,4)=TotDis8(3); total_disp.val(8,5)=TotDis8(4);total_disp.val(8,6)=TotDis8(5);
//    total_disp.val(8,7)=TotDis8(6);


    total_dispDu.val(9,1)=TotDis9(0);
    total_dispDu.val(9,2)=TotDis9(1);
    total_dispDu.val(9,3)=TotDis9(2);
//    total_disp.val(9,4)=TotDis9(3); total_disp.val(9,5)=TotDis9(4);total_disp.val(9,6)=TotDis9(5);
//    total_disp.val(9,7)=TotDis9(6);

    total_dispDu.val(10,1)=TotDis10(0);
    total_dispDu.val(10,2)=TotDis10(1);
    total_dispDu.val(10,3)=TotDis10(2);
//    total_disp.val(10,4)=TotDis10(3); total_disp.val(10,5)=TotDis10(4);total_disp.val(10,6)=TotDis10(5);
//    total_disp.val(10,7)=TotDis10(6);

    total_dispDu.val(11,1)=TotDis11(0);
    total_dispDu.val(11,2)=TotDis11(1);
    total_dispDu.val(11,3)=TotDis11(2);
//    total_disp.val(11,4)=TotDis11(3); total_disp.val(11,5)=TotDis11(4);total_disp.val(11,6)=TotDis11(5);
//    total_disp.val(11,7)=TotDis11(6);

    total_dispDu.val(12,1)=TotDis12(0);
    total_dispDu.val(12,2)=TotDis12(1);
    total_dispDu.val(12,3)=TotDis12(2);
//    total_disp.val(12,4)=TotDis12(3); total_disp.val(12,5)=TotDis12(4);total_disp.val(12,6)=TotDis12(5);
//    total_disp.val(12,7)=TotDis12(6);


    total_dispDu.val(13,1)=TotDis13(0);
    total_dispDu.val(13,2)=TotDis13(1);
    total_dispDu.val(13,3)=TotDis13(2);
//    total_disp.val(13,4)=TotDis13(3); total_disp.val(13,5)=TotDis13(4);total_disp.val(13,6)=TotDis13(5);
//    total_disp.val(13,7)=TotDis13(6);

    total_dispDu.val(14,1)=TotDis14(0);
    total_dispDu.val(14,2)=TotDis14(1);
    total_dispDu.val(14,3)=TotDis14(2);
//    total_disp.val(14,4)=TotDis14(3); total_disp.val(14,5)=TotDis14(4);total_disp.val(14,6)=TotDis14(5);
//    total_disp.val(14,7)=TotDis14(6);

    total_dispDu.val(15,1)=TotDis15(0);
    total_dispDu.val(15,2)=TotDis15(1);
    total_dispDu.val(15,3)=TotDis15(2);
//    total_disp.val(15,4)=TotDis15(3); total_disp.val(15,5)=TotDis15(4);total_disp.val(15,6)=TotDis15(5);
//    total_disp.val(15,7)=TotDis15(6);

    total_dispDu.val(16,1)=TotDis16(0);
    total_dispDu.val(16,2)=TotDis16(1);
    total_dispDu.val(16,3)=TotDis16(2);
//    total_disp.val(16,4)=TotDis16(3); total_disp.val(16,5)=TotDis16(4);total_disp.val(16,6)=TotDis16(5);
//    total_disp.val(16,7)=TotDis16(6);

    total_dispDu.val(17,1)=TotDis17(0);
    total_dispDu.val(17,2)=TotDis17(1);
    total_dispDu.val(17,3)=TotDis17(2);
//    total_disp.val(17,4)=TotDis17(3); total_disp.val(17,5)=TotDis17(4);total_disp.val(17,6)=TotDis17(5);
//    total_disp.val(17,7)=TotDis17(6);

    total_dispDu.val(18,1)=TotDis18(0);
    total_dispDu.val(18,2)=TotDis18(1);
    total_dispDu.val(18,3)=TotDis18(2);
//    total_disp.val(18,4)=TotDis18(3); total_disp.val(18,5)=TotDis18(4);total_disp.val(18,6)=TotDis18(5);
//    total_disp.val(18,7)=TotDis18(6);

    total_dispDu.val(19,1)=TotDis19(0);
    total_dispDu.val(19,2)=TotDis19(1);
    total_dispDu.val(19,3)=TotDis19(2);
//    total_disp.val(19,4)=TotDis19(3); total_disp.val(19,5)=TotDis19(4);total_disp.val(19,6)=TotDis19(5);
//    total_disp.val(19,7)=TotDis19(6);

    total_dispDu.val(20,1)=TotDis20(0);
    total_dispDu.val(20,2)=TotDis20(1);
    total_dispDu.val(20,3)=TotDis20(2);
//    total_disp.val(20,4)=TotDis20(3); total_disp.val(20,5)=TotDis20(4);total_disp.val(20,6)=TotDis20(5);
//    total_disp.val(20,7)=TotDis20(6);


    return total_dispDu;
  }

////#############################################################################
//=========================================================================
// Total displacement of fluid part--U     02/10/2002
//=========================================================================
tensor TwentyNodeBrick_u_p_U::total_dispDU()
  {
    const int dimensions[] = {20,3};
    tensor total_dispDU(2, dimensions, 0.0);

    //Zhaohui using node pointers, which come from the Domain
    const Vector &TotDis1 = nd1Ptr->getTrialDisp();
    opserr<<"\ntot node " << nd1Ptr->getTag()
      <<" Ux "<< TotDis1(4) << " Uy "<< TotDis1(5)<< " Uz "<< TotDis1(6)<< endln;

    const Vector &TotDis2 = nd2Ptr->getTrialDisp();
    opserr << "tot node " << nd2Ptr->getTag()
          <<" Ux "<< TotDis2(4) << " Uy "<< TotDis2(5)<< " Uz "<< TotDis2(6)<< endln;

    const Vector &TotDis3 = nd3Ptr->getTrialDisp();
    opserr << "tot node " << nd3Ptr->getTag()
      <<" Ux "<< TotDis3(4) << " Uy "<< TotDis3(5)<< " Uz "<< TotDis3(6)<< endln;

    const Vector &TotDis4 = nd4Ptr->getTrialDisp();
    opserr << "tot node " << nd4Ptr->getTag()
      <<" Ux "<< TotDis4(4) << " Uy "<< TotDis4(5)<< " Uz "<< TotDis4(6)<< endln;

    const Vector &TotDis5 = nd5Ptr->getTrialDisp();
    opserr << "tot node " << nd5Ptr->getTag()
      <<" Ux "<< TotDis5(4) << " Uy "<< TotDis5(5)<< " Uz "<< TotDis5(6)<< endln;

    const Vector &TotDis6 = nd6Ptr->getTrialDisp();
    opserr << "tot node " << nd6Ptr->getTag()
      <<" Ux "<< TotDis6(4) << " Uy "<< TotDis6(5)<< " Uz "<< TotDis6(6)<< endln;

    const Vector &TotDis7 = nd7Ptr->getTrialDisp();
    opserr << "tot node " << nd7Ptr->getTag()
      <<" Ux "<< TotDis7(4) << " Uy "<< TotDis7(5)<< " Uz "<< TotDis7(6)<< endln;

    const Vector &TotDis8 = nd8Ptr->getTrialDisp();
    opserr << "tot node " << nd8Ptr->getTag()
      <<" Ux "<< TotDis8(4) << " Uy "<< TotDis8(5)<< " Uz "<< TotDis8(6)<< endln;

    const Vector &TotDis9 = nd9Ptr->getTrialDisp();
    opserr << "tot node " << nd9Ptr->getTag()
      <<" Ux "<< TotDis9(4) << " Uy "<< TotDis9(5)<< " Uz "<< TotDis9(6)<< endln;

    const Vector &TotDis10 = nd10Ptr->getTrialDisp();
    opserr << "tot node " << nd10Ptr->getTag()
      <<" Ux "<< TotDis10(4) << " Uy "<< TotDis10(5)<< " Uz "<< TotDis10(6)<< endln;

    const Vector &TotDis11 = nd11Ptr->getTrialDisp();
    opserr << "tot node " << nd11Ptr->getTag()
      <<" Ux "<< TotDis11(4) << " Uy "<< TotDis11(5)<< " Uz "<< TotDis11(6)<< endln;

    const Vector &TotDis12 = nd12Ptr->getTrialDisp();
    opserr << "tot node " << nd12Ptr->getTag()
      <<" Ux "<< TotDis12(4) << " Uy "<< TotDis12(5)<< " Uz "<< TotDis12(6)<< endln;


    const Vector &TotDis13 = nd13Ptr->getTrialDisp();
    opserr << "tot node " << nd13Ptr->getTag()
      <<" Ux "<< TotDis13(4) << " Uy "<< TotDis13(5)<< " Uz "<< TotDis13(6)<< endln;

    const Vector &TotDis14 = nd14Ptr->getTrialDisp();
    opserr << "tot node " << nd14Ptr->getTag()
      <<" Ux "<< TotDis14(4) << " Uy "<< TotDis14(5)<< " Uz "<< TotDis14(6)<< endln;

    const Vector &TotDis15 = nd15Ptr->getTrialDisp();
    opserr << "tot node " << nd15Ptr->getTag()
      <<" Ux "<< TotDis15(4) << " Uy "<< TotDis15(5)<< " Uz "<< TotDis15(6)<< endln;

    const Vector &TotDis16 = nd16Ptr->getTrialDisp();
    opserr << "tot node " << nd16Ptr->getTag()
      <<" Ux "<< TotDis16(4) << " Uy "<< TotDis16(5)<< " Uz "<< TotDis16(6)<< endln;

    const Vector &TotDis17 = nd17Ptr->getTrialDisp();
    opserr << "tot node " << nd17Ptr->getTag()
      <<" Ux "<< TotDis17(4) << " Uy "<< TotDis17(5)<< " Uz "<< TotDis17(6)<< endln;

    const Vector &TotDis18 = nd18Ptr->getTrialDisp();
    opserr << "tot node " << nd18Ptr->getTag()
      <<" Ux "<< TotDis18(4) << " Uy "<< TotDis18(5)<< " Uz "<< TotDis18(6)<< endln;

    const Vector &TotDis19 = nd19Ptr->getTrialDisp();
    opserr << "tot node " << nd19Ptr->getTag()
      <<" Ux "<< TotDis19(4) << " Uy "<< TotDis19(5)<< " Uz "<< TotDis19(6)<< endln;

    const Vector &TotDis20 = nd20Ptr->getTrialDisp();
    opserr << "tot node " << nd20Ptr->getTag()
      <<" Ux "<< TotDis20(4) << " Uy "<< TotDis20(5)<< " Uz "<< TotDis20(6)<< endln;



    total_dispDU.val(1,1)=TotDis1(4);
    total_dispDU.val(1,2)=TotDis1(5);
    total_dispDU.val(1,3)=TotDis1(6);

    total_dispDU.val(2,1)=TotDis2(4);
    total_dispDU.val(2,2)=TotDis2(5);
    total_dispDU.val(2,3)=TotDis2(6);

    total_dispDU.val(3,1)=TotDis3(4);
    total_dispDU.val(3,2)=TotDis3(5);
    total_dispDU.val(3,3)=TotDis3(6);

    total_dispDU.val(4,1)=TotDis4(4);
    total_dispDU.val(4,2)=TotDis4(5);
    total_dispDU.val(4,3)=TotDis4(6);

    total_dispDU.val(5,1)=TotDis5(4);
    total_dispDU.val(5,2)=TotDis5(5);
    total_dispDU.val(5,3)=TotDis5(6);

    total_dispDU.val(6,1)=TotDis6(4);
    total_dispDU.val(6,2)=TotDis6(5);
    total_dispDU.val(6,3)=TotDis6(6);

    total_dispDU.val(7,1)=TotDis7(4);
    total_dispDU.val(7,2)=TotDis7(5);
    total_dispDU.val(7,3)=TotDis7(6);

    total_dispDU.val(8,1)=TotDis8(4);total_dispDU.val(8,2)=TotDis8(5);
    total_dispDU.val(8,3)=TotDis8(6);


    total_dispDU.val(9,1)=TotDis9(4);total_dispDU.val(9,2)=TotDis9(5);
    total_dispDU.val(9,3)=TotDis9(6);

    total_dispDU.val(10,1)=TotDis10(4);total_dispDU.val(10,2)=TotDis10(5);
    total_dispDU.val(10,3)=TotDis10(6);

    total_dispDU.val(11,1)=TotDis11(4);total_dispDU.val(11,2)=TotDis11(5);
    total_dispDU.val(11,3)=TotDis11(6);

    total_dispDU.val(12,1)=TotDis12(4);total_dispDU.val(12,2)=TotDis12(5);
    total_dispDU.val(12,3)=TotDis12(6);


    total_dispDU.val(13,1)=TotDis13(4);total_dispDU.val(13,2)=TotDis13(5);
    total_dispDU.val(13,3)=TotDis13(6);

    total_dispDU.val(14,1)=TotDis14(4);total_dispDU.val(14,2)=TotDis14(5);
    total_dispDU.val(14,3)=TotDis14(6);

    total_dispDU.val(15,1)=TotDis15(4);total_dispDU.val(15,2)=TotDis15(5);
    total_dispDU.val(15,3)=TotDis15(6);

    total_dispDU.val(16,1)=TotDis16(4);total_dispDU.val(16,2)=TotDis16(5);
    total_dispDU.val(16,3)=TotDis16(6);

    total_dispDU.val(17,1)=TotDis17(4);total_dispDU.val(17,2)=TotDis17(5);
    total_dispDU.val(17,3)=TotDis17(6);

    total_dispDU.val(18,1)=TotDis18(4);total_dispDU.val(18,2)=TotDis18(5);
    total_dispDU.val(18,3)=TotDis18(6);

    total_dispDU.val(19,1)=TotDis19(4);total_dispDU.val(19,2)=TotDis19(5);
    total_dispDU.val(19,3)=TotDis19(6);

    total_dispDU.val(20,1)=TotDis20(4);total_dispDU.val(20,2)=TotDis20(5);
    total_dispDU.val(20,3)=TotDis20(6);


    return total_dispDU;
  }


////#############################################################################

tensor TwentyNodeBrick_u_p_U::total_dispDu(FILE *fp, double * u)
  {
    const int dimensions[] = {20,3};
    tensor total_dispDu(2, dimensions, 0.0);

    const Vector &TotDis1  = nd1Ptr->getTrialDisp();
    const Vector &TotDis2  = nd2Ptr->getTrialDisp();
    const Vector &TotDis3  = nd3Ptr->getTrialDisp();
    const Vector &TotDis4  = nd4Ptr->getTrialDisp();
    const Vector &TotDis5  = nd5Ptr->getTrialDisp();
    const Vector &TotDis6  = nd6Ptr->getTrialDisp();
    const Vector &TotDis7  = nd7Ptr->getTrialDisp();
    const Vector &TotDis8  = nd8Ptr->getTrialDisp();
    const Vector &TotDis9  = nd9Ptr->getTrialDisp();
    const Vector &TotDis10 = nd10Ptr->getTrialDisp();
    const Vector &TotDis11 = nd11Ptr->getTrialDisp();
    const Vector &TotDis12 = nd12Ptr->getTrialDisp();
    const Vector &TotDis13 = nd13Ptr->getTrialDisp();
    const Vector &TotDis14 = nd14Ptr->getTrialDisp();
    const Vector &TotDis15 = nd15Ptr->getTrialDisp();
    const Vector &TotDis16 = nd16Ptr->getTrialDisp();
    const Vector &TotDis17 = nd17Ptr->getTrialDisp();
    const Vector &TotDis18 = nd18Ptr->getTrialDisp();
    const Vector &TotDis19 = nd19Ptr->getTrialDisp();
    const Vector &TotDis20 = nd20Ptr->getTrialDisp();

// Get first three total displacement for solid part
    total_dispDu.val(1,1)=TotDis1(0); total_dispDu.val(1,2)=TotDis1(1);total_dispDu.val(1,3)=TotDis1(2);
//    total_disp.val(1,4)=TotDis1(3); total_disp.val(1,5)=TotDis1(4);total_disp.val(1,6)=TotDis1(5);
//    total_disp.val(1,7)=TotDis1(6);

    total_dispDu.val(2,1)=TotDis2(0); total_dispDu.val(2,2)=TotDis2(1);total_dispDu.val(2,3)=TotDis2(2);
//    total_disp.val(2,4)=TotDis2(3); total_disp.val(2,5)=TotDis2(4);total_disp.val(2,6)=TotDis2(5);
//    total_disp.val(2,7)=TotDis2(6);

    total_dispDu.val(3,1)=TotDis3(0); total_dispDu.val(3,2)=TotDis3(1);total_dispDu.val(3,3)=TotDis3(2);
//    total_disp.val(3,4)=TotDis3(3); total_disp.val(3,5)=TotDis3(4);total_disp.val(3,6)=TotDis3(5);
//    total_disp.val(3,7)=TotDis3(6);

    total_dispDu.val(4,1)=TotDis4(0); total_dispDu.val(4,2)=TotDis4(1);total_dispDu.val(4,3)=TotDis4(2);
//    total_disp.val(4,4)=TotDis4(3); total_disp.val(4,5)=TotDis4(4);total_disp.val(4,6)=TotDis4(5);
//    total_disp.val(4,7)=TotDis4(6);

    total_dispDu.val(5,1)=TotDis5(0); total_dispDu.val(5,2)=TotDis5(1);total_dispDu.val(5,3)=TotDis5(2);
//    total_disp.val(5,4)=TotDis5(3); total_disp.val(5,5)=TotDis5(4);total_disp.val(5,6)=TotDis5(5);
//    total_disp.val(5,7)=TotDis5(6);

    total_dispDu.val(6,1)=TotDis6(0); total_dispDu.val(6,2)=TotDis6(1);total_dispDu.val(6,3)=TotDis6(2);
//    total_disp.val(6,4)=TotDis6(3); total_disp.val(6,5)=TotDis6(4);total_disp.val(6,6)=TotDis6(5);
//    total_disp.val(6,7)=TotDis6(6);

    total_dispDu.val(7,1)=TotDis7(0); total_dispDu.val(7,2)=TotDis7(1);total_dispDu.val(7,3)=TotDis7(2);
//    total_disp.val(7,4)=TotDis7(3); total_disp.val(7,5)=TotDis7(4);total_disp.val(7,6)=TotDis7(5);
//    total_disp.val(7,7)=TotDis7(6);

    total_dispDu.val(8,1)=TotDis8(0); total_dispDu.val(8,2)=TotDis8(1);total_dispDu.val(8,3)=TotDis8(2);
//    total_disp.val(8,4)=TotDis8(3); total_disp.val(8,5)=TotDis8(4);total_disp.val(8,6)=TotDis8(5);
//    total_disp.val(8,7)=TotDis8(6);


    total_dispDu.val(9,1)=TotDis9(0); total_dispDu.val(9,2)=TotDis9(1);total_dispDu.val(9,3)=TotDis9(2);
//    total_disp.val(9,4)=TotDis9(3); total_disp.val(9,5)=TotDis9(4);total_disp.val(9,6)=TotDis9(5);
//    total_disp.val(9,7)=TotDis9(6);

    total_dispDu.val(10,1)=TotDis10(0); total_dispDu.val(10,2)=TotDis10(1);total_dispDu.val(10,3)=TotDis10(2);
//    total_disp.val(10,4)=TotDis10(3); total_disp.val(10,5)=TotDis10(4);total_disp.val(10,6)=TotDis10(5);
//    total_disp.val(10,7)=TotDis10(6);

    total_dispDu.val(11,1)=TotDis11(0); total_dispDu.val(11,2)=TotDis11(1);total_dispDu.val(11,3)=TotDis11(2);
//    total_disp.val(11,4)=TotDis11(3); total_disp.val(11,5)=TotDis11(4);total_disp.val(11,6)=TotDis11(5);
//    total_disp.val(11,7)=TotDis11(6);

    total_dispDu.val(12,1)=TotDis12(0); total_dispDu.val(12,2)=TotDis12(1);total_dispDu.val(12,3)=TotDis12(2);
//    total_disp.val(12,4)=TotDis12(3); total_disp.val(12,5)=TotDis12(4);total_disp.val(12,6)=TotDis12(5);
//    total_disp.val(12,7)=TotDis12(6);


    total_dispDu.val(13,1)=TotDis13(0); total_dispDu.val(13,2)=TotDis13(1);total_dispDu.val(13,3)=TotDis13(2);
//    total_disp.val(13,4)=TotDis13(3); total_disp.val(13,5)=TotDis13(4);total_disp.val(13,6)=TotDis13(5);
//    total_disp.val(13,7)=TotDis13(6);

    total_dispDu.val(14,1)=TotDis14(0); total_dispDu.val(14,2)=TotDis14(1);total_dispDu.val(14,3)=TotDis14(2);
//    total_disp.val(14,4)=TotDis14(3); total_disp.val(14,5)=TotDis14(4);total_disp.val(14,6)=TotDis14(5);
//    total_disp.val(14,7)=TotDis14(6);

    total_dispDu.val(15,1)=TotDis15(0); total_dispDu.val(15,2)=TotDis15(1);total_dispDu.val(15,3)=TotDis15(2);
//    total_disp.val(15,4)=TotDis15(3); total_disp.val(15,5)=TotDis15(4);total_disp.val(15,6)=TotDis15(5);
//    total_disp.val(15,7)=TotDis15(6);

    total_dispDu.val(16,1)=TotDis16(0); total_dispDu.val(16,2)=TotDis16(1);total_dispDu.val(16,3)=TotDis16(2);
//    total_disp.val(16,4)=TotDis16(3); total_disp.val(16,5)=TotDis16(4);total_disp.val(16,6)=TotDis16(5);
//    total_disp.val(16,7)=TotDis16(6);

    total_dispDu.val(17,1)=TotDis17(0); total_dispDu.val(17,2)=TotDis17(1);total_dispDu.val(17,3)=TotDis17(2);
//    total_disp.val(17,4)=TotDis17(3); total_disp.val(17,5)=TotDis17(4);total_disp.val(17,6)=TotDis17(5);
//    total_disp.val(17,7)=TotDis17(6);

    total_dispDu.val(18,1)=TotDis18(0); total_dispDu.val(18,2)=TotDis18(1);total_dispDu.val(18,3)=TotDis18(2);
//    total_disp.val(18,4)=TotDis18(3); total_disp.val(18,5)=TotDis18(4);total_disp.val(18,6)=TotDis18(5);
//    total_disp.val(18,7)=TotDis18(6);

    total_dispDu.val(19,1)=TotDis19(0); total_dispDu.val(19,2)=TotDis19(1);total_dispDu.val(19,3)=TotDis19(2);
//    total_disp.val(19,4)=TotDis19(3); total_disp.val(19,5)=TotDis19(4);total_disp.val(19,6)=TotDis19(5);
//    total_disp.val(19,7)=TotDis19(6);

    total_dispDu.val(20,1)=TotDis20(0); total_dispDu.val(20,2)=TotDis20(1);total_dispDu.val(20,3)=TotDis20(2);
//    total_disp.val(20,4)=TotDis20(3); total_disp.val(20,5)=TotDis20(4);total_disp.val(20,6)=TotDis20(5);
//    total_disp.val(20,7)=TotDis20(6);


    return total_dispDu;

  }
////#############################################################################

tensor TwentyNodeBrick_u_p_U::total_dispDU(FILE *fp, double * u)
  {
    const int dimensions[] = {20,3};
    tensor total_dispDU(2, dimensions, 0.0);

    const Vector &TotDis1  = nd1Ptr->getTrialDisp();
    const Vector &TotDis2  = nd2Ptr->getTrialDisp();
    const Vector &TotDis3  = nd3Ptr->getTrialDisp();
    const Vector &TotDis4  = nd4Ptr->getTrialDisp();
    const Vector &TotDis5  = nd5Ptr->getTrialDisp();
    const Vector &TotDis6  = nd6Ptr->getTrialDisp();
    const Vector &TotDis7  = nd7Ptr->getTrialDisp();
    const Vector &TotDis8  = nd8Ptr->getTrialDisp();
    const Vector &TotDis9  = nd9Ptr->getTrialDisp();
    const Vector &TotDis10 = nd10Ptr->getTrialDisp();
    const Vector &TotDis11 = nd11Ptr->getTrialDisp();
    const Vector &TotDis12 = nd12Ptr->getTrialDisp();
    const Vector &TotDis13 = nd13Ptr->getTrialDisp();
    const Vector &TotDis14 = nd14Ptr->getTrialDisp();
    const Vector &TotDis15 = nd15Ptr->getTrialDisp();
    const Vector &TotDis16 = nd16Ptr->getTrialDisp();
    const Vector &TotDis17 = nd17Ptr->getTrialDisp();
    const Vector &TotDis18 = nd18Ptr->getTrialDisp();
    const Vector &TotDis19 = nd19Ptr->getTrialDisp();
    const Vector &TotDis20 = nd20Ptr->getTrialDisp();

//  Get last three total displacement for fluid part.
    total_dispDU.val(1,1)=TotDis1(4);total_dispDU.val(1,2)=TotDis1(5);
    total_dispDU.val(1,3)=TotDis1(6);

    total_dispDU.val(2,1)=TotDis2(4);total_dispDU.val(2,2)=TotDis2(5);
    total_dispDU.val(2,3)=TotDis2(6);

    total_dispDU.val(3,1)=TotDis3(4);total_dispDU.val(3,2)=TotDis3(5);
    total_dispDU.val(3,3)=TotDis3(6);

    total_dispDU.val(4,1)=TotDis4(4);total_dispDU.val(4,2)=TotDis4(5);
    total_dispDU.val(4,3)=TotDis4(6);

    total_dispDU.val(5,1)=TotDis5(4);total_dispDU.val(5,2)=TotDis5(5);
    total_dispDU.val(5,3)=TotDis5(6);

    total_dispDU.val(6,1)=TotDis6(4);total_dispDU.val(6,2)=TotDis6(5);
    total_dispDU.val(6,3)=TotDis6(6);

    total_dispDU.val(7,1)=TotDis7(4);total_dispDU.val(7,2)=TotDis7(5);
    total_dispDU.val(7,3)=TotDis7(6);

    total_dispDU.val(8,1)=TotDis8(4);total_dispDU.val(8,2)=TotDis8(5);
    total_dispDU.val(8,3)=TotDis8(6);


    total_dispDU.val(9,1)=TotDis9(4);total_dispDU.val(9,2)=TotDis9(5);
    total_dispDU.val(9,3)=TotDis9(6);

    total_dispDU.val(10,1)=TotDis10(4);total_dispDU.val(10,2)=TotDis10(5);
    total_dispDU.val(10,3)=TotDis10(6);

    total_dispDU.val(11,1)=TotDis11(4);total_dispDU.val(11,2)=TotDis11(5);
    total_dispDU.val(11,3)=TotDis11(6);

    total_dispDU.val(12,1)=TotDis12(4);total_dispDU.val(12,2)=TotDis12(5);
    total_dispDU.val(12,3)=TotDis12(6);


    total_dispDU.val(13,1)=TotDis13(4);total_dispDU.val(13,2)=TotDis13(5);
    total_dispDU.val(13,3)=TotDis13(6);

    total_dispDU.val(14,1)=TotDis14(4);total_dispDU.val(14,2)=TotDis14(5);
    total_dispDU.val(14,3)=TotDis14(6);

    total_dispDU.val(15,1)=TotDis15(4);total_dispDU.val(15,2)=TotDis15(5);
    total_dispDU.val(15,3)=TotDis15(6);

    total_dispDU.val(16,1)=TotDis16(4);total_dispDU.val(16,2)=TotDis16(5);
    total_dispDU.val(16,3)=TotDis16(6);

    total_dispDU.val(17,1)=TotDis17(4);total_dispDU.val(17,2)=TotDis17(5);
    total_dispDU.val(17,3)=TotDis17(6);

    total_dispDU.val(18,1)=TotDis18(4);total_dispDU.val(18,2)=TotDis18(5);
    total_dispDU.val(18,3)=TotDis18(6);

    total_dispDU.val(19,1)=TotDis19(4);total_dispDU.val(19,2)=TotDis19(5);
    total_dispDU.val(19,3)=TotDis19(6);

    total_dispDU.val(20,1)=TotDis20(4);total_dispDU.val(20,2)=TotDis20(5);
    total_dispDU.val(20,3)=TotDis20(6);


    return total_dispDU;
}

//====================================================================
// This function is not used. Just keep here and add another function
// incremental_UpdateDU for furture use. // Xiaoyan 02/08/2002
void TwentyNodeBrick_u_p_U::incremental_UpdateDu()
  {
    double r  = 0.0;
    // double rw = 0.0;
    double s  = 0.0;
    // double sw = 0.0;
    double t  = 0.0;
    // double tw = 0.0;

    short where = 0;
    //,,,,,    double weight = 0.0;

    //double this_one_PP = (matpoint)->operator[](where).IS_Perfect_Plastic();

    int dh_dim[] = {20,3};
    tensor dh(2, dh_dim, 0.0);

//    tensor Constitutive( 4, def_dim_4, 0.0);

    //    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {20,3};
    tensor incremental_displacements(2,disp_dim,0.0);

    straintensor incremental_strain;
//    straintensor total_strain_at_GP;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    //....    int number_of_subincrements = 1;
    //....    double this_one_PP = 1.0; // if set to 0.0 -> perfectly plastic
    //....                              // if set to 1.0 -> elasto plastic

//    stresstensor final_stress_after_integration;

    ///    stresstensor incremental_stress;
    // tensor of incremental displacements taken from node objects
    incremental_displacements = incr_dispDu();

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        //--        rw = get_Gauss_p_w( r_integration_order, GP_c_r );
        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            //--            sw = get_Gauss_p_w( s_integration_order, GP_c_s );
            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
            {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                //--                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                   ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordiantes with respect to local coordiantes
                dh = dh_drst_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //....                Jacobian.print("J");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");
                // determinant of Jacobian tensor ( matrix )
                //--                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");
                //....                dhGlobal.print("dh","dhGlobal");
                //weight
                //                weight = rw * sw * tw * det_of_Jacobian;
                //::::::   ::printf("\n\nIN THE STIFFNESS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                //::::::   ::printf(" void TwentyNodeBrick_u_p_U::incremental_Update()\n");
                //::::::   ::printf(" GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d    --->>>  where = %d \n",
                //::::::                      GP_c_r,GP_c_s,GP_c_t,where);
                //::::::   ::printf("WEIGHT = %f", weight);
                //::::::   ::printf("determinant of Jacobian = %f", determinant_of_Jacobian);
                //::::::   matpoint[where].report("Gauss Point\n");
                // incremental straines at this Gauss point
                // now in Update we know the incremental displacements so let's find
                // the incremental strain
                incremental_strain =
                    (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
                incremental_strain.null_indices();
                //incremental_strain.reportshort("\n incremental_strain tensor at GAUSS point\n");

                // here comes the final_stress calculation actually on only needs to copy stresses
                // from the iterative data . . .
                //(GPstress+where)->reportshortpqtheta("\n stress START GAUSS \n");

    // Getting final_stress_after_integration is  Done inside CDriver on EPState____ZHaohui
    //final_stress_after_integration = GPiterative_stress[where];
    //(matpoint)->operator[](where).kappa_set(final_stress_after_integration,
                //                                 GPq_ast_iterative[where]);

    //....         final_stress_after_integration =
                //....           (matpoint)->operator[](where).FinalStress(*(GPstress+where),
                //....                                                     incremental_strain,
                //....                                                     (matpoint)->operator[](where),
                //....                                                     number_of_subincrements,
                //....                                                     this_one_PP);
                //....//final_stress_after_integration.reportshortpqtheta("\n final_stress_after_integration GAUSS \n");
                // calculate the constitutive tensor

                // We do not need: final_stress_after_integration


          //Constitutive =
                //  (matpoint)->operator[](where).ConstitutiveTensor(final_stress_after_integration,
                //                                                   *(GPstress+where),
                //                                                   incremental_strain,
                //                                                   (matpoint)->operator[](where),
                //                                                   this_one_PP);

    // ZHaohui modified __09-29-2000

          // Now no EPState  but NDMaterial for each MatPoint
    //EPState *tmp_eps = (matpoint[where]).getEPS();
          //NDMaterial *tmp_ndm = (matpoint[where]).getNDMat();

    //if ( tmp_eps ) { //if there is an EPState for the MatPoint3D
    //  mmodel->setEPS( *tmp_eps );

    if ( ! ( (matpoint[where]->matmodel)->setTrialStrainIncr( incremental_strain)) )
      opserr << "TwentyNodeBrick_u_p_U::incremental_Update (tag: " << this->getTag() << "), not converged\n";

    //matpoint[where].setEPS( mmodel->getEPS() );
    //}

    //else if ( tmp_ndm )
    //  (matpoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
    //else {
                 //   g3ErrorHandler->fatal("TwentyNodeBrick_u_p_U::incremental_Update (tag: %d), no strain or stress state vars", this->getTag());
    //   exit(1);
    //}

    //Constitutive = trialEPS.getEep();

    //::::::                   Constitutive.print("C","\n\n C tensor \n");
          // this is update of constitutive tensor at this Gauss point

    // All done in EPState when calling setTrialStrainIncr and converged
    //GPtangent_E[where].Initialize(Constitutive);
          //GPtangent_E[where].print("\n tangent E at GAUSS point \n");

                //total_strain_at_GP.Initialize(*(GPstrain+where));
                //total_strain_at_GP.reportshort("\n total_strain tensor at GAUSS point \n");
                //total_strain_at_GP = total_strain_at_GP + incremental_strain;
                //total_strain_at_GP.reportshort("\n total_strain tensor at GAUSS point AFTER\n");
                //GPstress[where].Initialize(final_stress_after_integration);
                //GPstress[where].reportshortpqtheta("\n stress at GAUSS point \n");

    //GPstrain[where].Initialize(total_strain_at_GP);

    //GPstrain[where].reportshort("\n strain at GAUSS point \n");
            }
          }
      }
  }

//============================================================================

void TwentyNodeBrick_u_p_U::incremental_UpdateDU()
  {
    double r  = 0.0;
    // double rw = 0.0;
    double s  = 0.0;
    // double sw = 0.0;
    double t  = 0.0;
    // double tw = 0.0;

    short where = 0;
    //,,,,,    double weight = 0.0;

    //double this_one_PP = (matpoint)->operator[](where).IS_Perfect_Plastic();

    int dh_dim[] = {20,3};
    tensor dh(2, dh_dim, 0.0);

//    tensor Constitutive( 4, def_dim_4, 0.0);

    //    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {20,3};
    tensor incremental_displacements(2,disp_dim,0.0);

    straintensor incremental_strain;
//    straintensor total_strain_at_GP;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    //....    int number_of_subincrements = 1;
    //....    double this_one_PP = 1.0; // if set to 0.0 -> perfectly plastic
    //....                              // if set to 1.0 -> elasto plastic

//    stresstensor final_stress_after_integration;

    ///    stresstensor incremental_stress;
    // tensor of incremental displacements taken from node objects
    incremental_displacements = incr_dispDU();

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        //--        rw = get_Gauss_p_w( r_integration_order, GP_c_r );
        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            //--            sw = get_Gauss_p_w( s_integration_order, GP_c_s );
            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
            {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                //--                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                   ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordiantes with respect to local coordiantes
                dh = dh_drst_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //....                Jacobian.print("J");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");
                // determinant of Jacobian tensor ( matrix )
                //--                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");
                //....                dhGlobal.print("dh","dhGlobal");
                //weight
                //                weight = rw * sw * tw * det_of_Jacobian;
                //::::::   ::printf("\n\nIN THE STIFFNESS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                //::::::   ::printf(" void TwentyNodeBrick_u_p_U::incremental_Update()\n");
                //::::::   ::printf(" GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d    --->>>  where = %d \n",
                //::::::                      GP_c_r,GP_c_s,GP_c_t,where);
                //::::::   ::printf("WEIGHT = %f", weight);
                //::::::   ::printf("determinant of Jacobian = %f", determinant_of_Jacobian);
                //::::::   matpoint[where].report("Gauss Point\n");
                // incremental straines at this Gauss point
                // now in Update we know the incremental displacements so let's find
                // the incremental strain
                incremental_strain =
                    (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
                incremental_strain.null_indices();
                //incremental_strain.reportshort("\n incremental_strain tensor at GAUSS point\n");

                // here comes the final_stress calculation actually on only needs to copy stresses
                // from the iterative data . . .
                //(GPstress+where)->reportshortpqtheta("\n stress START GAUSS \n");

    // Getting final_stress_after_integration is  Done inside CDriver on EPState____ZHaohui
    //final_stress_after_integration = GPiterative_stress[where];
    //(matpoint)->operator[](where).kappa_set(final_stress_after_integration,
                //                                 GPq_ast_iterative[where]);

    //....         final_stress_after_integration =
                //....           (matpoint)->operator[](where).FinalStress(*(GPstress+where),
                //....                                                     incremental_strain,
                //....                                                     (matpoint)->operator[](where),
                //....                                                     number_of_subincrements,
                //....                                                     this_one_PP);
                //....//final_stress_after_integration.reportshortpqtheta("\n final_stress_after_integration GAUSS \n");
                // calculate the constitutive tensor

                // We do not need: final_stress_after_integration


          //Constitutive =
                //  (matpoint)->operator[](where).ConstitutiveTensor(final_stress_after_integration,
                //                                                   *(GPstress+where),
                //                                                   incremental_strain,
                //                                                   (matpoint)->operator[](where),
                //                                                   this_one_PP);

    // ZHaohui modified __09-29-2000

          // Now no EPState  but NDMaterial for each MatPoint
    //EPState *tmp_eps = (matpoint[where]).getEPS();
          //NDMaterial *tmp_ndm = (matpoint[where]).getNDMat();

    //if ( tmp_eps ) { //if there is an EPState for the MatPoint3D
    //  mmodel->setEPS( *tmp_eps );

    if ( ! ( (matpoint[where]->matmodel)->setTrialStrainIncr( incremental_strain)) )
      opserr << "TwentyNodeBrick_u_p_U::incremental_Update (tag: " << this->getTag() << "), not converged\n";

    //matpoint[where].setEPS( mmodel->getEPS() );
    //}

    //else if ( tmp_ndm )
    //  (matpoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
    //else {
                 //   g3ErrorHandler->fatal("TwentyNodeBrick_u_p_U::incremental_Update (tag: %d), no strain or stress state vars", this->getTag());
    //   exit(1);
    //}

    //Constitutive = trialEPS.getEep();

    //::::::                   Constitutive.print("C","\n\n C tensor \n");
          // this is update of constitutive tensor at this Gauss point

    // All done in EPState when calling setTrialStrainIncr and converged
    //GPtangent_E[where].Initialize(Constitutive);
          //GPtangent_E[where].print("\n tangent E at GAUSS point \n");

                //total_strain_at_GP.Initialize(*(GPstrain+where));
                //total_strain_at_GP.reportshort("\n total_strain tensor at GAUSS point \n");
                //total_strain_at_GP = total_strain_at_GP + incremental_strain;
                //total_strain_at_GP.reportshort("\n total_strain tensor at GAUSS point AFTER\n");
                //GPstress[where].Initialize(final_stress_after_integration);
                //GPstress[where].reportshortpqtheta("\n stress at GAUSS point \n");

    //GPstrain[where].Initialize(total_strain_at_GP);

    //GPstrain[where].reportshort("\n strain at GAUSS point \n");
            }
          }
      }
  }

//#############################################################################
// This function is not called now. Keep here for future use.  Xiaoyan 02/108/2002
void TwentyNodeBrick_u_p_U::set_strain_stress_tensorDu(FILE *fp, double * u)
  {
    int dh_dim[] = {20,3};
    tensor dh(2, dh_dim, 0.0);

//    tensor Constitutive( 4, def_dim_4, 0.0);
    tensor Constitutive;
    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;
    int where = 0;

    double det_of_Jacobian;

    straintensor strain;
    stresstensor stress;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;


    static int disp_dim[] = {20,3};
    tensor total_displacements(2,disp_dim,0.0); //

    total_displacements = total_dispDu(fp, u);

    ::printf("\n  displacement(x-y-z) at GAUSS pt %d \n\n", where+1);
    for (int ii=1; ii<=8;ii++)
     {
      ::printf("Global# %d Local#%d  %+0.5e %+0.5e %+0.5e\n",
                     //G_N_numbs[ii-1],
         connectedExternalNodes(ii-1),
         ii,total_displacements.val(ii,1),
               total_displacements.val(ii,2),
         total_displacements.val(ii,3));
     }
    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");
                //weight
                // straines at this Gauss point from displacement
                strain = (dhGlobal("ib")*total_displacements("ia")).symmetrize11();
                strain.null_indices();
                // here comes the final_stress calculation
                // at this Gauss point.

                //Constitutive =  GPtangent_E[where];
                //Constitutive =  (matpoint->getEPS() )->getEep();
                // if set total displ, then it should be elstic material
    Constitutive =  ( matpoint[where]->matmodel)->getTangentTensor();

           stress = Constitutive("ijkl") * strain("kl");   //<<<<<<<<<<<<<<<
                stress.null_indices();

                ::printf("\n  strain tensor at GAUSS point %d \n", where+1);
                strain.reportshort("");
                ::printf("\n  stress tensor at GAUSS point %d \n", where+1);
                stress.reportshort("");


              }
          }
      }
  }


////#############################################################################

//#############################################################################
// This function is not called now. Keep here for future use.  Xiaoyan 02/108/2002
void TwentyNodeBrick_u_p_U::set_strain_stress_tensorDU(FILE *fp, double * u)
  {
    int dh_dim[] = {20,3};
    tensor dh(2, dh_dim, 0.0);

//    tensor Constitutive( 4, def_dim_4, 0.0);
    tensor Constitutive;
    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;
    int where = 0;

    double det_of_Jacobian;

    straintensor strain;
    stresstensor stress;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;


    static int disp_dim[] = {20,3};
    tensor total_displacements(2,disp_dim,0.0); //

    total_displacements = total_dispDU(fp, u);

    ::printf("\n  displacement(x-y-z) at GAUSS pt %d \n\n", where+1);
    for (int ii=1; ii<=8;ii++)
     {
      ::printf("Global# %d Local#%d  %+0.5e %+0.5e %+0.5e\n",
                     //G_N_numbs[ii-1],
         connectedExternalNodes(ii-1),
         ii,total_displacements.val(ii,1),
               total_displacements.val(ii,2),
         total_displacements.val(ii,3));
     }
    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");
                //weight
                // straines at this Gauss point from displacement
                strain = (dhGlobal("ib")*total_displacements("ia")).symmetrize11();
                strain.null_indices();
                // here comes the final_stress calculation
                // at this Gauss point.

                //Constitutive =  GPtangent_E[where];
                //Constitutive =  (matpoint->getEPS() )->getEep();
                // if set total displ, then it should be elstic material
    Constitutive =  ( matpoint[where]->matmodel)->getTangentTensor();

           stress = Constitutive("ijkl") * strain("kl");   //<<<<<<<<<<<<<<<
                stress.null_indices();

                ::printf("\n  strain tensor at GAUSS point %d \n", where+1);
                strain.reportshort("");
                ::printf("\n  stress tensor at GAUSS point %d \n", where+1);
                stress.reportshort("");


              }
          }
      }
  }


////#############################################################################

////#############################################################################
TwentyNodeBrick_u_p_U & TwentyNodeBrick_u_p_U::operator[](int subscript)
  {
    return ( *(this+subscript) );
  }

//Finite_Element & TwentyNodeBrick_u_p_U::operator[](short subscript)
//  {
//    return ( *(this+subscript) );
//  }

//Finite_Element & TwentyNodeBrick_u_p_U::operator[](unsigned subscript)
//  {
//    return ( *(this+subscript) );
//  }


////#############################################################################
int TwentyNodeBrick_u_p_U::commitState ()
{
    // int order = theQuadRule->getOrder();     // Commented by Xiaoyan

    int i;
    //int j, k;      // Xiaoyan added k for three dimension
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "TwentyNodeBrick_u_p_U::commitState () - failed in base class";
    }

    // Loop over the integration points and commit the material states
    int count  = r_integration_order* s_integration_order * t_integration_order;

    //for (i = 0; i < r_integration_order; i++)        // Xiaoyan chaneged order to
    //  for (j = 0; j < s_integration_order; j++)      // r_integration_order,
    //                  // s_integration_order, and
    //      for (k = 0; k < t_integration_order; k++)      // added t_integration_order,
    //         retVal += (GaussPtheMaterial[i][j][k]).commitState();

    Vector pp = getResistingForce();

    //if ( this->getTag() == 1 || this->getTag() == 700)
    //{
      for (i = 0; i < count; i++)  // use this instead of the below wxy 09/28/2001
      //for (i = 0; i < 8; i++)
      {
         retVal += matpoint[i]->commitState();
         //if (i == 4 && strcmp(matpoint[i]->matmodel->getType(),"Template3Dep") == 0)
         stresstensor st;
   stresstensor prin;
         straintensor stn;
         straintensor stnprin;

         st = matpoint[i]->getStressTensor();
          prin = st.principal();
         stn = matpoint[i]->getStrainTensor();
          stnprin = stn.principal();
         /*
   opserr << "\nGauss Point: " << i << endln;
   opserr << "sigma11: "<< st.cval(1, 1) << " "<< st.cval(1, 2) << " " << st.cval(1, 3) << endln;
   opserr << "sigma21: "<< st.cval(2, 1) << " "<< st.cval(2, 2) << " " << st.cval(2, 3) << endln;
    opserr << "sigma31: "<< st.cval(3, 1) << " "<< st.cval(3, 2) << " " << st.cval(3, 3) << endln << endln;
   */
   //opserr << "strain11: "<< stn.cval(1, 1) << " "<< stn.cval(1, 2) << " " << stn.cval(1, 3) << endln;
   //opserr << "strain21: "<< stn.cval(2, 1) << " "<< stn.cval(2, 2) << " " << stn.cval(2, 3) << endln;
    //opserr << "strain31: "<< stn.cval(3, 1) << " "<< stn.cval(3, 2) << " " << stn.cval(3, 3) << endln;

//   double  p = -1*( prin.cval(1, 1)+ prin.cval(2, 2) +prin.cval(3, 3) )/3.0;          commented by Xiaoyan
//   double  ev = -1*( stnprin.cval(1, 1)+ stnprin.cval(2, 2) + stnprin.cval(3, 3) )/3.0;   commented by Xiaoyan
                                                                                                // 09/24/2001.
                         // These are unused variable

   //opserr << "   " << p;

   //if (p < 0)
   //  opserr  << "gs pnt:" << i << "  p="<< p;


   double q;
   //if ( fabs(prin.cval(1, 1) - prin.cval(2, 2) ) <=  0.0001 )
         if ( fabs(prin.cval(1, 1) - prin.cval(2, 2) ) <=  0.001 )
         {
             q = prin.cval(1, 1) - prin.cval(3, 3);
             //opserr << "1 = 2";
         }
         else
             q = prin.cval(3, 3) - prin.cval(1, 1);

   //Triaxial compr.  fabs
         //opserr << "     " << st.cval(2, 3); //tau_yz
   //opserr << "     " << q;
   ////----opserr << "     " << fabs(q);

         //opserr << "     " << ev << endln;

//out22Jan2001   if (strcmp(matpoint[i]->matmodel->getType(),"Template3Dep") == 0)
//out22Jan2001          {
//out22Jan2001           st = ( ((Template3Dep *)(matpoint[i]->matmodel))->getEPS())->getStress();
//out22Jan2001           prin = st.principal();
//out22Jan2001    }
//out22Jan2001    else
//out22Jan2001    {
//out22Jan2001            st = matpoint[i]->getStressTensor();
//out22Jan2001           prin = st.principal();
//out22Jan2001
//out22Jan2001    }

    //double  p = st.p_hydrostatic();
    //double  p = -1*( prin.cval(1, 1)+ prin.cval(2, 2) +prin.cval(3, 3) )/3.0;
          //opserr << "\n " << prin.cval(1, 1) << "   " << prin.cval(2, 2) << "  " <<  prin.cval(3, 3) << endln;
          //if ( getTag() == 960)
          //opserr << " El= " << getTag() << " , p    " << p << endln;

    //printf(stderr, " Gauss Point i = %d ", (i+1));
    //printf(stderr, " Gauss Point i = %d ", (i+1));


          //if ( p < 0 )
    //{
    //  opserr << getTag();
    //  opserr << " ***p  =    " << p << endln;
    //}
          //J2D
          //opserr << "        " << st.q_deviatoric();

          //double q;
          //if ( fabs(prin.cval(1, 1) - prin.cval(2, 2) ) <=  0.0001 )
          //{
          //    q = prin.cval(1, 1) - prin.cval(3, 3);
          //    //opserr << "1 = 2";
          //}
          //else
          //    q = prin.cval(3, 3) - prin.cval(1, 1);

          //Triaxial compr.
          //opserr << "        " << q;
         //}
      }

      //opserr << " at elements " << this->getTag() << endln;


      //output nodal force
      //opserr << "    " << pp(2) << endln;
    //}
    return retVal;
}

////#############################################################################
int TwentyNodeBrick_u_p_U::revertToLastCommit ()
{
  //  int order = theQuadRule->getOrder();  // Commented by Xiaoyan
    int i;
    //int j, k;     // Xiaoyan added k for three dimension
    int retVal = 0;

    // Loop over the integration points and revert to last committed material states
    int count  = r_integration_order* s_integration_order * t_integration_order;
    //for (i = 0; i < r_integration_order; i++)       // Xiaoyan chaneged order to
    //  for (j = 0; j < s_integration_order; j++)     // r_integration_order,
    //      for (k = 0; k < t_integration_order; k++)     // s_integration_order, and
                       // added t_integration_order,
      //retVal += (theMaterial[i][j][k]).revertToLastCommit();

    for (i = 0; i < count; i++)
       retVal += matpoint[i]->revertToLastCommit();


    return retVal;
}

////#############################################################################
//=============================================================================
int TwentyNodeBrick_u_p_U::revertToStart ()
{
    int i;     // Xiaoyan added k for three dimension
    int retVal = 0;

    // Loop over the integration points and revert to last committed material states
    //for (i = 0; i < r_integration_order; i++)       // Xiaoyan chaneged order to
    //  for (j = 0; j < s_integration_order; j++)     // r_integration_order,
    //      for (k = 0; k < t_integration_order; k++)     // s_integration_order, and
                 // added t_integration_order,
    //      retVal += (theMaterial[i][j][k]).revertToLastCommit();

    int count  = r_integration_order* s_integration_order * t_integration_order;

    for (i = 0; i < count; i++)
       retVal += matpoint[i]->revertToStart();


    return retVal;

    // Loop over the integration points and revert to initial material states
   }
////#############################################################################

//=============================================================================
//  The following are come from FourNodeQuad.cc   Xiaoyan 07/06/00
//  The following are come from FourNodeQuad.cc   Xiaoyan 07/06/00
//  The following are come from FourNodeQuad.cc   Xiaoyan 07/06/00
//=============================================================================

//=============================================================================
int TwentyNodeBrick_u_p_U::getNumExternalNodes () const
{
    return nodes_in_brick;  //changed from 4 to 8 Xiaoyan 07/06/00
}


//=============================================================================
const ID& TwentyNodeBrick_u_p_U::getExternalNodes ()
{
    return connectedExternalNodes;
}

Node **
TwentyNodeBrick_u_p_U::getNodePtrs(void)
{
  theNodes[0] = nd1Ptr;
  theNodes[1] = nd2Ptr;
  theNodes[2] = nd3Ptr;
  theNodes[3] = nd4Ptr;
  theNodes[4] = nd5Ptr;
  theNodes[5] = nd6Ptr;
  theNodes[6] = nd7Ptr;
  theNodes[7] = nd8Ptr;
  theNodes[8] = nd9Ptr;
  theNodes[9] = nd10Ptr;
  theNodes[10] = nd11Ptr;
  theNodes[11] = nd12Ptr;
  theNodes[12] = nd13Ptr;
  theNodes[13] = nd14Ptr;
  theNodes[14] = nd15Ptr;
  theNodes[15] = nd16Ptr;
  theNodes[16] = nd17Ptr;
  theNodes[17] = nd18Ptr;
  theNodes[18] = nd19Ptr;
  theNodes[19] = nd20Ptr;

  return theNodes;
}

//=============================================================================
int TwentyNodeBrick_u_p_U::getNumDOF ()
{
    return 7*nodes_in_brick;       // 3*20=60 Xiaoyan 09/28/2001
}

//=============================================================================
void TwentyNodeBrick_u_p_U::setDomain (Domain *theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
  nd1Ptr = 0;
  nd2Ptr = 0;
  nd3Ptr = 0;
  nd4Ptr = 0;
  //Xiaoyan added 5-8  07/06/00
  nd5Ptr = 0;
  nd6Ptr = 0;
  nd7Ptr = 0;
  nd8Ptr = 0;

  nd9Ptr  = 0;
        nd10Ptr = 0;
        nd11Ptr = 0;
        nd12Ptr = 0;

        nd13Ptr = 0;
        nd14Ptr = 0;
        nd15Ptr = 0;
        nd16Ptr = 0;

        nd17Ptr = 0;
        nd18Ptr = 0;
        nd19Ptr = 0;
        nd20Ptr = 0;
    }
    //Added if-else for found a bug when trying removeElement from theDomain  07-19-2001 Zhaohui
    else {
      int Nd1 = connectedExternalNodes(0);
      int Nd2 = connectedExternalNodes(1);
      int Nd3 = connectedExternalNodes(2);
      int Nd4 = connectedExternalNodes(3);
      //Xiaoyan added 5-8  07/06/00

      int Nd5 = connectedExternalNodes(4);
      int Nd6 = connectedExternalNodes(5);
      int Nd7 = connectedExternalNodes(6);
      int Nd8 = connectedExternalNodes(7);

      int Nd9  = connectedExternalNodes( 8);
      int Nd10 = connectedExternalNodes( 9);
      int Nd11 = connectedExternalNodes(10);
      int Nd12 = connectedExternalNodes(11);

      int Nd13 = connectedExternalNodes(12);
      int Nd14 = connectedExternalNodes(13);
      int Nd15 = connectedExternalNodes(14);
      int Nd16 = connectedExternalNodes(15);

      int Nd17 = connectedExternalNodes(16);
      int Nd18 = connectedExternalNodes(17);
      int Nd19 = connectedExternalNodes(18);
      int Nd20 = connectedExternalNodes(19);

      nd1Ptr = theDomain->getNode(Nd1);
      nd2Ptr = theDomain->getNode(Nd2);
      nd3Ptr = theDomain->getNode(Nd3);
      nd4Ptr = theDomain->getNode(Nd4);

      //Xiaoyan added 5-8  07/06/00
      nd5Ptr = theDomain->getNode(Nd5);
      nd6Ptr = theDomain->getNode(Nd6);
      nd7Ptr = theDomain->getNode(Nd7);
      nd8Ptr = theDomain->getNode(Nd8);

      nd9Ptr = theDomain->getNode(Nd9);
      nd10Ptr = theDomain->getNode(Nd10);
      nd11Ptr = theDomain->getNode(Nd11);
      nd12Ptr = theDomain->getNode(Nd12);

      nd13Ptr = theDomain->getNode(Nd13);
      nd14Ptr = theDomain->getNode(Nd14);
      nd15Ptr = theDomain->getNode(Nd15);
      nd16Ptr = theDomain->getNode(Nd16);

      nd17Ptr = theDomain->getNode(Nd17);
      nd18Ptr = theDomain->getNode(Nd18);
      nd19Ptr = theDomain->getNode(Nd19);
      nd20Ptr = theDomain->getNode(Nd20);

      if (nd1Ptr  == 0 || nd2Ptr  == 0 || nd3Ptr  == 0 || nd4Ptr  == 0 ||
          nd5Ptr  == 0 || nd6Ptr  == 0 || nd7Ptr  == 0 || nd8Ptr  == 0 ||
    nd9Ptr  == 0 || nd10Ptr == 0 || nd11Ptr == 0 || nd12Ptr == 0 ||
          nd13Ptr == 0 || nd14Ptr == 0 || nd15Ptr == 0 || nd16Ptr == 0 ||
          nd17Ptr == 0 || nd18Ptr == 0 || nd19Ptr == 0 || nd20Ptr == 0 ) {

        opserr << "FATAL ERROR TwentyNodeBrick_u_p_U (tag: " << this->getTag() << "), node not found in domain\n";
        exit(-1);
      }

      int dofNd1 = nd1Ptr->getNumberDOF();
      int dofNd2 = nd2Ptr->getNumberDOF();
      int dofNd3 = nd3Ptr->getNumberDOF();
      int dofNd4 = nd4Ptr->getNumberDOF();

      int dofNd5 = nd5Ptr->getNumberDOF();
      int dofNd6 = nd6Ptr->getNumberDOF();
      int dofNd7 = nd7Ptr->getNumberDOF();
      int dofNd8 = nd8Ptr->getNumberDOF();

      int dofNd9 = nd9Ptr->getNumberDOF();
      int dofNd10 = nd10Ptr->getNumberDOF();
      int dofNd11 = nd11Ptr->getNumberDOF();
      int dofNd12 = nd12Ptr->getNumberDOF();

      int dofNd13 = nd13Ptr->getNumberDOF();
      int dofNd14 = nd14Ptr->getNumberDOF();
      int dofNd15 = nd15Ptr->getNumberDOF();
      int dofNd16 = nd16Ptr->getNumberDOF();

      int dofNd17 = nd17Ptr->getNumberDOF();
      int dofNd18 = nd18Ptr->getNumberDOF();
      int dofNd19 = nd19Ptr->getNumberDOF();
      int dofNd20 = nd20Ptr->getNumberDOF();

      if (dofNd1  != 7 || dofNd2  != 7 || dofNd3  != 7 || dofNd4  != 7 ||
          dofNd5  != 7 || dofNd6  != 7 || dofNd7  != 7 || dofNd8  != 7 ||
          dofNd9  != 7 || dofNd10 != 7 || dofNd11 != 7 || dofNd12 != 7 ||
          dofNd13 != 7 || dofNd14 != 7 || dofNd15 != 7 || dofNd16 != 7 ||
          dofNd17 != 7 || dofNd18 != 7 || dofNd19 != 7 || dofNd20 != 7 ) {
        opserr << "FATAL ERROR TwentyNodeBrick_u_p_U (tag: " << this->getTag() << "), has differing number of DOFs at its nodes\n";
  exit(-1);
      }
      this->DomainComponent::setDomain(theDomain);
    }
}
//=============================================================================
/*const Matrix &TwentyNodeBrick_u_p_U::getTangentStiff ()
{
     opserr<<"stiffness matrix\n\n";
     for (int i=0;i<60;i++)
       for ( int j=0;j<60;j++)
        opserr<<K(i,j);
     opserr<<endln;

     return K;
}    wxy commented 01/15/2002  */
//=============================================================================
const Matrix &TwentyNodeBrick_u_p_U::getInitialStiff ()
{
  if (Ki == 0)
    Ki = new Matrix(this->getTangentStiff());

  if (Ki == 0) {
    opserr << "FATAL TwentyNodeBrick_u_p_U::getInitialStiff() -";
    opserr << "ran out of memory\n";
    exit(-1);
  }

  return *Ki;
}
//=============================================================================
/* const Matrix &TwentyNodeBrick_u_p_U::getDamp ()
{
     return C;
}      */
//=============================================================================
/*const Matrix &TwentyNodeBrick_u_p_U::getMass ()
{
      return M;
}    */
//=============================================================================
void TwentyNodeBrick_u_p_U::zeroLoad()
{
     Q.Zero();
}

//=============================================================================
int
TwentyNodeBrick_u_p_U::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "TwentyNodeBrick_u_p_U::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;

  return -1;
}

//=============================================================================
int TwentyNodeBrick_u_p_U::addInertiaLoadToUnbalance(const Vector &accel)
{

  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = nd1Ptr->getRV(accel);
  const Vector &Raccel2 = nd2Ptr->getRV(accel);
  const Vector &Raccel3 = nd3Ptr->getRV(accel);
  const Vector &Raccel4 = nd4Ptr->getRV(accel);
  const Vector &Raccel5 = nd5Ptr->getRV(accel);
  const Vector &Raccel6 = nd6Ptr->getRV(accel);
  const Vector &Raccel7 = nd7Ptr->getRV(accel);
  const Vector &Raccel8 = nd8Ptr->getRV(accel);
  const Vector &Raccel9 = nd9Ptr->getRV(accel);
  const Vector &Raccel10 = nd10Ptr->getRV(accel);
  const Vector &Raccel11 = nd11Ptr->getRV(accel);
  const Vector &Raccel12 = nd12Ptr->getRV(accel);
  const Vector &Raccel13 = nd13Ptr->getRV(accel);
  const Vector &Raccel14 = nd14Ptr->getRV(accel);
  const Vector &Raccel15 = nd15Ptr->getRV(accel);
  const Vector &Raccel16 = nd16Ptr->getRV(accel);
  const Vector &Raccel17 = nd17Ptr->getRV(accel);
  const Vector &Raccel18 = nd18Ptr->getRV(accel);
  const Vector &Raccel19 = nd19Ptr->getRV(accel);
  const Vector &Raccel20 = nd20Ptr->getRV(accel);
// change the size of Raccel??.Size() to 7   Xiaoyan 02/08/2002
    if (7 != Raccel1.Size()  || 7 != Raccel2.Size()  || 7 != Raccel7.Size()  || 7 != Raccel4.Size() ||
        7 != Raccel5.Size()  || 7 != Raccel6.Size()  || 7 != Raccel7.Size()  || 7 != Raccel8.Size() ||
        7 != Raccel9.Size()  || 7 != Raccel10.Size() || 7 != Raccel11.Size() || 7 != Raccel12.Size()||
        7 != Raccel17.Size() || 7 != Raccel14.Size() || 7 != Raccel15.Size() || 7 != Raccel16.Size()||
        7 != Raccel17.Size() || 7 != Raccel18.Size() || 7 != Raccel19.Size() || 7 != Raccel20.Size()   ){
  // Xiaoyan changed 2 to 3 and added Racce15-18  09/27/00
      opserr << "TwentyNodeBrick_u_p_U::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
      return -1;
    }

  static Vector ra(140);  // change from ra(60) to ra(140)  Xiaoyan 02/01/02

  ra( 0) = Raccel1(0);
  ra( 1) = Raccel1(1);
  ra( 2) = Raccel1(2);
  ra( 3) = 0.0;
  ra( 4) = Raccel1(4);
  ra( 5) = Raccel1(5);
  ra( 6) = Raccel1(6);

  ra( 7) = Raccel2(0);
  ra( 8) = Raccel2(1);
  ra( 9) = Raccel2(2);
  ra( 10) = 0.0;
  ra( 11) = Raccel2(4);
  ra( 12) = Raccel2(5);
  ra( 13) = Raccel2(6);

  ra(14) = Raccel3(0);
  ra(15) = Raccel3(1);
  ra(16) = Raccel3(2);
      ra(17) = 0.0;
  ra(18) = Raccel3(4);
  ra(19) = Raccel3(5);
  ra(20) = Raccel3(6);

  ra(21) = Raccel4(0);
  ra(22) = Raccel4(1);
  ra(23) = Raccel4(2);
  ra(24) = 0.0;
   ra(25) = Raccel4(4);
  ra(26) = Raccel4(5);
  ra(27) = Raccel4(6);

  ra(28) = Raccel5(0);
  ra(29) = Raccel5(1);
  ra(30) = Raccel5(2);
  ra(31) = 0.0;
  ra(32) = Raccel5(4);
  ra(33) = Raccel5(5);
  ra(34) = Raccel5(6);

  ra(35) = Raccel6(0);
  ra(36) = Raccel6(1);
  ra(37) = Raccel6(2);
  ra(38) = 0.0;
  ra(39) = Raccel6(4);
  ra(40) = Raccel6(5);
  ra(41) = Raccel6(6);

  ra(42) = Raccel7(0);
  ra(43) = Raccel7(1);
  ra(44) = Raccel7(2);
  ra(45) = 0.0;
  ra(46) = Raccel7(4);
  ra(47) = Raccel7(5);
  ra(48) = Raccel7(6);

  ra(49) = Raccel8(0);
  ra(50) = Raccel8(1);
  ra(51) = Raccel8(2);
  ra(52) = 0.0;
  ra(53) = Raccel8(4);
  ra(54) = Raccel8(5);
  ra(55) = Raccel8(6);

  ra(56) = Raccel9(0);
  ra(57) = Raccel9(1);
  ra(58) = Raccel9(2);
  ra(59) = 0.0;
  ra(60) = Raccel9(4);
  ra(61) = Raccel9(5);
  ra(62) = Raccel9(6);

   ra(63) = Raccel10(0);
  ra(64) = Raccel10(1);
  ra(65) = Raccel10(2);
  ra(66) = 0.0;
  ra(67) = Raccel10(4);
  ra(68) = Raccel10(5);
  ra(69) = Raccel10(6);

   ra(70) = Raccel11(0);
  ra(71) = Raccel11(1);
  ra(72) = Raccel11(2);
  ra(73) = 0.0;
  ra(74) = Raccel11(4);
  ra(75) = Raccel11(5);
  ra(76) = Raccel11(6);

   ra(77) = Raccel12(0);
  ra(78) = Raccel12(1);
  ra(79) = Raccel12(2);
  ra(80) = 0.0;
  ra(81) = Raccel12(4);
  ra(82) = Raccel12(5);
  ra(83) = Raccel12(6);

  ra( 84) = Raccel13(0);
  ra( 85) = Raccel13(1);
  ra( 86) = Raccel13(2);
  ra( 87) = 0.0;
  ra( 88) = Raccel13(4);
  ra( 89) = Raccel13(5);
  ra( 90) = Raccel13(6);

  ra(91) = Raccel14(0);
  ra(92) = Raccel14(1);
      ra(92) = Raccel14(2);
  ra(94) = 0.0;
  ra(95) = Raccel14(4);
  ra(96) = Raccel14(5);
  ra(97) = Raccel14(6);

  ra(98) = Raccel15(0);
  ra(99) = Raccel15(1);
  ra(100) = Raccel15(2);
  ra(101) = 0.0;
  ra(102) = Raccel15(4);
  ra(103) = Raccel15(5);
  ra(104) = Raccel15(6);

  ra(105) = Raccel16(0);
   ra(106) = Raccel16(1);
  ra(107) = Raccel16(2);
  ra(108) = 0.0;
  ra(109) = Raccel16(4);
  ra(110) = Raccel16(5);
  ra(111) = Raccel16(6);

  ra(112) = Raccel17(0);
  ra(113) = Raccel17(1);
  ra(114) = Raccel17(2);
  ra(115) = 0.0;
  ra(116) = Raccel17(4);
  ra(117) = Raccel17(5);
  ra(118) = Raccel17(6);

  ra(119) = Raccel18(0);
  ra(120) = Raccel18(1);
  ra(121) = Raccel18(2);
  ra(122) = 0.0;
  ra(123) = Raccel18(4);
  ra(124) = Raccel18(5);
  ra(125) = Raccel18(6);

  ra(126) = Raccel19(0);
  ra(127) = Raccel19(1);
  ra(128) = Raccel19(2);
  ra(129) = 0.0;
  ra(130) = Raccel19(4);
  ra(131) = Raccel19(5);
  ra(132) = Raccel19(6);

  ra(133) = Raccel20(0);
  ra(134) = Raccel20(1);
  ra(135) = Raccel20(2);
  ra(136) = 0.0;
  ra(137) = Raccel20(4);
  ra(138) = Raccel20(5);
  ra(139) = Raccel20(6);
//out 06/03/2002      int i;
//out 06/03/2002      for ( i=0; i<140; i++)
//out 06/03/2002       {
//out 06/03/2002         if(i%10==0) opserr<<endln;
//out 06/03/2002          opserr<<" "<<ra(i);
//out 06/03/2002       }     // added to check the value of ra      02/18/2002

    // Want to add ( - fact * M R * accel ) to unbalance
    // Take advantage of lumped mass matrix
    // Mass matrix is computed in setDomain()

    //double column_mass = 0;
    //for (int i = 0; i < 24; i++)
    //   column_mass += M(1,i);
    //column_mass = column_mass/3.0;

    //opserr << " addInerti... column_mass " << column_mass << endln;

//    for ( i = 0; i < nodes_in_brick*7; i++)
//       {
//          Q(i)    += -M(i,i)*ra(i);
//    opserr<<Q(i);  // 05/02/2002
//       }
    Q.addMatrixVector(1.0, M, ra, -1.0);
    return 0;
}

//=============================================================================
const Vector TwentyNodeBrick_u_p_U::FormEquiBodyForce()
{
    Vector bforce(140);

    // Check for a quick return
    //opserr << "rho " << rho << endln;
    if (rho == 0.0)
      return bforce;

    Vector ba(140);
// The body force(gravity g) are same for solid part and fluid part
//     Solid part         Fluid part
    ba( 0) = bf(0);
    ba( 1) = bf(1);
    ba( 2) = bf(2);
    ba( 3)=  0.0;
    ba( 4) = bf(0);
    ba( 5) = bf(1);
    ba( 6) = bf(2);

    ba( 7) = bf(0);
    ba( 8) = bf(1);
    ba( 9) = bf(2);
    ba(10) = 0.0;
    ba(11) = bf(0);
    ba(12) = bf(1);
    ba(13) = bf(2);

    ba(14) = bf(0);
    ba(15) = bf(1);
    ba(16) = bf(2);
    ba(17) = 0.0;
    ba(18) = bf(0);
    ba(19) = bf(1);
    ba(20) = bf(2);

    ba(21) = bf(0);
    ba(22) = bf(1);
    ba(23) = bf(2);
    ba(24) = 0.0;
    ba(25) = bf(0);
    ba(26) = bf(1);
    ba(27) = bf(2);

    ba(28) = bf(0);
    ba(29) = bf(1);
    ba(30) = bf(2);
    ba(31) = 0.0;
    ba(32) = bf(0);
    ba(33) = bf(1);
    ba(34) = bf(2);

    ba(35) = bf(0);
    ba(36) = bf(1);
    ba(37) = bf(2);
    ba(38) = 0.0;
    ba(39) = bf(0);
    ba(40) = bf(1);
    ba(41) = bf(2);

    ba(42) = bf(0);
    ba(43) = bf(1);
    ba(44) = bf(2);
    ba(45) = 0.0;
    ba(46) = bf(0);
    ba(47) = bf(1);
    ba(48) = bf(2);

    ba(49) = bf(0);
    ba(50) = bf(1);
    ba(51) = bf(2);
    ba(52) = 0.0;
    ba(53) = bf(0);
    ba(54) = bf(1);
    ba(55) = bf(2);

    ba(56) = bf(0);
    ba(57) = bf(1);
    ba(58) = bf(2);
    ba(59) = 0.0;
    ba(60) = bf(0);
    ba(61) = bf(1);
    ba(62) = bf(2);

    ba(63) = bf(0);
    ba(64) = bf(1);
    ba(65) = bf(2);
    ba(66) = 0.0;
    ba(67) = bf(0);
    ba(68) = bf(1);
    ba(69) = bf(2);

    ba(70) = bf(0);
    ba(71) = bf(1);
    ba(72) = bf(2);
    ba(73) = 0.0;
    ba(74) = bf(0);
    ba(75) = bf(1);
    ba(76) = bf(2);

    ba(77) = bf(0);
    ba(78) = bf(1);
    ba(79) = bf(2);
    ba(80) = 0.0;
    ba(81)  =bf(0);
    ba(82)  =bf(1);
    ba(83)  =bf(2);

    ba(84)  =bf(0);
    ba(85)  =bf(1);
    ba(86)  =bf(2);
    ba(87)  =0.0;
    ba(88)  =bf(0);
    ba(89)  =bf(1);
    ba(90)  =bf(2);

    ba(91)  =bf(0);
    ba(92)  =bf(1);
    ba(93)  =bf(2);
    ba(94)  =0.0;
    ba(95)  =bf(0);
    ba(96)  =bf(1);
    ba(97)  =bf(2);

    ba(98)  =  bf(0);
    ba(99)  =  bf(1);
    ba(100) =  bf(2);
    ba(101)  = 0.0;
    ba(102)  = bf(0);
    ba(103)  = bf(1);
    ba(104)  = bf(2);

    ba(105)  = bf(0);
    ba(106)  = bf(1);
    ba(107)  = bf(2);
    ba(108)  = 0.0;
    ba(109)  = bf(0);
    ba(110)  = bf(1);
    ba(111)  = bf(2);

    ba(112)  = bf(0);
    ba(113)  = bf(1);
    ba(114)  = bf(2);
    ba(115)  = 0.0;
    ba(116)  = bf(0);
    ba(117)  = bf(1);
    ba(118)  = bf(2);

    ba(119)  = bf(0);
    ba(120)  = bf(1);
    ba(121)  = bf(2);
    ba(122)  = 0.0;
    ba(123)  = bf(0);
    ba(124)  = bf(1);
    ba(125)  = bf(2);

    ba(126)  = bf(0);
    ba(127)  = bf(1);
    ba(128)  = bf(2);
    ba(129)  = 0.0;
    ba(130)  = bf(0);
    ba(131)  = bf(1);
    ba(132)  = bf(2);

    ba(133)  = bf(0);
    ba(134)  = bf(1);
    ba(135)  = bf(2);
    ba(136)  = 0.0;
    ba(137)  = bf(0);
    ba(138)  = bf(1);
    ba(139)  = bf(2);


    //Form equivalent body force
    bforce.addMatrixVector(0.0, M, ba, 1.0);
    //opserr << " ba " << ba;
    //opserr << " M " << M;
    //if (getTag() == 886)
//    opserr << " @@@@@ FormEquiBodyForce  " << bforce;

    return bforce;
}

//=============================================================================
// Setting initial E according to the initial pressure p
//void TwentyNodeBrick_u_p_U::setInitE(void)
//{
//    //Get the coors of each node
//
//    const Vector &nd1Crds = nd1Ptr->getCrds();
//    const Vector &nd2Crds = nd2Ptr->getCrds();
//    const Vector &nd3Crds = nd3Ptr->getCrds();
//    const Vector &nd4Crds = nd4Ptr->getCrds();
//    const Vector &nd5Crds = nd5Ptr->getCrds();
//    const Vector &nd6Crds = nd6Ptr->getCrds();
//    const Vector &nd7Crds = nd7Ptr->getCrds();
//    const Vector &nd8Crds = nd8Ptr->getCrds();
//
//    //dir is the ID for vertial direction, e.g. 1 means x-dir is vertical...
//    double Zavg = nd1Crds( dir-1)+
//           nd2Crds( dir-1)+
//           nd3Crds( dir-1)+
//           nd4Crds( dir-1)+
//           nd5Crds( dir-1)+
//           nd6Crds( dir-1)+
//           nd7Crds( dir-1)+
//           nd8Crds( dir-1);
//    Zavg = Zavg / 8;
//
//    //Estimate the pressure at that depth
//    double sigma_v = (Zavg - surflevel) * rho * 9.81; //units in SI system
//    double ko = 0.5;
//    double p_est = sigma_v*( 2.0*ko+1.0)/3.0;
//    //opserr << " Initial P " << p_est << endln;
//
//    int i;
//
//    // Loop over the integration points and set the initial material state
//    int count  = r_integration_order* s_integration_order * t_integration_order;
//
//    //For elastic-isotropic material
//    if (strcmp(matpoint[i]->matmodel->getType(),"ElasticIsotropic3D") == 0)
//    {
//       for (i = 0; i < count; i++)
//           (matpoint[i]->matmodel)->setElasticStiffness( p_est );
//    }
//
//    //return ;
//}


//=============================================================================
const Vector &TwentyNodeBrick_u_p_U::getResistingForce ()
{
    int force_dim[] = {20,7};
    tensor nodalforces(2,force_dim,0.0);   // The index of tensor starts from 1
                  // See SRC/nDarray/BJtensor.cpp
    nodalforces = nodal_forces();          // The index of Vector starts from 0
                   // SRC/matrix/Vector.h  Xiaoyan 02/10/20002
    //converting nodalforce tensor to vector
    for (int i = 0; i<nodes_in_brick; i++)
      for (int j = 0; j < 7; j++)
  p(i *7 + j) = nodalforces.cval(i+1, j+1);    // p--Vector--start 0
                 // nodalforces--tensor start 1
   // opserr << "p" << p;
   // opserr << "Q" << Q;
/*   for ( int i=0; i<140; i++)
    {
      printf("p(%d) = %lf    ",i, p(i));
      printf("Q(%d) = %lf    ",i, Q(i));
      p(i)=p(i)-Q(i);
//      printf("           --->>>  p(%d) = %lf    ",i, p(i));
      printf("p(%d)-Q(%d) = %lf   \n ",i, i, p(i));

    }
*/
    p = p - Q;
//    for ( int i=0; i<140; i++)
//      printf("p(%d) = %lf   ",i, p(i));   // check
//    //opserr << "p-Q" << p;
//    fflush(stdout);
    return p;
}

//=============================================================================
const Vector &TwentyNodeBrick_u_p_U::getResistingForceIncInertia ()
{
  // Check for a quick return
  if (rho == 0.0)
    return this->getResistingForce();

  //opserr << "Node555 trialDisp " << nd1Ptr->getTrialDisp();

  const Vector &accel1 = nd1Ptr->getTrialAccel();
        //opserr << "\nnode accel " << nd1Ptr->getTag() << " ux " << accel1(0) <<" uy "<< accel1(1) << " uz "<< accel1(2) << endln;

  const Vector &accel2 = nd2Ptr->getTrialAccel();
        //opserr << "node accel " << nd2Ptr->getTag() << " ux " << accel2(0) <<" uy "<< accel2(1) << " uz "<< accel2(2) << endln;

  const Vector &accel3 = nd3Ptr->getTrialAccel();
        //opserr << "node accel " << nd3Ptr->getTag() << " ux " << accel3(0) <<" uy "<< accel3(1) << " uz "<< accel3(2) << endln;

  const Vector &accel4 = nd4Ptr->getTrialAccel();
        //opserr << "node accel " << nd4Ptr->getTag() << " ux " << accel4(0) <<" uy "<< accel4(1) << " uz "<< accel4(2) << endln;

        // Xiaoyan added the following four 09/27/00
  const Vector &accel5 = nd5Ptr->getTrialAccel();
        //opserr << "node accel " << nd5Ptr->getTag() << " ux " << accel5(0) <<" uy "<< accel5(1) << " uz "<< accel5(2) << endln;

  const Vector &accel6 = nd6Ptr->getTrialAccel();
        //opserr << "node accel " << nd6Ptr->getTag() << " ux " << accel6(0) <<" uy "<< accel6(1) << " uz "<< accel6(2) << endln;

  const Vector &accel7 = nd7Ptr->getTrialAccel();
        //opserr << "node accel " << nd7Ptr->getTag() << " ux " << accel7(0) <<" uy "<< accel7(1) << " uz "<< accel7(2) << endln;

  const Vector &accel8 = nd8Ptr->getTrialAccel();
        //opserr << "node accel " << nd8Ptr->getTag() << " ux " << accel8(0) <<" uy "<< accel8(1) << " uz "<< accel8(2) << endln;

  const Vector &accel9 = nd9Ptr->getTrialAccel();
        //opserr << "node accel " << nd9Ptr->getTag() << " ux " << accel9(0) <<" uy "<< accel9(1) << " uz "<< accel9(2) << endln;

  const Vector &accel10 = nd10Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel11 = nd11Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel12 = nd12Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel13 = nd13Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel14 = nd14Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel15 = nd15Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel16 = nd16Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel17 = nd17Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel18 = nd18Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel19 = nd19Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;

  const Vector &accel20 = nd20Ptr->getTrialAccel();
        //opserr << "node accel " << nd10Ptr->getTag() << " ux " << accel10(0) <<" uy "<< accel10(1) << " uz "<< accel10(2) << endln;


  static Vector a(140);  // originally 8

// Modified to 3-solid 1-p 3-fluid        Xiaoyan 03/14/2002
  a( 0) = accel1(0);
  a( 1) = accel1(1);
  a( 2) = accel1(2);
  a( 3) = 0.0;
  a( 4) = accel1(4);
  a( 5) = accel1(5);
  a( 6) = accel1(6);

  a( 7) = accel2(0);
  a( 8) = accel2(1);
  a( 9) = accel2(2);
  a(10) = 0.0;
  a(11) = accel2(4);
      a(12) = accel2(5);
  a(13) = accel2(6);

  a(14) = accel3(0);
  a(15) = accel3(1);
  a(16) = accel3(2);
  a(17) = 0.0;
  a(18) = accel3(4);
  a(19) = accel3(5);
  a(20) = accel3(6);

  a(21) = accel4(0);
  a(22) = accel4(1);
  a(23) = accel4(2);
  a(24) = 0.0;
   a(25) = accel4(4);
  a(26) = accel4(5);
  a(27) = accel4(6);

  a(28) = accel5(0);
  a(29) = accel5(1);
  a(30) = accel5(2);
  a(31) = 0.0;
  a(32) = accel5(4);
  a(33) = accel5(5);
  a(34) = accel5(6);

  a(35) = accel6(0);
  a(36) = accel6(1);
  a(37) = accel6(2);
  a(38) = 0.0;
  a(39) = accel6(4);
  a(40) = accel6(5);
  a(41) = accel6(6);

  a(42) = accel7(0);
  a(43) = accel7(1);
  a(44) = accel7(2);
  a(45) = 0.0;
  a(46) = accel7(4);
  a(47) = accel7(5);
  a(48) = accel7(6);

  a(49) = accel8(0);
  a(50) = accel8(1);
  a(51) = accel8(2);
  a(52) = 0.0;
  a(53) = accel8(4);
  a(54) = accel8(5);
  a(55) = accel8(6);

  a(56) = accel9(0);
  a(57) = accel9(1);
  a(58) = accel9(2);
  a(59) = 0.0;
  a(60) = accel9(4);
  a(61) = accel9(5);
  a(62) = accel9(6);

  a(63) = accel10(0);
  a(64) = accel10(1);
  a(65) = accel10(2);
  a(66) = 0.0;
  a(67) = accel10(4);
  a(68) = accel10(5);
  a(69) = accel10(6);

   a(70) = accel11(0);
  a(71) = accel11(1);
  a(72) = accel11(2);
  a(73) = 0.0;
  a(74) = accel11(4);
  a(75) = accel11(5);
  a(76) = accel11(6);

   a(77) = accel12(0);
  a(78) = accel12(1);
  a(79) = accel12(2);
  a(80) = 0.0;
  a(81) = accel12(4);
  a(82) = accel12(5);
  a(83) = accel12(6);

  a( 84) = accel13(0);
  a( 85) = accel13(1);
  a( 86) = accel13(2);
  a( 87) = 0.0;
  a( 88) = accel13(4);
  a( 89) = accel13(5);
  a( 90) = accel13(6);

   a(91) = accel14(1);
  a(92) = accel14(1);
      a(93) = accel14(2);
  a(94) = 0.0;
  a(95) = accel14(4);
  a(96) = accel14(5);
  a(97) = accel14(6);

  a(98)  = accel15(0);
  a(99)  = accel15(1);
  a(100)  = accel15(2);
  a(101) = 0.0;
  a(102) = accel15(4);
  a(103) = accel15(5);
  a(104) = accel15(6);

  a(105) = accel16(0);
   a(106) = accel16(1);
  a(107) = accel16(2);
  a(108) = 0.0;;
  a(109) = accel16(4);
  a(110) = accel16(5);
  a(111) = accel16(6);

  a(112) = accel17(0);
  a(113) = accel17(1);
  a(114) = accel17(2);
  a(115) = 0.0;
  a(116) = accel17(4);
  a(117) = accel17(5);
  a(118) = accel17(6);

  a(119) = accel18(0);
  a(120) = accel18(1);
  a(121) = accel18(2);
  a(122) = 0.0;
  a(123) = accel18(4);
  a(124) = accel18(5);
  a(125) = accel18(6);

  a(126) = accel19(0);
  a(127) = accel19(1);
  a(128) = accel19(2);
  a(129) = 0.0;
  a(130) = accel19(4);
  a(131) = accel19(5);
  a(132) = accel19(6);

  a(133) = accel20(0);
  a(134) = accel20(1);
  a(135) = accel20(2);
  a(136) = 0.0;
  a(137) = accel20(4);
  a(138) = accel20(5);
  a(139) = accel20(6);



  // Compute the current resisting force
  this->getResistingForce();

  // Take advantage of lumped mass matrix
  // Mass matrix is computed in setDomain()
  //opserr << " M_ii \n";

        //double column_mass = 0;
        //for (int i = 0; i < 24; i++)
        //   column_mass += M(1,i);
        //column_mass = column_mass/3.0;

/*  for (int i = 0; i < 140; i++)
  {
     p(i ) += M(i,i)*a(i);
     //opserr << " " << M(i, i);
  }
  //opserr << endln;
  //opserr << "P+=Ma" << P<< endln;
*/
  p.addMatrixVector(1.0, M, a, 1.0);
  return p;
}


////#############################################################################
int TwentyNodeBrick_u_p_U::get_global_number_of_node(int local_node_number)
{
  //return G_N_numbs[local_node_number];
  return connectedExternalNodes(local_node_number);
}

////#############################################################################
int  TwentyNodeBrick_u_p_U::get_Brick_Number()
{
  //return elem_numb;
  return this->getTag();
}

////#############################################################################
//=============================================================================
double TwentyNodeBrick_u_p_U::get_Gauss_p_c(short order, short point_numb)
  {
//  Abscissae coefficient of the Gaussian quadrature formula
// starting from 1 not from 0
    static double Gauss_coordinates[7][7];

    Gauss_coordinates[1][1] = 0.0 ;
    Gauss_coordinates[2][1] = -0.577350269189626;
    Gauss_coordinates[2][2] = -Gauss_coordinates[2][1];
    Gauss_coordinates[3][1] = -0.774596669241483;
    Gauss_coordinates[3][2] = 0.0;
    Gauss_coordinates[3][3] = -Gauss_coordinates[3][1];
    Gauss_coordinates[4][1] = -0.861136311594053;
    Gauss_coordinates[4][2] = -0.339981043584856;
    Gauss_coordinates[4][3] = -Gauss_coordinates[4][2];
    Gauss_coordinates[4][4] = -Gauss_coordinates[4][1];
    Gauss_coordinates[5][1] = -0.906179845938664;
    Gauss_coordinates[5][2] = -0.538469310105683;
    Gauss_coordinates[5][3] = 0.0;
    Gauss_coordinates[5][4] = -Gauss_coordinates[5][2];
    Gauss_coordinates[5][5] = -Gauss_coordinates[5][1];
    Gauss_coordinates[6][1] = -0.932469514203152;
    Gauss_coordinates[6][2] = -0.661209386466265;
    Gauss_coordinates[6][3] = -0.238619186083197;
    Gauss_coordinates[6][4] = -Gauss_coordinates[6][3];
    Gauss_coordinates[6][5] = -Gauss_coordinates[6][2];
    Gauss_coordinates[6][6] = -Gauss_coordinates[6][1];

    return Gauss_coordinates[order][point_numb];
 }
////#############################################################################

double TwentyNodeBrick_u_p_U::get_Gauss_p_w(short order, short point_numb)
  {
//  Weight coefficient of the Gaussian quadrature formula
// starting from 1 not from 0
    static double Gauss_weights[7][7]; // static data ??

    Gauss_weights[1][1] = 2.0;
    Gauss_weights[2][1] = 1.0;
    Gauss_weights[2][2] = 1.0;
    Gauss_weights[3][1] = 0.555555555555556;
    Gauss_weights[3][2] = 0.888888888888889;
    Gauss_weights[3][3] = Gauss_weights[3][1];
    Gauss_weights[4][1] = 0.347854845137454;
    Gauss_weights[4][2] = 0.652145154862546;
    Gauss_weights[4][3] = Gauss_weights[4][2];
    Gauss_weights[4][4] = Gauss_weights[4][1];
    Gauss_weights[5][1] = 0.236926885056189;
    Gauss_weights[5][2] = 0.478628670499366;
    Gauss_weights[5][3] = 0.568888888888889;
    Gauss_weights[5][4] = Gauss_weights[5][2];
    Gauss_weights[5][5] = Gauss_weights[5][1];
    Gauss_weights[6][1] = 0.171324492379170;
    Gauss_weights[6][2] = 0.360761573048139;
    Gauss_weights[6][3] = 0.467913934572691;
    Gauss_weights[6][4] = Gauss_weights[6][3];
    Gauss_weights[6][5] = Gauss_weights[6][2];
    Gauss_weights[6][6] = Gauss_weights[6][1];

    return Gauss_weights[order][point_numb];
  }


////#############################################################################
int * TwentyNodeBrick_u_p_U::get_LM()
  {
    return LM;
  }

/////#############################################################################
//void TwentyNodeBrick_u_p_U::set_LM(Node * node)
//  {
////    unsigned int BrickNumber = this->get_Brick_Number();
////    this->reportshort("");
//// for element numbered BrickNumber create LM array (see Bathe pp984
////    for (int LocalNodeNumber = 1 ; LocalNodeNumber<=20 ; LocalNodeNumber++ )
//    for (int LocalNodeNumber = 1 ; LocalNodeNumber<=8 ; LocalNodeNumber++ )// for 8noded brick
//      {
////        unsigned int global_node_number = b3d[BrickNumber-1].get_global_number_of_node(LocalNodeNumber-1);
//        unsigned int global_node_number = this->get_global_number_of_node(LocalNodeNumber-1);
//        LM[3*LocalNodeNumber-3] = node[global_node_number].eqn_tx();
//        LM[3*LocalNodeNumber-2] = node[global_node_number].eqn_ty();
//        LM[3*LocalNodeNumber-1] = node[global_node_number].eqn_tz();
//      }
//
//      // ::printf("\n\n");
//
////===   this->reportLM("LM");
////   for (int count01=1;count01<=8;count01++)
////     {
////       ::printf("element %4d localNode %4d Globalnode %4d  LM   %4d   %4d   %4d\n", BrickNumber, count01,this->get_global_number_of_node(count01-1),  LM[count01*3-3], LM[count01*3-2], LM[count01*3-1] );
////     }
//
//  }


////#############################################################################
// returns nodal forces of solid part in an element
tensor TwentyNodeBrick_u_p_U::nodal_forcesFu()
  {
    int force_dim[] = {20,3};

    tensor nodal_forcesFu(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dh(2, dh_dim, 0.0);

    stresstensor stress_at_GP(0.0);

    double det_of_Jacobian = 0.0;

    straintensor incremental_strain;

    static int disp_dim[] = {20,3};
    tensor incremental_displacements(2,disp_dim,0.0); // \Delta u

    incremental_displacements = incr_dispDu();

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );

                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;

                // derivatives of local coordiantes with respect to local coordiantes
                dh = dh_drst_at(r,s,t);

                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //....                Jacobian.print("J");

                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");

                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );

                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

    incremental_strain =
                     (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
//    if (where == 0)
//       //opserr << " In nodal_force delta_incremental_strain tag "<< getTag() <<"  " <<incremental_strain << endln;
////    opserr << " el tag = "<< getTag();
//
    int err = ( matpoint[where]->matmodel )->setTrialStrainIncr( incremental_strain);
    if ( err)
      opserr << "TwentyNodeBrick_u_p_U::getStiffnessTensor (tag: " << this->getTag() << "), not converged\n";



    //char *test = matpoint[where]->matmodel->getType();
    // fmk - changing if so if into else block must be Template3Dep
//    if (strcmp(matpoint[where]->matmodel->getType(),"Template3Dep") != 0)
       stress_at_GP = matpoint[where]->getStressTensor();


                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forcesFu =
                     nodal_forcesFu + dhGlobal("ib")*stress_at_GP("ab")*weight;
                //nodal_forces.print("nf","\n\n Nodal Forces \n");

              }
          }
      }

//    opserr << "\n element no. " << getTag() << endln;
    nodal_forcesFu.print("Fu","\n Nodal Forces \n");
    return nodal_forcesFu;

  }
////#############################################################################
// returns nodal forces of fluid part in an element
tensor TwentyNodeBrick_u_p_U::nodal_forcesFU()
  {
 // This is a function for a generelized force for the fluid component...
 // For fluid there is only pressure not shear force, so the force for
 // fluid will be p*delta_{ij}      xiaoyan 05/20/2002
    int force_dim[] = {20,3};

    tensor nodal_forcesFU(2,force_dim,0.0);
    tensor I2("I", 2, def_dim_2);     // delta_{ij} added 05/20/2002

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};
    int h_dim[] = {20};     // added for shape function h  05/21/2002

    tensor dh(2, dh_dim, 0.0);
    tensor h(1, h_dim, 0.0);

    stresstensor stress_at_GP(0.0);
//    double  p_at_GP = 0;
    tensor p_at_GP;      // added 05/21/2002
    double pp;        // added 05/21/2002

    double det_of_Jacobian = 0.0;

//    straintensor incremental_p;

    static int disp_dim[] = {20};     // changed from {20,1}  05/20/2002
    tensor incremental_p(1,disp_dim,0.0); // \Delta U

    incremental_p = incr_dispDp();

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );

                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;

                // derivatives of local coordiantes with respect to local coordiantes
                dh = dh_drst_at(r,s,t);
    h=interp_poli_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //....                Jacobian.print("J");

                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");

                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );

                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

//    incremental_strain =
//                     (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
//    if (where == 0)
//       //opserr << " In nodal_force delta_incremental_strain tag "<< getTag() <<"  " <<incremental_strain << endln;
////    opserr << " el tag = "<< getTag();
//
//    int err = ( matpoint[where]->matmodel )->setTrialStrainIncr( incremental_strain);
//    if ( err)
//                    opserr << "TwentyNodeBrick_u_p_U::getStiffnessTensor (tag: %d), not converged",
//             this->getTag());



    //char *test = matpoint[where]->matmodel->getType();
    // fmk - changing if so if into else block must be Template3Dep
//    if (strcmp(matpoint[where]->matmodel->getType(),"Template3Dep") != 0)
//       stress_at_GP = matpoint[where]->getStressTensor();
       p_at_GP=h("i")*incremental_p("i");
       pp=p_at_GP.val(1);
/*       p_at_GP = h.val(1)*incremental_p.val(1)+h.val(2)*incremental_p.val(2)
                +h.val(3)*incremental_p.val(3)+h.val(4)*incremental_p.val(4)
          +h.val(5)*incremental_p.val(5)+h.val(6)*incremental_p.val(6)
          +h.val(7)*incremental_p.val(7)+h.val(8)*incremental_p.val(8)
          +h.val(9)*incremental_p.val(9)+h.val(10)*incremental_p.val(10)
          +h.val(11)*incremental_p.val(11)+h.val(12)*incremental_p.val(12)
          +h.val(13)*incremental_p.val(13)+h.val(14)*incremental_p.val(14)
          +h.val(15)*incremental_p.val(15)+h.val(16)*incremental_p.val(16)
          +h.val(17)*incremental_p.val(17)+h.val(18)*incremental_p.val(18)
          +h.val(19)*incremental_p.val(19)+h.val(20)*incremental_p.val(20);
 */
//          stress_at_GP = p_at_GP*I2("ab");

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forcesFU =
                     nodal_forcesFU + dhGlobal("ib")*I2("ab")*pp*weight;
                //nodal_forces.print("nf","\n\n Nodal Forces \n");

              }
          }
      }

//    opserr << "\n element no. " << getTag() << endln;
    nodal_forcesFU.print("FU","\n Nodal Forces \n");

    return nodal_forcesFU;

  }
////#############################################################################
// returns nodal forces for given stress field in an element
tensor TwentyNodeBrick_u_p_U::nodal_forces()
  {
    int force_dim[] = {20,7};

    tensor nodalForceFu = nodal_forcesFu();
    tensor nodalForceFU = nodal_forcesFU();

    tensor nodal_FORCES(2,force_dim,0.0);

    for( int i=1; i<=20; i++)
      {
    // part for u
          nodal_FORCES.val(i,1) = nodalForceFu.val(i,1);
          nodal_FORCES.val(i,2) = nodalForceFu.val(i,2);
          nodal_FORCES.val(i,3) = nodalForceFu.val(i,3);
    // part for p
          nodal_FORCES.val(i,4) = 0.0;
    // part for U
          nodal_FORCES.val(i,5) = nodalForceFU.val(i,1);
          nodal_FORCES.val(i,6) = nodalForceFU.val(i,2);
          nodal_FORCES.val(i,7) = nodalForceFU.val(i,3);
      }
// the following two lines are added  on 05/02/2002
//     opserr << "\n element no. " << getTag() << endln;
     nodal_FORCES.print("nf","\n Nodal Forces cd\n");
     return nodal_FORCES;
  }


////#############################################################################
// returns nodal forces for given ITERATIVE stress field in an element
// Solid part  Xiaoyan 02/08/2002
tensor TwentyNodeBrick_u_p_U::iterative_nodal_forcesFu()
  {
    int force_dim[] = {20,3};

    tensor nodal_forcesFu(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};

    tensor dh(2, dh_dim, 0.0);

    stresstensor stress_at_GP(0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );

                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                //.....
                //.....::printf("TwentyNodeBrick_u_p_U::iterative_nodal_forces()  ----**************** where = %d \n", where);
                //.....::printf("UPDATE ");
                //.....::printf("   GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //.....                           GP_c_r,GP_c_s,GP_c_t);
                // derivatives of local coordiantes with respect to local coordiantes
                dh = dh_drst_at(r,s,t);

                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //....                Jacobian.print("J");

                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");

                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );

                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

                //                   stress_at_GP = (GPstress)->operator[](where);
                //stress_at_GP = GPiterative_stress[where];

    //stress_at_GP = ( matpoint[where].getTrialEPS() )->getStress();
                stress_at_GP = matpoint[where]->getStressTensor();
                stress_at_GP.reportshortpqtheta("\n iterative_stress at GAUSS point in iterative_nodal_force\n");

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forcesFu =
                  nodal_forcesFu + dhGlobal("ib")*stress_at_GP("ab")*weight;
                //nodal_forces.print("nf","\n TwentyNodeBrick_u_p_U::iterative_nodal_forces Nodal Forces ~~~~\n");

              }
          }
      }


    return nodal_forcesFu;

  }
////#############################################################################
// returns nodal forces for given ITERATIVE stress field in an element
// Fluid part  Xiaoyan 02/08/2002
tensor TwentyNodeBrick_u_p_U::iterative_nodal_forcesFU()
  {
    int force_dim[] = {20,3};
    tensor nodal_forcesFU(2,force_dim,0.0);
    tensor I2("I", 2, def_dim_2);     // delta_{ij} added 05/20/2002

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};
    int h_dim[] = {20};     // added for shape function h  05/21/2002

    tensor dh(2, dh_dim, 0.0);
    tensor h(1, h_dim, 0.0);

    stresstensor stress_at_GP(0.0);
    tensor p_at_GP;      // added 05/21/2002
    double pp;        // added 05/21/2002

    double det_of_Jacobian = 0.0;
    static int disp_dim[] = {20};     // changed from {20,1}  05/20/2002
    tensor incremental_p(1,disp_dim,0.0); // \Delta U

    incremental_p = incr_dispDp();

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );

                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                //.....
                //.....::printf("TwentyNodeBrick_u_p_U::iterative_nodal_forces()  ----**************** where = %d \n", where);
                //.....::printf("UPDATE ");
                //.....::printf("   GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //.....                           GP_c_r,GP_c_s,GP_c_t);
                // derivatives of local coordiantes with respect to local coordiantes
                dh = dh_drst_at(r,s,t);
    h=interp_poli_at(r,s,t);

                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //....                Jacobian.print("J");

                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");

                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );

                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

                //                   stress_at_GP = (GPstress)->operator[](where);
                //stress_at_GP = GPiterative_stress[where];

    //stress_at_GP = ( matpoint[where].getTrialEPS() )->getStress();
//                stress_at_GP = matpoint[where]->getStressTensor();
//                stress_at_GP.reportshortpqtheta("\n iterative_stress at GAUSS point in iterative_nodal_force\n");
        p_at_GP=h("i")*incremental_p("i");
       pp=p_at_GP.val(1);

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forcesFU =
                  nodal_forcesFU + dhGlobal("ib")*I2("ab")*pp*weight;
                //nodal_forces.print("nf","\n TwentyNodeBrick_u_p_U::iterative_nodal_forces Nodal Forces ~~~~\n");

              }
          }
      }


    return nodal_forcesFU;

  }


////#############################################################################
// returns nodal forces for given constant stress field in the element
// This function is not called. Xiaoyan 02/08/2002
tensor TwentyNodeBrick_u_p_U::nodal_forces_from_stress(stresstensor & stress)
  {
    int force_dim[] = {20,3};

    tensor nodal_forces(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    double weight = 0.0;

    int dh_dim[] = {20,3}; // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor dh(2, dh_dim, 0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );

                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                //--                where =
                //--                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                //.....
                //.....::printf("TwentyNodeBrick_u_p_U::iterative_nodal_forces()  ----**************** where = %d \n", where);
                //.....::printf("UPDATE ");
                //.....::printf("   GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //.....                           GP_c_r,GP_c_s,GP_c_t);
                // derivatives of local coordiantes with respect to local coordiantes
                dh = dh_drst_at(r,s,t);

                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //....                Jacobian.print("J");

                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");

                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );

                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

                //                   stress_at_GP = (GPstress)->operator[](where);
                //                stress_at_GP = GPiterative_stress[where];
                //GPiterative_stress[where].reportshortpqtheta("\n iterative_stress at GAUSS point in iterative_nodal_force\n");

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forces =
                  nodal_forces + dhGlobal("ib")*stress("ab")*weight;
                //nodal_forces.print("nf","\n TwentyNodeBrick_u_p_U::iterative_nodal_forces Nodal Forces ~~~~\n");

              }
          }
      }

    return nodal_forces;

  }


////#############################################################################
// returns nodal forces for given incremental strain field in an element
// by using the linearized constitutive tensor from the begining of the step !
// This function is not called  Xiaoyan 02/08/2002
tensor TwentyNodeBrick_u_p_U::linearized_nodal_forces()
  {
    int force_dim[] = {20,3};  // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor linearized_nodal_forces(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {20,3};  // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor dh(2, dh_dim, 0.0);

    tensor Constitutive( 4, def_dim_4, 0.0);

    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {20,3};  // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor incremental_displacements(2,disp_dim,0.0);

    straintensor incremental_strain;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    stresstensor final_linearized_stress;
    //    stresstensor incremental_stress;
    // tensor of incremental displacements taken from node objects for this element !
    incremental_displacements = incr_dispDu();
    //incremental_displacements.print("disp","\n incremental_displacements tensor at GAUSS point in iterative_Update\n");

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                tw = get_Gauss_p_w( t_integration_order, GP_c_t );

                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;

                // derivatives of local coordiantes with respect to local coordiantes
                dh = dh_drst_at(r,s,t);

                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                //....                Jacobian.print("J");

                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");

                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );

                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;
                //..::printf("\n\nIN THE nodal forces ----**************** where = %d \n", where);
                //..::printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //..                           GP_c_r,GP_c_s,GP_c_t);
                //..::printf("WEIGHT = %f", weight);
                //..::printf("determinant of Jacobian = %f", det_of_Jacobian);
                // incremental straines at this Gauss point
                // now in Update we know the incremental displacements so let's find
                // the incremental strain
                incremental_strain =
                 (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
                incremental_strain.null_indices();
                //incremental_strain.reportshort("\n incremental_strain tensor at GAUSS point in iterative_Update\n");

                //Constitutive = GPtangent_E[where];

          //EPState *tmp_eps = (matpoint[where]).getEPS();
          //NDMaterial *tmp_ndm = (matpoint[where]).getNDMat();

    //if ( tmp_eps ) {     //Elasto-plastic case
    //    mmodel->setEPS( *tmp_eps );
    if ( ! (matpoint[where]->matmodel)->setTrialStrainIncr( incremental_strain)  )
      opserr << "TwentyNodeBrick_u_p_U::linearized_nodal_forces (tag: " << this->getTag() << "), not converged\n";
    Constitutive = (matpoint[where]->matmodel)->getTangentTensor();
          //    matpoint[where].setEPS( mmodel->getEPS() ); //Set the new EPState back
    //}
    //else if ( tmp_ndm ) { //Elastic case
    //    (matpoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
    //    Constitutive = (matpoint[where].p_matmodel)->getTangentTensor();
    //}
    //else {
                 //   g3ErrorHandler->fatal("TwentyNodeBrick_u_p_U::incremental_Update (tag: %d), could not getTangentTensor", this->getTag());
    //   exit(1);
    //}

    //Constitutive = ( matpoint[where].getEPS() )->getEep();
                //..//GPtangent_E[where].print("\n tangent E at GAUSS point \n");

                final_linearized_stress =
                  Constitutive("ijkl") * incremental_strain("kl");

                // nodal forces See Zienkievicz part 1 pp 108
                linearized_nodal_forces = linearized_nodal_forces +
                          dhGlobal("ib")*final_linearized_stress("ab")*weight;
                //::::::                   nodal_forces.print("nf","\n\n Nodal Forces \n");

              }
          }
      }


    return linearized_nodal_forces;

  }

////#############################################################################
//double TwentyNodeBrick_u_p_U::get_first_etacone()
//  {
//    double ret = matpoint[0].etacone();
//
//    return ret;
//
//  }
//

//=============================================================================
int TwentyNodeBrick_u_p_U::sendSelf (int commitTag, Channel &theChannel)
{
     // Not implemtented yet
     return 0;
}

//=============================================================================
int TwentyNodeBrick_u_p_U::recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
     // Not implemtented yet
     return 0;
}


//=============================================================================
int TwentyNodeBrick_u_p_U::displaySelf (Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = nd1Ptr->getCrds();
    const Vector &end2Crd = nd2Ptr->getCrds();
    const Vector &end3Crd = nd3Ptr->getCrds();
    const Vector &end4Crd = nd4Ptr->getCrds();
    const Vector &end5Crd = nd5Ptr->getCrds();
    const Vector &end6Crd = nd6Ptr->getCrds();
    const Vector &end7Crd = nd7Ptr->getCrds();
    const Vector &end8Crd = nd8Ptr->getCrds();
/*    const Vector &end9Crd = nd9Ptr->getCrds();
    const Vector &end10Crd = nd10Ptr->getCrds();
    const Vector &end11Crd = nd11Ptr->getCrds();
    const Vector &end12Crd = nd12Ptr->getCrds();
    const Vector &end13Crd = nd13Ptr->getCrds();
    const Vector &end14Crd = nd14Ptr->getCrds();
    const Vector &end15Crd = nd15Ptr->getCrds();
    const Vector &end16Crd = nd16Ptr->getCrds();
    const Vector &end17Crd = nd17Ptr->getCrds();
    const Vector &end18Crd = nd18Ptr->getCrds();
    const Vector &end19Crd = nd19Ptr->getCrds();
    const Vector &end20Crd = nd20Ptr->getCrds();  */

    const Vector &end1Disp = nd1Ptr->getDisp();
    const Vector &end2Disp = nd2Ptr->getDisp();
    const Vector &end3Disp = nd3Ptr->getDisp();
    const Vector &end4Disp = nd4Ptr->getDisp();
    const Vector &end5Disp = nd5Ptr->getDisp();
    const Vector &end6Disp = nd6Ptr->getDisp();
    const Vector &end7Disp = nd7Ptr->getDisp();
    const Vector &end8Disp = nd8Ptr->getDisp();
/*    const Vector &end9Disp = nd9Ptr->getDisp();
    const Vector &end10Disp = nd10Ptr->getDisp();
    const Vector &end11Disp = nd11Ptr->getDisp();
    const Vector &end12Disp = nd12Ptr->getDisp();
    const Vector &end13Disp = nd13Ptr->getDisp();
    const Vector &end14Disp = nd14Ptr->getDisp();
    const Vector &end15Disp = nd15Ptr->getDisp();
    const Vector &end16Disp = nd16Ptr->getDisp();
    const Vector &end17Disp = nd17Ptr->getDisp();
    const Vector &end18Disp = nd18Ptr->getDisp();
    const Vector &end19Disp = nd19Ptr->getDisp();
    const Vector &end20Disp = nd20Ptr->getDisp();  */

    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    static Vector v5(3);
    static Vector v6(3);
    static Vector v7(3);
    static Vector v8(3);
/*    static Vector v9(3);
    static Vector v10(3);
    static Vector v11(3);
    static Vector v12(3);
    static Vector v13(3);
    static Vector v14(3);
    static Vector v15(3);
    static Vector v16(3);
    static Vector v17(3);
    static Vector v18(3);
    static Vector v19(3);
    static Vector v20(3); */

    for (int i = 0; i < 2; i++)
    {
      v1(i) = end1Crd(i) + end1Disp(i)*fact;
      v2(i) = end2Crd(i) + end2Disp(i)*fact;
      v3(i) = end3Crd(i) + end3Disp(i)*fact;
      v4(i) = end4Crd(i) + end4Disp(i)*fact;
      v5(i) = end5Crd(i) + end5Disp(i)*fact;
      v6(i) = end6Crd(i) + end6Disp(i)*fact;
      v7(i) = end7Crd(i) + end7Disp(i)*fact;
      v8(i) = end8Crd(i) + end8Disp(i)*fact;
/*      v9(i) = end9Crd(i) + end9Disp(i)*fact;
      v10(i) = end10Crd(i) + end10Disp(i)*fact;
      v11(i) = end11Crd(i) + end11Disp(i)*fact;
      v12(i) = end12Crd(i) + end12Disp(i)*fact;
      v13(i) = end13Crd(i) + end13Disp(i)*fact;
      v14(i) = end14Crd(i) + end14Disp(i)*fact;
      v15(i) = end15Crd(i) + end15Disp(i)*fact;
      v16(i) = end16Crd(i) + end16Disp(i)*fact;
      v17(i) = end17Crd(i) + end17Disp(i)*fact;
      v18(i) = end18Crd(i) + end18Disp(i)*fact;
      v19(i) = end19Crd(i) + end19Disp(i)*fact;
      v20(i) = end20Crd(i) + end10Disp(i)*fact;     */
    }

    int error = 0;

    error += theViewer.drawLine (v1, v2, 1.0, 1.0);
    error += theViewer.drawLine (v2, v3, 1.0, 1.0);
    error += theViewer.drawLine (v3, v4, 1.0, 1.0);
    error += theViewer.drawLine (v4, v1, 1.0, 1.0);

    error += theViewer.drawLine (v5, v6, 1.0, 1.0);
    error += theViewer.drawLine (v6, v7, 1.0, 1.0);
    error += theViewer.drawLine (v7, v8, 1.0, 1.0);
    error += theViewer.drawLine (v8, v5, 1.0, 1.0);

    error += theViewer.drawLine (v1, v5, 1.0, 1.0);
    error += theViewer.drawLine (v2, v6, 1.0, 1.0);
    error += theViewer.drawLine (v3, v7, 1.0, 1.0);
    error += theViewer.drawLine (v4, v8, 1.0, 1.0);

    return error;

}

//=============================================================================
void TwentyNodeBrick_u_p_U::Print(OPS_Stream &s, int flag)
{
    //report(" TwentyNodeBrick_u_p_U ");
    s << "TwentyNodeBrick_u_p_U, element id:  " << this->getTag() << endln;
    s << "Connected external nodes:  " << connectedExternalNodes;

    int total_number_of_Gauss_points = r_integration_order*
                                       s_integration_order*
                                       t_integration_order;
    if ( total_number_of_Gauss_points != 0 )
      {
     nd1Ptr->Print(s);
     nd2Ptr->Print(s);
     nd3Ptr->Print(s);
     nd4Ptr->Print(s);

     nd5Ptr->Print(s);
     nd6Ptr->Print(s);
           nd7Ptr->Print(s);
     nd8Ptr->Print(s);

     nd9Ptr->Print(s);
     nd10Ptr->Print(s);
     nd11Ptr->Print(s);
     nd12Ptr->Print(s);

     nd13Ptr->Print(s);
     nd14Ptr->Print(s);
     nd15Ptr->Print(s);
     nd16Ptr->Print(s);

     nd17Ptr->Print(s);
     nd18Ptr->Print(s);
     nd19Ptr->Print(s);
     nd20Ptr->Print(s);
    }
    s << "Element mass density:  " << rho << endln << endln;
    s << "Material model:  " << endln;

    for( int GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
    {
      for( int GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
      {
        for( int GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
        {
           // this short routine is supposed to calculate position of
           // Gauss point from 3D array of short's
           short where =
           ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;

           s << "\n where = " << where << endln;
           s << " GP_c_r= " << GP_c_r << "GP_c_s = " << GP_c_s << " GP_c_t = " << GP_c_t << endln;
           matpoint[where]->report("Material Point\n");
           //GPstress[where].reportshort("stress at Gauss Point");
           //GPstrain[where].reportshort("strain at Gauss Point");
           //matpoint[where].report("Material model  at Gauss Point");
        }
      }
    }

}

//=============================================================================
Response * TwentyNodeBrick_u_p_U::setResponse (const char **argv, int argc, Information &eleInformation)
{
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new ElementResponse(this, 1, p);

    else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
    return new ElementResponse(this, 2, K);

  /*else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4)
      return theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, eleInfo);
      else
      return 0;
  }*/

    // otherwise response quantity is unknown for the quad class
    else
   return 0;
}
//=============================================================================

int TwentyNodeBrick_u_p_U::getResponse (int responseID, Information &eleInfo)
{
       switch (responseID) {

     case 1:
       return eleInfo.setVector(this->getResistingForce());

     /*case 2:
       return eleInfo.setMatrix(this->getTangentStiff());
      */
     default:
       return -1;
  }
     //return 0;
}




//#############################################################################
void TwentyNodeBrick_u_p_U::report(char * msg)
  {
    if ( msg ) ::printf("** %s",msg);
    ::printf("\n Element Number = %d\n", this->getTag() );
    ::printf("\n Number of nodes in a EightNodebrick_u_p_U = %d\n",
                                              nodes_in_brick);
    ::printf("\n Determinant of Jacobian (! ==0 before comp.) = %f\n",
                                                  determinant_of_Jacobian);

    ::printf("Node numbers \n");
    ::printf(
".....1.....2.....3.....4.....5.....6.....7.....8.....9.....0.....1.....2\n");
           for ( int i=0 ; i<20 ; i++ )
      //::printf("%6d",G_N_numbs[i]);
      ::printf("%6d",connectedExternalNodes(i));
    ::printf("\n");
    //           for ( int j=8 ; j<20 ; j++ )
    //             ::printf("%6d",G_N_numbs[j]);     // Commented by Xiaoyan
    ::printf("\n\n");

    //    ::printf("Node existance array \n");
    //           for ( int k=0 ; k<12 ; k++ )
    //             ::printf("%6d",node_existance[k]);
    //           ::printf("\n\n");          // Commented by Xiaoyan


    int total_number_of_Gauss_points = r_integration_order*
                                       s_integration_order*
                                       t_integration_order;
    if ( total_number_of_Gauss_points != 0 )
      {
           // report from Node class
           //for ( int in=0 ; in<8 ; in++ )
           //             (nodes[G_N_numbs[in]]).report("nodes from within element (first 8)\n");
           //Xiaoyan changed .report to . Print in above line 09/27/00
     //  (nodes[G_N_numbs[in]]).Print(opserr);

     nd1Ptr->Print(opserr);
     nd2Ptr->Print(opserr);
     nd3Ptr->Print(opserr);
     nd4Ptr->Print(opserr);

     nd5Ptr->Print(opserr);
     nd6Ptr->Print(opserr);
           nd7Ptr->Print(opserr);
     nd8Ptr->Print(opserr);

     nd9Ptr->Print(opserr);
     nd10Ptr->Print(opserr);
     nd11Ptr->Print(opserr);
      nd12Ptr->Print(opserr);

      nd13Ptr->Print(opserr);
     nd14Ptr->Print(opserr);
     nd15Ptr->Print(opserr);
     nd16Ptr->Print(opserr);

     nd17Ptr->Print(opserr);
     nd18Ptr->Print(opserr);
     nd19Ptr->Print(opserr);
     nd20Ptr->Print(opserr);

     //           for ( int jn=8 ; jn<20 ; jn++ )
           //             (nodes[G_N_numbs[jn]]).report("nodes from within element (last 12)\n");
           // Commented by Xiaoyan
      }

    ::printf("\n\nGauss-Legendre integration order\n");
    ::printf("Integration order in r direction = %d\n",r_integration_order);
    ::printf("Integration order in s direction = %d\n",s_integration_order);
    ::printf("Integration order in t direction = %d\n\n",t_integration_order);



    for( int GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        for( int GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            for( int GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                short where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;

                ::printf("\n\n----------------**************** where = %d \n", where);
                ::printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                            GP_c_r,GP_c_s,GP_c_t);
                matpoint[where]->report("Material Point\n");
                //GPstress[where].reportshort("stress at Gauss Point");
                //GPstrain[where].reportshort("strain at Gauss Point");
                //matpoint[where].report("Material model  at Gauss Point");
              }
          }
      }

  }


//#############################################################################
void TwentyNodeBrick_u_p_U::reportshort(char * msg)
  {
    if ( msg ) ::printf("** %s",msg);
    ::printf("\n Element Number = %d\n", this->getTag() );
    ::printf("\n Number of nodes in a TwentyNodeBrick_u_p_U = %d\n",
                                              nodes_in_brick);
    ::printf("\n Determinant of Jacobian (! ==0 before comp.) = %f\n",
                                                  determinant_of_Jacobian);

    ::printf("Node numbers \n");
    ::printf(
".....1.....2.....3.....4.....5.....6.....7.....8.....9.....0.....1.....2\n");
           for ( int i=0 ; i<20 ; i++ )
             //::printf("%6d",G_N_numbs[i]);
             ::printf( "%6d",connectedExternalNodes(i) );

     ::printf("\n");
           //           for ( int j=8 ; j<20 ; j++ )
           //             ::printf("%6d",G_N_numbs[j]);   //// Commented by Xiaoyan
           ::printf("\n\n");

           //    ::printf("Node existance array \n");
           //           for ( int k=0 ; k<12 ; k++ )
           //             ::printf("%6d",node_existance[k]);     // Commented by Xiaoyan
           ::printf("\n\n");

  }

//#############################################################################
void TwentyNodeBrick_u_p_U::reportPAK(char * msg)
  {
    if ( msg ) ::printf("%s",msg);
    ::printf("%10d   ",  this->getTag());
           for ( int i=0 ; i<8 ; i++ )
             ::printf( "%6d",connectedExternalNodes(i) );
             //::printf("%6d",G_N_numbs[i]);

    printf("\n");
  }

//#############################################################################
void TwentyNodeBrick_u_p_U::reportpqtheta(int GP_numb)
  {
    short where = GP_numb-1;
    matpoint[where]->reportpqtheta("");
  }

//#############################################################################
void TwentyNodeBrick_u_p_U::reportLM(char * msg)
  {
    if ( msg ) ::printf("%s",msg);
    ::printf("Element # %d, LM->", this->get_Brick_Number());
    for (int count = 0 ; count < 60 ; count++)
      {
        ::printf(" %d", LM[count]);
      }
    ::printf("\n");

  }

//#############################################################################
void TwentyNodeBrick_u_p_U::reportTensor(char * msg)
  {
        if ( msg ) ::printf("** %s\n",msg);

    // special case for 8 nodes only
    // special case for 8 nodes only
    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;

    short where = 0;

    // special case for 8 nodes only
    static const int dim[] = {3, 20}; // static-> see ARM pp289-290
    tensor NodalCoord(2, dim, 0.0);
    tensor matpointCoord(2, dim, 0.0);
    int h_dim[] = {60,3};   // Xiaoyan changed from {60,3} to {24,3} for 8 nodes
    tensor H(2, h_dim, 0.0);

    //for (int ncount = 1 ; ncount <= 8 ; ncount++ )
    ////  for (int ncount = 0 ; ncount <= 7 ; ncount++ )
    //  {
    //  //int global_node_number = get_global_number_of_node(ncount-1);
    //  // printf("global node num %d",global_node_number);
    //
    //    //   NodalCoord.val(1,ncount) = nodes[global_node_number].x_coordinate();
    //    //   NodalCoord.val(2,ncount) = nodes[global_node_number].y_coordinate();
    //    //   NodalCoord.val(3,ncount) = nodes[global_node_number].z_coordinate();
    //    // Xiaoyan changed to the following:  09/27/00
    //  Vector Coordinates = nodes[global_node_number].getCrds();
    //
    //    NodalCoord.val(1,ncount) = Coordinates(0);
    //    NodalCoord.val(2,ncount) = Coordinates(1);
    //    NodalCoord.val(3,ncount) = Coordinates(2);
    //printf("global point %d     x=%+.6e   y=%+.6e   z=%+.6e \n ", global_node_number,
    //                                                      NodalCoord.val(1,ncount),
    //                  NodalCoord.val(2,ncount),
    //                  NodalCoord.val(3,ncount));
    //}

    //Zhaohui using node pointers, which come from the Domain
    const Vector &nd1Crds = nd1Ptr->getCrds();
    const Vector &nd2Crds = nd2Ptr->getCrds();
    const Vector &nd3Crds = nd3Ptr->getCrds();
    const Vector &nd4Crds = nd4Ptr->getCrds();

    const Vector &nd5Crds = nd5Ptr->getCrds();
    const Vector &nd6Crds = nd6Ptr->getCrds();
    const Vector &nd7Crds = nd7Ptr->getCrds();
    const Vector &nd8Crds = nd8Ptr->getCrds();

    const Vector &nd9Crds = nd9Ptr->getCrds();
    const Vector &nd10Crds = nd10Ptr->getCrds();
    const Vector &nd11Crds = nd11Ptr->getCrds();
    const Vector &nd12Crds = nd12Ptr->getCrds();

    const Vector &nd13Crds = nd13Ptr->getCrds();
    const Vector &nd14Crds = nd14Ptr->getCrds();
    const Vector &nd15Crds = nd15Ptr->getCrds();
    const Vector &nd16Crds = nd16Ptr->getCrds();

    const Vector &nd17Crds = nd17Ptr->getCrds();
    const Vector &nd18Crds = nd18Ptr->getCrds();
    const Vector &nd19Crds = nd19Ptr->getCrds();
    const Vector &nd20Crds = nd20Ptr->getCrds();

    NodalCoord.val(1,1)=nd1Crds(0); NodalCoord.val(2,1)=nd1Crds(1); NodalCoord.val(3,1)=nd1Crds(2);
    NodalCoord.val(1,2)=nd2Crds(0); NodalCoord.val(2,2)=nd2Crds(1); NodalCoord.val(3,2)=nd2Crds(2);
    NodalCoord.val(1,3)=nd3Crds(0); NodalCoord.val(2,3)=nd3Crds(1); NodalCoord.val(3,3)=nd3Crds(2);
    NodalCoord.val(1,4)=nd4Crds(0); NodalCoord.val(2,4)=nd4Crds(1); NodalCoord.val(3,4)=nd4Crds(2);

    NodalCoord.val(1,5)=nd5Crds(0); NodalCoord.val(2,5)=nd5Crds(1); NodalCoord.val(3,5)=nd5Crds(2);
    NodalCoord.val(1,6)=nd6Crds(0); NodalCoord.val(2,6)=nd6Crds(1); NodalCoord.val(3,6)=nd6Crds(2);
    NodalCoord.val(1,7)=nd7Crds(0); NodalCoord.val(2,7)=nd7Crds(1); NodalCoord.val(3,7)=nd7Crds(2);
    NodalCoord.val(1,8)=nd8Crds(0); NodalCoord.val(2,8)=nd8Crds(1); NodalCoord.val(3,8)=nd8Crds(2);

    NodalCoord.val(1,9)=nd9Crds(0); NodalCoord.val(2,9)=nd8Crds(1); NodalCoord.val(3,9)=nd9Crds(2);
    NodalCoord.val(1,10)=nd10Crds(0); NodalCoord.val(2,10)=nd10Crds(1); NodalCoord.val(3,10)=nd10Crds(2);
    NodalCoord.val(1,11)=nd11Crds(0); NodalCoord.val(2,11)=nd11Crds(1); NodalCoord.val(3,11)=nd11Crds(2);
    NodalCoord.val(1,12)=nd12Crds(0); NodalCoord.val(2,12)=nd12Crds(1); NodalCoord.val(3,12)=nd12Crds(2);

    NodalCoord.val(1,13)=nd13Crds(0); NodalCoord.val(2,13)=nd13Crds(1); NodalCoord.val(3,13)=nd13Crds(2);
    NodalCoord.val(1,14)=nd14Crds(0); NodalCoord.val(2,14)=nd14Crds(1); NodalCoord.val(3,14)=nd14Crds(2);
    NodalCoord.val(1,15)=nd15Crds(0); NodalCoord.val(2,15)=nd15Crds(1); NodalCoord.val(3,15)=nd15Crds(2);
    NodalCoord.val(1,16)=nd16Crds(0); NodalCoord.val(2,16)=nd16Crds(1); NodalCoord.val(3,16)=nd16Crds(2);

    NodalCoord.val(1,17)=nd17Crds(0); NodalCoord.val(2,17)=nd17Crds(1); NodalCoord.val(3,17)=nd17Crds(2);
    NodalCoord.val(1,18)=nd18Crds(0); NodalCoord.val(2,18)=nd18Crds(1); NodalCoord.val(3,18)=nd18Crds(2);
    NodalCoord.val(1,19)=nd19Crds(0); NodalCoord.val(2,19)=nd19Crds(1); NodalCoord.val(3,19)=nd19Crds(2);
    NodalCoord.val(1,20)=nd20Crds(0); NodalCoord.val(2,20)=nd20Crds(1); NodalCoord.val(3,20)=nd20Crds(2);


    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates

               H = H_3D(r,s,t);

          for (int encount=1 ; encount <= 8 ; encount++ )
                //         for (int encount=0 ; encount <= 7 ; encount++ )
           {
                  //  matpointCoord.val(1,where+1) =+NodalCoord.val(1,where+1) * H.val(encount*3-2,1);
                  //  matpointCoord.val(2,where+1) =+NodalCoord.val(2,where+1) * H.val(encount*3-1,2);
                  //  matpointCoord.val(3,where+1) =+NodalCoord.val(3,where+1) * H.val(encount*3-0,3);
                  matpointCoord.val(1,where+1) +=NodalCoord.val(1,encount) * H.val(encount*3-2,1);
                  //::printf("-- NO nodal, H_val :%d %+.2e %+.2e %+.5e\n", encount,NodalCoord.val(1,encount),H.val(encount*3-2,1),matpointCoord.val(1,where+1) );
                  matpointCoord.val(2,where+1) +=NodalCoord.val(2,encount) * H.val(encount*3-1,2);
                  matpointCoord.val(3,where+1) +=NodalCoord.val(3,encount) * H.val(encount*3-0,3);

      }

    ::printf("gauss point# %d   %+.6e %+.6e %+.6e \n", where+1,
                                                       matpointCoord.val(1,where+1),
                                                       matpointCoord.val(2,where+1),
                                                       matpointCoord.val(3,where+1));

    //matpoint[where].reportTensor("");


              }
          }
      }

 }


////#############################################################################

//#############################################################################
//void TwentyNodeBrick_u_p_U::reportTensor(char * msg)
// ZHaohui added to print gauss point coord. to file fp

void TwentyNodeBrick_u_p_U::reportTensorF(FILE * fp)
  {
    //if ( msg ) ::printf("** %s\n",msg);

    // special case for 8 nodes only
    // special case for 8 nodes only
    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;

    short where = 0;

    // special case for 8 nodes only
    static const int dim[] = {3, 20}; // static-> see ARM pp289-290
    tensor NodalCoord(2, dim, 0.0);
    tensor matpointCoord(2, dim, 0.0);
    int h_dim[] = {60,3};  // Xiaoyan changed from {60,3} to {24,3} for 8 nodes

    tensor H(2, h_dim, 0.0);

    //for (int ncount = 1 ; ncount <= 8 ; ncount++ )
    //  // for (int ncount = 0 ; ncount <= 7 ; ncount++ )
    //  {
    //  int global_node_number = get_global_number_of_node(ncount-1);
    //  // printf("global node num %d",global_node_number);
    //
    //    //        NodalCoord.val(1,ncount) = nodes[global_node_number].x_coordinate();
    //    //        NodalCoord.val(2,ncount) = nodes[global_node_number].y_coordinate();
    //    //        NodalCoord.val(3,ncount) = nodes[global_node_number].z_coordinate();
    //    // Xiaoyan changed to the following:  09/27/00
    //  Vector Coordinates = nodes[global_node_number].getCrds();
    //    NodalCoord.val(1,ncount) = Coordinates(0);
    //    NodalCoord.val(2,ncount) = Coordinates(1);
    //    NodalCoord.val(3,ncount) = Coordinates(2);
    //printf("global point %d     x=%+.6e   y=%+.6e   z=%+.6e \n ", global_node_number,
    //                                                      NodalCoord.val(1,ncount),
    //                  NodalCoord.val(2,ncount),
    //                  NodalCoord.val(3,ncount));
    //  }

    //Zhaohui using node pointers, which come from the Domain
    const Vector &nd1Crds = nd1Ptr->getCrds();
    const Vector &nd2Crds = nd2Ptr->getCrds();
    const Vector &nd3Crds = nd3Ptr->getCrds();
    const Vector &nd4Crds = nd4Ptr->getCrds();

    const Vector &nd5Crds = nd5Ptr->getCrds();
    const Vector &nd6Crds = nd6Ptr->getCrds();
    const Vector &nd7Crds = nd7Ptr->getCrds();
    const Vector &nd8Crds = nd8Ptr->getCrds();

    const Vector &nd9Crds  =  nd9Ptr->getCrds();
    const Vector &nd10Crds = nd10Ptr->getCrds();
    const Vector &nd11Crds = nd11Ptr->getCrds();
    const Vector &nd12Crds = nd12Ptr->getCrds();

    const Vector &nd13Crds = nd13Ptr->getCrds();
    const Vector &nd14Crds = nd14Ptr->getCrds();
    const Vector &nd15Crds = nd15Ptr->getCrds();
    const Vector &nd16Crds = nd16Ptr->getCrds();

    const Vector &nd17Crds = nd17Ptr->getCrds();
    const Vector &nd18Crds = nd18Ptr->getCrds();
    const Vector &nd19Crds = nd19Ptr->getCrds();
    const Vector &nd20Crds = nd20Ptr->getCrds();

    NodalCoord.val(1,1)=nd1Crds(0); NodalCoord.val(2,1)=nd1Crds(1); NodalCoord.val(3,1)=nd1Crds(2);
    NodalCoord.val(1,2)=nd2Crds(0); NodalCoord.val(2,2)=nd2Crds(1); NodalCoord.val(3,2)=nd2Crds(2);
    NodalCoord.val(1,3)=nd3Crds(0); NodalCoord.val(2,3)=nd3Crds(1); NodalCoord.val(3,3)=nd3Crds(2);
    NodalCoord.val(1,4)=nd4Crds(0); NodalCoord.val(2,4)=nd4Crds(1); NodalCoord.val(3,4)=nd4Crds(2);

    NodalCoord.val(1,5)=nd5Crds(0); NodalCoord.val(2,5)=nd5Crds(1); NodalCoord.val(3,5)=nd5Crds(2);
    NodalCoord.val(1,6)=nd6Crds(0); NodalCoord.val(2,6)=nd6Crds(1); NodalCoord.val(3,6)=nd6Crds(2);
    NodalCoord.val(1,7)=nd7Crds(0); NodalCoord.val(2,7)=nd7Crds(1); NodalCoord.val(3,7)=nd7Crds(2);
    NodalCoord.val(1,8)=nd8Crds(0); NodalCoord.val(2,8)=nd8Crds(1); NodalCoord.val(3,8)=nd8Crds(2);

    NodalCoord.val(1,9)=nd9Crds(0); NodalCoord.val(2,9)=nd9Crds(1); NodalCoord.val(3,9)=nd9Crds(2);
    NodalCoord.val(1,10)=nd10Crds(0); NodalCoord.val(2,10)=nd10Crds(1); NodalCoord.val(3,10)=nd10Crds(2);
    NodalCoord.val(1,11)=nd11Crds(0); NodalCoord.val(2,11)=nd11Crds(1); NodalCoord.val(3,11)=nd11Crds(2);
    NodalCoord.val(1,12)=nd12Crds(0); NodalCoord.val(2,12)=nd12Crds(1); NodalCoord.val(3,12)=nd12Crds(2);

    NodalCoord.val(1,13)=nd13Crds(0); NodalCoord.val(2,13)=nd13Crds(1); NodalCoord.val(3,13)=nd13Crds(2);
    NodalCoord.val(1,14)=nd14Crds(0); NodalCoord.val(2,14)=nd14Crds(1); NodalCoord.val(3,14)=nd14Crds(2);
    NodalCoord.val(1,15)=nd15Crds(0); NodalCoord.val(2,15)=nd15Crds(1); NodalCoord.val(3,15)=nd15Crds(2);
    NodalCoord.val(1,16)=nd16Crds(0); NodalCoord.val(2,16)=nd16Crds(1); NodalCoord.val(3,16)=nd16Crds(2);

    NodalCoord.val(1,17)=nd17Crds(0); NodalCoord.val(2,17)=nd17Crds(1); NodalCoord.val(3,17)=nd17Crds(2);
    NodalCoord.val(1,18)=nd18Crds(0); NodalCoord.val(2,18)=nd18Crds(1); NodalCoord.val(3,18)=nd18Crds(2);
    NodalCoord.val(1,19)=nd19Crds(0); NodalCoord.val(2,19)=nd19Crds(1); NodalCoord.val(3,19)=nd19Crds(2);
    NodalCoord.val(1,20)=nd20Crds(0); NodalCoord.val(2,20)=nd20Crds(1); NodalCoord.val(3,20)=nd20Crds(2);


    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates

               H = H_3D(r,s,t);

          for (int encount=1 ; encount <= 8 ; encount++ )
                //         for (int encount=0 ; encount <= 7 ; encount++ )
         {
                  //  matpointCoord.val(1,where+1) =+NodalCoord.val(1,where+1) * H.val(encount*3-2,1);
                  //  matpointCoord.val(2,where+1) =+NodalCoord.val(2,where+1) * H.val(encount*3-1,2);
                  //  matpointCoord.val(3,where+1) =+NodalCoord.val(3,where+1) * H.val(encount*3-0,3);
                  matpointCoord.val(1,where+1) +=NodalCoord.val(1,encount) * H.val(encount*3-2,1);
                  //::printf("-- NO nodal, H_val :%d %+.2e %+.2e %+.5e\n", encount,NodalCoord.val(1,encount),H.val(encount*3-2,1),matpointCoord.val(1,where+1) );
                  matpointCoord.val(2,where+1) +=NodalCoord.val(2,encount) * H.val(encount*3-1,2);
                  matpointCoord.val(3,where+1) +=NodalCoord.val(3,encount) * H.val(encount*3-0,3);

         }

    fprintf(fp, "gauss point# %d   %+.6e %+.6e %+.6e \n", where+1,
                                                          matpointCoord.val(1,where+1),
                                                          matpointCoord.val(2,where+1),
                                                          matpointCoord.val(3,where+1));

    //matpoint[where].reportTensor("");


              }
          }
      }

 }


#endif
