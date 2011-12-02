///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              EightNodeBrick_u_p_U.cpp
// CLASS:             EightNodeBrick_u_p_U
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
// DATE:              Aug. 2001
// UPDATE HISTORY:    Modified from EightNodeBrick.cpp  reorganized a lot by Xiaoyan
//		      01/24/2002    Xiaoyan
//		      Add the permeability tensor and ks, kf to the constructor  Xiaoyan
//
//  "Coupled system": Solid and fluid coexist.
//                    u-- Solid displacement
//                    p-- Pore pressure
//                    U-- Absolute fluid displacement
//
//
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef EIGHTNODEBRICK_U_P_U_CPP
#define EIGHTNODEBRICK_U_P_U_CPP

#include <EightNodeBrick_u_p_U.h>
#include <ElementalLoad.h>	     //06/03/2002

#define FixedOrder 2

// Changed to static data members on 01/16/2002

Matrix EightNodeBrick_u_p_U::K(56, 56);      
Matrix EightNodeBrick_u_p_U::C(56, 56);      
Matrix EightNodeBrick_u_p_U::M(56, 56);      
Vector EightNodeBrick_u_p_U::p(56);	   
tensor EightNodeBrick_u_p_U::k(2,def_dim_2,0.0);
Node * EightNodeBrick_u_p_U::theNodes[8];

//=========================================================================
// Constructor. The dimension of K, C and M are(56,56)       Wxy 08/28/2001
//=========================================================================

EightNodeBrick_u_p_U::EightNodeBrick_u_p_U(int element_number,
                               int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
                               int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
                               NDMaterial * Globalmmodel, double b1, double b2,double b3,
			       double nn, double alf, double rs, double rf,
			       double permb_x,double permb_y,double permb_z, 
			       double kks, double kkf, double pp)	  
			       // wxy added rs and rf for the solid and fluid density    08/28/2001

			       //, EPState *InitEPS)  const char * type,
                               // Set it to 3 //int r_int_order, //int s_int_order, //int t_int_order,
			       //tensor * IN_tangent_E,  //stresstensor * INstress, //stresstensor * INiterative_stress, //double * IN_q_ast_iterative, //straintensor * INstrain):  __ZHaohui 09-29-2000
		               
  :Element(element_number, ELE_TAG_EightNodeBrick_u_p_U ),
   connectedExternalNodes(8), Q(56), bf(3), Ki(0),
  n(nn), alpha(alf), rho_s(rs), rho_f(rf), ks(kks), kf(kkf), pressure(pp)
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
//    k.print("k","\n test the permeability tensor k \n"); out 06/26/2002

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

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        double r = get_Gauss_p_c( r_integration_order, GP_c_r );
        double rw = get_Gauss_p_w( r_integration_order, GP_c_r );

        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            double s = get_Gauss_p_c( s_integration_order, GP_c_s );
            double sw = get_Gauss_p_w( s_integration_order, GP_c_s );

            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                double t = get_Gauss_p_c( t_integration_order, GP_c_t );
                double tw = get_Gauss_p_w( t_integration_order, GP_c_t );

                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;

 
                //cout << "where = " << where << endlnn;
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

      nodes_in_brick = 8;
      
}

////#############################################################################
//=========================================================================
// Default Constructor. The dimension of K, C and M are(56,56)       Wxy 08/28/2001
//=========================================================================

EightNodeBrick_u_p_U::EightNodeBrick_u_p_U ()
  :Element(0, ELE_TAG_EightNodeBrick_u_p_U ),
   connectedExternalNodes(8), Q(56), bf(3), Ki(0),
   n(0), alpha(1), rho_s(0.0),rho_f(0.0), ks(0.0), kf(0.0), pressure(0.0), mmodel(0)
{
     matpoint = 0;
}   

////#############################################################################
//=========================================================================
// Destructor.                                                             
//=========================================================================

EightNodeBrick_u_p_U::~EightNodeBrick_u_p_U ()
{
    
    int total_number_of_Gauss_points = r_integration_order*s_integration_order*t_integration_order;
    
    for (int i = 0; i < total_number_of_Gauss_points; i++) 
    {
	// Delete the NDMaterials at each integration point
	if (matpoint[i])
	    delete matpoint[i];    	
    }	
    
    // Delete the array of pointers to NDMaterial pointer arrays
    if (matpoint)
    	delete [] matpoint;
    
    //if (mmodel)
    //	delete [] mmodel;
    
    if (Ki != 0)
      delete Ki;

}
//=========================================================================
// Shape functions in "element golbal level". dimension are(24,3)          
// Since we define the mass Matrix Mf or Ms as four order tensor this 
// function will not used.  Just keep here. Wxy 09/26/2001
//=========================================================================
tensor EightNodeBrick_u_p_U::H_3D(double r1, double r2, double r3)
  {

    int dimension[] = {24,3}; // Xiaoyan changed from {60,3} to {24,3} for 8 nodes
                              // 3*8=24  07/12/00
    tensor H(2, dimension, 0.0);

    // influence of the node number 20
    //    H.val(58,1)=node_existance[20-1-8]*(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    H.val(59,2)=node_existance[20-1-8]*(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    H.val(60,3)=node_existance[20-1-8]*(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 19
    //    H.val(55,1)=node_existance[19-1-8]*(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    H.val(56,2)=node_existance[19-1-8]*(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    H.val(57,3)=node_existance[19-1-8]*(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 18
    //    H.val(52,1)=node_existance[18-1-8]*(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    H.val(53,2)=node_existance[18-1-8]*(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    H.val(54,3)=node_existance[18-1-8]*(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 17
    //    H.val(49,1)=node_existance[17-1-8]*(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    H.val(50,2)=node_existance[17-1-8]*(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    H.val(51,3)=node_existance[17-1-8]*(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;

    // influence of the node number 16
    //    H.val(46,1)=node_existance[16-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    H.val(47,2)=node_existance[16-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    H.val(48,3)=node_existance[16-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 15
    //    H.val(43,1)=node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    //    H.val(44,2)=node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    //    H.val(45,3)=node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    // influence of the node number 14
    //    H.val(40,1)=node_existance[14-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    H.val(41,2)=node_existance[14-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    H.val(42,3)=node_existance[14-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 13
    //    H.val(37,1)=node_existance[13-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;
    //    H.val(38,2)=node_existance[13-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;
    //    H.val(39,3)=node_existance[13-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;

    // influence of the node number 12
    //    H.val(34,1)=node_existance[12-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    H.val(35,2)=node_existance[12-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    H.val(36,3)=node_existance[12-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 11
    //    H.val(31,1)=node_existance[11-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    //    H.val(32,2)=node_existance[11-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    //    H.val(33,3)=node_existance[11-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    // influence of the node number 10
    //    H.val(28,1)=node_existance[10-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    H.val(29,2)=node_existance[10-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    H.val(30,3)=node_existance[10-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 9
    //    H.val(25,1)=node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;
    //    H.val(26,2)=node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;
    //    H.val(27,3)=node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;
    //

    // 9-20 nodes commented by Xiaoyan  07/12/00
    
    // influence of the node number 8
    //    H.val(22,1)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(15)+H.val(16)+H.val(20))/2.0;
    //    H.val(23,2)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(15)+H.val(16)+H.val(20))/2.0;
    //    H.val(24,3)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (H.val(15)+H.val(16)+H.val(20))/2.0;
    H.val(22,1)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0;// - (H.val(43,1)+H.val(48,3)+H.val(60,3))/2.0;
    H.val(23,2)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0;// - (H.val(43,1)+H.val(48,3)+H.val(60,3))/2.0;
    H.val(24,3)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0;// - (H.val(43,1)+H.val(48,3)+H.val(60,3))/2.0;
    // influence of the node number 7
    H.val(19,1)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0;// - (H.val(42,3)+H.val(43,1)+H.val(57,3))/2.0;
    H.val(20,2)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0;// - (H.val(42,3)+H.val(43,1)+H.val(57,3))/2.0;
    H.val(21,3)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0;// - (H.val(42,3)+H.val(43,1)+H.val(57,3))/2.0;
    // influence of the node number 6
    H.val(16,1)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0 ;//- (H.val(39,3)+H.val(42,3)+H.val(54,3))/2.0;
    H.val(17,2)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0 ;//- (H.val(39,3)+H.val(42,3)+H.val(54,3))/2.0;
    H.val(18,3)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0 ;//- (H.val(39,3)+H.val(42,3)+H.val(54,3))/2.0;
    // influence of the node number 5
    H.val(13,1)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0 ;//- (H.val(39,3)+H.val(48,3)+H.val(51,3))/2.0;
    H.val(14,2)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0 ;//- (H.val(39,3)+H.val(48,3)+H.val(51,3))/2.0;
    H.val(15,3)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0 ;//- (H.val(39,3)+H.val(48,3)+H.val(51,3))/2.0;

    // influence of the node number 4
    H.val(10,1)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0 ;//- (H.val(33,3)+H.val(36,3)+H.val(60,3))/2.0;
    H.val(11,2)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0 ;//- (H.val(33,3)+H.val(36,3)+H.val(60,3))/2.0;
    H.val(12,3)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0 ;//- (H.val(33,3)+H.val(36,3)+H.val(60,3))/2.0;
    // influence of the node number 3		        
    H.val(7,1)=(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0 ;//- (H.val(30,3)+H.val(33,3)+H.val(57,3))/2.0;
    H.val(8,2)=(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0 ;//- (H.val(30,3)+H.val(33,3)+H.val(57,3))/2.0;
    H.val(9,3)=(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0 ;//- (H.val(30,3)+H.val(33,3)+H.val(57,3))/2.0;
    // influence of the node number 2
    H.val(4,1)=(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0 ;//- (H.val(30,3)+H.val(54,3)+H.val(27,3))/2.0;
    H.val(5,2)=(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0 ;//- (H.val(30,3)+H.val(54,3)+H.val(27,3))/2.0;
    H.val(6,3)=(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0 ;//- (H.val(30,3)+H.val(54,3)+H.val(27,3))/2.0;
    // influence of the node number 1
    H.val(1,1)=(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0 ;//- (H.val(36,3)+H.val(51,3)+H.val(27,3))/2.0;
    H.val(2,2)=(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0 ;//- (H.val(36,3)+H.val(51,3)+H.val(27,3))/2.0;
    H.val(3,3)=(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0 ;//- (H.val(36,3)+H.val(51,3)+H.val(27,3))/2.0;

					       // The second part were commented by Xiaoyan
    //         double sum = 0;
    // 
    // 	for (int i=1; i<=60 ; i++)
    //           {
    // //  	    sum+=H.cval(i,1);
    // 	    for (int j=1; j<= 1; j++)
    // 	       {
    //        	          sum+=H.cval(i,1);
    // 	          ::printf( "  %+9.2e", H.cval(i,j) );
    // 	        }
    //            // ::printf( "  %d \n", i);
    // 	   }
    // 	    ::printf( " \n sum= %+6.2e\n", sum );
    

    //    printf("r1 = %lf, r2 = %lf, r3 = %lf\n", r1, r2, r3);
    //    H.print("h");

    return H;
  }
//==============================================================================
// Derivative of Shape functions in "element golbal level". dimension are(24,3) 
// Xiaoyna added this function  08/28/2001                                      
// Since we define the mass Matrix G1 or G2 as four order tensor this           
// function will not used.  Just keep here. Wxy 09/26/2001                      
//==============================================================================
tensor EightNodeBrick_u_p_U::dH_drst_at(double r1, double r2, double r3)
  {

    int dimension[] = {24,3}; // Xiaoyan changed from {60,3} to {24,3} for 8 nodes
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
    // 	for (int i=1; i<=60 ; i++)
    //           {
    // //  	    sum+=dH.cval(i,1);
    // 	    for (int j=1; j<= 1; j++)
    // 	       {
    //        	          sum+=dH.cval(i,1);
    // 	          ::printf( "  %+9.2e", dH.cval(i,j) );
    // 	        }
    //            // ::printf( "  %d \n", i);
    // 	   }
    // 	    ::printf( " \n sum= %+6.2e\n", sum );
    

    //    printf("r1 = %lf, r2 = %lf, r3 = %lf\n", r1, r2, r3);
    //    dH.print("h");

    return dH;
  }


//=========================================================================
// Shape functions in "element local level". dimension are(8)     Wxy 09/26/2001
//=========================================================================
tensor EightNodeBrick_u_p_U::interp_poli_at(double r1, double r2, double r3)
  {

    int dimension[] = {8};  // Xiaoyan changed from {20} to {8} for 8 nodes 07/12
    tensor h(1, dimension, 0.0);


    // influence of the node number 20
    //    h.val(20)=node_existance[20-1-8]*(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 19
    //    h.val(19)=node_existance[19-1-8]*(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 18
    //    h.val(18)=node_existance[18-1-8]*(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
    // influence of the node number 17
    //    h.val(17)=node_existance[17-1-8]*(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;

    // influence of the node number 16
    //    h.val(16)=node_existance[16-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 15
    //    h.val(15)=node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
    // influence of the node number 14
    //    h.val(14)=node_existance[14-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
    // influence of the node number 13
    //    h.val(13)=node_existance[13-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;

    // influence of the node number 12
    //    h.val(12)=node_existance[12-1-8]*(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 11
    //    h.val(11)=node_existance[11-1-8]*(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
    // influence of the node number 10
    //    h.val(10)=node_existance[10-1-8]*(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
    // influence of the node number 9
    //    h.val(9)=node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;

    // Commented by Xiaoyan

    // influence of the node number 8
    h.val(8)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0;// - (h.val(15)+h.val(16)+h.val(20))/2.0;
    // influence of the node number 7
    h.val(7)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0;// - (h.val(14)+h.val(15)+h.val(19))/2.0;
    // influence of the node number 6
    h.val(6)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0;// - (h.val(13)+h.val(14)+h.val(18))/2.0;
    // influence of the node number 5
    h.val(5)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0;// - (h.val(13)+h.val(16)+h.val(17))/2.0;

    // influence of the node number 4
    h.val(4)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0;// - (h.val(11)+h.val(12)+h.val(20))/2.0;
    // influence of the node number 3
    h.val(3)=(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0;// - (h.val(10)+h.val(11)+h.val(19))/2.0;
    // influence of the node number 2
    h.val(2)=(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0;// - (h.val(10)+h.val(18)+h.val(9))/2.0;
    // influence of the node number 1
    h.val(1)=(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0;// - (h.val(12)+h.val(17)+h.val(9))/2.0;
					    // The second part were commented by Xiaoyan 
					    // for 8 nodes

    //    printf("r1 = %lf, r2 = %lf, r3 = %lf\n", r1, r2, r3);
    //    h.print("h");

    return h;
  }

//#############################################################################
//=========================================================================
// Derivarite of Shape functions  dimension are(8,3)     Wxy 09/26/2001
//=========================================================================

tensor EightNodeBrick_u_p_U::dh_drst_at(double r1, double r2, double r3)
  {

    int dimensions[] = {8,3};  // Changed from{20,3} to {8,3} Xiaoyan 07/12
    tensor dh(2, dimensions, 0.0);


    // influence of the node number 20
    //    dh.val(20,1) =   node_existance[20-1-8]*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    dh.val(20,2) = - node_existance[20-1-8]*(1.0+r1)*(1.0-r3*r3)/4.0;
    //    dh.val(20,3) = - node_existance[20-1-8]*(1.0+r1)*(1.0-r2)*r3/2.0;
    // influence of the node number 19
    //    dh.val(19,1) = - node_existance[19-1-8]*(1.0-r2)*(1.0-r3*r3)/4.0;
    //    dh.val(19,2) = - node_existance[19-1-8]*(1.0-r1)*(1.0-r3*r3)/4.0;
    //    dh.val(19,3) = - node_existance[19-1-8]*(1.0-r1)*(1.0-r2)*r3/2.0;
    // influence of the node number 18
    //    dh.val(18,1) = - node_existance[18-1-8]*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    dh.val(18,2) =   node_existance[18-1-8]*(1.0-r1)*(1.0-r3*r3)/4.0;
    //    dh.val(18,3) = - node_existance[18-1-8]*(1.0-r1)*(1.0+r2)*r3/2.0;
    // influence of the node number 17
    //    dh.val(17,1) =   node_existance[17-1-8]*(1.0+r2)*(1.0-r3*r3)/4.0;
    //    dh.val(17,2) =   node_existance[17-1-8]*(1.0+r1)*(1.0-r3*r3)/4.0;
    //    dh.val(17,3) = - node_existance[17-1-8]*(1.0+r1)*(1.0+r2)*r3/2.0;

    // influence of the node number 16
    //    dh.val(16,1) =   node_existance[16-1-8]*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    dh.val(16,2) = - node_existance[16-1-8]*(1.0+r1)*r2*(1.0-r3)/2.0;
    //    dh.val(16,3) = - node_existance[16-1-8]*(1.0+r1)*(1.0-r2*r2)/4.0;
    // influnce of the node number 15
    //    dh.val(15,1) = - node_existance[15-1-8]*r1*(1.0-r2)*(1.0-r3)/2.0;
    //    dh.val(15,2) = - node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r3)/4.0;
    //    dh.val(15,3) = - node_existance[15-1-8]*(1.0-r1*r1)*(1.0-r2)/4.0;
    // influence of the node number 14
    //    dh.val(14,1) = - node_existance[14-1-8]*(1.0-r2*r2)*(1.0-r3)/4.0;
    //    dh.val(14,2) = - node_existance[14-1-8]*(1.0-r1)*r2*(1.0-r3)/2.0;
    //    dh.val(14,3) = - node_existance[14-1-8]*(1.0-r1)*(1.0-r2*r2)/4.0;
    // influence of the node number 13
    //    dh.val(13,1) = - node_existance[13-1-8]*r1*(1.0+r2)*(1.0-r3)/2.0;
    //    dh.val(13,2) =   node_existance[13-1-8]*(1.0-r1*r1)*(1.0-r3)/4.0;
    //    dh.val(13,3) = - node_existance[13-1-8]*(1.0-r1*r1)*(1.0+r2)/4.0;

    // influence of the node number 12
    //    dh.val(12,1) =   node_existance[12-1-8]*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    dh.val(12,2) = - node_existance[12-1-8]*(1.0+r1)*r2*(1.0+r3)/2.0;
    //    dh.val(12,3) =   node_existance[12-1-8]*(1.0+r1)*(1.0-r2*r2)/4.0;
    // influence of the node number 11
    //    dh.val(11,1) = - node_existance[11-1-8]*r1*(1.0-r2)*(1.0+r3)/2.0;
    //    dh.val(11,2) = - node_existance[11-1-8]*(1.0-r1*r1)*(1.0+r3)/4.0; // bug discovered 01 aug '95 2.0 -> 4.0
    //    dh.val(11,3) =   node_existance[11-1-8]*(1.0-r1*r1)*(1.0-r2)/4.0;
    // influence of the node number 10
    //    dh.val(10,1) = - node_existance[10-1-8]*(1.0-r2*r2)*(1.0+r3)/4.0;
    //    dh.val(10,2) = - node_existance[10-1-8]*(1.0-r1)*r2*(1.0+r3)/2.0;
    //    dh.val(10,3) =   node_existance[10-1-8]*(1.0-r1)*(1.0-r2*r2)/4.0;
    // influence of the node number 9
    //    dh.val(9,1) = - node_existance[9-1-8]*r1*(1.0+r2)*(1.0+r3)/2.0;
    //    dh.val(9,2) =   node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r3)/4.0;
    //    dh.val(9,3) =   node_existance[9-1-8]*(1.0-r1*r1)*(1.0+r2)/4.0;

    //   Commented by Xiaoyan for 8 nodes

    // influence of the node number 8
    dh.val(8,1)= (1.0-r2)*(1.0-r3)/8.0;// - (dh.val(15,1)+dh.val(16,1)+dh.val(20,1))/2.0;
    dh.val(8,2)=-(1.0+r1)*(1.0-r3)/8.0;// - (dh.val(15,2)+dh.val(16,2)+dh.val(20,2))/2.0;
    dh.val(8,3)=-(1.0+r1)*(1.0-r2)/8.0;// - (dh.val(15,3)+dh.val(16,3)+dh.val(20,3))/2.0;
    // influence of the node number 7
    dh.val(7,1)=-(1.0-r2)*(1.0-r3)/8.0;// - (dh.val(14,1)+dh.val(15,1)+dh.val(19,1))/2.0;
    dh.val(7,2)=-(1.0-r1)*(1.0-r3)/8.0;// - (dh.val(14,2)+dh.val(15,2)+dh.val(19,2))/2.0;
    dh.val(7,3)=-(1.0-r1)*(1.0-r2)/8.0;// - (dh.val(14,3)+dh.val(15,3)+dh.val(19,3))/2.0;
    // influence of the node number 6
    dh.val(6,1)=-(1.0+r2)*(1.0-r3)/8.0;// - (dh.val(13,1)+dh.val(14,1)+dh.val(18,1))/2.0;
    dh.val(6,2)= (1.0-r1)*(1.0-r3)/8.0;// - (dh.val(13,2)+dh.val(14,2)+dh.val(18,2))/2.0;
    dh.val(6,3)=-(1.0-r1)*(1.0+r2)/8.0;//- (dh.val(13,3)+dh.val(14,3)+dh.val(18,3))/2.0;
    // influence of the node number 5
    dh.val(5,1)= (1.0+r2)*(1.0-r3)/8.0;// - (dh.val(13,1)+dh.val(16,1)+dh.val(17,1))/2.0;
    dh.val(5,2)= (1.0+r1)*(1.0-r3)/8.0;// - (dh.val(13,2)+dh.val(16,2)+dh.val(17,2))/2.0;
    dh.val(5,3)=-(1.0+r1)*(1.0+r2)/8.0;// - (dh.val(13,3)+dh.val(16,3)+dh.val(17,3))/2.0;

    // influence of the node number 4
    dh.val(4,1)= (1.0-r2)*(1.0+r3)/8.0;// - (dh.val(11,1)+dh.val(12,1)+dh.val(20,1))/2.0;
    dh.val(4,2)=-(1.0+r1)*(1.0+r3)/8.0;// - (dh.val(11,2)+dh.val(12,2)+dh.val(20,2))/2.0;
    dh.val(4,3)= (1.0+r1)*(1.0-r2)/8.0;// - (dh.val(11,3)+dh.val(12,3)+dh.val(20,3))/2.0;
    // influence of the node number 3
    dh.val(3,1)=-(1.0-r2)*(1.0+r3)/8.0;// - (dh.val(10,1)+dh.val(11,1)+dh.val(19,1))/2.0;
    dh.val(3,2)=-(1.0-r1)*(1.0+r3)/8.0;// - (dh.val(10,2)+dh.val(11,2)+dh.val(19,2))/2.0;
    dh.val(3,3)= (1.0-r1)*(1.0-r2)/8.0;// - (dh.val(10,3)+dh.val(11,3)+dh.val(19,3))/2.0;
    // influence of the node number 2
    dh.val(2,1)=-(1.0+r2)*(1.0+r3)/8.0;// - (dh.val(10,1)+dh.val(18,1)+dh.val(9,1))/2.0;
    dh.val(2,2)= (1.0-r1)*(1.0+r3)/8.0;// - (dh.val(10,2)+dh.val(18,2)+dh.val(9,2))/2.0;
    dh.val(2,3)= (1.0-r1)*(1.0+r2)/8.0;// - (dh.val(10,3)+dh.val(18,3)+dh.val(9,3))/2.0;
    // influence of the node number 1
    dh.val(1,1)= (1.0+r2)*(1.0+r3)/8.0;// - (dh.val(12,1)+dh.val(17,1)+dh.val(9,1))/2.0;
    dh.val(1,2)= (1.0+r1)*(1.0+r3)/8.0;// - (dh.val(12,2)+dh.val(17,2)+dh.val(9,2))/2.0;
    dh.val(1,3)= (1.0+r1)*(1.0+r2)/8.0;//- (dh.val(12,3)+dh.val(17,3)+dh.val(9,3))/2.0;
				       // Commented by Xiaoyan
    return dh;
  }

//CE Dynamic Allocation for brick3d


//=========================================================================
// Permeability tensor  dimension are(3,3)     Wxy 09/26/2001              
//=========================================================================

// Xiaoyan added this function just want to test program. the values of k are not correct. 08/28/2001
/*tensor EightNodeBrick_u_p_U::k_at(double r1, double r2, double r3)
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
//Finite_Element * EightNodeBrick_u_p_U::new_el(int total)
//  {
//    EightNodeBrick_u_p_U *el_p;
//    el_p = new EightNodeBrick_u_p_U[total];
//    //DB//-------------------------------------------
//    //DB    for ( int i=0 ; i<total ; i++ )
//    //DB      {
//    //DB        el_p[i].report("derived EightNodeBrick_u_p_U\n");
//    //DB      }
//    //DB//-------------------------------------------
//    return el_p;
//  }
//=========================================================================
// Definition of Stiffness tensor Kep(8,3,3,8)     Wxy 09/26/2001          
//=========================================================================

tensor EightNodeBrick_u_p_U::getStiffnessTensorKep()
  {
    int K_dim[] = {8,3,3,8};  
			      
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

    int dh_dim[] = {8,3};  // Xiaoyan changed from {20,3} to {8,3}
    tensor dh(2, dh_dim, 0.0);

    //    tensor Constitutive( 4, def_dim_4, 0.0);
    tensor Constitutive;

    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {8,3};   // Xiaoyan changed from {20,3} to {8,3}
    tensor incremental_displacements(2,disp_dim,0.0); // \Delta u

    straintensor incremental_strain;
//    straintensor total_strain_at_GP;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
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
                dhGlobal = dh("ij") * JacobianINV("kj");
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
// Definition of Stiffness tensor G1(8,3,8)     Wxy 09/26/2001             
//=========================================================================
tensor EightNodeBrick_u_p_U::getStiffnessTensorG1()  //(double rho_s, double n,)
  {
    //int M_dim[] = {8,3,3,8}; 
    int G_dim[] = {8,3,8};    
    tensor G1(3,G_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}

    tensor dh(2, dh_dim, 0.0);

    int hp_dim[] = {8};   	// Xiaoyan changed from {60,3} to {24,3}
   tensor hp(1, hp_dim,0.0);     // Shape function.

    //int h_dim[] = {8,3};	// Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
//    tensor dH(2, h_dim, 0.0);
//    tensor Hp(1, Hp_dim,0.0);    // Xiaoyan added 08/27/2001

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
    tensor dhGlobal;

//    double RHO;
//    RHO= rho; 	 //global
//    double RHO_F=rho_f;
    double N=n;
    double ALPHA=alpha;

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
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
		hp= interp_poli_at(r,s,t);  // Assume the shape function of p (pressure) is the same as u. 
		                            // Xiaoyan 09/20/01
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // 		Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // 		printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates

	  
                //weight
                weight = rw * sw * tw * det_of_Jacobian;
  	                       
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

		G1 = G1 + dhGlobal("Ki")*hp("L")*weight*(ALPHA-N);
	       //	printf("\n +++++++++++++++++++++++++ \n\n");
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
// Definition of Stiffness tensor G2(8,3,8)     Wxy 09/26/2001             
//=========================================================================
tensor EightNodeBrick_u_p_U::getStiffnessTensorG2()  //(double rho_s, double n,)
  {
    //int M_dim[] = {8,3,3,8}; 
    int G_dim[] = {8,3,8};    
    tensor G2(3,G_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}

    tensor dhU(2, dh_dim, 0.0);

    int hp_dim[] = {8};   	// Xiaoyan changed from {60,3} to {24,3}
   tensor hp(1, hp_dim,0.0);     // Shape function.

    //int h_dim[] = {8,3};	// Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
//    tensor dH(2, h_dim, 0.0);
//    tensor Hp(1, Hp_dim,0.0);    // Xiaoyan added 08/27/2001

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
    tensor dhGlobal;

//    double RHO;
//    RHO= rho; 	 //global
//    double RHO_F=rho_f;
    double N=n;
//    double ALPHA=alpha;

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
                // derivatives of local coordinates with respect to local coordinates
                dhU = dh_drst_at(r,s,t);  // Assume  the shape function of u and U are same. Xiaoyan 09/20/01
		hp= interp_poli_at(r,s,t); // Assume  the shape function of p is same as the u's. Xiaoyan 09/20/01
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dhU);
                // 		Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dhU);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // 		printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dhU("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates

                //weight
                weight = rw * sw * tw * det_of_Jacobian;
  	        //	printf("weight = %6.2e \n",weight);

		//M.print("M","BEFORE");
                
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

		G2 = G2 + dhGlobal("Ki")*hp("L")*weight * N ;
	       //	printf("\n +++++++++++++++++++++++++ \n\n");
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
// Definition of Stiffness tensor P(8,8)     Wxy 09/26/2001                
//=========================================================================
tensor EightNodeBrick_u_p_U::getStiffnessTensorP()  
  {
    //int M_dim[] = {8,3,3,8}; 
    int G_dim[] = {8,8};    
    tensor P(2,G_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}

    tensor dhU(2, dh_dim, 0.0);

    int hp_dim[] = {8};   	
   tensor hp(1, hp_dim,0.0);     // Shape function of pore presssure.

    //int h_dim[] = {8,3};	// Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
//    tensor dH(2, h_dim, 0.0);
//    tensor Hp(1, Hp_dim,0.0);    // Xiaoyan added 08/27/2001

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
    tensor dhGlobal;

//    double RHO;
//    RHO= rho; 	 //global
//    double RHO_F=rho_f;
    double N=n;
    double ALPHA=alpha;
    double KS=ks;
    double KF=kf;
    double QQ = N/KF+(ALPHA-N)/KS;
    double e;
    double nu;

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
                // derivatives of local coordinates with respect to local coordinates
                dhU = dh_drst_at(r,s,t);  // Assume  the shape function of u and U are same. Xiaoyan 09/20/01
		hp= interp_poli_at(r,s,t); // Assume  the shape function of p is same as the u's. Xiaoyan 09/20/01
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dhU);
                // 		Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dhU);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // 		printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dhU("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates

                //weight
                weight = rw * sw * tw * det_of_Jacobian;
  	        //	printf("weight = %6.2e \n",weight);

		//M.print("M","BEFORE");
                
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

//		e = (matpoint[where]->matmodel)->getE();
//		nu = (matpoint[where]->matmodel)->getnu();
//		KS=e/(3*(1-2*nu));
//		cout<<"KS="<<KS<<endln;
//		QQ= N/KF+(ALPHA-N)/KS;
		P = P + hp("K") * hp("L")* QQ*weight;
	       //	printf("\n +++++++++++++++++++++++++ \n\n");
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
// Definition of Mass tensor Ms(8,3,3,8)     Wxy 09/26/2001                
//=========================================================================
tensor EightNodeBrick_u_p_U::getMassTensorMs()  //(double rho_s, double n,)
  {
// add some lines to test whether H^T*H =? dh^t*\delta_{ij}*dh? 
// xiaoyan 06/02/2002

    int M_dim[] = {8,3,3,8}; 
    tensor Ms(4,M_dim,0.0);
//out06/03/2002    tensor M83(4,M_dim,0.0);          // 06/02/2002-------

//out06/03/2002    int M24_dim[]={24,24};	      // 06/02/2002-------
//out06/03/2002    tensor M24(2,M24_dim,0.0);	      // 06/02/2002-------



    tensor I2("I", 2, def_dim_2);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {8};	
    tensor H(1, h_dim, 0.0);
// add the following two lines
//out06/03/2002     int h24_dim[]={24,3};	       	// 06/02/2002-------
//out06/03/2002     tensor H24(2,h24_dim,0.0);		// 06/02/2002-------
    double det_of_Jacobian = 0.0;	

    tensor Jacobian;

//    rho=(1-n)*rho_s+n*rho_f;
//    double RHO;
//    RHO= rho; 	 //global
    double RHO_S=rho_s;
    double N=n;

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
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // 		Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // 		printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                //                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates


                // printf("\n\nIN THE MASS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                // printf("  Mass_Tensor \n");
                // printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //                             GP_c_r,GP_c_s,GP_c_t);
                // 
                //weight
                weight = rw * sw * tw * det_of_Jacobian;
                H = interp_poli_at(r,s,t);
		
//out06/03/2002 		H24=H_3D(r,s,t);             // 06/02/2002-------
//out06/03/2002 		M24=M24+H24("ib")*H24("kb")*((1-N)*RHO_S *weight); // 06/02/2002-------
//out06/03/2002 //		tensor tem=H("K")*I2("ij");  // 06/02/2002-------
//out06/03/2002 //		Ms=Ms+tem("Kij")*H("L")*((1-N)*RHO_S *weight);   // 06/02/2002-------
//out06/03/2002 //                //weight
//out06/03/2002 //                weight = rw * sw * tw * det_of_Jacobian;
  	        //	printf("weight = %6.2e \n",weight);

		//M.print("M","BEFORE");
                
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );
		 tensor t=H("K") * I2("ij");
                 Ms = Ms + t("Kij") * H("L") * ((1-N)*RHO_S *weight);
//		 M24 = M24 +H24("ib")*H24("kb") * ((1-N)*RHO_S *weight);
		 //Ms.printshort("M");
              }
          }
      }
    //Ms.printshort("M");
//out06/03/2002 //    M83.printshort("M83_short");	                       // 06/02/2002-------
//out06/03/2002     M24.printshort("M24_short");	                       // 06/02/2002-------
//out06/03/2002 //    matrix mM83 = mass_matrixMs(M83);                          // 06/02/2002-------
//out06/03/2002 //    mM83.write_standard("mM83", "4th tensor to mass matrix");  // 06/02/2002-------
//out06/03/2002     matrix mMs=mass_matrixMs(Ms);                              // 06/02/2002-------
//out06/03/2002     mMs.write_standard("mMs", "4th tensor to mass matrix");    // 06/02/2002-------
//out06/03/2002     matrix mM24=mass_matrix(M24);				// 06/02/2002-------
//out06/03/2002     mM24.write_standard("mM24", "4th tensor to mass matrix");    // 06/02/2002-------

    return Ms;				
  }

////#############################################################################
//=========================================================================
// Definition of Mass tensor Mf(8,3,3,8)     Wxy 09/26/2001                
//=========================================================================
tensor EightNodeBrick_u_p_U::getMassTensorMf()  //(double rho_s, double n,)
  {
    //int M_dim[] = {8,3,3,8}; 
    int M_dim[] = {8,3,3,8}; 
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

    int dh_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {8};	
    tensor H(1, h_dim, 0.0);	 // Shape fun. of u
    tensor HU(1, h_dim, 0.0);	 // Shape fun. of U

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

//    rho=(1-n)*rho_s+n*rho_f;
//    double RHO;
//    RHO= rho; 	 //global
//    double RHO_S=rho_s;
    double RHO_F=rho_f;
    double N=n;


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
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t);
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // 		Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // 		printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                //                dhGlobal = dh("ij") * JacobianINV("kj");
                // derivatives of local coordinates with respect to local coordinates


                HU = interp_poli_at(r,s,t);   // now assume HU_3D(r,s,t) is the same as H_3D(r,s,t), 
		                              // wxy 08/28/2001

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

  	        //	printf("weight = %6.2e \n",weight);

		//M.print("M","BEFORE");
                
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );
		tensor temp=HU("K")* I2("ij");
		Mf = Mf + temp("Kij") * HU("L")* N * RHO_F * weight;
	       //	printf("\n +++++++++++++++++++++++++ \n\n");
	      	//Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
//    matrix mMf=mass_matrixMs(Mf);                              // 06/02/2002-------
//    mMf.write_standard("mMf", "4th tensor to mass matrix");    // 06/02/2002-------

    return Mf;
  }




//=========================================================================
// Definition of Damping tensor C1(8,3,3,8)     Wxy 09/26/2001
//=========================================================================
tensor EightNodeBrick_u_p_U::getDampTensorC1()  //(double rho_s, double n,)
  {
    int C_dim[] = {8,3,3,8};    
    tensor C1(4,C_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {8};	
    //int h_dim[] = {8,3};	// Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
    tensor H(1, h_dim, 0.0);
//    int k_dim[]={3,3};	       // Xiaoyan added for permeability tensor 08/27/2001
//    tensor k(2,k_dim,0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

    double N=n;

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
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t); 
//		k=k_at(r,s,t);    // wxy added for permeability.
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // 		Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // 		printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
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

  	        //	printf("weight = %6.2e \n",weight);

                
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

		tensor k_inverse=k("mn").inverse();
		tensor temp=   H("K")* k_inverse("ij");
		C1 = C1 + temp("Kij") *H("L")*weight * N*N;
	       //	printf("\n +++++++++++++++++++++++++ \n\n");
	      	//Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return C1;
  }


//=========================================================================
// Definition of Damping tensor C2(8,3,3,8)     Wxy 09/26/2001
//=========================================================================
tensor EightNodeBrick_u_p_U::getDampTensorC2()  //(double rho_s, double n,)
  {
    //int M_dim[] = {8,3,3,8}; 
    int C_dim[] = {8,3,3,8};    
    tensor C2(4,C_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {8};	
    //int h_dim[] = {8,3};	// Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
    tensor H(1, h_dim, 0.0);
    tensor HU(1, h_dim, 0.0);
//   int k_dim[]={3,3};	       // Xiaoyan added for permeability tensor 08/27/2001
//    tensor k(2,k_dim,0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

    double N=n;


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
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t); 
//		k=k_at(r,s,t);    // wxy added for permeability.
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // 		Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // 		printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
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
  	        //	printf("weight = %6.2e \n",weight);

		//M.print("M","BEFORE");
                
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

		tensor k_inverse=k("mn").inverse();
		tensor temp=  H("L")* k_inverse("ij") ;
		C2 = C2 + temp("Lij") *HU("K")*weight * N*N;
	       //	printf("\n +++++++++++++++++++++++++ \n\n");
	      	//Mf.printshort("M");
              }
          }
      }
    //M = Mf;
    //Mf.printshort("M");
    return C2;
  }

 
//=========================================================================
// Definition of Damping tensor C3(8,3,3,8)     Wxy 09/26/2001             
//=========================================================================
tensor EightNodeBrick_u_p_U::getDampTensorC3()  //(double rho_s, double n,)
  {
    //int M_dim[] = {8,3,3,8}; 
    int C_dim[] = {8,3,3,8};    
    tensor C3(4,C_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}

    tensor dh(2, dh_dim, 0.0);

    int h_dim[] = {8};	
    //int h_dim[] = {8,3};	// Xiaoyan changed from {60,3} to {24,3}
    //int h_dim[] = {20,3};
    tensor HU(1, h_dim, 0.0);
//    int k_dim[]={3,3};	       // Xiaoyan added for permeability tensor 08/27/2001
//    tensor k(2,k_dim,0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

    double N=n;


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
                // derivatives of local coordinates with respect to local coordinates
                dh = dh_drst_at(r,s,t); 
//		k=k_at(r,s,t);    // wxy added for permeability.
                // Jacobian tensor ( matrix )
                Jacobian = Jacobian_3D(dh);
                // 		Jacobian.print("J","Jacobian");
                // Inverse of Jacobian tensor ( matrix )
                //                JacobianINV = Jacobian_3Dinv(dh);
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                // 		printf("det_of_Jacobian = %6.2e \n",det_of_Jacobian);
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
  	        //	printf("weight = %6.2e \n",weight);

		//M.print("M","BEFORE");
                
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

		tensor k_inverse=k("mn").inverse();
		tensor temp =  HU("L")* k_inverse("ij");
		C3 = C3 + temp("Lij") *HU("K")*weight * N*N;
	       //	printf("\n +++++++++++++++++++++++++ \n\n");
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
// Converting stiffness tensor to stiffness matrix K^ep 24*24              
//=========================================================================

matrix EightNodeBrick_u_p_U::stiffness_matrixKep(const tensor  Kep) 
						       
  {						  
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    matrix Kepmatrix(24,24,0.0);	  

    int Ki=0;
    int Kj=0;
    int i;
    for ( i=1 ; i<=8 ; i++ )  // i<=8 for 8 nodes
      {
        int j;
        for ( j=1 ; j<=8 ; j++ )  // i<=8 for 8 nodes
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
// Converting stiffness tensor to stiffness matrix G1 24*8                 
//=========================================================================

matrix EightNodeBrick_u_p_U::stiffness_matrixG1(const tensor  G1) 
						       
  {						  
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    matrix G1matrix(24,8,0.0);	  

    int Ki=0;
    int Kj=0;

    for ( int i=1 ; i<=8 ; i++ )  
      {
        for ( int k=1 ; k<=3 ; k++ )
          {
            for ( int l=1 ; l<=8 ; l++ )
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
// Converting stiffness tensor to stiffness matrix G2 24*8                 
//=========================================================================


matrix EightNodeBrick_u_p_U::stiffness_matrixG2(const tensor  G2) 
						       
  {						  
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    matrix G2matrix(24,8,0.0);	  

    int Ki=0;
    int Kj=0;

    for ( int i=1 ; i<=8 ; i++ )  // i<=8 for 8 nodes
      {
        for ( int k=1 ; k<=3 ; k++ )
          {
            for ( int l=1 ; l<=8 ; l++ )
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
// Converting stiffness tensor to stiffness matrix P 8*8                   
//=========================================================================

matrix EightNodeBrick_u_p_U::stiffness_matrixP(const tensor P) 
						       
  {						  
    matrix Pmatrix(8,8,0.0);	  

    for ( int i=1 ; i<=8 ; i++ )  // i<=8 for 8 nodes
      {
        for ( int j=1 ; j<=8 ; j++ )
          {
             Pmatrix.val( i , j ) = P.cval(i,j);
            
          }						 
      }
    return Pmatrix;
  }

//=====================================================================================
// Constructing the whole system K (including Kep, G1, G2 and P 56*56)   Wxy 09/26/2001
//=====================================================================================

const Matrix &EightNodeBrick_u_p_U::getTangentStiff () 
						       
  {						  
    tensor tKep = getStiffnessTensorKep();
    tensor tG1  = getStiffnessTensorG1();
    tensor tG2  = getStiffnessTensorG2();
    tensor tP   = getStiffnessTensorP();
    matrix Kep = stiffness_matrixKep(tKep);
    matrix G1  = stiffness_matrixG1(tG1);
    matrix G2  = stiffness_matrixG2(tG2);
    matrix P   = stiffness_matrixP(tP);
    matrix G1t=G1.transpose();
    matrix G2t=G2.transpose();
//    Kep.write_standard("Kep_8upu.dat", "stiffness part of KupU");
//    G1.write_standard("G1_8upu.dat", "stiffness part of G1upU");
//    G2.write_standard("G2_8upu.dat", "stiffness part of G2upU");
//    P.write_standard("P_8upu.dat", "stiffness part of PupU");

    int i;
			 for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
          {
	    for(int n=0; n<3; n++)
	      {
	        for(int m=0; m<3; m++)
		  {
 		    K(i*7+n ,  j*7+m)=Kep.val(3*i+n+1, 3*j+1+m);     // Add Kep to K
		  }
	      }
	  }
      }

    for ( i=0 ; i<8 ; i++ )  		
      {
        for(int n=0; n<3; n++)
          {
            for ( int j=0 ; j<8 ; j++ )
	      {
	         K(i*7+n ,  j*7+3)=-G1.val(3*i+n+1, j+1);     // Add -G1 to K
	      }
	  }
      }

    for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
           {
              for(int n=0; n<3; n++)
   	      {
   	         K(i*7+3 ,  j*7+n)=-G1t.val(i+1, 3*j+1+n);     // Add -G1T to K
	      }
	  }
      }


    for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
	  {
	     K(i*7+3 ,  j*7+3)=-P.val(i+1, j+1);     // Add -P to K
	  }
      }


    for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
           {
              for(int n=0; n<3; n++)
   	      {
   	         K(i*7+3 ,  j*7+n+4)=-G2t.val(i+1, 3*j+1+n);     // Add -G2T to K
	      }
	  }
      }

    for ( i=0 ; i<8 ; i++ )  		
      {
        for(int n=0; n<3; n++)
          {
            for ( int j=0 ; j<8 ; j++ )
	      {
	         K(i*7+n+4 ,  j*7+3)=-G2.val(3*i+n+1, j+1);     // Add -G2 to K
	      }
	  }
      }

//     ofstream outK("K8.dat");	// want to check the whole K
//     K.Output(outK);

     return K;
  }

 	   
//=========================================================================
// Converting damping tensor to damping matrix C1(24*24)                   
//=========================================================================

matrix EightNodeBrick_u_p_U::damping_matrixC1(const tensor  C1)
						       
  {
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    matrix C1matrix(24,24,0.0);	  

    int Ki=0;
    int Kj=0;

    for ( int i=1 ; i<=8 ; i++ )  //  i<=8 for 8 nodes
      {
        for ( int j=1 ; j<=8 ; j++ )  //  i<=8 for 8 nodes
          {
            for ( int k=1 ; k<=3 ; k++ )
              {
                for ( int l=1 ; l<=3 ; l++ )
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
// Converting damping tensor to damping matrix  C2(24*240                  
//=========================================================================

matrix EightNodeBrick_u_p_U::damping_matrixC2(const tensor  C2) 
  {						   
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    matrix C2matrix(24,24,0.0);	  

    int Ki=0;
    int Kj=0;

    for ( int i=1 ; i<=8 ; i++ )  //  i<=8 for 8 nodes
      {
        for ( int j=1 ; j<=8 ; j++ )  //  i<=8 for 8 nodes
          {
            for ( int k=1 ; k<=3 ; k++ )
              {
                for ( int l=1 ; l<=3 ; l++ )
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
// Converting damping tensor to damping matrix C3(24*24)                   
//=========================================================================

matrix EightNodeBrick_u_p_U::damping_matrixC3(const tensor  C3)   
  {						     
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    matrix C3matrix(24,24,0.0);	  

    int Ki=0;
    int Kj=0;

    for ( int i=1 ; i<=8 ; i++ )  //  i<=8 for 8 nodes
      {
        for ( int j=1 ; j<=8 ; j++ )  //  i<=8 for 8 nodes
          {
            for ( int k=1 ; k<=3 ; k++ )
              {
                for ( int l=1 ; l<=3 ; l++ )
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
// Constructing Damping matrix of the whole system C(56*56) (including C1 C2 and C3) 
//=========================================================================

const Matrix &EightNodeBrick_u_p_U::getDamp () 
						       
  {						  

    tensor tC1  = getDampTensorC1();
    tensor tC2  = getDampTensorC2();
    tensor tC3  = getDampTensorC3();
    matrix C1 = damping_matrixC1(tC1);
    matrix C2 = damping_matrixC2(tC2);
    matrix C3 = damping_matrixC3(tC3);
//    C1.write_standard("C1", " damping  matrix");    // 06/04/2002
//    C2.write_standard("C2", " damping  matrix");    // 06/04/2002
//    C3.write_standard("C3", " damping  matrix");    // 06/04/2002

    matrix C2t=C2.transpose();

    int i;
    for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
          {
	    for(int n=0; n<3; n++)
	      {
	        for(int m=0; m<3; m++)
		  {
 		    C(i*7+n ,  j*7+m)=C1.val(3*i+n+1, 3*j+1+m);     // Add C1 to C
		  }
	      }
	  }
      }

     for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
          {
	    for(int n=0; n<3; n++)
	      {
	        for(int m=0; m<3; m++)
		  {
 		    C(i*7+n ,  j*7+m+4)=-C2.val(3*i+n+1, 3*j+1+m);     // Add -C2 to C
		  }
	      }
	  }
      }

    for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
          {
	    for(int n=0; n<3; n++)
	      {
	        for(int m=0; m<3; m++)
		  {
 		    C(i*7+n+4 ,  j*7+m)=-C2t.val(3*i+n+1, 3*j+1+m);     // Add -C2t to C
		  }
	      }
	  }
      }

    for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
          {
	    for(int n=0; n<3; n++)
	      {
	        for(int m=0; m<3; m++)
		  {
 		    C(i*7+n+4,  j*7+m+4)=C3.val(3*i+n+1, 3*j+1+m);     // Add C3 to C
		  }
	      }
	  }
      }

    return C;
  }

//=========================================================================
// Converting mass tensor(8,3,3,8) to mass matrix Ms ___Xiaoyan 08/27/2001          
//=========================================================================
matrix EightNodeBrick_u_p_U::mass_matrixMs(const tensor  Ms)
  {
    matrix Msmatrix(24,24,0.0);	  

    int Ki=0;
    int Kj=0;

    for ( int i=1 ; i<=8 ; i++ )  //  i<=8 for 8 nodes
      {
        for ( int j=1 ; j<=8 ; j++ )  // i<=8 for 8 nodes
          {
            for ( int k=1 ; k<=3 ; k++ )
              {
                for ( int l=1 ; l<=3 ; l++ )
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
// Converting mass tensor(8,3,3,8) to mass matrix Mf___Xiaoyan 08/27/2001           
//=========================================================================

matrix EightNodeBrick_u_p_U::mass_matrixMf(const tensor  Mf)
  {
    matrix Mfmatrix(24,24,0.0);	  

    int Ki=0;
    int Kj=0;

    for ( int i=1 ; i<=8 ; i++ )  // i<=8 for 8 nodes
      {
        for ( int j=1 ; j<=8 ; j++ )  //  i<=8 for 8 nodes
          {
            for ( int k=1 ; k<=3 ; k++ )
              {
                for ( int l=1 ; l<=3 ; l++ )
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
//=========================================================================
// Converting mass tensor(24,24) to mass matrix Mf___Xiaoyan 08/27/2001           
//=========================================================================
 
matrix EightNodeBrick_u_p_U::mass_matrix(const tensor  M)
  {
    //    int K_dim[] = {20,3,3,20};
    //    tensor K(4,K_dim,0.0);
    matrix Mmatrix(24,24,0.0);       
				     
    for ( int i=1 ; i<=24 ; i++ )    
      {				     
        for ( int j=1 ; j<=24 ; j++ )
          {
             Mmatrix.val( i , j ) = M.cval(i,j);
          }
      }
    return Mmatrix;
  }


//#############################################################################  
//=========================================================================
// Constructing Mass matrix of the whole system M (including Ms and Mf 56*56    
//=========================================================================

const Matrix &EightNodeBrick_u_p_U::getMass () 
 {

    tensor tMs  = getMassTensorMs();
    tensor tMf  = getMassTensorMf();
    matrix  Ms =  mass_matrixMs(tMs);
    matrix  Mf =  mass_matrixMf(tMf);

    int i;
			 for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
          {
	    for(int n=0; n<3; n++)
	      {
	        for(int m=0; m<3; m++)
		  {
 		    M(i*7+n ,  j*7+m)=Ms.val(3*i+n+1, 3*j+1+m);     // Add Ms to M
		  }
	      }
	  }
      }
    for ( i=0 ; i<8 ; i++ )  		
      {
        for ( int j=0 ; j<8 ; j++ )
          {
	    for(int n=0; n<3; n++)
	      {
	        for(int m=0; m<3; m++)
		  {
 		    M(i*7+n+4 ,  j*7+m+4)=Mf.val(3*i+n+1, 3*j+1+m);     // Add Mf to M
		  }
	      }
	  }
      }

//     ofstream outM("M8.dat");	// want to check the whole M
//     M.Output(outM);

     return M;
  }

//#############################################################################  
 

//=========================================================================
// Jacobian tensor J = dx/dr=dN/dr*x_i=dh*coordinate  wxy 09/26/2001          
//=========================================================================
tensor EightNodeBrick_u_p_U::Jacobian_3D(tensor dh)
  {                       
     //       dh ( 20*3)  // dh(8*3) Xiaoyan
     tensor N_C = Nodal_Coordinates(); // 20*3	  // 8*3 Xiaoyan
     tensor Jacobian_3D = dh("ij") * N_C("ik");
     return Jacobian_3D;
  }

//#############################################################################
//=========================================================================
// Jacobian tensor J^(-1)  wxy 09/26/2001          
//=========================================================================
tensor EightNodeBrick_u_p_U::Jacobian_3Dinv(tensor dh)
  {                       
     //       dh ( 20*3)	  // dh(8*3) Xiaoyan  
     tensor N_C = Nodal_Coordinates(); // 20*3	  	  // 8*3 Xiaoyan   
     tensor Jacobian_3Dinv = (dh("ij") * N_C("ik")).inverse();
     return Jacobian_3Dinv;
  }


////#############################################################################
tensor EightNodeBrick_u_p_U::Nodal_Coordinates()
  {
    const int dimensions[] = {8,3};  // Xiaoyan changed from {20,3} to {8,3} for 8 nodes
    tensor N_coord(2, dimensions, 0.0);

    //for ( int i=0 ; i<8 ; i++ )	  // Xiaoyan changed from 20 to 8 for 8 nodes
    //  {
    //    //        N_coord.val(i+1,1) = nodes[ G_N_numbs[i] ].x_coordinate();
    //    //        N_coord.val(i+1,2) = nodes[ G_N_numbs[i] ].y_coordinate();
    //    //        N_coord.val(i+1,3) = nodes[ G_N_numbs[i] ].z_coordinate();
    //    // Xiaoyan changed to the following:  09/27/00
    //    /// LOOK WITH DDD
    //	Vector Coordinates = nodes[ G_N_numbs[i] ].getCrds();
    //    N_coord.val(i+1,1) = Coordinates(0);
    //    N_coord.val(i+1,2) = Coordinates(1);
    //    N_coord.val(i+1,3) = Coordinates(2);
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
    
    N_coord.val(1,1)=nd1Crds(0); N_coord.val(1,2)=nd1Crds(1); N_coord.val(1,3)=nd1Crds(2);
    N_coord.val(2,1)=nd2Crds(0); N_coord.val(2,2)=nd2Crds(1); N_coord.val(2,3)=nd2Crds(2);
    N_coord.val(3,1)=nd3Crds(0); N_coord.val(3,2)=nd3Crds(1); N_coord.val(3,3)=nd3Crds(2);
    N_coord.val(4,1)=nd4Crds(0); N_coord.val(4,2)=nd4Crds(1); N_coord.val(4,3)=nd4Crds(2);
    N_coord.val(5,1)=nd5Crds(0); N_coord.val(5,2)=nd5Crds(1); N_coord.val(5,3)=nd5Crds(2);
    N_coord.val(6,1)=nd6Crds(0); N_coord.val(6,2)=nd6Crds(1); N_coord.val(6,3)=nd6Crds(2);
    N_coord.val(7,1)=nd7Crds(0); N_coord.val(7,2)=nd7Crds(1); N_coord.val(7,3)=nd7Crds(2);
    N_coord.val(8,1)=nd8Crds(0); N_coord.val(8,2)=nd8Crds(1); N_coord.val(8,3)=nd8Crds(2);

    return N_coord;
  }

////#############################################################################
tensor EightNodeBrick_u_p_U::incr_dispDu()
  {
    const int dimensions[] = {8,3};  // Xiaoyan changed from {20,3} to {8,3} for 8 nodes
    tensor increment_dispDu(2, dimensions, 0.0);

    //for ( int i=0 ; i<8 ; i++ )	 // Xiaoyan changed from 20 to 8 for 8 nodes
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

    //if ( getTag() == 486 || getTag() == 566 || getTag() == 956)
    //{
    //opserr <<" \n\n element " << getTag() << endlnn;
    //opserr<<"\n tot node " << nd1Ptr->getTag() <<" x "<< TotDis1(0) <<" y "<< TotDis1(1) << " z "<< TotDis1(2);
    //
    //opserr<<"\n incr node " << nd1Ptr->getTag() <<" x "<< incrdelDis1(0) <<" y "<< incrdelDis1(1) << " z "<< incrdelDis1(2);
    //
    //opserr << "\nincr del node " << nd1Ptr->getTag() << " x " << IncrDis1(0) <<" y "<< IncrDis1(1) << " z "<< IncrDis1(2) << endlnn;
    //
    //opserr << " incr del node " << nd2Ptr->getTag() << " x " << IncrDis2(0) <<" y "<< IncrDis2(1) << " z "<< IncrDis2(2) << endlnn;
    //
    //opserr << " incr del node " << nd3Ptr->getTag() << " x " << IncrDis3(0) <<" y "<< IncrDis3(1) << " z "<< IncrDis3(2) << endlnn;
    //
    //opserr << " incr del node " << nd4Ptr->getTag() << " x " << IncrDis4(0) <<" y "<< IncrDis4(1) << " z "<< IncrDis4(2) << endlnn;
    //
    //opserr << " incr del node " << nd5Ptr->getTag() << " x " << IncrDis5(0) <<" y "<< IncrDis5(1) << " z "<< IncrDis5(2) << endlnn;
    //opserr << " incr del node " << nd6Ptr->getTag() << " x " << IncrDis6(0) <<" y "<< IncrDis6(1) << " z "<< IncrDis6(2) << endlnn;
    //opserr << " incr del node " << nd7Ptr->getTag() << " x " << IncrDis7(0) <<" y "<< IncrDis7(1) << " z "<< IncrDis7(2) << endlnn;
    //opserr << " incr del node " << nd8Ptr->getTag() << " x " << IncrDis8(0) <<" y "<< IncrDis8(1) << " z "<< IncrDis8(2) << endlnn;
    //}

    increment_dispDu.val(1,1)=IncrDis1(0); increment_dispDu.val(1,2)=IncrDis1(1);increment_dispDu.val(1,3)=IncrDis1(2);
    increment_dispDu.val(2,1)=IncrDis2(0); increment_dispDu.val(2,2)=IncrDis2(1);increment_dispDu.val(2,3)=IncrDis2(2);
    increment_dispDu.val(3,1)=IncrDis3(0); increment_dispDu.val(3,2)=IncrDis3(1);increment_dispDu.val(3,3)=IncrDis3(2);
    increment_dispDu.val(4,1)=IncrDis4(0); increment_dispDu.val(4,2)=IncrDis4(1);increment_dispDu.val(4,3)=IncrDis4(2);
    increment_dispDu.val(5,1)=IncrDis5(0); increment_dispDu.val(5,2)=IncrDis5(1);increment_dispDu.val(5,3)=IncrDis5(2);
    increment_dispDu.val(6,1)=IncrDis6(0); increment_dispDu.val(6,2)=IncrDis6(1);increment_dispDu.val(6,3)=IncrDis6(2);
    increment_dispDu.val(7,1)=IncrDis7(0); increment_dispDu.val(7,2)=IncrDis7(1);increment_dispDu.val(7,3)=IncrDis7(2);
    increment_dispDu.val(8,1)=IncrDis8(0); increment_dispDu.val(8,2)=IncrDis8(1);increment_dispDu.val(8,3)=IncrDis8(2);
    
    return increment_dispDu;
  }

////#############################################################################
//=========================================================================
// Incremental displacement of pressure part-- p        05/20/2002
//=========================================================================
tensor EightNodeBrick_u_p_U::incr_dispDp()
  {
    const int dimensions[] = {8};  
    tensor increment_dispDp(1, dimensions, 0.0);

     
    const Vector &IncrDis1 = nd1Ptr->getIncrDeltaDisp();
    const Vector &IncrDis2 = nd2Ptr->getIncrDeltaDisp();
    const Vector &IncrDis3 = nd3Ptr->getIncrDeltaDisp();
    const Vector &IncrDis4 = nd4Ptr->getIncrDeltaDisp();

    const Vector &IncrDis5 = nd5Ptr->getIncrDeltaDisp();
    const Vector &IncrDis6 = nd6Ptr->getIncrDeltaDisp();
    const Vector &IncrDis7 = nd7Ptr->getIncrDeltaDisp();
    const Vector &IncrDis8 = nd8Ptr->getIncrDeltaDisp();

 // Get the fourth incremental displacement for pressure part. Xiaoyan 05/20/2002
    increment_dispDp.val(1)=IncrDis1(3);
    increment_dispDp.val(2)=IncrDis2(3);
    increment_dispDp.val(3)=IncrDis3(3);
    increment_dispDp.val(4)=IncrDis4(3);
    increment_dispDp.val(5)=IncrDis5(3);
    increment_dispDp.val(6)=IncrDis6(3);
    increment_dispDp.val(7)=IncrDis7(3);
    increment_dispDp.val(8)=IncrDis8(3);

    return increment_dispDp;
  }
////#############################################################################
//=========================================================================
// Incremental displacement of fluid part-- U     02/10/2002
//=========================================================================
tensor EightNodeBrick_u_p_U::incr_dispDU()
  {
    const int dimensions[] = {8,3};  
    tensor increment_dispDU(2, dimensions, 0.0);

     
    const Vector &IncrDis1 = nd1Ptr->getIncrDeltaDisp();
    const Vector &IncrDis2 = nd2Ptr->getIncrDeltaDisp();
    const Vector &IncrDis3 = nd3Ptr->getIncrDeltaDisp();
    const Vector &IncrDis4 = nd4Ptr->getIncrDeltaDisp();

    const Vector &IncrDis5 = nd5Ptr->getIncrDeltaDisp();
    const Vector &IncrDis6 = nd6Ptr->getIncrDeltaDisp();
    const Vector &IncrDis7 = nd7Ptr->getIncrDeltaDisp();
    const Vector &IncrDis8 = nd8Ptr->getIncrDeltaDisp();

 // Get the Last three incremental displacement for fluid part. Xiaoyan 02/08/2002
 
    increment_dispDU.val(1,1)=IncrDis1(4);increment_dispDU.val(1,2)=IncrDis1(5); 
    increment_dispDU.val(1,3)=IncrDis1(6);

    increment_dispDU.val(2,1)=IncrDis2(4);increment_dispDU.val(2,2)=IncrDis2(5); 
    increment_dispDU.val(2,3)=IncrDis2(6);

    increment_dispDU.val(3,1)=IncrDis3(4);increment_dispDU.val(3,2)=IncrDis3(5); 
    increment_dispDU.val(3,3)=IncrDis3(6);

    increment_dispDU.val(4,1)=IncrDis4(4);increment_dispDU.val(4,2)=IncrDis4(5); 
    increment_dispDU.val(4,3)=IncrDis4(6);
 
    increment_dispDU.val(5,1)=IncrDis5(4);increment_dispDU.val(5,2)=IncrDis5(5); 
    increment_dispDU.val(5,3)=IncrDis5(6);
 
    increment_dispDU.val(6,1)=IncrDis6(4);increment_dispDU.val(4,2)=IncrDis6(5); 
    increment_dispDU.val(6,3)=IncrDis6(6);
 
    increment_dispDU.val(7,1)=IncrDis7(4);increment_dispDU.val(7,2)=IncrDis7(5); 
    increment_dispDU.val(7,3)=IncrDis7(6);
 
    increment_dispDU.val(8,1)=IncrDis8(4);increment_dispDU.val(8,2)=IncrDis8(5); 
    increment_dispDU.val(8,3)=IncrDis8(6);

    return increment_dispDU;
  }
//=========================================================================
// Total displacement of solid part --u    03/22/2002
//=========================================================================
tensor EightNodeBrick_u_p_U::total_dispDu()
  {
    const int dimensions[] = {8,3};  
    tensor total_dispDu(2, dimensions, 0.0);
      
    // Get first three total displacement for solid part. Xiaoyan 03/22/2002
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
 
    total_dispDu.val(1,1)=TotDis1(0); total_dispDu.val(1,2)=TotDis1(1);total_dispDu.val(1,3)=TotDis1(2);
    total_dispDu.val(2,1)=TotDis2(0); total_dispDu.val(2,2)=TotDis2(1);total_dispDu.val(2,3)=TotDis2(2);
    total_dispDu.val(3,1)=TotDis3(0); total_dispDu.val(3,2)=TotDis3(1);total_dispDu.val(3,3)=TotDis3(2);
    total_dispDu.val(4,1)=TotDis4(0); total_dispDu.val(4,2)=TotDis4(1);total_dispDu.val(4,3)=TotDis4(2);
    total_dispDu.val(5,1)=TotDis5(0); total_dispDu.val(5,2)=TotDis5(1);total_dispDu.val(5,3)=TotDis5(2);
    total_dispDu.val(6,1)=TotDis6(0); total_dispDu.val(6,2)=TotDis6(1);total_dispDu.val(6,3)=TotDis6(2);
    total_dispDu.val(7,1)=TotDis7(0); total_dispDu.val(7,2)=TotDis7(1);total_dispDu.val(7,3)=TotDis7(2);
    total_dispDu.val(8,1)=TotDis8(0); total_dispDu.val(8,2)=TotDis8(1);total_dispDu.val(8,3)=TotDis8(2);

    return total_dispDu;
  }

////#############################################################################
//=========================================================================
// Total displacement of fluid part--U     03/22/2002
//=========================================================================
tensor EightNodeBrick_u_p_U::total_dispDU()
  {
    const int dimensions[] = {8,3};  
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


 
    return total_dispDU;
  }


////#############################################################################



////#############################################################################
tensor EightNodeBrick_u_p_U::total_dispDu(FILE *fp, double * u)
  {
    const int dimensions[] = {8,3};  
    tensor total_dispDu(2, dimensions, 0.0);
    //    double totalx, totaly, totalz;
    //    totalx=0;
    //    totaly=0;
    //    totalz=0;

    //for ( int i=0 ; i<8 ; i++ )  // Xiaoyan changed from 20 to 8 for 8 nodes
    //
    //  {
    //    // total_disp.val(i+1,1) = nodes[ G_N_numbs[i] ].total_translation_x(u);
    //    // total_disp.val(i+1,2) = nodes[ G_N_numbs[i] ].total_translation_y(u);
    //    // total_disp.val(i+1,3) = nodes[ G_N_numbs[i] ].total_translation_z(u);
    //    // Xiaoyan changed to the following 09/27/00
    //    Vector TotalTranDis = nodes[ G_N_numbs[i] ].getDisp();
    //
    //    total_disp.val(i+1,1) = TotalTranDis(0);
    //	total_disp.val(i+1,2) = TotalTranDis(1);
    //    total_disp.val(i+1,3) = TotalTranDis(2);
    //
    //  }
      
    //Zhaohui using node pointers, which come from the Domain
    const Vector &TotDis1 = nd1Ptr->getTrialDisp();
    const Vector &TotDis2 = nd2Ptr->getTrialDisp();
    const Vector &TotDis3 = nd3Ptr->getTrialDisp();
    const Vector &TotDis4 = nd4Ptr->getTrialDisp();
    const Vector &TotDis5 = nd5Ptr->getTrialDisp();
    const Vector &TotDis6 = nd6Ptr->getTrialDisp();
    const Vector &TotDis7 = nd7Ptr->getTrialDisp();
    const Vector &TotDis8 = nd8Ptr->getTrialDisp();
    
    total_dispDu.val(1,1)=TotDis1(0); total_dispDu.val(1,2)=TotDis1(1);total_dispDu.val(1,3)=TotDis1(2);
    total_dispDu.val(2,1)=TotDis2(0); total_dispDu.val(2,2)=TotDis2(1);total_dispDu.val(2,3)=TotDis2(2);
    total_dispDu.val(3,1)=TotDis3(0); total_dispDu.val(3,2)=TotDis3(1);total_dispDu.val(3,3)=TotDis3(2);
    total_dispDu.val(4,1)=TotDis4(0); total_dispDu.val(4,2)=TotDis4(1);total_dispDu.val(4,3)=TotDis4(2);
    total_dispDu.val(5,1)=TotDis5(0); total_dispDu.val(5,2)=TotDis5(1);total_dispDu.val(5,3)=TotDis5(2);
    total_dispDu.val(6,1)=TotDis6(0); total_dispDu.val(6,2)=TotDis6(1);total_dispDu.val(6,3)=TotDis6(2);
    total_dispDu.val(7,1)=TotDis7(0); total_dispDu.val(7,2)=TotDis7(1);total_dispDu.val(7,3)=TotDis7(2);
    total_dispDu.val(8,1)=TotDis8(0); total_dispDu.val(8,2)=TotDis8(1);total_dispDu.val(8,3)=TotDis8(2);

    return total_dispDu;
  }
//====================================================================
////#############################################################################

tensor EightNodeBrick_u_p_U::total_dispDU(FILE *fp, double * u)
  {
    const int dimensions[] = {8,3};  
    tensor total_dispDU(2, dimensions, 0.0);
  
    const Vector &TotDis1  = nd1Ptr->getTrialDisp();
    const Vector &TotDis2  = nd2Ptr->getTrialDisp();
    const Vector &TotDis3  = nd3Ptr->getTrialDisp();
    const Vector &TotDis4  = nd4Ptr->getTrialDisp();
    const Vector &TotDis5  = nd5Ptr->getTrialDisp();
    const Vector &TotDis6  = nd6Ptr->getTrialDisp();
    const Vector &TotDis7  = nd7Ptr->getTrialDisp();
    const Vector &TotDis8  = nd8Ptr->getTrialDisp();

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

    return total_dispDU;
}
//====================================================================


void EightNodeBrick_u_p_U::incremental_UpdateDu()
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

    int dh_dim[] = {8,3};   //Xiaoyan changed from {20,3} to {8,3}  07/12/00
    tensor dh(2, dh_dim, 0.0);

//    tensor Constitutive( 4, def_dim_4, 0.0);

    //    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {8,3}; //Xiaoyan changed from {20,3} to {8,3}  07/12/00
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
                //::::::   ::printf(" void EightNodeBrick_u_p_U::incremental_Update()\n");
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
		  opserr << "EightNodeBrick_u_p_U::incremental_Update (tag: " << this->getTag() << "), not converged\n";
		//matpoint[where].setEPS( mmodel->getEPS() );
		//}
		
		//else if ( tmp_ndm ) 
		//  (matpoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
		//else {
               	//   g3ErrorHandler->fatal("EightNodeBrick_u_p_U::incremental_Update (tag: %d), no strain or stress state vars", this->getTag());
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

void EightNodeBrick_u_p_U::incremental_UpdateDU()
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

    int dh_dim[] = {8,3};   
    tensor dh(2, dh_dim, 0.0);

//    tensor Constitutive( 4, def_dim_4, 0.0);

    //    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {8,3}; 
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
                //::::::   ::printf(" void EightNodeBrick_u_p_U::incremental_Update()\n");
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
		  opserr << "EightNodeBrick_u_p_U::incremental_Update (tag: " << this->getTag() << "), not converged\n";

		//matpoint[where].setEPS( mmodel->getEPS() );
		//}
		
		//else if ( tmp_ndm ) 
		//  (matpoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
		//else {
               	//   g3ErrorHandler->fatal("EightNodeBrick_u_p_U::incremental_Update (tag: " << this->getTag() << "), no strain or stress state vars", this->getTag());
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

//#############################################################################
void EightNodeBrick_u_p_U::set_strain_stress_tensorDu(FILE *fp, double * u)
  {
    int dh_dim[] = {8,3};   // Xiaoyan changed from {20,3} to {8,3}
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


    static int disp_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}
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
void EightNodeBrick_u_p_U::set_strain_stress_tensorDU(FILE *fp, double * u)
  {
    int dh_dim[] = {8,3};   
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


    static int disp_dim[] = {8,3};    
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
EightNodeBrick_u_p_U & EightNodeBrick_u_p_U::operator[](int subscript)
  {
    return ( *(this+subscript) );
  }

//Finite_Element & EightNodeBrick_u_p_U::operator[](short subscript)
//  {
//    return ( *(this+subscript) );
//  }

//Finite_Element & EightNodeBrick_u_p_U::operator[](unsigned subscript)
//  {
//    return ( *(this+subscript) );
//  }


////#############################################################################
int EightNodeBrick_u_p_U::commitState ()
{
    // int order = theQuadRule->getOrder();     // Commented by Xiaoyan

    int i;
    //int j, k;      // Xiaoyan added k for three dimension		       
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "EightNodeBrick-u_p_U::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
//    int count  = r_integration_order* s_integration_order * t_integration_order;    commented by Xiaoyan 
// 09/24/2001.  This is a unused variable

    //for (i = 0; i < r_integration_order; i++)		    // Xiaoyan chaneged order to
    //  for (j = 0; j < s_integration_order; j++)	    // r_integration_order,
    //							    // s_integration_order, and
    //	    for (k = 0; k < t_integration_order; k++)	    // added t_integration_order,
    //         retVal += (GaussPtheMaterial[i][j][k]).commitState();  

    Vector pp = getResistingForce();

    //if ( this->getTag() == 1 || this->getTag() == 700)
    //{
      //for (i = 0; i < count; i++)
      for (i = 0; i < 8; i++)
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
	 
//	 double  p = -1*( prin.cval(1, 1)+ prin.cval(2, 2) +prin.cval(3, 3) )/3.0;	        commented by Xiaoyan
//	 double  ev = -1*( stnprin.cval(1, 1)+ stnprin.cval(2, 2) + stnprin.cval(3, 3) )/3.0;   commented by Xiaoyan 
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

//out22Jan2001	 if (strcmp(matpoint[i]->matmodel->getType(),"Template3Dep") == 0)
//out22Jan2001          {
//out22Jan2001       	  st = ( ((Template3Dep *)(matpoint[i]->matmodel))->getEPS())->getStress();
//out22Jan2001       	  prin = st.principal();
//out22Jan2001 	 }
//out22Jan2001 	 else
//out22Jan2001 	 {
//out22Jan2001        	  st = matpoint[i]->getStressTensor();
//out22Jan2001       	  prin = st.principal();
//out22Jan2001 	 
//out22Jan2001 	 }
      
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
int EightNodeBrick_u_p_U::revertToLastCommit ()
{
  //  int order = theQuadRule->getOrder();	// Commented by Xiaoyan
    int i;
    //int j, k;     // Xiaoyan added k for three dimension	
    int retVal = 0;

    // Loop over the integration points and revert to last committed material states
    int count  = r_integration_order* s_integration_order * t_integration_order;
    //for (i = 0; i < r_integration_order; i++)		   // Xiaoyan chaneged order to 
    //	for (j = 0; j < s_integration_order; j++)	   // r_integration_order,      
    //	    for (k = 0; k < t_integration_order; k++)	   // s_integration_order, and  
		      					   // added t_integration_order,
	    //retVal += (theMaterial[i][j][k]).revertToLastCommit();

    for (i = 0; i < count; i++)	    
       retVal += matpoint[i]->revertToLastCommit(); 
   

    return retVal;
}

////#############################################################################
//=============================================================================
int EightNodeBrick_u_p_U::revertToStart () 
{
    int i;     // Xiaoyan added k for three dimension	
    int retVal = 0;

    // Loop over the integration points and revert to last committed material states
    //for (i = 0; i < r_integration_order; i++)		   // Xiaoyan chaneged order to 
    //	for (j = 0; j < s_integration_order; j++)	   // r_integration_order,      
    //	    for (k = 0; k < t_integration_order; k++)	   // s_integration_order, and  
							   // added t_integration_order,
    //	    retVal += (theMaterial[i][j][k]).revertToLastCommit();

    int count  = r_integration_order* s_integration_order * t_integration_order;
    	     
    for (i = 0; i < count; i++)
       retVal += matpoint[i]->revertToStart(); 
    
    
    return retVal;

    // Loop over the integration points and revert to initial material states
   } 
////#############################################################################

//=============================================================================
//  The following are come from FourNodeQuad.cc	 Xiaoyan 07/06/00
//  The following are come from FourNodeQuad.cc	 Xiaoyan 07/06/00
//  The following are come from FourNodeQuad.cc	 Xiaoyan 07/06/00
//=============================================================================

//=============================================================================
int EightNodeBrick_u_p_U::getNumExternalNodes () const
{
    return 8;  //changed from 4 to 8 Xiaoyan 07/06/00
}


//=============================================================================
const ID& EightNodeBrick_u_p_U::getExternalNodes ()
{
    return connectedExternalNodes;
}

Node **
EightNodeBrick_u_p_U::getNodePtrs (void)
{
  theNodes[0] = nd1Ptr;
  theNodes[1] = nd2Ptr;
  theNodes[2] = nd3Ptr;
  theNodes[3] = nd4Ptr;
  theNodes[4] = nd5Ptr;
  theNodes[5] = nd6Ptr;
  theNodes[6] = nd7Ptr;
  theNodes[7] = nd8Ptr;

  return theNodes;
}


//=============================================================================
int EightNodeBrick_u_p_U::getNumDOF ()
{
    return 7*8;	     // one node has 7 DOF. Xiaoyan 01/13/2002
}

//=============================================================================
void EightNodeBrick_u_p_U::setDomain (Domain *theDomain)
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
      
      nd1Ptr = theDomain->getNode(Nd1);
      nd2Ptr = theDomain->getNode(Nd2);
      nd3Ptr = theDomain->getNode(Nd3);
      nd4Ptr = theDomain->getNode(Nd4);
              
      //Xiaoyan added 5-8  07/06/00
      nd5Ptr = theDomain->getNode(Nd5);
      nd6Ptr = theDomain->getNode(Nd6);
      nd7Ptr = theDomain->getNode(Nd7);
      nd8Ptr = theDomain->getNode(Nd8);
      
      if (nd1Ptr == 0 || nd2Ptr == 0 || nd3Ptr == 0 || nd4Ptr == 0||
          nd5Ptr == 0 || nd6Ptr == 0 || nd7Ptr == 0 || nd8Ptr == 0 ) { 
      	//Xiaoyan added 5-8  07/06/00
      
      	opserr << "FATAL ERROR EightNodeBrick_u_p_U (tag: " << this->getTag() << "), node not found in domain\n";
	exit(-1);
      }
      
      int dofNd1 = nd1Ptr->getNumberDOF();
      int dofNd2 = nd2Ptr->getNumberDOF();
      int dofNd3 = nd3Ptr->getNumberDOF();
      int dofNd4 = nd4Ptr->getNumberDOF();
      
      //Xiaoyan added 5-8  07/06/00
      int dofNd5 = nd5Ptr->getNumberDOF();
      int dofNd6 = nd6Ptr->getNumberDOF();
      int dofNd7 = nd7Ptr->getNumberDOF();
      int dofNd8 = nd8Ptr->getNumberDOF();
      								      
      if (dofNd1 != 7 || dofNd2 != 7 || dofNd3 != 7 || dofNd4 != 7 ||  // Changed 3 to 7 Xiaoyan
          dofNd5 != 7 || dofNd6 != 7 || dofNd7 != 7 || dofNd8 != 7 ) {
      	opserr << "FATAL ERROR EightNodeBrick_u_p_U (tag: " << this->getTag() << "), has differing number of DOFs at its nodes\n";
	exit(-1);
      }
      this->DomainComponent::setDomain(theDomain);
    }   
}
//=============================================================================
//=============================================================================
const Matrix &EightNodeBrick_u_p_U::getInitialStiff () 
{
  if (Ki == 0)
    Ki = new Matrix(this->getTangentStiff());

  if (Ki == 0) {
    opserr << "FATAL EightNodeBrick_u_p_U::getInitialStiff() -";
    opserr << "ran out of memory\n";
    exit(-1);
  }  
    
  return *Ki;
}
//=============================================================================
void EightNodeBrick_u_p_U::zeroLoad()
{
     Q.Zero();
}

//=============================================================================
int 
EightNodeBrick_u_p_U::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
/*  opserr << "EightNodeBrick_u_p_U::addLoad - load type unknown for ele with tag: " << this->getTag() << "\n",
			  this->getTag());
  
  return -1;
*/
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  
  if (type == LOAD_TAG_BrickSelfWeight) {
  
    Vector bforce(56);  
    // Check for a quick return
    //opserr << "rho " << rho << endln;
    if (rho == 0.0)
    	return 0;
  
    Vector ba(56), bfx(3);  
    bfx(0) = bf(0) * loadFactor;
    bfx(1) = bf(1) * loadFactor;
    bfx(2) = bf(2) * loadFactor;
  
    ba( 0) = bfx(0);    
    ba( 1) = bfx(1);    
    ba( 2) = bfx(2);
    ba( 3)=  0.0;
    ba( 4) = bfx(0);    
    ba( 5) = bfx(1);    
    ba( 6) = bfx(2);    

    ba( 7) = bfx(0);    
    ba( 8) = bfx(1);    
    ba( 9) = bfx(2);    
    ba(10) = 0.0;    
    ba(11) = bfx(0);    
    ba(12) = bfx(1);    
    ba(13) = bfx(2); 
       
    ba(14) = bfx(0);    
    ba(15) = bfx(1);    
    ba(16) = bfx(2);    
    ba(17) = 0.0;      
    ba(18) = bfx(0);    
    ba(19) = bfx(1);    
    ba(20) = bfx(2);    

    ba(21) = bfx(0);    
    ba(22) = bfx(1);    
    ba(23) = bfx(2);    
    ba(24) = 0.0;      
    ba(25) = bfx(0);    
    ba(26) = bfx(1);    
    ba(27) = bfx(2);    

    ba(28) = bfx(0);    
    ba(29) = bfx(1);    
    ba(30) = bfx(2);    
    ba(31) = 0.0;      
    ba(32) = bfx(0);    
    ba(33) = bfx(1);    
    ba(34) = bfx(2);    

    ba(35) = bfx(0);    
    ba(36) = bfx(1);    
    ba(37) = bfx(2);    
    ba(38) = 0.0;      
    ba(39) = bfx(0);    
    ba(40) = bfx(1);    
    ba(41) = bfx(2);    

    ba(42) = bfx(0);    
    ba(43) = bfx(1);    
    ba(44) = bfx(2);    
    ba(45) = 0.0;      
    ba(46) = bfx(0);    
    ba(47) = bfx(1);    
    ba(48) = bfx(2);    

    ba(49) = bfx(0);    
    ba(50) = bfx(1);    
    ba(51) = bfx(2);    
    ba(52) = 0.0;      
    ba(53) = bfx(0);    
    ba(54) = bfx(1);    
    ba(55) = bfx(2);

    //Form equivalent body force
    this->getMass();
    bforce.addMatrixVector(0.0, M, ba, 1.0);
    Q.addVector(1.0, bforce, 1.0);       
  } else  {
    opserr << "EightNodeBrick_u_p_U::addLoad() - 8NodeBrick_u_p_U " << this->getTag() << ",load type " << type << " unknown\n";
    return -1;
  }
  return 0;
}


//=============================================================================
int EightNodeBrick_u_p_U::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0) 
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = nd1Ptr->getRV(accel);
	const Vector &Raccel2 = nd2Ptr->getRV(accel);
	const Vector &Raccel3 = nd3Ptr->getRV(accel);
	const Vector &Raccel4 = nd4Ptr->getRV(accel);
        // Xiaoyan added the following four 09/27/00	
	const Vector &Raccel5 = nd5Ptr->getRV(accel);
	const Vector &Raccel6 = nd6Ptr->getRV(accel);
	const Vector &Raccel7 = nd7Ptr->getRV(accel);
	const Vector &Raccel8 = nd8Ptr->getRV(accel);

    if (7 != Raccel1.Size() || 7 != Raccel2.Size() || 7 != Raccel3.Size() ||
       	7 != Raccel4.Size() || 7 != Raccel5.Size() || 7 != Raccel6.Size() || 
	7 != Raccel7.Size() || 7 != Raccel8.Size()  ){	
		opserr << "EightNodeBrick_u_p_U::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
		return -1;
    }

	static Vector ra(56);  // Changed form 24 to 56  Xiaoyan 03/22/2002

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

    // Want to add ( - fact * M R * accel ) to unbalance
    // Take advantage of lumped mass matrix
    // Mass matrix is computed in setDomain()
    
    //double column_mass = 0;
    //for (int i = 0; i < 24; i++)   
    //   column_mass += M(1,i);
    //column_mass = column_mass/3.0;

    //opserr << " addInerti... column_mass " << column_mass << endln;

//    for (int i = 0; i < 56; i++)   
//	    	Q(i) += -M(i,i)*ra(i);
    Q.addMatrixVector(1.0,M,ra,-1.0);
    return 0;
}

//=============================================================================
const Vector EightNodeBrick_u_p_U::FormEquiBodyForce()
{
    Vector bforce(56);  

    // Check for a quick return
    //opserr << "rho " << rho << endln;
    if (rho == 0.0) 
    	return bforce;

    Vector ba(56);  

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

    //Form equivalent body force
    bforce.addMatrixVector(0.0, M, ba, 1.0);
    //opserr << " ba " << ba;
    //opserr << " M " << M;
    //if (getTag() == 886)
    //opserr << " @@@@@ FormEquiBodyForce  " << bforce;
    
    return bforce;
}

//=============================================================================
// Setting initial E according to the initial pressure p
//void EightNodeBrick_u_p_U::setInitE(void)
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
//    		   nd2Crds( dir-1)+
//    		   nd3Crds( dir-1)+
//    		   nd4Crds( dir-1)+
//    		   nd5Crds( dir-1)+
//    		   nd6Crds( dir-1)+
//    		   nd7Crds( dir-1)+
//    		   nd8Crds( dir-1);
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
const Vector &EightNodeBrick_u_p_U::getResistingForce () 
{     
    int force_dim[] = {8,7}; 
    tensor nodalforces(2,force_dim,0.0);
     
    nodalforces = nodal_forces();

    //converting nodalforce tensor to vector
    for (int i = 0; i< 8; i++)
      for (int j = 0; j < 7; j++)
	p(i *7 + j) = nodalforces.cval(i+1, j+1);

    //opserr << "p" << p;
    //opserr << "Q" << Q;

//    p = p - Q;
    p.addVector(1.0, Q, -1.0);	     // 06/03/2002

    //opserr << "p-Q" << p;
    return p;
}

//=============================================================================
const Vector &EightNodeBrick_u_p_U::getResistingForceIncInertia () 
{
	// Check for a quick return
	if (rho == 0.0)
		return this->getResistingForce();

	//opserr << "Node555 trialDisp " << nd1Ptr->getTrialDisp();

	const Vector &accel1 = nd1Ptr->getTrialAccel();
        //opserr << "\nnode accel " << nd1Ptr->getTag() << " x " << accel1(0) <<" y "<< accel1(1) << " z "<< accel1(2) << endln;

	const Vector &accel2 = nd2Ptr->getTrialAccel();
        //opserr << "node accel " << nd2Ptr->getTag() << " x " << accel2(0) <<" y "<< accel2(1) << " z "<< accel2(2) << endln;
							               		    
	const Vector &accel3 = nd3Ptr->getTrialAccel();	               		    
        //opserr << "node accel " << nd3Ptr->getTag() << " x " << accel3(0) <<" y "<< accel3(1) << " z "<< accel3(2) << endln;
							               		    
	const Vector &accel4 = nd4Ptr->getTrialAccel();	               		    
        //opserr << "node accel " << nd4Ptr->getTag() << " x " << accel4(0) <<" y "<< accel4(1) << " z "<< accel4(2) << endln;
							               		    
        // Xiaoyan added the following four 09/27/00
	const Vector &accel5 = nd5Ptr->getTrialAccel();	               		    
        //opserr << "node accel " << nd5Ptr->getTag() << " x " << accel5(0) <<" y "<< accel5(1) << " z "<< accel5(2) << endln;
							               		    
	const Vector &accel6 = nd6Ptr->getTrialAccel();	               		    
        //opserr << "node accel " << nd6Ptr->getTag() << " x " << accel6(0) <<" y "<< accel6(1) << " z "<< accel6(2) << endln;
							               		    
	const Vector &accel7 = nd7Ptr->getTrialAccel();	               		    
        //opserr << "node accel " << nd7Ptr->getTag() << " x " << accel7(0) <<" y "<< accel7(1) << " z "<< accel7(2) << endln;
							               		    
	const Vector &accel8 = nd8Ptr->getTrialAccel();	               		    
        //opserr << "node accel " << nd8Ptr->getTag() << " x " << accel8(0) <<" y "<< accel8(1) << " z "<< accel8(2) << endln;
	
	static Vector a(56);  // originally 8

// Modified to 3-solid 1-p 3-fluid        Xioayan 03/14/2002
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
   
	// Compute the current resisting force
	this->getResistingForce();

	// Take advantage of lumped mass matrix
	// Mass matrix is computed in setDomain()
	//opserr << " M_ii \n";

        //double column_mass = 0;
        //for (int i = 0; i < 24; i++)   
        //   column_mass += M(1,i);
        //column_mass = column_mass/3.0;

/*	for (int i = 0; i < 56; i++)
	{   
	   p(i) += M(i,i)*a(i);
	   //opserr << " " << M(i, i);
	}
	//opserr << endln;
	//opserr << "P+=Ma" << P<< endln;
*/
	p.addMatrixVector(1.0, M, a,1.0);
	return p;
} 


////#############################################################################
int EightNodeBrick_u_p_U::get_global_number_of_node(int local_node_number)
{
  //return G_N_numbs[local_node_number];	
  return connectedExternalNodes(local_node_number);
}

////#############################################################################
int  EightNodeBrick_u_p_U::get_Brick_Number()
{
  //return elem_numb;
  return this->getTag();
}

////#############################################################################
//=============================================================================
double EightNodeBrick_u_p_U::get_Gauss_p_c(short order, short point_numb)
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

double EightNodeBrick_u_p_U::get_Gauss_p_w(short order, short point_numb)
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
int * EightNodeBrick_u_p_U::get_LM()
  {
    return LM;
  }

/////#############################################################################
//void EightNodeBrick_u_p_U::set_LM(Node * node)
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
// returns nodal forces for given stress field in an element
tensor EightNodeBrick_u_p_U::nodal_forcesFu()
  {
    int force_dim[] = {8,3};  

    tensor nodal_forcesFu(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};  

    tensor dh(2, dh_dim, 0.0);

    stresstensor stress_at_GP(0.0);

    double det_of_Jacobian = 0.0;

    straintensor incremental_strain;

    static int disp_dim[] = {8,3};   // Xiaoyan changed from {20,3} to {8,3}
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
//		if (where == 0)
//   		//opserr << " In nodal_force delta_incremental_strain tag "<< getTag() <<"  " <<incremental_strain << endln;
////		opserr << " el tag = "<< getTag();
//		
		int err = ( matpoint[where]->matmodel )->setTrialStrainIncr( incremental_strain);
		if ( err)
		  opserr << "EightNodeBrick_u_p_U::getStiffnessTensor (tag: " << this->getTag() << "), not converged\n";



		//char *test = matpoint[where]->matmodel->getType();
		// fmk - changing if so if into else block must be Template3Dep
//		if (strcmp(matpoint[where]->matmodel->getType(),"Template3Dep") != 0)
		   stress_at_GP = matpoint[where]->getStressTensor();


                //stress_at_GP.report("\n stress_at_GPtensor at GAUSS point for nodal forces \n");

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forcesFu =
                     nodal_forcesFu + dhGlobal("ib")*stress_at_GP("ab")*weight;
                //nodal_forces.print("nf","\n\n Nodal Forces \n");
 
              }
          }
      }

    //opserr << "\n element no. " << getTag() << endln;
    //nodal_forces.print("nf","\n Nodal Forces \n");
    return nodal_forcesFu;

  }
////#############################################################################
// returns nodal forces of fluid part in an element
tensor EightNodeBrick_u_p_U::nodal_forcesFU()
  {
 // This is a function for a generelized force for the fluid component...  
 // For fluid there is only pressure not shear force, so the force for 
 // fluid will be p*delta_{ij}      xiaoyan 05/20/2002
    int force_dim[] = {8,3};  

    tensor nodal_forcesFU(2,force_dim,0.0);
    tensor I2("I", 2, def_dim_2);	   // delta_{ij} added 05/20/2002

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};  
    int h_dim[] = {8};     // added for shape function h  05/21/2002

    tensor dh(2, dh_dim, 0.0);
    tensor h(1, h_dim, 0.0);

    stresstensor stress_at_GP(0.0);
    tensor p_at_GP;		  // added 05/21/2002
    double pp;			  // added 05/21/2002

    double det_of_Jacobian = 0.0;

//    straintensor incremental_strain;

    static int disp_dim[] = {8};     // changed from {8,1}
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
                
//		incremental_strain =
//                     (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
//		if (where == 0)
//   		//opserr << " In nodal_force delta_incremental_strain tag "<< getTag() <<"  " <<incremental_strain << endln;
////		opserr << " el tag = "<< getTag();
//		
//		int err = ( matpoint[where]->matmodel )->setTrialStrainIncr( incremental_strain);
//		if ( err)
//               	   opserr << "EightNodeBrick_u_p_U::getStiffnessTensor (tag: " << this->getTag() << "), not converged",
//		    		 this->getTag());



		//char *test = matpoint[where]->matmodel->getType();
		// fmk - changing if so if into else block must be Template3Dep
//		if (strcmp(matpoint[where]->matmodel->getType(),"Template3Dep") != 0)
		   stress_at_GP = matpoint[where]->getStressTensor();
		   p_at_GP=h("i")*incremental_p("i");
		   pp=p_at_GP.val(1);

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forcesFU =
                     nodal_forcesFU + dhGlobal("ib")*I2("ab")*pp*weight;
                //nodal_forces.print("nf","\n\n Nodal Forces \n");
 
              }
          }
      }

    //opserr << "\n element no. " << getTag() << endln;
    //nodal_forces.print("nf","\n Nodal Forces \n");
    return nodal_forcesFU;

  }
////#############################################################################
////#############################################################################
// returns nodal forces for given stress field in an element
tensor EightNodeBrick_u_p_U::nodal_forces()
  {
    int force_dim[] = {8,7};  

    tensor nodalForceFu = nodal_forcesFu();
    tensor nodalForceFU = nodal_forcesFU();

    tensor nodal_FORCES(2,force_dim,0.0);

    for( int i=1; i<=8; i++)
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
//     nodal_FORCES.print("nf","\n Nodal Forces \n");

     return nodal_FORCES;
  }


////#############################################################################
////#############################################################################
// returns nodal forces for given ITERATIVE stress field in an element
// Solid part  Xiaoyan 02/08/2002
tensor EightNodeBrick_u_p_U::iterative_nodal_forcesFu()
  {
    int force_dim[] = {8,3}; 

    tensor nodal_forcesFu(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};   

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
                //.....::printf("EightNodeBrick_u_p_U::iterative_nodal_forces()  ----**************** where = %d \n", where);
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
                //nodal_forces.print("nf","\n EightNodeBrick_u_p_U::iterative_nodal_forces Nodal Forces ~~~~\n");

              }
          }
      }


    return nodal_forcesFu;

  }
////#############################################################################
// returns nodal forces for given ITERATIVE stress field in an element
// Fluid part  Xiaoyan 02/08/2002
tensor EightNodeBrick_u_p_U::iterative_nodal_forcesFU()
  {
    int force_dim[] = {8,3}; 

    tensor nodal_forcesFU(2,force_dim,0.0);
    tensor I2("I", 2, def_dim_2);	   // delta_{ij} added 05/20/2002

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};   
    int h_dim[] = {8};     // added for shape function h  05/21/2002

    tensor dh(2, dh_dim, 0.0);
    tensor h(1, h_dim, 0.0);

    stresstensor stress_at_GP(0.0);
    tensor p_at_GP;		  // added 05/21/2002
    double pp;			  // added 05/21/2002

    double det_of_Jacobian = 0.0;
    static int disp_dim[] = {8};   	// changed from {20,1}  05/20/2002
    tensor incremental_p(1,disp_dim,0.0); // \Delta U

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
                //.....::printf("EightNodeBrick_u_p_U::iterative_nodal_forces()  ----**************** where = %d \n", where);
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
                //nodal_forces.print("nf","\n EightNodeBrick_u_p_U::iterative_nodal_forces Nodal Forces ~~~~\n");

              }
          }
      }


    return nodal_forcesFU;

  }

////#############################################################################
// returns nodal forces for given constant stress field in the element
tensor EightNodeBrick_u_p_U::nodal_forces_from_stress(stresstensor & stress)
  {
    int force_dim[] = {8,3};  // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor nodal_forces(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    double weight = 0.0;

    int dh_dim[] = {8,3}; // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

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
                //.....::printf("EightNodeBrick_u_p_U::iterative_nodal_forces()  ----**************** where = %d \n", where);
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
                //nodal_forces.print("nf","\n EightNodeBrick_u_p_U::iterative_nodal_forces Nodal Forces ~~~~\n");

              }
          }
      }

    return nodal_forces;

  }


////#############################################################################
// returns nodal forces for given incremental strain field in an element
// by using the linearized constitutive tensor from the begining of the step !
tensor EightNodeBrick_u_p_U::linearized_nodal_forces()
  {
    int force_dim[] = {8,3};  // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor linearized_nodal_forces(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};  // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor dh(2, dh_dim, 0.0);

    tensor Constitutive( 4, def_dim_4, 0.0);

    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {8,3};  // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

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
               	   opserr << "EightNodeBrick_u_p_U::linearized_nodal_forces (tag: %d), not converged",
			     
		Constitutive = (matpoint[where]->matmodel)->getTangentTensor();
      		//    matpoint[where].setEPS( mmodel->getEPS() ); //Set the new EPState back
		//}
		//else if ( tmp_ndm ) { //Elastic case
		//    (matpoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
		//    Constitutive = (matpoint[where].p_matmodel)->getTangentTensor();
		//}
		//else {
               	//   g3ErrorHandler->fatal("EightNodeBrick_u_p_U::incremental_Update (tag: " << this->getTag() << "), could not getTangentTensor", this->getTag());
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
//double EightNodeBrick_u_p_U::get_first_etacone()
//  {
//    double ret = matpoint[0].etacone();
//
//    return ret;
//
//  }
//

//=============================================================================
int EightNodeBrick_u_p_U::sendSelf (int commitTag, Channel &theChannel) 
{ 
     // Not implemtented yet
     return 0;
}

//=============================================================================
int EightNodeBrick_u_p_U::recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) 
{   
     // Not implemtented yet
     return 0;
}


//=============================================================================
int EightNodeBrick_u_p_U::displaySelf (Renderer &theViewer, int displayMode, float fact)
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

    const Vector &end1Disp = nd1Ptr->getDisp();
    const Vector &end2Disp = nd2Ptr->getDisp();
    const Vector &end3Disp = nd3Ptr->getDisp();
    const Vector &end4Disp = nd4Ptr->getDisp();
    const Vector &end5Disp = nd5Ptr->getDisp();
    const Vector &end6Disp = nd6Ptr->getDisp();
    const Vector &end7Disp = nd7Ptr->getDisp();
    const Vector &end8Disp = nd8Ptr->getDisp();

    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    static Vector v5(3);
    static Vector v6(3);
    static Vector v7(3);
    static Vector v8(3);

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
void EightNodeBrick_u_p_U::Print(OPS_Stream &s, int flag)
{
    //report(" EightNodeBrick_u_p_U ");
    s << "EightNodeBrick_u_p_U, element id:  " << this->getTag() << endln;
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
Response * EightNodeBrick_u_p_U::setResponse (const char **argv, int argc, Information &eleInformation) 
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

int EightNodeBrick_u_p_U::getResponse (int responseID, Information &eleInfo)
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
void EightNodeBrick_u_p_U::report(char * msg)
  {
    if ( msg ) ::printf("** %s",msg);
    ::printf("\n Element Number = %d\n", this->getTag() );
    ::printf("\n Number of nodes in a EightNodebrick = %d\n",
                                              nodes_in_brick);
    ::printf("\n Determinant of Jacobian (! ==0 before comp.) = %f\n",
                                                  determinant_of_Jacobian);

    ::printf("Node numbers \n");
    ::printf(
".....1.....2.....3.....4.....5.....6.....7.....8.....9.....0.....1.....2\n");
           for ( int i=0 ; i<8 ; i++ ) 
	    //::printf("%6d",G_N_numbs[i]);
	    ::printf("%6d",connectedExternalNodes(i));
    ::printf("\n");
    //           for ( int j=8 ; j<20 ; j++ )
    //             ::printf("%6d",G_N_numbs[j]);	   // Commented by Xiaoyan
    ::printf("\n\n");

    //    ::printf("Node existance array \n");
    //           for ( int k=0 ; k<12 ; k++ )
    //             ::printf("%6d",node_existance[k]);
    //           ::printf("\n\n");			    // Commented by Xiaoyan


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
void EightNodeBrick_u_p_U::reportshort(char * msg)
  {
    if ( msg ) ::printf("** %s",msg);
    ::printf("\n Element Number = %d\n", this->getTag() );
    ::printf("\n Number of nodes in a EightNodeBrick_u_p_U = %d\n",
                                              nodes_in_brick);
    ::printf("\n Determinant of Jacobian (! ==0 before comp.) = %f\n",
                                                  determinant_of_Jacobian);

    ::printf("Node numbers \n");
    ::printf(
".....1.....2.....3.....4.....5.....6.....7.....8.....9.....0.....1.....2\n");
           for ( int i=0 ; i<8 ; i++ )
             //::printf("%6d",G_N_numbs[i]);
             ::printf( "%6d",connectedExternalNodes(i) );
           
	   ::printf("\n");
           //           for ( int j=8 ; j<20 ; j++ )
           //             ::printf("%6d",G_N_numbs[j]);   //// Commented by Xiaoyan
           ::printf("\n\n");

           //    ::printf("Node existance array \n");
           //           for ( int k=0 ; k<12 ; k++ )
           //             ::printf("%6d",node_existance[k]);	   // Commented by Xiaoyan
           ::printf("\n\n");
						 
  }

//#############################################################################
void EightNodeBrick_u_p_U::reportPAK(char * msg)
  {
    if ( msg ) ::printf("%s",msg);
    ::printf("%10d   ",  this->getTag());
           for ( int i=0 ; i<8 ; i++ )
             ::printf( "%6d",connectedExternalNodes(i) );
             //::printf("%6d",G_N_numbs[i]);

    printf("\n");
  }

//#############################################################################
void EightNodeBrick_u_p_U::reportpqtheta(int GP_numb)
  {
    short where = GP_numb-1;
    matpoint[where]->reportpqtheta("");
  }

//#############################################################################
void EightNodeBrick_u_p_U::reportLM(char * msg) 
  {
    if ( msg ) ::printf("%s",msg);
    ::printf("Element # %d, LM->", this->get_Brick_Number());
    for (int count = 0 ; count < 24 ; count++)
      {
        ::printf(" %d", LM[count]);
      }
    ::printf("\n");

  }

//#############################################################################
void EightNodeBrick_u_p_U::reportTensor(char * msg)
  {
    //    if ( msg ) ::printf("** %s\n",msg);
    
    // special case for 8 nodes only
    // special case for 8 nodes only
    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;

    short where = 0;

    // special case for 8 nodes only
    static const int dim[] = {3, 8}; // static-> see ARM pp289-290
    tensor NodalCoord(2, dim, 0.0);
    tensor matpointCoord(2, dim, 0.0);
    int h_dim[] = {24,3};   // Xiaoyan changed from {60,3} to {24,3} for 8 nodes
    tensor H(2, h_dim, 0.0);

    //for (int ncount = 1 ; ncount <= 8 ; ncount++ )
    ////  for (int ncount = 0 ; ncount <= 7 ; ncount++ )
    //  { 
    //	//int global_node_number = get_global_number_of_node(ncount-1);
    //	// printf("global node num %d",global_node_number);
    //
    //    //   NodalCoord.val(1,ncount) = nodes[global_node_number].x_coordinate();
    //    //   NodalCoord.val(2,ncount) = nodes[global_node_number].y_coordinate();
    //    //   NodalCoord.val(3,ncount) = nodes[global_node_number].z_coordinate();
    //    // Xiaoyan changed to the following:  09/27/00
    //	Vector Coordinates = nodes[global_node_number].getCrds();
    //    
    //    NodalCoord.val(1,ncount) = Coordinates(0);
    //    NodalCoord.val(2,ncount) = Coordinates(1);
    //    NodalCoord.val(3,ncount) = Coordinates(2);
    //printf("global point %d     x=%+.6e   y=%+.6e   z=%+.6e \n ", global_node_number, 
    //                                                      NodalCoord.val(1,ncount),
    //						      NodalCoord.val(2,ncount),
    //						      NodalCoord.val(3,ncount));
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
    
    NodalCoord.val(1,1)=nd1Crds(0); NodalCoord.val(2,1)=nd1Crds(1); NodalCoord.val(3,1)=nd1Crds(2);
    NodalCoord.val(1,2)=nd2Crds(0); NodalCoord.val(2,2)=nd2Crds(1); NodalCoord.val(3,2)=nd2Crds(2);
    NodalCoord.val(1,3)=nd3Crds(0); NodalCoord.val(2,3)=nd3Crds(1); NodalCoord.val(3,3)=nd3Crds(2);
    NodalCoord.val(1,4)=nd4Crds(0); NodalCoord.val(2,4)=nd4Crds(1); NodalCoord.val(3,4)=nd4Crds(2);
    NodalCoord.val(1,5)=nd5Crds(0); NodalCoord.val(2,5)=nd5Crds(1); NodalCoord.val(3,5)=nd5Crds(2);
    NodalCoord.val(1,6)=nd6Crds(0); NodalCoord.val(2,6)=nd6Crds(1); NodalCoord.val(3,6)=nd6Crds(2);
    NodalCoord.val(1,7)=nd7Crds(0); NodalCoord.val(2,7)=nd7Crds(1); NodalCoord.val(3,7)=nd7Crds(2);
    NodalCoord.val(1,8)=nd8Crds(0); NodalCoord.val(2,8)=nd8Crds(1); NodalCoord.val(3,8)=nd8Crds(2);


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
                //	       for (int encount=0 ; encount <= 7 ; encount++ ) 
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
//void EightNodeBrick_u_p_U::reportTensor(char * msg)
// ZHaohui added to print gauss point coord. to file fp

void EightNodeBrick_u_p_U::reportTensorF(FILE * fp)
  {
    //if ( msg ) ::printf("** %s\n",msg);

    // special case for 8 nodes only
    // special case for 8 nodes only
    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;

    short where = 0;

    // special case for 8 nodes only
    static const int dim[] = {3, 8}; // static-> see ARM pp289-290
    tensor NodalCoord(2, dim, 0.0);
    tensor matpointCoord(2, dim, 0.0);
    int h_dim[] = {24,3};  // Xiaoyan changed from {60,3} to {24,3} for 8 nodes

    tensor H(2, h_dim, 0.0);

    //for (int ncount = 1 ; ncount <= 8 ; ncount++ )
    //  // for (int ncount = 0 ; ncount <= 7 ; ncount++ )
    //  { 
    //	int global_node_number = get_global_number_of_node(ncount-1);
    //	// printf("global node num %d",global_node_number);
    //
    //    //        NodalCoord.val(1,ncount) = nodes[global_node_number].x_coordinate();
    //    //        NodalCoord.val(2,ncount) = nodes[global_node_number].y_coordinate();
    //    //        NodalCoord.val(3,ncount) = nodes[global_node_number].z_coordinate();
    //    // Xiaoyan changed to the following:  09/27/00
    //	Vector Coordinates = nodes[global_node_number].getCrds();
    //    NodalCoord.val(1,ncount) = Coordinates(0); 
    //    NodalCoord.val(2,ncount) = Coordinates(1); 
    //    NodalCoord.val(3,ncount) = Coordinates(2); 
    //printf("global point %d     x=%+.6e   y=%+.6e   z=%+.6e \n ", global_node_number, 
    //                                                      NodalCoord.val(1,ncount),
    //						      NodalCoord.val(2,ncount),
    //						      NodalCoord.val(3,ncount));
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
    
    NodalCoord.val(1,1)=nd1Crds(0); NodalCoord.val(2,1)=nd1Crds(1); NodalCoord.val(3,1)=nd1Crds(2);
    NodalCoord.val(1,2)=nd2Crds(0); NodalCoord.val(2,2)=nd2Crds(1); NodalCoord.val(3,2)=nd2Crds(2);
    NodalCoord.val(1,3)=nd3Crds(0); NodalCoord.val(2,3)=nd3Crds(1); NodalCoord.val(3,3)=nd3Crds(2);
    NodalCoord.val(1,4)=nd4Crds(0); NodalCoord.val(2,4)=nd4Crds(1); NodalCoord.val(3,4)=nd4Crds(2);
    NodalCoord.val(1,5)=nd5Crds(0); NodalCoord.val(2,5)=nd5Crds(1); NodalCoord.val(3,5)=nd5Crds(2);
    NodalCoord.val(1,6)=nd6Crds(0); NodalCoord.val(2,6)=nd6Crds(1); NodalCoord.val(3,6)=nd6Crds(2);
    NodalCoord.val(1,7)=nd7Crds(0); NodalCoord.val(2,7)=nd7Crds(1); NodalCoord.val(3,7)=nd7Crds(2);
    NodalCoord.val(1,8)=nd8Crds(0); NodalCoord.val(2,8)=nd8Crds(1); NodalCoord.val(3,8)=nd8Crds(2);
      
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
                //	       for (int encount=0 ; encount <= 7 ; encount++ ) 
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
