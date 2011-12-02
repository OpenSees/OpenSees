///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              EightNodeBrick.cpp
// CLASS:             EightNodeBrick
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Zhaohui Yang and Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Zhaohui Yang  and Xiaoyan Wu
// DATE:              Aug. 2000
// UPDATE HISTORY:			 Modified from Brick3D and FourNodeQuad.hh  07/06/00
//																			 Sept. - Oct 2000 connected to OpenSees by Zhaohui
//
//
// CONTACT:           jeremic@ucdavis.edu
///////////////////////////////////////////////////////////////////////////////
//

#include <EightNodeBrick.h>
#define FixedOrder 3


//====================================================================
//Reorganized constructor ____ Zhaohui 02-10-2000
//====================================================================

EightNodeBrick::EightNodeBrick(int element_number,
                               int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
                               int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
                               NDMaterial * Globalmmodel, const char * type, double b1, double b2,
			       double p, double r, EPState *InitEPS)
                               // Set it to 3 //int r_int_order, //int s_int_order, //int t_int_order,
			       //tensor * IN_tangent_E,  //stresstensor * INstress, //stresstensor * INiterative_stress, //double * IN_q_ast_iterative, //straintensor * INstrain):  __ZHaohui 09-29-2000
		               
  :Element(element_number, ELE_TAG_EightNodeBrick ),
  connectedExternalNodes(8), K(24, 24), C(24, 24), M(24, 24), P(24),Q(24), b(2), rho(r), pressure(p)
  {
    //elem_numb = element_number;
    b(0) = b1;
    b(1) = b2;    
    determinant_of_Jacobian = 0.0;
    //r_integration_order = r_int_order; 
    //s_integration_order = s_int_order; 
    //t_integration_order = t_int_order; 
    r_integration_order = FixedOrder; // Gauss-Legendre integration order in r direction
    s_integration_order = FixedOrder; // Gauss-Legendre integration order in s direction
    t_integration_order = FixedOrder; // Gauss-Legendre integration order in t direction

    mmodel = Globalmmodel->getCopy( type ); // One global mat model

    int total_number_of_Gauss_points = r_integration_order*s_integration_order*t_integration_order;
    
    // according to ARM pp.61 default constructor will be called!!
    //IntegrationPoint * MatPoint = new IntegrationPoint[total_number_of_Gauss_points];
    //prebaci sve u jednodimenzioni niz jer samo prvi stepen pointera moze da se pokriva
    //sa onim stosom derived * ->> base * !!

    if ( total_number_of_Gauss_points != 0 )
      {
        //MatPoint = Globalmmodel->new_mm(total_number_of_Gauss_points);
	
        // Zhaohui -- All the following members are moved into gauss-point EPState: gpEPS
	//GPstress = new stresstensor[total_number_of_Gauss_points];
        //GPiterative_stress = new stresstensor[total_number_of_Gauss_points];
        //GPq_ast_iterative  = new double[total_number_of_Gauss_points];
        //GPstrain = new straintensor[total_number_of_Gauss_points];
        //GPtangent_E = new tensor[total_number_of_Gauss_points];//default constructor called
        
	MatPoint  = new IntegrationPoint[total_number_of_Gauss_points];

      }
    else
      {
	//GPstress = 0;//GPiterative_stress = 0;//GPq_ast_iterative  = 0; //GPstrain = 0;//GPtangent_E = 0;
        MatPoint  = 0;
    }
    ////////////////////////////////////////////////////////////////////
    //dakle posto:
    //// according to ARM pp.61 default constructor will be called!!
    //onda oni vec postoje u memoriji i samo im treba dodeliti prave
    //vrednosti iz onog modela koji je prenesen unutra. Znaci onInitializxe funkcija
    //sa modelom kao argumentom Initialize(taj_Model).
    //za stresstensor i straintensor isto to napravi
    //tu Initialize funkciju pa da uzima argument i da samo
    //koristi njegove vrednosti pa stavi u vec postojece mesto
    //u memoriji ove vrednsoti. Takodje unutar brick3d napravi ih
    //( te stresstensor , straintensor i mmodel ) static da ostanu tu
    //da se ne pozove destructor . . .

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

                //DB::printf("\n\nBefore Initialization **************** where = %d \n",where);
                //DB::printf("GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //DB            GP_c_r,GP_c_s,GP_c_t);
                //DB
                //DBGPstress[where].reportshort("stress within before Initialization");
                //DBGPstrain[where].reportshort("strain within before Initialization");
                //DB
                //DB// I suspect that [] should be overloaded so that compiler knows which
                //DB// material model is returning a pointer and fot the purpose
                //DB//MatPoint[where].report("mmodel within before Initialization");
                //DB//MatPoint[where].report("mmodel within before Initialization"); // operator[] overloaded
                //DB(MatPoint)->operator[](where).report("mmodel within before Initialization"); // operator[] overloaded
                //DB                                                               // for NDMaterial and
                //DB                                                               // derived types!

		// Zhaohui  09-30-2000
                //GPtangent_E[where].Initialize_all(*IN_tangent_E);
                //GPstress[where].Initialize(*INstress);
                //GPiterative_stress[where].Initialize(*INiterative_stress);
                //GPq_ast_iterative[where] = IN_q_ast_iterative[where];
                //GPstrain[where].Initialize(*INstrain);
                
                // Already assigned the copy of Globalmmodel to each element of MatPoint
		//(MatPoint)->operator[](where).Initialize(*Globalmmodel);

                NDMaterial *NMD = Globalmmodel;
		if (strcmp( Globalmmodel->getType(), "Template3Dep" ) == 0)
 		   NMD = 0;

                MatPoint[where].Initialize(GP_c_r,
                                         GP_c_s,
                                         GP_c_t,
                                         r, s, t,
                                         rw, sw, tw,
                                         InitEPS, 
					 NMD);
					 //&( GPstress[where] ), //&( GPiterative_stress[where] ), //IN_q_ast_iterative[where] ,//&( GPstrain[where] ),  //&( GPtangent_E[where] ),
                                         //&( (MatPoint)->operator[](where) )
                                          // ugly syntax but it works! Still don't know what's wrong   // with the old style MatPoint[where]
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

      //nodes = GlobalNodes;

      // loop nodes 9-20 and :
      //  1) put the right coordinates on place,
      //  2) calculate the total number of nodes

      nodes_in_brick = 8;
      
}

//====================================================================
EightNodeBrick::EightNodeBrick ():Element(0, ELE_TAG_EightNodeBrick ),
connectedExternalNodes(8), K(24, 24), C(24, 24), M(24, 24), P(24),Q(24), b(2), rho(0.0), pressure(0.0), mmodel(0)
{
     MatPoint = 0;
}   


// Initialize
// default constructor
// this one takes only one material model and forwards it to all
// Gauss points.
//#############################################################################

int EightNodeBrick::Initialize(int element_number,
                                int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
                                int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
                                NDMaterial * Globalmmodel, const char * type, double b1, double b2,
                                double p, double r, EPState * InitEPS
			        //tensor       * IN_tangent_E,//stresstensor * INstress, //stresstensor * INiterative_stress, //double  * IN_q_ast_iterative, //straintensor * INstrain
			        )   
  {
    //elem_numb = element_number;
    b(0) = b1;
    b(1) = b2;
    rho = r;
    pressure = p;
    // r_integration_order = r_int_order; 
    // s_integration_order = s_int_order; 
    // t_integration_order = t_int_order; 
    r_integration_order = FixedOrder; // Gauss-Legendre integration order in r direction
    s_integration_order = FixedOrder; // Gauss-Legendre integration order in s direction
    t_integration_order = FixedOrder; // Gauss-Legendre integration order in t direction

    mmodel = Globalmmodel->getCopy(); // One global mat model

    int total_number_of_Gauss_points = r_integration_order*
                                       s_integration_order*
                                       t_integration_order;

    //MatPoint = Globalmmodel->new_mm(total_number_of_Gauss_points);
    //MatPoint = Globalmmodel->getCopy(total_number_of_Gauss_points);
    
    //GPstress = new stresstensor[total_number_of_Gauss_points];
    //GPiterative_stress = new stresstensor[total_number_of_Gauss_points];
    //GPq_ast_iterative  = new double[total_number_of_Gauss_points];
    //GPstrain = new straintensor[total_number_of_Gauss_points];
    //GPtangent_E = new tensor[total_number_of_Gauss_points];
    //MatPoint  = new IntegrationPoint[total_number_of_Gauss_points];

    MatPoint = new IntegrationPoint[total_number_of_Gauss_points];

    ////////////////////////////////////////////////////////////////////
    //dakle posto:
    //// according to ARM pp.61 default constructor will be called!!
    //onda oni vec postoje u memoriji i samo im treba dodeliti prave
    //vrednosti iz onog modela koji je prenesen unutra. Znaci onInitializxe funkcija
    //sa modelom kao argumentom Initialize(taj_Model).
    //za stresstensor i straintensor isto to napravi
    //tu Initialize funkciju pa da uzima argument i da samo
    //koristi njegove vrednosti pa stavi u vec postojece mesto
    //u memoriji ove vrednsoti. Takodje unutar brick3d napravi ih
    //( te stresstensor , straintensor i mmodel ) static da ostanu tu
    //da se ne pozove destructor . . .

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

                //DB::printf("\n\nBefore Initialization **************** where = %d \n",where);
                //DB::printf("GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //DB            GP_c_r,GP_c_s,GP_c_t);
                //DB
                //DBGPstress[where].reportshort("stress within before Initialization");
                //DBGPstrain[where].reportshort("strain within before Initialization");
                //DB
                //DB// I suspect that [] should be overloaded so that compiler knows which
                //DB// material model is returning a pointer and fot the purpose
                //DB//MatPoint[where].report("mmodel within before Initialization");
                //DB//MatPoint[where].report("mmodel within before Initialization"); // operator[] overloaded
                //DB(MatPoint)->operator[](where).report("mmodel within before Initialization"); // operator[] overloaded
                //DB                                                               // for NDMaterial and
                //DB                                                               // derived types!
                
                //.MatPoint[where].Initialize(GP_c_r,
                //.                         GP_c_s,
                //.                         GP_c_t,
                //.                         r, s, t,
                //.                         rw, sw, tw,
                //.                         &( GPstress[where] ),
                //.                         &( GPstrain[where] ),
                //.                         &( (MatPoint)->operator[](where) )
                //.                        );      // ugly syntax but it works!
                //.                                // still don't know what's wrong
                //.                                // with the old style MatPoint[where]
                // Initialize it with the elastic stiffness tensor at the begining
                //````double Ey = (MatPoint)->operator[](where).E();
                //````double nu = (MatPoint)->operator[](where).nu();
                //````tensor Constitutive = (MatPoint)->operator[](where).StiffnessTensorE(Ey,nu);
                //````
                //````GPtangent_E[where].Initialize_all(Constitutive);
                
		//GPtangent_E[where].Initialize_all(*IN_tangent_E);                
                //GPstress[where].Initialize(*INstress);                
                //GPiterative_stress[where].Initialize(*INiterative_stress);
                //GPq_ast_iterative[where] = IN_q_ast_iterative[where];	  
                //GPstrain[where].Initialize(*INstrain);
                
		//(MatPoint)->operator[](where).Initialize(*Globalmmodel);
		// Already initialized using Globalmmodel in getCopy() --Zhaohui 09-29-2000

                NDMaterial *NMD = Globalmmodel;
		if (strcmp( Globalmmodel->getType(), "Template3Dep" ) == 0)
 		   NMD = 0;

                MatPoint[where].Initialize(GP_c_r,
                                           GP_c_s,
                                           GP_c_t,
                                           r, s, t,
                                           rw, sw, tw,
                                           InitEPS,
					   NMD );
                                           // &( GPstress[where] ), //&( GPiterative_stress[where] ), // IN_q_ast_iterative[where],
                                           // &( GPstrain[where] ),  //&( GPtangent_E[where] ),  // ZHaohui 09-29-2000
                                           // &( (MatPoint)->operator[](where) )
                                           // ugly syntax but it works! // still don't know what's wrong   // with the old style MatPoint[where]
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
      
      //G_N_numbs[0] = node_numb_1;
      //G_N_numbs[1] = node_numb_2;
      //G_N_numbs[2] = node_numb_3;
      //G_N_numbs[3] = node_numb_4;
      //G_N_numbs[4] = node_numb_5;
      //G_N_numbs[5] = node_numb_6;
      //G_N_numbs[6] = node_numb_7;
      //G_N_numbs[7] = node_numb_8;
      
      //nodes = GlobalNodes;

      // loop nodes 9-20 and :
      //  1) put the right coordinates on place,
      //  2) calculate the total number of nodes
      // loop nodes 9-20 and :
      //  1) put the right coordinates on place,
      //  2) calculate the total number of nodes
      nodes_in_brick = 8;
      //    for ( int i=8 ; i<20 ; i++ )
      //      {
      //        if (G_N_numbs[i] == 0 )
      //          {
      //                   node_existance[i-8] = 0;
      //          }
      //        else
      //          {
      //            nodes_in_brick++;
      //            node_existance[i-8] = 1;
      //          }
      //      }
      //    Node * N = new Node[nodes_in_brick];
      //  // Xiaoyan comment for only 8 nodes  07/11/00
      return 0;
  }

//#############################################################################



EightNodeBrick::~EightNodeBrick ()
{

    //int order = theQuadRule->getOrder();
    //
    //for (int i = 0; i < order; i++) {
    //	for (int j = 0; j < order; j++)
    //	    // Delete the NDMaterials at each integration point
    //	    if (theMaterial[i][j])
    //		delete theMaterial[i][j];
    //	
    //	// Delete each array of NDMaterial pointers
    //	if (theMaterial[i])
    //	    delete [] theMaterial[i];
    //}	
    //
    //// Delete the array of pointers to NDMaterial pointer arrays
    //if (theMaterial)
    //	delete [] theMaterial;
    //
    //// Delete the quadrature rule
    //if (theQuadRule)
    //	delete theQuadRule;

}
     // comment by Xiaoyan 07/11/00


void EightNodeBrick::incremental_Update()
  {
    double r  = 0.0;
    // double rw = 0.0;
    double s  = 0.0;
    // double sw = 0.0;
    double t  = 0.0;
    // double tw = 0.0;

    short where = 0;
    //,,,,,    double weight = 0.0;
    
    //double this_one_PP = (MatPoint)->operator[](where).IS_Perfect_Plastic();

    int dh_dim[] = {8,3};   //Xiaoyan changed from {20,3} to {8,3}  07/12/00
    tensor dh(2, dh_dim, 0.0);

    tensor Constitutive( 4, def_dim_4, 0.0);

    //    double det_of_Jacobian = 0.0;

    static int disp_dim[] = {8,3}; //Xiaoyan changed from {20,3} to {8,3}  07/12/00
    tensor incremental_displacements(2,disp_dim,0.0);

    straintensor incremental_strain;
    straintensor total_strain_at_GP;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    //....    int number_of_subincrements = 1;
    //....    double this_one_PP = 1.0; // if set to 0.0 -> perfectly plastic
    //....                              // if set to 1.0 -> elasto plastic

    stresstensor final_stress_after_integration;
    
    ///    stresstensor incremental_stress;
    // tensor of incremental displacements taken from node objects
    incremental_displacements = incr_disp();

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
                dhGlobal = dh("ij") * JacobianINV("jk");
                //....                dhGlobal.print("dh","dhGlobal");
                //weight
                //                weight = rw * sw * tw * det_of_Jacobian;
                //::::::   ::printf("\n\nIN THE STIFFNESS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                //::::::   ::printf(" void EightNodeBrick::incremental_Update()\n");
                //::::::   ::printf(" GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d    --->>>  where = %d \n",
                //::::::                      GP_c_r,GP_c_s,GP_c_t,where);
                //::::::   ::printf("WEIGHT = %f", weight);
                //::::::   ::printf("determinant of Jacobian = %f", determinant_of_Jacobian);
                //::::::   MatPoint[where].report("Gauss Point\n");
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
		//(MatPoint)->operator[](where).kappa_set(final_stress_after_integration,
                //                                 GPq_ast_iterative[where]);
                
		//....         final_stress_after_integration =
                //....           (MatPoint)->operator[](where).FinalStress(*(GPstress+where),
                //....                                                     incremental_strain,
                //....                                                     (MatPoint)->operator[](where),
                //....                                                     number_of_subincrements,
                //....                                                     this_one_PP);
                //....//final_stress_after_integration.reportshortpqtheta("\n final_stress_after_integration GAUSS \n");
                // calculate the constitutive tensor

                // We do not need: final_stress_after_integration

	        
	        //Constitutive =
                //  (MatPoint)->operator[](where).ConstitutiveTensor(final_stress_after_integration,
                //                                                   *(GPstress+where),
                //                                                   incremental_strain,
                //                                                   (MatPoint)->operator[](where),
                //                                                   this_one_PP);
	        // ZHaohui modified __09-29-2000
		
		//EPState trialEPS = *(mmodel->getEPS());
		
	        EPState *tmp_eps = (MatPoint[where]).getEPS();
	        NDMaterial *tmp_ndm = (MatPoint[where]).getNDMat();

		if ( tmp_eps ) { //if there is an EPState for the MatPoint
		  mmodel->setEPS( *tmp_eps );
		  if ( ! (mmodel->setTrialStrainIncr( incremental_strain)) )
               	     g3ErrorHandler->warning("EightNodeBrick::incremental_Update (tag: %d), not converged",
					 this->getTag());
		  MatPoint[where].setEPS( mmodel->getEPS() );
		}
		else if ( tmp_ndm ) 
		  (MatPoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
		else {
               	   g3ErrorHandler->fatal("EightNodeBrick::incremental_Update (tag: %d), no strain or stress state vars", this->getTag());
		   exit(1);
		}
		   	        
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

//Zhaohui: We never need this one now for it's combined in FE(BE)EPStates 02-10-2000

////#############################################################################
//// update iterative_stresses during iterative process
//// so that function iterative_nodal_forces can return right internal forces
//// (iterative in this case ) !!!!
////#############################################################################
//void EightNodeBrick::iterative_Update()
//  {
//    double r  = 0.0;
////    double rw = 0.0;
//    double s  = 0.0;
////    double sw = 0.0;
//    double t  = 0.0;
////    double tw = 0.0;
//
//    short where = 0;
////,,,,,    double weight = 0.0;
//
//    int dh_dim[] = {8,3}; // Xiaoyan changed from {20,3} to {8,3} for 8 nodes brick
//    tensor dh(2, dh_dim, 0.0);
//
////    tensor Constitutive( 4, def_dim_4, 0.0);
//
////    double det_of_Jacobian = 0.0;
//
//    static int disp_dim[] = {8,3};  // Xiaoyan changed from {20,3} to {8,3} for 8 nodes brick
//
//    tensor incremental_displacements(2,disp_dim,0.0);
//
//    straintensor incremental_strain;
////    straintensor total_strain_at_GP;
//
//    tensor Jacobian;
//    tensor JacobianINV;
//    tensor dhGlobal;
//
//    int number_of_subincrements = 1;
//    //double this_one_PP = (MatPoint)->operator[](where).IS_Perfect_Plastic();
//
//    stresstensor final_stress_after_integration;
//    //    stresstensor incremental_stress;
//    // tensor of incremental displacements taken from node objects for this element !
//    incremental_displacements = incr_disp();
//    //incremental_displacements.print("disp","\n incremental_displacements tensor at GAUSS point in iterative_Update\n");
//
//    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
//      {
//        r = get_Gauss_p_c( r_integration_order, GP_c_r );
//        //--        rw = get_Gauss_p_w( r_integration_order, GP_c_r );
//        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
//          {
//            s = get_Gauss_p_c( s_integration_order, GP_c_s );
//            //--            sw = get_Gauss_p_w( s_integration_order, GP_c_s );
//            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
//              {
//                t = get_Gauss_p_c( t_integration_order, GP_c_t );
//                //--                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
//                // this short routine is supposed to calculate position of
//                // Gauss point from 3D array of short's
//                where =
//                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
//                // derivatives of local coordiantes with respect to local coordiantes
//                dh = dh_drst_at(r,s,t);
//                // Jacobian tensor ( matrix )
//                Jacobian = Jacobian_3D(dh);
//                //Jacobian.print("J");
//                // Inverse of Jacobian tensor ( matrix )
//                JacobianINV = Jacobian_3Dinv(dh);
//                //JacobianINV.print("JINV");
//                // determinant of Jacobian tensor ( matrix )
//                //--                det_of_Jacobian  = Jacobian.determinant();
//                //::printf("determinant of Jacobian is %f\n",Jacobian_determinant );
//                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
//                dhGlobal = dh("ij") * JacobianINV("jk");
//                //....                dhGlobal.print("dh","dhGlobal");
//                //weight
//                //                weight = rw * sw * tw * det_of_Jacobian;
//                //::::::   ::printf(" void EightNodeBrick::iterative_Update() \n ");
//                //::::::   ::printf(" GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d    --->>>  where = %d \n",
//                //::::::                      GP_c_r,GP_c_s,GP_c_t,where);
//                //::::::   ::printf("WEIGHT = %f", weight);
//                //::::::   ::printf("determinant of Jacobian = %f", determinant_of_Jacobian);
//                //::::::   MatPoint[where].report("Gauss Point\n");
//                incremental_strain =
//                   (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
//                incremental_strain.null_indices();
//                //incremental_strain.reportshort("incremental_strain in void EightNodeBrick::iterative_Update()");
//
//                //..   dakle ovde posalji strain_increment jer se stari stress cuva u okviru svake
//                //..   Gauss tacke a samo saljes strain_increment koji ce da se prenese
//                //..   u integracionu rutinu pa ce ta da vrati krajnji napon i onda moze da
//                //..   se pravi ConstitutiveStiffnessTensor.
//                // here comes the final_stress calculation
//                //(GPstress+where)->reportshortpqtheta("\n stress START GAUSS  in iterative_Update\n");
//         
//	        // Zhaohui commented out_____02-10-2000
//		final_stress_after_integration =
//                    (MatPoint)->operator[](where).FinalStress(*(GPstress+where),
//                                                     incremental_strain,
//                                                     (MatPoint)->operator[](where),
//                                                     number_of_subincrements,
//                                                     this_one_PP);
//                
//		final_stress_after_integration.reportshort("\n final_stress_after_integration in void EightNodeBrick::iterative_Update()\n");
//                
//		//----// intialize total strain with the strain at this Gauss point before
//                //----// adding this increments strains!
//                //----                total_strain_at_GP.Initialize(*(GPstrain+where));
//                //----//total_strain_at_GP.reportshort("\n total_strain tensor at GAUSS point \n");
//                //----// this is the addition of incremental strains to the previous strain state at
//                //----// this Gauss point
//                //----                total_strain_at_GP = total_strain_at_GP + incremental_strain;
//                //----//total_strain_at_GP.reportshort("\n total_strain tensor at GAUSS point AFTER\n");
//                // calculate the constitutive tensor
//            
//	        // this is update of iterative_stress state at this Gauss point
//                GPiterative_stress[where].Initialize(final_stress_after_integration);
//                GPq_ast_iterative[where] =
//                  (MatPoint)->operator[](where).kappa_get(final_stress_after_integration);
//                //             GPstress[where].Initialize(final_stress_after_integration);
//                //GPiterative_stress[where].reportshortpqtheta("\n iterative_stress at GAUSS point in iterative_Update\n");
//                GPiterative_stress[where].reportshort("\n iterative_stress at GAUSS point in iterative_Update\n");
//                //----    // this is update of strain state at this Gauss point
//                //----                GPstrain[where].Initialize(total_strain_at_GP);
//                //GPstrain[where].reportshort("\n strain at GAUSS point \n");
//              }
//          }
//      }
//  }




//#############################################################################
//#############################################################################
//***************************************************************
tensor EightNodeBrick::H_3D(double r1, double r2, double r3)
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

//#############################################################################
//***************************************************************
tensor EightNodeBrick::interp_poli_at(double r1, double r2, double r3)
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



tensor EightNodeBrick::dh_drst_at(double r1, double r2, double r3)
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


////#############################################################################
//Finite_Element * EightNodeBrick::new_el(int total)
//  {
//    EightNodeBrick *el_p;
//    el_p = new EightNodeBrick[total];
//    //DB//-------------------------------------------
//    //DB    for ( int i=0 ; i<total ; i++ )
//    //DB      {
//    //DB        el_p[i].report("derived EightNodeBrick\n");
//    //DB      }
//    //DB//-------------------------------------------
//    return el_p;
//  }

////#############################################################################
Finite_Element & EightNodeBrick::operator[](int subscript)
  {
    return ( *(this+subscript) );
  }

//Finite_Element & EightNodeBrick::operator[](short subscript)
//  {
//    return ( *(this+subscript) );
//  }

//Finite_Element & EightNodeBrick::operator[](unsigned subscript)
//  {
//    return ( *(this+subscript) );
//  }


////#############################################################################
tensor EightNodeBrick::getStiffnessTensor(void)
  {
    int K_dim[] = {8,3,3,8};  // Xiaoyan changed from {20,3,3,20} to {8,3,3,8} 
			      // for 8 nodes      07/12
    tensor Kk(4,K_dim,0.0);

    //debugging
    matrix Kmat;

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
    straintensor total_strain_at_GP;

    tensor Jacobian;
    tensor JacobianINV;
    tensor JacobianINVtemp;
    tensor dhGlobal;

    //int number_of_subincrements = 1;
    //double this_one_PP = 1.0; // if set to 0.0 -> perfectly plastic
                              // if set to 1.0 -> elasto plastic
    //    tensor Ktemp(4,K_dim,0.0);
    //    char * integr_type = 0;

    stresstensor final_stress_after_integration;
    stresstensor incremental_stress;
    // tensor of incremental displacements taken from node objects
    incremental_displacements = incr_disp();
    //..    incremental_displacements.print("incr_disp");

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
                dhGlobal = dh("ij") * JacobianINV("jk");
                //        ::fprintf(stdout," # %d \n\n\n\n\n\n\n\n", El_count);
		//dhGlobal.print("dhGlobal");
                //weight
                weight = rw * sw * tw * det_of_Jacobian;
                //::::::
                //printf("\n\nIN THE STIFFNESS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                //printf("  Stifness_Tensor \n");
                //printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //                            GP_c_r,GP_c_s,GP_c_t);
                //printf("WEIGHT = %f", weight);
                //Jacobian.print("J");
                //JacobianINV.print("JINV");
                //JacobianINVtemp.print("JINVtemp");
                //tensor I = JacobianINV("ij")*Jacobian("jk");
                //I.print("I");

                //printf("determinant of Jacobian = %.6e", det_of_Jacobian);
                //MatPoint[where].report("Gauss Point\n");

                // incremental straines at this Gauss point
                //GPstress[where].reportshortpqtheta("\n stress at GAUSS point in stiffness_tensor1\n");
                
		incremental_strain =
                     (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
                //incremental_strain.null_indices();
                //incremental_strain.report("\n incremental_strain tensor at GAUSS point\n");
                
		// incremental_strain.reportshort("\n incremental_strain tensor at GAUSS point\n");
                //----   GPstress[where].reportshortpqtheta("\n stress at GAUSS point in stiffness_tensor2\n");
                // intialize total strain with the strain at this Gauss point before
                // adding this increments strains!
                //                total_strain_at_GP.Initialize(*(GPstrain+where));
                //total_strain_at_GP.reportshort("\n total_strain tensor at GAUSS point BEFORE\n");
                // this is the addition of incremental strains to the previous strain state at
                // this Gauss point
                //                total_strain_at_GP = total_strain_at_GP + incremental_strain;
                //total_strain_at_GP.reportshort("\n total_strain tensor at GAUSS point AFTER\n");
                //..   dakle ovde posalji strain_increment jer se stari stress cuva u okviru svake
                //..   Gauss tacke a samo saljes strain_increment koji ce da se prenese
                //..   u integracionu rutinu pa ce ta da vrati krajnji napon i onda moze da
                //..   se pravi ConstitutiveStiffnessTensor.
                // here comes the final_stress calculation
                // this final stress after integration is used ONLY to obtain Constitutive tensor
                // at this Gauss point.
                
                //=====================================================================================
                // Need to change:  to use new Template3Dep 
                // ZHaohui 09-27-2000   xxxxxxxxxxxxxxxxxxxxxxxxx

                //final_stress_after_integration =
                //    (MatPoint)->operator[](where).FinalStress(*(GPstress+where),
                //                                 incremental_strain,
                //                                 (MatPoint)->operator[](where),
                //                                 number_of_subincrements,
                //                                 this_one_PP);
                //final_stress_after_integration.reportshortpqtheta("\n final_stress_after_integration in stiffness_tensor5\n");

                //----   GPstress[where].reportshortpqtheta("\n stress at GAUSS point in stiffness_tensor3\n");
                //final_stress_after_integration.reportshortpqtheta("\n final_stress_after_integration GAUSS \n");
                //----   GPstress[where].reportshortpqtheta("\n stress at GAUSS point in stiffness_tensor4\n");

                // this final stress after integration is used ONLY to obtain Constitutive tensor
                // at this Gauss point AND to set up the iterative_stress that is used in calculting
                // internal forces during iterations!!

                //GPiterative_stress[where].Initialize(final_stress_after_integration);
                //GPiterative_stress[where].reportshortpqtheta("\n iterative_stress at GAUSS point in stiffness_tensor5\n");


                // Stress state at Gauss point will be updated ( in the
                // sense of updating stresses ( integrating them ) ) in function Update (...) ! ! ! ! !
                // calculate the constitutive tensor
                //......         Constitutive =  GPtangent_E[where];

		//Constitutive =  (MatPoint)->operator[](where).ConstitutiveTensor(final_stress_after_integration,
                //                                         *(GPstress+where),
                //                                          incremental_strain,
                //                                          (MatPoint)->operator[](where),
                //                                          this_one_PP);
                //Constitutive.print("C","\n\n C tensor \n");
                
	        EPState *tmp_eps = (MatPoint[where]).getEPS();
	        NDMaterial *tmp_ndm = (MatPoint[where]).getNDMat();

		if ( tmp_eps ) {     //Elasto-plastic case
		    mmodel->setEPS( *tmp_eps );
		    if ( ! (mmodel->setTrialStrainIncr( incremental_strain)) )
               	       g3ErrorHandler->warning("EightNodeBrick::incremental_Update (tag: %d), not converged",
					 this->getTag());
		    Constitutive = mmodel->getTangentTensor();
      		    MatPoint[where].setEPS( mmodel->getEPS() );
		}
		else if ( tmp_ndm ) { //Elastic case
		    (MatPoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
		    Constitutive = (MatPoint[where].p_matmodel)->getTangentTensor();
		}
		else {
               	   g3ErrorHandler->fatal("EightNodeBrick::incremental_Update (tag: %d), could not getTangentTensor", this->getTag());
		   exit(1);
		}
			  
		//printf("Constitutive.trace = %12.6e\n", Constitutive.trace());
                //Kmat = this->stiffness_matrix(Constitutive);
                //printf("Constitutive tensor max:= %10.3e\n", Kmat.mmax());

                //----   GPstress[where].reportshortpqtheta("\n stress at GAUSS point in stiffness_tensor5\n");
                // this is update of constitutive tensor at this Gauss point
                //GPtangent_E[where].Initialize(Constitutive);
                //....GPtangent_E[where].print(" tangent E at GAUSS point");

                //GPstress[where].reportshortpqtheta("\n stress at GAUSS point in stiffness_tensor6\n");

                //K = K + temp2;
                Kk = Kk + dhGlobal("ib")*Constitutive("abcd")*dhGlobal("jd")*weight;
                //....K.print("K","\n\n K tensor \n"); 
                
		//Kmat = this->stiffness_matrix(Kk);
                //printf("K tensor max= %10.3e\n", Kmat.mmax());

                //convert constitutive and K to matrix and find min and max and print!



              }
          }
      }
    //K = Kk;
    return Kk;
  }


//#############################################################################

void EightNodeBrick::set_strain_stress_tensor(FILE *fp, double * u)
  {
    int dh_dim[] = {8,3};   // Xiaoyan changed from {20,3} to {8,3}
    tensor dh(2, dh_dim, 0.0);

//    tensor Constitutive( 4, def_dim_4, 0.0);
    tensor Constitutive;
    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;
    int where = 0;

    double det_of_Jacobian = 0.0;

    straintensor strain;
    stresstensor stress;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;


    static int disp_dim[] = {8,3};    // Xiaoyan changed from {20,3} to {8,3}
    tensor total_displacements(2,disp_dim,0.0); //

    total_displacements = total_disp(fp, u);

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
                dhGlobal = dh("ij") * JacobianINV("jk");
                //weight
                // straines at this Gauss point from displacement
                strain = (dhGlobal("ib")*total_displacements("ia")).symmetrize11();
                strain.null_indices();
                // here comes the final_stress calculation
                // at this Gauss point.

                //Constitutive =  GPtangent_E[where];
                //Constitutive =  (MatPoint->getEPS() )->getEep();
                // if set total displ, then it should be elstic material
		Constitutive =  ( MatPoint->getNDMat() )->getTangentTensor();

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
//  tensor EightNodeBrick::mass_tensor(Elastic  mmodel)
tensor EightNodeBrick::getMassTensor(void)
  {
    int M_dim[] = {24,24};    // Xiaoyan changed from {60,60} to {24,24}
    tensor Mm(2,M_dim,0.0);

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

    int h_dim[] = {24,3};	// Xiaoyan changed from {60,3} to {24,3}
    //    int h_dim[] = {20,3};
    tensor H(2, h_dim, 0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;

    double RHO;
    RHO= rho; 

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
                //                dhGlobal = dh("ij") * JacobianINV("jk");
                // derivatives of local coordinates with respect to local coordinates


                // printf("\n\nIN THE MASS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                // printf("  Mass_Tensor \n");
                // printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //                             GP_c_r,GP_c_s,GP_c_t);
                // 
                H = H_3D(r,s,t);

                //	double sum = 0.0;
                //	for (int i=1; i<=60 ; i++)
                //           {
                // //  	    sum+=H.cval(i,1);
                // 	    for (int j=1; j<= 3; j++)
                // 	       {
                //        	          sum+=H.cval(i,j);
                // 	          ::printf( "  %+9.2e", H.cval(i,j) );
                // 	        }
                //             ::printf( "  %d \n", i);
                // 	   }
                // 	    ::printf( " \n sum= %+6.2e\n", sum );


     

                // IntegrationPoint GaPo = IntegrationPoint::GP()+where;

		// if it's elasto-plastic, then use the rho at each integration point
		if ( MatPoint[where].getEPS() ) 
		   RHO= (MatPoint[where].getEPS())->getrho();  

		//printf("RHO = %10.3e",RHO);
                //		getchar();

                //weight
                weight = rw * sw * tw * RHO * det_of_Jacobian;
  	        //	printf("weight = %6.2e \n",weight);

		//M.print("M","BEFORE");
                
	        //	tensor temp = H("ib")*H("kb");
		//temp.print("t","temporary tensor H(\"ib\")*H(\"kb\") \n\n" );

		Mm = Mm + H("ib")*H("kb")*weight;
	       //	printf("\n +++++++++++++++++++++++++ \n\n");
	      //	M.print("M","AFTER");
              }
          }
      }
    //M = Mm;
    return Mm;
  }


////#############################################################################
double EightNodeBrick::Potential_Energy(void)
  {
//    double Potential_Energy_Estimate = 0.0;
    double Delta_Potential_Energy_Estimate = 0.0;

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};	// Xiaoyan changed from {20,3} to {8,3}
    tensor dh(2, dh_dim, 0.0);

    static int disp_dim[] = {8,3};	// Xiaoyan changed from {20,3} to {8,3}
    tensor incremental_displacements(2,disp_dim,0.0); // \Delta u

    double det_of_Jacobian = 0.0;

    straintensor incremental_strain;
//    straintensor total_strain_at_GP;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;

    stresstensor stress_sum;
// tensor of incremental displacements taken from node objects
    incremental_displacements = incr_disp();

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
                //....                Jacobian.print("J");
                // Inverse of Jacobian tensor ( matrix )
                JacobianINV = Jacobian_3Dinv(dh);
                //....                JacobianINV.print("JINV");
                // determinant of Jacobian tensor ( matrix )
                det_of_Jacobian  = Jacobian.determinant();
                //....  ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );
                // Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
                dhGlobal = dh("ij") * JacobianINV("jk");
                //weight
                weight = rw * sw * tw * det_of_Jacobian;
                //::::::   ::printf("\n\nIN THE STIFFNESS TENSOR INTEGRATOR ----**************** where = %d \n", where);
                //::fprintf(stdout," Potential_Energy\n");
                //::fprintf(stdout,"                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //                                      GP_c_r,GP_c_s,GP_c_t);
                //::fprintf(stdout,"WEIGHT = %f", weight);
                //::fprintf(stdout,"determinant of Jacobian = %f", determinant_of_Jacobian);
                //::::::   MatPoint[where].report("Gauss Point\n");
                //::::::
                // incremental straines at this Gauss point
                //----   GPstress[where].reportshortpqtheta("\n stress at GAUSS point in stiffness_tensor1\n");
                incremental_strain =
                 (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
                incremental_strain.null_indices();
                //GPstress[where].reportshortpqtheta("\n stress at GAUSS point in Potential_Energy\n");
                //GPiterative_stress[where].reportshortpqtheta("\n ITERATIVE stress at GAUSS point in Potential_Energy\n");

                //stress_sum = GPstress[where] + GPiterative_stress[where];
                stress_sum = (MatPoint[where].getEPS())->getStress() 
		                + (MatPoint[where].getEPS())->getStress_commit();
                
		//stress_sum.reportshortpqtheta("\n ITERATIVE stress_SUM at GAUSS point in Potential_Energy\n");

                Delta_Potential_Energy_Estimate += weight * (stress_sum("ij")*incremental_strain("ij")).trace();
                //         Potential_Energy_Estimate =

              }
          }
      }
//::fprintf(stdout,"Delta_Potential_Energy_Estimate = %.20e \n",Delta_Potential_Energy_Estimate);
    return Delta_Potential_Energy_Estimate;
  }


////#############################################################################

tensor EightNodeBrick::stiffness_matrix(const tensor & K)
  {
//    int K_dim[] = {20,3,3,20};
//    tensor K(4,K_dim,0.0);
    matrix Kmatrix(24,24,0.0);	 // Xiaoyan changed from (60,60,0,0) to (24,24,0,0)

    int Ki=0;
    int Kj=0;

    for ( int i=1 ; i<=8 ; i++ )  // Xiaoyan changed from i<=20 to i<=8 for 8 nodes
      {
        for ( int j=1 ; j<=8 ; j++ )  // Xiaoyan changed from i<=20 to i<=8 for 8 nodes
          {
            for ( int k=1 ; k<=3 ; k++ )
              {
                for ( int l=1 ; l<=3 ; l++ )
                  {
                    Ki = k+3*(i-1);
                    Kj = l+3*(j-1);
                    //::printf("i=%d k=%d  Ki=%d       j=%d l=%d  Kj=%d\n",i,k,Ki,j,l,Kj);

                    Kmatrix.val( Ki , Kj ) = K.cval(i,k,l,j);
                  }
              }
          }
      }
    return Kmatrix;
  }

//#############################################################################

////#############################################################################
// Constructing mass matrix from mass tensor ___Zhaohui 07-05-99
tensor EightNodeBrick::mass_matrix(const tensor & M)
  {
    //    int K_dim[] = {20,3,3,20};
    //    tensor K(4,K_dim,0.0);
    matrix Mmatrix(24,24,0.0);	// Xiaoyan changed from (60,60,0,0) to (24,24,0,0) for 8 nodes

    for ( int i=1 ; i<=24 ; i++ )   // Xiaoyan changed from i<=60 to i<=24 for 8 nodes
      {
        for ( int j=1 ; j<=24 ; j++ ) // Xiaoyan changed from i<=60 to i<=24 for 8 nodes
          {
             Mmatrix.val( i , j ) = M.cval(i,j);
             //  ::printf("Mi Mj %d %d %+6.2e ",Mi,Mj,Mmatrix.val( Mi , Mj ) );
          }
      }
    return Mmatrix;
  }
////#############################################################################

////#############################################################################
tensor EightNodeBrick::Jacobian_3D(tensor dh)
  {                       
     //       dh ( 20*3)  // dh(8*3) Xiaoyan
     tensor N_C = Nodal_Coordinates(); // 20*3	  // 8*3 Xiaoyan
     tensor Jacobian_3D = dh("ij") * N_C("ik");
     return Jacobian_3D;
  }

//#############################################################################
tensor EightNodeBrick::Jacobian_3Dinv(tensor dh)
  {                       
     //       dh ( 20*3)	  // dh(8*3) Xiaoyan  
     tensor N_C = Nodal_Coordinates(); // 20*3	  	  // 8*3 Xiaoyan   
     tensor Jacobian_3Dinv = (dh("ij") * N_C("ik")).inverse();
     return Jacobian_3Dinv;
  }


////#############################################################################
tensor EightNodeBrick::Nodal_Coordinates(void)
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
tensor EightNodeBrick::incr_disp(void)
  {
    const int dimensions[] = {8,3};  // Xiaoyan changed from {20,3} to {8,3} for 8 nodes
    tensor increment_disp(2, dimensions, 0.0);

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
    const Vector &IncrDis1 = nd1Ptr->getIncrDisp();
    const Vector &IncrDis2 = nd2Ptr->getIncrDisp();
    const Vector &IncrDis3 = nd3Ptr->getIncrDisp();
    const Vector &IncrDis4 = nd4Ptr->getIncrDisp();
    const Vector &IncrDis5 = nd5Ptr->getIncrDisp();
    const Vector &IncrDis6 = nd6Ptr->getIncrDisp();
    const Vector &IncrDis7 = nd7Ptr->getIncrDisp();
    const Vector &IncrDis8 = nd8Ptr->getIncrDisp();
    
    increment_disp.val(1,1)=IncrDis1(0); increment_disp.val(1,2)=IncrDis1(1);increment_disp.val(1,3)=IncrDis1(2);
    increment_disp.val(2,1)=IncrDis2(0); increment_disp.val(2,2)=IncrDis2(1);increment_disp.val(2,3)=IncrDis2(2);
    increment_disp.val(3,1)=IncrDis3(0); increment_disp.val(3,2)=IncrDis3(1);increment_disp.val(3,3)=IncrDis3(2);
    increment_disp.val(4,1)=IncrDis4(0); increment_disp.val(4,2)=IncrDis4(1);increment_disp.val(4,3)=IncrDis4(2);
    increment_disp.val(5,1)=IncrDis5(0); increment_disp.val(5,2)=IncrDis5(1);increment_disp.val(5,3)=IncrDis5(2);
    increment_disp.val(6,1)=IncrDis6(0); increment_disp.val(6,2)=IncrDis6(1);increment_disp.val(6,3)=IncrDis6(2);
    increment_disp.val(7,1)=IncrDis7(0); increment_disp.val(7,2)=IncrDis7(1);increment_disp.val(7,3)=IncrDis7(2);
    increment_disp.val(8,1)=IncrDis8(0); increment_disp.val(8,2)=IncrDis8(1);increment_disp.val(8,3)=IncrDis8(2);
    
    return increment_disp;
  }

////#############################################################################
tensor EightNodeBrick::total_disp(FILE *fp, double * u)
  {
    const int dimensions[] = {8,3};  // Xiaoyan changed from {20,3} to {8,3} for 8 nodes
    tensor total_disp(2, dimensions, 0.0);
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
    
    total_disp.val(1,1)=TotDis1(0); total_disp.val(1,2)=TotDis1(1);total_disp.val(1,3)=TotDis1(2);
    total_disp.val(2,1)=TotDis2(0); total_disp.val(2,2)=TotDis2(1);total_disp.val(2,3)=TotDis2(2);
    total_disp.val(3,1)=TotDis3(0); total_disp.val(3,2)=TotDis3(1);total_disp.val(3,3)=TotDis3(2);
    total_disp.val(4,1)=TotDis4(0); total_disp.val(4,2)=TotDis4(1);total_disp.val(4,3)=TotDis4(2);
    total_disp.val(5,1)=TotDis5(0); total_disp.val(5,2)=TotDis5(1);total_disp.val(5,3)=TotDis5(2);
    total_disp.val(6,1)=TotDis6(0); total_disp.val(6,2)=TotDis6(1);total_disp.val(6,3)=TotDis6(2);
    total_disp.val(7,1)=TotDis7(0); total_disp.val(7,2)=TotDis7(1);total_disp.val(7,3)=TotDis7(2);
    total_disp.val(8,1)=TotDis8(0); total_disp.val(8,2)=TotDis8(1);total_disp.val(8,3)=TotDis8(2);

    return total_disp;
  }

////#############################################################################
int EightNodeBrick::get_global_number_of_node(int local_node_number)
{
  //return G_N_numbs[local_node_number];	
  return connectedExternalNodes(local_node_number);
}

////#############################################################################
int  EightNodeBrick::get_Brick_Number(void)
{
  //return elem_numb;
  return this->getTag();
}

////#############################################################################
int * EightNodeBrick::get_LM(void)
  {
    return LM;
  }

//Commented out Zhaohui 09-27-2000

//////#############################################################################
//void EightNodeBrick::set_LM(Node * node)
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
tensor EightNodeBrick::nodal_forces(void)
  {
    int force_dim[] = {8,3};  // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor nodal_forces(2,force_dim,0.0);

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
                dhGlobal = dh("ij") * JacobianINV("jk");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;
                //..::printf("\n\nIN THE nodal forces ----**************** where = %d \n", where);
                //..::printf("                    GP_c_r = %d,  GP_c_s = %d,  GP_c_t = %d\n",
                //..                           GP_c_r,GP_c_s,GP_c_t);
                //..::printf("WEIGHT = %f", weight);
                //..::printf("determinant of Jacobian = %f", det_of_Jacobian);
                //..MatPoint[where].report("Gauss Point\n");

                //..   // samo jos odredi ovaj tensor E i to za svaku gauss tacku drugaciji !!!!!!!!!!!!
                //..   ovde negde bi trebalo da se na osnovu inkrementalnih pomeranja
                //..   nadje inkrementalna deformacija ( strain_increment ) pa sa njom dalje:
                //..
                //// tensor of incremental displacements taken from node objects
                //                incremental_displacements = incr_disp();
                //
                //// incremental straines at this Gauss point
                //                incremental_strain =
                //                  (dhGlobal("ib")*incremental_displacements("ia")).symmetrize11();
                //
                //                incremental_strain.null_indices();
                ////incremental_strain.reportshort("\n incremental_strain tensor at GAUSS point\n");
                //
                ////                integr_type = (MatPoint)->operator[](where).integration_type();
                ////                if ( !strcmp(integr_type,"BakcwardEuler")

                //..   dakle ovde posalji strain_increment jer se stari stress cuva u okviru svake
                //..   Gauss tacke a samo saljes strain_increment koji ce da se prenese
                //..   u integracionu rutinu pa ce ta da vrati krajnji napon i onda moze da
                //..   se pravi ConstitutiveStiffnessTensor.
                //.. Ustvari posalji sve sto imas ( incremental_strain, start_stress,
                //.. number_of_subincrements . . . u ovu Constitutive_tensor funkciju
                //.. pa ona nek ide, u zavisnosti od modela koji se koristi i neka
                //.. onda tamo u svakoj posebnoj modelskoj funkciji vrati sta treba
                //.. ( recimo Elastic odmah vraca Eelastic a recimo MRS_Lade prvo
                //.. pita koji nacin integracije da koristi pa onda u zvisnosti od toga
                //.. zove funkcuju koja integrali za taj algoritam ( ForwardEuler, BakcwardEuler,
                //.. SemiBackwardEuler, . . . ) i onda kada funkcija vrati napon onda
                //.. se opet pita koji je tip integracije bio u pitanju pa pravi odgovarajuci
                //.. ConstitutiveTensor i vraca ga nazad!

                //                   stress_at_GP = (GPstress)->operator[](where);
                //stress_at_GP = GPstress[where];

	        EPState *tmp_eps = (MatPoint[where]).getEPS();
	        NDMaterial *tmp_ndm = (MatPoint[where]).getNDMat();

		if ( tmp_eps ) {     //Elasto-plastic case
 		    stress_at_GP = (MatPoint[where].getEPS())->getStress();
		}
		else if ( tmp_ndm ) { //Elastic case
             	    stress_at_GP = (MatPoint[where].getNDMat())->getStressTensor();
		}
		else {
               	   g3ErrorHandler->fatal("EightNodeBrick::nodal_forces (tag: %d), could not getStress", this->getTag());
		   exit(1);
		}

                //..stress_at_GP.reportshort("\n stress_at_GPtensor at GAUSS point for nodal forces \n");

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forces =
                     nodal_forces + dhGlobal("ib")*stress_at_GP("ab")*weight;
                //::::::  nodal_forces.print("nf","\n\n Nodal Forces \n");
 
              }
          }
      }


    return nodal_forces;

  }

////#############################################################################
// returns nodal forces for given ITERATIVE stress field in an element
tensor EightNodeBrick::iterative_nodal_forces(void)
  {
    int force_dim[] = {8,3}; // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

    tensor nodal_forces(2,force_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {8,3};   // Xiaoyan changed from {20,3 to {8,3} for 8 nodes

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
                //.....::printf("EightNodeBrick::iterative_nodal_forces(void)  ----**************** where = %d \n", where);
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
                dhGlobal = dh("ij") * JacobianINV("jk");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

                //                   stress_at_GP = (GPstress)->operator[](where);
                //stress_at_GP = GPiterative_stress[where];
                
		//stress_at_GP = ( MatPoint[where].getTrialEPS() )->getStress();
                stress_at_GP = ( MatPoint[where].getEPS() )->getStress();
                stress_at_GP.reportshortpqtheta("\n iterative_stress at GAUSS point in iterative_nodal_force\n");

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forces =
                  nodal_forces + dhGlobal("ib")*stress_at_GP("ab")*weight;
                //nodal_forces.print("nf","\n EightNodeBrick::iterative_nodal_forces Nodal Forces ~~~~\n");

              }
          }
      }


    return nodal_forces;

  }

////#############################################################################
// returns nodal forces for given constant stress field in the element
tensor EightNodeBrick::nodal_forces_from_stress(stresstensor & stress)
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
                //.....::printf("EightNodeBrick::iterative_nodal_forces(void)  ----**************** where = %d \n", where);
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
                dhGlobal = dh("ij") * JacobianINV("jk");

                //weight
                weight = rw * sw * tw * det_of_Jacobian;

                //                   stress_at_GP = (GPstress)->operator[](where);
                //                stress_at_GP = GPiterative_stress[where];
                //GPiterative_stress[where].reportshortpqtheta("\n iterative_stress at GAUSS point in iterative_nodal_force\n");

                // nodal forces See Zienkievicz part 1 pp 108
                nodal_forces =
                  nodal_forces + dhGlobal("ib")*stress("ab")*weight;
                //nodal_forces.print("nf","\n EightNodeBrick::iterative_nodal_forces Nodal Forces ~~~~\n");

              }
          }
      }

    return nodal_forces;

  }


////#############################################################################
// returns nodal forces for given incremental strain field in an element
// by using the linearized constitutive tensor from the begining of the step !
tensor EightNodeBrick::linearized_nodal_forces(void)
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
    incremental_displacements = incr_disp();
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
                dhGlobal = dh("ij") * JacobianINV("jk");

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
                
	        EPState *tmp_eps = (MatPoint[where]).getEPS();
	        NDMaterial *tmp_ndm = (MatPoint[where]).getNDMat();

		if ( tmp_eps ) {     //Elasto-plastic case
		    mmodel->setEPS( *tmp_eps );
		    if ( ! (mmodel->setTrialStrainIncr( incremental_strain)) )
               	       g3ErrorHandler->warning("EightNodeBrick::incremental_Update (tag: %d), not converged",
					 this->getTag());
		    Constitutive = mmodel->getTangentTensor();
      		    MatPoint[where].setEPS( mmodel->getEPS() ); //Set the new EPState back
		}
		else if ( tmp_ndm ) { //Elastic case
		    (MatPoint[where].p_matmodel)->setTrialStrainIncr( incremental_strain );
		    Constitutive = (MatPoint[where].p_matmodel)->getTangentTensor();
		}
		else {
               	   g3ErrorHandler->fatal("EightNodeBrick::incremental_Update (tag: %d), could not getTangentTensor", this->getTag());
		   exit(1);
		}
		
		//Constitutive = ( MatPoint[where].getEPS() )->getEep();
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

//....////#############################################################################
//....// updates Gauss point stresses and strains from given displacements
//....void EightNodeBrick::update_stress_strain(tensor & displacementsT)
//....  {
//....//    int force_dim[] = {20,3};
//....//    tensor nodal_forces(2,force_dim,0.0);
//....
//....    double r  = 0.0;
//....    double rw = 0.0;
//....    double s  = 0.0;
//....    double sw = 0.0;
//....    double t  = 0.0;
//....    double tw = 0.0;
//....
//....    short where = 0;
//....    double weight = 0.0;
//....
//....    int dh_dim[] = {20,3};
//....    tensor dh(2, dh_dim, 0.0);
//....
//....    stresstensor stress_at_GP(0.0);
//....    straintensor strain_at_GP(0.0);
//....
//....    double det_of_Jacobian = 0.0;
//....
//....    tensor Jacobian;
//....    tensor JacobianINV;
//....    tensor dhGlobal;
//....
//....    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
//....      {
//....        r = get_Gauss_p_c( r_integration_order, GP_c_r );
//....        rw = get_Gauss_p_w( r_integration_order, GP_c_r );
//....
//....        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
//....          {
//....            s = get_Gauss_p_c( s_integration_order, GP_c_s );
//....            sw = get_Gauss_p_w( s_integration_order, GP_c_s );
//....
//....            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
//....              {
//....                t = get_Gauss_p_c( t_integration_order, GP_c_t );
//....                tw = get_Gauss_p_w( t_integration_order, GP_c_t );
//....
//....// this short routine is supposed to calculate position of
//....// Gauss point from 3D array of short's
//....                where =
//....                ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
//....
//....//........................................................
//....//........................................................
//....// interpolation functions
//....                tensor h = b3darray[0].interp_poli_at(r,s,t);
//....                ::printf("\n\n r = %f, s = %f, t = %f\n", r, s, t);
//....//  h.print("h");
//....
//....// displacements
//....//....   tensor disp_at_rst = h("i")*displacementsT("ia");
//....//....   disp_at_rst.print("disp");
//....
//....// derivatives of interpolation functions
//....                dh = dh_drst_at(r,s,t);
//....//                ::printf("\n\n r = %f, s = %f, t = %f\n", r, s, t);
//....//  dh.print("dh");
//....
//....                Jacobian = b3darray[0].Jacobian_3D(dh);
//....//                Jacobian.print("J");
//....
//....                JacobianINV = b3darray[0].Jacobian_3Dinv(dh);
//....//                JacobianINV.print("JINV");
//....
//....//                det_of_Jacobian = Jacobian.determinant();
//....//                ::printf("determinant of Jacobian is %f\n",Jacobian_determinant );
//....
//....// Derivatives of local coordinates multiplied with inverse of Jacobian (see Bathe p-202)
//....                dhGlobal = dh("ij") * JacobianINV("jk");
//....// straines
//....//  strain = (dh("ib")*displacements("ia")).symmetrize11();
//....                strain = (dhGlobal("ib")*displacementsT("ia")).symmetrize11();
//....//  straintensor strain = dh("ib")*displacements("ia");
//....                strain.reportshort("\n strain tensor\n");
//....                strain.null_indices();
//....
//....//                tensor E = mmElastic.ElasticStiffness();
//....
//....//stresses
//....                stress = E("ijkl") * strain("kl");
//....                stress.reportshort("\n\n stress tensor \n");
//....//...
//....//........................................................
//....//........................................................
//....//........................................................
//....//........................................................
//....//........................................................
//....//........................................................
//....//........................................................
//....
//....
//....              }
//....          }
//....      }
//....
//....  }

////#############################################################################
////#############################################################################
//double EightNodeBrick::get_first_q_ast(void)
//  {
//    double ret = MatPoint[0].kappa_cone_get();
//
//    return ret;
//
//  }
////#############################################################################
//double EightNodeBrick::get_first_etacone(void)
//  {
//    double ret = MatPoint[0].etacone();
//
//    return ret;
//
//  }
//

//#############################################################################
void EightNodeBrick::report(char * msg)
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
	   //  (nodes[G_N_numbs[in]]).Print(cout);

	   nd1Ptr->Print(cout);
	   nd2Ptr->Print(cout);
	   nd3Ptr->Print(cout);
	   nd4Ptr->Print(cout);
	   nd5Ptr->Print(cout);
	   nd6Ptr->Print(cout);
           nd7Ptr->Print(cout);
	   nd8Ptr->Print(cout);

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
                MatPoint[where].report("Integration Point\n");
                //GPstress[where].reportshort("stress at Gauss Point");
                //GPstrain[where].reportshort("strain at Gauss Point");
                //MatPoint[where].report("Material model  at Gauss Point");
              }
          }
      }

  }


//#############################################################################
void EightNodeBrick::reportshort(char * msg)
  {
    if ( msg ) ::printf("** %s",msg);
    ::printf("\n Element Number = %d\n", this->getTag() );
    ::printf("\n Number of nodes in a EightNodeBrick = %d\n",
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
void EightNodeBrick::reportPAK(char * msg)
  {
    if ( msg ) ::printf("%s",msg);
    ::printf("%10d   ",  this->getTag());
           for ( int i=0 ; i<8 ; i++ )
             ::printf( "%6d",connectedExternalNodes(i) );
             //::printf("%6d",G_N_numbs[i]);

    printf("\n");
  }


//#############################################################################
void EightNodeBrick::reportpqtheta(int GP_numb)
  {
    short where = GP_numb-1;
    MatPoint[where].reportpqtheta("");
  }

//#############################################################################
void EightNodeBrick::reportLM(char * msg) 
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
void EightNodeBrick::reportTensor(char * msg)
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
    tensor IntegrationPointCoord(2, dim, 0.0);
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
                  //  IntegrationPointCoord.val(1,where+1) =+NodalCoord.val(1,where+1) * H.val(encount*3-2,1);
                  //  IntegrationPointCoord.val(2,where+1) =+NodalCoord.val(2,where+1) * H.val(encount*3-1,2);
                  //  IntegrationPointCoord.val(3,where+1) =+NodalCoord.val(3,where+1) * H.val(encount*3-0,3);
                  IntegrationPointCoord.val(1,where+1) +=NodalCoord.val(1,encount) * H.val(encount*3-2,1);
                  //::printf("-- NO nodal, H_val :%d %+.2e %+.2e %+.5e\n", encount,NodalCoord.val(1,encount),H.val(encount*3-2,1),IntegrationPointCoord.val(1,where+1) );
                  IntegrationPointCoord.val(2,where+1) +=NodalCoord.val(2,encount) * H.val(encount*3-1,2);
                  IntegrationPointCoord.val(3,where+1) +=NodalCoord.val(3,encount) * H.val(encount*3-0,3);

		  }
		  			  
    ::printf("gauss point# %d   %+.6e %+.6e %+.6e \n", where+1, 
                                                       IntegrationPointCoord.val(1,where+1),
                                                       IntegrationPointCoord.val(2,where+1),
                                                       IntegrationPointCoord.val(3,where+1));

    //MatPoint[where].reportTensor("");

		
              }
          }
      }

 }


////#############################################################################

//#############################################################################
//void EightNodeBrick::reportTensor(char * msg)
// ZHaohui added to print gauss point coord. to file fp

void EightNodeBrick::reportTensorF(FILE * fp)
  {
//    if ( msg ) ::printf("** %s\n",msg);

// special case for 8 nodes only
// special case for 8 nodes only
// special case for 8 nodes only
// special case for 8 nodes only
// special case for 8 nodes only
    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;

    short where = 0;

    // special case for 8 nodes only
    static const int dim[] = {3, 8}; // static-> see ARM pp289-290
    tensor NodalCoord(2, dim, 0.0);
    tensor IntegrationPointCoord(2, dim, 0.0);
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
                  //  IntegrationPointCoord.val(1,where+1) =+NodalCoord.val(1,where+1) * H.val(encount*3-2,1);
                  //  IntegrationPointCoord.val(2,where+1) =+NodalCoord.val(2,where+1) * H.val(encount*3-1,2);
                  //  IntegrationPointCoord.val(3,where+1) =+NodalCoord.val(3,where+1) * H.val(encount*3-0,3);
                  IntegrationPointCoord.val(1,where+1) +=NodalCoord.val(1,encount) * H.val(encount*3-2,1);
                  //::printf("-- NO nodal, H_val :%d %+.2e %+.2e %+.5e\n", encount,NodalCoord.val(1,encount),H.val(encount*3-2,1),IntegrationPointCoord.val(1,where+1) );
                  IntegrationPointCoord.val(2,where+1) +=NodalCoord.val(2,encount) * H.val(encount*3-1,2);
                  IntegrationPointCoord.val(3,where+1) +=NodalCoord.val(3,encount) * H.val(encount*3-0,3);

	       }
		  			  
    fprintf(fp, "gauss point# %d   %+.6e %+.6e %+.6e \n", where+1, 
                                                          IntegrationPointCoord.val(1,where+1),
                                                          IntegrationPointCoord.val(2,where+1),
                                                          IntegrationPointCoord.val(3,where+1));

    //MatPoint[where].reportTensor("");

		
              }
          }
      }

 }

//=============================================================================
//  The following are come from FourNodeQuad.cc	 Xiaoyan 07/06/00
//  The following are come from FourNodeQuad.cc	 Xiaoyan 07/06/00
//  The following are come from FourNodeQuad.cc	 Xiaoyan 07/06/00
//=============================================================================


//=============================================================================
int EightNodeBrick::getNumExternalNodes () const
{
    return 8;  //changed from 4 to 8 Xiaoyan 07/06/00
}


//=============================================================================
const ID& EightNodeBrick::getExternalNodes ()
{
    return connectedExternalNodes;
}

//=============================================================================
int EightNodeBrick::getNumDOF ()
{
    return 24;	     //Changed from 2*4=8 to 3*8=24 Xiaoyan 07/06/00
}

//=============================================================================
void EightNodeBrick::setDomain (Domain *theDomain)
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

	g3ErrorHandler->fatal("FATAL ERROR EightNodeBrick (tag: %d), node not found in domain",
		this->getTag());
	
	return;
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
								      
    if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3 ||  // Changed 2 to 3 Xiaoyan
        dofNd5 != 3 || dofNd6 != 3 || dofNd7 != 3 || dofNd8 != 3 ) {
	g3ErrorHandler->fatal("FATAL ERROR EightNodeBrick (tag: %d), has differing number of DOFs at its nodes",
 		this->getTag());
	
	return;
    }
    this->DomainComponent::setDomain(theDomain);
}

//=============================================================================
int EightNodeBrick::commitState ()
{
    // int order = theQuadRule->getOrder();     // Commented by Xiaoyan

    int i;
    //int j, k;      // Xiaoyan added k for three dimension		       
    int retVal = 0;

    // Loop over the integration points and commit the material states
    int count  = r_integration_order* s_integration_order * t_integration_order;
    //for (i = 0; i < r_integration_order; i++)		    // Xiaoyan chaneged order to
    //  for (j = 0; j < s_integration_order; j++)	    // r_integration_order,
    //							    // s_integration_order, and
    //	    for (k = 0; k < t_integration_order; k++)	    // added t_integration_order,
    //         retVal += (GaussPtheMaterial[i][j][k]).commitState();  

    EPState *tmp_eps;
    NDMaterial *tmp_ndm; 
    
    for (i = 0; i < count; i++)	{
    
       tmp_eps = (MatPoint[i]).getEPS();
       tmp_ndm = (MatPoint[i]).getNDMat();
       
       if ( tmp_eps ) {     //Elasto-plastic case
           retVal += (MatPoint[i].gpEPS)->commitState(); //commit the State vars
       }
       else if ( tmp_ndm ) { //Elastic case
           retVal += (MatPoint[i].p_matmodel)->commitState();
       }
       else {
          g3ErrorHandler->warning("EightNodeBrick::commit (tag: %d), No state params", this->getTag());
          //exit(1);
       }
    }
     
    return retVal;
}

//=============================================================================
int EightNodeBrick::revertToLastCommit ()
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

    EPState *tmp_eps;
    NDMaterial *tmp_ndm; 
    
    for (i = 0; i < count; i++)	{
    
       tmp_eps = (MatPoint[i]).getEPS();
       tmp_ndm = (MatPoint[i]).getNDMat();
       
       if ( tmp_eps ) {     //Elasto-plastic case
           retVal += (MatPoint[i].gpEPS)->revertToLastCommit(); 
       }
       else if ( tmp_ndm ) { //Elastic case
           retVal += (MatPoint[i].p_matmodel)->revertToLastCommit();
       }
       else {
          g3ErrorHandler->warning("EightNodeBrick::revertToLastCommit (tag: %d), No state params", this->getTag());
          //exit(1);
       }
    }

    return retVal;
}

//=============================================================================
int EightNodeBrick::revertToStart () 
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
    	     
    EPState *tmp_eps;
    NDMaterial *tmp_ndm; 
    
    for (i = 0; i < count; i++)	{
    
       tmp_eps = (MatPoint[i]).getEPS();
       tmp_ndm = (MatPoint[i]).getNDMat();
       
       if ( tmp_eps ) {     //Elasto-plastic case
           retVal += (MatPoint[i].gpEPS)->revertToStart(); //commit the State vars
       }
       else if ( tmp_ndm ) { //Elastic case
           retVal += (MatPoint[i].p_matmodel)->revertToStart();
       }
       else {
          g3ErrorHandler->warning("EightNodeBrick::revertToStart (tag: %d), no state parameters", this->getTag());
          //exit(1);
       }
    }
    
    return retVal;

    // Loop over the integration points and revert to initial material states
   } 


//=============================================================================
const Matrix &EightNodeBrick::getTangentStiff () 
{ 
     tensor stifftensor = getStiffnessTensor();
     tensor stiffmatrix = stiffness_matrix ( stifftensor );
     for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 24; j++) {
	    K(i, j) = stiffmatrix.val(i+1, j+1);
	}
     }
     return K;
}

//=============================================================================
const Matrix &EightNodeBrick::getSecantStiff () 
{
     return K;
}

//=============================================================================
const Matrix &EightNodeBrick::getDamp () 
{    
     return C;
}

//=============================================================================
const Matrix &EightNodeBrick::getMass () 
{
     return M;
}

//=============================================================================
void EightNodeBrick::zeroLoad(void)
{
     Q.Zero();
}


//=============================================================================
int  EightNodeBrick::addLoad(const Vector &addLoad)
{
     if (addLoad.Size() != 24) {
     	g3ErrorHandler->warning("EightNodeBrick::addLoad %s\n",
     			"Vector not of correct size");
     	return -1;
     }

     // Add to the external nodal loads
     Q += addLoad;

     return 0;
}

//=============================================================================
int EightNodeBrick::addInertiaLoadToUnbalance(const Vector &accel)
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

    if (3 != Raccel1.Size() || 3 != Raccel2.Size() || 3 != Raccel3.Size() ||
       	3 != Raccel4.Size() || 3 != Raccel5.Size() || 3 != Raccel6.Size() || 
	3 != Raccel7.Size() || 3 != Raccel8.Size()  ){	
	// Xiaoyan changed 2 to 3 and added Racce15-18  09/27/00
		g3ErrorHandler->warning("FourNodeQuad::addInertiaLoadToUnbalance %s\n",
				"matrix and vector sizes are incompatable");
		return -1;
    }

	static Vector ra(24);  // Changed form 8 to 24(3*8)  Xiaoyan 09/27/00

	ra(0) =  Raccel1(0);
	ra(1) =  Raccel1(1);
	ra(2) =  Raccel1(2);
	ra(3) =  Raccel2(0);
	ra(4) =  Raccel2(1);
	ra(5) =  Raccel2(2);
	ra(6) =  Raccel3(0);
	ra(7) =  Raccel3(1);
	ra(8) =  Raccel3(2);
	ra(9) =  Raccel4(0);
	ra(10) = Raccel4(1);
	ra(11) = Raccel4(2);
    	// Xiaoyan added the following 09/27/00
    	ra(12) = Raccel5(0);
	ra(13) = Raccel5(1);
	ra(14) = Raccel5(2);
	ra(15) = Raccel6(0);
	ra(16) = Raccel6(1);
	ra(17) = Raccel6(2);
	ra(18) = Raccel7(0);
	ra(19) = Raccel7(1);
	ra(20) = Raccel7(2);
	ra(21) = Raccel8(0);
	ra(22) = Raccel8(1);
	ra(23) = Raccel8(2);

    // Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	// Mass matrix is computed in setDomain()
    for (int i = 0; i < 24; i++)   // Changed form 8 to 24  Xiaoyan 09/27/00
		Q(i) += -M(i,i)*ra(i);

    return 0;
}



//=============================================================================
const Vector &EightNodeBrick::getResistingForce () 
{     
     return P;
}

//=============================================================================
const Vector &EightNodeBrick::getResistingForceIncInertia () 
{
	// Check for a quick return
	if (rho == 0.0)
		return this->getResistingForce();

	const Vector &accel1 = nd1Ptr->getTrialAccel();
	const Vector &accel2 = nd2Ptr->getTrialAccel();
	const Vector &accel3 = nd3Ptr->getTrialAccel();
	const Vector &accel4 = nd4Ptr->getTrialAccel();
        // Xiaoyan added the following four 09/27/00	
	const Vector &accel5 = nd5Ptr->getTrialAccel();
	const Vector &accel6 = nd6Ptr->getTrialAccel();
	const Vector &accel7 = nd7Ptr->getTrialAccel();
	const Vector &accel8 = nd8Ptr->getTrialAccel();
	
	static Vector a(24);  // originally 8

	a(0) =  accel1(0);
	a(1) =  accel1(1);
	a(2) =  accel1(2);
	a(3) =  accel2(0);
	a(4) =  accel2(1);
	a(5) =  accel2(2);
	a(6) =  accel3(0);
	a(7) =  accel3(1);
	a(8) =  accel3(2);
	a(9) =  accel4(0);
	a(10) = accel4(1);
	a(11) = accel4(2);
    	// Xiaoyn added the following 09/27/00
    	a(12) = accel5(0);
	a(13) = accel5(1);
	a(14) = accel5(2);
	a(15) = accel6(0);
	a(16) = accel6(1);
	a(17) = accel6(2);
	a(18) = accel7(0);
	a(19) = accel7(1);
	a(20) = accel7(2);
	a(21) = accel8(0);
	a(22) = accel8(1);
	a(23) = accel8(2);
		   		   
	// Compute the current resisting force
	this->getResistingForce();

	// Take advantage of lumped mass matrix
	// Mass matrix is computed in setDomain()
	for (int i = 0; i < 24; i++)
		P(i) += M(i,i)*a(i);

	return P;
} 

//=============================================================================
int EightNodeBrick::sendSelf (int commitTag, Channel &theChannel) 
{ 
     // Not implemtented yet
     return 0;
}

//=============================================================================
int EightNodeBrick::recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) 
{   
     // Not implemtented yet
     return 0;
}


//=============================================================================
int EightNodeBrick::displaySelf (Renderer &theViewer, int displayMode, float fact)
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
void EightNodeBrick::Print(ostream &s, int flag =0)
{
     // Not implemtented yet
     //return 0;
}

//=============================================================================
int EightNodeBrick::setResponse (char **argv, int argc, Information &eleInformation) 
{
     return 0;
}
//=============================================================================

int EightNodeBrick::getResponse (int responseID, Information &eleInformation)
{
     return 0;
}
//=============================================================================


//const Matrix&
//EightNodeBrick::getTangentStiff ()
//{
//	int order = theQuadRule->getOrder();
//	const Vector &intPt = theQuadRule->getIntegrPointCoords();
//	const Vector &intWt = theQuadRule->getIntegrPointWeights();
//
//	const Vector &disp1 = nd1Ptr->getTrialDisp();
//        const Vector &disp2 = nd2Ptr->getTrialDisp();
//	const Vector &disp3 = nd3Ptr->getTrialDisp();
//        const Vector &disp4 = nd4Ptr->getTrialDisp();
//       // Xiaoyan added 5-8 07/06/00
//        const Vector &disp5 = nd5Ptr->getTrialDisp();
//        const Vector &disp6 = nd6Ptr->getTrialDisp();
//	const Vector &disp7 = nd7Ptr->getTrialDisp();
//        const Vector &disp8 = nd8Ptr->getTrialDisp();
//
//	static Vector u(24);	    //Changed from u(8) to u(24) Xiaoyn 07/06/00
//
//	u(0) = disp1(0);
//	u(1) = disp1(1);
//        u(2) = disp1(2);
//	u(3) = disp2(0);
//	u(4) = disp2(1);
//	u(5) = disp2(2);
//        u(6) = disp3(0);	     
//	u(7) = disp3(1);
//	u(8) = disp3(2);
//	u(9) = disp4(0);
//	u(10) = disp4(1);
//	u(11) = disp4(2);
//	u(12) = disp5(0);
//	u(13) = disp5(1);
//	u(14) = disp5(2);
//	u(15) = disp6(0);
//	u(16) = disp6(1);
//	u(17) = disp6(2);
//	u(18) = disp7(0);
//	u(19) = disp7(1);
//	u(20) = disp7(2);
//	u(21) = disp8(0);
//	u(22) = disp8(1);
//	u(23) = disp8(2);


//	static Vector eps (6);		  // Changed eps(3) to eps(6) Xiaoyan 07/06/00

//	K.Zero();

//	// Loop over the integration points
//	for (int i = 0; i < order; i++)
//	{
//		for (int j = 0; j < order; j++)
//		{
//
//			// Determine Jacobian for this integration point
//			this->setJacobian (intPt(i), intPt(j));
//
//			// Interpolate strains
//			this->formBMatrix (intPt(i), intPt(j));
//			eps = B*u;
//
//			// Set the material strain
//			(theMaterial[i][j])->setTrialStrain (eps);
//
//			// Get the material tangent
//			const Matrix &D = (theMaterial[i][j])->getTangent();
//
//			// Form the Jacobian of the coordinate transformation
//			double detJ = this->formDetJ (intPt(i), intPt(j));
//
//			// Perform numerical integration
//			K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
//		}
//	}
//
//	K = K * thickness;
//
//	return K;
//}

//const Matrix&
//EightNodeBrick::getSecantStiff ()
//{
//	return K;
//}

//Commented by Xiaoyan     Use the form like Brick3d
//const Matrix & EightNodeBrick::getDamp ()
//{
//	return C;
//}
// Commented by Xiaoyan 08/04/00

//const Matrix&
//EightNodeBrick::getMass ()
//{
//	int order = theQuadRule->getOrder();
//	const Vector &intPt = theQuadRule->getIntegrPointCoords();
//	const Vector &intWt = theQuadRule->getIntegrPointWeights();
//
//	M.Zero();
//
//	int i, j;
//
//	// Loop over the integration points
//	for (i = 0; i < order; i++)
//	{
//		for (j = 0; j < order; j++)
//		{
//			// Determine Jacobian for this integration point
//			this->setJacobian (intPt(i), intPt(j));
//
//			// Interpolate strains
//			this->formNMatrix (intPt(i), intPt(j));
//
//			// Form the Jacobian of the coordinate transformation
//			double detJ = this->formDetJ (intPt(i), intPt(j));
//
//			// Perform numerical integration
//			M = M + (N^ N) * intWt(i)*intWt(j) * detJ;
//		}
//	}
//
//	M = M * thickness * rho;
//
//	// Lumped mass ... can be optional
//	for (i = 0; i < 24; i++)	     // Changed 8 to 24  Xiaoyan 07/06/00
//	{
//		double sum = 0.0;
//		for (j = 0; j < 24; j++)    // Changed 8 to 24  Xiaoyan 07/06/00
//		{			    
//			sum += M(i,j);
//			M(i,j) = 0.0;
//		}
//		M(i,i) = sum;
//	}
//
//	return M;
//}
//
//const Vector&
//EightNodeBrick::getResistingForce ()
//{
//	int order = theQuadRule->getOrder();
//	const Vector &intPt = theQuadRule->getIntegrPointCoords();
//	const Vector &intWt = theQuadRule->getIntegrPointWeights();
//	
//	const Vector &disp1 = nd1Ptr->getTrialDisp();
//        const Vector &disp2 = nd2Ptr->getTrialDisp();
//	const Vector &disp3 = nd3Ptr->getTrialDisp();
//        const Vector &disp4 = nd4Ptr->getTrialDisp();
//	//6-8 added by Xiaoyan 07/06/00
//	const Vector &disp5 = nd5Ptr->getTrialDisp();
//        const Vector &disp6 = nd6Ptr->getTrialDisp();
//	const Vector &disp7 = nd7Ptr->getTrialDisp();
//        const Vector &disp8 = nd8Ptr->getTrialDisp();
//
//
//	static Vector u(24);	    //Changed from u(8) to u(24) Xiaoyn 07/06/00
//
//	u(0) = disp1(0);
//	u(1) = disp1(1);
//        u(2) = disp1(2);
//	u(3) = disp2(0);
//	u(4) = disp2(1);
//	u(5) = disp2(2);
//        u(6) = disp3(0);	     
//	u(7) = disp3(1);
//	u(8) = disp3(2);
//	u(9) = disp4(0);
//	u(10) = disp4(1);
//	u(11) = disp4(2);
//	u(12) = disp5(0);
//	u(13) = disp5(1);
//	u(14) = disp5(2);
//	u(15) = disp6(0);
//	u(16) = disp6(1);
//	u(17) = disp6(2);
//	u(18) = disp7(0);
//	u(19) = disp7(1);
//	u(20) = disp7(2);
//	u(21) = disp8(0);
//	u(22) = disp8(1);
//	u(23) = disp8(2);
//
//	eps (6);      //Changed eps(3) to eps(6) Xiaoyan 07/06/00
//
//	P.Zero();
//
//	// Loop over the integration points
//	for (int i = 0; i < order; i++)
//	{
//		for (int j = 0; j < order; j++)
//		{
//			// Determine Jacobian for this integration point
//			this->setJacobian (intPt(i), intPt(j));
//
//			// Interpolate strains
//			this->formBMatrix (intPt(i), intPt(j));
//			eps = B*u;
//
//			// Set the material strain
//			(theMaterial[i][j])->setTrialStrain (eps);
//
//			// Get material stress response
//			const Vector &sigma = (theMaterial[i][j])->getStress();
//
//			// Form the Jacobian of the coordinate transformation
//			double detJ = this->formDetJ (intPt(i), intPt(j));
//
//			// Perform numerical integration
//			P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
//		}
//	}
//
//	P = P * thickness * -1;
//
//	return P;
//}
//
//const Vector&
//EightNodeBrick::getResistingForceIncInertia ()
//{
//	// Yet to implement
//	return P;
//}
//
//
//
//void
//EightNodeBrick::Print (ostream &s, int flag)
//{
//	s << "EightNodeBrick, element id:  " << this->getTag() << endl;
//	s << "Connected external nodes:  " << connectedExternalNodes;
//	s << "Material model:  " << theMaterial[0][0]->getType() << endl;
//	s << "Element thickness:  " << thickness << endl;
//	s << "Element mass density:  " << rho << endl << endl;
//}
//
//
//int
//EightNodeBrick::displaySelf (Renderer &theViewer, int displayMode, float fact)
//{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
//        const Vector &end1Crd = nd1Ptr->getCrds();
//        const Vector &end2Crd = nd2Ptr->getCrds();	
//	const Vector &end3Crd = nd3Ptr->getCrds();	
//	const Vector &end4Crd = nd4Ptr->getCrds();	
//	// 5-8 were added by Xiaoyan
//        const Vector &end5Crd = nd5Ptr->getCrds();
//        const Vector &end6Crd = nd6Ptr->getCrds();	
//	const Vector &end7Crd = nd7Ptr->getCrds();	
//	const Vector &end8Crd = nd8Ptr->getCrds();	
////---------------------------------------------------------------
//    	const Vector &end1Disp = nd1Ptr->getDisp();
//	const Vector &end2Disp = nd2Ptr->getDisp();
//	const Vector &end3Disp = nd3Ptr->getDisp();
//	const Vector &end4Disp = nd4Ptr->getDisp();
//
	// 5-8 were added by Xiaoyan
//        const Vector &end5Disp = nd5Ptr->getDisp();
//	const Vector &end6Disp = nd6Ptr->getDisp();
//	const Vector &end7Disp = nd7Ptr->getDisp();
//	const Vector &end8Disp = nd8Ptr->getDisp();
//
//	Vector v1(3);
//	Vector v2(3);
//	Vector v3(3);
//	Vector v4(3);
//	//5-8 added by Xiaoyan 07/06/00
//	Vector v5(3);
//	Vector v6(3);
//	Vector v7(3);
//	Vector v8(3);
//
//	for (int i = 0; i < 3; i++)	    //Changed from i<2 to i<3, Xiaonyan 07/06/00
//	{
//		v1(i) = end1Crd(i) + end1Disp(i)*fact;
//		v2(i) = end2Crd(i) + end2Disp(i)*fact;    
//		v3(i) = end3Crd(i) + end3Disp(i)*fact;    
//		v4(i) = end4Crd(i) + end4Disp(i)*fact; 
//	
//		//5-8 added by Xiaoyan 07/06/00
//   		v5(i) = end5Crd(i) + end1Disp(i)*fact;
//		v6(i) = end6Crd(i) + end2Disp(i)*fact;    
//		v7(i) = end7Crd(i) + end3Disp(i)*fact;    
//		v8(i) = end8Crd(i) + end4Disp(i)*fact;    
//	}
//	int error = 0;
//	
//	error += theViewer.drawLine (v1, v2, 1.0, 1.0);
//	error += theViewer.drawLine (v2, v3, 1.0, 1.0);
//	error += theViewer.drawLine (v3, v4, 1.0, 1.0);
//	error += theViewer.drawLine (v4, v5, 1.0, 1.0);   // 5-8 added by Xiaoyan 07/06/00
//	error += theViewer.drawLine (v5, v6, 1.0, 1.0);
//	error += theViewer.drawLine (v6, v7, 1.0, 1.0);
//	error += theViewer.drawLine (v7, v8, 1.0, 1.0);
//	error += theViewer.drawLine (v8, v1, 1.0, 1.0);
//
//	return error;
//}
// The following are all commented by  Xiaoyan. We use the Brick3D to form these

//
//void
//EightNodeBrick::formNMatrix (double r, double s,double t) 
////Changed xi, eta to r,s and added t Xiaoyan  07/06/00
//{
//	N.Zero();
//
////	N(0,0) = N(1,1) = 0.25*(1.0-xi)*(1.0-eta);		// N_1
//// 	N(0,2) = N(1,3) = 0.25*(1.0+xi)*(1.0-eta);		// N_2
////	N(0,4) = N(1,5) = 0.25*(1.0+xi)*(1.0+eta);		// N_3
////	N(0,6) = N(1,7) = 0.25*(1.0-xi)*(1.0+eta);		// N_4
//
////	Changed by Xiaoyan 07/06/00
// The shape functions have been changed from N(2,8) to N(3,24)
// I take the node order according to Bathe's book p344-345. Xiaoyan
//        N(0,0)=N(1,1)=N(2,2)=1/8.*(1.0+r)*(1.0+s)*(1.0+t);
//	N(0,3)=N(1,4)=N(2,5)=1/8.*(1.0-r)*(1.0+s)*(1.0+t);
//	N(0,6)=N(1,7)=N(2,8)=1/8.*(1.0-r)*(1.0-s)*(1.0+t);
//	N(0,9)=N(1,10)=N(2,11)=1/8.*(1.0+r)*(1.0-s)*(1.0+t);
//	N(0,12)=N(1,13)=N(2,14)=1/8.*(1.0+r)*(1.0+s)*(1.0-t);
//	N(0,15)=N(1,16)=N(2,17)=1/8.*(1.0-r)*(1.0+s)*(1.0-t);
//	N(0,18)=N(1,19)=N(2,20)=1/8.*(1.0-r)*(1.0-s)*(1.0-t);
//	N(0,21)=N(1,22)=N(2,23)=1/8.*(1.0+r)*(1.0-s)*(1.0-t);
// }
//void
//EightNodeBrick::setJacobian (double r, double s, double t)
////Changed xi, eta to r,s and added t Xiaoyan 07/06/00
//{
//	const Vector &nd1Crds = nd1Ptr->getCrds();
//	const Vector &nd2Crds = nd2Ptr->getCrds();
//	const Vector &nd3Crds = nd3Ptr->getCrds();
//	const Vector &nd4Crds = nd4Ptr->getCrds();
//	// Xiaoyan added 5-8 07/06/00
//	const Vector &nd5Crds = nd5Ptr->getCrds();
//	const Vector &nd6Crds = nd6Ptr->getCrds();
//	const Vector &nd7Crds = nd7Ptr->getCrds();
//	const Vector &nd8Crds = nd8Ptr->getCrds();
//
////	J(0,0) = -nd1Crds(0)*(1.0-eta) + nd2Crds(0)*(1.0-eta) +
////				nd3Crds(0)*(1.0+eta) - nd4Crds(0)*(1.0+eta);
////
//	J(0,1) = -nd1Crds(0)*(1.0-xi) - nd2Crds(0)*(1.0+xi) +
//				nd3Crds(0)*(1.0+xi) + nd4Crds(0)*(1.0-xi);
//
//	J(1,0) = -nd1Crds(1)*(1.0-eta) + nd2Crds(1)*(1.0-eta) +
//				nd3Crds(1)*(1.0+eta) - nd4Crds(1)*(1.0+eta);
//
//	J(1,1) = -nd1Crds(1)*(1.0-xi) - nd2Crds(1)*(1.0+xi) +
//				nd3Crds(1)*(1.0+xi) + nd4Crds(1)*(1.0-xi);
//    	J = J * 0.25;
//
//	// For 3D problem Jacobi Matrix changed from J(2,2) to J(3,3)
//	// Xiaoyan  changed 07/06/00
//
//
//	J(0,0) = nd1Crds(0)*(1.0+s)*(1.0+t) - nd2Crds(0)*(1.0+s)*(1.0+t) -              
//		 nd3Crds(0)*(1.0-s)*(1.0+t) + nd4Crds(0)*(1.0-s)*(1.0+t) +
//		 nd5Crds(0)*(1.0+s)*(1.0-t) - nd6Crds(0)*(1.0+s)*(1.0-t) -              
//		 nd7Crds(0)*(1.0-s)*(1.0-t) + nd8Crds(0)*(1.0-s)*(1.0-t);
//	
//	J(0,1) = nd1Crds(1)*(1.0+s)*(1.0+t) - nd2Crds(1)*(1.0+s)*(1.0+t) -              
//		 nd3Crds(1)*(1.0-s)*(1.0+t) + nd4Crds(1)*(1.0-s)*(1.0+t) +
//		 nd5Crds(1)*(1.0+s)*(1.0-t) - nd6Crds(1)*(1.0+s)*(1.0-t) -              
//		 nd7Crds(1)*(1.0-s)*(1.0-t) + nd8Crds(1)*(1.0-s)*(1.0-t);
//	
//	J(0,2) = nd1Crds(2)*(1.0+s)*(1.0+t) - nd2Crds(2)*(1.0+s)*(1.0+t) -              
//		 nd3Crds(2)*(1.0-s)*(1.0+t) + nd4Crds(2)*(1.0-s)*(1.0+t) +
//		 nd5Crds(2)*(1.0+s)*(1.0-t) - nd6Crds(2)*(1.0+s)*(1.0-t) -              
//		 nd7Crds(2)*(1.0-s)*(1.0-t) + nd8Crds(2)*(1.0-s)*(1.0-t);
//                       
//	J(1,0) = nd1Crds(0)*(1.0+r)*(1.0+t) + nd2Crds(0)*(1.0-r)*(1.0+t) -              
//		 nd3Crds(0)*(1.0-r)*(1.0+t) - nd4Crds(0)*(1.0+r)*(1.0+t) +
//		 nd5Crds(0)*(1.0+r)*(1.0-t) + nd6Crds(0)*(1.0-r)*(1.0-t) -              
//		 nd7Crds(0)*(1.0-r)*(1.0-t) - nd8Crds(0)*(1.0+r)*(1.0-t);
//                                                               
//	J(1,1) = nd1Crds(1)*(1.0+r)*(1.0+t) + nd2Crds(1)*(1.0-r)*(1.0+t) -              
//		 nd3Crds(1)*(1.0-r)*(1.0+t) - nd4Crds(1)*(1.0+r)*(1.0+t) +
//		 nd5Crds(1)*(1.0+r)*(1.0-t) + nd6Crds(1)*(1.0-r)*(1.0-t) -              
//		 nd7Crds(1)*(1.0-r)*(1.0-t) - nd8Crds(1)*(1.0+r)*(1.0-t);
//       
//        J(1,2) = nd1Crds(2)*(1.0+r)*(1.0+t) + nd2Crds(2)*(1.0-r)*(1.0+t) -              
//		 nd3Crds(2)*(1.0-r)*(1.0+t) - nd4Crds(2)*(1.0+r)*(1.0+t) +
//		 nd5Crds(2)*(1.0+r)*(1.0-t) + nd6Crds(2)*(1.0-r)*(1.0-t) -              
//		 nd7Crds(2)*(1.0-r)*(1.0-t) - nd8Crds(2)*(1.0+r)*(1.0-t);
//
//	J(2,0) = nd1Crds(0)*(1.0+r)*(1.0+s) + nd2Crds(0)*(1.0-r)*(1.0+s) +              
//		 nd3Crds(0)*(1.0-r)*(1.0-s) + nd4Crds(0)*(1.0+r)*(1.0-s) -
//		 nd5Crds(0)*(1.0+r)*(1.0+s) - nd6Crds(0)*(1.0-r)*(1.0+s) -              
//		 nd7Crds(0)*(1.0-r)*(1.0-s) - nd8Crds(0)*(1.0+r)*(1.0-s);
//	 	
//	J(2,1) = nd1Crds(1)*(1.0+r)*(1.0+s) + nd2Crds(1)*(1.0-r)*(1.0+s) +              
//		 nd3Crds(1)*(1.0-r)*(1.0-s) + nd4Crds(1)*(1.0+r)*(1.0-s) -
//		 nd5Crds(1)*(1.0+r)*(1.0+s) - nd6Crds(1)*(1.0-r)*(1.0+s) -              
//		 nd7Crds(1)*(1.0-r)*(1.0-s) - nd8Crds(1)*(1.0+r)*(1.0-s);
//
//	J(2,2) = nd1Crds(2)*(1.0+r)*(1.0+s) + nd2Crds(2)*(1.0-r)*(1.0+s) +              
//		 nd3Crds(2)*(1.0-r)*(1.0-s) + nd4Crds(2)*(1.0+r)*(1.0-s) -
//		 nd5Crds(2)*(1.0+r)*(1.0+s) - nd6Crds(2)*(1.0-r)*(1.0+s) -              
//		 nd7Crds(2)*(1.0-r)*(1.0-s) - nd8Crds(2)*(1.0+r)*(1.0-s);
//				            
//	 J=J*0.125
//
//		// L = inv(J)  Changed from L(2,2) to L(3,3)  07/07/00
//
//	L(0,0)=-J(1,2)*J(2,1) + J(1,1)*J(2,2);
//	L(0.1)= J(0,2)*J(2,1) - J(0,1)*J(2,2);
//	L(0,3)=-J(0,2)*J(1,1) + J(0,1)*J(1,2);
//	L(1,0)= J(1,2)*J(2,0) - J(1,0)*J(2,2);
//	L(1,1)=-J(0,2)*J(2,0) + J(0,0)*J(2.2);
//	L(1,2)= J(0,2)*J(1,0) - J(0,0)*J(1,2);
//	L(2,0)=-J(1,1)*J(2,0) + J(1,0)*J(2,1);
//	L(2,1)= J(0,1)*J(2,0) - J(0,0)*J(2,1);
//	L(2,2)=-J(0,1)*J(1,0) + J(0,0)*J(1,1);
//	L=L/formDetJ(r,s,t)
//	    
//	L(0,0) = J(1,1);
//	L(1,0) = -J(0,1);
//	L(0,1) = -J(1,0);
//	L(1,1) = J(0,0);

//	L = L / formDetJ (xi, eta);
//}
//
//void
//EightNodeBrick::formBMatrix (double r, double s, double t)
////Changed xi, eta to r,s and added t Xiaoyan  07/06/00
//{
//    B.Zero();
//
//    //Changed by Xiaoyan 07/06/00
//    double L00 = L(0,0);
//    double L01 = L(0,1);
//    double L02 = L(0,1);
//    double L10 = L(1,0);
//    double L11 = L(1,1);
//    double L12 = L(1,2);
//    double L20 = L(2,0);
//    double L21 = L(2,1);
//    double L22 = L(2,2);
//
//    // See Cook, Malkus, Plesha p. 169 for the derivation of these terms
//    B(0,0) = L00*-0.25*(1.0-eta) + L01*-0.25*(1.0-xi);		// N_1,1
//    B(0,2) = L00*0.25*(1.0-eta) + L01*-0.25*(1.0+xi);		// N_2,1
//    B(0,4) = L00*0.25*(1.0+eta) + L01*0.25*(1.0+xi);		// N_3,1
//    B(0,6) = L00*-0.25*(1.0+eta) + L01*0.25*(1.0-xi);		// N_4,1
//	
//    B(1,1) = L10*-0.25*(1.0-eta) + L11*-0.25*(1.0-xi);  	// N_1,2
//    B(1,3) = L10*0.25*(1.0-eta) + L11*-0.25*(1.0+xi);		// N_2,2
//    B(1,5) = L10*0.25*(1.0+eta) + L11*0.25*(1.0+xi);		// N_3,2
//    B(1,7) = L10*-0.25*(1.0+eta) + L11*0.25*(1.0-xi);		// N_4,2
//
//    B(2,0) = B(1,1);
//    B(2,1) = B(0,0);
//    B(2,2) = B(1,3);
//    B(2,3) = B(0,2);
//    B(2,4) = B(1,5);
//    B(2,5) = B(0,4);
//    B(2,6) = B(1,7);
//    B(2,7) = B(0,6);
//}
//
//
//
////The derivative  of shape function to r,s,t respectly.
//// For example dh1dr means dh1/dr etc. Xiaoyan  07/07/00
//double  dh1dr=0.125*(1+s)*(1+t);
//double  dh1ds=0.125*(1+r)*(1+t);
//double  dh1dt=0.125*(1+r)*(1+s);
//
//double  dh2dr=-0.125*(1+s)*(1+t);
//double  dh2ds=0.125*(1-r)*(1+t);
//double  dh2dt=0.125*(1-r)*(1+s);
//
//double  dh3dr=-0.125*(1-s)*(1+t);
//double  dh3ds=-0.125*(1-r)*(1+t);
//double  dh3dt=0.125*(1-r)*(1-s);
//
//double  dh4dr=0.125*(1-s)*(1+t);
//double  dh4ds=-0.125*(1+r)*(1+t);
//double  dh4dt=0.125*(1+r)*(1-s);
//
//double  dh5dr=0.125*(1+s)*(1-t);
//double  dh5ds=0.125*(1+r)*(1-t);
//double  dh5dt=-0.125*(1+r)*(1+s);
//
//double  dh6dr=-0.125*(1+s)*(1-t);
//double  dh6ds=0.125*(1-r)*(1-t);
//double  dh6dt=-0.125*(1-r)*(1+s);
//
//double  dh7dr=-0.125*(1-s)*(1-t);
//double  dh7ds=-0.125*(1-r)*(1-t);
//double  dh7dt=-0.125*(1-r)*(1-s);
//
//double  dh8dr=0.125*(1-s)*(1-t);
//double  dh8ds=-0.125*(1+r)*(1-t);
//double  dh8dt=-0.125*(1+r)*(1-s);
//
//// B Matrix B(6,24)
//B(0,0)=L00*dh1dr+L01*dh1ds+L02*dh1dt;
//B(0,3)=L00*dh2dr+L01*dh2ds+L02*dh2dt;
//B(0,6)=L00*dh3dr+L01*dh3ds+L02*dh3dt;
//B(0,9)=L00*dh4dr+L01*dh4ds+L02*dh4dt;
//B(0,12)=L00*dh5dr+L01*dh5ds+L02*dh5dt;
//B(0,15)=L00*dh6dr+L01*dh6ds+L02*dh6dt;
//B(0,18)=L00*dh7dr+L01*dh7ds+L02*dh7dt;
//B(0,21)=L00*dh8dr+L01*dh8ds+L02*dh8dt;
//
//B(1,1)=L10*dh1dr+L11*dh1ds+L12*dh1dt;
//B(1,4)=L10*dh2dr+L11*dh2ds+L12*dh2dt;
//B(1,7)=L10*dh3dr+L11*dh3ds+L12*dh3dt;
//B(1,10)=L10*dh4dr+L11*dh4ds+L12*dh4dt;
//B(1,13)=L10*dh5dr+L11*dh5ds+L12*dh5dt;
//B(1,16)=L10*dh6dr+L11*dh6ds+L12*dh6dt;
//B(1,19)=L10*dh7dr+L11*dh7ds+L12*dh7dt;
//B(1,22)=L10*dh8dr+L11*dh8ds+L12*dh8dt;
//
//B(2,2)=L20d*h1dr+L21*dh1ds+L22*dh1dt;
//B(2,5)=L20d*h2dr+L21*dh2ds+L22*dh2dt;
//B(2,8)=L20d*h3dr+L21*dh3ds+L22*dh3dt;
//B(2,11)=L20*dh4dr+L21*dh4ds+L22*dh4dt;
//B(2,14)=L20*dh5dr+L21*dh5ds+L22*dh5dt;
//B(2,17)=L20*dh6dr+L21*dh6ds+L22*dh6dt;
//B(2,20)=L20*dh7dr+L21*dh7ds+L22*dh7dt;
//B(2,23)=L20*dh8dr+L21*dh8ds+L22*dh8dt;
//
//B(3,0)=B(1,1);
//B(3,1)=B(0,0);
//B(3,3)=B(1,4);
//B(3,4)=B(0,3);
//B(3,6)=B(1,7);
//B(3,7)=B(0,6);
//B(3,9)=B(1,10);
//B(3,10)=B(0,9);
//B(3,12)=B(1,13);
//B(3,13)=B(0,12);
//B(3,15)=B(1,16);
//B(3,16)=B(0,15);
//B(3,18)=B(1,19);
//B(3,19)=B(0,18);
//B(3,21)=B(1,22);
//B(3,22)=B(0,21);
//
//B(4,1)=B(2,2);
//B(4,2)=B(1,1);
//B(4,4)=B(2,5);
//B(4,5)=B(1,4);
//B(4,7)=B(2,8);
//B(4,8)=B(1,7);
//B(4,10)=B(2,11);
//B(4,11)=B(1,10);
//B(4,13)=B(2,14);
//B(4,14)=B(1,13);
//B(4,16)=B(2,17);
//B(4,17)=B(1,16);
//B(4,19)=B(2,20);
//B(4,20)=B(1,19);
//B(4,21)=B(2,23);
//B(4,23)=B(1,22);
//
//B(5,0)=B(2,2);
//B(5,2)=B(0,0);
//B(5,3)=B(2,5);
//B(5,5)=B(0,3);
//B(5,6)=B(2,8);
//B(5,8)=B(0,6);
//B(5,9)=B(2,11);
//B(5,11)=B(0,9);
//B(5,12)=B(2,14);
//B(5,14)=B(0,12);
//B(5,15)=B(2,17);
//B(5,17)=B(0,15);
//B(5,18)=B(2,20);
//B(5,20)=B(2,18);
//B(5,21)=B(0,23);
//B(5,23)=B(0,21);
//
//B(3,3)= L00*dh2dr+L01*dh2ds+L02*dh2dt;
//B(3,6)= L00*dh3dr+L01*dh3ds+L02*dh3dt;
//B(3,9)= L00*dh4dr+L01*dh4ds+L02*dh4dt;
//B(3,12)=L00*dh5dr+L01*dh5ds+L02*dh5dt;
//B(3,15)=L00*dh6dr+L01*dh6ds+L02*dh6dt;
//B(3,18)=L00*dh7dr+L01*dh7ds+L02*dh7dt;
//B(3,21)=L00*dh8dr+L01*dh8ds+L02*dh8dt;
//double
//EightNodeBrick::formDetJ (double r, double s, double t)
//{
//    return J(0,0)*J(1,1)*J(2,2)+J(1,0)*J(2,1)*J(0,2)+J(2,0)*J(0,1)*J(1,2)
//         - J(2,0)*J(1,1)*J(0,2)-J(0,0)*J(2,1)*J(1,2)-J(0,1)*J(1,0)*J(2,2);
//}
//
