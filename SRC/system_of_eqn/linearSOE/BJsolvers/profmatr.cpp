//############################################################################
//#                                                                          #
//#             /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/~~\              #
//#            |                                          |____|             #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |                 B A S E                  |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |              C L A S S E S               |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |          C + +     S O U R C E           |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#         /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/   |                  #
//#        |                                         |    |                  #
//#         \_________________________________________\__/                   #
//#                                                                          #
//#                                                                          #
//############################################################################
//
//   "C makes it easy to shoot yourself in the foot, C++ makes it harder,
//   but when you do, it blows away your whole leg" -- Bjarne Stroustrup
//
//##############################################################################
#ifndef PROFMATR_CC
#define PROFMATR_CC

#include "profmatr.h"


//##############################################################################
//# COPYRIGHT (C):     :-))                                                    #
//# PROJECT:           Object Oriented Finite Element Program                  #
//# PURPOSE:                                                                   #
//# CLASS:             profmatrix                                               #
//#                                                                            #
//# VERSION:                                                                   #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.10, SUN C++ ver=2.1 )#
//# TARGET OS:         DOS || UNIX || . . .                                    #
//# PROGRAMMER(S):     Boris Jeremic                                           #
//#                                                                            #
//#                                                                            #
//# DATE:              September 12-13 '94 looking for the solver for symmetric#
//#                                        and unsymmetric sparse matrices.    #
//#                                        One solution is Taylor's profile    #
//#                                        solver ( see FEM4ed by O.Z. and R.T.#
//#                    September 27 '94 profile solver for symmetric and       #
//#                                     Nonsymmetric systems works!            #
//#                                     (from FEM4ed by O.Z. and R.T.)         #
//#                                                                            #
//#                                                                            #
//##############################################################################
// this one is based on Zienkiewicz's and Taylor's book!

//##########################################################################


//tempout//##########################################################################
//tempout// create the structure and initialization of profilematrix!___Zhaohui 07-06-99
//tempoutprofilematrix::profilematrix(FEModelData & FEMD, Brick3D * b3d, Node *  node)
//tempout//skymatrix::skymatrix(FEModelData & FEMD, Finite_Element & FE)
//tempout  {
//tempout// create the structure:
//tempout     pc_profilematrix_rep = new profilematrix_rep; // this 'new' is overloaded
//tempout
//tempout
//tempout    int Number_of_DOFs = FEMD.get_number_of_DOFs();
//tempout//    pc_profilematrix_rep->square_dim = Number_of_DOFs;
//tempout// create ColumnHeight _________________________________________________
//tempout// column height starts from 1 to  Number_of_DOFs.
//tempout    pc_profilematrix_rep->columnheight = new int [Number_of_DOFs+1];
//tempout      if (!pc_profilematrix_rep->columnheight)
//tempout        {
//tempout          ::printf("\a\nInsufficient memory for pc_skymatrix_rep->columnheight\n");
//tempout          ::exit(1);
//tempout        }
//tempout    for (int count = 1 ; count <= Number_of_DOFs ; count++) 
//tempout
//tempout      //**** Changing 1 to 0 to set initial value of colheight. 
//tempout      pc_profilematrix_rep->columnheight[count]=0;  
//tempout
//tempout    int Number_of_Elements = FEMD.get_number_of_Bricks();
//tempout
//tempout    int II = 0;
//tempout    int *LM = NULL;
//tempout    int ME = 0;
//tempout
//tempout    int Number_of_DOFs_for_brick = 24;
//tempout
//tempout    for (int el_count = 1 ; el_count <= Number_of_Elements ; el_count++)
//tempout      {
//tempout        b3d[el_count].set_LM(node);
//tempout//        b3d[el_count].reportLM("LM");
//tempout      }
//tempout
//tempout// KJBathe: Page 1002 sub COLHT
//tempout    int LS = 1000000;
//tempout    int temp=0;
//tempout
//tempout    for (int el_cnt =1 ; el_cnt <= Number_of_Elements ; el_cnt++)
//tempout      {
//tempout        LS = 1000000; 
//tempout        LM = b3d[el_cnt].get_LM();
//tempout        for (int count01 = 0 ; count01 < Number_of_DOFs_for_brick ; count01++)
//tempout          {
//tempout            if (LM[count01] !=0 )
//tempout              {
//tempout//              ::printf("+++++++ LM[%d]= %d  LS = %d \n", count01, LM[count01],  LS);
//tempout                temp =  LS;
//tempout//              ::printf("temp = %d  LM[%d]= %d \n", temp, count01, LM[count01]);
//tempout                temp = LM[count01] - temp;
//tempout//              ::printf("temp = %d  LM[%d]= %d  LS = %d \n", temp, count01, LM[count01],  LS);
//tempout                if ( temp < 0 ) 
//tempout                  {
//tempout                     LS = LM[count01];
//tempout//                 ::printf("(LM[count01]-LS)<0->> LM[%d]= %d  LS = %d \n", count01, LM[count01],  LS);
//tempout                  }
//tempout
//tempout              }
//tempout          // getchar();
//tempout          
//tempout          }
//tempout      
//tempout      
//tempout        for (int count02 = 0 ; count02 < Number_of_DOFs_for_brick ; count02++)
//tempout          {
//tempout            II = LM[count02];
//tempout          // ::printf("II= %d  LM[%d] = %d \n", II,count02, LM[count02] );
//tempout
//tempout            if (II!=0) 
//tempout              {
//tempout                ME = II-LS;
//tempout         //  ::printf("ME= %d  II = %d  LS = %d \n", ME, II, LS );
//tempout
//tempout         //  ::printf("ME= %d  II = %d  ColumnHeight[%d] = %d \n", ME, II, II, ColumnHeight[II]);
//tempout              if (ME > pc_profilematrix_rep->columnheight[II]) 
//tempout                {
//tempout                  pc_profilematrix_rep->columnheight[II] = ME;
//tempout           //  ::printf("ME > ColumnHeight[II] ->> ME= %d  ColumnHeight[%d] = %d \n", ME, II, ColumnHeight[II]);
//tempout           //  getchar();
//tempout                  }
//tempout              }
//tempout         }
//tempout      
//tempout//        ::printf("IN colheigth = \n");
//tempout
//tempout//    for (int count03 = 1 ; count03 <= Number_of_DOFs ; count03++)
//tempout//      {
//tempout//        ::printf(" %d ",  pc_profilematrix_rep->columnheight[count03] );
//tempout//      }
//tempout//      getchar();
//tempout      
//tempout
//tempout      }
//tempout
//tempout
//tempout
//tempout
//tempout//        ::printf(" \n ------+------ \n ");
//tempout
//tempout//--    //KJBathe pp 1002
//tempout
//tempout//// Create the pc_profilematrix->jp vector -----Zhaohui
//tempout// jp[i], i starts from 0 to Number_of_DOFs - 1 !!
//tempout// column height [i], i starts from 1 to Number_of_DOFs !!
//tempout    pc_profilematrix_rep->jp = new int [Number_of_DOFs];
//tempout    if (!pc_profilematrix_rep->jp)
//tempout      {
//tempout        ::printf("\a\nInsufficient memory for JP vector\n");
//tempout        ::exit(1);
//tempout      }
//tempout    pc_profilematrix_rep->jp[ 0 ] = 0;
//tempout    for (int count04 = 1 ; count04 < Number_of_DOFs ; count04++)
//tempout      {
//tempout	pc_profilematrix_rep->jp[count04]  = 
//tempout	  pc_profilematrix_rep->jp[count04 - 1] + pc_profilematrix_rep->columnheight[count04+1] ;
//tempout      }
//tempout    
//tempout    int Numb_of_DOF = FEMD.get_number_of_DOFs(); 
//tempout
//tempout// set value for member of pc_profilematrix_rep---Zhaohui
//tempout    pc_profilematrix_rep->neq = Numb_of_DOF;
//tempout    pc_profilematrix_rep->total_numb = pc_profilematrix_rep->jp[Numb_of_DOF-1];
//tempout    pc_profilematrix_rep->flag = 'N';
//tempout
//tempout//    ::printf("\n\n Profile____JP = \n");
//tempout//    for (int count05 = 0 ; count05 < Number_of_DOFs ; count05++)
//tempout//      {
//tempout//       ::printf(" %d ",  pc_profilematrix_rep->jp[count05] );
//tempout//      }
//tempout
//tempout// Initialization and allocation memory of profilematrix____zhaohui
//tempout//  void profilematrix::profile_init( FEModelData      & FEMD)
//tempout//    {
//tempout    int length_of_AL = pc_profilematrix_rep->total_numb;
//tempout
//tempout// allocate memory for the actual profilematrix as profilematrix
//tempout    pc_profilematrix_rep->al = new double [(size_t) length_of_AL];
//tempout      if (!pc_profilematrix_rep->al)
//tempout        { printf("\a\nInsufficient memory for AL\n"); exit(1);}
//tempout    pc_profilematrix_rep->ad = new double [(size_t) Numb_of_DOF];
//tempout      if (!pc_profilematrix_rep->ad)
//tempout        { printf("\a\nInsufficient memory for AD\n"); exit(1);}
//tempout    pc_profilematrix_rep->au = new double [(size_t) length_of_AL];
//tempout      if (!pc_profilematrix_rep->au)
//tempout        { printf("\a\nInsufficient memory for AU\n"); exit(1);}
//tempout
//tempout    for (int I = 0 ; I < length_of_AL ; I++)   // Initializing AL, AU.
//tempout      {
//tempout	 pc_profilematrix_rep->al[I] = 0.0;
//tempout	 pc_profilematrix_rep->au[I] = 0.0;
//tempout       }
//tempout    for (int I = 0 ; I < Numb_of_DOF ; I++)   // Initializing AD.
//tempout      {
//tempout	 pc_profilematrix_rep->ad[I] = 0.0;
//tempout       }
//tempout
//tempout
//tempout    pc_profilematrix_rep->n = 1;  // so far, there's one reference
//tempout
//tempout
//tempout// Initialization and allocation memory of profilematrix____zhaohui
//tempout
//tempout//--
//tempout//--   int Total_K_length = pc_profilematrix_rep->maxa[Number_of_DOFs+1];
//tempout//--
//tempout//--    pc_profilematrix_rep->data = new double [Total_K_length+1];
//tempout//--
//tempout//--      if (!pc_profilematrix_rep->data)
//tempout//--        {
//tempout//--          ::printf("\a\nInsufficient memory for pc_profilematrix_rep->data\n");
//tempout//--          ::exit(1);
//tempout//--        }
//tempout//--    for (int count08 = 1 ; count08 <= Number_of_DOFs ; count08++) 
//tempout//--      pc_profilematrix_rep->data[count08]=0;
//tempout//--
//tempout
//tempout
//tempout  }

////#############################################################################



profilematrix::profilematrix(int matrix_order,
                             int *jp_init,
                             char flag_init,
                             double *au_init_val,
                             double *al_init_val,
                             double *ad_init_val)
  {
 // create the structure:
    pc_profilematrix_rep = new profilematrix_rep; // this 'new' is overloaded

    pc_profilematrix_rep->neq = matrix_order;

// get space for the JP vector
    pc_profilematrix_rep->jp = new int[matrix_order];
// put all jp_init's in the jp
    for ( int j=0 ; j<pc_profilematrix_rep->neq ; j++ )
      pc_profilematrix_rep->jp[j] = jp_init[j];

    pc_profilematrix_rep->total_numb = pc_profilematrix_rep->jp[matrix_order-1];

// set flag for Symmetry of Nonsymmetry
    pc_profilematrix_rep->flag = flag_init;

// allocate memory for the actual profilematrix as profilematrix
    pc_profilematrix_rep->au = new double [(size_t) pc_profilematrix_rep->total_numb];
      if (!pc_profilematrix_rep->au)
        {
          ::printf("\a\nInsufficient memory for pc_profilematrix_rep->au\n");
          ::exit(1);
        }

    if ( pc_profilematrix_rep->flag == 'N' )
      {
        pc_profilematrix_rep->al = new double [(size_t) pc_profilematrix_rep->total_numb];
          if (!pc_profilematrix_rep->al)
            {
              ::printf("\a\nInsufficient memory for pc_profilematrix_rep->al\n");
              ::exit(1);
            }
      }
    else if ( pc_profilematrix_rep->flag == 'S' )
      {
        pc_profilematrix_rep->al = pc_profilematrix_rep->au;
      }
    else
      {
        ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
      }

    pc_profilematrix_rep->ad = new double [(size_t) pc_profilematrix_rep->neq];
      if (!pc_profilematrix_rep->ad)
        {
          ::printf("\a\nInsufficient memory for pc_profilematrix_rep->ad\n");
          ::exit(1);
        }


    pc_profilematrix_rep->n = 1;  // so far, there's one reference

// this constructors is just for testing purposes.
// The real one is to be initialized from 0.0 init_vals and then
// class stiffness matrix those values are to be altered!!!!!!!!!!!!!!!!!!!!!!!!
    for ( int iau=0 ; iau<pc_profilematrix_rep->total_numb ; iau++ )
      pc_profilematrix_rep->au[iau] = au_init_val[iau];

    for ( int iad=0 ; iad<pc_profilematrix_rep->neq ; iad++ )
      pc_profilematrix_rep->ad[iad] = ad_init_val[iad];

    if ( pc_profilematrix_rep->flag == 'N' )
      {
        for ( int ial=0 ; ial<pc_profilematrix_rep->total_numb ; ial++ )
          pc_profilematrix_rep->al[ial] = al_init_val[ial];
      }


  }

//##########################################################################
// overloaded. Different calling methods
profilematrix::profilematrix(int matrix_order,
                             int *jp_init,
                             char flag_init,
                             double au_init_val,
                             double al_init_val,
                             double ad_init_val)
  {
 // create the structure:
    pc_profilematrix_rep = new profilematrix_rep; // this 'new' is overloaded

    pc_profilematrix_rep->neq = matrix_order;

// get space for the JP vector
    pc_profilematrix_rep->jp = new int[matrix_order];
// put all jp_init's in the jp
    for ( int j=0 ; j<pc_profilematrix_rep->neq ; j++ )
      pc_profilematrix_rep->jp[j] = jp_init[j];

    pc_profilematrix_rep->total_numb = pc_profilematrix_rep->jp[matrix_order-1];

// set flag for Symmetry of Nonsymmetry
    pc_profilematrix_rep->flag = flag_init;

// allocate memory for the actual profilematrix as profilematrix
    pc_profilematrix_rep->au = new double [(size_t) pc_profilematrix_rep->total_numb];
      if (!pc_profilematrix_rep->au)
        {
          ::printf("\a\nInsufficient memory for pc_profilematrix_rep->au\n");
          ::exit(1);
        }

    if ( pc_profilematrix_rep->flag == 'N' )
      {
        pc_profilematrix_rep->al = new double [(size_t) pc_profilematrix_rep->total_numb];
          if (!pc_profilematrix_rep->al)
            {
              ::printf("\a\nInsufficient memory for pc_profilematrix_rep->al\n");
              ::exit(1);
            }
      }
    else if ( pc_profilematrix_rep->flag == 'S' )
      {
        pc_profilematrix_rep->al = pc_profilematrix_rep->au;
      }
    else
      {
        ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
      }

    pc_profilematrix_rep->ad = new double [(size_t) pc_profilematrix_rep->neq];
      if (!pc_profilematrix_rep->ad)
        {
          ::printf("\a\nInsufficient memory for pc_profilematrix_rep->ad\n");
          ::exit(1);
        }


    pc_profilematrix_rep->n = 1;  // so far, there's one reference

// The real one is to be initialized from 0.0 init_vals and then
// class stiffness matrix those values are to be altered!!!!!!!!!!!!!!!!!!!!!!!!
    for ( int iau=0 ; iau<pc_profilematrix_rep->total_numb ; iau++ )
      pc_profilematrix_rep->au[iau] = au_init_val;

    for ( int iad=0 ; iad<pc_profilematrix_rep->neq ; iad++ )
      pc_profilematrix_rep->ad[iad] = ad_init_val;

    if ( pc_profilematrix_rep->flag == 'N' )
      {
        for ( int ial=0 ; ial<pc_profilematrix_rep->total_numb ; ial++ )
          pc_profilematrix_rep->al[ial] = al_init_val;
      }


  }

//##############################################################################
profilematrix::profilematrix(const profilematrix & x)   // copy initializer
 {
  x.pc_profilematrix_rep->n++; // we're adding another reference.
  pc_profilematrix_rep = x.pc_profilematrix_rep;  // point to the new profilematrix_rep.
 }



//##########################################################################
profilematrix::~profilematrix()
  {
    if (--pc_profilematrix_rep->n == 0)  // if reference count  goes to 0
  {
// DEallocate memory of the actual nDarray
//    delete [pc_nDarray_rep->pc_nDarray_rep->total_numb] pc_nDarray_rep->pd_nDdata;
//  see ELLIS & STROUSTRUP $18.3
//  and note on the p.65($5.3.4)
//  and the page 276 ($12.4)
    delete [] pc_profilematrix_rep->au;
    delete [] pc_profilematrix_rep->ad;
    if ( pc_profilematrix_rep->flag == 'N' )
      {
        delete [] pc_profilematrix_rep->al;
      }
    delete [] pc_profilematrix_rep->jp;
    //delete pc_profilematrix_rep;
  }
}

////##########################################################################
profilematrix & profilematrix::operator=(const profilematrix & rval)
  {
    rval.pc_profilematrix_rep->n++; // we're adding another reference.
//    rval.reference_count(+1);  // tell the rval it has another reference
// It is important to increment the reference_counter in the new
// tensor before decrementing the reference_counter in the
// old tensor_rep to ensure proper operation when assigning a
// tensor_rep to itself ( after ARKoenig JOOP May/June '90 )
// clean up current value;
    if( reference_count(-1) == 0)  // if nobody else is referencing us.
      {
        delete [] pc_profilematrix_rep->au;
        delete [] pc_profilematrix_rep->ad;
        if ( pc_profilematrix_rep->flag == 'N' )
          {
            delete [] pc_profilematrix_rep->al;
          }
        delete [] pc_profilematrix_rep->jp;
        delete pc_profilematrix_rep;
      }
// connect to new value
    pc_profilematrix_rep = rval.pc_profilematrix_rep;// point at the rval profilematrix_rep

    return *this;
  }

//##########################################################################
int profilematrix::dimension_of_profile_M(void) const // dimension of profile matrix
  {
    return pc_profilematrix_rep->neq;
  }

//##########################################################################
int * profilematrix::get_jp(void) const  // get pointer to array of
  {                                      // Locations of Diagonals
    return pc_profilematrix_rep->jp;
  }

//##########################################################################
//Zhaohui  added to set the max and min element
double profilematrix::get_max(void) // 
 {                                                 
    double max_element = -1000000.0;
    for (int i =0; i< pc_profilematrix_rep->neq; i++)
     //{
        if (pc_profilematrix_rep->ad[i] > max_element)
           max_element = pc_profilematrix_rep->ad[i];
//	 else if (pc_profilematrix_rep->ad[i] < min_element)
//           min_element = pc_profilematrix_rep->ad[i];

      //}
    //for ( int i=0 ; i<pc_profilematrix_rep->total_numb ; i++ )
    // {
    //    if (pc_profilematrix_rep->au[i] > pc_profilematrix_rep->max_element)
    //       pc_profilematrix_rep->max_element = pc_profilematrix_rep->au[i];
    //	 else if (pc_profilematrix_rep->au[i] < pc_profilematrix_rep->min_element)
    //       pc_profilematrix_rep->min_element = pc_profilematrix_rep->au[i];
    //
    //  }
    //for ( int i=0 ; i<pc_profilematrix_rep->total_numb ; i++ )
    // {
    //    if (pc_profilematrix_rep->al[i] > pc_profilematrix_rep->max_element)
    //       pc_profilematrix_rep->max_element = pc_profilematrix_rep->al[i];
    //	 else if (pc_profilematrix_rep->al[i] < pc_profilematrix_rep->min_element)
    //       pc_profilematrix_rep->min_element = pc_profilematrix_rep->al[i];
    //
    //  }

    return max_element;
 }


//##########################################################################
//Zhaohui  added to set the max-element in diagonal of eqn_no_shake!
void profilematrix::set_penalty_element(int *eqn_no_shake, int no_of_shake, double max_element) 
 {                                                  
 //double scale = 100000.0;
 for (int count00 = 0 ; count00 < no_of_shake; count00++)
  {
//     ::printf("mass_pen before :  %+10.3e, %d ", pc_profilematrix_rep->ad[eqn_no_shake[count00]-1] , eqn_no_shake[count00]);
//     pc_profilematrix_rep->ad[eqn_no_shake[count00]] = max_element;
     // correct for eq_no =1, 2 ,..., but ad starts from 0.
     pc_profilematrix_rep->ad[eqn_no_shake[count00]-1] = max_element;
//     ::printf("mass_pen after  :  %+10.3e, %d ", pc_profilematrix_rep->ad[eqn_no_shake[count00]-1] , eqn_no_shake[count00]);
  }
}

//##########################################################################
//double  profilematrix::get_max(void)        // get pointer to array of
//  {                                      // Locations of Diagonals
//    return pc_profilematrix_rep->max_element;
//  }
//##########################################################################
//double  profilematrix::get_min(void)        // get pointer to array of
//  {                                      // Locations of Diagonals
//    return pc_profilematrix_rep->min_element;
//  }
//
////##############################################################################
// profilematrix addition
profilematrix& profilematrix::operator+=(const profilematrix & rval)
  {
    long int this_total_numb = this->pc_profilematrix_rep->total_numb;
    long int rval_total_numb =  rval.pc_profilematrix_rep->total_numb;
//    this_total_numb = 696;
//    rval_total_numb = 696;

    if(this_total_numb != rval_total_numb)
      {
        ::printf("\a\nprofilematrixs of different sizes: += not YET possible\n");
        ::exit ( 1 );
      }

// Copy *this if necessary
    if ( this->pc_profilematrix_rep->n > 1 )// see ARK in JOOP may/june '90
      {                                    // "Letter From a Newcomer"
// create the structure:
        profilematrix_rep * New_pc_profilematrix_rep = new profilematrix_rep; // this 'new' is overloaded

        New_pc_profilematrix_rep->neq = this->pc_profilematrix_rep->neq ;

// get space for the JP vector
        New_pc_profilematrix_rep->jp = new int[New_pc_profilematrix_rep->neq];
// put all jp_init's in the jp
        for ( int j=0 ; j<this->pc_profilematrix_rep->neq ; j++ )
          New_pc_profilematrix_rep->jp[j] = this->pc_profilematrix_rep->jp[j];

        New_pc_profilematrix_rep->total_numb = this->pc_profilematrix_rep->total_numb;

// set flag for Symmetry of Nonsymmetry
        New_pc_profilematrix_rep->flag = this->pc_profilematrix_rep->flag;

// allocate memory for the actual profilematrix as profilematrix
        New_pc_profilematrix_rep->au = new double [(size_t) New_pc_profilematrix_rep->total_numb];
          if (!New_pc_profilematrix_rep->au)
            {
              ::printf("\a\nInsufficient memory for New_pc_profilematrix_rep->au\n");
              ::exit(1);
            }

        if ( New_pc_profilematrix_rep->flag == 'N' )
          {
            New_pc_profilematrix_rep->al = new double [(size_t) New_pc_profilematrix_rep->total_numb];
              if (!New_pc_profilematrix_rep->al)
                {
                  ::printf("\a\nInsufficient memory for New_pc_profilematrix_rep->al\n");
                  ::exit(1);
                }
          }
        else if ( New_pc_profilematrix_rep->flag == 'S' )
          {
            New_pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
          }
        else
          {
            ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
          }

        New_pc_profilematrix_rep->ad = new double [(size_t) New_pc_profilematrix_rep->neq];
          if (!New_pc_profilematrix_rep->ad)
            {
              ::printf("\a\nInsufficient memory for New_pc_profilematrix_rep->ad\n");
              ::exit(1);
            }

        New_pc_profilematrix_rep->n = 1;  // so far, there's one reference

// The real one is to be initialized from 0.0 init_vals and then
// class stiffness matrix those values are to be altered!!!!!!!!!!!!!!!!!!!!!!!!
//        for ( int iau=0 ; iau<New_pc_profilematrix_rep->total_numb ; iau++ )
        for ( int iau=0 ; iau<this_total_numb ; iau++ )
          New_pc_profilematrix_rep->au[iau] = this->pc_profilematrix_rep->au[iau];

        for ( int iad=0 ; iad<New_pc_profilematrix_rep->neq ; iad++ )
          New_pc_profilematrix_rep->ad[iad] = this->pc_profilematrix_rep->ad[iad];

        if ( New_pc_profilematrix_rep->flag == 'N' )
          {
            for ( int ial=0 ; ial<this_total_numb ; ial++ )
              New_pc_profilematrix_rep->al[ial] = this->pc_profilematrix_rep->al[ial];
          }
        else if ( New_pc_profilematrix_rep->flag == 'S' )
          {
            New_pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
          }
        else
          {
            ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
          }


//??????????????????
//        this->pc_profilematrix_rep->total_numb--;
        this->pc_profilematrix_rep = New_pc_profilematrix_rep;
//..............................................................................
      }
// I can add this two profilematrices just as a simple vectors:
//    for (int kau=0 ; kau<this->pc_profilematrix_rep->total_numb ; kau++)
    for (int kau=0 ; kau<this_total_numb ; kau++)
      this->pc_profilematrix_rep->au[kau] += rval.pc_profilematrix_rep->au[kau];
    for (int kad=0 ; kad<this->pc_profilematrix_rep->neq ; kad++)
      this->pc_profilematrix_rep->ad[kad] += rval.pc_profilematrix_rep->ad[kad];
    if ( this->pc_profilematrix_rep->flag == 'N' )
      {
        for ( int kal=0 ; kal<this_total_numb ; kal++ )
         this->pc_profilematrix_rep->al[kal] += rval.pc_profilematrix_rep->al[kal];
      }
//..// this is OK
//..    else if ( this->pc_profilematrix_rep->flag == 'S' )
//..      {
//..        this->pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
//..      }
//..    else
//..      {
//..        ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
//..      }
//..
    return *this;
  }


//##############################################################################
// profilematrix addition
profilematrix operator+(const profilematrix & lval, const profilematrix & rval)
  {
    profilematrix result(lval);
    result += rval;
    return result;
  }


//##############################################################################
// scalar addition
profilematrix profilematrix::operator+( double rval)
  {
// construct profilematrix using the same control numbers as for the
// original one.
    profilematrix add(this->pc_profilematrix_rep->neq,
                      this->pc_profilematrix_rep->jp,
                      this->pc_profilematrix_rep->flag,
                      this->pc_profilematrix_rep->al,
                      this->pc_profilematrix_rep->au,
                      this->pc_profilematrix_rep->ad);

// I can add this two profilematrices just as a simple vectors:
    for (int kau=0 ; kau<this->pc_profilematrix_rep->total_numb ; kau++)
      add.pc_profilematrix_rep->au[kau] += rval;
    for (int kad=0 ; kad<this->pc_profilematrix_rep->neq ; kad++)
      add.pc_profilematrix_rep->ad[kad] += rval;
    if ( this->pc_profilematrix_rep->flag == 'N' )
      {
        for ( int kal=0 ; kal<this->pc_profilematrix_rep->total_numb ; kal++ )
         add.pc_profilematrix_rep->al[kal] += rval;
      }
//..// this is OK!
//..    else if ( this->pc_profilematrix_rep->flag == 'S' )
//..      {
//..        this->pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
//..      }
//..    else
//..      {
//..        ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
//..      }

    return add;
 }

//##############################################################################
// profilematrix substraction
profilematrix& profilematrix::operator-=(const profilematrix & rval)
  {
    long int this_total_numb = this->pc_profilematrix_rep->total_numb;
    long int rval_total_numb =  rval.pc_profilematrix_rep->total_numb;

    if(this_total_numb != rval_total_numb)
      {
        ::printf("\a\nprofilematrixs of different sizes: -= not YET possible\n");
        ::exit ( 1 );
      }

// Copy *this if necessary
    if ( this->pc_profilematrix_rep->n > 1 )// see ARK in JOOP may/june '90
      {                                    // "Letter From a Newcomer"
// create the structure:
        profilematrix_rep * New_pc_profilematrix_rep = new profilematrix_rep; // this 'new' is overloaded

        New_pc_profilematrix_rep->neq = this->pc_profilematrix_rep->neq ;

// get space for the JP vector
        New_pc_profilematrix_rep->jp = new int[New_pc_profilematrix_rep->neq];
// put all jp_init's in the jp
        for ( int j=0 ; j<this->pc_profilematrix_rep->neq ; j++ )
          New_pc_profilematrix_rep->jp[j] = this->pc_profilematrix_rep->jp[j];

        New_pc_profilematrix_rep->total_numb = New_pc_profilematrix_rep->total_numb;

// set flag for Symmetry of Nonsymmetry
        New_pc_profilematrix_rep->flag = this->pc_profilematrix_rep->flag;

// allocate memory for the actual profilematrix as profilematrix
        New_pc_profilematrix_rep->au = new double [(size_t) New_pc_profilematrix_rep->total_numb];
          if (!New_pc_profilematrix_rep->au)
            {
              ::printf("\a\nInsufficient memory for New_pc_profilematrix_rep->au\n");
              ::exit(1);
            }

        if ( New_pc_profilematrix_rep->flag == 'N' )
          {
            New_pc_profilematrix_rep->al = new double [(size_t) New_pc_profilematrix_rep->total_numb];
              if (!New_pc_profilematrix_rep->al)
                {
                  ::printf("\a\nInsufficient memory for New_pc_profilematrix_rep->al\n");
                  ::exit(1);
                }
          }
        else if ( New_pc_profilematrix_rep->flag == 'S' )
          {
            New_pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
          }
        else
          {
            ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
          }

        New_pc_profilematrix_rep->ad = new double [(size_t) New_pc_profilematrix_rep->neq];
          if (!New_pc_profilematrix_rep->ad)
            {
              ::printf("\a\nInsufficient memory for New_pc_profilematrix_rep->ad\n");
              ::exit(1);
            }

        New_pc_profilematrix_rep->n = 1;  // so far, there's one reference

// The real one is to be initialized from 0.0 init_vals and then
// class stiffness matrix those values are to be altered!!!!!!!!!!!!!!!!!!!!!!!!
        for ( int iau=0 ; iau<New_pc_profilematrix_rep->total_numb ; iau++ )
          New_pc_profilematrix_rep->au[iau] = this->pc_profilematrix_rep->au[iau];

        for ( int iad=0 ; iad<New_pc_profilematrix_rep->neq ; iad++ )
          New_pc_profilematrix_rep->au[iad] = this->pc_profilematrix_rep->ad[iad];

        if ( New_pc_profilematrix_rep->flag == 'N' )
          {
            for ( int ial=0 ; ial<New_pc_profilematrix_rep->total_numb ; ial++ )
              New_pc_profilematrix_rep->al[ial] = this->pc_profilematrix_rep->al[ial];
          }
        else if ( New_pc_profilematrix_rep->flag == 'S' )
          {
            New_pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
          }
        else
          {
            ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
          }

        this->pc_profilematrix_rep->total_numb--;
        this->pc_profilematrix_rep = New_pc_profilematrix_rep;
//..............................................................................
      }
// I can add this two profilematrices just as a simple vectors:
    for (int kau=0 ; kau<this->pc_profilematrix_rep->total_numb ; kau++)
      this->pc_profilematrix_rep->au[kau] -= rval.pc_profilematrix_rep->au[kau];
    for (int kad=0 ; kad<this->pc_profilematrix_rep->neq ; kad++)
      this->pc_profilematrix_rep->ad[kad] -= rval.pc_profilematrix_rep->ad[kad];
    if ( this->pc_profilematrix_rep->flag == 'N' )
      {
        for ( int kal=0 ; kal<this->pc_profilematrix_rep->total_numb ; kal++ )
         this->pc_profilematrix_rep->al[kal] -= rval.pc_profilematrix_rep->al[kal];
      }
//..// this is OK!
//..    else if ( this->pc_profilematrix_rep->flag == 'S' )
//..      {
//..        this->pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
//..      }
//..    else
//..      {
//..        ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
//..      }

    return *this;
  }


//##############################################################################
// profilematrix substraction
profilematrix operator-(const profilematrix & lval, const profilematrix & rval)
  {
    profilematrix result(lval);
    result -= rval;
    return result;
  }


//##############################################################################
// scalar substraction
profilematrix profilematrix::operator-( double rval)
  {
// construct profilematrix using the same control numbers as for the
// original one.
    profilematrix sub(this->pc_profilematrix_rep->neq,
                      this->pc_profilematrix_rep->jp,
                      this->pc_profilematrix_rep->flag,
                      this->pc_profilematrix_rep->al,
                      this->pc_profilematrix_rep->au,
                      this->pc_profilematrix_rep->ad);

// I can add this two profilematrices just as a simple vectors:
    for (int kau=0 ; kau<this->pc_profilematrix_rep->total_numb ; kau++)
      sub.pc_profilematrix_rep->au[kau] -= rval;
    for (int kad=0 ; kad<this->pc_profilematrix_rep->neq ; kad++)
      sub.pc_profilematrix_rep->ad[kad] -= rval;
    if ( this->pc_profilematrix_rep->flag == 'N' )
      {
        for ( int kal=0 ; kal<this->pc_profilematrix_rep->total_numb ; kal++ )
         sub.pc_profilematrix_rep->al[kal] -= rval;
      }
//..// This is OK!
//..    else if ( this->pc_profilematrix_rep->flag == 'S' )
//..      {
//..        this->pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
//..      }
//..    else
//..      {
//..        ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
//..      }

    return sub;
 }

//##############################################################################
// profilematrix multiply by a vector___Zhaohui__07-07-1999                     
// vector arg(i), i from 1 to neq, but return[i], i from 0 to neq-1!!!!
// bug found in July 27,99__Zhaohui
 double * profilematrix::operator*( BJvector & arg)
  {
    int irow, iau, colht, jp, neq;
    neq = this->pc_profilematrix_rep->neq;

    if( neq != arg.rows())
      error("# rows of second mat must equal "
               "# cols of first for multiply#");
    double  *result;
    result = new double [ neq];
    if (! result)
        {
          ::printf("\a\nInsufficient memory for result\n");
          ::exit(1);
        }
    // initializing result!
    for (int init=0; init<neq; init++)
       result[init] = 0.0;

//    this->profile_al_print();
//    for( int i=0 ; i<neq ; i++ )
//       ::printf(" \n jp=%d",this->pc_profilematrix_rep->jp[i] );

    for( int i=0 ; i<neq ; i++ )
     {
      result[i] += this->pc_profilematrix_rep->ad[i] * arg.val(i+1);
//       ::printf(" \n\nI%d %+6.2e%+6.2e%+6.2e\n\n",i,this->pc_profilematrix_rep->ad[i],arg.val(i+1),result[i+1] );
      jp = this->pc_profilematrix_rep->jp[i];
      if ( i==0 )
         {
	  colht = 0;
	 }
        else
         {
	  colht = this->pc_profilematrix_rep->jp[i]-this->pc_profilematrix_rep->jp[i-1];
	 }

      for( int j=1 ; j<=colht ; j++ )
       {
	  // remember that au[i], i starts from 0!!!
 	  //iau = jp - j +1;
	  iau = jp - j ;
	  irow = i  -j ;
	  result[irow] +=this->pc_profilematrix_rep->au[iau] * arg.val(i+1);
	  result[i]  +=this->pc_profilematrix_rep->al[iau] * arg.val(irow+1);
//       ::printf(" \n\nJ%d %+6.2e%+6.2e%+6.2e\n\n",j,this->pc_profilematrix_rep->au[iau],arg.val(i+1),result[irow] );
       }
     
//      ::printf(" \n\n  %+6.2e\n", result[i+1] );
      
     }
//   for( int i=0 ; i<neq ; i++ )
//       ::printf(" \n result=%+6.2e arg=%+6.2e\n", result[i+1], arg.val(i+1) );
  return result;
  }


//..//##############################################################################
// scalar multiplication
profilematrix  profilematrix::operator*( double rval)
 {
// construct profilematrix using the same control numbers as for the
// original one.
    profilematrix mult(this->pc_profilematrix_rep->neq,
                       this->pc_profilematrix_rep->jp,
                       this->pc_profilematrix_rep->flag,
                       this->pc_profilematrix_rep->al,
                       this->pc_profilematrix_rep->au,
                       this->pc_profilematrix_rep->ad);

  for (int kau=0 ; kau<pc_profilematrix_rep->total_numb ; kau++)
     mult.pc_profilematrix_rep->au[kau] *= rval;
   for (int kad=0 ; kad<this->pc_profilematrix_rep->neq ; kad++)
     mult.pc_profilematrix_rep->ad[kad] *= rval;
  if ( this->pc_profilematrix_rep->flag == 'N' )
    {
      for ( int kal=0 ; kal<this->pc_profilematrix_rep->total_numb ; kal++ )
       mult.pc_profilematrix_rep->al[kal] *= rval;
    }
//    for (int kau=0 ; kau<this->pc_profilematrix_rep->total_numb ; kau++)
//      this->pc_profilematrix_rep->au[kau] *= rval;
//    for (int kad=0 ; kad<this->pc_profilematrix_rep->neq ; kad++)
//      this->pc_profilematrix_rep->ad[kad] *= rval;
//    if ( this->pc_profilematrix_rep->flag == 'N' )
//      {
//        for ( int kal=0 ; kal<this->pc_profilematrix_rep->total_numb ; kal++ )
//         this->pc_profilematrix_rep->al[kal] *= rval;
//      }
//..// This is OK
//..    else if ( this->pc_profilematrix_rep->flag == 'S' )
//..      {
//..        this->pc_profilematrix_rep->al = New_pc_profilematrix_rep->au;
//..      }
//..    else
//..      {
//..        ::printf("flags not set properly, flag=='S' symmetric; flag=='N' NON(un)symmetric\n");
//..      }

    return  mult;
 }



// only for the members inside
// the profile !!!!!!!!!!!
// so row and col must be within profile !!!!!!!!!!!!!!!
//##########################################################################
// this one is the same as mval except that it is more convenient  //  ________
// to overload operator (row,col).                                 // |  rows
double & profilematrix::operator( )(int row, int col) const        // |c
  {                                                                // |o
    if( row==col )                                                 // |l
      {                                                            // |u
        return this->pc_profilematrix_rep->ad[row-1];              // |m
      }                                                            // |n
    else if( row<col ) // upper part AU                                 // |s
      {
        int  next_to_diag_memb = *(this->pc_profilematrix_rep->jp+col-1);//-1:starts from 0
        int member_of_profile_n = next_to_diag_memb - (col-1) + row;
        double *member_of_profile = &(pc_profilematrix_rep->au[member_of_profile_n-1]);//-1:starts from 0
        return (*member_of_profile);
      }
    else // if( col<row ) // lower part AL
      {
        int  next_to_diag_memb = *(this->pc_profilematrix_rep->jp+row-1);//-1:starts from 0
        int member_of_profile_n = next_to_diag_memb - (row-1) + col;
        double *member_of_profile = &(pc_profilematrix_rep->al[member_of_profile_n-1]);//-1:starts from 0
        return (*member_of_profile);
      }
//..    else
//..      {
//..        double zero=0.0;
//..        return zero;
//..      }
  }

//##########################################################################
double & profilematrix::mval(int row, int col) const // only for the members inside
  {                                                  // the profile !!!!!!!!!!!
    if( row==col )                                   // so row and col must be
      {                                              // within profile !!!!!!!!!!!!!!!
        return this->pc_profilematrix_rep->ad[row-1];
      }
    else if( row<col ) // upper part AU
      {
        int  next_to_diag_memb = *(this->pc_profilematrix_rep->jp+col-1);//-1:starts from 0
        int member_of_profile_n = next_to_diag_memb - (col-1) + row;
        double *member_of_profile = &(pc_profilematrix_rep->au[member_of_profile_n-1]);//-1:starts from 0
        return (*member_of_profile);
      }
    else //if( col<row ) // lower part AL
      {
        int  next_to_diag_memb = *(this->pc_profilematrix_rep->jp+row-1);//-1:starts from 0
        int member_of_profile_n = next_to_diag_memb - (row-1) + col;
        double *member_of_profile = &(pc_profilematrix_rep->al[member_of_profile_n-1]);//-1:starts from 0
        return (*member_of_profile);
      }
//..    else
//..      {
//..        double zero=0.0;
//..        return zero;
//..      }
  }



//##########################################################################
// full_val function allows you to treat profilematrix as if it is full
// float matrix. The function will calculate position inside profile matrix
// and return appropriate number if row and col are below/right profile line or
// return zero (0) if row and col are above/left profile line
double profilematrix::full_val(int row, int col) const
  {
    if( row==col )
      {
        return this->pc_profilematrix_rep->ad[row-1];
      }
    else if( row<col ) // upper part AU
      {
        int total_column_depth = col-1;
        int actual_column_position = row;  // - total_column_depth + 2;
        int real_column_depth = pc_profilematrix_rep->jp[col-1] -
                                 pc_profilematrix_rep->jp[col-1-1];
        int zeros_above_profile =
          total_column_depth-real_column_depth;
        int how_much_above_profile =
          zeros_above_profile-actual_column_position+1;
        if ( how_much_above_profile > 0 )
          {
            double back = 0.0;    ////    // not very smart ???????
            return ( back );
          }
        else
          {
            int  next_to_diag_memb =
              *(this->pc_profilematrix_rep->jp+col-1);//-1:starts from 0
            int member_of_profile_n =
              next_to_diag_memb - total_column_depth + row - 1;
            double *member_of_profile =
              &(pc_profilematrix_rep->au[member_of_profile_n]);
            return (*member_of_profile);
          }
      }
    else if( col<row ) // lower part AL
      {
        int total_row_width = row-1;
        int actual_row_position = col;  // - total_row_width + 2;
        int real_row_width = pc_profilematrix_rep->jp[row-1] -
                             pc_profilematrix_rep->jp[row-1-1];
        int zeros_left_of_profile =
          total_row_width-real_row_width;
        int how_much_left_of_profile =
          zeros_left_of_profile-actual_row_position+1;
        if ( how_much_left_of_profile > 0 )
          {
            double back = 0.0;    ////    // not very smart ???????
            return ( back );
          }
        else
          {
            int  next_to_diag_memb =
              *(this->pc_profilematrix_rep->jp+row-1);//-1:starts from 0
            int member_of_profile_n =
              next_to_diag_memb - total_row_width + col - 1;
            double *member_of_profile =
              &(pc_profilematrix_rep->al[member_of_profile_n]);
            return (*member_of_profile);
          }
      }
    else return 0.0;
  }
     // used by profmatrix functions which KNOW they aren't
     // exceeding the boundaries


//##########################################################################
double & profilematrix::val(int row, int col)
  {
//    double zero=0.0;
    if ( row<=0 && row>=dimension_of_profile_M() && col>=dimension_of_profile_M() )
      {
        error("index out of range");
//        return zero;
      }
//    else
      return (mval(row,col));
  }

//##########################################################################
double profilematrix::cval(int row, int col) const
  {
    double zero=0.0;
    if ( row<=0 && row>=dimension_of_profile_M() && col>=dimension_of_profile_M() )
      {
        error("index out of range");
        return zero;
      }
//    else
      return (mval(row,col));
  }


//##########################################################################
double profilematrix::mmin()
  {
    double temp;
    if ( dimension_of_profile_M()<=0 )
      error("bad profilematrix size for min ()");
    double minimum = full_val(1,1);
    for ( int row=1 ; row<=dimension_of_profile_M() ; row++ )
      for ( int col=1 ; col<=dimension_of_profile_M() ; col++ )
        if ( (temp=full_val(row,col)) < minimum )
          minimum = temp;
    return minimum;
  }

//##########################################################################
double profilematrix::mmax()
  {
    double temp=0.0;
    if( dimension_of_profile_M()<=0 )
      error("bad profilematrix size for max()");
    double maximum = full_val(1,1);
    for ( int row=1 ; row<=dimension_of_profile_M() ; row++ )
      for ( int col=1 ; col<=dimension_of_profile_M() ; col++ )
        {
          if ( (temp=full_val(row,col)) > maximum )
          maximum = temp;
        }
    return maximum;
  }

//##########################################################################
double profilematrix::mean()
  {
    int col = 0;
    int row = 1; 
    double sum = 0;
    for ( row=1 ; row<=dimension_of_profile_M() ; row++ )
      for ( col=1 ; col<=dimension_of_profile_M() ; col++ )
        sum += fabs(full_val(row,col));
    return sum/(row*col);
  }


//..//##########################################################################
//..void profilematrix::lower_print(char *msg)
//..  {
//..    if (*msg) printf("%s\n",msg);
//..    for ( int row=1 ; row<=dimension_of_prof_M() ; row++ )
//..      {
//..        int total_column_height = row;
//..        int real_column_height = pc_profilematrix_rep->p_maxa[row] - pc_profilematrix_rep->p_maxa[row-1];
//..        int numb_of_voids = total_column_height - real_column_height;
//..        int n_of_voids_to_reach_number = numb_of_voids;
//..        for ( int col=1 ; col<=row ; col++ )
//..          {
//..            if( n_of_voids_to_reach_number > 0 )
//..              {
//..                 for(int void_count=1; void_count<=numb_of_voids; void_count++)
//..                  {
//..                    printf("   **   ");
//..                    col++;
//..                    n_of_voids_to_reach_number--;
//..                  }
//..              }
//..            printf( "%7.4f ", mval(row,col) );
//..          }
//..        printf("\n");
//..      }
//..  }
//..
//..//##########################################################################
//..void profilematrix::upper_print(char *msg)
//..  {
//..    if (*msg) printf("%s\n",msg);
//..    for ( int row=1 ; row<=dimension_of_prof_M() ; row++ )
//..      {
//..        for ( int voids=0 ; voids<row-1 ; voids++ ) printf("        ");
//..        for ( int col=row ; col<=dimension_of_prof_M() ; col++ )
//..          {
//..            int total_column_height = col;
//..            int actual_column_position = total_column_height - row + 1;
//..            int real_column_height = pc_profilematrix_rep->p_maxa[col] - pc_profilematrix_rep->p_maxa[col-1];
//..            int how_much_above_profline=actual_column_position-real_column_height;
//..            if ( how_much_above_profline > 0 )
//..              {
//..                printf("   **   ");
//..              }
//..            else
//..              {
//..                printf( "%7.4f ", mval(col,row) );
//..              }
//..          }
//..        printf("\n");
//..      }
//..  }


//tempout////#############################################################################
//tempout// Assemble global stiffness matrix from local contributions 
//tempout// and store it in porfile format!_____Zhaohui 06-30-99
//tempoutprofilematrix & profilematrix::AssembleBricksInProfMatrix(stiffness_matrix & Ke,
//tempout                                                          Brick3D          * b3d,
//tempout                                                          Node             * node,
//tempout		       				          FEModelData      & FEMD,
//tempout                                                          float              scale // used to add damping matrix C to Kbar Zhaohui 
//tempout							                         )
//tempout  {
//tempout
//tempout//    ::printf("\n\n");
//tempout
//tempout    int *lm = b3d->get_LM();
//tempout
//tempout
//tempout// KJBathe pp1003 ADDBAN
//tempout
//tempout    int DOF_count_for_Brick = 24;
//tempout    int II = 0;
//tempout    int MI = 0;
//tempout    int MJ = 0;
//tempout    int JJ = 0;
//tempout    int IJ = 0;
//tempout    int KK = 0;
//tempout
//tempout//    int Numb_of_DOF = FEMD.get_number_of_DOFs(); 
//tempout//    int length_of_AL = pc_profilematrix_rep->jp[Numb_of_DOF-1];
//tempout//
//tempout//// allocate memory for the actual profilematrix as profilematrix
//tempout//    pc_profilematrix_rep->al = new double [(size_t) length_of_AL];
//tempout//      if (!pc_profilematrix_rep->al)
//tempout//        { printf("\a\nInsufficient memory for AL\n"); exit(1);}
//tempout//    pc_profilematrix_rep->ad = new double [(size_t) Numb_of_DOF];
//tempout//      if (!pc_profilematrix_rep->ad)
//tempout//        { printf("\a\nInsufficient memory for AD\n"); exit(1);}
//tempout//    pc_profilematrix_rep->au = new double [(size_t) length_of_AL];
//tempout//      if (!pc_profilematrix_rep->au)
//tempout//        { printf("\a\nInsufficient memory for AU\n"); exit(1);}
//tempout//
//tempout//    for (int I = 0 ; I < length_of_AL ; I++)   // Initialization.
//tempout//      {
//tempout//	 pc_profilematrix_rep->al[I] = 0.0;
//tempout//	 pc_profilematrix_rep->au[I] = 0.0;
//tempout//       }
//tempout//    for (int I = 0 ; I < Numb_of_DOF ; I++)   // Initialization.
//tempout//      {
//tempout//	 pc_profilematrix_rep->ad[I] = 0.0;
//tempout//       }
//tempout//  THIS PART HAS BEEN MOVED OUTSIDE AND BECOME INDEPENDENT SUBROUTINE!!
//tempout
//tempout    for (int I = 1 ; I <= DOF_count_for_Brick ; I++)
//tempout      {
//tempout        II =  lm[I-1];
//tempout        if (II > 0) 
//tempout          {
//tempout            MI = pc_profilematrix_rep->jp[II - 1];  
//tempout            for (int J = 1 ; J <= DOF_count_for_Brick ; J++)
//tempout              {
//tempout                JJ = lm[J-1];
//tempout                if (JJ > 0)
//tempout                  {
//tempout                    MJ = pc_profilematrix_rep->jp[JJ - 1];//new variable corresponding to MI
//tempout                    IJ = II - JJ;
//tempout                    
//tempout		    if (IJ > 0 ) // Lower triangular!
//tempout                      { 
//tempout                        KK = MI - IJ + 1;
//tempout                        pc_profilematrix_rep->al[KK-1] = 
//tempout			  pc_profilematrix_rep->al[KK-1] + Ke.val(I,J)*scale;
//tempout                      }
//tempout                    
//tempout		    else if (IJ == 0 ) // Diagonal elements
//tempout                      { 
//tempout                        KK = II ;
//tempout                        pc_profilematrix_rep->ad[KK-1] = 
//tempout			  pc_profilematrix_rep->ad[KK-1] + Ke.val(I,J)*scale;
//tempout                      }
//tempout                    
//tempout		    else if (IJ < 0 )  // Upper triangular!
//tempout                      { 
//tempout                        KK = MJ + IJ + 1;
//tempout                        pc_profilematrix_rep->au[KK-1] = 
//tempout		          pc_profilematrix_rep->au[KK-1] + Ke.val(I,J)*scale;
//tempout                      }
//tempout
//tempout                  }
//tempout
//tempout              }
//tempout
//tempout          }
//tempout      }
//tempout  
//tempout
//tempout//// Print AL, AD, AU for chencking purpose.
//tempout//    ::printf("\n\nAD= ");
//tempout//      for (int count05 = 0 ; count05 < Numb_of_DOF ; count05++)
//tempout//      {
//tempout//       ::printf("  ad %6.2e ",  pc_profilematrix_rep->ad[count05] );
//tempout////       ::printf("  assem jp %d \n",  pc_profilematrix_rep->jp[count05] );
//tempout// }
//tempout//
//tempout//       ::printf(" \n\n AL= ");
//tempout//      for (int count05 = 0 ; count05 < length_of_AL ; count05++)
//tempout//      {
//tempout//       ::printf(" %6.2e ",  pc_profilematrix_rep->al[count05] );
//tempout//      }
//tempout//       ::printf(" \n\nAU= ");
//tempout//      for (int count05 = 0 ; count05 < length_of_AL ; count05++)
//tempout//      {
//tempout//       ::printf(" %6.2e ",  pc_profilematrix_rep->au[count05] );
//tempout//      }
//tempout//
//tempout//      ::printf("\n\n=============\n\n");
//tempout//
//tempout//
//tempout    return * this;
//tempout }
//tempout

//##########################################################################
// print jp of profilematrix
void profilematrix::profile_jp_print( void )
 {
   ::printf("\n jp \n\n " );
   for (int count05 = 0 ; count05 < pc_profilematrix_rep->neq ; count05++)
      ::printf("  jp %d ",  pc_profilematrix_rep->jp[count05] );
 }

//##########################################################################
//##########################################################################
// print pc_profilematrix_rep->au--Zhaohui
void profilematrix::profile_al_print( void )
{ ::printf("\n al \n\n " );
 for (int count05 = 0 ; count05 <pc_profilematrix_rep->jp[pc_profilematrix_rep->neq-1] ; count05++)
      {
       ::printf(" al %10.3e ",  pc_profilematrix_rep->al[count05] );
      }
}

//##########################################################################
// print pc_profilematrix_rep->ad--Zhaohui
void profilematrix::profile_ad_print( void )
{ ::printf("\n al \n\n " );
 for (int count05 = 0 ; count05 <pc_profilematrix_rep->neq ; count05++)
      {
       ::printf(" ad %10.3e\n ",  pc_profilematrix_rep->ad[count05] );
      }
}

//##########################################################################
void profilematrix::full_print(char *msg)
  {
    if (*msg) printf("%s\n",msg);
    for ( int row=1 ; row<=dimension_of_profile_M() ; row++ )
      {
        for ( int col=1 ; col<=dimension_of_profile_M() ; col++ )
          {
            printf( "%+8.3e ", full_val(row,col) );
          }
        printf("\n");
      }
  }


//##########################################################################
void profilematrix::error(char * msg1, char * msg2) const
  {
    ::fprintf(stderr,"profmatrix error: %s %s\n", msg1, msg2);
    exit( 1 );
  }




///////////////////////////////////////////////////////////////////////////////


//##############################################################################
// very private part
//##############################################################################
//..   double * profilematrix::data(void) const
//..    {
//..      return this->pc_profilematrix_rep->pd_nDdata;
//..    }
//..
//..    void profilematrix::set_data_pointer(double * data_pointer)
//..    {
//..      this->pc_profilematrix_rep->pd_nDdata = data_pointer;
//..    }
//
//    int profilematrix::rank(void) const
//    {
//      return this->pc_profilematrix_rep->profilematrix_rank;
//    }

//    void profilematrix::rank(int profilematrix_rank)
//    {
//      this->pc_profilematrix_rep->profilematrix_rank = profilematrix_rank;
//    }
//
    long int profilematrix::total_number(void) const
    {
      return this->pc_profilematrix_rep->total_numb;
    }

    void profilematrix::total_number(int number)
    {
      this->pc_profilematrix_rep->total_numb = number;
    }

//    int * profilematrix::dim(void) const
//    {
//      return this->pc_profilematrix_rep->dim;
//    }
//
//    int & profilematrix::get_dim_pointer(void) const
//    {
//      return this->pc_profilematrix_rep->dim[0];
//    }
//
//    void profilematrix::set_dim_pointer(int * dim_pointer)
//    {
//      this->pc_profilematrix_rep->dim = dim_pointer;
//    }
//
//    int profilematrix::dim(int which) const
//    {
//      return this->pc_profilematrix_rep->dim[which-1];
//    }
//
    int profilematrix::reference_count(int up_down)
    {
      this->pc_profilematrix_rep->n += up_down;
      return(this->pc_profilematrix_rep->n);
    }

    void profilematrix::set_reference_count(int ref_count)
    {
      this->pc_profilematrix_rep->n=ref_count;
    }


//##############################################################################
// TENSOR_REP_CC
// #######################
// memory manager part
// overloading operator new in profilematrix::profilematrix_rep class  ##################
void * profilematrix_rep::operator new(size_t s)
  {                                       // see C++ reference manual by
    void *void_pointer;                   // ELLIS and STROUSTRUP page 283.
    void_pointer = ::operator new(s);     // and ECKEL page 529.
//    ::printf("\nnew pointer %p of size %d\n",void_pointer,s);
    if (!void_pointer)
      {
        ::printf("\a\nInsufficient memory\n");
        ::exit(1);
      }
    return void_pointer;
  }

// overloading operator delete in profilematrix::profilematrix_rep class  ##################
void profilematrix_rep::operator delete(void *p)
  {                                       // see C++ reference manual by
                                          // ELLIS and STROUSTRUP page 283.
                                          // and ECKEL page 529.
//    ::printf("deleted pointer %p\n",p);
    ::operator delete(p);
  }



//##############################################################################
//##############################################################################
// DATRI triangularization of profile matrix stored in profile form.
// Both symmetric and unsymmetric profile matrices can be solved!!
// this is important for the nonassociated flow theories in plasticity!
// FORTRAN source code taken from the chapter 15 of FEM IV ed. by
// O. C. Zienkiewicz and R. L. Taylor. See also " Solution of linear
// equations by a profile solver" in Eng. Comp. 1985 Vol 2
// December pages: 344-350.
profilematrix & profilematrix::datri(void)
  {
//  c
//        subroutine datri(al,au,ad,jp,neq,flg)
//        implicit real*8 (a-h,o-z)
//  c.... triangular decomposition of a matrix stored in profile form
//        logical flg
//        integer*2 jp(1)
//        real*8 al(1),au(1),ad(1)
//        common /iofile/ ior,iow
//  c.... n.b.  tol should be set to approximate half-word precision.
//        data zero,one/0.0d0,1.0d0/, tol/0.5d-07/
    double zero = 0.0;
    double one  = 1.0;
    double tol  = sqrt(d_macheps());
//  c.... set initial values for conditioning check
//        dimx = zero
//        dimn = zero
    double dimx = zero;
    double dimn = zero;
//        do 50 j = 1,neq
//        dimn = max(dimn,abs(ad(j)))
//  50    continue
    int j = 1;
    for ( j=1 ; j<=pc_profilematrix_rep->neq ; j++ )
//      dimn = max(dimn,(fabs(pc_profilematrix_rep->ad[j-1])));
//..      dimn =
     if ( dimn < fabs(pc_profilematrix_rep->ad[j-1]) )
        dimn = fabs(pc_profilematrix_rep->ad[j-1]);
//        dfig = zero
    double dfig = zero;
//  c.... loop through the columns to perform the triangular decomposition
//        jd = 1
    int jd = 1;
//        do 200 j = 1,neq

    int jr = 0;
    int jh = 0;
    int is = 0;
    int ie = 0;
    int id = 0;
    int ih = 0;
    int jrh = 0;
    int idh = 0;
    double dd    = 0;
//    double dj    = 0;
//    double ud    = 0;
    double ifig  = 0;
    double daval = 0.0;


    for ( j=1 ; j<=pc_profilematrix_rep->neq ; j++ )
      {
//          jr = jd + 1
        jr = jd + 1;
//          jd = jp(j)
        jd = pc_profilematrix_rep->jp[j-1];
//          jh = jd - jr
        jh = jd -jr;
//          if(jh.gt.0) then
        if ( jh > 0 )
          {
//            is = j - jh
            is = j - jh;
//            ie = j - 1
            ie = j - 1;
//  c.... if diagonal is zero compute a norm for singularity test
//            if(ad(j).eq.zero) call Ndatest(au(jr),jh,daval)
            if ( pc_profilematrix_rep->ad[j-1] == zero )
              daval = datest((pc_profilematrix_rep->ad+jr-1), jh);
//            do 100 i = is,ie
            for ( int i=is ; i<=ie ; i++ )
              {
//              jr = jr + 1
                jr = jr + 1;
//              id = jp(i)
                id = pc_profilematrix_rep->jp[i-1];
//              ih = min(id-jp(i-1),i-is+1)
//                ih = min((id-pc_profilematrix_rep->jp[i-1-1]),(i-is+1));
                ih = id - pc_profilematrix_rep->jp[i-1-1]; // ____Zhaohui
		if (ih >(i-is+1))  ih = i-is+1; 	 // ____Zhaohui

                if ( ih > 0 )
                  {
//                jrh = jr - ih
                    jrh = jr - ih;
//                idh = id - ih + 1
                    idh = id -ih + 1;
//                au(jr) = au(jr) - dot(au(jrh),al(idh),ih)
                    pc_profilematrix_rep->au[jr-1] =
                      pc_profilematrix_rep->au[jr-1] -
                      dot((pc_profilematrix_rep->au+jrh-1),
                          (pc_profilematrix_rep->al+idh-1), ih);
//                if(flg) al(jr) = al(jr) - dot(al(jrh),au(idh),ih)
                    if ( pc_profilematrix_rep->flag == 'N' )
                    pc_profilematrix_rep->al[jr-1] =
                      pc_profilematrix_rep->al[jr-1] -
                      dot((pc_profilematrix_rep->al+jrh-1),
                          (pc_profilematrix_rep->au+idh-1), ih);
//              endif
                  }

//  100       continue
              }
//          endif
          }
//  c.... reduce the diagonal
//          if(jh.ge.0) then
        if ( jh >= 0 )
          {
//            dd = ad(j)
            dd = pc_profilematrix_rep->ad[j-1];
//            jr = jd - jh
            jr = jd - jh;
//            jrh = j - jh - 1
            jrh = j - jh - 1;
//            call dredu(al(jr),au(jr),ad(jrh),jh+1,flg  ,ad(j))
            dredu((pc_profilematrix_rep->al+jr-1),
                  (pc_profilematrix_rep->au+jr-1),
                  (pc_profilematrix_rep->ad+jrh-1),
                  (jh+1),
                  (pc_profilematrix_rep->flag),
                  (pc_profilematrix_rep->ad+j-1) );
//  c.... check for possible errors and print warnings
//            if(abs(ad(j)).lt.tol*abs(dd))  write(iow,2000) j
            if ( (fabs(pc_profilematrix_rep->ad[j-1])<(tol*fabs(dd))) )
::printf("***DATRI WARNING 1***\n Loss of at least N digits in reducing diagonal of equation %d\n",j);
//            if(dd.lt.zero.and.ad(j).gt.zero) write(iow,2001) j
            if ( (dd<zero) && (pc_profilematrix_rep->ad[j-1]>zero) )
::printf("***DATRI WARNING 2***\n Sign of diagonal changed when reducing equation %d\n",j);
//            if(dd.gt.zero.and.ad(j).lt.zero) write(iow,2001) j
            if ( (dd>zero) && (pc_profilematrix_rep->ad[j-1]<zero) )
::printf("***DATRI WARNING 2***\n Sign of diagonal changed when reducing equation %d\n",j);
//            if(ad(j) .eq.  zero)             write(iow,2002) j
            if ( (pc_profilematrix_rep->ad[j-1]==zero) ) // strange it is never == 0 but macheps!!!!!!!!!!!!
::printf("***DATRI WARNING 3***\n Reduced diagonal is zero for equation %d\n",j);
//            if(dd.eq.zero.and.jh.gt.0) then
            if ( (dd==zero) && (jh>0) ) // strange it is never == 0 but macheps!!!!!!!!!!!!
//              if(abs(ad(j)).lt.tol*daval)   write(iow,2003) j
              if ( (fabs(pc_profilematrix_rep->ad[j-1])<(tol*daval)) )
::printf("***DATRI WARNING 4***\n Rank failure for zero unreduced diagonal in equation %d\n",j);
//            endif
//          endif
          }
//  c.... store reciprocal of diagonal, compute condition checks
//          if(ad(j).ne.zero) then
        if ( pc_profilematrix_rep->ad[j-1] != zero )
          {
//            dimx  = max(dimx,abs(ad(j)))
//            dimx = max(dimx,fabs(pc_profilematrix_rep->ad[j-1]));
	      if (dimx < fabs(pc_profilematrix_rep->ad[j-1])) 
	         dimx = fabs(pc_profilematrix_rep->ad[j-1]);
//            dimn  = min(dimn,abs(ad(j)))
//            dimn = min(dimn,fabs(pc_profilematrix_rep->ad[j-1]));
	      if (dimn > fabs(pc_profilematrix_rep->ad[j-1])) 
	         dimn = fabs(pc_profilematrix_rep->ad[j-1]);
//            dfig  = max(dfig,abs(dd/ad(j)))
//            dfig = max(dfig,fabs(dd/pc_profilematrix_rep->ad[j-1]));
	      if (dfig <fabs(dd/pc_profilematrix_rep->ad[j-1])) 
	         dfig=fabs(dd/pc_profilematrix_rep->ad[j-1]); 
//            ad(j) = one/ad(j)
            pc_profilematrix_rep->ad[j-1] = one/pc_profilematrix_rep->ad[j-1];
//          endif
          }
//  200   continue
      }
//  c.... print conditioning information
//        dd = zero
    dd = zero;
//        if(dimn.ne.zero) dd = dimx/dimn
    if ( dimn != zero ) dd = dimx/dimn;
//        ifig = dlog10(dfig) + 0.6
    ifig = log (dfig) + 0.6;
//        write(iow,2004) dimx,dimn,dd,ifig
    ::printf("\n-------- Condition check:\n D-max %f\n D-min %f\n Ratio %f\n Maximum no. diagonal digits lost: %f\n",
              dimx , dimn, dd, ifig);
//        if(ior.lt.0) write(*,2004) dimx,dimn,dd,ifig
//        return
//        end
//

    return *this;
  }



//##############################################################################
//##############################################################################
// DASOL solution of the equations after triangular decomposition
// Both symmetric and unsymmetric profile matrices can be solved!!
// this is important for the nonassociated flow theories in plasticity!
// FORTRAN source code taken from the chapter 15 of FEM IV ed. by
// O. C. Zienkiewicz and R. L. Taylor. See also " Solution of linear
// equations by a profile solver" in Eng. Comp. 1985 Vol 2
// December pages: 344-350.
double * profilematrix::dasol(double * b)
  {
//c
//      subroutine dasol(al,au,ad,b,jp,neq, energy)
//      implicit real*8 (a-h,o-z)
//c.... solution of symmetric equations stored in profile form
//c.... coefficient matrix must be decomposed into its triangular
//c.... factors using datri before using dasol.
//      integer*2 jp(1)
//      real*8 al(1),au(1),ad(1),b(1)
//      common /iofile/ ior,iow
//      data zero/0.0d0/
    double zero = 0.0;
//c.... find the first nonzero entry in the right hand side
//      do 100 is = 1,neq
//        if(b(is).ne.zero) go to 200
//100   continue

    int jr = 0;
    int jh = 0;
    double bd  = 0.0;
    int is = 1;
    for ( is=1 ; is<=pc_profilematrix_rep->neq ; is++ )
      {
        if ( b[is-1] != zero ) break;
      }
//      write(iow,2000)
//      if(ior.lt.0) write(*,2000)
    if ( is > pc_profilematrix_rep->neq )
      {
        ::printf("***DASOL WARNING 1***\n Zero right-hand-side vector\n");
        return b;
      }
//      return
//200   if(is.lt.neq) then
    if ( is < pc_profilematrix_rep->neq )
      {

//c.... reduce the right hand side
//        do 300 j = is+1,neq
        for ( int j=(is+1) ; j<=pc_profilematrix_rep->neq ; j++ )
          {
//          jr = jp(j-1)
            jr = pc_profilematrix_rep->jp[j-1-1];
//          jh = jp(j) - jr
            jh = pc_profilematrix_rep->jp[j-1] - jr;
//          if(jh.gt.0) then
            if ( jh > 0 )
              {
//            b(j) = b(j) - dot(al(jr+1),b(j-jh),jh)
                b[j-1] = b[j-1] -
                         dot((pc_profilematrix_rep->al+jr+1-1),
                             (b+j-jh-1), jh);
              }
//          endif
//300     continue
          }
//      endif
      }
//c.... multiply by inverse of diagonal elements
//      energy = zero
    double energy = zero;
//      do 400 j = is,neq
    int j = is;
    for ( j=is ; j<=pc_profilematrix_rep->neq ; j++ )
      {
//        bd = b(j)
        bd = b[j-1];
//        b(j) = b(j)*ad(j)
        b[j-1] = b[j-1]*pc_profilematrix_rep->ad[j-1];
//        energy = energy + bd*b(j)
        energy = energy + bd*b[j-1];
//400   continue
      }
//c.... backsubstitution
//      if(neq.gt.1) then
    if ( pc_profilematrix_rep->neq > 1 )
      {
//        do 500 j = neq,2,-1
        for ( j=pc_profilematrix_rep->neq ; j>=2 ; j-- )
          {
//          jr = jp(j-1)
            jr = pc_profilematrix_rep->jp[j-1-1];
//          jh = jp(j) - jr
            jh = pc_profilematrix_rep->jp[j-1] - jr;
//          if(jh.gt.0) then
            if ( jh > 0 )
              {
//            call saxpb(au(jr+1),b(j-jh),-b(j),jh, b(j-jh))
                saxpb((pc_profilematrix_rep->au+jr+1-1),
                      (b+j-jh-1),
                      (-(b[j-1])),
                      jh,
                      (b+j-jh-1) );
//          endif
              }
//500     continue
          }
//      endif
      }
//      return
    return b;
//      end

  }


//      subroutine datest(au,jh,daval)
//      implicit real*8 (a-h,o-z)
//      real*8 au(jh)
//c.... test for rank
//      daval = 0.0d0
//      do 100 j = 1,jh
//         daval = daval + abs(au(j))
//100   continue
//      return
//      end
double profilematrix::datest(double* au, int jh)
  {
    double daval = 0.0;
    for ( int j=1 ; j<=jh ; j++ )
      daval += fabs(au[j-1]);
    return daval;
  }

double profilematrix::dot(double* a, double* b, int n)
  {
    double dot = 0.0;
    for ( int i=1 ; i<=n ; i++ )
      dot += a[i-1]*b[i-1];
    return dot;
  }


//..//            call dredu(al(jr),au(jr),ad(jrh),jh+1,flg  ,ad(j))
//..      subroutine dredu(al,au,ad,jh,flg  ,dj)
//..      implicit real*8 (a-h,o-z)
//..c.... reduce diagonal element in triangular decomposition
//..      logical flg
//..      real*8 al(jh),au(jh),ad(jh)
//..c.... computation of column for unsymmetric matrices
//..      if(flg) then
//..        do 100 j = 1,jh
//..          au(j) = au(j)*ad(j)
//..          dj    = dj - al(j)*au(j)
//..          al(j) = al(j)*ad(j)
//..100   continue
//..c.... computation of column for symmetric matrices
//..      else
//..        do 200 j = 1,jh
//..          ud    = au(j)*ad(j)
//..          dj    = dj - au(j)*ud
//..          au(j) = ud
//..200     continue
//..      endif
//..      return
//..      end
void profilematrix::dredu(double* al,
                          double* au,
                          double* ad,
                          int jh,
                          char flag,
                          double* dj)
  {
    double ud = 0.0;
    if ( flag == 'N' )
      {
        for ( int j=1 ; j<=jh ; j++ )
          {
            au[j-1] = au[j-1]*ad[j-1];
            *dj = *dj - al[j-1]*au[j-1];
            al[j-1] = al[j-1]*ad[j-1];
          }
      }
    else
      {
        for ( int j=1 ; j<=jh ; j++ )
          {
            ud = au[j-1]*ad[j-1];
            *dj = *dj - au[j-1]*ud;
            au[j-1] = ud;
          }
      }
  }

//            call saxpb(au(jr+1),b(j-jh),-b(j),jh, b(j-jh))
//c                saxpb(a,       b,      x,     n, c      )
//      subroutine saxpb (a,b,x,n,c)
//      real*8 a(1),b(1),c(1),x
//c... vector times scalar added to second vector
//      do 10 k=1,n
//        c(k) = a(k)*x +b(k)
//   10 continue
//      end
void profilematrix::saxpb(double* a,
                          double* b,
                          double x,
                          int jh,
                          double* c)
  {
   for ( int k=1 ; k<=jh ; k++ )
     c[k-1] = a[k-1] * x + b[k-1];
  }


#endif

