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
//# COPY-YES  (C):     :-))                                                    #
//# PROJECT:           Object Oriented Finite Element Program                  #
//# PURPOSE:                                                                   #
//# CLASS:             skymatrix class                                         #
//#                                                                            #
//# VERSION:                                                                   #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )#
//# TARGET OS:         DOS || UNIX || . . .                                    #
//# PROGRAMMER(S):     Boris Jeremic ( with parts of fortran COLSOL which      #
//#                                    is ported to C++ by BJ                  #
//#                                                                            #
//# DATE:              Nov. 12. 1992.                                          #
//# UPDATE HISTORY:    Nov. 14. 1992.  Improved operator=                      #
//#                    Nov. 16. 1992.  More print options ( upper, lower, full)#
//#                    ___. __. 199_.  derived class from nDarray class        #
//#                    August 22-29 '94 choped to separate files and worked on #
//#                                   const and & issues                       #
//#                    August 30-31 '94 added use_def_dim to full the CC       #
//#                                   resolved problem with temoraries for     #
//#                                   operators + and - ( +=, -= )             #
//#                                                                            #
//#                    September 09 '94 starting to rewrite things after a talk#
//#                    by Stephen Jonson " Objecting the Objects". The point is#
//#                    to forget about inheriting skymatrix from nDarray and   #
//#                    start from separate branch!                             #
//#                    September 11 '94 it works                               #
//#                    September 12-13 '94 looking for the solver for symmetric#
//#                                        and unsymmetric sparse matrices.    #
//#                                        One solution is Taylor's profile    #
//#                                        solver ( see FEM4ed by O.Z. and R.T.#
//#                    September 27 '94 profile solver for symmetric and       #
//#                                     Nonsymmetric systems works!            #
//#                                     (from FEM4ed by O.Z. and R.T.)         #
//#                                                                            #
//##############################################################################
#ifndef SKYMATRIX_CC
#define SKYMATRIX_CC

#include "skymatr.h"

// SYMSKYMATRIX.CC  Symmetric Skyline Sparse Matrix Class



//.. //##########################################################################
//.. // empty constructor ( with default values )
//.. skymatrix::skymatrix( int matrix_order, double init_val)
//.. {
//..  // create the structure:
//..    pc_skymatrix_rep = new skymatrix_rep; // this 'new' is overloaded
//.. 
//.. 
//..    pc_skymatrix_rep->square_dim = matrix_order;
//.. 
//.. // get space for the MAXA vector
//..    pc_skymatrix_rep->p_maxa = new int[matrix_order];
//.. // put all 1 in the MAXA ( somewhat superficial! )
//..    for ( int j=0 ; j<pc_skymatrix_rep->square_dim ; j++ )
//..      pc_skymatrix_rep->p_maxa[j] = 1;
//.. 
//..    pc_skymatrix_rep->total_numb = total_numb;
//.. 
//..  // allocate memory for the actual skymatrix as skymatrix
//..    pc_skymatrix_rep->pd_nDdata = new double [(size_t) pc_skymatrix_rep->total_numb];
//..      if (!pc_skymatrix_rep->pd_nDdata)
//..        {
//..          ::printf("\a\nInsufficient memory for skymatrix_rep\n");
//..          ::exit(1);
//..        }
//.. 
//..    pc_skymatrix_rep->n = 1;  // so far, there's one reference
//.. 
//..    for ( int i=0 ; i<pc_skymatrix_rep->total_numb ; i++ )
//..      pc_skymatrix_rep->pd_nDdata[i] = init_val;
//.. }
//.. 
//tempout//##########################################################################
//tempoutskymatrix::skymatrix(FEModelData & FEMD, Brick3D * b3d, Node *  node)
//tempout//skymatrix::skymatrix(FEModelData & FEMD, Finite_Element & FE)
//tempout  {
//tempout// create the structure:
//tempout     pc_skymatrix_rep = new skymatrix_rep; // this 'new' is overloaded
//tempout
//tempout
//tempout    int Number_of_DOFs = FEMD.get_number_of_DOFs();
//tempout    pc_skymatrix_rep->square_dim = Number_of_DOFs;
//tempout// create ColumnHeight _________________________________________________
//tempout    pc_skymatrix_rep->columnheight = new int [Number_of_DOFs+1];
//tempout      if (!pc_skymatrix_rep->columnheight)
//tempout        {
//tempout          ::printf("\a\nInsufficient memory for pc_skymatrix_rep->columnheight\n");
//tempout          ::exit(1);
//tempout        }
//tempout    for (int count = 1 ; count <= Number_of_DOFs ; count++) 
//tempout
//tempout      //**** Changing 1 to 0 to set initial value of colheight. 
//tempout      pc_skymatrix_rep->columnheight[count]=0;  
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
//tempout       // b3d[el_count].reportLM("LM");
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
//tempout              if (ME > pc_skymatrix_rep->columnheight[II]) 
//tempout                {
//tempout                  pc_skymatrix_rep->columnheight[II] = ME;
//tempout           //  ::printf("ME > ColumnHeight[II] ->> ME= %d  ColumnHeight[%d] = %d \n", ME, II, ColumnHeight[II]);
//tempout           //  getchar();
//tempout                  }
//tempout              }
//tempout         }
//tempout      
//tempout        //::printf("IN colheigth = \n");
//tempout
//tempout    // for (int count03 = 1 ; count03 <= Number_of_DOFs ; count03++)
//tempout    //   {
//tempout    //     ::printf(" %d ",  pc_skymatrix_rep->columnheight[count03] );
//tempout    //   }
//tempout//      getchar();
//tempout      
//tempout      }
//tempout
//tempout
//tempout
//tempout
//tempout        ::printf(" \n ------+------ \n ");
//tempout
//tempout// create MAXA ______________________________________________
//tempout//    pc_skymatrix_rep->maxa = new int [Number_of_DOFs+1];
//tempout    pc_skymatrix_rep->maxa = new int [Number_of_DOFs];
//tempout      if (!pc_skymatrix_rep->maxa)
//tempout        {
//tempout          ::printf("\a\nInsufficient memory for maxa\n");
//tempout          ::exit(1);                                     
//tempout        }
//tempout    for (int count05 = 0 ; count05 < Number_of_DOFs ; count05++) 
//tempout      pc_skymatrix_rep->maxa[count05]=0;
//tempout
//tempout//  For the convience of calculating full_val from skymatrix: maxa[col]-maxa[col-1], adding
//tempout//    pc_skymatrix_rep->maxa[0] = 1;
//tempout
//tempout//    pc_skymatrix_rep->maxa[1] = 1;
//tempout//    pc_skymatrix_rep->maxa[2] = 2;
//tempout    pc_skymatrix_rep->maxa[0] = 1;
//tempout    pc_skymatrix_rep->maxa[1] = 2;
//tempout    //KJBathe pp 1002
//tempout    int MK = 0;
//tempout
//tempout    for (int count07 = 1 ; count07 <=  Number_of_DOFs ; count07++) 
//tempout      {
//tempout        if (pc_skymatrix_rep->columnheight[count07] > MK )
//tempout          {
//tempout            MK = pc_skymatrix_rep->columnheight[count07];
//tempout          }
//tempout
//tempout//        pc_skymatrix_rep->maxa[count07+1] = 
//tempout//           pc_skymatrix_rep->maxa[count07]+pc_skymatrix_rep->columnheight[count07] + 1;
//tempout        pc_skymatrix_rep->maxa[count07] = 
//tempout           pc_skymatrix_rep->maxa[count07-1]+pc_skymatrix_rep->columnheight[count07] + 1;
//tempout      }
//tempout
//tempout   //     ::printf("MAXA = \n");
//tempout   // for (int count = 0 ; count < Number_of_DOFs+1 ; count++)
//tempout   //   {
//tempout   //     ::printf("Col %d MAXA =  %d \n", count, pc_skymatrix_rep->maxa[count] );
//tempout   //   }
//tempout
//tempout
//tempout
//tempout
//tempout
//tempout   int Total_K_length = pc_skymatrix_rep->maxa[Number_of_DOFs+1];
//tempout
//tempout    pc_skymatrix_rep->data = new double [Total_K_length+1];
//tempout
//tempout      if (!pc_skymatrix_rep->data)
//tempout        {
//tempout          ::printf("\a\nInsufficient memory for pc_skymatrix_rep->data\n");
//tempout          ::exit(1);
//tempout        }
//tempout  //ZHaohui Feb. 11, 2000 SOmething wrong with this statememt!
//tempout    for (int count08 = 1 ; count08 <= Number_of_DOFs ; count08++) 
//tempout      pc_skymatrix_rep->data[count08]=0;
//tempout
//tempout
//tempout
//tempout  }
//tempout
//##########################################################################
//##########################################################################
skymatrix::skymatrix(int Number_of_DOFs, int *MAXA, double *initval)
  {
 // create the structure:
    pc_skymatrix_rep = new skymatrix_rep; // this 'new' is overloaded

    pc_skymatrix_rep->square_dim = Number_of_DOFs;

// get space for the MAXA vector
    pc_skymatrix_rep->maxa = new int[Number_of_DOFs+1];
// put all maxa's in the MAXA
    for ( int j=0 ; j<=Number_of_DOFs ; j++ )  
//      pc_skymatrix_rep->maxa[0] = 1;
//    for ( int j=1 ; j<=Number_of_DOFs+1 ; j++ )  
      pc_skymatrix_rep->maxa[j] = MAXA[j];

   int Total_K_length = pc_skymatrix_rep->maxa[Number_of_DOFs];
 // allocate memory for the actual skymatrix as skymatrix
    pc_skymatrix_rep->data = new double [(size_t) Total_K_length];
      if (!pc_skymatrix_rep->data)
        {
          ::printf("\a\nInsufficient memory for skymatrix_rep\n");
          ::exit(1);
        }

    for ( int i=0 ; i<Total_K_length ; i++ )
      pc_skymatrix_rep->data[i] = initval[i];
  }

//outOLD//##############################################################################
//outOLDskymatrix::skymatrix(const skymatrix & x)   // copy initializer
//outOLD {
//outOLD  x.pc_skymatrix_rep->n++; // we're adding another reference.
//outOLD  pc_skymatrix_rep = x.pc_skymatrix_rep;  // point to the new skymatrix_rep.
//outOLD }



//oldDestructor//##########################################################################
//oldDestructorskymatrix::~skymatrix()
//oldDestructor  {
//oldDestructor    if (--pc_skymatrix_rep->n == 0)  // if reference count  goes to 0
//oldDestructor  {
//oldDestructor// DEallocate memory of the actual nDarray
//oldDestructor//    delete [pc_nDarray_rep->pc_nDarray_rep->total_numb] pc_nDarray_rep->pd_nDdata;
//oldDestructor//  see ELLIS & STROUSTRUP $18.3
//oldDestructor//  and note on the p.65($5.3.4)
//oldDestructor//  and the page 276 ($12.4)
//oldDestructor    delete [] pc_skymatrix_rep->pd_nDdata;
//oldDestructor    delete [] pc_skymatrix_rep->p_maxa;
//oldDestructor    delete pc_skymatrix_rep;
//oldDestructor  }
//oldDestructor}

//##########################################################################
skymatrix::~skymatrix()
  {
        ::printf(" ------------  skymatrix::~skymatrix()   \n ");
         delete [] pc_skymatrix_rep->columnheight;
         delete [] pc_skymatrix_rep->maxa;
         delete [] pc_skymatrix_rep->data;
  }
//outOLD////##########################################################################
//outOLDskymatrix & skymatrix::operator=(const skymatrix & rval)
//outOLD  {
//outOLD    rval.pc_skymatrix_rep->n++; // we're adding another reference.
//outOLD//    rval.reference_count(+1);  // tell the rval it has another reference
//outOLD// It is important to increment the reference_counter in the new
//outOLD// tensor before decrementing the reference_counter in the
//outOLD// old tensor_rep to ensure proper operation when assigning a
//outOLD// tensor_rep to itself ( after ARKoenig JOOP May/June '90 )
//outOLD// clean up current value;
//outOLD    if( reference_count(-1) == 0)  // if nobody else is referencing us.
//outOLD      {
//outOLD        delete [] pc_skymatrix_rep->pd_nDdata;
//outOLD        delete [] pc_skymatrix_rep->p_maxa;
//outOLD        delete pc_skymatrix_rep;
//outOLD      }
//outOLD// connect to new value
//outOLD    pc_skymatrix_rep = rval.pc_skymatrix_rep;// point at the rval skymatrix_rep
//outOLD
//outOLD    return *this;
//outOLD  }
//outOLD
//outOLD


////#############################################################################
// Assemble global stiffness matrix from local contributions

//        ::printf("ColHeigth = \n");
//    for (int count = 1 ; count <= Number_of_DOFs ; count++)
//      {
//        ::printf("Col %d height =  %d \n", count, ColumnHeight[count] );
//      }

//##########################################################################
int skymatrix::dimension_of_sky_M(void) const   // dimension of  sky matrix
  {
    return pc_skymatrix_rep->square_dim;
  }

//##########################################################################
int * skymatrix::get_MAXA(void) const  // get pointer to array of
  {                                 // Locations of Diagonals
    return pc_skymatrix_rep->maxa;
  }


//outOLD//##############################################################################
//outOLD// skymatrix addition
//outOLDskymatrix& skymatrix::operator+=(const skymatrix & rval)
//outOLD  {
//outOLD    long int this_total_numb = this->pc_skymatrix_rep->total_numb;
//outOLD    long int rval_total_numb =  rval.pc_skymatrix_rep->total_numb;
//outOLD
//outOLD    if(this_total_numb != rval_total_numb)
//outOLD      {
//outOLD        ::printf("\a\nskymatrixs of different sizes: += not YET possible\n");
//outOLD        ::exit ( 1 );
//outOLD      }
//outOLD
//outOLD// Copy *this if necessary
//outOLD    if ( this->pc_skymatrix_rep->n > 1 )// see ARK in JOOP may/june '90
//outOLD      {                                    // "Letter From a Newcomer"
//outOLD// create the structure:
//outOLD        skymatrix_rep * New_pc_skymatrix_rep = new skymatrix_rep; // this 'new' is overloaded
//outOLD
//outOLD        New_pc_skymatrix_rep->square_dim = this->pc_skymatrix_rep->square_dim;
//outOLD
//outOLD// get space for the MAXA vector
//outOLD        New_pc_skymatrix_rep->p_maxa = new int[this->pc_skymatrix_rep->square_dim];
//outOLD// put all maxa's in the MAXA
//outOLD        for ( int j=0 ; j<this->pc_skymatrix_rep->square_dim ; j++ )
//outOLD          New_pc_skymatrix_rep->p_maxa[j] = this->pc_skymatrix_rep->p_maxa[j];
//outOLD
//outOLD        New_pc_skymatrix_rep->total_numb = this->pc_skymatrix_rep->total_numb;
//outOLD
//outOLD // allocate memory for the actual skymatrix as skymatrix
//outOLD        New_pc_skymatrix_rep->pd_nDdata = new double [(size_t) pc_skymatrix_rep->total_numb];
//outOLD          if (!New_pc_skymatrix_rep->pd_nDdata)
//outOLD            {
//outOLD              ::printf("\a\nInsufficient memory for New_skymatrix_rep\n");
//outOLD              ::exit(1);
//outOLD            }
//outOLD
//outOLD        New_pc_skymatrix_rep->n = 1;  // so far, there's one reference
//outOLD
//outOLD        for ( int i=0 ; i<New_pc_skymatrix_rep->total_numb ; i++ )
//outOLD          New_pc_skymatrix_rep->pd_nDdata[i] = this->pc_skymatrix_rep->pd_nDdata[i];
//outOLD
//outOLD        this->pc_skymatrix_rep->total_numb--;
//outOLD        this->pc_skymatrix_rep = New_pc_skymatrix_rep;
//outOLD//..............................................................................
//outOLD      }
//outOLD// I can add this two skymatrices just as a simple vectors:
//outOLD    for (int k=0 ; k<this->pc_skymatrix_rep->total_numb ; k++)
//outOLD      this->pc_skymatrix_rep->pd_nDdata[k] += rval.pc_skymatrix_rep->pd_nDdata[k];
//outOLD
//outOLD    return *this;
//outOLD  }
//outOLD
//outOLD
//outOLD//##############################################################################
//outOLD// skymatrix addition
//outOLDskymatrix operator+(const skymatrix & lval, const skymatrix & rval)
//outOLD  {
//outOLD    skymatrix result(lval);
//outOLD    result += rval;
//outOLD    return result;
//outOLD  }
//outOLD
//outOLD
//outOLD//##############################################################################
//outOLD// scalar addition
//outOLDskymatrix skymatrix::operator+( double rval)
//outOLD {
//outOLD// construct skymatrix using the same control numbers as for the
//outOLD// original one.
//outOLD      skymatrix add(pc_skymatrix_rep->square_dim,
//outOLD                       pc_skymatrix_rep->p_maxa,
//outOLD                       pc_skymatrix_rep->pd_nDdata);
//outOLD
//outOLD      for ( int i1=0 ; i1<this->pc_skymatrix_rep->total_numb ; i1++ )
//outOLD        {
//outOLD          add.pc_skymatrix_rep->pd_nDdata[i1] += rval;
//outOLD        }
//outOLD    return add;
//outOLD }
//outOLD
//outOLD
//outOLD
//outOLD//##############################################################################
//outOLD// skymatrix substraction
//outOLDskymatrix& skymatrix::operator-=(const skymatrix & rval)
//outOLD  {
//outOLD    long int this_total_numb = this->pc_skymatrix_rep->total_numb;
//outOLD    long int rval_total_numb =  rval.pc_skymatrix_rep->total_numb;
//outOLD
//outOLD    if(this_total_numb != rval_total_numb)
//outOLD      {
//outOLD        ::printf("\a\nskymatrixs of different sizes: -= not YET possible\n");
//outOLD        ::exit ( 1 );
//outOLD      }
//outOLD
//outOLD// Copy *this if necessary
//outOLD    if ( this->pc_skymatrix_rep->n > 1 )// see ARK in JOOP may/june '90
//outOLD      {                                    // "Letter From a Newcomer"
//outOLD//..............................................................................
//outOLD      // create the structure:
//outOLD // create the structure:
//outOLD        skymatrix_rep * New_pc_skymatrix_rep = new skymatrix_rep; // this 'new' is overloaded
//outOLD
//outOLD        New_pc_skymatrix_rep->square_dim = this->pc_skymatrix_rep->square_dim;
//outOLD
//outOLD// get space for the MAXA vector
//outOLD        New_pc_skymatrix_rep->p_maxa = new int[this->pc_skymatrix_rep->square_dim];
//outOLD// put all maxa's in the MAXA
//outOLD        for ( int j=0 ; j<this->pc_skymatrix_rep->square_dim ; j++ )
//outOLD          New_pc_skymatrix_rep->p_maxa[j] = this->pc_skymatrix_rep->p_maxa[j];
//outOLD
//outOLD        New_pc_skymatrix_rep->total_numb = this->pc_skymatrix_rep->total_numb;
//outOLD
//outOLD // allocate memory for the actual skymatrix as skymatrix
//outOLD        New_pc_skymatrix_rep->pd_nDdata = new double [(size_t) pc_skymatrix_rep->total_numb];
//outOLD          if (!New_pc_skymatrix_rep->pd_nDdata)
//outOLD            {
//outOLD              ::printf("\a\nInsufficient memory for New_skymatrix_rep\n");
//outOLD              ::exit(1);
//outOLD            }
//outOLD
//outOLD        New_pc_skymatrix_rep->n = 1;  // so far, there's one reference
//outOLD
//outOLD        for ( int i=0 ; i<New_pc_skymatrix_rep->total_numb ; i++ )
//outOLD          New_pc_skymatrix_rep->pd_nDdata[i] = this->pc_skymatrix_rep->pd_nDdata[i];
//outOLD
//outOLD        this->pc_skymatrix_rep->total_numb--;
//outOLD        this->pc_skymatrix_rep = New_pc_skymatrix_rep;
//outOLD//..............................................................................
//outOLD      }
//outOLD// I can add this two skymatrices just as a simple vectors:
//outOLD    for (int k=0 ; k<this->pc_skymatrix_rep->total_numb ; k++)
//outOLD      this->pc_skymatrix_rep->pd_nDdata[k] -= rval.pc_skymatrix_rep->pd_nDdata[k];
//outOLD
//outOLD    return *this;
//outOLD  }
//outOLD
//outOLD
//outOLD//##############################################################################
//outOLD// skymatrix substraction
//outOLDskymatrix operator-(const skymatrix & lval, const skymatrix & rval)
//outOLD  {
//outOLD    skymatrix result(lval);
//outOLD    result -= rval;
//outOLD    return result;
//outOLD  }
//outOLD
//outOLD
//outOLD//##############################################################################
//outOLD// scalar substraction
//outOLDskymatrix skymatrix::operator-( double rval)
//outOLD {
//outOLD// construct skymatrix using the same control numbers as for the
//outOLD// original one.
//outOLD      skymatrix substract(pc_skymatrix_rep->square_dim,
//outOLD                             pc_skymatrix_rep->p_maxa,
//outOLD                             pc_skymatrix_rep->pd_nDdata);
//outOLD
//outOLD      for ( int i1=0 ; i1<this->pc_skymatrix_rep->total_numb ; i1++ )
//outOLD        {
//outOLD          substract.pc_skymatrix_rep->pd_nDdata[i1] -= rval;
//outOLD        }
//outOLD    return substract;
//outOLD }
//outOLD
//outOLD
//outOLD//##############################################################################
//outOLD// scalar multiplication
//outOLDskymatrix skymatrix::operator*( double rval)
//outOLD {
//outOLD// construct skymatrix using the same control numbers as for the
//outOLD// original one.
//outOLD      skymatrix mult(pc_skymatrix_rep->square_dim,
//outOLD                        pc_skymatrix_rep->p_maxa,
//outOLD                        pc_skymatrix_rep->pd_nDdata);
//outOLD
//outOLD      for ( int i1=0 ; i1<this->pc_skymatrix_rep->total_numb ; i1++ )
//outOLD        {
//outOLD          mult.pc_skymatrix_rep->pd_nDdata[i1] *= rval;
//outOLD        }
//outOLD    return mult;
//outOLD }
//outOLD
//outOLD
//outOLD
//outOLD
//outOLD
//outOLD
//outOLD
//outOLD
//outOLD
//outOLD//##########################################################################
//outOLD// this one is the same as mval except that it is more convinient
//outOLD// to overload operator (row,col).
//outOLDdouble & skymatrix::operator( )(int row, int col) const
//outOLD  {
//outOLD    if( row>col )         // Now we can call the matrix
//outOLD      {                   // as if it is full matrix.
//outOLD        int temp = row;   // This makes small overhead
//outOLD        row = col;        // and it will be removed latter on.
//outOLD        col = temp;
//outOLD      }
//outOLD    int diagonal_member_n = *(pc_skymatrix_rep->p_maxa+col-1);  //-1:starts from 0
//outOLD    int member_of_sky_n = diagonal_member_n + col - row -1; //-1:starts from 0
//outOLD    double * member_of_sky = &(pc_skymatrix_rep->pd_nDdata[member_of_sky_n]) ;
//outOLD    return( * member_of_sky );
//outOLD  }

//##########################################################################
double & skymatrix::mval(int row, int col) const
  {
    if( row>col )         // Now we can call the matrix
      {                   // as if it is full matrix.
        int temp = row;   // This makes small overhead
        row = col;        // and it will be removed latter on.
        col = temp;
      }
    int diagonal_member_n = *(pc_skymatrix_rep->maxa+col-1);  //-1:starts from 0
   int member_of_sky_n = diagonal_member_n + col - row -1; //-1:starts from 0
    double * member_of_sky = &(pc_skymatrix_rep->data[member_of_sky_n]) ;
    return( * member_of_sky );
// put this latter after debugging
//          return(*( ps_sky_m_rep->pd_nDdata[*(ps_sky_m_rep->p_maxa+col)+col-row-1]));
  }



//##########################################################################
// full_val function allows you to treat skymatrix as if it is full
// float matrix. The function will calculate position inside sky matrix
// and return appropriate number if row and col are bellow skyline or
// return zero (0) if row and col are above sky line
double skymatrix::full_val(int row, int col) const
  {
    if( row>col )
      {
        int temp = row;
        row = col;
        col = temp;
      }
  int total_column_height = col;
  int actual_column_position = total_column_height - row + 1;
  int real_column_height = pc_skymatrix_rep->maxa[col]-pc_skymatrix_rep->maxa[col-1]; 
  int how_much_above_skyline=actual_column_position-real_column_height; //  adding -1 -Zhaohui
  if ( how_much_above_skyline > 0 )
    {
      double back = 0.0;    ////    // not very smart ???????
      return ( back );
    }
  else
    {
      int diagonal_member_n = *(pc_skymatrix_rep->maxa+col-1);  //-1:starts from 0 // deleting -1 --Zhaohui
      int member_of_sky_n = diagonal_member_n + col - row-1 ; //-1:starts from 0 // deleting -1  --Zhaohui
      double member_of_sky = pc_skymatrix_rep->data[member_of_sky_n] ;
      return( member_of_sky );
// put this latter after debugging
//        return(*( ps_sky_m_rep->pd_nDdata[*(ps_sky_m_rep->p_maxa+col)+col-row-1]));
    }
  }
     // used by skymatrix functions which KNOW they aren't
     // exceeding the boundaries





//##########################################################################
double & skymatrix::val(int row, int col)
  {
//    double zero=0.0;
    if ( row<=0 && row>=dimension_of_sky_M() && col>=dimension_of_sky_M() )
      {
        error("index out of range");
//        return zero;
      }
//    else
      return (mval(row,col));
  }

//##########################################################################
double skymatrix::cval(int row, int col) const
  {
    double zero=0.0;
    if ( row<=0 && row>=dimension_of_sky_M() && col>=dimension_of_sky_M() )
      {
        error("index out of range");
        return zero;
      }
//    else
      return (mval(row,col));
  }


//##########################################################################
double skymatrix::mmin()
  {
    double temp = 0.0;
    if ( dimension_of_sky_M()<=0 )
      {
        error("bad skymatrix size for min ()");
        return 0.0;
      }
    double minimum = mval(1,1);
    for ( int row=1 ; row<=dimension_of_sky_M() ; row++ )
      for ( int col=1 ; col<=dimension_of_sky_M() ; col++ )
        if ( (temp=mval(row,col)) < minimum )
          minimum = temp;
    return minimum;
  }

//##########################################################################
double skymatrix::mmax()
  {
    double temp = 0.0;
    if( dimension_of_sky_M()<=0 )
      {
        error("bad skymatrix size for max()");
        double zero=0.0;
        return zero;
      }
    double maximum = mval(1,1);
    for ( int row=1 ; row<=dimension_of_sky_M() ; row++ )
      for ( int col=1 ; col<=dimension_of_sky_M() ; col++ )
        {
          if ( (temp=mval(row,col)) > maximum )
          maximum = temp;
        }
    return maximum;
  }

//outOLD//##########################################################################
//outOLDdouble skymatrix::mean()
//outOLD  {
//outOLD    int col = 0;
//outOLD    int row = 1;
//outOLD    double sum = 0;
//outOLD    for ( row=1 ; row<=dimension_of_sky_M() ; row++ )
//outOLD      for ( col=1 ; col<=dimension_of_sky_M() ; col++ )
//outOLD        sum += fabs(mval(row,col));
//outOLD    return sum/(row*col);
//outOLD  }
//outOLD
//outOLD
//##########################################################################
void skymatrix::lower_print(char *msg)
  {
    if (*msg) printf("%s\n",msg);
    for ( int row=1 ; row<=dimension_of_sky_M() ; row++ )
      {
        int total_column_height = row;
        int real_column_height = pc_skymatrix_rep->maxa[row] - pc_skymatrix_rep->maxa[row-1];
        int numb_of_voids = total_column_height - real_column_height;
        int n_of_voids_to_reach_number = numb_of_voids;
        for ( int col=1 ; col<=row ; col++ )
          {
            if( n_of_voids_to_reach_number > 0 )
              {
                 for(int void_count=1; void_count<=numb_of_voids; void_count++)
                  {
                    printf("********* ");
                    col++;
                    n_of_voids_to_reach_number--;
                  }
              }
            printf( "%+6.2e ", cval(row,col) );
          }
        printf("\n");
      }
  }

//##########################################################################
void skymatrix::upper_print(char *msg)
  {
    if (*msg) printf("%s\n",msg);
    for ( int row=1 ; row<=dimension_of_sky_M() ; row++ )
      {
        for ( int voids=0 ; voids<row-1 ; voids++ ) printf("        ");
        for ( int col=row ; col<=dimension_of_sky_M() ; col++ )
          {
            int total_column_height = col;
            int actual_column_position = total_column_height - row + 1;
            int real_column_height = pc_skymatrix_rep->maxa[col] - pc_skymatrix_rep->maxa[col-1];
            int how_much_above_skyline=actual_column_position-real_column_height;
            if ( how_much_above_skyline > 0 )
              {
                printf("********* ");
              }
            else
              {
                printf( "%+6.2e ", cval(col,row) );
              }
          }
        printf("\n");
      }
  }


//##########################################################################
void skymatrix::full_print(char *msg)
  {
    if (*msg) printf("%s\n",msg);
    for ( int row=1 ; row<=dimension_of_sky_M() ; row++ )
      {
        for ( int col=1 ; col<=dimension_of_sky_M() ; col++ )
          {
            printf( "%+6.2e ", full_val(row,col) );
          }
        printf("\n");
      }
  }


//##########################################################################
void skymatrix::error(char * msg1, char * msg2) const
  {
    ::fprintf(stderr,"skymatrix error: %s %s\n", msg1, msg2);
    exit( 1 );
  }



//outOLD


//tempout//#############################################################################
//tempout// Assemble global stiffness matrix from local contributions
//tempoutskymatrix & skymatrix::AssembleBricksInSkyMatrix(stiffness_matrix & Ke,
//tempout                                                 Brick3D          * b3d,
//tempout                                                 Node             * node
//tempout                                                )
//tempout  {
//tempout//    int BrickNumber = b3d->get_Brick_Number();
//tempout//    b3d->reportshort("");
//tempout
//tempout// for element numbered BrickNumber create LM array (see Bathe pp984
//tempout//    for (int LocalNodeNumber = 1 ; LocalNodeNumber<=20 ; LocalNodeNumber++ )
//tempout
//tempout
//tempout//--    for (int LocalNodeNumber = 1 ; LocalNodeNumber<=8 ; LocalNodeNumber++ )// for 8noded brick
//tempout//--      {
//tempout//--//        int global_node_number = b3d[BrickNumber-1].get_global_number_of_node(LocalNodeNumber-1);
//tempout//--        int global_node_number = b3d->get_global_number_of_node(LocalNodeNumber-1);
//tempout//--        LM[3*LocalNodeNumber-3] = node[global_node_number].eqn_tx();
//tempout//--        LM[3*LocalNodeNumber-2] = node[global_node_number].eqn_ty();
//tempout//--        LM[3*LocalNodeNumber-1] = node[global_node_number].eqn_tz();
//tempout//--      }
//tempout
//tempout       ::printf("\n\n");
//tempout//..for (int count01=1;count01<=8;count01++)
//tempout//..  {
//tempout//..::printf("element %4d localNode %4d Globalnode %4d  LM   %4d   %4d   %4d\n", BrickNumber, count01,b3d->get_global_number_of_node(count01-1),  LM[count01*3-3], LM[count01*3-2], LM[count01*3-1] );
//tempout//..  }
//tempout
//tempout        int *lm = b3d->get_LM();
//tempout
//tempout
//tempout// KJBathe pp1003 ADDBAN
//tempout
//tempout    int DOF_count_for_Brick = 24;
//tempout    int ND = DOF_count_for_Brick;
//tempout    int NDI = 0;
//tempout    int II = 0;
//tempout    int MI = 0;
//tempout    int KS = 0;
//tempout    int KSS = 0;
//tempout    int JJ = 0;
//tempout    int IJ = 0;
//tempout    int KK = 0;
//tempout    for (int I = 1 ; I <= DOF_count_for_Brick ; I++)
//tempout      {
//tempout        II =  lm[I-1];
//tempout        if (II > 0) 
//tempout          {
//tempout            MI = pc_skymatrix_rep->maxa[II-1];   //II-1--> II
//tempout
//tempout            KS = I;
//tempout            for (int J = 1 ; J <= DOF_count_for_Brick ; J++)
//tempout              {
//tempout                JJ = lm[J-1];
//tempout                if (JJ > 0)
//tempout                  {
//tempout                    IJ = II - JJ;
//tempout                    if (IJ >= 0 ) 
//tempout                      { 
//tempout                        KK = MI + IJ;
//tempout                        KSS = KS;
//tempout                        if ( J >= I ) 
//tempout                          {
//tempout                            KSS = J + NDI;
//tempout                          }
//tempout                        pc_skymatrix_rep->data[KK-1] = 
//tempout                          pc_skymatrix_rep->data[KK-1] + Ke.val(I,J); // J-->KSS
//tempout                      }
//tempout                  }
//tempout                KS = KS + ND -J;
//tempout
//tempout              }
//tempout
//tempout          }
//tempout        NDI = NDI + ND - I;
//tempout      }
//tempout
//tempout
//tempout    return * this;
//tempout  }
//tempout
//tempout


///****************************************************************************
//* COPY-YES  (C):   1990,1991,1992
//* PROJECT:
//* FILE:
//* FUNCTION:
//* PURPOSE:
//* VERSION
//* LANGUAGE:        Microsoft C 6.0 , Borland C++ 3.1
//* TARGET OS:       DOS
//* PROGRAMMER:      Jeremic Boris
//* DATE:
//* UPDATE HISTORY:  Nov. 12. 1992. C++ ver
//****************************************************************************/
///*....................................................................*/
///*.                                                                  .*/
///*.                    C  O  L  S  O  L                              .*/
///*.                                                                  .*/
///*. Program to solve Finite Element Static Equilibrium Equations     .*/
///*. in Core, using compacted and column reduction scheme             .*/
///*.                                                                  .*/
///*.                                                                  .*/
///*. Input Variables                                                  .*/
///*.                                                                  .*/
///*.       a[nwk]      = stiffness matrix stored in compacted form    .*/
///*.       v[nn]       = right_hand_side load vector                  .*/
///*.       maxa[nn+1]  = vector containing addresses of diagonal      .*/
///*.                     elements of stiffness matrix in a            .*/
///*.       nn          = number of equations                          .*/
///*.       nwk         = number of elements below skyline of matrix   .*/
///*.                                                                  .*/
///*. Output                                                           .*/
///*.                                                                  .*/
///*.       a[nwk]      = D and L - factors of stiffness matrix        .*/
///*.       v[nn]       = displacement vector                          .*/
///*.                                                                  .*/
///*....................................................................*/
///*....................................................................*/
///*.      Fortran source code was taken from the book:                .*/
///*.                                                                  .*/
///*.      Klaus-Jürgen Bathe ;                                        .*/
///*.                                                                  .*/
///*.          Finite Element Procedures In Engineering Analysis       .*/
///*.                                                                  .*/
///*....................................................................*/
///*.      Rearanged for C language by Jeremic Boris                   .*/
///*.      Bekhme dam site, August & September 1990.                   .*/
///*.      Beograd, home,   January & February 1991.                   .*/
///*....................................................................*/
///*.      Rearanged for C++ language by Jeremic Boris                 .*/
///*.      Boulder, home,   Nov. 12. 1992.                             .*/
///*....................................................................*/
///*....................................................................*/
///* void function v_ldl_factorize will                                 */
///* factorize matrix a in L_D_L factors                                */
///* arguments are : a    - pointer (double) to vector of name a */
///*                 maxa - pointer (int) to vector maxa         */
///*                 nn   - pointer (int) to number nn           */
///*                                                                    */
///* 03. august.   1990. Bekhme Dam Site          0. revision           */
///* 04. september 1990. Bekhme Dam Site          1. revision           */
///* 14. september 1990. Bekhme Dam Site          2. revision           */
///* 01. february  1991. E.P. Beograd             3. revision           */
///* 12. November  1992. Boulder CU ( home)       4. revision           */
///* 28. May       1999. Clarkson, Potsdam)       5. revision           */
///*....................................................................*/
skymatrix & skymatrix::v_ldl_factorize()
{
  int kn  = 0;
  int kl  = 0;
  int ku  = 0;
  int kh  = 0;
  int k   = 0;
  int ic  = 0;
  int klt = 0;
  int ki  = 0;
  int nd  = 0;
  int kk  = 0;
  int n   = 0;
  int j   = 0;
  int l   = 0;
  double c = 0.0;
  double b = 0.0;
  ::printf(" \n\n* * * Equations to factorize : ");
  for ( n=1 ; n<=pc_skymatrix_rep->square_dim ; n++ )
    {
      printf(" %5d\b\b\b\b\b", pc_skymatrix_rep->square_dim - n);
      kn=*(pc_skymatrix_rep->maxa-1+n);
      kl=kn+1;
      ku=*(pc_skymatrix_rep->maxa-1+n+1)-1;
      kh=ku-kl;         // changes ######## from colsol.c
      if ( kh>0 )       // *(pd_ldl_a. . . ) --> *(pc_skymatrix_rep->pd_nDdata-1. . . )
        {               // *(pi_ldl_maxa. . . ) --> *(pc_skymatrix_rep->p_maxa-1. . .)
          k=n-kh;       // *(pi_ldl_nn) --> pc_skymatrix_rep->square_dim
          ic=0;
          klt=ku;
          for ( j=1 ; j<=kh ; j++ )
            {
              ic=ic+1;
              klt=klt-1;
              ki=*(pc_skymatrix_rep->maxa-1+k);
              nd=*(pc_skymatrix_rep->maxa-1+k+1)-ki-1;
              if ( nd>0 )
                {
                  kk=( (ic<nd) ? ic : nd );
                  c=0.0;
                  for ( l=1 ; l<=kk ; l++ )
                    c=c+(*(pc_skymatrix_rep->data-1+ki+l))*(*(pc_skymatrix_rep->data-1+klt+l));
                  *(pc_skymatrix_rep->data-1+klt)=*(pc_skymatrix_rep->data-1+klt)-c;
                }
              k=k+1;
            }
        }
      if ( kh>=0 )
        {
          k=n;
          b=0.0;
          for ( kk=kl ; kk<=ku ; kk++ )
            {
              k=k-1;
              ki=*(pc_skymatrix_rep->maxa-1+k);
              c=(*(pc_skymatrix_rep->data-1+kk))/(*(pc_skymatrix_rep->data-1+ki));
              b=b+c*(*(pc_skymatrix_rep->data-1+kk));
              *(pc_skymatrix_rep->data-1+kk)=c;
            }
          *(pc_skymatrix_rep->data-1+kn)=*(pc_skymatrix_rep->data-1+kn)-b;
        }
      if ( *(pc_skymatrix_rep->data-1+kn)<=0 )
        {
          printf("\n Colsol Stoped - Stiffness Matrix not positive definite \n");
          printf(" non positive pivot for equation, %d\n ", n);
          printf(" pivot, %.12e \n", *(pc_skymatrix_rep->data-1+kn) );
          exit(1);
        }
//  printf("--------------------------  %d\n",n);
//  for( int i=0 ; i<=11 ; i++ )
//    {
//      printf("pc_skymatrix_rep->pd_nDdata[%d] = %8.4f\n",i, pc_skymatrix_rep->pd_nDdata[i]);
//    }
//  getch();


  }
  printf("\n");
  return(*this);
}

///****************************************************************************
//* COPY-YES  (C):   1990,1991   Jeremic Boris,
//* PROJECT:
//* FILE:
//* FUNCTION:
//* PURPOSE:
//* VERSION
//* LANGUAGE:        Microsoft C 6.0
//* TARGET OS:       DOS
//* PROGRAMMER:      Jeremic Boris
//* DATE:
//* UPDATE HISTORY:
//****************************************************************************/
///*....................................................................*/
///* double function  d_reduce_r_h_s_l_v                                  */
///* will reduce right hand side load vector v                          */
///* arguments are : a    - pointer (double) to vector of name a */
///*                 v    - pointer (double) to vector of name v */
///*                 maxa - pointer (int) to vector maxa         */
///*                 nn   - pointer (int) to number nn           */
///*                                                                    */
///* 10. august.   1990. Bekhme Dam Site      0. revision               */
///* 04. september 1990. Bekhme Dam Site      1. revision               */
///* 15. september 1990. Bekhme Dam Site      2. revision               */
///* 01. february  1991. E.P. Beograd         3. revision               */
///* 12. November  1992. Boulder CU ( home)   4. revision               */
///*....................................................................*/
double * skymatrix::d_reduce_r_h_s_l_v ( double *pd_rhs )
{
  int kl = 0;
  int ku = 0;
  int k  = 0;
  int n  = 0;
  int kk = 0;
  double c = 0.0;
  printf("\n * * * Right hand side loads to reduce :");
  for ( n=1 ; n<=pc_skymatrix_rep->square_dim ; n++ )
    {
      printf(" %5d\b\b\b\b\b", pc_skymatrix_rep->square_dim - n);
      kl=*(pc_skymatrix_rep->maxa-1+n)+1;
      ku=*(pc_skymatrix_rep->maxa-1+n+1)-1;  // changes ######## from colsol.c
      if ( ku>=kl )             // *(pd_red_a. . . ) --> *(pc_skymatrix_rep->pd_nDdata-1. . . )
        {                       // *(pi_red_maxa. . . ) --> *(pc_skymatrix_rep->p_maxa-1. . .)
          k=n;                  // *(pi_red_nn) --> pc_skymatrix_rep->square_dim
          c=0.0;                // *(pd_rhs. . .) --> *(pd_rhs-1 . . . )
          for ( kk=kl ; kk<=ku ; kk++ )
            {
              k=k-1;
              c=c+(*(pc_skymatrix_rep->data-1+kk))*(*(pd_rhs-1+k));
            }
          *(pd_rhs-1+n)=*(pd_rhs-1+n)-c;   //pd_rhs[-1+n]==*(pd_rhs-1+n)
        }
    }
  printf("\n");
  return( pd_rhs );
}

///****************************************************************************
//* COPY-YES  (C):   1990,1991     Jeremic Boris,
//* PROJECT:
//* FILE:
//* FUNCTION:
//* PURPOSE:
//* VERSION
//* LANGUAGE:        MicroSoft C 6.0
//* TARGET OS:       DOS
//* PROGRAMMER:      Jeremic Boris  (J_P_S_B)
//* DATE:
//* UPDATE HISTORY:
//*
//*
//*
//****************************************************************************/
///*....................................................................*/
///* void function v_back_substitute will                               */
///* back substitute factorized matrix a and solve for unknown vector v */
///* arguments are : pd_bac_a    - pointer (double) to vector of name a */
///*                 pd_bac_v    - pointer (double) to vector of name v */
///*                 pi_bac_maxa - pointer (int) to vector maxa         */
///*                 pi_bac_nn   - pointer (int) to number nn           */
///*                                                                    */
///* 10. august.   1990. Bekhme dam site     0.  revision               */
///* 04. september 1990. Bekhme Dam Site     1.  revision               */
///* 17. septembar 1990. Bekhme Dam Site     2.  revision               */
///* 01. february  1991. E.P. Beograd        3.  revision               */
///* 12. November  1992. Boulder CU ( home)  4. revision                */
///*....................................................................*/
double * skymatrix::d_back_substitute ( double *pd_rhs )
{
  int k  = 0;
  int n  = 0;
  int kl = 0;
  int ku = 0;
  int l  = 0;
  int kk = 0;
  printf(" \n * * * Equations to back substitute :");
  for ( n=1 ; n<=pc_skymatrix_rep->square_dim ; n++ )
    {
      printf(" %5d\b\b\b\b\b", pc_skymatrix_rep->square_dim - n);
      k=*(pc_skymatrix_rep->maxa-1+n);
      *(pd_rhs-1+n)=(*(pd_rhs-1+n))/(*(pc_skymatrix_rep->data-1+k));
    }
  if ( pc_skymatrix_rep->square_dim==1 ) return(pd_rhs);
  n=pc_skymatrix_rep->square_dim ;
  for ( l=2 ; l<=pc_skymatrix_rep->square_dim ; l++ )
    {                            // changes ######## from colsol.c
      kl=*(pc_skymatrix_rep->maxa-1+n)+1;     // *(pd_bac_a. . . ) --> *(pc_skymatrix_rep->pd_nDdata-1. . . )
      ku=*(pc_skymatrix_rep->maxa-1+n+1)-1;   // *(pi_bac_maxa. . . ) --> *(pc_skymatrix_rep->p_maxa-1. . .)
      if( ku>=kl )               // *(pi_bac_nn) --> pc_skymatrix_rep->square_dim
        {                        // *(pd_rhs. . .) --> *(pd_rhs-1 . . . )
          k=n;
          for ( kk=kl ; kk<=ku ; kk++ )
            {
              k=k-1;
              *(pd_rhs-1+k)=*(pd_rhs-1+k)-(*(pc_skymatrix_rep->data-1+kk))*(*(pd_rhs-1+n));
            }
          n=n-1;
        }
    }
  printf("\n");
  return (pd_rhs);
}

///////////////////////////////////////////////////////////////////////////////


//outOLD//##############################################################################
//outOLD// very private part
//outOLD//##############################################################################
//outOLD   double * skymatrix::data(void) const
//outOLD    {
//outOLD      return this->pc_skymatrix_rep->pd_nDdata;
//outOLD    }
//outOLD
//outOLD    void skymatrix::set_data_pointer(double * data_pointer)
//outOLD    {
//outOLD      this->pc_skymatrix_rep->pd_nDdata = data_pointer;
//outOLD    }
//outOLD//
//outOLD//    int skymatrix::rank(void) const
//outOLD//    {
//outOLD//      return this->pc_skymatrix_rep->skymatrix_rank;
//outOLD//    }
//outOLD
//outOLD//    void skymatrix::rank(int skymatrix_rank)
//outOLD//    {
//outOLD//      this->pc_skymatrix_rep->skymatrix_rank = skymatrix_rank;
//outOLD//    }
//outOLD//
//outOLD    long int skymatrix::total_number(void) const
//outOLD    {
//outOLD      return this->pc_skymatrix_rep->total_numb;
//outOLD    }
//outOLD
//outOLD    void skymatrix::total_number(long int number)
//outOLD    {
//outOLD      this->pc_skymatrix_rep->total_numb = number;
//outOLD    }
//outOLD
//outOLD//    int * skymatrix::dim(void) const
//outOLD//    {
//outOLD//      return this->pc_skymatrix_rep->dim;
//outOLD//    }
//outOLD//
//outOLD//    int & skymatrix::get_dim_pointer(void) const
//outOLD//    {
//outOLD//      return this->pc_skymatrix_rep->dim[0];
//outOLD//    }
//outOLD//
//outOLD//    void skymatrix::set_dim_pointer(int * dim_pointer)
//outOLD//    {
//outOLD//      this->pc_skymatrix_rep->dim = dim_pointer;
//outOLD//    }
//outOLD//
//outOLD//    int skymatrix::dim(int which) const
//outOLD//    {
//outOLD//      return this->pc_skymatrix_rep->dim[which-1];
//outOLD//    }
//outOLD//
//outOLD    int skymatrix::reference_count(int up_down)
//outOLD    {
//outOLD      this->pc_skymatrix_rep->n += up_down;
//outOLD      return(this->pc_skymatrix_rep->n);
//outOLD    }
//outOLD
//outOLD    void skymatrix::set_reference_count(int ref_count)
//outOLD    {
//outOLD      this->pc_skymatrix_rep->n=ref_count;
//outOLD    }
//outOLD
//outOLD
//outOLD//##############################################################################
//outOLD// TENSOR_REP_CC
//outOLD// #######################
//outOLD// memory manager part
//outOLD// overloading operator new in skymatrix::skymatrix_rep class  ##################
//outOLDvoid * skymatrix_rep::operator new(size_t s)
//outOLD  {                                       // see C++ reference manual by
//outOLD    void *void_pointer;                   // ELLIS and STROUSTRUP page 283.
//outOLD    void_pointer = ::operator new(s);     // and ECKEL page 529.
//outOLD//    ::printf("\nnew pointer %p of size %d\n",void_pointer,s);
//outOLD    if (!void_pointer)
//outOLD      {
//outOLD        ::printf("\a\nInsufficient memory\n");
//outOLD        ::exit(1);
//outOLD      }
//outOLD    return void_pointer;
//outOLD  }
//outOLD
//outOLD// overloading operator delete in skymatrix::skymatrix_rep class  ##################
//outOLDvoid skymatrix_rep::operator delete(void *p)
//outOLD  {                                       // see C++ reference manual by
//outOLD                                          // ELLIS and STROUSTRUP page 283.
//outOLD                                          // and ECKEL page 529.
//outOLD//    ::printf("deleted pointer %p\n",p);
//outOLD    ::operator delete(p);
//outOLD  }
//outOLD

#endif

