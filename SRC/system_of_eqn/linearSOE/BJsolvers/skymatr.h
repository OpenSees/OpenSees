//#############################################################################
//                                                                            #
//                                                                            #
//             /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/~~\                #
//            |                                          |____|               #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |        B A S E   C L A S S E S           |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |          C + +     H E A D E R           |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//         /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/   |                    #
//         \_________________________________________\__/                     #
//                                                                            #
//                                                                            #
//#############################################################################
//#############################################################################
//##############################################################################
//# COPYRIGHT (C):     :-))                                                    #
//# PROJECT:           Object Oriented Finite Element Program                  #
//# PURPOSE:                                                                   #
//# CLASS:             skymatrix                                               #
//#                                                                            #
//# VERSION:                                                                   #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.10, SUN C++ ver=2.1 )#
//# TARGET OS:         DOS || UNIX || . . .                                    #
//# PROGRAMMER(S):     Boris Jeremic                                           #
//#                                                                            #
//#                                                                            #
//# DATE:              November '92                                            #
//# UPDATE HISTORY:    05 - __ avgust '93.  redefined as derived class from    #
//#                                 nDarray class                              #
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
//#                                                                            #
//#                    September 11 '94 it works                               #
//#                    September 12-13 '94 looking for the solver for symmetric#
//#                                        and unsymmetric sparse matrices.    #
//#                                        One solution is Taylor's profile    #
//#                                        solver ( see FEM4ed by O.Z. and R.T.#
//#                    September 27 '94 profile solver for symmetric and       #
//#                                     Nonsymmetric systems works!            #
//#                                     (from FEM4ed by O.Z. and R.T.)         #
//#                                                                            #
//#                                                                            #
//##############################################################################

#ifndef SKYMATRIX_HH
#define SKYMATRIX_HH

#include <basics.h>


//#include <femdata.h>
//#include <brick3d.h>
//#include <node.h>
//#include <stifmat.h>


//outOLD// SKYMATRIX_CC Skyline Sparse Matrix Class
//outOLD// this one is based on Bathe's book!
//outOLDclass skymatrix_rep
//outOLD  {
//outOLD    public:
//outOLD      friend class skymatrix;
//outOLD    private:
//outOLD      double *pd_nDdata;   // skymatrix as 1D array
//outOLD
//outOLD      int *p_maxa;         // array of diagonal location pointers (MAXA)
//outOLD
//outOLD      long int total_numb;      // total number of elements in skymatrix
//outOLD
//outOLD      int square_dim;      // dimension of a system ( as square matrix )
//outOLD
//outOLD      int n;               // reference count
//outOLD
//outOLD    public:
//outOLD// overloading operator new and delete in skymatrix_rep class  ########
//outOLD      void * operator new(size_t s); // see C++ reference manual by
//outOLD      void operator delete(void *);  // by ELLIS and STROUSTRUP page 283.
//outOLD                                   // and ECKEL page 529.
//outOLD  };


class skymatrix
  {
//    private:
      struct skymatrix_rep
        {
          int * columnheight;
          int * maxa;
          double * data;
          int square_dim;
        } *pc_skymatrix_rep;
//outOLD      skymatrix_rep * pc_skymatrix_rep;
    public:
//..       skymatrix(int matrix_order=1, double init_value=0.0 ); // default constructor
//tempout      skymatrix(FEModelData & FEMD );
      skymatrix(int order_n, int *maxa, double *initval);
//outOLD//      skymatrix(char *flag, int dimension ); // create an ident skymatrix
//outOLD      skymatrix(const skymatrix & x); // copy-initializer
//tempout      skymatrix(FEModelData & FEMD, Brick3D * b3d, Node * node);
//      skymatrix(FEModelData & FEMD, Finite_Element & );
      ~skymatrix();


//outOLD//      void ColumnHeights(Brick3D * b3d, Node * node, FEModelData & FEMD);
//outOLD//      void create_MAXA(FEModelData  & );
//outOLD
      int dimension_of_sky_M(void ) const; // dimension of  sky matrix
      int *get_MAXA(void) const; // get pointer to array of
//outOLD                              // Locations of Diagonals
//outOLD//outOLD      skymatrix & operator=(const skymatrix & rval); // skymatrix assignment
//outOLD
//outOLD//....// This is from JOOP May/June 1990 after ARKoenig
//outOLD      skymatrix& operator +=( const skymatrix & ); // skymatrix addition
//outOLD      friend skymatrix operator+(const skymatrix & , const skymatrix & ); // skymatrix addition
//outOLD//....// This is from JOOP May/June 1990 after ARKoenig
//outOLD      skymatrix& operator -=( const skymatrix & ); // skymatrix subtraction
//outOLD      friend skymatrix operator-(const skymatrix & , const skymatrix & ); // skymatrix subtraction
//outOLD
//outOLD      skymatrix operator+( double   rval);  // scalar addition
//outOLD      skymatrix operator-( double   rval);  // scalar subtraction
//outOLD      skymatrix operator*( double   rval);  // scalar multiplication
//outOLD
//outOLD//......    nDarray deep_copy( ); // make an image
//outOLD
//outOLD//      int operator==( skymatrix & rval);  // skymatrix comparison
//outOLD//                                             // returns 1 if they are same
//outOLD//                                             // returns 0 if they are not
//outOLD
      double & val(int row,  int col); // element selection;
      double cval(int row,  int col) const; // element selection;
//outOLD// can be used to read or write an element.
      double mmin( ); // find minimum element in the skymatrix
      double mmax( ); // find maximum element in the skymatrix
//outOLD      double mean( );  // average all the elements of the skymatrix
      void lower_print(char *msg = ""); // print lower part of
                                        // skymatrix with a message
      void upper_print(char *msg = ""); // print upper part of
                                        // skymatrix with a message
      void full_print(char *msg =  ""); // print sky matrix
                                        // as a full matrix with a message



//tempout      skymatrix  & AssembleBricksInSkyMatrix(stiffness_matrix   & Ke,  
//tempout                                             Brick3D            * b3d,
//tempout                                             Node               * node
//tempout                                            );

      skymatrix & v_ldl_factorize( ); // ldl factorizing sky matrix
      double * d_reduce_r_h_s_l_v ( double * );
      double * d_back_substitute ( double * );


    private:
      void error(char *msg1, char *msg2 = "") const; // private function
//outOLD// this one is the same as mval except that it is more convenient
//outOLD// to overload operator (row,col).
//outOLD      double & operator( )(int , int ) const;
      double & mval(int , int ) const;
// full_val inline function allows you to treat skymatrix as if it is full
// float matrix. The function will calculate position inside sky matrix
// and return appropriate number if row and col are below skyline or
// return zero (0) if row and col are above sky line
      double full_val (int , int ) const;

  private:
//outOLD    double* data(void) const;
//outOLD    void set_data_pointer(double* );
//outOLD//    int rank(void) const;
//outOLD//    void rank(int );
//outOLD    long int total_number(void) const ;
//outOLD    void total_number(long int );
//outOLD//    int* dim(void) const ;
//outOLD//    int& get_dim_pointer(void) const ;
//outOLD//    void set_dim_pointer(int* );
//outOLD//    int dim(int which) const;
//outOLD    int reference_count(int );
//outOLD    void set_reference_count(int );

  };


#endif
