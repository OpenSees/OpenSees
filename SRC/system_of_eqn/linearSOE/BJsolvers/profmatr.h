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
//##############################################################################
//# COPYRIGHT (C):     :-))                                                    #
//# PROJECT:           Object Oriented Finite Element Program                  #
//# PURPOSE:                                                                   #
//# CLASS:             profilematrix                                               #
//#                                                                            #
//# VERSION:                                                                   #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.10, SUN C++ ver=2.1 )#
//# TARGET OS:         DOS || UNIX || . . .                                    #
//# PROGRAMMER(S):     Boris Jeremic                                           #
//#                                                                            #
//#                                                                            #
//# DATE:              September 13 '94 Unsymmetric Solver -> profile solver   #
//#                                     ( see FEM4ed by O.Z. and R.T.)         #
//#                                                                            #
//# UPDATE HISTORY:    September 27 '94 profile solver for symmetric and       #
//#                                     Nonsymmetric systems works!            #
//#                                     (from FEM4ed by O.Z. and R.T.)         #
//#                                                                            #
//#                                                                            #
//##############################################################################
// Profile Sparse Matrix Class
// this one is based on Zienkiewicz's and Taylor's book!

#ifndef PROFMATR_HH
#define PROFMATR_HH

#include <basics.h>

#include <BJvector.h>

//tempout#include <femdata.h>
//tempout#include <brick3d.h>
//tempout#include <node.h>
//tempout#include <stifmat.h>



class profilematrix_rep
  {
    public:
      friend class profilematrix;
    private:
      double *al;   // lower triangular part of matrix
      double *au;   // upper triangular part of matrix
      double *ad;   // diagonals of triangular of matrix
      //double max_element;   // max ele. of  matrix
      //double min_element;   // min ele. of  matrix

      int *jp;      // pointers to bottom of columns of al and au arrays
      
      int * columnheight;    // pointers to the height of each column may be replaced by jp __Zhaohui

      int neq;      // dimension of a system ( as square matrix )

      int total_numb; // total number of elements in al and/or au

      char flag;    // if flag=='S' symmetric; if flag=='N' NON(un)symmetric

      int n;        // reference count

    public:
// overloading operator new and delete in profilematrix_rep class  ########
      void * operator new(size_t s); // see C++ reference manual by
      void operator delete(void *);  // by ELLIS and STROUSTRUP page 283.
                                   // and ECKEL page 529.
  };


class profilematrix
  {
    private:
      profilematrix_rep * pc_profilematrix_rep;
    public:
// I don't need default constructor, this is not a very common data type!!!!!!!!!!!
//      profilematrix(int matrix_order=1, double init_value=0.0 ); // default constructor
// this constructor just for testing purposes!
      
//tempout      profilematrix(FEModelData & FEMD, 
//tempout                         Brick3D * b3d, 
//tempout			 Node *  node);

      profilematrix(int matrix_order,
                    int *jp_init,
                    char flag_init,
                    double *au_init_val,
                    double *al_init_val,
                    double *ad_init_val);
// The real constructor. This one will be initialized to the original size
// with the init_val values and then in stiffness matrix class those
// default inti_val values will be altered ( to the right stiffness ).
      profilematrix(int matrix_order,
                    int *jp_init,
                    char flag_init,
                    double au_init_val,
                    double al_init_val,
                    double ad_init_val);
//      skymatrix(char *flag, int dimension ); // create an ident skymatrix
      profilematrix(const profilematrix & x); // copy-initializer
      ~profilematrix();


      int dimension_of_profile_M(void ) const; // dimension of profile matrix
      int *get_jp(void) const; // get pointer to array of
                              // Locations of Diagonals
      profilematrix  &operator=(const profilematrix & rval); // profilematrix assignment

//....// This is from JOOP May/June 1990 after ARKoenig
      profilematrix& operator +=( const profilematrix & ); // profilematrix addition
      friend profilematrix operator+(const profilematrix & , const profilematrix & ); // profilematrix addition
//....// This is from JOOP May/June 1990 after ARKoenig
      profilematrix& operator -=( const profilematrix & ); // profilematrix subtraction
      friend profilematrix operator-(const profilematrix & , const profilematrix & ); // profilematrix subtraction

      profilematrix operator+( double   rval);  // scalar addition
      profilematrix operator-( double   rval);  // scalar subtraction
      profilematrix operator*( double   rval);  // scalar multiplication
      double  *     operator*( BJvector  &arg);  // profile multiplied by vector__Zhaohui 

//......    nDarray deep_copy( ); // make an image

//      int operator==( profilematrix & rval);  // profilematrix comparison
//                                             // returns 1 if they are same
//                                             // returns 0 if they are not

      double & val(int row,  int col); // element selection;
      double cval(int row,  int col) const; // element selection;
// can be used to read or write an element.
      double mmin( ); // find minimum element in the skymatrix
      double mmax( ); // find maximum element in the skymatrix
      double mean( );  // average all the elements of the skymatrix
      double get_max(void);        // get max element
      //double get_min(void);        // get max element
      //void  set_mm(void);        // set max and min element

     //Zhaohui  added to set the max-element in diagonal of eqn_no_shake!
     void  set_penalty_element(int *eqn_no_shake, int no_of_shake, double max_element);

      void lower_print(char *msg = ""); // print lower part of
                                        // skymatrix with a message
      void upper_print(char *msg = ""); // print upper part of
                                        // skymatrix with a message

//      void profile_init( FEModelData      & FEMD);	        // initialize profile
// This sub has been merged with initializor No.1                    //Zhaohui added   

      void profile_jp_print( );         // print jp vector
                                        //Zhaohui added  

      void profile_ad_print( void );
      void profile_al_print( );         // print al vector
                                        //Zhaohui added  

      void full_print(char *msg =  ""); // print sky matrix
                                        // as a full matrix with a message

//tempout      profilematrix & AssembleBricksInProfMatrix(stiffness_matrix & Ke,
//tempout                                            Brick3D          * b3d,
//tempout                                            Node             * node,
//tempout			                                         FEModelData      & FEMD,
//tempout                                            float              scale // used to add damping matrix C to Kbar Zhaohui 
//tempout	                                   );
//tempout

    private:
      void error(char *msg1, char *msg2 = "") const; // private function
// this one is the same as mval except that it is more convenient
// to overload operator (row,col).
      double & operator( )(int , int ) const;
      double & mval(int , int ) const;
// full_val inline function allows you to treat profilematrix as if it is full
// float matrix. The function will calculate position inside sky matrix
// and return appropriate number if row and col are below skyline or
// return zero (0) if row and col are above sky line
      double full_val (int , int ) const;

  private:
//..    double* data(void) const;
//..    void set_data_pointer(double* );
//    int rank(void) const;
//    void rank(int );
    long int total_number(void) const ;
    void total_number(int );
//    int* dim(void) const ;
//    int& get_dim_pointer(void) const ;
//    void set_dim_pointer(int* );
//    int dim(int which) const;
    int reference_count(int );
    void set_reference_count(int );

// profile solvers  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  public:
//    void datri(); // triangular decomposition of profile symmetric/unsymmetric matrix!
    profilematrix & datri(void); // triangular decomposition of profile symmetric/unsymmetric matrix!
    double * dasol(double *); // backsubstitution of profile symmetric/unsymmetric matrix!
    double datest(double* , int );
    double dot(double* , double* , int );
    void dredu(double* , double* , double* , int , char , double* );
    void saxpb(double* , double* , double , int , double* );



  };

#endif 
