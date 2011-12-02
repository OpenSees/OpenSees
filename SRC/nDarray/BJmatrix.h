// $Revision: 1.1 $                                                              
// $Date: 2001-08-23 16:45:50 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/nDarray/BJmatrix.h,v $                                                                

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
///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:                                                                     #
//# CLASS:             BJmatrix                                                    #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.10, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              November '92                                              #
//# UPDATE HISTORY:    05 - 07 avgust '93.  redefined as derived class from      #
//#                                 nDarray class                                #
//#                    january 06 '93  added BJmatrix2tensor_1, BJmatrix2tensor_2    #
//#                    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/
//

// MATRIX.hhH: fully functional BJmatrix class based on
// design in chapter 9.   ( Bruce Eckel: " Using C++ " )
// improved a lot by Boris Jeremic
#ifndef MATRIX_HH
#define MATRIX_HH

#include "nDarray.h"
#include "BJtensor.h"

//class vector;


class BJmatrix : public nDarray
  {
    friend class BJvector; // explanation why this one should be a friend
                         // instead of inheriting all data through protected
                         // construct -> see in J. Coplien "Advanced C++ ..."
                         // page 96.


    private:
      void error(char *msg1, char *msg2 = ""); // private function

    public:
      BJmatrix(int mrows = 1, int columns = 1, double initval = 0.0);
      BJmatrix(int mrows, int columns, double *initvalues);
 //special for vector
      BJmatrix(int rank, int mrows, int columns, double *initvalues);
      BJmatrix(int rank, int mrows, int columns, double initvalues);

      BJmatrix(char *flag, int dimension ); // create an ident BJmatrix
      BJmatrix(char *matfile); // read from a "standard" BJmatrix file
      BJmatrix(char *matfile, char *outfile); // read from a flat BJmatrix file
                                            // and write test output to another mat file
      BJmatrix(const BJmatrix & x);  // copy-initializer
      BJmatrix(const nDarray & x);

//--      ~BJmatrix( );

      int rows( void ) const;  // rows in BJmatrix
      int cols( void ) const;  // cols in BJmatrix

      BJmatrix & operator=(const BJmatrix & rval); // BJmatrix assignment

// Write a "standard" BJmatrix file:
      void write_standard(char *filename, char *msg = "");


      BJmatrix operator*( BJmatrix &); // BJmatrix multiplication
      BJmatrix operator*( double rval); // scalar multiplication
//....      vector operator*( vector &); // vector multiplication
//      vector operator*( double ); // vector multiplication

//-----//forwarding definition
//-----      virtual vector operator*( double ); // scalar multiplication 

// this is COUPLING between ordinary BJmatrix and SKYMATRIX
// it will be usefull in multiplying vectors with skyBJmatrix
//####      BJmatrix operator*(const skyBJmatrix & rval);  // BJmatrix/skyBJmatrix multiplication
//####      double & val(int row,  int col); // element selection;

//####// this one is the same as mval except that it is more convinient
//####// to overload operator (row,col).
//####// THE ROW AND COL COUNTER STARTS FROM 1 ( NOT FROM 0 )
//####      double & operator( )(int row, int col);
// can be used to read or write an element.

      BJmatrix transpose( );       // transpose a square BJmatrix
      double determinant( );
      BJmatrix inverse( );
      double mmin( );            // find minimum element in the BJmatrix
      double mmax( );            // find maximum element in the BJmatrix
      double mean( );            // average all the elements of the BJmatrix
      double sum( );             // sum all the elements in BJmatrix
      double variance( );        // statistical variance of all elements

      BJtensor BJmatrix2BJtensor_1( );     // convert BJmatrix back to 4th order BJtensor
                                   // I_ijkl scheme
      BJtensor BJmatrix2BJtensor_2( );     // convert BJmatrix back to 4th order BJtensor
                                   // I_ikjl scheme
      BJtensor BJmatrix2BJtensor_22( );  // convert BJmatrix back to 2th order BJtensor

      BJtensor BJmatrix2BJtensor_3( );     // convert BJmatrix back to 4th order BJtensor
                                   // I_iljk scheme

//####
//####// Compress columns and rows
//####      BJmatrix compress_col(int col1, int col2, int to_col);
//####      BJmatrix compress_row(int row1, int row2, int to_row);
//####

    private: // functions used by inverse() and determinant()
      void switch_columns(int col1, int co12);
      void copy_column(BJmatrix & m, int from_col, int to_col);
      BJmatrix scale( ); // scale a BJmatrix (used in L-U decomposition)
      void deepcopy(BJmatrix & from, BJmatrix & to); // make an image
      BJmatrix lu_decompose(BJmatrix & indx, int & d );
           // Returns the L-U decomposition of a BJmatrix
      void lu_back_subst(BJmatrix & indx, BJmatrix & b);
           // Uses L-U decomposition for BJmatrix inverse

      double & mval (int row, int col);  // I am still keeping mval
                                         // operator for compatibility
                                         // with old BJmatrix class members
                                         // and they start from 0 ###
          // used by BJmatrix functions which KNOW they aren't
          // exceeding the boundaries

		 public:
// Tiejun Li Jan 2000
     double *  BJmatrixtoarray(int &);
 

// // prebacen u nDarray 14 oktobra 1996
//     public:
//       BJvector eigenvalues(void);
//       BJmatrix eigenvectors(void);

// // from Numerical recipes in C
//     private:
//       void tqli(double * d, double * e, int n, double ** z);
//       void tred2(double ** a, int n, double * d, double * e);
//       void eigsrt(double * d, double ** v, int n);
          
};
#endif 
