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
//# CLASS:             BJtensor                                                    #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.10, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              May 28. - July  __  '93                                   #
//# UPDATE HISTORY:    july 8. '93. BJtensor02 - BJtensor multiplication             #
//#                                 inner and outer products                     #
//#                    august 17-19 '93 fixed default constructor that wasted    #
//#                                     memory ###                               #
//#                    october 11 '93 added transpose0110, transpose0101,        #
//#                                   transpose0111 so the creation of           #
//#                                   isotropic BJtensor is much easer and         #
//#                                   understandable!                            #
//#                    januar 06 '93  added BJtensor2BJmatrix_1, BJtensor2BJmatrix_2     #
//#                                   BJtensor2BJmatrix_3, inverse_1, inverse_2,     #
//#                                   inverse_3                                  #
//#                    januar 20 '93  added inverse  TRUE ONE                    #
//#                    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                    Oktobar 15 '96  added  trasposeoverbar trasposeunderbar   #
//#                                    ( the same as transpose0110 and           #
//#                                    transpose0101) but just to make it        #
//#                                    easier to work with Runesson's notation   #
//#                    28June2004     added val4 for efficiency still            #
//#                                   to be worked on                            #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

#ifndef TENSOR_HH
#define TENSOR_HH


#include "nDarray.h"
#include "BJmatrix.h"


class BJtensor : public nDarray
{
  friend class stresstensor;
  friend class straintensor;
//  private:
  public:
    char *indices1;   //  character array of indices
    char *indices2;    //  they should be checked for matching and
                       //  then operator will be identified
  public: // just send appropriate arguments to the base constructor
    BJtensor(int rank_of_BJtensor=1, double initval=0); // default constructor
    BJtensor(int rank_of_BJtensor, const int *pdim, double *values);
    BJtensor(int rank_of_BJtensor, const int *pdim, double initvalue);
    BJtensor(char *flag, int rank_of_BJtensor, const int *pdim);  // create a unit nDarray
    BJtensor(char *flag);   //this one used to send "NO" message
    BJtensor(const BJtensor & x);  // instead of: "BJtensor(nDarray & x):"   :-)
    BJtensor(const nDarray & x);

    //~BJtensor( );

    BJtensor&  operator=( const BJtensor & rval);       // BJtensor assignment

    friend BJtensor operator+(const BJtensor & , const BJtensor & ); // BJtensor addition
    friend BJtensor operator-(const BJtensor & , const BJtensor & ); // BJtensor substraction

    BJtensor operator*( const double   rval) const; //Guanzhou   // scalar multiplication
    BJtensor operator*( BJtensor & rval);       // inner/outter product
    BJtensor operator/( BJtensor & rval);       // BJtensor division rval MUST BE BJtensor of
                                            // order 0 ( i.e. scalar in BJtensor form )

    BJtensor transpose0110( ) const;      // transpose ijkl  -->> ikjl
    BJtensor transposeoverbar( ) const;   // transpose ijkl  -->> ikjl
    BJtensor transpose0101( ) const;      // transpose ijkl  -->> ilkj
    BJtensor transpose0111( ) const;      // transpose ijkl  -->> iljk
    BJtensor transposeunderbar( ) const;  // transpose ijkl  -->> iljk
    BJtensor transpose1100( ) const;      // transpose ijkl  -->> jikl    First minor symm
    BJtensor transpose0011( ) const;      // transpose ijkl  -->> ijlk    Second minor symm
    BJtensor transpose1001( ) const;      // transpose ijkl  -->> ljki
    BJtensor transpose11( ) const;        // transpose ij  -->> ji

    BJtensor symmetrize11( ) const;   // symmetrize with respect to ij
    double determinant( ) const;  // determinant of a BJtensor
    BJmatrix BJtensor2BJmatrix_1( ) const;// convert BJtensor of even order to BJmatrix
                              // to be used in inversion process
                              // I_ijkl scheme
    BJmatrix BJtensor2BJmatrix_2( ) const;// convert BJtensor of even order to BJmatrix
                              // to be used in inversion process
                              // I_ikjl scheme
    BJmatrix BJtensor2BJmatrix_3( ) const;// convert BJtensor of even order to BJmatrix
                              // to be used in inversion process
                               // I_iljk scheme
    BJtensor inverse( ) const;  // inverse of a BJtensor
                                 // I_ikjl scheme TRUE ONE         ____
                                                         //           |
    BJtensor inverse_2( )     const;    // inverse of a BJtensor   <--|
                               // I_ikjl scheme

    BJtensor & operator()(char *indices_from_user);

    void null_indices( );
//K    char *BJtensor::f_indices1( ) const;
    char *f_indices1( ) const;
//K    char *BJtensor::f_indices2( ) const;
    char *f_indices2( ) const;

  private:
//K    int BJtensor::contracted_ind(char *, char *, int *, int *, int , int );
//K    int BJtensor::uncontracted_ind(int *, int *, int);
    int contracted_ind(char *, char *, int *, int *, int , int );
    int uncontracted_ind(int *, int *, int);



};
// GLOBAL
BJtensor operator*( const double lval, const BJtensor & rval);  // REVIEWER global *

#endif




