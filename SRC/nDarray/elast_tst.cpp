///*
//################################################################################
//# COPY-YES  (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:                                                                     #
//# CLASS:                                                                       #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0                                            #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              April 2001                                                #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/
#include "time.h"
#include "limits.h"
#include "math.h"
#include "float.h"
// nDarray tools
#include "BJtensor.h"
#include "stresst.h"
#include "straint.h"
#include "BJvector.h"
#include "BJmatrix.h"
#include "iostream.h"

int main(int argc, char *argv[])
{
  if (argc != 2)
    {
      puts("\a\n usage: elast_tst test01.std\n");
      exit( 1 );
    }
::printf("\n\n-------------  MACHINE (CPU) DEPENDENT THINGS  --------------\n\n");

// defining machine epsilon for different built in data types supported by C++
float        float_eps       = f_macheps();  // or use float.h values FLT_EPSILON
double       double_eps      = d_macheps();  // or use float.h values DBL_EPSILON
//long double  long_double_eps = ld_macheps(); // or use float.h values LDBL_EPSILON

::printf("\n float macheps = %.20e \n",float_eps);
::printf(" double macheps = %.20e \n",double_eps);
//::printf(" long double macheps = %.20Le \n\n\n",long_double_eps);


// Usual tolerance is defined as square root of machine epsilon
float        sqrt_f_eps       = sqrt(f_macheps());
double       sqrt_d_eps       = sqrt(d_macheps());
//long double  sqrt_ld_eps      = sqrt(ld_macheps());

::printf("\n sqrt float macheps = %.20e \n",sqrt_f_eps);
::printf(" sqrt double macheps = %.20e \n",sqrt_d_eps);
//::printf(" sqrt long double macheps = %.20Le \n",sqrt_ld_eps);


double one = 1.0;
double onepmc = 1.0+double_eps;
::printf("\n one = %.32e \n",one);
::printf("onepmc = %.32e \n",onepmc);
if (one == onepmc) 
  {
    printf("1.0 == 1.0+double_eps -->>  OK \n");
  }
else 
  {
    printf("1.0 == 1.0+double_eps -->>  NOT OK \n");
  }   



double one2 = 10e20;
double onepmc2 = 10e20*(1.0+double_eps);
::printf("\n one2 = %.32e \n",one2);
::printf("onepmc2 = %.32e \n",onepmc2);
if (one2 == onepmc2) 
  {
    printf("10e20 == 10e20*(1.0+double_eps) -->>  OK \n");
  }
else 
  {
    printf("10e20 == 10e20*(1.0+double_eps) -->>  NOT OK \n");
  }   



double one3 = 1.0;
double onepmc3 = 1.0+double_eps*10000.0;
::printf("\n one3 = %.32e \n",one3);
::printf("onepmc3 = %.32e\n",onepmc3);
if (one3 == onepmc3) 
  {
    printf("1.0 == 1.0+double_eps*10.0 -->>  OK \n");
  }
else 
  {
    printf("1.0 == 1.0+double_eps*10.0 -->>  NOT OK \n");
  }   


::printf("\n\n----------------  MATRIX  -----------------\n\n");

// New data type MATRIX is defined in such a way that it can be regarded
// as concrete data type.
// default constructor:
  BJmatrix matrix1;
  matrix1.print("m1","matrix1");

  static double initvalue[] = {1.0, 2.0,
                               3.0, 4.0};
// constructor with matrix initialization
  BJmatrix matrix2(2,2,initvalue);
  matrix2.print("m2","matrix2");

// unit matrix initialization ( third order )
  BJmatrix matrix3("I",3);
  matrix3.print("m3","matrix3");

// matrix initialization from the file that is in predefined format!
  BJmatrix  m(argv[1]);
  m.print("m","m is");
//====
// matrix assignment
  BJmatrix first_matrix = m;
  first_matrix.print("fm","first_matrix is");


// ---------------------------------
::printf("\n\n--------EIGENVECTORS and EIGENVALUES-----------\n");
  static double symm_initvalue[] = {1.0, 0.0,
                                    0.0, 0.0};
// constructor with matrix initialization
  BJmatrix sym_matrix(2,2,symm_initvalue);
  sym_matrix.print("m","Sym_matrix");

// Eigenvalues
  BJvector eval = sym_matrix.eigenvalues();
  eval.print("ev","Eigenvalues of sym_matrix");
// Eigenvectors
  BJmatrix evect = sym_matrix.eigenvectors();
  evect.print("ev","Eigenvectors of sym_matrix");

// matrix to file "test.dat"
  m.write_standard("test.dat", "this is a test");

// equivalence of two matrices ( within tolerance sqrt(macheps) )
  if ( first_matrix == m ) ::printf("\n first_matrix == m  TRUE \n");
  else ::printf("\a\n first_matrix == m  NOTTRUE \n");

// matrix addition
  BJmatrix m2 =m+m;
  m2.print("m2","\nm+m = \n");

// equivalence of two matrices ( within tolerance sqrt(macheps) )
  if ( m2 == m ) ::printf("\a\n m2 == m  TRUE \n");
  else ::printf("\n m2 == m  NOTTRUE \n");
// matrix substraction
  BJmatrix m3 =m-m;
  m3.print("m3","\nm-m = \n");

// matrix multiplication
  BJmatrix m4 = m*m;
  m4.print("m4","\n m*m = \n");

// matrix transposition
  BJmatrix mtran = m.transpose();
  mtran.print("mt","transpose of m is");

// determinant of matrix

  printf("determinant = %6.6f \n",m.determinant());

//  exit(0);
// matrix inversion
  BJmatrix minv = m.inverse();
  minv.print("mi","inverse of m is");
// check matrix inversion
  ( m * minv ).print("1"," m * m.inverse() ");

// minima, maxima, mean and variance values for matrix:
  printf("min = %6.6f \n", m.mmin() );
  printf("max = %6.6f \n", m.mmax() );
  printf("mean = %6.6f \n", m.mean() );
  printf("variance = %6.6f \n", m.variance() );

//====//===printf("\n\n-----------------  SKY MATRIX ----------------------\n\n");
//====//===
//====//===// skymatrix as defined by K.J. Bathe in "FE procedures in Engineering Analysis"
//====//===// The example is from that book too!
//====//===
//====//===int maxa[] = { 1,   2,         4,          6,          8,                  13};
//====//===double a[] = {2.0, 3.0, -2.0, 5.0, -2.0, 10.0, -3.0, 10.0, 4.0, 0.0, 0.0, -1.0};
//====//===double rhs[] = { 0.0, 1.0, 0.0, 0.0, 0.0 };
//====//===
//====//===// skymatrix constructor
//====//===  skymatrix first( 5, maxa, a);
//====//===
//====//===// simple print functions:
//====//===  first.full_print("\nfirst Sparse Skyline Matrix as FULL :\n ");
//====//===  first.lower_print("\nfirst Sparse Skyline Matrix :\n ");
//====//===  first.upper_print("\nfirst Sparse Skyline Matrix :\n ");
//====//===
//====//===// skymatrix assignment
//====//===  skymatrix second = first;
//====//===  second.full_print("\nsecond Sparse Skyline Matrix as FULL :\n ");
//====//===
//====//===//  skymatrix third; // not defined !
//====//===// assignments:
//====//===  skymatrix third = second = first;
//====//===
//====//===// LDL factorization ( see K.J. Bathe's book! )
//====//===  first = first.v_ldl_factorize();
//====//===
//====//===// simple print functions:
//====//===  first.lower_print("\n lower first but factorized :\n ");
//====//===  first.full_print("\nfirst but factorized Sparse Skyline Matrix as FULL :\n ");
//====//===  first.upper_print("\n upper first but factorized :\n ");
//====//===
//====//===// right hand side reduction ( see K.J. Bathe's book! ) ( vetor of unknowns is rhs )
//====//===  first.d_reduce_r_h_s_l_v( rhs );
//====//===  int i=0;
//====//===  for( i=0 ; i<5 ; i++ )
//====//===    {
//====//===      printf("intermediate_rhs[%d] = %8.4f\n",i, rhs[i]);
//====//===    }
//====//===
//====//===// backsubstitution
//====//===  first.d_back_substitute ( rhs );
//====//===  for( i=0 ; i<5 ; i++ )
//====//===    {
//====//===      printf("solution_rhs[%d] = %8.4f\n",i, rhs[i]);
//====//===    }
//====//===
//====//===// minima, maxima and mean values for skymatrix:
//====//===  printf("min = %6.6f \n", first.mmin() );
//====//===  printf("max = %6.6f \n", first.mmax() );
//====//===  printf("mean = %6.6f \n", first.mean() );
//====//===//
//====//===//
//====//===
//====//===printf("\n\n-----------------  PROFILE MATRIX ----------------------\n");
//====//===printf("-----------------  EXAMPLE 1 ----------------------\n\n");
//====//===
//====//===// profilematrix as defined by O.C. Zienkiewicz and R.L. Taylor
//====//===// in "The FEM  - fourth edition"
//====//===// The example is from: K.J. Bathe: "FE procedures in Engineering Analysis"
//====//===
//====//===int jp[]    = { 0,    1,    2,    3,                   7};
//====//===double au[] = {    -2.0, -2.0, -3.0, -1.0, 0.0, 0.0, 4.0};
//====//===double ad[] = { 2.0, 3.0, 5.0, 10.0, 10.0 };
//====//===double* al  = au;
//====//===double b[] = { 0.0, 1.0, 0.0, 0.0, 0.0 };
//====//===
//====//===// profile matrix constructor
//====//===  profilematrix first_profile( 5, jp, 'S', au, al, ad);
//====//===
//====//===// simple print function
//====//===  first_profile.full_print("\nfirst_profile Sparse Profile Matrix as FULL :\n ");
//====//===
//====//===// asignment operetor
//====//===  profilematrix second_profile = first_profile;
//====//===  second_profile.full_print("\nsecond_profile Sparse Profile Matrix as FULL :\n ");
//====//===
//====//===// triangularization ( Crout's method ); see: O.C. Zienkiewicz and R.L. Taylor:
//====//===// "The FEM - fourth edition"
//====//===  first_profile = first_profile.datri();
//====//===  first_profile.full_print("\nfirst but factorized Sparse Profile Matrix as FULL :\n ");
//====//===  ::printf("\n\n");
//====//===
//====//===// right hand side reduction and backsubstitution ( vector of unknowns is b )
//====//===  first_profile.dasol(b);
//====//===  for( i=0 ; i<5 ; i++ )
//====//===    {
//====//===      printf("solution_b[%d] = %8.4f\n",i, b[i]);
//====//===    }
//====//===
//====//===  printf("\n\nmin = %6.6f \n", first_profile.mmin() );
//====//===  printf("max = %6.6f \n", first_profile.mmax() );
//====//===  printf("mean = %6.6f \n", first_profile.mean() );
//====//===
//====//===
//====//===printf("\n-----------------  EXAMPLE 2 ----------------------\n\n");
//====//===
//====//===// The example is from: O.C. Zienkiewicz and R.L. Taylor:
//====//===// "The FEM - fourth edition"
//====//===int jp2[]    = { 0,    1,    3};
//====//===double au2[] = {    2.0, 1.0, 2.0};
//====//===double ad2[] = { 4.0, 4.0, 4.0};
//====//===double* al2  = au2;
//====//===double b2[] = { 0.0, 0.0, 1.0 };
//====//===
//====//===// profile matrix constructor
//====//===  profilematrix first_profile2( 3, jp2, 'S', au2, al2, ad2);
//====//===
//====//===// simple print function
//====//===  first_profile2.full_print("\nfirst_profile Sparse Profile Matrix as FULL :\n ");
//====//===
//====//===// asignment operetor
//====//===  profilematrix second_profile2 = first_profile2;
//====//===  second_profile2.full_print("\nsecond_profile Sparse Profile Matrix as FULL :\n ");
//====//===
//====//===// triangularization ( Crout's method ); see: O.C. Zienkiewicz and R.L. Taylor:
//====//===// "The FEM - fourth edition"
//====//===  profilematrix first_profile2TRI = first_profile2.datri();
//====//===  first_profile2TRI.full_print("\nfirst but factorized Sparse Profile Matrix as FULL :\n ");
//====//===  ::printf("\n\n");
//====//===
//====//===// right hand side reduction and backsubstitution ( vector of unknowns is b2 )
//====//===  first_profile2TRI.dasol( b2 );
//====//===  for( i=0 ; i<3 ; i++ )
//====//===    {
//====//===      printf("solution_b[%d] = %8.4f\n",i, b2[i]);
//====//===    }
//====//===
//====//===  printf("\n\nmin = %6.6f \n", first_profile2.mmin() );
//====//===  printf("max = %6.6f \n", first_profile2.mmax() );
//====//===  printf("mean = %6.6f \n", first_profile2.mean() );
//====//===
//====//===
//====//===printf("\n-----------------  EXAMPLE 3 ----------------------\n\n");
//====//===
//====//===// The main advantage of profile solver as described in
//====//===// O.C. Zienkiewicz and R.L. Taylor: "The FEM - fourth edition"
//====//===// is that it supports nonsymetric matrices stored in profile format.
//====//===
//====//===int jp3[]    = { 0,    1,    3};
//====//===double au3[] = { 2.0, 3.0, 2.0};
//====//===double ad3[] = { 4.0, 4.0, 4.0};
//====//===double al3[] = { 2.0, 1.0, 2.0};
//====//===double  b3[] = { 1.0, 1.0, 1.0 };
//====//===
//====//===// nonsymetric profile matrix constructor
//====//===  profilematrix first_profile3( 3, jp3, 'N', au3, al3, ad3);
//====//===
//====//===// print out
//====//===  first_profile3.full_print("\nfirst_profile Sparse Profile Matrix as FULL :\n ");
//====//===
//====//===// triangularization ( Crout's method ); see: O.C. Zienkiewicz and R.L. Taylor:
//====//===// "The FEM - fourth edition"
//====//===  profilematrix first_profile3TRI = first_profile3.datri();
//====//===  first_profile3TRI.full_print("\nfirst but factorized \
//====//===  Sparse Profile Matrix as FULL :\n ");
//====//===  ::printf("\n\n");
//====//===
//====//===// right hand side reduction and backsubstitution ( vector of unknowns is b3 )
//====//===  first_profile3TRI.dasol( b3 );
//====//===  for( i=0 ; i<3 ; i++ )
//====//===    {
//====//===      printf("solution_b[%d] = %8.4f\n",i, b3
//====//===      [i]);
//====//===    }
//====//===
//====//===  printf("\n\nmin = %6.6f \n", first_profile3.mmin() );
//====//===  printf("max = %6.6f \n", first_profile3.mmax() );
//====//===  printf("mean = %6.6f \n", first_profile3.mean() );
//====//===
//====//===
//====printf("\n\n------------  VECTOR  -------------------\n\n");

static double vect1[] = {2.0, 3.0, -2.0, 5.0, -2.0};
static double vect2[] = {1.0, 1.0, 1.0, 1.0, 1.0};

// vector constructor
  BJvector first_vector( 5, vect1);
  first_vector.print("f","\n\nfirst_vector vector\n");

// vector constructor
  BJvector second_vector ( 5, vect2);
  second_vector.print("s","\n\nsecond_vector vector\n");

// default constructor for vector
  BJvector third_vector;
// vector assignment
  third_vector = first_vector + second_vector;
  third_vector.print("t","\nfirst_vector + second_vector = \n");

  BJvector dot_product = first_vector.transpose() * second_vector; 
//  vector dot product;
//  dot_product = first_vector.transpose() * second_vector;
//====  dot_product.print("d","\nfirst_vector * second_vector \n");
//====
      
printf("\n\n------------VECTOR AND MATRIX -----------------\n\n");
static double vect3[] = {1.0, 2.0, 3.0, 4.0, 5.0};

//default construcotrs
BJmatrix res1;
BJmatrix res2;

// matrix constructor
BJmatrix vect_as_matrix(5, 1, vect3);
BJvector vect(5, vect3);

first_matrix.print("f","\n first matrix");

// vector and matrix are quite similar !
vect_as_matrix.print("vm","\n vector as matrix");
vect.print("v","\n vector");

// this should give the same answer! ( vector multiplication with matrix
res1 = vect_as_matrix.transpose() * first_matrix;
res2 = vect.transpose() * first_matrix;

res1.print("r1","\n result");
res2.print("r2","\n result");



::printf("\n\n----------------- TENSOR --------------------\n\n");
  
//  Tensor t0();
//  Tensor *ta = new Tensor;
  
// tensor default construcotr
  tensor t1;
  t1.print("t1","\n tensor t1");

  static double t2values[] = {  1,2,3,
                                4,5,6,
                                7,8,9  };
// tensor constructor with values assignments ( order 2 ; like matrix )
  tensor t2( 2, def_dim_2, t2values);
  t2.print("t2","\ntensor t2 (2nd-order with values assignement)");

  cerr << "rank of t2: " << t2.rank() << endl;
  for (int ii=1; ii<=t2.rank(); ii++) {
    cerr << "dimension of t2 in " << ii << " is " << t2.dim(ii) << endl;
  }


  tensor tst( 2, def_dim_2, t2values); 
//  t2 = t2*2.0;
//  t2.print("t2x2","\ntensor t2 x 2.0");
  t2.print("t2","\noriginal tensor t2 ");
  tst.print("tst","\ntst ");

  tensor t2xt2 = t2("ij")*t2("jk");
  t2xt2.print("t2xt2","\ntensor t2 x t2");

  tensor t2xtst = t2("ij")*tst("jk");
  t2xtst.print("t2xtst","\ntensor t2 x tst");
//  exit(0);
  t2.print("t2","\noriginal tensor t2 ");

// tensor constructor with 0.0 assignment ( order 4 ; 4D array )
  tensor ZERO(4,def_dim_4,0.0);
  ZERO.print("z","\ntensor ZERO (4-th order with value assignment)");

  static double t4val[] =
{ 11.100,  12.300,  13.100,   15.230,  16.450,   14.450,  17.450,  18.657,  19.340,
 110.020, 111.000, 112.030,  114.034, 115.043,  113.450, 116.670, 117.009, 118.060,
 128.400, 129.500, 130.060,  132.560, 133.000,  131.680, 134.100, 135.030, 136.700,
 119.003, 120.400, 121.004,  123.045, 124.023,  122.546, 125.023, 126.043, 127.000,
 137.500, 138.600, 139.000,  141.067, 142.230,  140.120, 143.250, 144.030, 145.900,
 146.060, 147.700, 148.990,  150.870, 151.000,  149.780, 152.310, 153.200, 154.700,
 164.050, 165.900, 166.003,  168.450, 169.023,  167.067, 170.500, 171.567, 172.500,
 155.070, 156.800, 157.067,  159.000, 160.230,  158.234, 161.470, 162.234, 163.800,
 173.090, 174.030, 175.340,  177.060, 178.500,  176.000, 179.070, 180.088, 181.200};

// tensor constructor with values assignment ( order 4 ; 4D array )
  tensor t4(4, def_dim_4, t4val);
  t4.print("t4","\ntensor t4 (4th-order with values assignment)");

//  Some relevant notes from the Dec. 96 Draft (as pointed out by KAI folks):
//  
//  "A binary operator shall be implemented either by a  non-static  member
//    function (_class.mfct_) with one parameter or by a non-member function
//    with two parameters.  Thus, for any binary  operator  @,  x@y  can  be
//    interpreted as either x.operator@(y) or operator@(x,y)." [13.5.2]
//  
//  "The order of evaluation of arguments [of a function call] is
//  unspecified."  [5.2.2 paragraph 9]
//  
//  
// SO THAT IS WHY I NEED THIS temp
// It is not evident which operator will be called up first: operator() or operator* .
// 
//% ////////////////////////////////////////////////////////////////////////
// testing the multi products using indicial notation:
  tensor tst1;
  tensor temp = t2("ij") *  t4("ijkl");
  tst1 = temp("kl")  * t4("klpq")  * t2("pq");
  tst1.print("tst1","\ntensor tst1 = t2(\"ij\")*t4(\"ijkl\")*t4(\"klpq\")*t2(\"pq\");");


// testing inverse function for 4. order tensor. Inverse for 4. order tensor
// is done by converting that tensor to matrix, inverting that matrix
// and finally converting matrix back to tensor
  tensor t4inv_2 = t4.inverse_2();
  t4inv_2.print("t4i","\ntensor t4 inverted_2");
  tensor unity_2 = t4("ijkl")*t4inv_2("klpq");
  unity_2.print("uno2","\ntensor unity_2");

//....................................................................
// There  are several built in tensor types. One of them is:
// Levi-Civita permutation tensor
  tensor e("e",3,def_dim_3);
  e.print("e","\nLevi-Civita permutation tensor e");

// The other one is:
// Kronecker delta tensor
  tensor I2("I", 2, def_dim_2);
  I2.print("I2","\ntensor I2 (2nd order unit tensor: Kronecker Delta -> d_ij)");



	 tensor Ipmnq = I2("pm") * I2("nq");
  tensor Imnpq = I2("mn") * I2("pq");



////  tensor Apqmn = alpha("pq") * I2("mn");
//  tensor X  = (1.0/3.0) * Imnpq ;// - Apqmn * (1.0/3.0);
//	 X.print();


// tensor multiplications 
  tensor I2I2 = I2("ij")*I2("ij");
  I2I2.print("I2I2","\ntensor I2 * I2");

// tensor multiplications ( from right)
  tensor I2PI = I2*PI;
  I2PI.print("I2PI","\ntensor I2 * PI");

//// tensor multiplications ( from left)
//  tensor PII2 = PI*I2;
//  PII2.print("PII2","\ntensor PI * I2");

//// Unary minus
//  tensor mPII2 = -PII2;
//  mPII2.print("mPII2","\ntensor -PI * I2");

// trace function
  double deltatrace = I2.trace();
  ::printf("\ntrace of Kronecker delta is %f \n",deltatrace);

// determinant of 2. order tensor only!
  double deltadet = I2.determinant();
  ::printf("\ndeterminant of Kronecker delta is %f \n",deltadet);

// comparison  of tensor ( again within sqrt(macheps) tolerance )
  tensor I2again = I2;
  if ( I2again == I2 ) ::printf("\n I2again == I2  TRUE (OK)\n");
  else ::printf("\a\n\n\n I2again == I2  NOTTRUE \n\n\n");

//....................................................................
// some of the 4. order tensor used in conituum mechanics:
// the most general representation of the fourth order isotropic tensor
// includes the following fourth orther unit isotropic tensors:
  tensor I_ijkl( 4, def_dim_4, 0.0 );
  I_ijkl = I2("ij")*I2("kl");
//  I_ijkl.print("ijkl","\ntensor I_ijkl = d_ij d_kl)");
  tensor I_ikjl( 4, def_dim_4, 0.0 );
//  I_ikjl = I2("ij")*I2("kl");
  I_ikjl = I_ijkl.transpose0110();
//  I_ikjl.print("ikjl","\ntensor I_ikjl) ");
  tensor I_ikjl_inv_2 = I_ikjl.inverse_2();
//  I_ikjl_inv_2.print("I_ikjl","\ntensor I_ikjl_inverse_2)");
  if ( I_ikjl == I_ikjl_inv_2 ) ::printf("\n\n I_ikjl == I_ikjl_inv_2 !\n");
  else ::printf("\a\n\n           I_ikjl != I_ikjl_inv_2 !\n");

  tensor I_iljk( 4, def_dim_4, 0.0 );
//  I_ikjl = I2("ij")*I2("kl");
//..  I_iljk = I_ikjl.transpose0101();
  I_iljk = I_ijkl.transpose0111();
//  I_iljk.print("iljk","\ntensor I_iljk) ");
// To check out this three fourth-order isotropic tensor please see
// W.Michael Lai, David Rubin, Erhard Krempl
// " Introduction to Continuum Mechanics"
// QA808.2
// ISBN 0-08-022699-X

//....................................................................
// Multiplications between 4. order tensors
  ::printf("\n\n intermultiplications \n");
//  tensor I_ijkl_I_ikjl = I_ijkl("ijkl")*I_ikjl("ikjl");
  tensor I_ijkl_I_ikjl = I_ijkl("ijkl")*I_ikjl("ijkl");
  I_ijkl_I_ikjl.print("i","\n\n I_ijkl*I_ikjl \n");

//  tensor I_ijkl_I_iljk = I_ijkl("ijkl")*I_iljk("iljk");
  tensor I_ijkl_I_iljk = I_ijkl("ijkl")*I_iljk("ijkl");
  I_ijkl_I_iljk.print("i","\n\n I_ijkl*I_iljk \n");

//  tensor I_ikjl_I_iljk = I_ikjl("ikjl")*I_iljk("iljk");
  tensor I_ikjl_I_iljk = I_ikjl("ijkl")*I_iljk("ijkl");
 I_ikjl_I_iljk.print("i","\n\n I_ikjl*I_iljk \n");

//....................................................................
// More multiplications between 4. order tensors
  ::printf("\n\n   intermultiplications same with same \n");
  tensor I_ijkl_I_ijkl = I_ijkl("ijkl")*I_ijkl("ijkl");
  I_ijkl_I_ijkl.print("i","\n\n I_ijkl*I_ijkl \n");

  tensor I_ikjl_I_ikjl = I_ikjl("ikjl")*I_ikjl("ikjl");
  I_ikjl_I_ikjl.print("i","\n\n I_ikjl*I_ikjl \n");

  tensor I_iljk_I_iljk = I_iljk("iljk")*I_iljk("iljk");
  I_iljk_I_iljk.print("i","\n\n I_iljk*I_iljk \n");


// tensor additions ans scalar multiplications
// symmetric part of fourth oder unit isotropic tensor
  tensor I4s = (I_ikjl+I_iljk)*0.5;
//  tensor I4s = 0.5*(I_ikjl+I_iljk);
  I4s.print("I4s","\n symmetric tensor I4s = (I_ikjl+I_iljk)*0.5 ");
// skew-symmetric part of fourth oder unit isotropic tensor
  tensor I4sk = (I_ikjl-I_iljk)*0.5;
  I4sk.print("I4sk","\n skew-symmetric tensor I4sk = (I_ikjl-I_iljk)*0.5 ");

// equvalence check ( they should be different )
  if ( I4s == I4sk ) ::printf("\n\n I4s == I4sk !\n");
  else ::printf("\a\n\n           I4s != I4sk !\n");


//....................................................................
// e-delta identity  ( see Lubliner ( QA931 .L939 1990) page 2.
// and Lai ( QA808.2 L3 1978) page 12 )
  tensor ee = e("ijm")*e("klm");
  ee.print("ee","\ntensor e_ijm*e_klm");
  tensor id = ee - (I_ikjl - I_iljk);
  if ( id == ZERO ) ::printf("\n\n e-delta identity HOLDS !! \n\n");
  tensor ee1 = e("ilm")*e("jlm");
  ee1.print("ee1","\n\n e(\"ilm\")*e(\"jlm\") \n");
  tensor ee2 = e("ijk")*e("ijk");
  ee2.print("ee2","\n\n e(\"ijk\")*e(\"ijk\") \n");









//....................................................................
// Linear Isotropic Elasticity Tensor
// building elasticity tensor
  double Ey = 20000;
  double nu = 0.2;
// building stiffness tensor in one command line
  tensor E1 = (I_ijkl*((Ey*nu*2.0)/(2.0*(1.0+nu)*(1-2.0*nu)))) +
              ( (I_ikjl+I_iljk)*(Ey/(2.0*(1.0+nu))));
// building compliance tensor in one command line
  tensor D1= (I_ijkl * (-nu/Ey)) + ( (I_ikjl+I_iljk) * ((1.0+nu)/(2.0*Ey)));

// multiplication between them
  tensor test = E1("ijkl")*D1("klpq");
  test.print("t","\n\n\n test =  E1(\"ijkl\")*D1(\"klpq\") \n");
// and this should be equal to the symmetric part of 4. order unit tensor
  if ( test == I4s ) ::printf("\n test == I4s  TRUE  ( up to sqrt(macheps())) \n");
  else ::printf("\a\n\n\n   test == I4s    NOTTRUE  ( up to sqrt(macheps())) \n\n\n");

//............Different Way...........................................
// Linear Isotropic Elasticity Tensor
  double lambda = nu * Ey / (1. + nu) / (1. - 2. * nu);
  double mu = Ey / (2. * (1. + nu));
  ::printf("\n  lambda = %.12e\n",lambda);
  ::printf("      mu = %.12e\n",mu);

  ::printf("\nYoung Modulus = %.4e",Ey);
  ::printf("\nPoisson ratio = %.4e\n",nu);
  ::printf(  "    lambda + 2 mu = %.4e\n", (lambda + 2 * mu));

// building elasticity tensor
  tensor E2 = (I_ijkl * lambda) + (I4s * (2. * mu));
  E2.print("E2","tensor E2: elastic moduli Tensor = lambda*I2*I2+2*mu*I4s");

  tensor E2inv = E2.inverse();

// building compliance tensor
  tensor D2= (I_ijkl * (-nu/Ey)) + (I4s * (1./(2.*mu)));
  D2.print("D2","tensor D2: compliance Tensor (I_ijkl * (-nu/Ey)) + (I4s * (1./2./mu))");

  if ( E2inv == D2 ) ::printf("\n E2inv == D2   TRUE \n");
  else ::printf("\a\n\n\n  E2inv == D2   NOTTRUE \n\n\n");

  if ( E1 == E2 ) ::printf("\n E1 == E2  TRUE \n");
  else ::printf("\a\n\n\n  E1 == E2  NOTTRUE \n\n\n");

  if ( D1 == D2 ) ::printf("\n D1 == D2  TRUE \n");
  else ::printf("\a\n\n\n  D1 == D2  NOTTRUE \n\n\n");


 
//::printf("stack length is : %u \n",_stklen);
double PIo3 = PI/3.0;
//::printf("Pi/3.0 = %f \n",PIo3);
//-------------------------------------------
  double p_start     = 5000.000;
  double q_start     = 2000.0;
  double theta_start = PIo3;
//..//  double theta_start = 0.4;

stresstensor stress = stress.pqtheta2stress(p_start, q_start, theta_start);
stress.report("stress\n");


straintensor strain = D2("ijkl") * stress("kl");
strain.report("strain = D2(\"ijkl\") * stress(\"kl\") \n");






::printf("\n\n\n\n\n");
::printf("\n\n----------------- STRESS TENSOR --------(Kaspar Willam tests)--------\n\n");
  
	 static double stressvalues[] = {  1.0, 2.0, 3.0,
                                    2.0, 2.0, 5.0,
                                    3.0, 5.0, 9.0  };

  stresstensor sigma(stressvalues);
  sigma.report("\n\n stress tensor sigma (2nd-order with values assignment)");


	 stresstensor s = sigma.deviator();
	 s.report("\n\n stress deviator s \n\n"); 

	 stresstensor test01 = s("ij") * s("ij");
	 test01.reportshort("\n\n\n test01 \n\n\n ");

	 stresstensor KWtest01 =  s("ij") * s("jk");
	 stresstensor KWtest02 =  KWtest01("ik") * I2("ik");
	 KWtest02.reportshort("\n\n\n KWtest02 \n\n\n "); 


//  static double stressvector[] = {-3.0, -2.0, 5.0, 2.0, 3.0, 5.0};
  static double stressvector[] = {-3.0, -2.0, 5.0, 2.0*sqrt(2.0), 3.0*sqrt(2.0), 5.0*sqrt(2.0)};

// vector constructor
  BJvector stressvector01( 6, stressvector);
  stressvector01.print("sigma","\n\n stressvector01\n");

	 BJvector product = stressvector01.transpose() * stressvector01;
  product.print("product","\n\n product = stressvector01*stressvector01 \n");


::printf("\n\n\n");



// --------------
::printf("\n\n\n --------------------------TENSOR SCALAR MULTIPLICATION ----\n\n");
  tensor EEEE1 = E2;
  tensor EEEE2 = EEEE1*2.0;
  EEEE1.print("E1","tensor EEEE1");
  EEEE2.print("E2","tensor EEEE2");


//###printf("\n\n------------VECTOR AND SKYMATRIX -----------------\n\n");
//###// 1) use of member function full_val is of crucial importance here
//###// and it enable us to use skymatrix as if it is full matrix.
//###// There is however small time overhead but that can be optimized . . .
//###double vect4[] = {1.0, 2.0, 3.0, 4.0, 5.0};
//###
//###matrix res2;
//###vector vector2(5, vect4);
//###second.full_print();
//###vector2.print();
//###res2 = vector2.transpose() * second;
//###res2.print();
//###
//###
//###printf("\n\n------------ MATRIX AND SKYMATRIX -----------------\n\n");
//###// 1) use of member function full_val is of crucial importance here
//###// and it enable us to use skymatrix as if it is full matrix.
//###// There is however small time overhead but that can be optimized . . .
//###
//###matrix big_result;
//###big_result = first_matrix * second;
//###big_result.print();




  exit(1);
//  return 1;
}
