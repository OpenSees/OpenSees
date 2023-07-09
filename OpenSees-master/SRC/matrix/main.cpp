/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.3 $
// $Date: 2005-12-14 23:49:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/main.cpp,v $
                                                                        
                                                                        
#include "Vector.h"
#include "ID.h"
#include "Matrix.h"
#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

int main()
{

  ID data(0);

  opserr << data << endln;
  opserr << data.insert(1) << endln;
  opserr << data << endln;
  opserr << data.insert(2) << endln;
  opserr << data << endln;
  opserr << data.insert(4) << endln;
  opserr << data << endln;
  opserr << data.insert(5) << endln;
  opserr << data << endln;
  opserr << data.insert(6) << endln;
  opserr << data << endln;
  opserr << data.insert(3) << endln;
  opserr << data << endln;
  opserr << data.insert(-1) << endln;
  opserr << data << endln;

  opserr << data.insert(3) << endln;
  opserr << data << endln;
  exit(0);

  opserr.setPrecision(14);
  opserr.setFloatField(SCIENTIFIC);

  Matrix a(3,3);
  a(0,0)= 1.35691200000000e+08; a(0,1)= 1.49035621973463e+07; a(0,2) =  -3.72529029846191e-09;
  a(1,0) = 1.49035621973463e+07; a(1,1)= 1.63692388430620e+06; a(1,2)= -4.65661287307739e-10;
  a(2,0) =-3.72529029846191e-09; a(2,1)=  -4.65661287307739e-10; a(2,2)= 3.19209024960000e+07;

  opserr << a;

  Matrix b(3,3);

  a.Invert(b);
  opserr << b;

  opserr << a*b;

  Matrix I(3,3);
  for (int i=0; i<3; i++)
    I(i,i) = 1.0;

  opserr << a.Solve(I,b);
  
  opserr << b;
  opserr << a*b;
  opserr << I;

  /*
    // Vector  TEST
    cout << "\n***** Vector TEST *****  \n";
    Vector a;
    cout << "Vector a;     a -> " << a;
    Vector b(5);
    cout << "Vector b(5);  b -> " << b;    
    Vector c(b);
    cout << "Vector c(b);  c -> " << c;        
    Vector d = a;
    cout << "Vector d = a; d -> " << d;        
    d = c;
    cout << "       d = c; d -> " << d;            
    d+=5.0;
    cout << "       d +=5.0; d -> " << d;                
    cout << "                c -> " << c;
    d-=3.0;
    cout << "       d -=3.0; d -> " << d;                    
    d*=6.0;
    cout << "       d *=6.0; d -> " << d;                    
    d/=4.0;
    cout << "       d /=4.0; d -> " << d;                        
    c = d + 2;
    cout << "     c = d + 2; c -> " << c;                            
    cout << "                d -> " << d;

    c = c * 2;
    cout << "     c = c * 2; c -> " << c;                            
    c = c - 2;
    cout << "     c = c - 2; c -> " << c;                            
    c = c / 2;
    cout << "     c = c / 2; c -> " << c;                            
    double x = c^c;
    cout << "double x = c^c; x -> " << x << "\n";
    x = c^d;
    cout << "       x = c^d; x -> " << x << "\n";                                
    a = c * (c^d);
    cout << "   a = c*(c^d); a -> " << a;

    cout << "\n\n\n***** ID TEST ***** \n";    	    
    ID id1(2);
    cout << "     ID id1(2): id1 -> " << id1;
    ID id2(3);
    cout << "     ID id2(3): id2 -> " << id2;    

    id2(0) = 2; id2(1) = 0; id2(2) = 3;
    cout << "      revised : id2 -> " << id2;        
    Vector y(5); y(0) = 1.0; y(1) = 2.0; y(2) = 3.0; y(3) = 4.0; y(4) = 5.0;
    cout << "              : y ->   " << y;
    
    Vector z = y(id2);
    cout << "Vector &z=y(id2): z ->   " << z;

  */

  /*
    Matrix ZA = A(id2,id2);
    cout << "Matrix &ZA=A(id2,id2): ZA -> " << ZA;
    Matrix &ZB = A(id1,id2);
    cout << "Matrix &ZB=A(id1,id2): ZB -> " << ZB;    
    ZB = A(id2,id1);
    cout << "        ZB=A(id2,id1): ZB -> " << ZB;    
    ZB = A(id1,id1);
    cout << "        ZB=A(id1,id1): ZB -> " << ZB;        
    */

  /*
    cout << "\n\n\n***** Matrix TEST ***** \n";    	
    Matrix A(2,2);
    cout << "      Matrix A; A -> " << A;
    Matrix B(3,2);
    cout << " Matrix B(3,2); B -> " << B;    
    Matrix C(B);
    cout << "   Matrix C(B); C -> " << C;        
    
    C+=2.0;
    cout << "       C +=2.0; C -> " << C;
    C*=6.0;
    cout << "       C *=6.0; C -> " << C;
    C/=4.0;
    cout << "       C /=4.0; C -> " << C;
    C-=1.0;
    cout << "       C -=1.0; C -> " << C;
    B = C *2.0;
    cout << "     B = C*2.0; B -> " << B;    
    B = C /2.0;
    cout << "     B = C/2.0; B -> " << B;    
    B = C -2.0;
    cout << "     B = C-2.0; B -> " << B;    
    B = C +2.0;
    cout << "     B = C+2.0; B -> " << B;    
    A = C^B;
    cout << "       A = C^B; A -> " << A;        
    A(1,1) += 6;
    cout << "     A(1,1)+=6; A -> " << A;            
    
    Vector e(2); e(0) = 72; e(1) = 84;
    cout << "                e -> " << e;
    //    a = e / A;
    //    cout << "      a= e / A: a -> " << a;
    cout << "                A -> " << A;
    
    Matrix E(2,5);
    E += 2.0;
    E(1,3) = 5.0;
    cout << "               E -> " << E;

    cout << " A * E -> " << A*E;

    Matrix G(2,5);
    cout << " G -> " << G;
    G.addMatrixProduct(1.0,A,E,-1.0);
    cout << " G.addMatrixProduct(A,E,-1.0) -> " << G;
    G.addMatrixProduct(1.0,A,E,-1.0);
    cout << " G.addMatrixProduct(A,E,-1.0) -> " << G;

    Matrix F(5,5);
    F = E^A*E;
    cout << "    F = E^A*E: F -> " << F;    

    F.addMatrixTripleProduct(0.0,E,A,1.0);
    cout << "    F.addMatrixTripleProduct(0.0,E,A,1.0): F -> " << F;    
    F.addMatrixTripleProduct(1.0,E,A,1.0);
    cout << "    F.addMatrixTripleProduct(E,A,1.0): F -> " << F;    

    Matrix AA(2,2); AA(0,0) = 1; AA(0,1) = 2; AA(1,0) = 4; AA(1,1) = 9;
    Vector b(2), x(2); b(0) = 1; b(1) = 2;

    opserr << " AA -> " << AA;
    opserr << "  b -> " << b;

    AA.Solve(b,x);

    opserr << "AA.solve(b,x) x-> " << x;
    opserr << "AA * x -> " << AA*x;



    Matrix bb(2,3); bb(0,0) = 1; bb(1,0) = 2;
    bb(0,1) = 2; bb(1,1) = 4;
    bb(0,2) = 11; bb(1,2) = 12;
    Matrix xx(2,3);

    opserr << " AA -> " << AA;
    opserr << " bb -> " << bb;

    AA.Solve(bb,xx);

    opserr << "AA.solve(bb,xx) xx-> " << xx;
    opserr << "AA * xx -> " << AA*xx;
  */
}



