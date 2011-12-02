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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/main.cpp,v $
                                                                        
                                                                        

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>

main() 
{
    Matrix a(2,2);
    a(0,0) = 2.0;
    a(0,1) = 1.0;
    a(1,0) = 1.0;
    a(1,1) = 2.0;
    
    Matrix b(1,1);
    b(0,0) = 10;
    
    Matrix c(3,3);
    c(0,0) = 9.0;
    c(0,1) = 2.0;
    c(0,2) = 2.0;
    c(2,0) = 2.0;
    c(2,1) = 2.0;
    c(2,2) = 9.0;
    c(1,0) = 2.0;
    c(1,1) = 9.0;
    c(1,2) = 2.0;		      
    

    cout << a;
    cout << c;
    ID d(1);    d(0) =  1;
    ID e(1);    e(0) =  2;
    ID f(2);    f(0) =  1; f(1) = 0;
    ID g(2);    g(0) = -1; g(1) = 2;
    ID h(3);    h(0) =  2; h(1) = 1; h(2) = 3;
    ID m(3);    m(0) =  5; m(1) = 4; m(2) = 3;    


    Vector x(3); x(0) = 2; x(1) = 3; x(2) = 4.5;
    

    BandGenLinLapackSolver theSolver;
    BandGenLinSOE sys(6,2,2,theSolver);
    sys.addA(a,f);
    sys.addA(a,g);
    sys.addA(c,h);
    sys.addA(c,m,2.0);    

    sys.addB(x,h);
    sys.addB(x,m,-1);

    cout << "INFO FLAG: " << sys.solve() << endl;
    cout << sys.getX();

    cout << "INFO FLAG: " << sys.solve() << endl;
    cout << sys.getX();
    
    
    
}
