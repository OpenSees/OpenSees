      program  TESTLSRT39
C     All inputs built-in; execute "as is"     
      double precision   x(3), y(3),  Emat(3,3)
      double precision   Em, nu, h
      double precision   Kb(9,9), Kh(9,9), Kt(9,9)
      double precision   e(9), v(9,9), z(9)
      real               tbeg, tend
      integer            n,m, ls(9), NT,nesec
      character*40       status
      data    x /0.,4.08,3.4/, y /0.,-3.44,1.14/
      data    Em /120./, nu /.25/, h /0.125/
      data    m  /9/
      data    ls /1,2,3,4,5,6,7,8,9/
      data    NT /0/
C    
      print '("x=",3F10.4)', x
      print '("y=",3F10.4)', y
      call  GETEM (Em, nu, Emat)
C     Override GETEM output to test anisotropic material
      Emat(1,1)=1703./10
      Emat(1,2)= 382./10
      Emat(2,1)= 382./10
      Emat(1,3)=-420./10
      Emat(3,1)=-420./10
      Emat(2,2)= 800./10
      Emat(2,3)=   5./10
      Emat(3,2)=   5./10
      Emat(3,3)=1056./10
C
      print'('' Constitutive matrix Emat:'')'
      call   MATRIXPRINT (Emat, 3, 3)
      call   LST39RMembMatStiffANDESTemp 
     &      ('OPT', 'BC', x, y, Emat, h, ls, Kb, m, status)
      print '(A)', status
      print '('' Basic membrane stiffness:'')'
      call   MATRIXPRINT (Kb, m, m)
      call   JACOBI  (Kb, m, m, e, .false., v, z)
      print '('' Eigenvalues of Kb:'')'
      call   MATRIXPRINT (e, 1, m)
      call   LST39RMembMatStiffANDESTemp 
     &      ('OPT', 'HC', x, y, Emat, h, ls, Kh, m, status)
      print '(A)', status
      print '('' HO membrane stiffness:'')'
      call   MATRIXPRINT (Kh, m, m)
      call   JACOBI  (Kh, m, m, e, .false., v, z)
      print '('' Eigenvalues of Kh:'')'
      call   MATRIXPRINT (e, 1, m)
      call   LST39RMembMatStiffANDESTemp 
     &      ('OPT', ' C', x, y, Emat, h, ls, Kt, m, status)
      print '('' Total membrane stiffness:'')'
      call   MATRIXPRINT (Kt, m, m)
      call   JACOBI  (Kt, m, m, e, .false., v, z)
      print '('' Eigenvalues of Kt:'')'
      call   MATRIXPRINT (e, 1, m)
C     Timing test.  It is run only if NT>0. SECNDS is system dep     
      if (NT .le. 0)   stop
      call   SECNDS (tbeg)
      do  n = 1,NT
        call   LST39RMembMatStiffANDESTemp 
     &        ('OPT', 'TC', x, y, Emat, h, ls, Kt, m, status)
      end do
      call   SECNDS (tend)
      if (tbeg .eq. tend)                   tend = tbeg + 1.0E-6
      nesec = float(NT)/(tend-tbeg)
      print '("Elements/sec=",I10)', nesec
      end
      subroutine   SECNDS (tim)
      real         t(3), tim
      call etime(t)
      tim =        t(1)
      return
      end
C=DECK MATRIXPRINT
C=PURPOSE Print (m x n) matrix as 2D array, 9 items/line
C=AUTHOR C. A. Felippa, Jun 1968
C=VERSION Dec 1985
C=BLOCK FORTRAN
      subroutine   MATRIXPRINT (a, m, n)
      integer           m, n, i, j, jref
      double precision  a(m,n)
C
      do 2000  jref = 0,n-1,9
        print '(3X,9I8)',(j,j=jref+1,min(jref+9,n))
        do 1500  i = 1,m
          print '(I5,9F8.4)',i,(a(i,j),j=jref+1,min(jref+9,n))
 1500     continue
 2000   continue
      return
      end
C=END FORTRAN
C=DECK JACOBI
C=PURPOSE Jacobi eigensolver for a real symmetric matrix
C=AUTHOR P. S. Jensen, Jan 1975
C=BLOCK USAGE
C
C     Fortran implementation of the Jacobi eigenvalue-vector Algol
C     program by H. Rutishauser,  Handbook for Automatic Computing,
C     Vol II, Springer-Verlag, New York (1971) pp.202-211
C
C     Input Arguments:
C
C       A    symmetric N by N matrix with Fortran dimension NA by NA.
C       N    eigensystem order.
C       NA   Fortran dimension of A (and V if EIGEC is .true.) in the
C            calling programn.
C     EIVEC  logical variable set to .true. if the eigenvectors are
C            to be calculated.
C       Z    scratch vector of length N for rotation accumulation
C
C     Output Arguments:
C
C       A    The upper triangular portion of A is modified by JACOBI.
C       D    a vector of dimension N in which the computed eigenvalues,
C            ordered by decreasing algebraic values, are stored.
C       V    N by N matrix whose columns are set to the computed  eigen-
C            vectors if EIVEC is .true.  The Fortran dimension of V is NA
C            by NA.  These vectors will be orthogonal to working precision.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    JACOBI  (a, na, n, d, eivec, v, z)
      integer   n,na
      double precision  a(na,na), d(n), v(na,na), z(n)
      double precision  dsum,tfac,tol,dp,dq,t,s,c,h,g,th,sm,tresh
      logical eivec
      integer i,j,nm1,k,kk,rot,p,q,pp1
      equivalence (sm,c,th)
      if  (n.le.0)                           return
C-------- initialize
      if (eivec)                             then
C     --- set v to the identity matrix
      do 105  i = 1,n
        do 100  j = 1,n
          v(j,i) = 0.0
  100     continue
        v(i,i) = 1.0
  105   continue
      end if
C     --- set d to the diagonal of a
      dsum = 0.0
      do 115  i = 1,n
        d(i) = a(i,i)
        dsum = dsum + abs(d(i))
 115    continue
      if  (n.eq.1)                           return
C
      rot  = 0
      nm1  = n-1
      tfac = 0.2/n**2
      tol =  dsum*1.0e-12*n
C-------- iteration loop
      do 250 i = 1,50
        do 118  j = 1,n
          z(j) = 0.0
 118      continue
C        --- determine the sum of the abs. values of the off diagonals
         sm = 0.
         do 120 p=1,nm1
            pp1 = p+1
            do 120 q=pp1,n
               sm = abs(a(p,q)) + sm
  120       continue
C        if (sm)                             260,260,130
         if (sm-tol)                         260,260,130
  130    tresh = 0.
         if (i.lt.4)  tresh=tfac*sm
C     ------ loop for one sweep over the matrix ------
         do 230 p=1,nm1
            pp1 = p+1
            do 220 q=pp1,n
               dp = d(p)+z(p)
               dq = d(q)+z(q)
               g = 100.*abs(a(p,q))
               if (i-4)                      170,170,140
C              --- rotation angle size check
  140          if (abs(dp)+g - abs(dp))      170,150,170
  150          if (abs(dq)+g - abs(dq))      170,160,170
  160          a(p,q) = 0.
                                             go to 220
  170          if (abs(a(p,q))-tresh)        220,220,180
C              --- perform a plane rotation
  180          h = dq-dp
               if (abs(h)+g - abs(h))        200,190,200
  190          t = a(p,q)/h
                                             go to 210
C              --- calculate  t = tan(rotation angle, phi)
  200          t = 0.5*h/a(p,q)
               t = 1.0/(t + sign(sqrt(1.0+t**2),t))
C              --- calculate the sin and cos and modify the diagonals
  210          c = 1.0/sqrt(1.0+t**2)
               s = t*c
               h = t*a(p,q)
               z(p)= z(p)-h
               z(q)= z(q)+h
               a(p,q)= 0.
C              --- apply the rotation to a
C                  use th = tan(phi/2) for better accuracy
               th= s/(1.0+c)
               if (p.gt.1)  call SJ2(p-1,  a(1,p),1,   a(1,q),1,   s,th)
               if (pp1.lt.q)call SJ2(q-pp1,a(p,pp1),na,a(pp1,q),1, s,th)
               if (q.lt.n)  call SJ2(n-q,  a(p,q+1),na,a(q,q+1),na,s,th)
C              --- also apply the rotation to v if desired
               if (eivec)   call SJ2(n,    v(1,p),1,   v(1,q),1,   s,th)
C
               rot = 1+rot
  220       continue
  230    continue
         do 240 p=1,n
  240       d(p) = d(p)+z(p)
  250 continue
      print *, '*** Remark: Jacobi eigensolution failed to converge.'
C ---  order eigenvalues by decreasing algebraic value
  260 continue
      do 400  kk = 2,n
        k =   n+2-kk
        do 350  j = 2,k
          if (d(j-1).ge.d(j))          go to 350
          t =        d(j-1)
          d(j-1) =   d(j)
          d(j) =     t
          if (.not.eivec)              go to 350
          do 320  i = 1,n
            t =        v(i,j-1)
            v(i,j-1) = v(i,j)
            v(i,j) =   t
  320       continue
  350     continue
  400   continue
  500 continue
      return
      end
      subroutine SJ2 (n, x, ix, y, iy, s,tau)
C-----------------------------------------------------------------------
C     Apply a jacobi plane rotation of angle theta, where  s = sin(theta)
C     and  tau = s/(1+cos(theta))  are given, to vector pair (x,y)**t.
C     == Restriction ==  iy must be positive.
C-----------------------------------------------------------------------
      double precision  x(*),y(*),g,h,s,tau
      integer  n,ix,iy,jx,jy,niy
C
      if (n .le. 0)             return
  100 jx = 1
      niy= n*iy
      do 110 jy=1,niy,iy
         g = x(jx)
         h = y(jy)
         x(jx) = g-s*(h+g*tau)
         y(jy) = h+s*(g-h*tau)
         jx = ix+jx
  110 continue
      return
      end
C=END FORTRAN
C=DECK GETEM GETEM FORTRAN
C=PURPOSE Form isotropic membrane constitutive matrix Emat
C=AUTHOR C. A. Felippa, April 1966
C=VERSION July 1988
C=BLOCK FORTRAN
      subroutine  GETEM (Em, nu, Emat)
      double precision       Em, nu, Emat(3,3)
      double precision       c
      c =          Em/(1.-nu**2)
      Emat(1,1) =    c
      Emat(1,2) =    nu*c
      Emat(2,1) =    nu*c
      Emat(2,2) =    c
      Emat(3,3) =    0.5d0*(1.-nu)*c
      Emat(1,3) =    0.0
      Emat(2,3) =    0.0
      Emat(3,1) =    0.0
      Emat(3,2) =    0.0
      return
      end
C=END FORTRAN
C=DECK LST39RMembMatStiffANDESTemp 
C=PURPOSE Form material stiffness of 9-DOF memb triangle ANDES template
C=AUTHOR C. A. Felippa, April 2002
C=VERSION January 2003
C=EQUIPMENT Machine independent
C=KEYWORDS finite element
C=KEYWORDS material stiffness matrix ANDES 2002
C=KEYWORDS optimal triangle membrane drilling freedoms
C=BLOCK ABSTRACT
C
C     LST39RMembMatStiffTemp forms the material stiffness matrix 
C     of a 9-dof membrane triangle based on ANDES template.  
C     This is a reformulation of the 1991 Militello-Felippa 
C     element of the same configuration (C. A. Felippa and 
C     C. Militello, Membrane triangles with corner drilling 
C     freedoms: Part II. The ANDES element, FEAD, Vol 12, 
C     189-201, 1992.)
C
C     Emphasis placed on formation speed by unrolling loops.
C     Clocked at 350000 elems/sec on a 3GHz P4 with the
C     optimized Intel f90 compiler.
C
C     Reference for this implementation: CMAME, Vol 192, 
C     2125-2168, 2003.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C     call   LST39RMembMatStiffANDESTemp 
C           (name, options, x, y, Em, h, ls, Km, m, status)
C
C     The inputs are:
C    
C      name     Name of template (see LST39RSignature)
C      options  2-letter character string process options
C               options(1:1)='B'  return only basic stiffness
C               options(1:1)='H'  return only higher order stiffness
C               options(1:1)=' '  return total stiffness
C               options(2:2)='C'  clear full array Km on entry
C                    
C     x,y       (3 x 1) arrays of x,y coordinates of triangle corners
C     Em        (3 x 3) elasticity matrix 
C     h         Element tchickness
C     ls        (9 x 1) array of stiffness location pointers.  To form
C               Km in internal DOF order, set ls={1,2,3,4,5,6,7,8,9}  
C     Km        Input stiffness array dimensioned (m x m).  Cleared
C               on entry if options(2:2)='C'.  If not cleared, stiffness
C               is added to incoming matrix (useful for shell assembly).
C     m         First dimension of Km in calling program.
C
C     The outputs are:
C
C     Km        Output stiffness array with membrane stiffness
C               coefficients added in.  The (i,j)-th entry of the
C               (9 by 9) membrane stiffness is added to Km(k,l), 
C               where k=ls(i) and l=ls(j).  Internal DOF order is
C               ux1,uy1,theta1,ux2,uy2,theta2,ux3,uy3,theta3. If
C               options(2:2)='C' they are added over a cleared array.
C     status    Status character variable.  Blank if no error
C               detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine  LST39RMembMatStiffANDESTemp
     $           (name, options, x, y, Em, h, ls, Km, m, status)
C
C                   A R G U M E N T S
C
      character*(*)  name, options, status
      integer       ls(9), m
      double precision  x(3), y(3), Em(3,3), h, Km(m,m) 
C
C                   L O C A L   V A R I A B L E S
C
      double precision     x21,x12,x32,x23,x13,x31
      double precision     y21,y12,y32,y23,y13,y31
      double precision     Eh(3,3), Ehn(3,3), Te(3,3)
      double precision     Q4(3,3), Q5(3,3), Q6(3,3)
      double precision     Lt(9,3), Kb(9,9), Kh(9,9)
      double precision     S11,S12,S13,S22,S23,S33
      double precision     Ssum1,Ssum2,Ssum3,Ssum123
      double precision     area2, area, area4    
      double precision     a1,a2,a3, b1,b2,b3, c1,c2,c3
      double precision     a4,a5,a6, b4,b5,b6, c4,c5,c6
      double precision     cL,cL3,cL6, tfac,sfac,fK, s1,s2,s3
      double precision     fp(11),alphab,beta0,beta1,beta2,beta3
      double precision     beta4,beta5,beta6,beta7,beta8,beta9
      double precision     ll21,ll32,ll13, f21,f32,f13

      integer              i, j, k, l  
      equivalence          (fp(1), alphab),(fp(2), beta0),(fp(3),beta1)
      equivalence          (fp(4), beta2), (fp(5), beta3),(fp(6),beta4)
      equivalence          (fp(7), beta5), (fp(8), beta6),(fp(9),beta7)
      equivalence          (fp(10),beta8), (fp(11),beta9) 
*23456789012345678901234567890123456789012345678901234567890123456789012
C
C                   L O G I C
C
      status =   ' '
      x21 =      x(2)-x(1)
      x12 =     -x21
      x32 =      x(3)-x(2)
      x23 =     -x32
      x13 =      x(1)-x(3)
      x31 =     -x13
      y21 =      y(2)-y(1)
      y12 =     -y21
      y32 =      y(3)-y(2)
      y23 =     -y32
      y13 =      y(1)-y(3)
      y31 =     -y13
      area2 =    y21*x13 - x21*y13
      if (area2 .le. 0.0d0)      then
            status = 'LST39RMembMatStiffTemplate: Negative area'
        if (area2 .eq. 0.0d0)   
     &      status = 'LST39RMembMatStiffTemplate: Zero area'
        return
      end if
      area =  0.5d0*area2
      area4 = 2.0d0*area2
      call       LST39RSignature (name, Em, fp, status)
*     print '("fpars:",11F7.3)', fp
      if (status(1:1) .ne. ' ')                 return
      do 1200  i = 1,9
         Lt(i,1) =  0.0d0
         Lt(i,2) =  0.0d0
         do 1100  j = 1,9
            Kb(i,j) =  0.0d0
            Kh(i,j) =  0.0d0
 1100       continue
 1200    continue
      do 1500  i = 1,3
         do 1500  j = 1,3
           Eh(i,j) =  h*Em(i,j)
 1500    continue
C
      if (options(1:1) .ne. 'H')                  then
        cL =        0.5d0/area
        cL3 =       alphab*cL/3.0d0
        cL6 =       alphab*cL/6.0d0
        Lt(1,1) =   y23*cL
        Lt(1,3) =   x32*cL
        Lt(2,2) =   x32*cL
        Lt(2,3) =   y23*cL
        Lt(4,1) =   y31*cL
        Lt(4,3) =   x13*cL
        Lt(5,2) =   x13*cL
        Lt(5,3) =   y31*cL
        Lt(7,1) =   y12*cL
        Lt(7,3) =   x21*cL
        Lt(8,2) =   x21*cL
        Lt(8,3) =   y12*cL
        Lt(3,1) =   y23*(y13-y21)   *cL6
        Lt(3,2) =   x32*(x31-x12)   *cL6
        Lt(3,3) =  (x31*y13-x12*y21)*cL3
        Lt(6,1) =   y31*(y21-y32)   *cL6
        Lt(6,2) =   x13*(x12-x23)   *cL6
        Lt(6,3) =  (x12*y21-x23*y32)*cL3
        Lt(9,1) =   y12*(y32-y13)   *cL6
        Lt(9,2) =   x21*(x23-x31)   *cL6
        Lt(9,3) =  (x23*y32-x31*y13)*cL3
        do 2000  j = 1,9
          s1 =   area*(Eh(1,1)*Lt(j,1)+Eh(1,2)*Lt(j,2)+Eh(1,3)*Lt(j,3))
          s2 =   area*(Eh(1,2)*Lt(j,1)+Eh(2,2)*Lt(j,2)+Eh(2,3)*Lt(j,3))
          s3 =   area*(Eh(1,3)*Lt(j,1)+Eh(2,3)*Lt(j,2)+Eh(3,3)*Lt(j,3))
          do 1800  i = 1,j
            Kb(i,j) =  Kb(i,j) + (s1*Lt(i,1) + s2*Lt(i,2) + s3*Lt(i,3))
            Kb(j,i) =  Kb(i,j)
 1800     continue
 2000   continue
      end if
C
      if (options(1:1) .ne. 'B')                  then
        ll21 =       x21*x21+y21*y21
        ll32 =       x32*x32+y32*y32
        ll13 =       x13*x13+y13*y13
        tfac =       1.0d0/(area2*area2)
        Te(1,1) =    tfac*y23*y13*ll21
        Te(1,2) =    tfac*y31*y21*ll32
        Te(1,3) =    tfac*y12*y32*ll13
        Te(2,1) =    tfac*x23*x13*ll21
        Te(2,2) =    tfac*x31*x21*ll32
        Te(2,3) =    tfac*x12*x32*ll13
        Te(3,1) =    tfac*(y23*x31+x32*y13)*ll21
        Te(3,2) =    tfac*(y31*x12+x13*y21)*ll32
        Te(3,3) =    tfac*(y12*x23+x21*y32)*ll13
*23456789012345678901234567890123456789012345678901234567890123456789012
        a1 =        Eh(1,1)*Te(1,1)+Eh(1,2)*Te(2,1)+Eh(1,3)*Te(3,1)
        b1 =        Eh(2,1)*Te(1,1)+Eh(2,2)*Te(2,1)+Eh(2,3)*Te(3,1)
        c1 =        Eh(3,1)*Te(1,1)+Eh(3,2)*Te(2,1)+Eh(3,3)*Te(3,1)
        Ehn(1,1) =  Te(1,1)*a1+Te(2,1)*b1+Te(3,1)*c1
        Ehn(2,1) =  Te(1,2)*a1+Te(2,2)*b1+Te(3,2)*c1
        Ehn(1,2) =  Ehn(2,1)
        Ehn(3,1) =  Te(1,3)*a1+Te(2,3)*b1+Te(3,3)*c1
        Ehn(1,3) =  Ehn(3,1)
        a2 =        Eh(1,1)*Te(1,2)+Eh(1,2)*Te(2,2)+Eh(1,3)*Te(3,2)
        b2 =        Eh(2,1)*Te(1,2)+Eh(2,2)*Te(2,2)+Eh(2,3)*Te(3,2)
        c2 =        Eh(3,1)*Te(1,2)+Eh(3,2)*Te(2,2)+Eh(3,3)*Te(3,2)
        Ehn(2,2) =  Te(1,2)*a2+Te(2,2)*b2+Te(3,2)*c2
        Ehn(3,2) =  Te(1,3)*a2+Te(2,3)*b2+Te(3,3)*c2
        Ehn(2,3) =  Ehn(3,2)
        a3 =        Eh(1,1)*Te(1,3)+Eh(1,2)*Te(2,3)+Eh(1,3)*Te(3,3)
        b3 =        Eh(2,1)*Te(1,3)+Eh(2,2)*Te(2,3)+Eh(2,3)*Te(3,3)
        c3 =        Eh(3,1)*Te(1,3)+Eh(3,2)*Te(2,3)+Eh(3,3)*Te(3,3)
        Ehn(3,3) =  Te(1,3)*a3+Te(2,3)*b3+Te(3,3)*c3
        f21 =     area/(3.d0*ll21)
        f32 =     area/(3.d0*ll32)
        f13 =     area/(3.d0*ll13)
        Q4(1,1) = (beta1+beta9)*f21
        Q4(1,2) = (beta2+beta7)*f21
        Q4(1,3) = (beta3+beta8)*f21
        Q4(2,1) = (beta3+beta4)*f32
        Q4(2,2) = (beta1+beta5)*f32
        Q4(2,3) = (beta2+beta6)*f32
        Q4(3,1) = (beta6+beta7)*f13
        Q4(3,2) = (beta4+beta8)*f13
        Q4(3,3) = (beta5+beta9)*f13
        Q5(1,1) = (beta5+beta9)*f21
        Q5(1,2) = (beta6+beta7)*f21
        Q5(1,3) = (beta4+beta8)*f21
        Q5(2,1) = (beta3+beta8)*f32
        Q5(2,2) = (beta1+beta9)*f32
        Q5(2,3) = (beta2+beta7)*f32
        Q5(3,1) = (beta2+beta6)*f13
        Q5(3,2) = (beta3+beta4)*f13
        Q5(3,3) = (beta1+beta5)*f13
        Q6(1,1) = (beta1+beta5)*f21
        Q6(1,2) = (beta2+beta6)*f21
        Q6(1,3) = (beta3+beta4)*f21
        Q6(2,1) = (beta4+beta8)*f32
        Q6(2,2) = (beta5+beta9)*f32 
        Q6(2,3) = (beta6+beta7)*f32
        Q6(3,1) = (beta2+beta7)*f13
        Q6(3,2) = (beta3+beta8)*f13
        Q6(3,3) = (beta1+beta9)*f13
        sfac =    0.75d0*beta0*area
        a4 =     Ehn(1,1)*Q4(1,1)+Ehn(1,2)*Q4(2,1)+Ehn(1,3)*Q4(3,1)
        b4 =     Ehn(2,1)*Q4(1,1)+Ehn(2,2)*Q4(2,1)+Ehn(2,3)*Q4(3,1)
        c4 =     Ehn(3,1)*Q4(1,1)+Ehn(3,2)*Q4(2,1)+Ehn(3,3)*Q4(3,1)
        a5 =     Ehn(1,1)*Q5(1,1)+Ehn(1,2)*Q5(2,1)+Ehn(1,3)*Q5(3,1)
        b5 =     Ehn(2,1)*Q5(1,1)+Ehn(2,2)*Q5(2,1)+Ehn(2,3)*Q5(3,1)
        c5 =     Ehn(3,1)*Q5(1,1)+Ehn(3,2)*Q5(2,1)+Ehn(3,3)*Q5(3,1)
        a6 =     Ehn(1,1)*Q6(1,1)+Ehn(1,2)*Q6(2,1)+Ehn(1,3)*Q6(3,1)
        b6 =     Ehn(2,1)*Q6(1,1)+Ehn(2,2)*Q6(2,1)+Ehn(2,3)*Q6(3,1)
        c6 =     Ehn(3,1)*Q6(1,1)+Ehn(3,2)*Q6(2,1)+Ehn(3,3)*Q6(3,1)
        S11 =   (Q4(1,1)*a4+Q4(2,1)*b4+Q4(3,1)*c4+Q5(1,1)*a5+Q5(2,1)*b5
     &           +Q5(3,1)*c5+Q6(1,1)*a6+Q6(2,1)*b6+Q6(3,1)*c6)*sfac
        S12 =   (Q4(1,2)*a4+Q4(2,2)*b4+Q4(3,2)*c4+Q5(1,2)*a5+Q5(2,2)*b5
     &           +Q5(3,2)*c5+Q6(1,2)*a6+Q6(2,2)*b6+Q6(3,2)*c6)*sfac
        S13 =   (Q4(1,3)*a4+Q4(2,3)*b4+Q4(3,3)*c4+Q5(1,3)*a5+Q5(2,3)*b5
     &           +Q5(3,3)*c5+Q6(1,3)*a6+Q6(2,3)*b6+Q6(3,3)*c6)*sfac
        a4 =     Ehn(1,1)*Q4(1,2)+Ehn(1,2)*Q4(2,2)+Ehn(1,3)*Q4(3,2)
        b4 =     Ehn(2,1)*Q4(1,2)+Ehn(2,2)*Q4(2,2)+Ehn(2,3)*Q4(3,2)
        c4 =     Ehn(3,1)*Q4(1,2)+Ehn(3,2)*Q4(2,2)+Ehn(3,3)*Q4(3,2)
        a5 =     Ehn(1,1)*Q5(1,2)+Ehn(1,2)*Q5(2,2)+Ehn(1,3)*Q5(3,2)
        b5 =     Ehn(2,1)*Q5(1,2)+Ehn(2,2)*Q5(2,2)+Ehn(2,3)*Q5(3,2)
        c5 =     Ehn(3,1)*Q5(1,2)+Ehn(3,2)*Q5(2,2)+Ehn(3,3)*Q5(3,2)
        a6 =     Ehn(1,1)*Q6(1,2)+Ehn(1,2)*Q6(2,2)+Ehn(1,3)*Q6(3,2)
        b6 =     Ehn(2,1)*Q6(1,2)+Ehn(2,2)*Q6(2,2)+Ehn(2,3)*Q6(3,2)
        c6 =     Ehn(3,1)*Q6(1,2)+Ehn(3,2)*Q6(2,2)+Ehn(3,3)*Q6(3,2)
        S22 =   (Q4(1,2)*a4+Q4(2,2)*b4+Q4(3,2)*c4+Q5(1,2)*a5+Q5(2,2)*b5
     &           +Q5(3,2)*c5+Q6(1,2)*a6+Q6(2,2)*b6+Q6(3,2)*c6)*sfac
        S23 =   (Q4(1,3)*a4+Q4(2,3)*b4+Q4(3,3)*c4+Q5(1,3)*a5+Q5(2,3)*b5
     &           +Q5(3,3)*c5+Q6(1,3)*a6+Q6(2,3)*b6+Q6(3,3)*c6)*sfac
        a4 =     Ehn(1,1)*Q4(1,3)+Ehn(1,2)*Q4(2,3)+Ehn(1,3)*Q4(3,3)
        b4 =     Ehn(2,1)*Q4(1,3)+Ehn(2,2)*Q4(2,3)+Ehn(2,3)*Q4(3,3)
        c4 =     Ehn(3,1)*Q4(1,3)+Ehn(3,2)*Q4(2,3)+Ehn(3,3)*Q4(3,3)
        a5 =     Ehn(1,1)*Q5(1,3)+Ehn(1,2)*Q5(2,3)+Ehn(1,3)*Q5(3,3)
        b5 =     Ehn(2,1)*Q5(1,3)+Ehn(2,2)*Q5(2,3)+Ehn(2,3)*Q5(3,3)
        c5 =     Ehn(3,1)*Q5(1,3)+Ehn(3,2)*Q5(2,3)+Ehn(3,3)*Q5(3,3)
        a6 =     Ehn(1,1)*Q6(1,3)+Ehn(1,2)*Q6(2,3)+Ehn(1,3)*Q6(3,3)
        b6 =     Ehn(2,1)*Q6(1,3)+Ehn(2,2)*Q6(2,3)+Ehn(2,3)*Q6(3,3)
        c6 =     Ehn(3,1)*Q6(1,3)+Ehn(3,2)*Q6(2,3)+Ehn(3,3)*Q6(3,3)
        S33 =   (Q4(1,3)*a4+Q4(2,3)*b4+Q4(3,3)*c4+Q5(1,3)*a5+Q5(2,3)*b5
     &           +Q5(3,3)*c5+Q6(1,3)*a6+Q6(2,3)*b6+Q6(3,3)*c6)*sfac
C
*23456789012345678901234567890123456789012345678901234567890123456789012
*       print'("S mtx:"/(3F10.5))',S11,S12,S13,S21,S22,S23,S31,S32,S33
        Ssum1 =    (S11+S12+S13)/area4
        Ssum2 =    (S12+S22+S23)/area4
        Ssum3 =    (S13+S23+S33)/area4
        Ssum123 =  (Ssum1+Ssum2+Ssum3)/area4
        Kh(1,1) =  Ssum123*x32*x32
        Kh(1,2) =  Ssum123*x32*y32
        Kh(1,3) =  Ssum1*x32
        Kh(1,4) =  Ssum123*x13*x32
        Kh(1,5) =  Ssum123*x32*y13
        Kh(1,6) =  Ssum2*x32
        Kh(1,7) =  Ssum123*x21*x32
        Kh(1,8) =  Ssum123*x32*y21
        Kh(1,9) =  Ssum3*x32
        Kh(2,2) =  Ssum123*y32*y32
        Kh(2,3) =  Ssum1*y32
        Kh(2,4) =  Ssum123*x13*y32
        Kh(2,5) =  Ssum123*y13*y32
        Kh(2,6) =  Ssum2*y32
        Kh(2,7) =  Ssum123*x21*y32
        Kh(2,8) =  Ssum123*y21*y32
        Kh(2,9) =  Ssum3*y32
        Kh(3,3) =  S11
        Kh(3,4) =  Ssum1*x13
        Kh(3,5) =  Ssum1*y13
        Kh(3,6) =  S12
        Kh(3,7) =  Ssum1*x21
        Kh(3,8) =  Ssum1*y21
        Kh(3,9) =  S13
        Kh(4,4) =  Ssum123*x13*x13
        Kh(4,5) =  Ssum123*x13*y13
        Kh(4,6) =  Ssum2*x13
        Kh(4,7) =  Ssum123*x13*x21
        Kh(4,8) =  Ssum123*x13*y21
        Kh(4,9) =  Ssum3*x13
        Kh(5,5) =  Ssum123*y13*y13
        Kh(5,6) =  Ssum2*y13
        Kh(5,7) =  Ssum123*x21*y13
        Kh(5,8) =  Ssum123*y13*y21
        Kh(5,9) =  Ssum3*y13
        Kh(6,6) =  S22
        Kh(6,7) =  Ssum2*x21
        Kh(6,8) =  Ssum2*y21
        Kh(6,9) =  S23
        Kh(7,7) =  Ssum123*x21*x21
        Kh(7,8) =  Ssum123*x21*y21
        Kh(7,9) =  Ssum3*x21
        Kh(8,8) =  Ssum123*y21*y21         
        Kh(8,9) =  Ssum3*y21
        Kh(9,9) =  S33
        do 3200  i = 2,9
           do 3100  j = 1,i-1
              Kh(i,j) = Kh(j,i)
 3100         continue
 3200      continue
      end if
C 
      fK =      1.0d0
      if (options(2:2) .eq. 'C')                 fK = 0.0d0
       do 4500  i = 1,9
          k = ls(i)
          do 4500  j = 1,i
             l = ls(j)
             Km(k,l) = fK*Km(k,l) + Kb(i,j) + Kh(i,j)
             Km(l,k) =    Km(k,l)
 4500     continue 
      return
      end
C=END FORTRAN
C=DECK LST39RSignature 
C=PURPOSE Return signature of 9-DOF memb triangle ANDES template
C=AUTHOR C. A. Felippa, April 2002
C=VERSION January 2003
C=EQUIPMENT Machine independent
C=KEYWORDS finite element
C=KEYWORDS material stiffness matrix ANDES template signature
C=KEYWORDS triangle membrane drilling freedoms
C=BLOCK ABSTRACT
C
C     LST39RSignature returns the ANDES template signature for a
C     LST-3/9 membrane triangle, given a template name.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C     call   LST39RMembMatStiffANDESTemp 
C           (name, Em, fpars, status)
C
C     The inputs are:
C    
C      name     Name of template.  Only 'OPT' implemented here
C               For other names see Mathemtica code.
C
C     The outputs are:
C
C     fpars     Free parameters making up the template signature
C     status    Status character variable.  Blank if no error
C               detected.
C
C=END USAGE
      subroutine  LST39RSignature
     $           (name, Em, fpars, status)
C
C                   A R G U M E N T S
C
      character*(*)  name, status
      double precision  Em(3,3), fpars(10)
C
C                   L O C A L   V A R I A B L E S
C
      double precision     fpOPT(11), beta0, Edet
      double precision     E11,E22,E33,E12,E13,E23,E11C11
      integer              i
      save                 fpOPT
      data                 fpOPT
     &  /1.5d0,0.d0,1.d0,2.d0,1.d0,0.d0,1.d0,-1.d0,-1.d0,-1.d0,-2.d0/ 
  
C
C                   L O G I C
C
*23456789012345678901234567890123456789012345678901234567890123456789012 
      status = ' '
      if (name(1:3) .eq. 'OPT')             then
        do 1500  i = 1,11
           fpars(i) =   fpOPT(i)
 1500      continue
        E11 =  Em(1,1)
        E22 =  Em(2,2)
        E33 =  Em(3,3)
        E12 =  Em(1,2)
        E13 =  Em(1,3)
        E23 =  Em(2,3)
        Edet=E11*E22*E33+2*E12*E13*E23-E11*E23**2-E22*E13**2-E33*E12**2
        E11C11=(-5*E11*E12**2-6*E12**3-3*E11*E13**2+14*E12*E13**2+
     &   5*E11**2*E22+6*E11*E12*E22-5*E12**2*E22-75*E13**2*E22+
     &   5*E11*E22**2-14*E11*E13*E23+92*E12*E13*E23-14*E13*E22*E23-
     &   75*E11*E23**2+14*E12*E23**2-3*E22*E23**2+(3*E11**2+82*E11*E22+
     &   3*E22**2-4*(6*E12**2+5*E13**2-6*E13*E23+5*E23**2))*E33+
     &   4*(5*E11-6*E12+5*E22)*E33**2)/(128*Edet)
        beta0 = max(2.d0/E11C11-1.5d0,0.01)
*       print '("Edet,E11C11,beta0=",3G16.7)', Edet,E11C11,beta0
        fpars(2) = beta0
        return
      end if
      status = 'LST39RSignature: Illegal template name'
      return
      end
C=END FORTRAN
