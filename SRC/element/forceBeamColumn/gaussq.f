c  To get dgamma,  "send dgamma from fnlib".
c  To get d1mach, mail netlib
c       send d1mach from core
c
      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
c
c           this set of routines computes the nodes t(j) and weights
c        w(j) for gaussian-type quadrature rules with pre-assigned
c        nodes.  these are used when one wishes to approximate
c
c                 integral (from a to b)  f(x) w(x) dx
c
c                              n
c        by                   sum w  f(t )
c                             j=1  j    j
c
c        (note w(x) and w(j) have no connection with each other.)
c        here w(x) is one of six possible non-negative weight
c        functions (listed below), and f(x) is the
c        function to be integrated.  gaussian quadrature is particularly
c        useful on infinite intervals (with appropriate weight
c        functions), since then other techniques often fail.
c
c           associated with each weight function w(x) is a set of
c        orthogonal polynomials.  the nodes t(j) are just the zeroes
c        of the proper n-th degree polynomial.
c
c     input parameters (all real numbers are in double precision)
c
c        kind     an integer between 1 and 6 giving the type of
c                 quadrature rule:
c
c        kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)
c        kind = 2:  chebyshev quadrature of the first kind
c                   w(x) = 1/sqrt(1 - x*x) on (-1, +1)
c        kind = 3:  chebyshev quadrature of the second kind
c                   w(x) = sqrt(1 - x*x) on (-1, 1)
c        kind = 4:  hermite quadrature, w(x) = exp(-x*x) on
c                   (-infinity, +infinity)
c        kind = 5:  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
c                   beta on (-1, 1), alpha, beta .gt. -1.
c                   note: kind=2 and 3 are a special case of this.
c        kind = 6:  generalized laguerre quadrature, w(x) = exp(-x)*
c                   x**alpha on (0, +infinity), alpha .gt. -1
c
c        n        the number of points used for the quadrature rule
c        alpha    real parameter used only for gauss-jacobi and gauss-
c                 laguerre quadrature (otherwise use 0.d0).
c        beta     real parameter used only for gauss-jacobi quadrature--
c                 (otherwise use 0.d0)
c        kpts     (integer) normally 0, unless the left or right end-
c                 point (or both) of the interval is required to be a
c                 node (this is called gauss-radau or gauss-lobatto
c                 quadrature).  then kpts is the number of fixed
c                 endpoints (1 or 2).
c        endpts   real array of length 2.  contains the values of
c                 any fixed endpoints, if kpts = 1 or 2.
c        b        real scratch array of length n
c
c     output parameters (both double precision arrays of length n)
c
c        t        will contain the desired nodes.
c        w        will contain the desired weights w(j).
c
c     underflow may sometimes occur, but is harmless.
c
c     references
c        1.  golub, g. h., and welsch, j. h., "calculation of gaussian
c            quadrature rules," mathematics of computation 23 (april,
c            1969), pp. 221-230.
c        2.  golub, g. h., "some modified matrix eigenvalue problems,"
c            siam review 15 (april, 1973), pp. 318-334 (section 7).
c        3.  stroud and secrest, gaussian quadrature formulas, prentice-
c            hall, englewood cliffs, n.j., 1966.
c
c        original version 20 jan 1975 from stanford
c        modified 21 dec 1983 by eric grosse
c          imtql2 => gausq2
c          hex constant => d1mach (from core library)
c          compute pi using datan
c          removed accuracy claims, description of method
c          added single precision version
c
      double precision b(n), t(n), w(n), endpts(2), muzero, t1,
     x gam, solve, dsqrt, alpha, beta
c
      call class (kind, n, alpha, beta, b, t, muzero)
c
c           the matrix of coefficients is assumed to be symmetric.
c           the array t contains the diagonal elements, the array
c           b the off-diagonal elements.
c           make appropriate changes in the lower right 2 by 2
c           submatrix.
c
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50
c
c           if kpts=1, only t(n) must be changed
c
      t(n) = solve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
c
c           if kpts=2, t(n) and b(n-1) must be recomputed
c
   50 gam = solve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(solve(endpts(2), n, t, b) - gam))
      b(n-1) = dsqrt(t1)
      t(n) = endpts(1) + gam*t1
c
c           note that the indices of the elements of b run from 1 to n-1
c           and thus the value of b(n) is arbitrary.
c           now compute the eigenvalues of the symmetric tridiagonal
c           matrix, which has been modified as necessary.
c           the method used is a ql-type method with origin shifting
c
  100 w(1) = 1.0d0
      do 105 i = 2, n
  105    w(i) = 0.0d0
c
      call gausq2 (n, t, b, w, ierr)
      do 110 i = 1, n
  110    w(i) = muzero * w(i) * w(i)
c
      return
      end
c
c
c
      double precision function solve(shift, n, a, b)
c
c       this procedure performs elimination to solve for the
c       n-th component of the solution delta to the equation
c
c             (jn - shift*identity) * delta  = en,
c
c       where en is the vector of all zeroes except for 1 in
c       the n-th position.
c
c       the matrix jn is symmetric tridiagonal, with diagonal
c       elements a(i), off-diagonal elements b(i).  this equation
c       must be solved to obtain the appropriate changes in the lower
c       2 by 2 submatrix of coefficients for orthogonal polynomials.
c
c
      double precision shift, a(n), b(n), alpha
c
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
   10    alpha = a(i) - shift - b(i-1)**2/alpha
      solve = 1.0d0/alpha
      return
      end
c
c
c
      subroutine class(kind, n, alpha, beta, b, a, muzero)
c
c           this procedure supplies the coefficients a(j), b(j) of the
c        recurrence relation
c
c             b p (x) = (x - a ) p   (x) - b   p   (x)
c              j j            j   j-1       j-1 j-2
c
c        for the various classical (normalized) orthogonal polynomials,
c        and the zero-th moment
c
c             muzero = integral w(x) dx
c
c        of the given polynomial's weight function w(x).  since the
c        polynomials are orthonormalized, the tridiagonal matrix is
c        guaranteed to be symmetric.
c
c           the input parameter alpha is used only for laguerre and
c        jacobi polynomials, and the parameter beta is used only for
c        jacobi polynomials.  the laguerre and jacobi polynomials
c        require the gamma function.
c
      double precision a(n), b(n), muzero, alpha, beta
      double precision abi, a2b2, pi, dsqrt, ab
c      double precision abi, a2b2, dgamma, pi, dsqrt, ab
c
      pi = 4.0d0 * datan(1.0d0)
      nm1 = n - 1
      go to (10, 20, 30, 40, 50, 60), kind
c
c              kind = 1:  legendre polynomials p(x)
c              on (-1, +1), w(x) = 1.
c
   10 muzero = 2.0d0
      do 11 i = 1, nm1
         a(i) = 0.0d0
         abi = i
   11    b(i) = abi/dsqrt(4*abi*abi - 1.0d0)
      a(n) = 0.0d0
      return
c
c              kind = 2:  chebyshev polynomials of the first kind t(x)
c              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
c
   20 muzero = pi
      do 21 i = 1, nm1
         a(i) = 0.0d0
   21    b(i) = 0.5d0
      b(1) = dsqrt(0.5d0)
      a(n) = 0.0d0
      return
c
c              kind = 3:  chebyshev polynomials of the second kind u(x)
c              on (-1, +1), w(x) = sqrt(1 - x*x)
c
   30 muzero = pi/2.0d0
      do 31 i = 1, nm1
         a(i) = 0.0d0
   31    b(i) = 0.5d0
      a(n) = 0.0d0
      return
c
c              kind = 4:  hermite polynomials h(x) on (-infinity,
c              +infinity), w(x) = exp(-x**2)
c
   40 muzero = dsqrt(pi)
      do 41 i = 1, nm1
         a(i) = 0.0d0
   41    b(i) = dsqrt(i/2.0d0)
      a(n) = 0.0d0
      return
c
c              kind = 5:  jacobi polynomials p(alpha, beta)(x) on
c              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
c              beta greater than -1
c
   50 ab = alpha + beta
c      abi = 2.0d0 + ab
c      muzero = 2.0d0 ** (ab + 1.0d0) * dgamma(alpha + 1.0d0) * dgamma(
c     x beta + 1.0d0) / dgamma(abi)
c      a(1) = (beta - alpha)/abi
c      b(1) = dsqrt(4.0d0*(1.0d0 + alpha)*(1.0d0 + beta)/((abi + 1.0d0)*
c     1  abi*abi))
c      a2b2 = beta*beta - alpha*alpha
c      do 51 i = 2, nm1
c         abi = 2.0d0*i + ab
c         a(i) = a2b2/((abi - 2.0d0)*abi)
c   51    b(i) = dsqrt (4.0d0*i*(i + alpha)*(i + beta)*(i + ab)/
c     1   ((abi*abi - 1)*abi*abi))
c      abi = 2.0d0*n + ab
c      a(n) = a2b2/((abi - 2.0d0)*abi)
      return
c
c              kind = 6:  laguerre polynomials l(alpha)(x) on
c              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
c              than -1.
c
 60   muzero = 1.0
c   60 muzero = dgamma(alpha + 1.0d0)
c      do 61 i = 1, nm1
c         a(i) = 2.0d0*i - 1.0d0 + alpha
c   61    b(i) = dsqrt(i*(i + alpha))
c      a(n) = 2.0d0*n - 1 + alpha
      return
      end
c
c
      subroutine gausq2(n, d, e, z, ierr)
c
c     this subroutine is a translation of an algol procedure,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c     this is a modified version of the 'eispack' routine imtql2.
c
c     this subroutine finds the eigenvalues and first components of the
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
c     method.
c
c     on input:
c
c        n is the order of the matrix;
c
c        d contains the diagonal elements of the input matrix;
c
c        e contains the subdiagonal elements of the input matrix
c          in its first n-1 positions.  e(n) is arbitrary;
c
c        z contains the first row of the identity matrix.
c
c      on output:
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1, 2, ..., ierr-1;
c
c        e has been destroyed;
c
c        z contains the first components of the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.  if an error exit is
c          made, z contains the eigenvectors associated with the stored
c          eigenvalues;
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     ------------------------------------------------------------------
c
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8 d(n), e(n), z(n), b, c, f, g, p, r, s, machep
      real*8 dsqrt, dabs, dsign, d1mach
c
      machep=d1mach(4)
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
c     :::::::::: look for small sub-diagonal element ::::::::::
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if (dabs(e(m)) .le. machep * (dabs(d(m)) + dabs(d(m+1))))
     x         go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = dsqrt(g*g+1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c
c     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (dabs(f) .lt. dabs(g)) go to 150
            c = g / f
            r = dsqrt(c*c+1.0d0)
            e(i+1) = f * r
            s = 1.0d0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = dsqrt(s*s+1.0d0)
            e(i+1) = g * r
            c = 1.0d0 / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     :::::::::: form first component of vector ::::::::::
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
  200       z(i) = c * z(i) - s * f
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
c
c     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         p = z(i)
         z(i) = z(k)
         z(k) = p
  300 continue
c
      go to 1001
c     :::::::::: set error -- no convergence to an
c                eigenvalue after 30 iterations ::::::::::
 1000 ierr = l
 1001 return
c     :::::::::: last card of gausq2 ::::::::::
      end
