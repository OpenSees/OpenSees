! Modify material properties
! Modify cp_ref for PML      
      
!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE PML_3D(MMATRX,CMATRX,KMATRX,GMATRX,HMATRX,
     1               NDOFEL,PROPS,COORDS,MCRD,NNODE,LFLAGS)
    !
    ! INCLUDE 'ABA_PARAM.INC'
    ! 
    !
      REAL *8 PROPS,COORDS,MMATRX, CMATRX, KMATRX, GMATRX, HMATRX
      INTEGER *4 NDOFEL,MCRD,NNODE,LFLAGS

      DIMENSION MMATRX(NDOFEL,NDOFEL),
     1          CMATRX(NDOFEL,NDOFEL),
     2          KMATRX(NDOFEL,NDOFEL),
     3          GMATRX(NDOFEL,NDOFEL),
     4          HMATRX(NDOFEL,NDOFEL),
     5          PROPS(*),COORDS(MCRD,NNODE),LFLAGS(*)

    ! Local Variables

      real*8 ZERO,ONE,HALF       
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0 )

      real*8 coef_alpha, coef_beta
      PARAMETER (coef_alpha = 1.d0/12.d0)
      PARAMETER (coef_beta = 1.d0/48.d0)

      integer      :: i,j,n_points,kint
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt spatial coords
    
      double precision  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  D(6,6)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(6,60)                           ! strain = B*(dof_total)
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  E, xnu                            ! Material properties
      double precision  ::  Phi(3,60)                         ! Matrix consists of shape functions [N1 N2 ... N20]
      
      double precision  ::  rho                               ! Density of the material

      double precision  ::  PML_L,afp                         ! Parameters of PML
      double precision  ::  alpha_0, beta_0                   ! Parameters of PML
      double precision  ::  x1_0, x2_0, x3_0                  ! Parameters of PML
      double precision  ::  n1, n2, n3                        ! Parameters of PML
      integer           ::  EleType_pos                       ! Parameters of PML

      double precision  ::  x1,x2,x3,PML_alpha_beta(2,3)      ! Coordinates x, y and z, alphas and betas    
      double precision  ::  coef_a, coef_b, coef_c, coef_d    ! Coefficients of the PML
      double precision  ::  coef_Le(3,3), coef_Lp(3,3), coef_Lw(3,3)    ! Coefficients of the PML

      double precision  ::  K_RD(24,24), Kxx(8,8), Kxy(8,8), Kxz(8,8)   ! Stiffness matrices
      double precision  ::  Kyx(8,8), Kyy(8,8), Kyz(8,8)                ! Stiffness matrices
      double precision  ::  Kzx(8,8), Kzy(8,8), Kzz(8,8)                ! Stiffness matrices
      double precision  ::  Phi_x(8), Phi_y(8), Phi_z(8)                ! Derivative of shape functions for displacement fields [N1,i N2,i ... N8,i]^T
      double precision  ::  lambda, mu                                  ! Lame' constants
      double precision  ::  M_RD(24,24)                                 ! Mass matrix for regular domain                              
      double precision  ::  C_RD(24,24)                                 ! Damping matrix for regular domain                              
      double precision  ::  M_a(24,24),M_b(24,24),M_c(24,24),M_d(24,24) ! Mass matrices
      double precision  ::  N_a(48,48),N_b(48,48),N_c(48,48),N_d(48,48) ! N matrices for PML
      double precision  ::  A_eu(24,48), A_wu(24,48), A_pu(24,48)       ! A_iu matrices for PML, where i = e, w, p
      double precision  ::  A_el(24,48), A_wl(24,48), A_pl(24,48)       ! A_il matrices for PML, where i = e, w, p
      double precision  ::  K_PML(72,72), M_PML(72,72), C_PML(72,72)    ! Stiffnee, mass and damping matrices for PML
      double precision  ::  G_PML(72,72)                                ! G matrices for PML
      double precision  ::  H_PML(72,72)                                ! H matrices for PML

      double precision  ::  d_bar_n(72), d_bar_n1(72)                   ! d_bar at previous and current step
      double precision  ::  d_bar_dot_n(72), d_bar_dot_n1(72)           ! d_bar_dot at previous and current step
      
      double precision  ::  U_n(72), V_n(72), A_n(72)                   ! U, V, A at previous step

      integer           ::  NodeDOF
      
      double precision  ::  Damp_alpha, Damp_beta                       ! Damping coefficients for Rayleigh Damping
      


      NodeDOF = NDOFEL/NNODE
      



      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 27               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex
          ! Local Variables
      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)


      MMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      KMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      CMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      GMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      HMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      

      E = PROPS(1)
      xnu = PROPS(2)
      rho = PROPS(3)
      EleType_pos = floor(PROPS(4))
      PML_L = PROPS(5)
      x1_0 = PROPS(6)
      x2_0 = PROPS(7)
      x3_0 = PROPS(8)
      n1 = PROPS(9)
      n2 = PROPS(10)
      n3 = PROPS(11)
      alpha_0 = PROPS(12)
      beta_0 = PROPS(13)
      Damp_alpha = PROPS(14)
      Damp_beta = PROPS(15)
      afp = PROPS(16)


      IF (afp .LT. 3.d0) THEN
        n_points = 27
      else
        n_points = 64
      endif
            
      lambda = xnu*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
      mu = 0.5D0*E/(1.d0+xnu)


      M_RD = 0.d0
      M_a = 0.d0
      M_b = 0.d0
      M_c = 0.d0
      M_d = 0.d0

      N_a = 0.d0
      N_b = 0.d0
      N_c = 0.d0
      N_d = 0.d0

      A_eu = 0.d0
      A_wu = 0.d0
      A_pu = 0.d0

      A_el = 0.d0
      A_wl = 0.d0
      A_pl = 0.d0

      K_RD = 0.d0
      C_RD = 0.d0

      K_PML = 0.d0
      C_PML = 0.d0
      M_PML = 0.d0
      G_PML = 0.d0
      H_PML = 0.d0

      d_bar_n = 0.d0
      d_bar_n1 = 0.d0
      
      d_bar_dot_n = 0.d0
      d_bar_dot_n1 = 0.d0

      U_n = 0.d0
      V_n = 0.d0
      A_n = 0.d0
      
    !     --  Loop over integration points
      do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
        B = 0.d0
        B(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
        B(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
        B(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
        B(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
        B(4,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
        B(5,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
        B(5,3:3*NNODE:3)   = dNdx(1:NNODE,1)
        B(6,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
        B(6,3:3*NNODE:3)   = dNdx(1:NNODE,2)

      !   strain = matmul(B(1:6,1:3*NNODE),U(1:3*NNODE))

        stress = matmul(D,strain)


         !  ALPHA = PARAMS(1)
         !  BETA  = PARAMS(2)
         !  GAMMA = PARAMS(3)

         !  DADU = ONE/(BETA*DTIME**2.d0)
         !  DVDU = GAMMA/(BETA*DTIME)

          Phi = 0.d0
          Phi(1,1:3*NNODE-2:3) = N(1:NNODE)
          Phi(2,2:3*NNODE-1:3) = N(1:NNODE)
          Phi(3,3:3*NNODE:3)   = N(1:NNODE)

          Phi_x = 0.d0
          Phi_y = 0.d0
          Phi_z = 0.d0
          Phi_x(1:NNODE) = dNdx(1:NNODE,1)
          Phi_y(1:NNODE) = dNdx(1:NNODE,2)
          Phi_z(1:NNODE) = dNdx(1:NNODE,3)


          x1 = 0.d0
          x2 = 0.d0
          x3 = 0.d0

          do i = 1,NNODE
            x1 = x1 + N(i)*coords(1,i)
            x2 = x2 + N(i)*coords(2,i)
            x3 = x3 + N(i)*coords(3,i)
          end do



         !  call PML_alpha_beta_function(PROPS,x1,x2,x3,PML_alpha_beta)
         !  PML_alpha_beta(1,1) = 1.d0 + alpha_0*((x1 -x1_0) * n1 /PML_L)**afp 
         !  PML_alpha_beta(1,2) = 1.d0 + alpha_0*((x2 -x2_0) * n2 /PML_L)**afp
         !  PML_alpha_beta(1,3) = 1.d0 + alpha_0*((x3 -x3_0) * n3 /PML_L)**afp
          PML_alpha_beta(1,1) = alpha_0*((x1 -x1_0) * n1 /PML_L)**afp
          PML_alpha_beta(1,2) = alpha_0*((x2 -x2_0) * n2 /PML_L)**afp
          PML_alpha_beta(1,3) = alpha_0*((x3 -x3_0) * n3 /PML_L)**afp
          PML_alpha_beta(1,1) = 1.d0 + PML_alpha_beta(1,1)
          PML_alpha_beta(1,2) = 1.d0 + PML_alpha_beta(1,2)
          PML_alpha_beta(1,3) = 1.d0 + PML_alpha_beta(1,3)
    
          PML_alpha_beta(2,1) = beta_0*((x1 -x1_0) * n1 /PML_L)**afp 
          PML_alpha_beta(2,2) = beta_0*((x2 -x2_0) * n2 /PML_L)**afp
          PML_alpha_beta(2,3) = beta_0*((x3 -x3_0) * n3 /PML_L)**afp

          coef_a = PML_alpha_beta(1,1)*PML_alpha_beta(1,2)
     1            *PML_alpha_beta(1,3)
         


          coef_b = PML_alpha_beta(1,1)*PML_alpha_beta(1,2)
     1            *PML_alpha_beta(2,3) + 
     2             PML_alpha_beta(1,1)*PML_alpha_beta(1,3)
     3            *PML_alpha_beta(2,2) + 
     4             PML_alpha_beta(1,2)*PML_alpha_beta(1,3)
     5            *PML_alpha_beta(2,1) 

          coef_c = PML_alpha_beta(1,1)*PML_alpha_beta(2,2)
     1            *PML_alpha_beta(2,3) + 
     2             PML_alpha_beta(1,2)*PML_alpha_beta(2,3)
     3            *PML_alpha_beta(2,1) + 
     4             PML_alpha_beta(1,3)*PML_alpha_beta(2,2)
     5            *PML_alpha_beta(2,1)

          coef_d = PML_alpha_beta(2,1)*PML_alpha_beta(2,2)
     1            *PML_alpha_beta(2,3) 

          coef_Le = 0.d0
          coef_Lp = 0.d0
          coef_Lw = 0.d0

          do i = 1,3
            do j = 1,3
              coef_Le(i,j) = PML_alpha_beta(1,i)*PML_alpha_beta(1,j)
              coef_Lp(i,j) = PML_alpha_beta(1,i)*PML_alpha_beta(2,j) +
     1                       PML_alpha_beta(2,i)*PML_alpha_beta(1,j) 
              coef_Lw(i,j) = PML_alpha_beta(2,i)*PML_alpha_beta(2,j)  
            end do
          end do
          
          
          do i = 1,NNODE
            do j = 1,NNODE
              Kxx(i,j) = (lambda + 2.d0*mu) * Phi_x(i)*Phi_x(j) + 
     1                    mu*(Phi_y(i)*Phi_y(j)+Phi_z(i)*Phi_z(j))     
              Kyy(i,j) = (lambda + 2.d0*mu) * Phi_y(i)*Phi_y(j) + 
     1                    mu*(Phi_x(i)*Phi_x(j)+Phi_z(i)*Phi_z(j)) 
              Kzz(i,j) = (lambda + 2.d0*mu) * Phi_z(i)*Phi_z(j) + 
     1                    mu*(Phi_x(i)*Phi_x(j)+Phi_y(i)*Phi_y(j))  

              Kxy(i,j) = lambda * Phi_x(i)*Phi_y(j) + 
     1                   mu * Phi_y(i)*Phi_x(j)

              Kxz(i,j) = lambda * Phi_x(i)*Phi_z(j) + 
     1                   mu * Phi_z(i)*Phi_x(j)

              Kyz(i,j) = lambda * Phi_y(i)*Phi_z(j) + 
     1                   mu * Phi_z(i)*Phi_y(j)

              M_RD(i,j) = M_RD(i,j) + rho*N(i)*N(j)*w(kint)*determinant 
              M_a(i,j) = M_a(i,j) + 
     1                         coef_a*rho*N(i)*N(j)*w(kint)*determinant



     
              M_b(i,j) = M_b(i,j) + 
     1                         coef_b*rho*N(i)*N(j)*w(kint)*determinant
              M_c(i,j) = M_c(i,j) + 
     1                         coef_c*rho*N(i)*N(j)*w(kint)*determinant
              M_d(i,j) = M_d(i,j) + 
     1                         coef_d*rho*N(i)*N(j)*w(kint)*determinant

              A_eu(i,j) = A_eu(i,j) + 
     1                   Phi_x(i)*N(j)*coef_Le(2,3)*w(kint)*determinant
              A_eu(i,j+NNODE*3) = A_eu(i,j+NNODE*3) + 
     1                   Phi_y(i)*N(j)*coef_Le(1,3)*w(kint)*determinant 
              A_eu(i,j+NNODE*4) = A_eu(i,j+NNODE*4) + 
     1                   Phi_z(i)*N(j)*coef_Le(1,2)*w(kint)*determinant 

              A_eu(i+NNODE,j+NNODE) = A_eu(i+NNODE,j+NNODE) + 
     1                   Phi_y(i)*N(j)*coef_Le(1,3)*w(kint)*determinant 
              A_eu(i+NNODE,j+NNODE*3) = A_eu(i+NNODE,j+NNODE*3) + 
     1                   Phi_x(i)*N(j)*coef_Le(2,3)*w(kint)*determinant 
              A_eu(i+NNODE,j+NNODE*5) = A_eu(i+NNODE,j+NNODE*5) + 
     1                   Phi_z(i)*N(j)*coef_Le(1,2)*w(kint)*determinant 

              A_eu(i+NNODE*2,j+NNODE*2) = A_eu(i+NNODE*2,j+NNODE*2) + 
     1                   Phi_z(i)*N(j)*coef_Le(1,2)*w(kint)*determinant 
              A_eu(i+NNODE*2,j+NNODE*4) = A_eu(i+NNODE*2,j+NNODE*4) + 
     1                   Phi_x(i)*N(j)*coef_Le(2,3)*w(kint)*determinant
              A_eu(i+NNODE*2,j+NNODE*5) = A_eu(i+NNODE*2,j+NNODE*5) + 
     1                   Phi_y(i)*N(j)*coef_Le(1,3)*w(kint)*determinant

              A_wu(i,j) = A_wu(i,j) + 
     1                   Phi_x(i)*N(j)*coef_Lw(2,3)*w(kint)*determinant
              A_wu(i,j+NNODE*3) = A_wu(i,j+NNODE*3) + 
     1                   Phi_y(i)*N(j)*coef_Lw(1,3)*w(kint)*determinant 
              A_wu(i,j+NNODE*4) = A_wu(i,j+NNODE*4) + 
     1                   Phi_z(i)*N(j)*coef_Lw(1,2)*w(kint)*determinant 

              A_wu(i+NNODE,j+NNODE) = A_wu(i+NNODE,j+NNODE) + 
     1                   Phi_y(i)*N(j)*coef_Lw(1,3)*w(kint)*determinant 
              A_wu(i+NNODE,j+NNODE*3) = A_wu(i+NNODE,j+NNODE*3) + 
     1                   Phi_x(i)*N(j)*coef_Lw(2,3)*w(kint)*determinant 
              A_wu(i+NNODE,j+NNODE*5) = A_wu(i+NNODE,j+NNODE*5) + 
     1                   Phi_z(i)*N(j)*coef_Lw(1,2)*w(kint)*determinant 

              A_wu(i+NNODE*2,j+NNODE*2) = A_wu(i+NNODE*2,j+NNODE*2) + 
     1                   Phi_z(i)*N(j)*coef_Lw(1,2)*w(kint)*determinant 
              A_wu(i+NNODE*2,j+NNODE*4) = A_wu(i+NNODE*2,j+NNODE*4) + 
     1                   Phi_x(i)*N(j)*coef_Lw(2,3)*w(kint)*determinant
              A_wu(i+NNODE*2,j+NNODE*5) = A_wu(i+NNODE*2,j+NNODE*5) + 
     1                   Phi_y(i)*N(j)*coef_Lw(1,3)*w(kint)*determinant

              A_pu(i,j) = A_pu(i,j) + 
     1                   Phi_x(i)*N(j)*coef_Lp(2,3)*w(kint)*determinant
              A_pu(i,j+NNODE*3) = A_pu(i,j+NNODE*3) + 
     1                   Phi_y(i)*N(j)*coef_Lp(1,3)*w(kint)*determinant 
              A_pu(i,j+NNODE*4) = A_pu(i,j+NNODE*4) + 
     1                   Phi_z(i)*N(j)*coef_Lp(1,2)*w(kint)*determinant 

              A_pu(i+NNODE,j+NNODE) = A_pu(i+NNODE,j+NNODE) + 
     1                   Phi_y(i)*N(j)*coef_Lp(1,3)*w(kint)*determinant 
              A_pu(i+NNODE,j+NNODE*3) = A_pu(i+NNODE,j+NNODE*3) + 
     1                   Phi_x(i)*N(j)*coef_Lp(2,3)*w(kint)*determinant 
              A_pu(i+NNODE,j+NNODE*5) = A_pu(i+NNODE,j+NNODE*5) + 
     1                   Phi_z(i)*N(j)*coef_Lp(1,2)*w(kint)*determinant 

              A_pu(i+NNODE*2,j+NNODE*2) = A_pu(i+NNODE*2,j+NNODE*2) + 
     1                   Phi_z(i)*N(j)*coef_Lp(1,2)*w(kint)*determinant 
              A_pu(i+NNODE*2,j+NNODE*4) = A_pu(i+NNODE*2,j+NNODE*4) + 
     1                   Phi_x(i)*N(j)*coef_Lp(2,3)*w(kint)*determinant
              A_pu(i+NNODE*2,j+NNODE*5) = A_pu(i+NNODE*2,j+NNODE*5) + 
     1                   Phi_y(i)*N(j)*coef_Lp(1,3)*w(kint)*determinant

                 end do
          end do


          Kyx(1:NNODE,1:NNODE) = transpose(Kxy(1:NNODE,1:NNODE)) 
          Kzx(1:NNODE,1:NNODE) = transpose(Kxz(1:NNODE,1:NNODE)) 
          Kzy(1:NNODE,1:NNODE) = transpose(Kyz(1:NNODE,1:NNODE)) 

          K_RD(1:NNODE,1:NNODE) = K_RD(1:NNODE,1:NNODE) + 
     1                          Kxx(1:NNODE,1:NNODE)*w(kint)*determinant 
          K_RD(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1            K_RD(NNODE+1:2*NNODE,NNODE+1:2*NNODE) + 
     2                          Kyy(1:NNODE,1:NNODE)*w(kint)*determinant 
          K_RD(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) = 
     1        K_RD(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) + 
     2                          Kzz(1:NNODE,1:NNODE)*w(kint)*determinant

          K_RD(1:NNODE,NNODE+1:2*NNODE) =
     1                    K_RD(1:NNODE,NNODE+1:2*NNODE) +  
     2                          Kxy(1:NNODE,1:NNODE)*w(kint)*determinant  
          K_RD(1:NNODE,2*NNODE+1:3*NNODE) =
     1                  K_RD(1:NNODE,2*NNODE+1:3*NNODE) +  
     2                          Kxz(1:NNODE,1:NNODE)*w(kint)*determinant     
          K_RD(NNODE+1:2*NNODE,2*NNODE+1:3*NNODE) =
     1          K_RD(NNODE+1:2*NNODE,2*NNODE+1:3*NNODE) +  
     2                          Kyz(1:NNODE,1:NNODE)*w(kint)*determinant

          K_RD(NNODE+1:2*NNODE,1:NNODE) =
     1                    K_RD(NNODE+1:2*NNODE,1:NNODE) +  
     2                          Kyx(1:NNODE,1:NNODE)*w(kint)*determinant  
          K_RD(2*NNODE+1:3*NNODE,1:NNODE) =
     1                  K_RD(2*NNODE+1:3*NNODE,1:NNODE) +  
     2                          Kzx(1:NNODE,1:NNODE)*w(kint)*determinant     
          K_RD(2*NNODE+1:3*NNODE,NNODE+1:2*NNODE) =
     1          K_RD(2*NNODE+1:3*NNODE,NNODE+1:2*NNODE) +  
     2                          Kzy(1:NNODE,1:NNODE)*w(kint)*determinant  


      end do

        

        M_RD(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_RD(1:NNODE,1:NNODE)
        M_RD(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) 
     1                                        = M_RD(1:NNODE,1:NNODE)
        
        C_RD(1:3*NNODE,1:3*NNODE) = Damp_alpha*M_RD(1:3*NNODE,1:3*NNODE)
     1                            + Damp_beta*K_RD(1:3*NNODE,1:3*NNODE)

        M_a(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_a(1:NNODE,1:NNODE)
        M_a(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) 
     1                                       = M_a(1:NNODE,1:NNODE)

        M_b(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_b(1:NNODE,1:NNODE)
        M_b(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) 
     1                                       = M_b(1:NNODE,1:NNODE)

        M_c(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_c(1:NNODE,1:NNODE)
        M_c(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) 
     1                                       = M_c(1:NNODE,1:NNODE)

        M_d(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_d(1:NNODE,1:NNODE)
        M_d(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) 
     1                                       = M_d(1:NNODE,1:NNODE)  


        N_a(1:NNODE,1:NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_a(1:NNODE,NNODE+1:2*NNODE) = 
     1  -M_a(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_a(1:NNODE,2*NNODE+1:3*NNODE) = 
     1  -M_a(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_a(NNODE+1:2*NNODE,1:NNODE) = 
     1  -M_a(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_a(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_a(NNODE+1:2*NNODE,2*NNODE+1:3*NNODE) = 
     1  -M_a(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_a(2*NNODE+1:3*NNODE,1:NNODE) = 
     1  -M_a(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_a(2*NNODE+1:3*NNODE,NNODE+1:2*NNODE) = 
     1  -M_a(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_a(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_a(3*NNODE+1:4*NNODE,3*NNODE+1:4*NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho/mu
        N_a(4*NNODE+1:5*NNODE,4*NNODE+1:5*NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho/mu
        N_a(5*NNODE+1:6*NNODE,5*NNODE+1:6*NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho/mu

        N_b(1:NNODE,1:NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_b(1:NNODE,NNODE+1:2*NNODE) = 
     1  -M_b(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_b(1:NNODE,2*NNODE+1:3*NNODE) = 
     1  -M_b(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_b(NNODE+1:2*NNODE,1:NNODE) = 
     1  -M_b(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_b(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_b(NNODE+1:2*NNODE,2*NNODE+1:3*NNODE) = 
     1  -M_b(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_b(2*NNODE+1:3*NNODE,1:NNODE) = 
     1  -M_b(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_b(2*NNODE+1:3*NNODE,NNODE+1:2*NNODE) = 
     1  -M_b(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_b(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_b(3*NNODE+1:4*NNODE,3*NNODE+1:4*NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho/mu
        N_b(4*NNODE+1:5*NNODE,4*NNODE+1:5*NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho/mu
        N_b(5*NNODE+1:6*NNODE,5*NNODE+1:6*NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho/mu

        N_c(1:NNODE,1:NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_c(1:NNODE,NNODE+1:2*NNODE) = 
     1  -M_c(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_c(1:NNODE,2*NNODE+1:3*NNODE) = 
     1  -M_c(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_c(NNODE+1:2*NNODE,1:NNODE) = 
     1  -M_c(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_c(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_c(NNODE+1:2*NNODE,2*NNODE+1:3*NNODE) = 
     1  -M_c(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_c(2*NNODE+1:3*NNODE,1:NNODE) = 
     1  -M_c(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_c(2*NNODE+1:3*NNODE,NNODE+1:2*NNODE) = 
     1  -M_c(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_c(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_c(3*NNODE+1:4*NNODE,3*NNODE+1:4*NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho/mu
        N_c(4*NNODE+1:5*NNODE,4*NNODE+1:5*NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho/mu
        N_c(5*NNODE+1:6*NNODE,5*NNODE+1:6*NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho/mu

        N_d(1:NNODE,1:NNODE) = 
     1   M_d(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_d(1:NNODE,NNODE+1:2*NNODE) = 
     1  -M_d(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_d(1:NNODE,2*NNODE+1:3*NNODE) = 
     1  -M_d(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_d(NNODE+1:2*NNODE,1:NNODE) = 
     1  -M_d(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_d(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1   M_d(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_d(NNODE+1:2*NNODE,2*NNODE+1:3*NNODE) = 
     1  -M_d(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_d(2*NNODE+1:3*NNODE,1:NNODE) = 
     1  -M_d(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_d(2*NNODE+1:3*NNODE,NNODE+1:2*NNODE) = 
     1  -M_d(1:NNODE,1:NNODE)/rho*(lambda)/mu/2.d0/(3.d0*lambda+2.d0*mu)
        N_d(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) = 
     1   M_d(1:NNODE,1:NNODE)/rho*(lambda+mu)/mu/(3.d0*lambda+2.d0*mu)
        N_d(3*NNODE+1:4*NNODE,3*NNODE+1:4*NNODE) = 
     1   M_d(1:NNODE,1:NNODE)/rho/mu
        N_d(4*NNODE+1:5*NNODE,4*NNODE+1:5*NNODE) = 
     1   M_d(1:NNODE,1:NNODE)/rho/mu
        N_d(5*NNODE+1:6*NNODE,5*NNODE+1:6*NNODE) = 
     1   M_d(1:NNODE,1:NNODE)/rho/mu          


      IF (EleType_pos .EQ. 1 .OR. LFLAGS(1).EQ.1 .OR. 
     &       LFLAGS(1).EQ.2) THEN

          M_PML(1:NNODE*3,1:NNODE*3) = M_RD(1:NNODE*3,1:NNODE*3)

          K_PML(1:NNODE*3,1:NNODE*3) = K_RD(1:NNODE*3,1:NNODE*3)
          write(*,*) 'this 2 if statement is working'



      else
          IF (EleType_pos .EQ. 1) THEN

          M_PML(1:NNODE*3,1:NNODE*3) = M_RD(1:NNODE*3,1:NNODE*3) + 
     1                             M_a(1:NNODE*3,1:NNODE*3)
          M_PML(NNODE*3+1:NNODE*9,NNODE*3+1:NNODE*9) = 
     1                                      -N_a(1:NNODE*6,1:NNODE*6)

          C_PML(1:NNODE*3,1:NNODE*3) = C_RD(1:NNODE*3,1:NNODE*3) + 
     1                             M_b(1:NNODE*3,1:NNODE*3) 
          C_PML(1:NNODE*3,NNODE*3+1:NNODE*9) = A_eu(1:NNODE*3,1:NNODE*6)
          C_PML(NNODE*3+1:NNODE*9,1:NNODE*3) = 
     1                         transpose(A_eu(1:NNODE*3,1:NNODE*6))
          C_PML(NNODE*3+1:NNODE*9,NNODE*3+1:NNODE*9) = 
     1                                      -N_b(1:NNODE*6,1:NNODE*6)

          K_PML(1:NNODE*3,1:NNODE*3) = K_RD(1:NNODE*3,1:NNODE*3) + 
     1                             M_c(1:NNODE*3,1:NNODE*3)
          K_PML(1:NNODE*3,NNODE*3+1:NNODE*9) = A_pu(1:NNODE*3,1:NNODE*6)
          K_PML(NNODE*3+1:NNODE*9,1:NNODE*3) = 
     1                         transpose(A_pu(1:NNODE*3,1:NNODE*6))
          K_PML(NNODE*3+1:NNODE*9,NNODE*3+1:NNODE*9) = 
     1                                      -N_c(1:NNODE*6,1:NNODE*6)

          G_PML(1:NNODE*3,1:NNODE*3) = M_d(1:NNODE*3,1:NNODE*3)
          G_PML(1:NNODE*3,NNODE*3+1:NNODE*9) = A_wu(1:NNODE*3,1:NNODE*6)
          G_PML(NNODE*3+1:NNODE*9,1:NNODE*3) = 
     1                         transpose(A_wu(1:NNODE*3,1:NNODE*6))
          G_PML(NNODE*3+1:NNODE*9,NNODE*3+1:NNODE*9) = 
     1                                      -N_d(1:NNODE*6,1:NNODE*6)
         !  write(*,*) 'this 1 if statement is working'

        else
         ! write(*,*) 'this if statement is working'
          M_PML(1:NNODE*3,1:NNODE*3) = 
     1                            M_a(1:NNODE*3,1:NNODE*3)
          M_PML(1:NNODE*3,NNODE*3+1:NNODE*9) = 
     1                            A_eu(1:NNODE*3,1:NNODE*6)*Damp_beta 
          M_PML(NNODE*3+1:NNODE*9,NNODE*3+1:NNODE*9) = 
     1                                      -N_a(1:NNODE*6,1:NNODE*6)

          C_PML(1:NNODE*3,1:NNODE*3) = 
     1    M_b(1:NNODE*3,1:NNODE*3) + 
     2    M_a(1:NNODE*3,1:NNODE*3)*Damp_alpha
          
          C_PML(1:NNODE*3,NNODE*3+1:NNODE*9) = 
     1                            A_pu(1:NNODE*3,1:NNODE*6)*Damp_beta +
     2                            A_eu(1:NNODE*3,1:NNODE*6)    
          C_PML(NNODE*3+1:NNODE*9,1:NNODE*3) = 
     1                         transpose(A_eu(1:NNODE*3,1:NNODE*6))
          C_PML(NNODE*3+1:NNODE*9,NNODE*3+1:NNODE*9) = 
     1                                      -N_b(1:NNODE*6,1:NNODE*6)

          K_PML(1:NNODE*3,1:NNODE*3) = 
     1    M_c(1:NNODE*3,1:NNODE*3) + M_b(1:NNODE*3,1:NNODE*3)*Damp_alpha
          K_PML(1:NNODE*3,NNODE*3+1:NNODE*9) = 
     1    A_pu(1:NNODE*3,1:NNODE*6)+A_wu(1:NNODE*3,1:NNODE*6)*Damp_beta
          K_PML(NNODE*3+1:NNODE*9,1:NNODE*3) = 
     1                         transpose(A_pu(1:NNODE*3,1:NNODE*6))
          K_PML(NNODE*3+1:NNODE*9,NNODE*3+1:NNODE*9) = 
     1                                      -N_c(1:NNODE*6,1:NNODE*6)

          G_PML(1:NNODE*3,1:NNODE*3) = 
     1    M_d(1:NNODE*3,1:NNODE*3) + M_c(1:NNODE*3,1:NNODE*3)*Damp_alpha
          G_PML(1:NNODE*3,NNODE*3+1:NNODE*9) = A_wu(1:NNODE*3,1:NNODE*6)
          G_PML(NNODE*3+1:NNODE*9,1:NNODE*3) = 
     1                         transpose(A_wu(1:NNODE*3,1:NNODE*6))
          G_PML(NNODE*3+1:NNODE*9,NNODE*3+1:NNODE*9) = 
     1                                      -N_d(1:NNODE*6,1:NNODE*6)
          
          H_PML(1:NNODE*3,1:NNODE*3) = 
     1     M_d(1:NNODE*3,1:NNODE*3)*Damp_alpha


        endif 


      endif 


        do i = 1,8
          do j = 1,8
            MMATRX((i-1)*9+1:i*9,(j-1)*9+1:j*9) = 
     1                        M_PML(i:NDOFEL-8+i:8,j:NDOFEL-8+j:8)

            CMATRX((i-1)*9+1:i*9,(j-1)*9+1:j*9) = 
     1                        C_PML(i:NDOFEL-8+i:8,j:NDOFEL-8+j:8)

            KMATRX((i-1)*9+1:i*9,(j-1)*9+1:j*9) = 
     1                        K_PML(i:NDOFEL-8+i:8,j:NDOFEL-8+j:8)

            GMATRX((i-1)*9+1:i*9,(j-1)*9+1:j*9) = 
     1                        G_PML(i:NDOFEL-8+i:8,j:NDOFEL-8+j:8)
            HMATRX((i-1)*9+1:i*9,(j-1)*9+1:j*9) = 
     1                        H_PML(i:NDOFEL-8+i:8,j:NDOFEL-8+j:8)
          end do
        end do

      return

      END SUBROUTINE PML_3D

   



      subroutine abq_UEL_3D_integrationpoints(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(3,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)

    !         Defines integration points and weights for 3D continuum elements

      if (n_nodes  == 4.or.n_nodes ==10 ) then   ! Tetrahedral elements
        if (n_points == 1) then
            xi(1,1) = 0.25D0
            xi(2,1) = 0.25D0
            xi(3,1) = 0.25D0
            w(1) = 1.D0/6.D0
        else if (n_points == 4) then
            xi(1,1) = 0.58541020
            xi(2,1) = 0.13819660
            xi(3,1) = xi(2,1)
            xi(1,2) = xi(2,1)
            xi(2,2) = xi(1,1)
            xi(3,2) = xi(2,1)
            xi(1,3) = xi(2,1)
            xi(2,3) = xi(2,1)
            xi(3,3) = xi(1,1)
            xi(1,4) = xi(2,1)
            xi(2,4) = xi(2,1)
            xi(3,4) = xi(2,1)
            w(1:4) = 1.D0/24.D0
        else if (n_points == 5) then
            xi(1,1) = 0.25d0
            xi(2,1) = 0.25d0
            xi(3,1) = 0.25d0
            xi(1,2) = 0.5d0
            xi(2,2) = 1.d0/6.d0
            xi(3,2) = 1.d0/6.d0
            xi(1,3) = 1.d0/6.d0
            xi(2,3) = 0.5d0
            xi(3,3) = 1.d0/6.d0
            xi(1,4) = 1.d0/6.d0
            xi(2,4) = 1.d0/6.d0
            xi(3,4) = 0.5d0
            xi(1,5) = 1.d0/6.d0
            xi(2,5) = 1.d0/6.d0
            xi(3,5) = 1.d0/6.d0
            w(1) = -4.d0/30.d0
            w(2:5) = 3.d0/40.d0
        else
            write(6,*) 'Incorrect # of int pts for tetrahedral element '
            write(6, *) ' called with ',n_points
            stop
        endif
      else if ( n_nodes == 8 .or. n_nodes == 20 ) then   ! 8 or 20 noded hexahedral elements
        if (n_points == 1) then
            xi(1,1) = 0.D0
            xi(2,1) = 0.D0
            xi(3,1) = 0.D0
            w(1) = 8.D0
        else if (n_points == 8) then
            x1D(1) = -0.5773502692
            x1D(2) =  0.5773502692
            do k = 1,2
                do j = 1,2
                    do i = 1,2
                        n = 4*(k-1) + 2*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                    end do
                end do
            end do
            w(1:8) = 1.D0
        else if (n_points == 27) then
            x1D(1) = -0.7745966692
            x1D(2) = 0.
            x1D(3) = 0.7745966692
            w1D(1) = 0.5555555555D0
            w1D(2) = 0.888888888D0
            w1D(3) = 0.55555555555D0
            do k = 1,3
                do j = 1,3
                    do i = 1,3
                        n = 9*(k-1) + 3*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                        w(n) = w1D(i)*w1D(j)*w1D(k)
                    end do
                end do
            end do
        else if (n_points == 64) then
            x1D(1) = .8611363115940526D+00
            x1D(2) = .3399810435848563D+00
            x1D(3) = -.3399810435848563D+00
            x1D(4) = -.8611363115940526D+00
            w1D(1) = .3478548451374538D+00
            w1D(2) = .6521451548625461D+00
            w1D(3) = .6521451548625461D+00
            w1D(4) = .3478548451374538D+00
            do k = 1,4
                do j = 1,4
                    do i = 1,4
                        n = 16*(k-1) + 4*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                        w(n) = w1D(i)*w1D(j)*w1D(k)
                    end do
                end do
            end do
        endif
      endif

      return

      end subroutine abq_UEL_3D_integrationpoints


      subroutine abq_UEL_3D_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(3)
      double precision, intent(out) :: f(20)
      double precision, intent(out) :: df(20,3)
      double precision xi4

!   Defines shape functions for 3D continuum elements

      if (n_nodes == 4) then
        f(1) = xi(1)
        f(2) = xi(2)
        f(3) = xi(3)
        f(4) = 1.-xi(1)-xi(2)-xi(3)
        df(1,1) = 1.
        df(2,2) = 1.
        df(3,3) = 1.
        df(4,1) = -1.
        df(4,2) = -1.
        df(4,3) = -1.
      else if (n_nodes == 10) then
        xi4 = 1.D0-xi(1)-xi(2)-xi(3)
        f(1) = (2.*xi(1)-1.)*xi(1)
        f(2) = (2.*xi(2)-1.)*xi(2)
        f(3) = (2.*xi(3)-1.)*xi(3)
        f(4) = (2.*xi4-1.)*xi4
        f(5) = 4.*xi(1)*xi(2)
        f(6) = 4.*xi(2)*xi(3)
        f(7) = 4.*xi(3)*xi(1)
        f(8) = 4.*xi(1)*xi4
        f(9) = 4.*xi(2)*xi4
        f(10) = 4.*xi(3)*xi4
        df(1,1) = (4.*xi(1)-1.)
        df(2,2) = (4.*xi(2)-1.)
        df(3,3) = (4.*xi(3)-1.)
        df(4,1) = -(4.*xi4-1.)
        df(4,2) = -(4.*xi4-1.)
        df(4,3) = -(4.*xi4-1.)
        df(5,1) = 4.*xi(2)
        df(5,2) = 4.*xi(1)
        df(6,2) = 4.*xi(3)
        df(6,3) = 4.*xi(2)
        df(7,1) = 4.*xi(3)
        df(7,3) = 4.*xi(1)
        df(8,1) = 4.*(xi4-xi(1))
        df(8,2) = -4.*xi(1)
        df(8,3) = -4.*xi(1)
        df(9,1) = -4.*xi(2)
        df(9,2) = 4.*(xi4-xi(2))
        df(9,3) = -4.*xi(2)
        df(10,1) = -4.*xi(3)*xi4
        df(10,2) = -4.*xi(3)
        df(10,3) = 4.*(xi4-xi(3))
      else if (n_nodes == 8) then
        f(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.
        f(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.
        f(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.
        f(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.
        f(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.
        f(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.
        f(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.
        f(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.
        df(1,1) = -(1.-xi(2))*(1.-xi(3))/8.
        df(1,2) = -(1.-xi(1))*(1.-xi(3))/8.
        df(1,3) = -(1.-xi(1))*(1.-xi(2))/8.
        df(2,1) = (1.-xi(2))*(1.-xi(3))/8.
        df(2,2) = -(1.+xi(1))*(1.-xi(3))/8.
        df(2,3) = -(1.+xi(1))*(1.-xi(2))/8.
        df(3,1) = (1.+xi(2))*(1.-xi(3))/8.
        df(3,2) = (1.+xi(1))*(1.-xi(3))/8.
        df(3,3) = -(1.+xi(1))*(1.+xi(2))/8.
        df(4,1) = -(1.+xi(2))*(1.-xi(3))/8.
        df(4,2) = (1.-xi(1))*(1.-xi(3))/8.
        df(4,3) = -(1.-xi(1))*(1.+xi(2))/8.
        df(5,1) = -(1.-xi(2))*(1.+xi(3))/8.
        df(5,2) = -(1.-xi(1))*(1.+xi(3))/8.
        df(5,3) = (1.-xi(1))*(1.-xi(2))/8.
        df(6,1) = (1.-xi(2))*(1.+xi(3))/8.
        df(6,2) = -(1.+xi(1))*(1.+xi(3))/8.
        df(6,3) = (1.+xi(1))*(1.-xi(2))/8.
        df(7,1) = (1.+xi(2))*(1.+xi(3))/8.
        df(7,2) = (1.+xi(1))*(1.+xi(3))/8.
        df(7,3) = (1.+xi(1))*(1.+xi(2))/8.
        df(8,1) = -(1.+xi(2))*(1.+xi(3))/8.
        df(8,2) = (1.-xi(1))*(1.+xi(3))/8.
        df(8,3) = (1.-xi(1))*(1.+xi(2))/8.
      else if (n_nodes == 20) then
        f(1)=(1.-xi(1))*(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)/8.
        f(2)=(1.+xi(1))*(1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)/8.
        f(3)=(1.+xi(1))*(1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)/8.
        f(4)=(1.-xi(1))*(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)/8.
        f(5)=(1.-xi(1))*(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)/8.
        f(6)=(1.+xi(1))*(1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)/8.
        f(7)=(1.+xi(1))*(1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)/8.
        f(8)=(1.-xi(1))*(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)/8.
        f(9) = (1.-xi(1)**2.)*(1.-xi(2))*(1.-xi(3))/4.
        f(10) = (1.+xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
        f(11) = (1.-xi(1)**2.)*(1.+xi(2))*(1.-xi(3))/4.
        f(12) = (1.-xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
        f(13) = (1.-xi(1)**2.)*(1.-xi(2))*(1.+xi(3))/4.
        f(14) = (1.+xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
        f(15) = (1.-xi(1)**2.)*(1.+xi(2))*(1.+xi(3))/4.
        f(16) = (1.-xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
        f(17) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
        f(18) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
        f(19) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
        f(20) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
        df(1,1) = (-(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(1,2) = (-(1.-xi(1))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(1,3) = (-(1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.

        df(2,1) = ((1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(2,2) = (-(1.+xi(1))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     1          -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(2,3) = (-(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-xi(3)-2.)
     1           -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.

        df(3,1) = ((1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(3,2) = ((1.+xi(1))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(3,3) = (-(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-xi(3)-2.)
     1           -(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.

        df(4,1) = (-(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(4,2) = ((1.-xi(1))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     1            +(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(4,3) = (-(1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(5,1) = (-(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(5,2) = (-(1.-xi(1))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(5,3) = ((1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,1) = ((1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,2) = (-(1.+xi(1))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     1           -(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,3) = ((1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(7,1) = ((1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(7,2) = ((1.+xi(1))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(7,3) = ((1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,1) = (-(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,2) = ((1.-xi(1))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,3) = ((1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(9,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
        df(9,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
        df(9,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
        df(10,1)  = (1.-xi(2)**2.)*(1.-xi(3))/4.
        df(10,2)  = -2.*xi(2)*(1.+xi(1))*(1.-xi(3))/4.
        df(10,3)  = -(1.-xi(2)**2.)*(1.+xi(1))/4.
        df(11,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
        df(11,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
        df(11,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
        df(12,1)  = -(1.-xi(2)**2.)*(1.-xi(3))/4.
        df(12,2)  = -2.*xi(2)*(1.-xi(1))*(1.-xi(3))/4.
        df(12,3)  = -(1.-xi(2)**2.)*(1.-xi(1))/4.
        df(13,1)  = -2.*xi(1)*(1.-xi(2))*(1.+xi(3))/4.
        df(13,2)  = -(1.-xi(1)**2.)*(1.+xi(3))/4.
        df(13,3)  = (1.-xi(1)**2.)*(1.-xi(2))/4.
        df(14,1)  = (1.-xi(2)**2.)*(1.+xi(3))/4.
        df(14,2)  = -2.*xi(2)*(1.+xi(1))*(1.+xi(3))/4.
        df(14,3)  = (1.-xi(2)**2.)*(1.+xi(1))/4.
        df(15,1)  = 2.*xi(1)*(1.+xi(2))*(1.+xi(3))/4.
        df(15,2)  = (1.-xi(1)**2.)*(1.+xi(3))/4.
        df(15,3)  = (1.-xi(1)**2.)*(1.+xi(2))/4.
        df(16,1)  = -(1.-xi(2)**2.)*(1.+xi(3))/4.
        df(16,2)  = -2.*xi(2)*(1.-xi(1))*(1.+xi(3))/4.
        df(16,3)  = (1.-xi(2)**2.)*(1.-xi(1))/4.
        df(17,1) = -(1.-xi(2))*(1.-xi(3)**2.)/4.
        df(17,2) = -(1.-xi(1))*(1.-xi(3)**2.)/4.
        df(17,3) = -xi(3)*(1.-xi(1))*(1.-xi(2))/2.
        df(18,1) = (1.-xi(2))*(1.-xi(3)**2.)/4.
        df(18,2) = -(1.+xi(1))*(1.-xi(3)**2.)/4.
        df(18,3) = -xi(3)*(1.+xi(1))*(1.-xi(2))/2.
        df(19,1) = (1.+xi(2))*(1.-xi(3)**2.)/4.
        df(19,2) = (1.+xi(1))*(1.-xi(3)**2.)/4.
        df(19,3) = -xi(3)*(1.+xi(1))*(1.+xi(2))/2.
        df(20,1) = -(1.+xi(2))*(1.-xi(3)**2.)/4.
        df(20,2) = (1.-xi(1))*(1.-xi(3)**2.)/4.
        df(20,3) = -xi(3)*(1.-xi(1))*(1.+xi(2))/2.
      endif


      end subroutine abq_UEL_3D_shapefunctions




      subroutine abq_UEL_invert3d(A,A_inverse,determinant)

      double precision, intent(in) :: A(3,3)
      double precision, intent(out) :: A_inverse(3,3)
      double precision, intent(out) :: determinant

      double precision COFACTOR(3,3)

!   Compute inverse and determinant of 3x3 matrix

      determinant =   A(1,1)*A(2,2)*A(3,3)
     1   - A(1,1)*A(2,3)*A(3,2)
     2   - A(1,2)*A(2,1)*A(3,3)
     3   + A(1,2)*A(2,3)*A(3,1)
     4   + A(1,3)*A(2,1)*A(3,2)
     5   - A(1,3)*A(2,2)*A(3,1)

      IF (determinant==0.d0) THEN
        write(6,*) ' Error in subroutine abq_UEL_inver3d'
        write(6,*) ' A 3x3 matrix has a zero determinant'
        stop
      endif
      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      A_inverse = transpose(COFACTOR) / determinant


      end subroutine abq_UEL_invert3d



!       subroutine PML_alpha_beta_function(PROPS,x1,x2,x3,PML_alpha_beta)

!       implicit none

!       double precision, intent(in) :: PROPS(*)
!       double precision, intent(in) :: x1, x2, x3
     
!       double precision, intent(out) :: PML_alpha_beta(2,3)

!       integer EleType_arg

!       double precision E, xnu, rho, PML_L, afp, PML_Rcoef
!       double precision x1_0, x2_0, x3_0, n1, n2, n3, cp_ref, PML_b
!       double precision alpha_0, beta_0
!       double precision RD_half_width_x, RD_half_width_y, RD_depth

!       E = PROPS(1)
!       xnu = PROPS(2)
!       rho = PROPS(3)
!       EleType_arg = floor(PROPS(4))
!       PML_L = PROPS(5)
!       afp = PROPS(6)
!       PML_Rcoef = PROPS(7)
!       RD_half_width_x = PROPS(8)
!       RD_half_width_y = PROPS(9)
!       RD_depth = PROPS(10)
      

!       cp_ref = SQRT(E *(1.d0-xnu)/rho/(1.d0+xnu)/(1.d0-2.d0*xnu))

!       IF(x2 < -RD_half_width_y) then
!           IF(x1 < -RD_half_width_x) then
!               IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 15
                  
!               else
                  
!                   EleType_arg = 6
                  
!               endif
              
!           ELSEIF(x1 < RD_half_width_x) then
!                IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 11
                  
!               else
                  
!                   EleType_arg = 2
                  
!               endif    
              
!           Else
!               IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 16
                  
!               else
                  
!                   EleType_arg = 7
                  
!               endif  
!           endif
          
!       ELSEIF(x2 < RD_half_width_y) then
!           IF(x1 < -RD_half_width_x) then
!               IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 14
                  
!               else
                  
!                   EleType_arg = 5
                  
!               endif
              
!           ELSEIF(x1 < RD_half_width_x) then
!                IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 10
                  
!               else
                  
!                   EleType_arg = 1
                  
!               endif    
              
!           Else
!               IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 12
                  
!               else
                  
!                   EleType_arg = 3
                  
!               endif  
!           endif  
          
!       Else
!           IF(x1 < -RD_half_width_x) then
!               IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 18
                  
!               else
                  
!                   EleType_arg = 9
                  
!               endif
              
!           ELSEIF(x1 < RD_half_width_x) then
!                IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 13
                  
!               else
                  
!                   EleType_arg = 4
                  
!               endif    
              
!           Else
!               IF(x3 < -RD_depth) then
                  
!                   EleType_arg = 17
                  
!               else
                  
!                   EleType_arg = 8
                  
!               endif  
!           endif  
          
!       endif
!       ! write(6,*) "EleType_arg ",EleType_arg
!       select case (EleType_arg) 
!       case(1) !Regular domain (do nothing) 
!           n1 = 0.d0 
!           n2 = 0.d0 
!           n3 = 0.d0 
!           x1_0 = 0
!           x2_0 = 0
!           x3_0 = 0

!       case(2) ! Top view, bottom surface PML
!           n1 = 0.d0
!           n2 = -1.d0 
!           n3 = 0.d0
!           x1_0 = 0
!           x2_0 = -1.d0* RD_half_width_y
!           x3_0 = 0

!       case(3) ! Top view, right surface PML
!           n1 = 1.d0
!           n2 = 0.d0 
!           n3 = 0.d0
!           x1_0 = 1.d0* RD_half_width_x
!           x2_0 = 0
!           x3_0 = 0

!       case(4) ! Top view, top surface PML
!           n1 = 0.d0
!           n2 = 1.d0 
!           n3 = 0.d0
!           x1_0 = 0
!           x2_0 = 1.d0* RD_half_width_y
!           x3_0 = 0

!       case(5) ! Top view, left surface PML
!           n1 = -1.d0
!           n2 = 0.d0 
!           n3 = 0.d0
!           x1_0 = -1.d0* RD_half_width_x
!           x2_0 = 0
!           x3_0 = 0          

!       case(6) ! Top view, left-bottom column PML
!           n1 = -1.d0
!           n2 = -1.d0 
!           n3 = 0.d0
!           x1_0 = -1.d0* RD_half_width_x
!           x2_0 = -1.d0* RD_half_width_y
!           x3_0 = 0

!       case(7) ! Top view, right-bottom column PML
!           n1 = 1.d0
!           n2 = -1.d0 
!           n3 = 0.d0
!           x1_0 = 1.d0* RD_half_width_x
!           x2_0 = -1.d0* RD_half_width_y
!           x3_0 = 0

!       case(8) ! Top view, right-top column PML
!           n1 = 1.d0
!           n2 = 1.d0 
!           n3 = 0.d0
!           x1_0 = 1.d0* RD_half_width_x
!           x2_0 = 1.d0* RD_half_width_y
!           x3_0 = 0

!       case(9) ! Top view, left-top column PML
!           n1 = -1.d0
!           n2 = 1.d0 
!           n3 = 0.d0
!           x1_0 = -1.d0* RD_half_width_x
!           x2_0 = 1.d0* RD_half_width_y
!           x3_0 = 0

!       case(10) ! Top view, bottom surface, surface PML
!           n1 = 0.d0
!           n2 = 0.d0 
!           n3 = -1.d0
!           x1_0 = 0
!           x2_0 = 0
!           x3_0 = -1.d0* RD_depth    

!       case(11) ! Top view, bottom surface, bottom column PML
!           n1 = 0.d0
!           n2 = -1.d0 
!           n3 = -1.d0
!           x1_0 = 0
!           x2_0 = -1.d0* RD_half_width_y
!           x3_0 = -1.d0* RD_depth    

!       case(12) ! Top view, bottom surface, right column PML
!           n1 = 1.d0
!           n2 = 0.d0 
!           n3 = -1.d0
!           x1_0 = 1.d0* RD_half_width_x
!           x2_0 = 0
!           x3_0 = -1.d0* RD_depth  

!       case(13) ! Top view, bottom surface, top column PML
!           n1 = 0.d0
!           n2 = 1.d0 
!           n3 = -1.d0
!           x1_0 = 0
!           x2_0 = 1.d0* RD_half_width_y
!           x3_0 = -1.d0* RD_depth 

!       case(14) ! Top view, bottom surface, left column PML
!           n1 = -1.d0
!           n2 = 0.d0 
!           n3 = -1.d0
!           x1_0 = -1.d0* RD_half_width_x
!           x2_0 = 0
!           x3_0 = -1.d0* RD_depth 

!       case(15) ! Top view, bottom surface, left-bottom corner PML
!           n1 = -1.d0
!           n2 = -1.d0 
!           n3 = -1.d0
!           x1_0 = -1.d0* RD_half_width_x
!           x2_0 = -1.d0* RD_half_width_y
!           x3_0 = -1.d0* RD_depth    

!       case(16) ! Top view, bottom surface, right-bottom corner PML
!           n1 = 1.d0
!           n2 = -1.d0 
!           n3 = -1.d0
!           x1_0 = 1.d0* RD_half_width_x
!           x2_0 = -1.d0* RD_half_width_y
!           x3_0 = -1.d0* RD_depth  

!       case(17) ! Top view, bottom surface, right-top corner PML
!           n1 = 1.d0
!           n2 = 1.d0 
!           n3 = -1.d0
!           x1_0 = 1.d0* RD_half_width_x
!           x2_0 = 1.d0* RD_half_width_y
!           x3_0 = -1.d0* RD_depth 

!       case(18) ! Top view, bottom surface, left-top corner PML
!           n1 = -1.d0
!           n2 = 1.d0 
!           n3 = -1.d0
!           x1_0 = -1.d0* RD_half_width_x
!           x2_0 = 1.d0* RD_half_width_y
!           x3_0 = -1.d0* RD_depth 

!       end select

!       PML_b = PML_L / 1.d0   !characteristic length (average element size in the PML domain) 
      
!       alpha_0 = ((afp+1)*PML_b) / (2.d0*PML_L )*LOG10(1.d0 / PML_Rcoef)
!       beta_0 = ((afp+1)*cp_ref) / (2.d0*PML_L )*LOG10(1.d0 / PML_Rcoef)
      
! !      alpha_0 = 5.d0
! !      beta_0 = 800.d0

!       write(6,*) , "alpha_____0", alpha_0
!       write(6,*) , "beta______0", beta_0

!       PML_alpha_beta(1,1) = 1.d0 + alpha_0*((x1 -x1_0) * n1 /PML_L)**afp 
!       PML_alpha_beta(1,2) = 1.d0 + alpha_0*((x2 -x2_0) * n2 /PML_L)**afp
!       PML_alpha_beta(1,3) = 1.d0 + alpha_0*((x3 -x3_0) * n3 /PML_L)**afp

!       PML_alpha_beta(2,1) = beta_0*((x1 -x1_0) * n1 /PML_L)**afp 
!       PML_alpha_beta(2,2) = beta_0*((x2 -x2_0) * n2 /PML_L)**afp
!       PML_alpha_beta(2,3) = beta_0*((x3 -x3_0) * n3 /PML_L)**afp

!       IF (EleType_arg .EQ.1) THEN
!         PML_alpha_beta(1:2,1:3) = 0.d0
!       end if

!       end subroutine PML_alpha_beta_function
