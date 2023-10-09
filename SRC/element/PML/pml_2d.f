c------------------------------------------------------------------------------
c     Complete UEL file for Perfectly-Matched-Layer (PML) (2D Plane Strain)
c------------------------------------------------------------------------------
c     Note: 1) Please do not distribute the material model routines, and if
c              anyone asks, please refer them to Prof. Ertugrul Taciroglu
c              (etacir@ucla.edu)
c           2) Please properly acknowledge our related publications when you
c              publish your results
c     Ref: Zhang W, Esmaeilzadeh Seylabi E, Taciroglu E. An ABAQUS toolbox for 
c          soil-structure interaction analysis. Computers and Geotechnics. 114 
c          (2019): 103143.
c------------------------------------------------------------------------------
c     Authors: Wenyang Zhang (zwyll@ucla.edu)
c              University of California, Los Angeles
c------------------------------------------------------------------------------     
c     Function inputs given by the user:
c
c     E = PROPS(1)                     --- Young's modulus
c     xnu = PROPS(2)                   --- Poisson's Ratio
c     rho = PROPS(3)                   --- Density
c     EleType_pos = floor(PROPS(4))    --- Element type, See line
c     PML_L = PROPS(5)                 --- Thickness of the PML
c     afp = PROPS(6)                   --- Coefficient m, typically m = 2
c     PML_Rcoef = PROPS(7)             --- Coefficient R, typically R = 1e-8
c     RD_half_width = PROPS(8)         --- Halfwidth of the regular domain
c     RD_depth = PROPS(9)              --- Depth of the regular domain
c     Damp_alpha = PROPS(10) 
c     Damp_beta = PROPS(11)
!=========================== ABAQUS format user element subroutine ===================
      

      SUBROUTINE PML_2D(KMATRX,CMATRX,MMATRX,GMATRX,NDOFEL,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE)

      DOUBLE PRECISION PROPS,COORDS,KMATRX,CMATRX,MMATRX,GMATRX
      INTEGER NDOFEL,NPROPS,MCRD,NNODE

    !
    !
    !
      DIMENSION KMATRX(NDOFEL,NDOFEL),CMATRX(NDOFEL,NDOFEL),
     1   MMATRX(NDOFEL,NDOFEL),PROPS(*),COORDS(MCRD,NNODE),
     2   GMATRX(NDOFEL,NDOFEL)

    !
    !       Variables that must be computed in this routine

    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node


    !
    ! Local Variables

      real*8 ZERO,ONE,HALF       
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0 )

      real*8 coef_alpha
      PARAMETER (coef_alpha = 1.d0/12.d0)

      integer      :: i,j,n_points,kint, nfacenodes, ipoin, ksize
      integer      :: face_node_list(3)                       ! List of nodes on an element face
    !
      double precision  ::  coords_new(2,4)                  ! Arranged coords

      double precision  ::  xi(2,9)                          ! Area integration points
      double precision  ::  w(9)                             ! Area integration weights
      double precision  ::  N(9)                             ! 2D shape functions
      double precision  ::  dNdxi(9,2)                       ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                        ! Spatial derivatives
      double precision  ::  dxdxi(2,2)                       ! Derivative of spatial coords wrt normalized coords

    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(2,3)                  ! Coords of nodes on an element face
      double precision  ::  xi1(6)                            ! 1D integration points
      double precision  ::  w1(6)                              ! Integration weights
      double precision  ::  N1(3)                             ! 1D shape functions
      double precision  ::  dN1dxi(3)                         ! 1D shape function derivatives
      double precision  ::  norm(2)                           ! Normal to an element face
      double precision  ::  dxdxi1(2)                         ! Derivative of 1D spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(4)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress(4)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  D(4,4)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(4,22)                           ! strain = B*(dof_total)
      double precision  ::  ktemp(22,22)                      ! Temporary stiffness (for incompatible mode elements)
      double precision  ::  dxidx(2,2), determinant, det0     ! Jacobian inverse and determinant
      double precision  ::  E, xnu, D44, D11, D12             ! Material properties

      double precision  ::  ALPHA, BETA, GAMMA, DADU, DVDU    ! Newmark Integration properties
      double precision  ::  Phi(2,18)                         ! Matrix consists of shape functions [N1 N2 ... N20]
c      double precision  ::  MMATRX(NDOFEL,NDOFEL)             ! Element mass matrix
c      double precision  ::  KMATRX(NDOFEL,NDOFEL)             ! Element stiffness matrix
c      double precision  ::  CMATRX(NDOFEL,NDOFEL)             ! Element damping matrix
c      double precision  ::  GMATRX(NDOFEL,NDOFEL)             ! Element G matrix
      double precision  ::  rho                               ! Density of the material

      double precision  ::  PML_L,afp,PML_Rcoef               ! Parameters of PML
      double precision  ::  RD_half_width                     ! Parameters of PML
      double precision  ::  RD_depth                          ! Parameters of PML
      integer           ::  EleType_pos                       ! Parameters of PML

      double precision  ::  x1,x2,x3,PML_alpha_beta(2,2)      ! Coordinates x and y, alphas and betas    
      double precision  ::  coef_a, coef_b, coef_c            ! Coefficients of the PML
      double precision  ::  coef_Le(2), coef_Lp(2)            ! Coefficients of the PML

      double precision  ::  K_RD(8,8), Kxx(4,4), Kxy(4,4)               ! Stiffness matrices
      double precision  ::  Kyx(4,4), Kyy(4,4)                          ! Stiffness matrices
      double precision  ::  Phi_x(4), Phi_y(4)                          ! Derivative of shape functions for displacement fields [N1,i N2,i ... N4,i]^T
      double precision  ::  lambda, mu                                  ! Lame' constants
      double precision  ::  C_RD(8,8)                                   ! Damping matrix for regular domain                              
      double precision  ::  M_RD(8,8)                                   ! Mass matrix for regular domain                              
      double precision  ::  M_a(8,8),M_b(8,8),M_c(8,8)                  ! Mass matrices
      double precision  ::  N_a(12,12),N_b(12,12),N_c(12,12)            ! N matrices for PML
      double precision  ::  A_eu(8,12), A_pu(8,12)                      ! A_iu matrices for PML, where i = e, w, p
      double precision  ::  A_el(8,12), A_pl(8,12)                      ! A_il matrices for PML, where i = e, w, p
      double precision  ::  K_PML(20,20), M_PML(20,20), C_PML(20,20)    ! Stiffnee, mass and damping matrices for PML
      double precision  ::  G_PML(20,20)                                ! G matrices for PML

      double precision  ::  d_bar_n(20), d_bar_n1(20)                   ! d_bar at previous and current step
      double precision  ::  U_n(20), V_n(20), A_n(20)                   ! U, V, A at previous step

      integer           ::  NodeDOF                                     ! Number of DOF per node

      integer           ::  FaceID                                      ! Face number 
      double precision  ::  IntegrationPointXY(2)                       ! Coordinates of the integration points 

      double precision  ::  Damp_alpha, Damp_beta                       ! Damping coefficients for Rayleigh Damping
      
      double precision  ::  Vs

      NodeDOF = NDOFEL/NNODE

    !
    !     Example ABAQUS PML2D implementing 2D linear elastic elements
    !     Includes option for incompatible mode elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio


      if (NNODE == 3) n_points = 1              ! Linear triangle
      if (NNODE == 4) n_points = 9              ! Linear rectangle
      if (NNODE == 6) n_points = 4              ! Quadratic triangle
      if (NNODE == 8) n_points = 9              ! Serendipity rectangle
      if (NNODE == 9) n_points = 9              ! Quadratic rect

      call abq_PML2D_integrationpoints(n_points, NNODE, xi, w)

      if (NDOFEL<2*NNODE) then
        write(6,*) ' Error in abaqus PML2D '
        write(6,*) ' Variable NDOFEL must exceed 2*NNODE'
        stop
      endif

      MMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      KMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      CMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      GMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      E = PROPS(1)
      xnu = PROPS(2)
      rho = PROPS(3)
      EleType_pos = floor(PROPS(4))
      PML_L = PROPS(5)
      afp = PROPS(6)
      PML_Rcoef = PROPS(7)
      RD_half_width = PROPS(8)
      RD_depth = PROPS(9)
      Damp_alpha = PROPS(10)
      Damp_beta = PROPS(11)

      lambda = xnu*E/( (1.d0+xnu)*(1.d0-2.D0*xnu) )
      mu = 0.5D0*E/(1.d0+xnu)
      
      Vs = sqrt(mu/rho)
!      write(6,*) 'Vs: ', Vs
!      write(6,*) 'Props: ', Props(1:11)

      M_RD = 0.d0
      M_a = 0.d0
      M_b = 0.d0
      M_c = 0.d0

      N_a = 0.d0
      N_b = 0.d0
      N_c = 0.d0

      A_eu = 0.d0
      A_pu = 0.d0

      A_el = 0.d0
      A_pl = 0.d0

      K_RD = 0.d0

      K_PML = 0.d0
      C_PML = 0.d0
      M_PML = 0.d0
      G_PML = 0.d0

      d_bar_n = 0.d0
      d_bar_n1 = 0.d0

      U_n = 0.d0
      V_n = 0.d0
      A_n = 0.d0

    !     --  Loop over integration points
!      write(6,*) 'Node1:', ' ', coords(1,1),' ',coords(2,1)
!      write(6,*) 'Node2:', ' ', coords(1,2),' ',coords(2,2)
!      write(6,*) 'Node3:', ' ', coords(1,3),' ',coords(2,3)
!      write(6,*) 'Node4:', ' ', coords(1,4),' ',coords(2,4)
      

!      write(6,*) 'nnode: ', NNODE, 'n_points', n_points

      do kint = 1, n_points
         call abq_PML2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
         dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))       

         determinant =   dxdxi(1,1)*dxdxi(2,2)
     1   - dxdxi(1,2)*dxdxi(2,1)
         
         dxidx = 0.d0
         dxidx(1,1) = dxdxi(2,2)/determinant
         dxidx(1,2) = -dxdxi(1,2)/determinant
         dxidx(2,1) = -dxdxi(2,1)/determinant
         dxidx(2,2) = dxdxi(1,1)/determinant

         dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)      
         
         


c         ALPHA = PARAMS(1)
c         BETA  = PARAMS(2)
c         GAMMA = PARAMS(3)
         
c         DADU = ONE/(BETA*DTIME**2)
c         DVDU = GAMMA/(BETA*DTIME)
         
         Phi = 0.d0
         Phi(1,1:2*NNODE-1:2) = N(1:NNODE)
         Phi(2,2:3*NNODE:2) = N(1:NNODE)

         Phi_x = 0.d0
         Phi_y = 0.d0
         Phi_x(1:NNODE) = dNdx(1:NNODE,1)
         Phi_y(1:NNODE) = dNdx(1:NNODE,2)
         
         x1 = 0.d0
         x2 = 0.d0
         
         do i = 1,NNODE
            x1 = x1 + N(i)*coords(1,i)
            x2 = x2 + N(i)*coords(2,i)
         end do

!         write(6,*) "x1 x2: ", x1, ' ' , x2

         call PML2D_alpha_beta_function(PROPS,x1,x2,PML_alpha_beta)
         
         coef_a = PML_alpha_beta(1,1)*PML_alpha_beta(1,2)
         
         coef_b = PML_alpha_beta(1,1)*PML_alpha_beta(2,2) + 
     1        PML_alpha_beta(2,1)*PML_alpha_beta(1,2)
         
         coef_c = PML_alpha_beta(2,1)*PML_alpha_beta(2,2)
         
         coef_Le = 0.d0
         coef_Lp = 0.d0
         
         coef_Le(1) = PML_alpha_beta(1,1)
         coef_Le(2) = PML_alpha_beta(1,2)
         
         coef_Lp(1) = PML_alpha_beta(2,1)
         coef_Lp(2) = PML_alpha_beta(2,2)

!         write(6,*) 'COEF:',coef_a,' ', coef_b,' ', coef_c
!         write(6,*) 'coef_Lp:',coef_Lp,' ', 'coef_Le:',' ', coef_Le
         
         do i = 1,NNODE
            do j = 1,NNODE
               Kxx(i,j) = (lambda + 2.d0*mu) * Phi_x(i)*Phi_x(j) + 
     1              mu*(Phi_y(i)*Phi_y(j))     
               Kyy(i,j) = (lambda + 2.d0*mu) * Phi_y(i)*Phi_y(j) + 
     1              mu*(Phi_x(i)*Phi_x(j))               
               
               Kxy(i,j) = lambda * Phi_x(i)*Phi_y(j) + 
     1              mu * Phi_y(i)*Phi_x(j)
               
               
               
               M_RD(i,j) = M_RD(i,j) + rho*N(i)*N(j)*w(kint)*determinant 
               M_a(i,j) = M_a(i,j) + 
     1              coef_a*rho*N(i)*N(j)*w(kint)*determinant 
               M_b(i,j) = M_b(i,j) + 
     1              coef_b*rho*N(i)*N(j)*w(kint)*determinant
               M_c(i,j) = M_c(i,j) + 
     1              coef_c*rho*N(i)*N(j)*w(kint)*determinant
               
               
               A_eu(i,j) = A_eu(i,j) + 
     1              Phi_x(i)*N(j)*coef_Le(2)*w(kint)*determinant
               A_eu(i,j+NNODE*2) = A_eu(i,j+NNODE*2) + 
     1              Phi_y(i)*N(j)*coef_Le(1)*w(kint)*determinant 
               
               A_eu(i+NNODE,j+NNODE) = A_eu(i+NNODE,j+NNODE) + 
     1              Phi_y(i)*N(j)*coef_Le(1)*w(kint)*determinant 
               A_eu(i+NNODE,j+NNODE*2) = A_eu(i+NNODE,j+NNODE*2) + 
     1              Phi_x(i)*N(j)*coef_Le(2)*w(kint)*determinant 
               
               A_pu(i,j) = A_pu(i,j) + 
     1              Phi_x(i)*N(j)*coef_Lp(2)*w(kint)*determinant
               A_pu(i,j+NNODE*2) = A_pu(i,j+NNODE*2) + 
     1              Phi_y(i)*N(j)*coef_Lp(1)*w(kint)*determinant 
               
               A_pu(i+NNODE,j+NNODE) = A_pu(i+NNODE,j+NNODE) + 
     1              Phi_y(i)*N(j)*coef_Lp(1)*w(kint)*determinant 
               A_pu(i+NNODE,j+NNODE*2) = A_pu(i+NNODE,j+NNODE*2) + 
     1              Phi_x(i)*N(j)*coef_Lp(2)*w(kint)*determinant 
               
            end do
         end do
         
          Kyx(1:NNODE,1:NNODE) = transpose(Kxy(1:NNODE,1:NNODE)) 

          K_RD(1:NNODE,1:NNODE) = K_RD(1:NNODE,1:NNODE) + 
     1                          Kxx(1:NNODE,1:NNODE)*w(kint)*determinant 
          K_RD(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1            K_RD(NNODE+1:2*NNODE,NNODE+1:2*NNODE) + 
     2                          Kyy(1:NNODE,1:NNODE)*w(kint)*determinant 

          K_RD(1:NNODE,NNODE+1:2*NNODE) =
     1                    K_RD(1:NNODE,NNODE+1:2*NNODE) +  
     2                          Kxy(1:NNODE,1:NNODE)*w(kint)*determinant  

          K_RD(NNODE+1:2*NNODE,1:NNODE) =
     1                    K_RD(NNODE+1:2*NNODE,1:NNODE) +  
     2                          Kyx(1:NNODE,1:NNODE)*w(kint)*determinant  

       end do

        M_RD(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_RD(1:NNODE,1:NNODE)

        C_RD(1:2*NNODE,1:2*NNODE) = Damp_alpha*M_RD(1:2*NNODE,1:2*NNODE)
     1                            + Damp_beta*K_RD(1:2*NNODE,1:2*NNODE)
        
        M_a(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_a(1:NNODE,1:NNODE)
        
        M_b(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_b(1:NNODE,1:NNODE)
        
        M_c(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = M_c(1:NNODE,1:NNODE)
       
        N_a(1:NNODE,1:NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho*(lambda+2.d0*mu)/mu/4.d0/(lambda+mu)
        N_a(1:NNODE,NNODE+1:2*NNODE) = 
     1  -M_a(1:NNODE,1:NNODE)/rho*(lambda)/mu/4.d0/(lambda+mu)
        
        N_a(NNODE+1:2*NNODE,1:NNODE) = 
     1  -M_a(1:NNODE,1:NNODE)/rho*(lambda)/mu/4.d0/(lambda+mu)
        N_a(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho*(lambda+2.d0*mu)/mu/4.d0/(lambda+mu)
        
        N_a(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) = 
     1   M_a(1:NNODE,1:NNODE)/rho/mu


        N_b(1:NNODE,1:NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho*(lambda+2.d0*mu)/mu/4.d0/(lambda+mu)
        N_b(1:NNODE,NNODE+1:2*NNODE) = 
     1  -M_b(1:NNODE,1:NNODE)/rho*(lambda)/mu/4.d0/(lambda+mu)
        
        N_b(NNODE+1:2*NNODE,1:NNODE) = 
     1  -M_b(1:NNODE,1:NNODE)/rho*(lambda)/mu/4.d0/(lambda+mu)
        N_b(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho*(lambda+2.d0*mu)/mu/4.d0/(lambda+mu)
        
        N_b(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) = 
     1   M_b(1:NNODE,1:NNODE)/rho/mu


        N_c(1:NNODE,1:NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho*(lambda+2.d0*mu)/mu/4.d0/(lambda+mu)
        N_c(1:NNODE,NNODE+1:2*NNODE) = 
     1  -M_c(1:NNODE,1:NNODE)/rho*(lambda)/mu/4.d0/(lambda+mu)
        
        N_c(NNODE+1:2*NNODE,1:NNODE) = 
     1  -M_c(1:NNODE,1:NNODE)/rho*(lambda)/mu/4.d0/(lambda+mu)
        N_c(NNODE+1:2*NNODE,NNODE+1:2*NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho*(lambda+2.d0*mu)/mu/4.d0/(lambda+mu)
        
        N_c(2*NNODE+1:3*NNODE,2*NNODE+1:3*NNODE) = 
     1   M_c(1:NNODE,1:NNODE)/rho/mu

        IF (EleType_pos .EQ. 1) THEN

          M_PML(1:NNODE*2,1:NNODE*2) = M_RD(1:NNODE*2,1:NNODE*2) + 
     1                             M_a(1:NNODE*2,1:NNODE*2)
          M_PML(NNODE*2+1:NNODE*5,NNODE*2+1:NNODE*5) = 
     1                                      -N_a(1:NNODE*3,1:NNODE*3)

          C_PML(1:NNODE*2,1:NNODE*2) = C_RD(1:NNODE*2,1:NNODE*2) + 
     1                             M_b(1:NNODE*2,1:NNODE*2)
          
          C_PML(1:NNODE*2,NNODE*2+1:NNODE*5) = A_eu(1:NNODE*2,1:NNODE*3)
          C_PML(NNODE*2+1:NNODE*5,1:NNODE*2) = 
     1                         transpose(A_eu(1:NNODE*2,1:NNODE*3))
          C_PML(NNODE*2+1:NNODE*5,NNODE*2+1:NNODE*5) = 
     1                                      -N_b(1:NNODE*3,1:NNODE*3)

          K_PML(1:NNODE*2,1:NNODE*2) = K_RD(1:NNODE*2,1:NNODE*2) + 
     1                             M_c(1:NNODE*2,1:NNODE*2)
          K_PML(1:NNODE*2,NNODE*2+1:NNODE*5) = A_pu(1:NNODE*2,1:NNODE*3)
          K_PML(NNODE*2+1:NNODE*5,1:NNODE*2) = 
     1                         transpose(A_pu(1:NNODE*2,1:NNODE*3))
          K_PML(NNODE*2+1:NNODE*5,NNODE*2+1:NNODE*5) = 
     1                                      -N_c(1:NNODE*3,1:NNODE*3)

          else

          M_PML(1:NNODE*2,1:NNODE*2) = M_a(1:NNODE*2,1:NNODE*2)
          M_PML(1:NNODE*2,NNODE*2+1:NNODE*5) = 
     1                            A_eu(1:NNODE*2,1:NNODE*3)*Damp_beta
          M_PML(NNODE*2+1:NNODE*5,NNODE*2+1:NNODE*5) = 
     1                                      -N_a(1:NNODE*3,1:NNODE*3)

          C_PML(1:NNODE*2,1:NNODE*2) = M_b(1:NNODE*2,1:NNODE*2) + 
     1                   M_a(1:NNODE*2,1:NNODE*2)*Damp_alpha     
          C_PML(1:NNODE*2,NNODE*2+1:NNODE*5) = A_eu(1:NNODE*2,1:NNODE*3)
     1       + A_pu(1:NNODE*2,1:NNODE*3)*Damp_beta     
          C_PML(NNODE*2+1:NNODE*5,1:NNODE*2) = 
     1                         transpose(A_eu(1:NNODE*2,1:NNODE*3))
          C_PML(NNODE*2+1:NNODE*5,NNODE*2+1:NNODE*5) = 
     1                                      -N_b(1:NNODE*3,1:NNODE*3)

          K_PML(1:NNODE*2,1:NNODE*2) = M_c(1:NNODE*2,1:NNODE*2) + 
     1                   M_b(1:NNODE*2,1:NNODE*2)*Damp_alpha
          K_PML(1:NNODE*2,NNODE*2+1:NNODE*5) = A_pu(1:NNODE*2,1:NNODE*3)
          K_PML(NNODE*2+1:NNODE*5,1:NNODE*2) = 
     1                         transpose(A_pu(1:NNODE*2,1:NNODE*3))
          K_PML(NNODE*2+1:NNODE*5,NNODE*2+1:NNODE*5) = 
     1                                      -N_c(1:NNODE*3,1:NNODE*3)
          
          G_PML(1:NNODE*2,1:NNODE*2)=M_c(1:NNODE*2,1:NNODE*2)*Damp_alpha


          endif 

        do i = 1,4
          do j = 1,4
            MMATRX((i-1)*5+1:i*5,(j-1)*5+1:j*5) = 
     1                        M_PML(i:NDOFEL-4+i:4,j:NDOFEL-4+j:4)

            CMATRX((i-1)*5+1:i*5,(j-1)*5+1:j*5) = 
     1                        C_PML(i:NDOFEL-4+i:4,j:NDOFEL-4+j:4)

            KMATRX((i-1)*5+1:i*5,(j-1)*5+1:j*5) = 
     1                        K_PML(i:NDOFEL-4+i:4,j:NDOFEL-4+j:4)

            GMATRX((i-1)*5+1:i*5,(j-1)*5+1:j*5) = 
     1                        G_PML(i:NDOFEL-4+i:4,j:NDOFEL-4+j:4)
          end do
        end do
        
        
!        MMATRX = M_PML
!        CMATRX = C_PML
!        KMATRX = K_PML


!      write(6,*) 'K11= ', KMATRX(1,1)
!      write(6,*) 'C11= ', CMATRX(1,1)
!      write(6,*) 'M11= ', MMATRX(1,1)

      return


      END SUBROUTINE PML_2D



      subroutine abq_PML2D_integrationpoints(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(2,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision :: cn,w1,w2,w11,w12,w22

    !         Defines integration points and weights for 2D continuum elements

      if ( n_points==1 ) then
        if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
            xi(1, 1) = 0.D0
            xi(2, 1) = 0.D0
            w(1) = 4.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = 1.D0/3.D0
            w(1) = 1.D0/2.D0
        end if
      else if ( n_points==3 ) then
        xi(1, 1) = 0.5D0
        xi(2, 1) = 0.5D0
        w(1) = 1.D0/6.D0
        xi(1, 2) = 0.D0
        xi(2, 2) = 0.5D0
        w(2) = w(1)
        xi(1, 3) = 0.5D0
        xi(2, 3) = 0.D0
        w(3) = w(1)
      else if ( n_points==4 ) then
        if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
            !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
            !     43
            !     12
            cn = 0.5773502691896260D0
            xi(1, 1) = -cn
            xi(1, 2) = cn
            xi(1, 3) = cn
            xi(1, 4) = -cn
            xi(2, 1) = -cn
            xi(2, 2) = -cn
            xi(2, 3) = cn
            xi(2, 4) = cn
            w(1) = 1.D0
            w(2) = 1.D0
            w(3) = 1.D0
            w(4) = 1.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then
            !     xi integration points for triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = xi(1, 1)
            w(1) = -27.D0/96.D0
            xi(1, 2) = 0.6D0
            xi(2, 2) = 0.2D0
            w(2) = 25.D0/96.D0
            xi(1, 3) = 0.2D0
            xi(2, 3) = 0.6D0
            w(3) = w(2)
            xi(1, 4) = 0.2D0
            xi(2, 4) = 0.2D0
            w(4) = w(2)
        end if

      else if ( n_points==7 ) then
        ! Quintic integration for triangle
        xi(1,1) = 1.d0/3.d0
        xi(2,1) = xi(1,1)
        w(1) = 0.1125d0
        xi(1,2) = 0.0597158717d0
        xi(2,2) = 0.4701420641d0
        w(2) = 0.0661970763d0
        xi(1,3) = xi(2,2)
        xi(2,3) = xi(1,2)
        w(3) = w(2)
        xi(1,4) = xi(2,2)
        xi(2,4) = xi(2,2)
        w(4) = w(2)
        xi(1,5) = 0.7974269853d0
        xi(2,5) = 0.1012865073d0
        w(5) = 0.0629695902d0
        xi(1,6) = xi(2,5)
        xi(2,6) = xi(1,5)
        w(6) = w(5)
        xi(1,7) = xi(2,5)
        xi(2,7) = xi(2,5)
        w(7) = w(5)
      else if ( n_points==9 ) then
        !     3X3 GAUSS INTEGRATION POINTS
        !     789
        !     456
        !     123
        cn = 0.7745966692414830D0
        xi(1, 1) = -cn
        xi(1, 2) = 0.D0
        xi(1, 3) = cn
        xi(1, 4) = -cn
        xi(1, 5) = 0.D0
        xi(1, 6) = cn
        xi(1, 7) = -cn
        xi(1, 8) = 0.D0
        xi(1, 9) = cn
        xi(2, 1) = -cn
        xi(2, 2) = -cn
        xi(2, 3) = -cn
        xi(2, 4) = 0.D0
        xi(2, 5) = 0.D0
        xi(2, 6) = 0.D0
        xi(2, 7) = cn
        xi(2, 8) = cn
        xi(2, 9) = cn
        w1 = 0.5555555555555560D0
        w2 = 0.8888888888888890D0
        w11 = w1*w1
        w12 = w1*w2
        w22 = w2*w2
        w(1) = w11
        w(2) = w12
        w(3) = w11
        w(4) = w12
        w(5) = w22
        w(6) = w12
        w(7) = w11
        w(8) = w12
        w(9) = w11
      end if

      return

      end subroutine abq_PML2D_integrationpoints



      subroutine abq_PML2D_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(2)
      double precision, intent(out) :: f(*)
      double precision, intent(out) :: df(9,2)
      double precision g1, g2, g3, dg1, dg2, dg3
      double precision h1, h2, h3, dh1, dh2, dh3
      double precision z,dzdp, dzdq

            if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = 1.D0 - xi(1) - xi(2)
                df(1, 1) = 1.D0
                df(1, 2) = 0.D0
                df(2, 1) = 0.D0
                df(2, 2) = 1.D0
                df(3, 1) = -1.D0
                df(3, 2) = -1.D0
            else if ( n_nodes==4 ) then
                !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
                !     43
                !     12
                g1 = 0.5D0*(1.D0 - xi(1))
                g2 = 0.5D0*(1.D0 + xi(1))
                h1 = 0.5D0*(1.D0 - xi(2))
                h2 = 0.5D0*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g2*h2
                f(4) = g1*h2
                dg1 = -0.5D0
                dg2 = 0.5D0
                dh1 = -0.5D0
                dh2 = 0.5D0
                df(1, 1) = dg1*h1
                df(2, 1) = dg2*h1
                df(3, 1) = dg2*h2
                df(4, 1) = dg1*h2
                df(1, 2) = g1*dh1
                df(2, 2) = g2*dh1
                df(3, 2) = g2*dh2
                df(4, 2) = g1*dh2

            else if ( n_nodes==6 ) then

                !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
                !          3

                !       6      5

                !     1    4     2

                !     P = L1
                !     Q = L2
                !     Z = 1 - P - Q = L3

                z = 1.D0 - xi(1) - xi(2)
                f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
                f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
                f(3) = (2.D0*z - 1.D0)*z
                f(4) = 4.D0*xi(1)*xi(2)
                f(5) = 4.D0*xi(2)*z
                f(6) = 4.D0*xi(1)*z
                dzdp = -1.D0
                dzdq = -1.D0
                df(1, 1) = 4.D0*xi(1) - 1.D0
                df(2, 1) = 0.D0
                df(3, 1) = 4.D0*z*dzdp - dzdp
                df(4, 1) = 4.D0*xi(2)
                df(5, 1) = 4.D0*xi(2)*dzdp
                df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
                df(1, 2) = 0.D0
                df(2, 2) = 4.D0*xi(2) - 1.D0
                df(3, 2) = 4.D0*z*dzdq - dzdq
                df(4, 2) = 4.D0*xi(1)
                df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
                df(6, 2) = 4.D0*xi(1)*dzdq

            else if ( n_nodes==8 ) then
                !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
                 f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
                 f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
                 f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
                 f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
                 f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
                 f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
                 f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
                 f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
                 df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
                 df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
                 df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
                 df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
                 df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
                 df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
                 df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
                 df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
                 df(5,1) = -xi(1)*(1.-xi(2));
                 df(5,2) = -0.5*(1.-xi(1)*xi(1));
                 df(6,1) = 0.5*(1.-xi(2)*xi(2));
                 df(6,2) = -(1.+xi(1))*xi(2);
                 df(7,1) = -xi(1)*(1.+xi(2));
                 df(7,2) = 0.5*(1.-xi(1)*xi(1));
                 df(8,1) = -0.5*(1.-xi(2)*xi(2));
                 df(8,2) = -(1.-xi(1))*xi(2);
            else if ( n_nodes==9 ) then
                !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
                !     789
                !     456
                !     123
                g1 = -.5D0*xi(1)*(1.D0 - xi(1))
                g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
                g3 = .5D0*xi(1)*(1.D0 + xi(1))
                h1 = -.5D0*xi(2)*(1.D0 - xi(2))
                h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
                h3 = .5D0*xi(2)*(1.D0 + xi(2))
                dg1 = xi(1) - 0.5d0
                dg2 = -2.d0*xi(1)
                dg3 = xi(1) + 0.5d0
                dh1 = xi(2)-0.5d0
                dh2 = -2.d0*xi(2)
                dh3 = xi(2) + 0.5d0
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g3*h1
                f(4) = g1*h2
                f(5) = g2*h2
                f(6) = g3*h2
                f(7) = g1*h3
                f(8) = g2*h3
                f(9) = g3*h3
                df(1,1) = dg1*h1
                df(1,2) = g1*dh1
                df(2,1) = dg2*h1
                df(2,2) = g2*dh1
                df(3,1) = dg3*h1
                df(3,2) = g3*dh1
                df(4,1) = dg1*h2
                df(4,2) = g1*dh2
                df(5,1) = dg2*h2
                df(5,2) = g2*dh2
                df(6,1) = dg3*h2
                df(6,2) = g3*dh2
                df(7,1) = dg1*h3
                df(7,2) = g1*dh3
                df(8,1) = dg2*h3
                df(8,2) = g2*dh3
                df(9,1) = dg3*h3
                df(9,2) = g3*dh3
            end if

      end subroutine abq_PML2D_shapefunctions


      subroutine abq_PML2D_1D_integrationpoints(n_points,n_nodes,xi,w)


      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)



      select case ( n_points )
        case (2)
            xi(1) = .5773502691896257D+00
            xi(2) = -.5773502691896257D+00
            w(1) = .1000000000000000D+01
            w(2) = .1000000000000000D+01
            return
        case (3)
            xi(1) = 0.7745966692414834D+00
            xi(2) = .0000000000000000D+00
            xi(3) = -.7745966692414834D+00
            w(1) = .5555555555555556D+00
            w(2) = .8888888888888888D+00
            w(3) = .5555555555555556D+00
            return
        case (4)
            xi(1) = .8611363115940526D+00
            xi(2) = .3399810435848563D+00
            xi(3) = -.3399810435848563D+00
            xi(4) = -.8611363115940526D+00
            w(1) = .3478548451374538D+00
            w(2) = .6521451548625461D+00
            w(3) = .6521451548625461D+00
            w(4) = .3478548451374538D+00
            return
        case (5)
            xi(1) = .9061798459386640D+00
            xi(2) = .5384693101056831D+00
            xi(3) = .0000000000000000D+00
            xi(4) = -.5384693101056831D+00
            xi(5) = -.9061798459386640D+00
            w(1) = .2369268850561891D+00
            w(2) = .4786286704993665D+00
            w(3) = .5688888888888889D+00
            w(4) = .4786286704993665D+00
            w(5) = .2369268850561891D+00
            return
        case (6)
            xi(1) = .9324695142031521D+00
            xi(2) = .6612093864662645D+00
            xi(3) = .2386191860831969D+00
            xi(4) = -.2386191860831969D+00
            xi(5) = -.6612093864662645D+00
            xi(6) = -.9324695142031521D+00
            w(1) = .1713244923791703D+00
            w(2) = .3607615730481386D+00
            w(3) = .4679139345726910D+00
            w(4) = .4679139345726910D+00
            w(5) = .3607615730481386D+00
            w(6) = .1713244923791703D+00
            return
        case DEFAULT
            write(6,*)'Error in routine abq_PML2D_1D_integrationpoints'
            write(6,*) ' Invalid number of integration points for 1D'
            write(6,*) ' n_points must be between 1 and 6'
            stop
      end select


      end subroutine ABQ_PML2D_1D_integrationpoints



      subroutine PML2D_alpha_beta_function(PROPS,x1,x2,PML_alpha_beta)

      implicit none

      double precision, intent(in) :: PROPS(*)
      double precision, intent(in) :: x1, x2
     
      double precision, intent(out) :: PML_alpha_beta(2,2)

      integer EleType_arg

      double precision E, xnu, rho, PML_L, afp, PML_Rcoef
      double precision x1_0, x2_0, n1, n2, cp_ref, PML_b
      double precision alpha_0, beta_0
      double precision RD_half_width, RD_depth

      E = PROPS(1)
      xnu = PROPS(2)
      rho = PROPS(3)
      EleType_arg = floor(PROPS(4))
      PML_L = PROPS(5)
      afp = PROPS(6)
      PML_Rcoef = PROPS(7)
      RD_half_width = PROPS(8)
      RD_depth = PROPS(9)
      

      cp_ref = SQRT(E *(1.d0-xnu)/rho/(1.d0+xnu)/(1.d0-2.d0*xnu))
      cp_ref = 750.d0
      
      
      IF(x2 < -RD_depth) then
          IF(x1 < -RD_half_width) then
              EleType_arg = 3
          ElSEIF(x1 < RD_half_width) then
              EleType_arg = 4
          ELSE
              EleType_arg = 5
          ENDIF
      ELSE
          IF(x1 < -RD_half_width) then
              EleType_arg = 2
          ElSEIF(x1 < RD_half_width) then
              EleType_arg = 1
          ELSE
              EleType_arg = 6
          ENDIF
      ENDIF
      
!      write(6,*) 'EleType_arg:',' ',EleType_arg


      select case (EleType_arg) 
      case(1) !Regular domain (do nothing) 
          n1 = 0.d0 
          n2 = 0.d0 
          x1_0 = 0
          x2_0 = 0

      case(2) ! Left column PML
          n1 = -1.d0
          n2 = 0.d0 
          x1_0 = -1.d0* RD_half_width
          x2_0 = 0

      case(3) ! Left-bottom corner PML
          n1 = -1.d0
          n2 = -1.d0
          x1_0 = -1.d0* RD_half_width
          x2_0 = -1.d0* RD_depth           

      case(4) ! Bottom row PML
          n1 = 0.d0
          n2 = -1.d0
          x1_0 = 0
          x2_0 = -1.d0* RD_depth 

      case(5) ! Right-bottom corner PML
          n1 = 1.d0
          n2 = -1.d0
          x1_0 = 1.d0* RD_half_width
          x2_0 = -1.d0* RD_depth 

      case(6) ! Right column PML
          n1 = 1.d0
          n2 = 0.d0 
          x1_0 = 1.d0* RD_half_width
          x2_0 = 0

      end select

      PML_b = PML_L / 1.d0   !characteristic length (average element size in the PML domain) 

      alpha_0 = ((afp+1)*PML_b) / (2.d0*PML_L )*LOG10(1.d0 / PML_Rcoef)
      beta_0 = ((afp+1)*cp_ref) / (2.d0*PML_L )*LOG10(1.d0 / PML_Rcoef)

      PML_alpha_beta(1,1) = 1.d0 + alpha_0*((x1 -x1_0) * n1 /PML_L)**afp
      PML_alpha_beta(1,2) = 1.d0 + alpha_0*((x2 -x2_0) * n2 /PML_L)**afp

      PML_alpha_beta(2,1) = beta_0*((x1 -x1_0) * n1 /PML_L)**afp 
      PML_alpha_beta(2,2) = beta_0*((x2 -x2_0) * n2 /PML_L)**afp

      IF (EleType_arg .EQ.1) THEN
        PML_alpha_beta(1:2,1:2) = 0.d0
      end if

      end subroutine PML2D_alpha_beta_function


