      Subroutine STEEL(Es,EpsLast,FpsLast,YpTanLast,EpsOld,Fy,Epy, ! Input
     *           EpSH,Epsu,Fpsu,Youngs,SHPower,                     ! Input
     *           Epr,Fpr,Epa,Fpa,Epo,EpoMax,EpsuSh,YoungsUn,        ! Changeable
     *           Power,BFlag,LMR,EprM,FprM,EpaM,FpaM,YpTanM,PowerM, ! Changeable
     *           Eps,Fps,Fs,YpTan,YTan,                             ! Output
     *           OmegFac)                           ! Temporary Variable
C
C==============================================================================
C  Written by L.L. Dodd and J. Restrepo  1991
C
C  This subroutine determines the cyclic stress history of reinforcing steel 
C  given the tensile skeleton curve properties and the strain history as 
C  described by thesis.
C
C  Upgraded by J. Restrepo 1994: Factor.dat is not longer required
C                                OmegaFac = 1.
C				 Units of output strain = units of input strain
C
C  Upgraded by J. Restrepo and M. Schoettler 2011:
C                                Bauschinger curve secant slope > minimum slope
C                                Updated maximum iterations from 5 to 10 in
C                                     subroutines Bausch and LinInterp
C                                Updated input values for "Between Yield 
C                                      Plateaus" and "MAJOR BAUSCHINGER CURVE 
C                                      moving toward yield plateau point"
C
C  Upgraded by J. Restrepo 2011:
C                               Changed the way OMEGA is calculated.
C
C  Upgraded by M. Schoettler 2011:
C                               Accounted for true stress in the yield
C                                     plateau in compression
C                               Corrected strain hardening curve's calculation
C                                     of FpSH.
C                               Corrected a Major Reversal's call to Bausch. when
C                                     within the yield plateau to use Fy*(1+Es),
C                                     instead of Fy.
C                               Corrected a small reversal to follow "reversal
C                                     to strain hardening" instead of "Elastic".
C
C  This subroutine determines the cyclic stress history of reinforcing steel 
C  given the tensile skeleton curve properties and the strain history as 
C  described by Chapter 2 of this thesis.
C
C (C) Copyright 2011, The Regents of the University of California    
C All Rights Reserved.                                               

C The Regents grants permission, without fee and without a written license 
C agreement, for (a) use, reproduction, modification, and distribution of 
C this software and its documentation by educational, research, and non-profit 
C entities for noncommercial purposes only; and (b) use, reproduction and 
C modification of this software by other entities for internal purposes 
C only. The above copyright notice, this paragraph and the following three 
C paragraphs must appear in all copies and modifications of the software 
C and/or documentation.
C
C Permission to incorporate this software into products for commercial 
C distribution may be obtained by contacting the University of California 
C Office of Technology Licensing 
C 2150 Shattuck Avenue #510, 
C Berkeley, CA 94720-1620, 
C (510) 643-7201.

C This software program and documentation are copyrighted by The Regents of
C the University of California. The Regents does not warrant that the operation 
C of the program will be uninterrupted or error-free. The end-user understands 
C that the program was developed for research purposes and is advised not to rely 
C exclusively on the program for any reason.

C IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, 
C INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
C USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF 
C THE POSSIBILITY OF SUCH DAMAGE.  REGENTS GRANTS NO EXPRESS OR IMPLIED LICENSE IN 
C ANY PATENT RIGHTS OF REGENTS BUT HAS IMPLEMENTED AN INDIVIDUAL CONTRIBUTOR LICENSE 
C AGREEMENT FOR THE OPENSEES PROJECT AT THE UNIVERSITY OF CALIFORNIA, BERKELEY TO 
C BENEFIT THE END USER.

C REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
C IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE 
C SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED 
C "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, 
C ENHANCEMENTS, OR MODIFICATIONS. 
C
C J. Retrepo UCSD.
C==============================================================================
C
C
C PASSED VARIABLES
C
      Implicit None
      Integer LMR   ! Last Major Reversal direction. Value of "s" after reversal
      Integer BFlag(2) ! Strain hardening (0) or Bauschinger (1) curve 
C
      Real*8 Epa(2)   ! Strain at end of linear branch (1=tension, 2=compression)
      Real*8 EpaM(2)  ! Major reversal Epa
      Real*8 Epo(2)   ! Maximum "natural" shift (1=compression, 2=tension)
      Real*8 EpoMax   ! The maximum magnitude of Epo(1) and Epo(2)
      Real*8 Eps      ! Natural strain
      Real*8 Epr(2)   ! Reversal strain (1=tension, 2=compression)
      Real*8 EprM(2)  ! Major reversal strain (1=tension, 2=compression)
      Real*8 EpSH     ! Natural coordinate strain hardening strain
      Real*8 EpsLast  ! Natural strain at last increment
      Real*8 EpsOld   ! Natural strain at second to last increment
      Real*8 Epsu     ! Natural coordinate "ultimate" strain
      Real*8 EpsuSh(2)! Shifted "ultimate" strain value (1=tension, 2=compression)
      Real*8 Epy      ! The yield strain, Fy/Youngs (a positive value)
      Real*8 Es       ! Engineering Strain
      Real*8 Fpa(2)   ! Stress at end of linear branch (1=tension, 2=compression)
      Real*8 FpaM(2)  ! Major reversal Fpa
      Real*8 Fpr(2)   ! Reversal stress (1=tension, 2=compression)
      Real*8 FprM(2)  ! Major reversal stress (1=tension, 2=compression)
      Real*8 Fps      ! True coordinates stress
      Real*8 FpsLast  ! True stress at last increment
      Real*8 Fpsu     ! True coordinate "ultimate" stress (slope at ultimate)
      Real*8 Fs       ! Engineering Stress
      Real*8 Fy       ! Yield Stress
      Real*8 OmegFac  ! Multiplication factor for Omega
      Real*8 Power(2) ! Exponent in normalised Bauschinger eq. (1=tens., 2=comp.)
      Real*8 PowerM(2)! Major reversal Power
      Real*8 SHPower  ! Exponent which governs the strain-hardening curve
      Real*8 Youngs   ! Youngs modulus
      Real*8 YTan     ! Tangential modulus
      Real*8 YpTan    ! True coordinates tangential modulus
      Real*8 YpTanM(2)! Tangential modulus at major reversals (1=tens, 2=comp)
      Real*8 YpTanLast! Tangential modulus at last increment
      Real*8 YoungsUn ! Unloading modulus
C
C INTERNAL VARIABLES
C
      Real*8 a        !
      Real*8 C1       ! Temporary constant
      Real*8 C2       ! Temporary constant
      Real*8 Delta    ! Strain change from previous increment
      Real*8 Epp      ! Abs((Epsush(K) - Epa(M))/Epsu)
      Real*8 FNorm    ! Abs(Fpp/Fpt)
      Real*8 Fpp      ! Fpsu*(s*1.-EpsuSh(K)+Epa(M)) - Fpa(M)
      Real*8 Fpt      ! Fpsu*(2-EpsuSh(1)+EpsuSh(2))
      Real*8 FpSH     ! Strain hardening natural strain (Fy*(1+Epsh))
      Real*8 Omega    ! Percent area term for Bauschinger curve
      Integer MaxFlag ! Flag to tell if reversal point is a new max 0-no,1-yes
      Integer S     ! Straining direction: -1 for compressing, 1 for tensioning
      Integer K     ! Index value           2                  1
      Integer M     ! Index value           1                  2
      Integer L     ! Index value: K for LMR*s.NE.-1,  M otherwise
C
C      write(*,*) 'DATA IN:'
C      write(*,*) Es
C      write(*,*) EpsLast
C      write(*,*) FpsLast
C      write(*,*) YpTanLast
C      write(*,*) EpsOld
C      write(*,*) Fy
C      write(*,*) Epy
C      write(*,*) EpSH
C      write(*,*) Epsu
C      write(*,*) Fpsu
C      write(*,*) Youngs
C      write(*,*) SHPower
C      write(*,*) Epr
C      write(*,*) Fpr
C      write(*,*) Epa
C      write(*,*) Fpa
C      write(*,*) Epo
C      write(*,*) EpoMax
C      write(*,*) 'EpsuSh',EpsuSh
C      write(*,*) YoungsUn
C      write(*,*) Power
C      write(*,*) BFlag
C      write(*,*) LMR
C      write(*,*) EprM
C      write(*,*) FprM
C      write(*,*) EpaM
C      write(*,*) FpaM
C      write(*,*) YpTanM
C      write(*,*) PowerM
C       write(*,*) Eps
C       write(*,*) Fps
C       write(*,*) Fs
C       write(*,*) YpTan
C       write(*,*) YTan
C       write(*,*) OmegFac

      Eps = log(1+Es)
      Delta = Eps - EpsLast
      If (Delta.Eq.0) Delta = EpsLast - EpsOld
      If (Delta.Gt.0.0) then                      ! tensioning
         M  = 2
         K  = 1
         S  = 1
      Else                                        ! compressing
         M  =  1
         K  =  2
         S  = -1
      End If
      If (Eps*s.Gt.EpsuSH(K)*s) then
         Fps  = 0.
         Fpsu = 0.
         Write (6,800)
                                                  RETURN
      End If
C
      If (((EpsLast-EpsOld)*Delta).Lt.0.) then    ! Reversal
C
C  =================
C   STRAIN REVERSAL
C  =================
C
        If ((Epo(1).Eq.0.).And.(Epo(2).Eq.0.).And.
     *                          (Abs(EpsLast).Lt.Epy)) then         !Elastic
          Continue
        Else If ((LMR*s.EQ.-1).And.((EpsLast-Epa(K))*s).GE.0.) then !Linear Rng
           Continue
        Else If ((LMR*s.EQ.1).And.(EpsLast-Epr(M))*s.GE.0.) then  !Linear Range
           Continue
        Else
C
          MaxFlag = 0
          If (s*Epo(K).GT.s*(EPsLast-FpsLast/YoungsUn)) then
C
C                                                 ! Max abs strain in direction
C
            MaxFlag  = 1
            Epo(K)    = EpsLast-FpsLast/YoungsUn
            EpsuSh(K) = Epsu*s+Epo(K)
            If (Abs(Epo(K)).Gt.EPoMax) then
              EpoMax   = Abs(Epo(K))                            ! New Max Strain
              YoungsUn = Youngs*(0.82 + 1./(5.55+1000.*EpoMax)) ! Unloading Mod.
            End If
          End If
          LMR    = s
          Epr(M) = EpsLast
          Fpr(M) = FpsLast
          Epa(M) = EpsLast + s*Fy/YoungsUn
          Fpa(M) = FpsLast + s*Fy

          If (((BFlag(K).EQ.0.OR.BFlag(M).EQ.0).And.MaxFlag.EQ.1) .OR.
     +           s*(FprM(K)-Fpr(M)).GT.2*Fy    .OR.
     +       ((Epr(M)-EpaM(K))/(EprM(M)-EpaM(K)).GT.1.0  .And.
     +                 BFlag(K).EQ.1))   then
C
C  ================
C   MAJOR REVERSAL
C  ================
C  Reversal from skeleton curve or more than 2Fy from previous major reversal
C  Reverse to a MAJOR BAUSCHINGER CURVE
C
C  Reversal from skeleton curve or more than 2Fy from previous major reversal
C  Reverse to a MAJOR BAUSCHINGER CURVE
C
            EprM(M)   = Epr(M)
            FprM(M)   = Fpr(M)
            EpaM(M)   = Epa(M)
            FpaM(M)   = Fpa(M)
            YpTanM(M) = YpTanLast
            Epr(K)    = Epr(M)
            EprM(K)   = Epr(M)
            If ((Epo(2)-Epo(1)).LT.(EpSh-Epy)) then    ! Between Yield Plateaus
              Power(K) = 0.35
            Else                                       ! Bauschinger Curve
              BFlag(K) = 1
C Changes were made on 10/05/11 by JR to the way Omega is computed
C   Tests on high elongation stainless steel shows Omega is not a 
C    function of Epsu
C   Values of Omega/OmegFac are limited to the 0.05 to 0.085 range

              Fpt      = Fpsu*(2-EpsuSh(1)+EpsuSh(2))
              Fpp      = Fpsu*(s*1.-EpsuSh(K)+Epa(M)) - Fpa(M)
C              Epp      = Abs((Epsush(K) - Epa(M))/Epsu) ! ORIGINAL
              Epp      = Abs((0.2*s+Epo(K) - Epa(M))/0.2) ! NEW

C              write(*,*)'s',s
C              write(*,*)'Epo(k)',Epo(k)
C              write(*,*)'Epa(M)',Epa(M)
C              write(*,*)'Epp',Epp

              FNorm = Abs(Fpp/Fpt)
              Omega = ((0.001+1.08E-3/(1.043-Epp))/0.18*(FNorm-0.69)
     +                +0.085)
              If (Omega.GT.0.085) Omega = 0.085
              If (Omega.LT.0.05) Omega = 0.05
              Omega = Omega*OmegFac
              Power(K) = 56.689*(Omega-0.077)**2-4.921*(Omega-0.077)+0.1
            End If
            PowerM(K) = Power(K)
          Else If ((Epr(M)-EprM(M))*s.LT.0) then
C
C  ================
C   MINOR REVERSAL
C  ================
C  Reversal from a major Bauschinger curve less than 2Fy from previous 
C  major reversal on 
C  Reverse to a MINOR BAUSCHINGER CURVE
C
            EprM(M)   = Epr(M)
            FprM(M)   = Fpr(M)
            YpTanM(M) = YpTanLast            
            Power(K)  = 0.35
          Else
C
C  Reversal from a minor Bauschinger curve
C
            Power(K)  = 0.35
          End If
        End If
      End If
C
C  ================================
C   GOVERNING STRESS-STRAIN CURVES
C  ================================
C
      If (LMR*s.NE.-1) Then
        L = M
      Else
        L = K
      End If

      If (BFlag(K).eq.1) then
        If (((LMR*s.Eq.-1).And.(Eps*s.Gt.Epr(K)*s)).Or.
     *      ((LMR*s.Eq.1).And.(Eps*s.Gt.Epa(M)*s))) then
C
C  Post Yield-Plateau Bausch1inger Curve
C
          If (s*(Eps-EprM(K)).GT.0.) then
C
C  MAJOR BAUSCHINGER CURVE moving toward ultimate point
C
      Call Bausch1 (Eps,EpaM(M),FpaM(M),EpsuSh(K),Fpsu*s,YoungsUn,     ! Input
     +             Fpsu,PowerM(K),                                    ! Input
     +             Fps,YpTan)                                         ! Output

          Else
C
C  MINOR BAUSCHINGER CURVE moving toward previous minor or major reversal point
C
      Call Bausch1 (Eps,Epa(M),Fpa(M),EprM(K),FprM(K),YoungsUn,        ! Input
     +             YpTanM(K),Power(K),                                ! Input
     +             Fps,YpTan)                                         ! Output
          End If
        Else
          Fps   = Fpr(L) + (Eps-Epr(L))*YoungsUn
          YpTan = YoungsUn
        End If
C
C  Elastic Branch
C
      Else If ((Epo(1).eq.0.).And.(Epo(2).eq.0.).And.
     *        (Abs(Eps).LE.Epy)) then
        Fps   = Eps*Youngs
        YpTan = Youngs

C
C  Skeleton Curve
C
      Else If (s*(Eps-Epo(M)-s*Fy/YoungsUn).GE.-1.E-5) Then
        If (s*(Eps-Epo(K)-s*Epy).LE.(EpSH-Epy)) then
          Fps   = Fy*s*(1+Es)                        ! Yield Plateau
          If (s.LT.0.0) then                         ! MODIFIED BY M. SCHOETTLER 11/8/11
             Fps = -Fy*(1.0/(1+Es))
          End If
          YpTan = Fy*exp(s*Eps)                      ! MODIFIED BY M. SCHOETTLER 11/8/11
        Else
C MODIFIED BY M. SCHOETTLER 11-9-11 FpSH computed incorrectly as Fy*(1+Epsh)
          FpSH  = Fy*exp(Epsh)                          ! Strain Hardening
          C1    = FpSH - Fpsu + Fpsu*(Epsu-EpSH)
          C2    = (Epsu-s*(Eps-Epo(K)))/(Epsu-EpSH)
          Fps   = s*C1*C2**SHPower - Fpsu*(s*Epsu-(Eps-Epo(K))) + s*Fpsu
          YpTan = -SHPower*(C1/(Epsu-EpSH))*C2**(SHPower-1) + FpSU
        End If
C
C Reloading to Strain Hardening
C
      Else If ((EprM(K)-Epo(K))*s.GE.Epsh) then
        If (((LMR*s.Eq.-1).And.(Eps*s.Gt.Epr(K)*s)).Or.
     *      ((LMR*s.Eq.1).And.(Eps*s.Gt.Epa(M)*s))) then
  
      Call Bausch1 (Eps,Epa(M),Fpa(M),EprM(K),FprM(K),YoungsUn,        ! Input
     +             YpTanM(K),Power(K),                                ! Input
     +             Fps,YpTan)                                         ! Output
        Else
          Fps   = Fpr(L) + (Eps-Epr(L))*YoungsUn
          YpTan = YoungsUn
        End If
C
C Between Yield Plateaus
C
      Else

C   Updated 10-5-2011 Epr(K) to EprM(K)
        If (((LMR*s.Eq.-1).And.(Eps*s.Gt.Epr(K)*s)).Or.
     *      ((LMR*s.Eq.1).And.(Eps*s.Gt.Epa(M)*s))) then
          If (s*(Eps-EprM(K)).GT.0.) then
C
C  MAJOR BAUSCHINGER CURVE moving toward yield plateau point
C
C    Updated 10-5-11:  Epa(M) and Fpa(M) to EpaM(M) and FpaM(M)
C    MODIFIED BY M. SCHOETTLER 11-9-11:  Fy*s was updated to Fy*s*(1+Es)
      Call Bausch1 (Eps,EpaM(M),FpaM(M),Epo(M)+s*Fy/YoungsUn,          ! Input
     +             Fy*s*(1+Es),YoungsUn,Fy,Power(K),                   ! Input
     +             Fps,YpTan)                                          ! Output
          Else
C
C  MINOR BAUSCHINGER CURVE moving toward previous minor or major reversal point
C
      Call Bausch1 (Eps,Epa(M),Fpa(M),EprM(K),FprM(K),YoungsUn,         ! Input
     +             YpTanM(K),Power(K),                                 ! Input
     +             Fps,YpTan)                                          ! Output
          End If
        Else
          Fps   = Fpr(L) + (Eps-Epr(L))*YoungsUn
          YpTan = YoungsUn
        End If        
      End If
      Fs   = Fps/(1+Es)
      YTan = (YpTan-fps)*exp(-2*Eps)

C      write(*,*) 'DATA OUT:'
C      write(*,*) Es
C      write(*,*) EpsLast
C      write(*,*) FpsLast
C      write(*,*) YpTanLast
C      write(*,*) EpsOld
C      write(*,*) Fy
C      write(*,*) Epy
C      write(*,*) EpSH
C      write(*,*) Epsu
C      write(*,*) Fpsu
C      write(*,*) Youngs
C      write(*,*) SHPower
C      write(*,*) Epr
C      write(*,*) Fpr
C      write(*,*) Epa
C      write(*,*) Fpa
C      write(*,*) Epo
C      write(*,*) EpoMax
C      write(*,*) EpsuSh
C      write(*,*) YoungsUn
C      write(*,*) 'Power',Power
C      write(*,*) BFlag
C      write(*,*) LMR
C      write(*,*) EprM
C      write(*,*) FprM
C      write(*,*) EpaM
C      write(*,*) FpaM
C      write(*,*) 'YPTan', YpTanM
C      write(*,*) 'PowerM',PowerM
C       write(*,*) Eps
C       write(*,*) Fps
C      write(*,*) Fs
C       write(*,*) YpTan
C       write(*,*) YTan
C       write(*,*) OmegFac

      Return
  800 Format (/' The peak strain has been exceeded, REBAR FRACTURE!'/)
      End
C
      Subroutine Bausch1 (Eps,E1,F1,E2,F2,Slope1,Slope2,Power,         ! Input
     +                   Fps,YpTan)                                   ! Output
C
C This subroutine calculates the stress for a given strain on the Bauschinger
C curve
C
      Integer ITest ! Convergence flag
      Integer I     ! Counter

      Real*8 C1       ! (Fpu-Epu*Slope1)/(Fpu-Epu*F2)
      Real*8 C2       ! Eps*(Slope1-F2)/(Fpu-F2*Epu)
      Real*8 C3       ! 1 - Eppn
      Real*8 C4       ! 1 - C3*C3 = 1-(1-Eppn)^2
      Real*8 C5       ! C4**Power-C1*Eppn-C2 (function for which Eppn is a root)
      Real*8 E1       ! Initial strain on Bauschinger curve
      Real*8 Eppn     ! Normalised Strain in the Bauschinger Branch
      Real*8 Eppn2    ! New estimate of Eppn
      Real*8 Eps      ! Natural strain
      Real*8 E2       ! Final strain on Bauschinger curve
      Real*8 Epu      ! Strain from linear branch to ultimate (E2 - E1)
      Real*8 F1       ! Initial stress on Bauschinger curve
      Real*8 Fps      ! True coordinates stress
      Real*8 Fpu      ! Stress from linear branch to ultimate (F2 - F1)
      Real*8 F2       ! Final stress on Bauschinger curve
      Real*8 Power    ! Exponent in the Normalised Bauschinger equation
      Real*8 Slope1   ! Initial Slope
      Real*8 Slope2   ! Final Slope
      Real*8 Slope3   ! Slope of the line passing through (E1,F1) and (E2,F3)
      Real*8 YpTan    ! Tangential modulus
C
      Slope3  = 4.0/5*(F2-F1)/(E2-E1) ! Updated 10-5-2011: Added check for Slope2
      Slope2  = MIN(Slope2, Slope3)
      fpu  = F2 - F1
      Epu  = E2 - E1
      C1   = (Fpu-Epu*Slope1)/(Fpu-Epu*Slope2)
      C2   = (Eps-E1)*(Slope1-Slope2)/(Fpu-Slope2*Epu)
      Eppn = (Eps-E1)/(E2-E1)
      ITest = 1
      I     = 0	

      Do While ((ITest.Eq.1).And.(I.LT.10))
         I     = I+1
         C3    = 1-Eppn
         C4    = 1-C3*C3
         C5    = C4**Power-C1*Eppn-C2
         Eppn2 = Eppn - C5/(2*Power*C4**(Power - 1)*C3-C1) ! Newton-Raphson

C	 If (Eppn2.LT.0) then
C	     Eppn2=0
C	 End If
         If (Eppn2.GT.0.02) then
            Eppn = Eppn2
            If (Abs(C5).LE.0.001) ITest = 0
         Else
            Call LinInterp(Eppn,C1,C2,Power)
            ITest = 0
         End If
      End Do

      Fps  = Eppn*(FPu-Epu*Slope1)+(Eps-E1)*Slope1+F1

      If (Eppn.LT.0.0001 .OR. (Slope1-Slope2)/slope1.LT.0.01) then
        YpTan = Slope1
      Else
        YpTan = 2*Power*(1-(1-Eppn)**2)**(Power-1.)*(1-Eppn) !Normal coordinates
        YpTan =YpTan*(Fpu-Slope2*Epu)/((Epu*Slope1-Fpu)/(Slope1-Slope2))
        YpTan = YpTan*Slope1/(YpTan+Slope1) + Slope2
      End If

      Return
      End
C
      Subroutine LinInterp (Eppn,C1,C2,Power)
C
C  Calculate Eppn using an iterative linear inperpolation
C         
      Integer ITest ! Convergence flag
      Integer I     ! counter
C
      Real*8 C1       ! (Fpu-Epu*Slope1)/(Fpu-Epu*Slope2)
      Real*8 C2       ! Eps*(Slope1-Slope2)/(Fpu-Slope2*Epu)
      Real*8 C3       ! 1 - Eppn
      Real*8 C4       ! 1 - C3*C3 = 1-(1-Eppn)^2
      Real*8 C5       ! C4**Power-C1*Eppn-C2 (function for which Eppn is a root)
      Real*8 C5L      ! Lower bound of C5
      Real*8 C5U      ! Upper bound of C5
      Real*8 Eppn     ! Normalised Strain in the Bauschinger Branch
      Real*8 EppnL    ! Lower bound of Eppn
      Real*8 EppnU    ! Upper bound of Eppn
      Real*8 Eps      ! Natural strain
      Real*8 Power    ! Exponent in the Normalised Bauschinger equation

      EppnU = Eppn
       C3    = 1-EppnU
       C4    = 1-C3*C3
       C5U   = C4**Power-C1*EppnU-C2
       EppnL = 0
       C5L   = -C2

C       ITest = 1
C       Do While (ITest.Eq.1)
       Do I=1,10
         Eppn = EppnL-C5L*(EppnU-EppnL)/(C5U-C5L)
            If (Eppn.Lt.0) Eppn = 0
          C3    = 1-Eppn
          C4    = 1-C3*C3
          C5    = C4**Power-C1*Eppn-C2
C          If (Abs(C5).LE.0.03) then
C             ITest = 0
C          Else 
             If (C5.Gt.0) then

                EppnU = Eppn
                C5U   = C5
             Else
                EppnL = Eppn
                C5L   = C5
             End If
C          End If
       End Do
       Return
       End

