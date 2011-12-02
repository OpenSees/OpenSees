c **********************************************************************
c  HEADER FILE:  infel00.h
c **********************************************************************
c  DRAIN-2DX ROTATIONAL SPRING ELEMENT - Bilinear Return Map
c  MHS
c  June 2001
c ---------------------------------------------------------------------
c  PURPOSE
c     This header file establishes the /infel/ common block
c     for the rotational spring element.
c ---------------------------------------------------------------------

      common /infel/ E,sigY,Hiso,Hkin,
     1               ep,alpha,kappa,
     2               epsP,sigP,tangP,
     3               tang

c ---------------------------------------------------------------------
c  VARIABLES
c
c     E           - elastic modulus
c     sigY        - yield stress
c     Hiso        - isotropic hardening modulus
c     Hkin        - kinematic hardening modulus
c     ep          - effective plastic strain
c     alpha       - internal hardening variables
c     kappa       - back stress
c     epsP        - committed strain
c     sigP        - committed stress
c     tangP       - committed tangent
c     tang        - trial tangent
c ---------------------------------------------------------------------















