! Minimal isotropic linear-elastic material in Abaqus UMAT format.
!
! Used by tests/test_umat_material.py to exercise the generic
! `nDMaterial UMAT` wrapper (SRC/material/nD/UmatMaterial.cpp): the test
! compiles this file into a shared library at run time, binds it with
! `-library`, and checks that a brick element driven by this UMAT matches
! the same brick driven by the built-in ElasticIsotropic material.
!
!   PROPS(1) = Young's modulus E
!   PROPS(2) = Poisson's ratio nu
!   NSTATV   = 0 (stateless)
!
! Three-dimensional only (NTENS = 6, NDI = 3, NSHR = 3), tension positive,
! engineering shear strains, Abaqus component order 11 22 33 12 13 23.
! Written in the standard UMAT argument list so the same file also builds
! under Abaqus for a side-by-side comparison.

subroutine umat_elastic(stress, statev, ddsdde, sse, spd, scd, &
                        rpl, ddsddt, drplde, drpldt, &
                        stran, dstran, time, dtime, temp, dtemp, &
                        predef, dpred, cmname, ndi, nshr, ntens, nstatv, &
                        props, nprops, coords, drot, pnewdt, celent, &
                        dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
  implicit none
  character(len=80) :: cmname
  integer :: ndi, nshr, ntens, nstatv, nprops, noel, npt, layer, kspt, kinc
  integer :: kstep(4)
  double precision :: stress(ntens), statev(nstatv), ddsdde(ntens, ntens)
  double precision :: sse, spd, scd, rpl, drpldt, dtime, temp, dtemp, pnewdt, celent
  double precision :: ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens)
  double precision :: time(2), predef(1), dpred(1), props(nprops)
  double precision :: coords(3), drot(3, 3), dfgrd0(3, 3), dfgrd1(3, 3)

  double precision :: e, nu, lam, mu
  integer :: i, j

  e  = props(1)
  nu = props(2)
  lam = e * nu / ((1.0d0 + nu) * (1.0d0 - 2.0d0 * nu))
  mu  = e / (2.0d0 * (1.0d0 + nu))

  ! Elastic tangent: lambda + 2mu on the normal diagonal, lambda off it,
  ! mu on the shear diagonal (engineering shear strains).
  ddsdde = 0.0d0
  do i = 1, ndi
    do j = 1, ndi
      ddsdde(i, j) = lam
    end do
    ddsdde(i, i) = lam + 2.0d0 * mu
  end do
  do i = ndi + 1, ntens
    ddsdde(i, i) = mu
  end do

  ! Linear elastic: stress increment is tangent times strain increment.
  do i = 1, ntens
    do j = 1, ntens
      stress(i) = stress(i) + ddsdde(i, j) * dstran(j)
    end do
  end do

  ! Elastic strain energy, for completeness (not used by the wrapper).
  sse = 0.0d0
  do i = 1, ntens
    sse = sse + 0.5d0 * stress(i) * (stran(i) + dstran(i))
  end do

  pnewdt = 1.0d0
end subroutine umat_elastic
