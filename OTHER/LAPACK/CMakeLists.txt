#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#------------------------------------------------------------------------------
# Claudio Perez

#project(LAPACK)

add_library(LAPACK)

target_sources(LAPACK PUBLIC
  lsame.f
  lsamen.f
  xerbla.f
  ilaenv.f
  xlaenv.f
  ieeeck.f

  dgeqr2.f
  dlabad.f
  dlacon.f
  dlacpy.f
  dladiv.f
  dlae2.f
  dlaev2.f
  dlaexc.f
  dlagtm.f
  dlahqr.f
  dlaln2.f
  dlamch.f
  dlange.f
  dlanhs.f
  dlanst.f
  dlanv2.f
  dlaptm.f
  dlapy2.f
  dlapy3.f
  dlaran.f
  dlarf.f
  dlarfg.f
  dlarfx.f
  dlarnd.f
  dlarnv.f
  dlartg.f
  dlaruv.f
  dlascl.f
  dlaset.f
  dlasr.f
  dlasrt.f
  dlassq.f
  dlasy2.f
  dorm2r.f
  dsteqr.f
  dtrevc.f
  dtrexc.f
  dtrsen.f
  dtrsyl.f
  dgbtf2.f
  dgbtrf.f
  dgbtrs.f
  dgttrf.f
  dgttrs.f
  dpttrf.f
  dpttrs.f
  dlaswp.f
  dgbsv.f
  dpbtrs.f
  dgesv.f
  dgetrs.f
  dpbsv.f
  dgetrf.f
  dgetf2.f
  dpbtrf.f
  dpotf2.f
  dpbtf2.f
  dgels.f
  dgeqrf.f
  dgelqf.f
  dormqr.f
  dormlq.f
  dgelq2.f
  dlarft.f
  dlarfb.f
  dorml2.f
  dlaebz.f
  dlagtf.f
  dlagts.f
  dlansb.f
  dlar2v.f
  dlargv.f
  dlartv.f
  dsbevx.f
  dsbtrd.f
  dstebz.f
  dstein.f
  dsterf.f
  dgetri.f
  dtrtri.f
  dtrti2.f
  dgebak.f
  dgebal.f
  dgeev.f
  dgehd2.f
  dgehrd.f
  dhseqr.f
  dlahrd.f  
  dorg2r.f
  dorghr.f
  dorgqr.f
  dpotrf.f
  dpotrs.f
  dtrtrs.f
  dggbak.f
  dggbal.f
  dggev.f
  dgghrd.f
  dhgeqz.f
  dlag2.f
  dlasv2.f
  dtgevc.f
  dgesvxx.f
  dla_gerfsx_extended.f
  dlascl2.f
  dlaqge.f
  dgeequb.f
  dgerfsx.f
  chla_transtype.f
  dla_gercond.f
  dla_wwaddw.f
  ilaprec.f
  ilatrans.f
  dgecon.f
  dla_geamv.f
  dlacn2.f
  dla_lin_berr.f
  drscl.f
  dlatrs.f
  dla_rpvgrw.f
  dla_gerpvgrw.f
)  

add_library(LAPACK_C)
target_sources(LAPACK_C PUBLIC
  BLAS_dgemv2_x.c 
  BLAS_dgemv2_x-f2c.c 
  BLAS_dgemv_x.c 
  BLAS_dgemv_x-f2c.c 
  BLAS_error.c 
  BLAS_error.c
)

target_link_libraries(LAPACK PUBLIC LAPACK_C BLAS)



