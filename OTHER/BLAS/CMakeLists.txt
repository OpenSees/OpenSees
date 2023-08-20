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

add_library(BLAS)

target_sources(BLAS PUBLIC
    dasum.f
    daxpy.f
    dcopy.f
    ddot.f
    dgbmv.f
    dgeev.f
    dgemm.f
    dgemv.f
    dger.f
    dnrm2.f
    dpotrf.f
    drot.f
    drotg.f
    dscal.f
    dswap.f
    dsymm.f
    dsymv.f
    dsyr.f
    dsyr2.f
    dsyr2k.f
    dsyrk.f
    dtbsv.f
    dtrmm.f
    dtrmv.f
    dtrsm.f
    dtrsv.f
    dtrtrs.f
    dzasum.f
    dznrm2.f
    icamax.f
    idamax.f
    isamax.f
    izamax.f
)