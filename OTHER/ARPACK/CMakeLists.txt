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


add_library(ARPACK STATIC
dgetv0.f
dlaqrb.f
dmout.f
dnaitr.f
dnapps.f
dnaup2.f
dnaupd.f
dnaupe.f
dnconv.f
dneigh.f
dneupd.f
dngets.f
dsaitr.f
dsapps.f
dsaup2.f
dsaupd.f
dsconv.f
dseigt.f
dsesrt.f
dseupd.f
dsgets.f
dsortc.f
dsortr.f
dstatn.f
dstats.f
dstqrb.f
dvout.f
dvout1.f
icnteq.f
icopy.f
iset.f
iswap.f
ivout.f
ivout1.f
second.f
)




