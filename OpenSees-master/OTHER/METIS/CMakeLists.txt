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

add_library(METIS)

target_sources(METIS PUBLIC
  coarsen.c
  fm.c
  initpart.c
  match.c
  ccgraph.c
  memory.c
  pmetis.c
  pqueue.c
  refine.c
  util.c
  timing.c
  debug.c
  bucketsort.c
  graph.c
  stat.c
  kmetis.c
  kwayrefine.c
  kwayfm.c
  balance.c
  srefine.c
  sfm.c
  separator.c
  mincover.c
  mmd.c
  mesh.c
  meshpart.c
  frename.c
  fortran.c
  myqsort.c
  compress.c
  parmetis.c
  estmem.c
  mpmetis.c
  mcoarsen.c
  mmatch.c
  minitpart.c
  mbalance.c
  mrefine.c
  mutil.c
  mfm.c
  mkmetis.c
  mkwayrefine.c
  mkwayfmh.c
  mrefine2.c
  minitpart2.c
  mbalance2.c
  mfm2.c
  kvmetis.c
  kwayvolrefine.c
  kwayvolfm.c
  subdomains.c
  ometis.c
)

