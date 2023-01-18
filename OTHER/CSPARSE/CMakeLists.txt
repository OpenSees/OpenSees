#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#------------------------------------------------------------------------------
#                       CSPARSE version 3.1.0, Jun 1 2012
#                Timothy A. Davis, http://www.suitesparse.com
#==============================================================================
# Claudio Perez

add_library(CSPARSE STATIC
    cs_add.c
    cs_amd.c
    cs_chol.c
    cs_cholsol.c
    cs_counts.c
    cs_cumsum.c
    cs_droptol.c
    cs_dropzeros.c
    cs_dupl.c
    cs_entry.c
    cs_etree.c
    cs_fkeep.c
    cs_gaxpy.c
    cs_happly.c
    cs_house.c
    cs_ipvec.c
    cs_lsolve.c
    cs_ltsolve.c
    cs_lu.c
    cs_lusol.c
    cs_util.c
    cs_multiply.c
    cs_permute.c
    cs_pinv.c
    cs_post.c
    cs_pvec.c
    cs_qr.c
    cs_qrsol.c
    cs_scatter.c
    cs_schol.c
    cs_sqr.c
    cs_symperm.c
    cs_tdfs.c
    cs_malloc.c
    cs_transpose.c
    cs_compress.c
    cs_usolve.c
    cs_utsolve.c
    cs_scc.c
    cs_maxtrans.c
    cs_dmperm.c
    cs_updown.c
    cs_print.c
    cs_norm.c
    cs_load.c
    cs_dfs.c
    cs_reach.c
    cs_spsolve.c
    cs_ereach.c
    cs_leaf.c
    cs_randperm.c
)

