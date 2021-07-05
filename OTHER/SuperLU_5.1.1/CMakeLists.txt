# makefile for sparse supernodal LU, implemented in ANSI C

#######################################################################
#  This is the makefile to create a library for SuperLU.
#  The files are organized as follows:
#
#       ALLAUX -- Auxiliary routines called from all precisions of SuperLU
#	SLUSRC -- Single precision real SuperLU routines
#       DLUSRC -- Double precision real SuperLU routines
#	CLUSRC -- Single precision complex SuperLU routines
#       ZLUSRC -- Double precision complex SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single double
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       superlu.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################

project(SUPERLU)

### return machine parameters
#SCAUX 	= smach.o
#DZAUX 	= dmach.o
### SuperLU 

target_sources(ALLAUX PUBLIC
    superlu_timer.c util.c memory.c get_perm_c.c mmd.c 
    sp_coletree.c sp_preorder.c sp_ienv.c relax_snode.c 
    heap_relax_snode.c colamd.c 
    ilu_relax_snode.c ilu_heap_relax_snode.c mark_relax.c 
    mc64ad.c qselect.c input_error.c smach.c dmach.c
)

target_sources(SLUSRC PUBLIC
	sgssv.c sgssvx.c 
	ssp_blas2.c ssp_blas3.c sgscon.c  
	slangs.c sgsequ.c slaqgs.c spivotgrowth.c 
	sgsrfs.c sgstrf.c sgstrs.c scopy_to_ucol.c 
	ssnode_dfs.c ssnode_bmod.c 
	spanel_dfs.c spanel_bmod.c 
	sreadhb.c sreadrb.c sreadtriple.c 
	scolumn_dfs.c scolumn_bmod.c spivotL.c spruneL.c 
	smemory.c sutil.c smyblas2.c 
	sgsisx.c sgsitrf.c sldperm.c 
	ilu_sdrop_row.c ilu_ssnode_dfs.c 
	ilu_scolumn_dfs.c ilu_spanel_dfs.c ilu_scopy_to_ucol.c 
	ilu_spivotL.c sdiagonal.c slacon2.c
)

target_sources(DLUSRC PUBLIC
	dgssv.c dgssvx.c 
	dsp_blas2.c dsp_blas3.c dgscon.c 
	dlangs.c dgsequ.c dlaqgs.c dpivotgrowth.c  
	dgsrfs.c dgstrf.c dgstrs.c dcopy_to_ucol.c 
	dsnode_dfs.c dsnode_bmod.c dpanel_dfs.c dpanel_bmod.c 
	dreadhb.c dreadrb.c dreadtriple.c 
	dcolumn_dfs.c dcolumn_bmod.c dpivotL.c dpruneL.c 
	dmemory.c dutil.c dmyblas2.c 
	dgsisx.c dgsitrf.c dldperm.c 
	ilu_ddrop_row.c ilu_dsnode_dfs.c 
	ilu_dcolumn_dfs.c ilu_dpanel_dfs.c ilu_dcopy_to_ucol.c 
	ilu_dpivotL.c ddiagonal.c dlacon2.c
        ## dgstrsL.c dgstrsU.c

target_sources(CLUSRC PUBLIC
	scomplex.c cgssv.c cgssvx.c csp_blas2.c csp_blas3.c cgscon.c
	clangs.c cgsequ.c claqgs.c cpivotgrowth.c 
	cgsrfs.c cgstrf.c cgstrs.c ccopy_to_ucol.c
	csnode_dfs.c csnode_bmod.c
	cpanel_dfs.c cpanel_bmod.c
	creadhb.c creadrb.c creadtriple.c 
	ccolumn_dfs.c ccolumn_bmod.c cpivotL.c cpruneL.c
	cmemory.c cutil.c cmyblas2.c
	cgsisx.c cgsitrf.c cldperm.c
	ilu_cdrop_row.c ilu_csnode_dfs.c
	ilu_ccolumn_dfs.c ilu_cpanel_dfs.c ilu_ccopy_to_ucol.c
	ilu_cpivotL.c cdiagonal.c clacon2.c scsum1.c icmax1.c
)

target_sources(ZLUSRC PUBLIC 
	dcomplex.c zgssv.c zgssvx.c zsp_blas2.c zsp_blas3.c zgscon.c
	zlangs.c zgsequ.c zlaqgs.c zpivotgrowth.c 
	zgsrfs.c zgstrf.c zgstrs.c zcopy_to_ucol.c
	zsnode_dfs.c zsnode_bmod.c
	zpanel_dfs.c zpanel_bmod.c
	zreadhb.c zreadrb.c zreadtriple.c
	zcolumn_dfs.c zcolumn_bmod.c zpivotL.c zpruneL.c
	zmemory.c zutil.c zmyblas2.c
	zgsisx.c zgsitrf.c zldperm.c
	ilu_zdrop_row.c ilu_zsnode_dfs.c
	ilu_zcolumn_dfs.c ilu_zpanel_dfs.c ilu_zcopy_to_ucol.c
	ilu_zpivotL.c zdiagonal.c zlacon2.c dzsum1.c izmax1.c
)

target_link_libraries(SUPERLU DLUSRC ALLAUX DZLAUX)

set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

export(TARGETS ARPACK
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/SUPERLU.cmake"
)

#all:    double 
#
#single: $(SLUSRC) $(ALLAUX) $(SCAUX)
#	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) $(SLUSRC) $(ALLAUX) $(SCAUX)
#	$(RANLIB) $(SUPERLULIB)
#
#double: $(DLUSRC) $(ALLAUX) $(DZLAUX)
#	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) $(DLUSRC) $(ALLAUX) $(DZLAUX)
#	$(RANLIB) $(SUPERLULIB)
#
#complex: $(CLUSRC) $(ALLAUX) $(SCAUX)
#	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) $(CLUSRC) $(ALLAUX) $(SCAUX)
#	$(RANLIB) $(SUPERLULIB)
#
#complex16: $(ZLUSRC) $(ALLAUX) $(DZLAUX)
#	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) $(ZLUSRC) $(ALLAUX) $(DZLAUX)
#	$(RANLIB) $(SUPERLULIB)
#
#
###################################
## Do not optimize these routines #
###################################
#smach.o: smach.c ; $(CC) -c -fPIC $(NOOPTS) $(CDEFS) $<
#dmach.o: dmach.c ; $(CC) -c -fPIC $(NOOPTS) $(CDEFS) $<
#superlu_timer.o:  superlu_timer.c ; $(CC) -c -fPIC $(NOOPTS) $<
###################################
#
#.c.o:
#	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)
#
#.f.o:
#	$(FORTRAN) $(FFLAGS) -c $<
#
#clean:	
#	rm -f *.o $(SUPERLULIB)

