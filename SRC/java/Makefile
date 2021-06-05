include ../../Makefile.def

SHARED = true

include ../Makefile 

OBJS       = OpenSeesEvaluator.o OpenSeesEvaluatorKepler.o

FE_LIBRARY      = $(HOME)/lib/libOpenSees.so 

all: $(OBJS)

java: $(OBJS) 
	g++ -shared $(ELE_LIBS) $(MATERIAL_LIBS) $(SECTION_LIBS) \
	$(DATABASE_LIBS) $(FEDEAS_LIBS) $(MODEL_BUILDER_LIBS) \
	$(DOMAIN_LIBS) $(MATMOD_BJ_LIBS) $(NDARRAY_BJ_LIBS) \
	$(ELEMENT_BJ_LIBS) $(FINITE_DEFORMATION_BJ_LIBS) \
	$(PY_SJB_RWB_BJ_LIBS) $(ANALYSIS_LIBS) $(SysOfEqn_LIBS) \
	$(GRAPH_LIBS) $(TAGGED_LIBS) $(UTILITY_LIBS) $(ACTOR_LIBS) \
	$(RECORDER_LIBS) $(LAW_LIBS) $(GGD_LIBS) $(RENDERER_LIBS) \
	$(REMO_LIBS) $(GLF_LIBS) $(ZHY_LIBS) $(ERRORHANDLER_LIBS) \
	$(DRAIN_LIBS)  $(SNAP_LIBS) $(FMK_LIBS) $(KJE_LIBS) $(MATRIX_LIBS) \
	$(ARPACK_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) \
	 /usr/local/lib/libg2c.a /usr/local/lib/libgfortran.a \
	$(SUPERLU_LIBS) $(UMFPACK_LIBS) $(METIS_LIBS) \
	OpenSeesEvaluator.o $(TCL_LIBS) -o libOpenSeesEvaluator.so

kepler: $(OBJS) 
	g++ -shared $(ELE_LIBS) $(MATERIAL_LIBS) $(SECTION_LIBS) \
	$(DATABASE_LIBS) $(FEDEAS_LIBS) $(MODEL_BUILDER_LIBS) \
	$(DOMAIN_LIBS) $(MATMOD_BJ_LIBS) $(NDARRAY_BJ_LIBS) \
	$(ELEMENT_BJ_LIBS) $(FINITE_DEFORMATION_BJ_LIBS) \
	$(PY_SJB_RWB_BJ_LIBS) $(ANALYSIS_LIBS) $(SysOfEqn_LIBS) \
	$(GRAPH_LIBS) $(TAGGED_LIBS) $(UTILITY_LIBS) $(ACTOR_LIBS) \
	$(RECORDER_LIBS) $(LAW_LIBS) $(GGD_LIBS) $(RENDERER_LIBS) \
	$(REMO_LIBS) $(GLF_LIBS) $(ZHY_LIBS) $(ERRORHANDLER_LIBS) \
	$(DRAIN_LIBS)  $(SNAP_LIBS) $(FMK_LIBS) $(KJE_LIBS) $(MATRIX_LIBS) \
	$(ARPACK_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)  \
	/usr/local/lib/libg2c.a /usr/local/lib/libgfortran.a \
	$(SUPERLU_LIBS) $(UMFPACK_LIBS) $(METIS_LIBS) OpenSeesEvaluator.o \
	$(TCL_LIBS) -o /usr/local/kepler/lib/libOpenSeesEvaluator.so



TCL_LIBS = $(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o \
	/usr/lib/libtcl.so

ARPACK_LIBS =  $(FE)/../OTHER/ARPACK/dgetv0.o \
	$(FE)/../OTHER/ARPACK/dlaqrb.o \
	$(FE)/../OTHER/ARPACK/dstqrb.o \
	$(FE)/../OTHER/ARPACK/dsortc.o \
	$(FE)/../OTHER/ARPACK/dsortr.o \
	$(FE)/../OTHER/ARPACK/dstatn.o \
	$(FE)/../OTHER/ARPACK/dstats.o \
	$(FE)/../OTHER/ARPACK/dnaitr.o \
	$(FE)/../OTHER/ARPACK/dnapps.o \
	$(FE)/../OTHER/ARPACK/dnaup2.o \
	$(FE)/../OTHER/ARPACK/dnaupd.o \
	$(FE)/../OTHER/ARPACK/dnconv.o \
	$(FE)/../OTHER/ARPACK/dneigh.o \
	$(FE)/../OTHER/ARPACK/dngets.o \
	$(FE)/../OTHER/ARPACK/dsaitr.o \
	$(FE)/../OTHER/ARPACK/dsapps.o \
	$(FE)/../OTHER/ARPACK/dsaup2.o \
	$(FE)/../OTHER/ARPACK/dsaupd.o \
	$(FE)/../OTHER/ARPACK/dsconv.o \
	$(FE)/../OTHER/ARPACK/dseigt.o \
	$(FE)/../OTHER/ARPACK/dsgets.o \
	$(FE)/../OTHER/ARPACK/dneupd.o \
	$(FE)/../OTHER/ARPACK/dseupd.o \
	$(FE)/../OTHER/ARPACK/dsesrt.o \
	$(FE)/../OTHER/ARPACK/dvout.o  \
	$(FE)/../OTHER/ARPACK/dmout.o \
	$(FE)/../OTHER/ARPACK/icnteq.o \
	$(FE)/../OTHER/ARPACK/icopy.o \
	$(FE)/../OTHER/ARPACK/iset.o \
	$(FE)/../OTHER/ARPACK/iswap.o \
	$(FE)/../OTHER/ARPACK/ivout.o \
	$(FE)/../OTHER/ARPACK/second.o


LAPACK_LIBS =  $(FE)/../OTHER/LAPACK/dgeqr2.o \
	$(FE)/../OTHER/LAPACK/dlabad.o \
	$(FE)/../OTHER/LAPACK/dlacpy.o \
	$(FE)/../OTHER/LAPACK/dladiv.o \
	$(FE)/../OTHER/LAPACK/dlae2.o \
	$(FE)/../OTHER/LAPACK/dlaev2.o \
	$(FE)/../OTHER/LAPACK/dlaexc.o \
	$(FE)/../OTHER/LAPACK/dlagtm.o \
	$(FE)/../OTHER/LAPACK/dlahqr.o \
	$(FE)/../OTHER/LAPACK/dlaln2.o \
	$(FE)/../OTHER/LAPACK/dlamch.o \
	$(FE)/../OTHER/LAPACK/dlange.o \
	$(FE)/../OTHER/LAPACK/dlanhs.o \
	$(FE)/../OTHER/LAPACK/dlanst.o \
	$(FE)/../OTHER/LAPACK/dlanv2.o \
	$(FE)/../OTHER/LAPACK/dlaptm.o \
	$(FE)/../OTHER/LAPACK/dlapy2.o \
	$(FE)/../OTHER/LAPACK/dlapy3.o \
	$(FE)/../OTHER/LAPACK/dlaran.o \
	$(FE)/../OTHER/LAPACK/dlarf.o \
	$(FE)/../OTHER/LAPACK/dlarfg.o \
	$(FE)/../OTHER/LAPACK/dlarfx.o \
	$(FE)/../OTHER/LAPACK/dlarnd.o \
	$(FE)/../OTHER/LAPACK/dlarnv.o \
	$(FE)/../OTHER/LAPACK/dlartg.o \
	$(FE)/../OTHER/LAPACK/dlaruv.o \
	$(FE)/../OTHER/LAPACK/dlascl.o \
	$(FE)/../OTHER/LAPACK/dlaset.o \
	$(FE)/../OTHER/LAPACK/dlasr.o \
	$(FE)/../OTHER/LAPACK/dlasrt.o \
	$(FE)/../OTHER/LAPACK/dlassq.o \
	$(FE)/../OTHER/LAPACK/dlasy2.o \
	$(FE)/../OTHER/LAPACK/dorm2r.o \
	$(FE)/../OTHER/LAPACK/dsteqr.o \
	$(FE)/../OTHER/LAPACK/dtrevc.o \
	$(FE)/../OTHER/LAPACK/dtrexc.o \
	$(FE)/../OTHER/LAPACK/dtrsen.o \
	$(FE)/../OTHER/LAPACK/dtrsyl.o \
	$(FE)/../OTHER/LAPACK/dgbtf2.o \
	$(FE)/../OTHER/LAPACK/dgbtrf.o \
	$(FE)/../OTHER/LAPACK/dgbtrs.o \
	$(FE)/../OTHER/LAPACK/dgttrf.o \
	$(FE)/../OTHER/LAPACK/dgttrs.o \
	$(FE)/../OTHER/LAPACK/dpttrf.o \
	$(FE)/../OTHER/LAPACK/dpttrs.o \
	$(FE)/../OTHER/LAPACK/dlaswp.o \
	$(FE)/../OTHER/LAPACK/dgbsv.o \
	$(FE)/../OTHER/LAPACK/dpbtrs.o \
	$(FE)/../OTHER/LAPACK/dgesv.o \
	$(FE)/../OTHER/LAPACK/dgetrs.o \
	$(FE)/../OTHER/LAPACK/dpbsv.o \
	$(FE)/../OTHER/LAPACK/dgetrf.o \
	$(FE)/../OTHER/LAPACK/dgetf2.o \
	$(FE)/../OTHER/LAPACK/dpbtrf.o \
	$(FE)/../OTHER/LAPACK/dpotf2.o \
	$(FE)/../OTHER/LAPACK/dpbtf2.o \
	$(FE)/../OTHER/LAPACK/dgels.o \
	$(FE)/../OTHER/LAPACK/dgeqrf.o \
	$(FE)/../OTHER/LAPACK/dgelqf.o \
	$(FE)/../OTHER/LAPACK/dormqr.o \
	$(FE)/../OTHER/LAPACK/dormlq.o \
	$(FE)/../OTHER/LAPACK/dgelq2.o \
	$(FE)/../OTHER/LAPACK/dlarft.o \
	$(FE)/../OTHER/LAPACK/dlarfb.o \
	$(FE)/../OTHER/LAPACK/dorml2.o \
	$(FE)/../OTHER/LAPACK/dlaebz.o \
	$(FE)/../OTHER/LAPACK/dlagtf.o \
	$(FE)/../OTHER/LAPACK/dlagts.o \
	$(FE)/../OTHER/LAPACK/dlansb.o \
	$(FE)/../OTHER/LAPACK/dlar2v.o \
	$(FE)/../OTHER/LAPACK/dlargv.o \
	$(FE)/../OTHER/LAPACK/dlartv.o \
	$(FE)/../OTHER/LAPACK/dsbevx.o \
	$(FE)/../OTHER/LAPACK/dsbtrd.o \
	$(FE)/../OTHER/LAPACK/dstebz.o \
	$(FE)/../OTHER/LAPACK/dstein.o \
	$(FE)/../OTHER/LAPACK/dsterf.o \
	$(FE)/../OTHER/LAPACK/dgetri.o \
	$(FE)/../OTHER/LAPACK/dtrtri.o \
	$(FE)/../OTHER/LAPACK/dtrti2.o \
	$(FE)/../OTHER/LAPACK/dgebak.o \
	$(FE)/../OTHER/LAPACK/dgebal.o \
	$(FE)/../OTHER/LAPACK/dgeev.o \
	$(FE)/../OTHER/LAPACK/dgehd2.o \
	$(FE)/../OTHER/LAPACK/dgehrd.o \
	$(FE)/../OTHER/LAPACK/dhseqr.o \
	$(FE)/../OTHER/LAPACK/dlahrd.o \
	$(FE)/../OTHER/LAPACK/dorg2r.o \
	$(FE)/../OTHER/LAPACK/dorghr.o \
	$(FE)/../OTHER/LAPACK/dorgqr.o \
	$(FE)/../OTHER/LAPACK/dpotrf.o \
	$(FE)/../OTHER/LAPACK/dpotrs.o \
	$(FE)/../OTHER/LAPACK/ilaenv.o \
	$(FE)/../OTHER/LAPACK/ieeeck.o

BLAS_LIBS = $(FE)/../OTHER/BLAS/idamax.o \
	$(FE)/../OTHER/BLAS/dasum.o \
	$(FE)/../OTHER/BLAS/daxpy.o \
	$(FE)/../OTHER/BLAS/dcopy.o \
	$(FE)/../OTHER/BLAS/ddot.o \
	$(FE)/../OTHER/BLAS/dnrm2.o \
	$(FE)/../OTHER/BLAS/drot.o \
	$(FE)/../OTHER/BLAS/drotg.o \
	$(FE)/../OTHER/BLAS/dscal.o \
	$(FE)/../OTHER/BLAS/dswap.o \
	$(FE)/../OTHER/BLAS/dgemv.o \
	$(FE)/../OTHER/BLAS/dger.o \
	$(FE)/../OTHER/BLAS/dsymv.o \
	$(FE)/../OTHER/BLAS/dsyr2.o \
	$(FE)/../OTHER/BLAS/dtrmm.o \
	$(FE)/../OTHER/BLAS/dgbmv.o \
	$(FE)/../OTHER/BLAS/dgemm.o \
	$(FE)/../OTHER/BLAS/dtrsm.o \
	$(FE)/../OTHER/BLAS/dtbsv.o \
	$(FE)/../OTHER/BLAS/dsyr.o \
	$(FE)/../OTHER/BLAS/dtrsv.o \
	$(FE)/../OTHER/BLAS/dsyrk.o \
	$(FE)/../OTHER/BLAS/dtrmv.o

SUPERLU_LIBS = $(FE)/../OTHER/SuperLU_3.0/SRC/superlu_timer.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/lsame.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/util.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/memory.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/get_perm_c.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/mmd.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/sp_coletree.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/sp_preorder.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/sp_ienv.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/relax_snode.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/heap_relax_snode.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/xerbla.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/colamd.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dgssv.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dgssvx.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dsp_blas2.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dsp_blas3.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dgscon.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dlacon.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dlangs.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dgsequ.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dlaqgs.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dpivotgrowth.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dgsrfs.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dgstrf.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dgstrs.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dcopy_to_ucol.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dsnode_dfs.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dsnode_bmod.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dpanel_dfs.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dpanel_bmod.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dreadhb.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dcolumn_dfs.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dcolumn_bmod.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dpivotL.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dpruneL.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dmemory.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dutil.o \
	$(FE)/../OTHER/SuperLU_3.0/SRC/dmyblas2.o


UMFPACK_LIBS = $(FE)/../OTHER/UMFPACK/umd2fb.o \
	$(FE)/../OTHER/UMFPACK/umd2co.o \
	$(FE)/../OTHER/UMFPACK/umd2fa.o \
	$(FE)/../OTHER/UMFPACK/umd2f0.o \
	$(FE)/../OTHER/UMFPACK/umd2f1.o \
	$(FE)/../OTHER/UMFPACK/umd2f2.o \
	$(FE)/../OTHER/UMFPACK/umd2fg.o \
	$(FE)/../OTHER/UMFPACK/umd2in.o \
	$(FE)/../OTHER/UMFPACK/umd2of.o \
	$(FE)/../OTHER/UMFPACK/umd2s2.o \
	$(FE)/../OTHER/UMFPACK/umd2sl.o \
	$(FE)/../OTHER/UMFPACK/umd2so.o \
	$(FE)/../OTHER/UMFPACK/umd2su.o \
	$(FE)/../OTHER/UMFPACK/umd2er.o \
	$(FE)/../OTHER/UMFPACK/umd2p1.o \
	$(FE)/../OTHER/UMFPACK/umd2p2.o \
	$(FE)/../OTHER/UMFPACK/umd2lt.o \
	$(FE)/../OTHER/UMFPACK/umd2ut.o \
	$(FE)/../OTHER/UMFPACK/umd2rf.o \
	$(FE)/../OTHER/UMFPACK/umd2ra.o \
	$(FE)/../OTHER/UMFPACK/umd2r0.o \
	$(FE)/../OTHER/UMFPACK/umd2r2.o \
	$(FE)/../OTHER/UMFPACK/umd2rg.o \
	$(FE)/../OTHER/UMFPACK/umd21i.o \
	$(FE)/../OTHER/UMFPACK/mc21b.o \
	$(FE)/../OTHER/UMFPACK/mc13e.o

METIS_LIBS	= $(FE)/../OTHER/METIS/io.o \
	$(FE)/../OTHER/METIS/util.o \
	$(FE)/../OTHER/METIS/match.o \
	$(FE)/../OTHER/METIS/coarsen.o \
	$(FE)/../OTHER/METIS/memory.o \
	$(FE)/../OTHER/METIS/fm.o  \
	$(FE)/../OTHER/METIS/GKlib.o \
	$(FE)/../OTHER/METIS/htable.o \
	$(FE)/../OTHER/METIS/list.o \
	$(FE)/../OTHER/METIS/refine.o \
	$(FE)/../OTHER/METIS/bucketlist.o \
	$(FE)/../OTHER/METIS/stat.o \
	$(FE)/../OTHER/METIS/mlevelpart.o \
	$(FE)/../OTHER/METIS/separator.o \
	$(FE)/../OTHER/METIS/mlevelorder.o \
	$(FE)/../OTHER/METIS/smbfactor.o \
	$(FE)/../OTHER/METIS/smbfct.o \
	$(FE)/../OTHER/METIS/mincover.o \
	$(FE)/../OTHER/METIS/initpart.o \
	$(FE)/../OTHER/METIS/lanczos.o \
	$(FE)/../OTHER/METIS/mmd.o \
	$(FE)/../OTHER/METIS/balpart.o \
	$(FE)/../OTHER/METIS/kwaypart.o \
	$(FE)/../OTHER/METIS/kwayfm.o \
	$(FE)/../OTHER/METIS/kwaygreedy.o \
	$(FE)/../OTHER/METIS/entrypoint.o 


# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.





