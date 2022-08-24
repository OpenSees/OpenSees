/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: protos.h
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99sep14
/
/ This file contains the prototypes of all non-static functions
/
******************************************************************************/

/* functions in lib/greg_pord.h */
PORD_INT 	greg_pord(PORD_INT, PORD_INT, PORD_INT *, PORD_INT *, PORD_INT *, PORD_INT *, PORD_INT *);

/* functions in lib/graph.c */
graph_t*	newGraph(PORD_INT, PORD_INT);
void		freeGraph(graph_t*);
void		printGraph(graph_t*);
void		randomizeGraph(graph_t*);
graph_t*	setupSubgraph(graph_t*, PORD_INT*, PORD_INT, PORD_INT*);
graph_t*	setupGraphFromMtx(inputMtx_t*);
graph_t*	setupGridGraph(PORD_INT, PORD_INT, PORD_INT);
PORD_INT       	connectedComponents(graph_t*);
graph_t*	compressGraph(graph_t*, PORD_INT*);

/* functions in lib/gbisect.c */
gbisect_t*      newGbisect(graph_t*);
void            freeGbisect(gbisect_t*);
void		printGbisect(gbisect_t*);
void            checkSeparator(gbisect_t*);
void            constructSeparator(gbisect_t*, options_t*, timings_t*);
PORD_INT       	smoothBy2Layers(gbisect_t*, PORD_INT*, PORD_INT*, PORD_INT, PORD_INT);
void            smoothSeparator(gbisect_t*, options_t*);

/* functions in lib/ddcreate.c */
domdec_t*	newDomainDecomposition(PORD_INT, PORD_INT);
void		freeDomainDecomposition(domdec_t*);
void		printDomainDecomposition(domdec_t*);
void		checkDomainDecomposition(domdec_t*);
void		buildInitialDomains(graph_t*, PORD_INT*, PORD_INT*, PORD_INT*);
void		mergeMultisecs(graph_t *G, PORD_INT*, PORD_INT*);
domdec_t*	initialDomainDecomposition(graph_t*, PORD_INT*, PORD_INT*, PORD_INT*);
domdec_t*	constructDomainDecomposition(graph_t*, PORD_INT*);
void		computePriorities(domdec_t*, PORD_INT*, PORD_INT*, PORD_INT);
void		eliminateMultisecs(domdec_t*, PORD_INT*, PORD_INT*);
void		findIndMultisecs(domdec_t*, PORD_INT*, PORD_INT*);
domdec_t*	coarserDomainDecomposition(domdec_t*, PORD_INT*);
void		shrinkDomainDecomposition(domdec_t*, PORD_INT);

/* functions in lib/ddbisect.c */
void		checkDDSep(domdec_t*);
PORD_INT       	findPseudoPeripheralDomain(domdec_t*, PORD_INT);
void		constructLevelSep(domdec_t*, PORD_INT);
void		initialDDSep(domdec_t*);
void		updateB2W(bucket_t*, bucket_t*, domdec_t*, PORD_INT, PORD_INT*,
			PORD_INT*, PORD_INT*, PORD_INT*);
void		updateW2B(bucket_t*, bucket_t*, domdec_t*, PORD_INT, PORD_INT*,
			PORD_INT*, PORD_INT*, PORD_INT*);
void		improveDDSep(domdec_t*);

/* functions in lib/gbipart.c */
gbipart_t*      newBipartiteGraph(PORD_INT, PORD_INT, PORD_INT);
void            freeBipartiteGraph(gbipart_t*);
void		printGbipart(gbipart_t*);
gbipart_t*      setupBipartiteGraph(graph_t*, PORD_INT*, PORD_INT, PORD_INT, PORD_INT*);
void            maximumMatching(gbipart_t*, PORD_INT*);
void            maximumFlow(gbipart_t*, PORD_INT*, PORD_INT*);
void            DMviaMatching(gbipart_t*, PORD_INT*, PORD_INT*, PORD_INT*);
void            DMviaFlow(gbipart_t*, PORD_INT*, PORD_INT*, PORD_INT*, PORD_INT*);

/* functions in lib/nestdiss.c */
nestdiss_t*	newNDnode(graph_t*, PORD_INT*, PORD_INT);
void		freeNDnode(nestdiss_t*);
nestdiss_t*	setupNDroot(graph_t*, PORD_INT*);
void		splitNDnode(nestdiss_t*, options_t*, timings_t*);
void		buildNDtree(nestdiss_t*, options_t*, timings_t*);
void		freeNDtree(nestdiss_t*);

/* functions in lib/multisector.c */
multisector_t*	newMultisector(graph_t*);
void		freeMultisector(multisector_t*);
multisector_t*	trivialMultisector(graph_t*);
multisector_t*	constructMultisector(graph_t*, options_t*, timings_t*);
multisector_t*	extractMS2stage(nestdiss_t*);
multisector_t*	extractMSmultistage(nestdiss_t*);

/* functions in lib/gelim.c */
gelim_t*	newElimGraph(PORD_INT, PORD_INT);
void		freeElimGraph(gelim_t*);
void		printElimGraph(gelim_t*);
gelim_t*	setupElimGraph(graph_t*);
PORD_INT       	crunchElimGraph(gelim_t*);
void		buildElement(gelim_t *Gelim, PORD_INT me);
void		updateAdjncy(gelim_t*, PORD_INT*, PORD_INT, PORD_INT*, PORD_INT*);
void		findIndNodes(gelim_t*, PORD_INT*, PORD_INT, PORD_INT*, PORD_INT*, PORD_INT*, PORD_INT*);
void		updateDegree(gelim_t*, PORD_INT*, PORD_INT, PORD_INT*);
void		updateScore(gelim_t*, PORD_INT*, PORD_INT, PORD_INT, PORD_INT*);
elimtree_t*	extractElimTree(gelim_t*);

/* functions in lib/bucket.c */
bucket_t*       newBucket(PORD_INT, PORD_INT, PORD_INT);
void            freeBucket(bucket_t*);
bucket_t*	setupBucket(PORD_INT, PORD_INT, PORD_INT);
PORD_INT        minBucket(bucket_t*);
void            insertBucket(bucket_t*, PORD_INT, PORD_INT);
void            removeBucket(bucket_t*, PORD_INT);

/* functions in lib/minpriority.c */
minprior_t*	newMinPriority(PORD_INT nvtx, PORD_INT nstages);
void		freeMinPriority(minprior_t*);
minprior_t*	setupMinPriority(multisector_t*);
elimtree_t*	orderMinPriority(minprior_t*, options_t*, timings_t*);
void		eliminateStage(minprior_t*, PORD_INT, PORD_INT, timings_t*);
PORD_INT       	eliminateStep(minprior_t*, PORD_INT, PORD_INT);

/* functions in lib/tree.c */
elimtree_t*     newElimTree(PORD_INT, PORD_INT);
void            freeElimTree(elimtree_t*);
void            printElimTree(elimtree_t *);
PORD_INT        firstPostorder(elimtree_t*);
PORD_INT        firstPostorder2(elimtree_t*, PORD_INT);
PORD_INT        nextPostorder(elimtree_t*, PORD_INT);
PORD_INT        firstPreorder(elimtree_t*);
PORD_INT	nextPreorder(elimtree_t*, PORD_INT);
elimtree_t*     setupElimTree(graph_t*, PORD_INT*, PORD_INT*);
void            initFchSilbRoot(elimtree_t*);
void		permFromElimTree(elimtree_t*, PORD_INT*);
elimtree_t*	expandElimTree(elimtree_t*, PORD_INT*, PORD_INT);
elimtree_t*	permuteElimTree(elimtree_t*, PORD_INT*);
elimtree_t*     fundamentalFronts(elimtree_t*);
elimtree_t*     mergeFronts(elimtree_t*, PORD_INT);
elimtree_t*     compressElimTree(elimtree_t*, PORD_INT*, PORD_INT);
PORD_INT        justifyFronts(elimtree_t*);
PORD_INT	nWorkspace(elimtree_t*);
PORD_INT	nFactorIndices(elimtree_t*);
PORD_INT	nFactorEntries(elimtree_t*);
FLOAT           nFactorOps(elimtree_t*);
void		subtreeFactorOps(elimtree_t*, FLOAT*);
FLOAT		nTriangularOps(elimtree_t*);

/* functions in lib/matrix.c */
inputMtx_t*	newInputMtx(PORD_INT, PORD_INT);
void 		freeInputMtx(inputMtx_t*);
void		printInputMtx(inputMtx_t*);
denseMtx_t*     newDenseMtx(workspace_t*, PORD_INT);
void		freeDenseMtx(denseMtx_t*);
void		printDenseMtx(denseMtx_t*);
void		checkDenseMtx(denseMtx_t*);
workspace_t*	initWorkspaceForDenseMtx(PORD_INT, PORD_INT);
FLOAT*		getWorkspaceForDenseMtx(workspace_t*, PORD_INT);
void		freeWorkspaceForDenseMtx(workspace_t*);
inputMtx_t*	setupInputMtxFromGraph(graph_t*);
inputMtx_t*	setupLaplaceMtx(PORD_INT, PORD_INT, PORD_INT);
inputMtx_t* 	permuteInputMtx(inputMtx_t*, PORD_INT*);

/* functions in lib/symbfac.c */
css_t*		newCSS(PORD_INT, PORD_INT, PORD_INT);
void		freeCSS(css_t*);
css_t*		setupCSSFromGraph(graph_t*, PORD_INT*, PORD_INT*);
css_t*		setupCSSFromFrontSubscripts(frontsub_t*);
frontsub_t*	newFrontSubscripts(elimtree_t*);
void		freeFrontSubscripts(frontsub_t*);
void		printFrontSubscripts(frontsub_t*);
frontsub_t*	setupFrontSubscripts(elimtree_t*, inputMtx_t*);
factorMtx_t*	newFactorMtx(PORD_INT);
void		freeFactorMtx(factorMtx_t*);
void		printFactorMtx(factorMtx_t*);
void		initFactorMtx(factorMtx_t *L, inputMtx_t*);
void		initFactorMtxNEW(factorMtx_t *L, inputMtx_t*);

/* functions in lib/numfac.c */
void		numfac(factorMtx_t *L, timings_t *cpus);
denseMtx_t*	setupFrontalMtx(workspace_t*, factorMtx_t*, PORD_INT);
void		initLocalIndices(denseMtx_t*, PORD_INT*, PORD_INT*);
denseMtx_t*	extendedAdd(denseMtx_t*, denseMtx_t*, PORD_INT*, PORD_INT*);
denseMtx_t*	setupUpdateMtxFromFrontalMtx(denseMtx_t*, factorMtx_t*);

/* functions in lib/kernel.c */
denseMtx_t*	factorize1x1Kernel(denseMtx_t*, PORD_INT);
denseMtx_t*	factorize2x2Kernel(denseMtx_t*, PORD_INT);
denseMtx_t*	factorize3x3Kernel(denseMtx_t*, PORD_INT);

/* functions in lib/triangular.c */
void		forwardSubst1x1(factorMtx_t*, FLOAT*);
void		backwardSubst1x1(factorMtx_t*, FLOAT*);
void		forwardSubst1x1NEW(factorMtx_t*, FLOAT*);
void		backwardSubst1x1NEW(factorMtx_t*, FLOAT*);

/* functions in lib/mapping.c */
mapping_t*	newMapping(elimtree_t*, PORD_INT);
void		freeMapping(mapping_t*);
void		printMapping(mapping_t*);
void		listing(mapping_t*, PORD_INT, PORD_INT, PORD_INT, FLOAT*, FLOAT*);
mapping_t*	setupMapping(elimtree_t*, PORD_INT, PORD_INT);
void		split(mapping_t*, PORD_INT, PORD_INT, PORD_INT, PORD_INT*, PORD_INT*, FLOAT*, PORD_INT);

/* functions in lib/interface.c */
elimtree_t*	SPACE_ordering(graph_t*, options_t*, timings_t*);
elimtree_t*	SPACE_transformElimTree(elimtree_t*, PORD_INT);
factorMtx_t*	SPACE_symbFac(elimtree_t*, inputMtx_t*);
void		SPACE_numFac(factorMtx_t*, timings_t*);
void		SPACE_solveTriangular(factorMtx_t *L, FLOAT *rhs, FLOAT *xvec);
void		SPACE_solve(inputMtx_t*, FLOAT*, FLOAT*, options_t*,
			timings_t*);
void		SPACE_solveWithPerm(inputMtx_t*, PORD_INT*, FLOAT*, FLOAT*,
			options_t*, timings_t*);
mapping_t*	SPACE_mapping(graph_t*, PORD_INT*, options_t*, timings_t*);

/* functions in lib/sort.c */
void 		insertUpInts(PORD_INT, PORD_INT*);
void		insertUpIntsWithStaticIntKeys(PORD_INT, PORD_INT*, PORD_INT*);
void		insertDownIntsWithStaticFloatKeys(PORD_INT, PORD_INT*, FLOAT*);
void 		insertUpFloatsWithIntKeys(PORD_INT, FLOAT*, PORD_INT*);
void 		qsortUpInts(PORD_INT, PORD_INT*, PORD_INT*);
void 		qsortUpFloatsWithIntKeys(PORD_INT, FLOAT*, PORD_INT*, PORD_INT*);
void		distributionCounting(PORD_INT, PORD_INT*, PORD_INT*);

/* functions in lib/read.c */
graph_t*	readChacoGraph(char*);
inputMtx_t*	readHarwellBoeingMtx(char*);

/* functions in libPAR/topology.c */
topology_t*	newTopology(PORD_INT);
void		freeTopology(topology_t*);
void		printTopology(topology_t*);
topology_t*	setupTopology(void);
void		recMapCube(topology_t*, PORD_INT, PORD_INT, PORD_INT, PORD_INT, PORD_INT, PORD_INT);
void		sendCube(topology_t*, void*, size_t, PORD_INT);
size_t		recvCube(topology_t*, void*, size_t, PORD_INT);
PORD_INT        myrank(void);

/* functions in libPAR/mask.c */
mask_t*		newMask(PORD_INT);
void		freeMask(mask_t*);
mask_t*		setupMask(PORD_INT, PORD_INT, PORD_INT);

/* functions in libPAR/broadcast.c */
void		broadcastInputMtx(topology_t*, inputMtx_t**);
void		broadcastElimTree(topology_t*, elimtree_t**);
void		broadcastArray(topology_t*, char*, size_t);

/* functions in libPAR/buffer.c */
buffer_t*	newBuffer(size_t);
void		freeBuffer(buffer_t*);
buffer_t*	exchangeBuffer(topology_t*, buffer_t*, PORD_INT);
buffer_t*	setupSymbFacBuffer(frontsub_t*, PORD_INT*);
void		readoutSymbFacBuffer(buffer_t*, frontsub_t*, PORD_INT*);
buffer_t*	setupNumFacBuffer(workspace_t*, mask_t*, PORD_INT);
void		readoutNumFacBuffer(workspace_t*, buffer_t*, denseMtx_t**);
buffer_t*	setupTriangularBuffer(frontsub_t*, PORD_INT*, FLOAT*);
void		readoutTriangularBuffer(buffer_t*, frontsub_t*, PORD_INT*, FLOAT*);

/* functions in libPAR/symbfacPAR.c */
frontsub_t*	newFrontSubscriptsPAR(mask_t*, mapping_t*, elimtree_t*);
frontsub_t*	setupFrontSubscriptsPAR(topology_t*, mask_t*, mapping_t*,
			elimtree_t*, inputMtx_t*);
css_t*		setupCSSFromFrontSubscriptsPAR(mask_t*, mapping_t*,
			frontsub_t*);
void		initFactorMtxPAR(mask_t*, mapping_t*, factorMtx_t*,
			inputMtx_t*);

/* functions in libPAR/numfacPAR.c */
void		numfacPAR(topology_t*, mask_t*, mapping_t*, factorMtx_t*,
			 PORD_INT msglvl, timings_t*);
denseMtx_t*	setupFrontalMtxPAR(mask_t*, PORD_INT, workspace_t*, factorMtx_t*,
			PORD_INT);
void		initLocalIndicesPAR(denseMtx_t*, PORD_INT*, PORD_INT*);
denseMtx_t*	extendedAddPAR(denseMtx_t*, denseMtx_t*, PORD_INT*, PORD_INT*);
denseMtx_t*	setupUpdateMtxFromFrontalMtxPAR(denseMtx_t*, factorMtx_t*);
denseMtx_t*	setupUpdateMtxFromBuffer(workspace_t*, FLOAT*);
void		splitDenseMtxColumnWise(denseMtx_t*, mask_t*, buffer_t*, PORD_INT);
void		splitDenseMtxRowWise(denseMtx_t*, mask_t*, buffer_t*, PORD_INT);

/* functions in libPAR/kernelPAR.c */
denseMtx_t*	factorize1x1KernelPAR(topology_t*, mask_t*, PORD_INT, denseMtx_t*,
			frontsub_t*, timings_t*);
denseMtx_t*	factorize2x2KernelPAR(topology_t*, mask_t*, PORD_INT, denseMtx_t*,
			frontsub_t*, timings_t*);
denseMtx_t*	factorize3x3KernelPAR(topology_t*, mask_t*, PORD_INT, denseMtx_t*,
			frontsub_t*, timings_t*);

/* functions in libPAR/triangularPAR.c */
void		forwardSubst1x1PAR(topology_t*, mask_t*, mapping_t*,
			factorMtx_t*, FLOAT*, FLOAT*);
void		backwardSubst1x1PAR(topology_t*, mask_t*, mapping_t*,
			factorMtx_t*, FLOAT*);
void		forwardSubst1x1KernelPAR(topology_t*, mask_t*, PORD_INT, PORD_INT,
			factorMtx_t*, FLOAT*, FLOAT*);
void		backwardSubst1x1KernelPAR(topology_t*, mask_t*, PORD_INT, PORD_INT,
			factorMtx_t*, FLOAT*);
void		accumulateVector(topology_t*, mask_t*, mapping_t*,
			factorMtx_t*, FLOAT*);

/* functions in libPAR/interfacePAR.c */
topology_t*	SPACE_setupTopology(void);
mask_t*		SPACE_setupMask(topology_t*, PORD_INT);
void		SPACE_cleanup(topology_t*, mask_t*);
factorMtx_t*	SPACE_symbFacPAR(topology_t*, mask_t*, mapping_t*, elimtree_t*,
			inputMtx_t*);
void		SPACE_numFacPAR(topology_t*, mask_t*, mapping_t*, factorMtx_t*,
			PORD_INT msglvl, timings_t*);
void		SPACE_solveTriangularPAR(topology_t*, mask_t*, mapping_t*,
			factorMtx_t*, FLOAT*, FLOAT*);
void		SPACE_solveWithPermPAR(topology_t *top, mask_t *mask,
			 inputMtx_t *A, PORD_INT *perm, FLOAT *rhs, FLOAT *xvec,
			 options_t *options, timings_t *cpus);

