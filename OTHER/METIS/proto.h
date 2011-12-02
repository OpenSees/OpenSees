/*
 * proto.h
 *
 * This file conatins function prototypes
 *
 * Started 8/27/94
 * George
 *
 * $Id: proto.h,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */


/* entrypoint.o */
int PMETIS(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
int KMETIS(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
int OMETIS(int *, int *, int *, int *, int *, int *, int *);
void ConvertGraph(CoarseGraphType *, int, int *, int *, int *, int *, int, int);
void CleanUpRootGraph(void);


/* balpart.c */
void FastInitBalance(CoarseGraphType *, int);
void FastBalance(CoarseGraphType *, int, int);
void FastBalance2(CoarseGraphType *graph, int);


/* bucketlist.c */
void initbucket(BucketListType *, int, int, int, int);
void resetbucket(BucketListType *);
void freebucket(BucketListType *);
int Add2Part(BucketListType *, int, int);
int UpdatePart(BucketListType *, int, int, int);
int GetMaxGainVtx(BucketListType *);
int SeeMaxGainVtx(BucketListType *);
void PrintPlusPart(char *, BucketListType *);
void PrintPartGains(BucketListType *, BucketListType *);


/* coarsen.c */
CoarseGraphType *Coarsen(CoarseGraphType *, int);
CoarseGraphType *KwayCoarsen(CoarseGraphType *, int);
int SelectMatching(CoarseGraphType *, int);
void KwaySelectMatching(CoarseGraphType *, int);


/* fm.c */
void FMR_Refine(CoarseGraphType *, int, int);
void BFMR_Refine(CoarseGraphType *, int, int, int);
void BFMR_Refine_Weighted(CoarseGraphType *, int, int, int);
void BFMR_Refine_EqWgt(CoarseGraphType *, int);
void Greedy_Refine(CoarseGraphType *, int);
void printwhere(CoarseGraphType *);
int CheckBndSize(CoarseGraphType *);


/* htable.c */
void CreateHTable(HTableType *, int);
int SelHTSize(int);
void AddHTable(HTableType *, int);
void DelHTable(HTableType *, int);
void IncreaseHTable(HTableType *);


/* initpart.c */
void InitPartition(CoarseGraphType *graph, int);
void GGPPartition(CoarseGraphType *, int, int);
void GGGPPartition(CoarseGraphType *, int);
void EigPartition(CoarseGraphType *, int);
int inccompeinz(const void *, const void *);
int ComputeCut(CoarseGraphType *);


/* io.c */
void readgraph(CoarseGraphType *, char *);
void ReadGraphSKIT(int *, int *, char *);
void WriteGraph(CoarseGraphType *, char *);
void printgraph(CoarseGraphType *, FILE *);
void WriteOrder(char *, int *, int);
void WritePartition(char *, int *, int, int);
void WriteCoarseGraph2Grid(char *, CoarseGraphType *);


/* kwayfm.c */
void KWay_RefineFM(CoarseGraphType *, int, int);
void KWayFMUpdateDegreesI(CoarseGraphType *, int, int, BucketListType *, int *, int);
void KWay_BalanceFM(CoarseGraphType *, int, int);
void KWayFMUpdateDegreesBal(CoarseGraphType *, int, int, BucketListType *, int *);
int GetMaxEwgtI(int, EdgeType *);


/* kwaygreedy.c */
void KWay_RefineGreedy(CoarseGraphType *, int, int);
void KWayUpdateDegrees(CoarseGraphType *, int, int);
void KWayUpdateVtxDegrees(CoarseGraphType *, int, int, int, int);


/* kwaypart.c */
int KWayPart(CoarseGraphType *, int, int, int, int, int, int, int *, int *);
void KWayRefine(CoarseGraphType *, CoarseGraphType *, int, int *);
void KWayComputePartitionParams(CoarseGraphType *, int, int);
void KWayProjectPartition(CoarseGraphType *);
int KWayCheckDegrees(CoarseGraphType *);


/* lanczos.c */
int lanczos(CoarseGraphType *, double *);
void givens(double *, double *, int, double, double *);
void pivot(double *, double *, int, double, double *);
void matvec(CoarseGraphType *, double *, double *);


/* list.c */
ListNodeType *addthisnode(ListNodeType *, ListNodeType *);
ListNodeType *delthisnode(ListNodeType *, int, ListNodeType **);


/* main.c */
int parsecmd(char *);
void PrintInfo(char *, CoarseGraphType, int);
void PrintOptions(void);
void InitTimers(void);
void PrintTimers(int);
int ispow2(int);


/* match.h */
int RM_Match(CoarseGraphType *);
int RM_Match_W(CoarseGraphType *);
int HEM_Match(CoarseGraphType *);
int HEM_Match_W(CoarseGraphType *);
int LEM_Match(CoarseGraphType *);
int LEM_Match_W(CoarseGraphType *);
int HCM_Match(CoarseGraphType *);
int HCM_Match_W(CoarseGraphType *);
int MHEM_Match(CoarseGraphType *);
int MHEM_Match_W(CoarseGraphType *);
int SRM_Match(CoarseGraphType *);
int SHEM_Match(CoarseGraphType *);
int SMHEM_Match(CoarseGraphType *);
int CreateCoarseGraph(CoarseGraphType *, int);
void mergevertices(VertexType **, int *, int, VertexType *, int, VertexType *);


/* memory.c */
void AllocatePools(CtrlType *);
void FreePools(CtrlType *);
CoarseGraphType *CreateGraph(void);
void FreeRootGraph(CoarseGraphType *);
void FreeGraph(CoarseGraphType *);
void InitGraph(CoarseGraphType *);
void ResetPools(void);
EdgeType *GetEdgePool(void);
int SetEdgePool(int);
void FreeEdgePool(int);
int EdgePoolSizeLeft(void);
EdgeType *GetnExtDegrees(int);
void ResetExtDegrees(void);
int *icoremalloc(int, char *, int);
void icorefree(int);


/* mincover.c */
void MinCover(int *, int *, int, int, int *, int *);
int MinCover_Augment(int *, int *, int, int *, int *, int *, int);
void MinCover_Decompose(int *, int *, int, int, int *, int *, int *);
void MinCover_ColDFS(int *, int *, int, int *, int *, int);
void MinCover_RowDFS(int *, int *, int, int *, int *, int);


/* mlevelorder.c */
void MultiLevelOrder(CoarseGraphType *, int, int, int, int, int, int *, SepNodeType *); 
void MLND(CoarseGraphType *, int *, int, SepNodeType *, int); 
void SplitGraphOrder(CoarseGraphType *, CoarseGraphType *, CoarseGraphType *, int *, int);
void SimpleOrder(CoarseGraphType *, int *, int); 
void MDOrder(CoarseGraphType *, int *, int, SepNodeType *, int); 


/* mlevelpart.c */
int MultiLevelPart(CoarseGraphType *, int, int, int, int, int, int, int, int *, int *, int *);
int RMLB(CoarseGraphType *, int, int *, int, int, int *, int *, int *); 
void SplitGraphPart(CoarseGraphType *, CoarseGraphType *, CoarseGraphType *);
void SplitGraphPart1_2(CoarseGraphType *, CoarseGraphType *);


/* mmd.c */
void genmmd(int, int *, int *, int *, int *, int,int *, int *, int *, int *, int , int *);


/* myrand48.c */
double drand48(void);
double erand48(unsigned short *);
void srand48(long);


/* refine.c */
void Refine(CoarseGraphType *, CoarseGraphType *, int);
void ComputePartitionParams(CoarseGraphType *);
void ProjectPartition(CoarseGraphType *, int);


/* separator.c */
void FindMinCovNodeSeparator(CoarseGraphType *, int *, int);


/* smbfactor.c */
void ComputeFillIn(char *, int *, int, int, int *);
void smbfactor(int *, int *, int, int *, int, int *, int *);
void ComputeElTree(int, int *, int *, int *, int *);


/* smbfct.c */
int smbfct(int, int *, int *, int *, int *, int *, int *, int *, int *, int *);


/* stat.c */
void PrintGraphMMM(CoarseGraphType *);
void PrintPartResults(char *, int, int *);
void PrintOrderResults(int, SepNodeType *, int *);
void CalcNodeOpc(SepNodeType *, int *, int);
double CalcParOpc(SepNodeType *, int, int);
float ComputePartBalance(CoarseGraphType *, int, int *);


/* util.c */
void InitRandom(void);
void PermuteGraphRandom(CoarseGraphType *);
void RandomPermute(int *, int, int);
float RelDiff(int, int);
int CheckDegrees(CoarseGraphType *);
void SortKeyValueNodesDec(KeyValueType *, int);
int DecKeyValueCmp(const void *, const void *);




