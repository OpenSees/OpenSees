/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * main.c
 *
 * This file contains the driving routine for multilevel method
 *
 * Started 8/28/94
 * George
 *
 * $Id: main.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"


/*************************************************************************
* Global Variables
**************************************************************************/
int dbglvl = DBG_TIME|DBG_OUTPUT;	/* Debuging Level */
int InitPartType; 			/* Type of Initial Partition */
int RefineType;				/* Type of Refinement Algorithm */
int MatchType;				/* Type of matching to use */
int CoarsenTo;				/* The number of nodes to coarsen down */
int OpType;				/* Operation type, order/partition */
int IsWeighted;				/* True if the graph is weighted */

/*************************************************************************
* Global Variables
**************************************************************************/
timer TotalTmr;		/* Times the entire algorithm */
timer CoarsenTmr;	/* Times total coarsening time */
timer GreedyTmr;        /* Times Total Greedy Time */
timer GreedyInitTmr;    /* Times initialization cost of Greedy */
timer GreedyIterTmr;    /* Times Iterative cost of Greedy */
timer GreedyWrapUpTmr;  /* Times the wrap up phase of Greedy */
timer MlevelTmr;	/* Times the entire multilevel algorithm */
timer InitPartTmr;	/* Times the initial partition phase */
timer ProjectTmr;	/* Times the projection of the partition */
timer SplitTmr;		/* Times the boundary creation */
timer BalanceTmr;	/* Times the time required for balancing */
timer IOTmr;		/* Times the file input time */
timer UncrsTmr;		/* Times the file input time */


/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int i;
  CoarseGraphType graph;
  char filename[256];
  int nparts;
  int totalcut;
  int *perm;		/* The array used for ordering */
  int *part, *kpwgts;		/* The array that stores the partitioning */
  int *cuts;		/* Stores the edges crossed at various levels */
  SepNodeType *stree;	/* Stores the separator tree */
  int *smbnz;		/* Stores the nonzeros of the sy,bolically factored matrix */
  int exectype;


  exectype = parsecmd(argv[0]);

  switch (exectype) {
    case ET_METIS:
      if (argc != 9) {
        printf("Usage: %s <GraphFile> <Nparts> <CoarsenTo> <Mtype> <Rtype> <IPtype> <OPtype> <Options> \n",argv[0]);
        PrintOptions();
        exit(0);
      }
    
      strcpy(filename, argv[1]);
      nparts = atoi(argv[2]);
      CoarsenTo = atoi(argv[3]);
      MatchType = atoi(argv[4]);
      RefineType = atoi(argv[5]);
      InitPartType = atoi(argv[6]); 
      OpType = atoi(argv[7]); 
      dbglvl = atoi(argv[8]);

      if (nparts <= 1)
        errexit("Nparts must be greater than 1!");
      break;
    case ET_METIS_P:
      if (argc != 3) {
        printf("Usage: %s <GraphFile> <Nparts>\n",argv[0]);
        exit(0);
      }
    
      strcpy(filename, argv[1]);
      nparts = atoi(argv[2]);
      OpType = OP_RMLB;

      CoarsenTo = PMETIS_CTO;
      MatchType = PMETIS_MTYPE;
      InitPartType = PMETIS_IPTYPE;
      RefineType = PMETIS_RTYPE;

      if (nparts <= 1)
        errexit("Nparts must be greater than 1!");
      break;
    case ET_METIS_P_K:
      if (argc != 3) {
        printf("Usage: %s <GraphFile> <Nparts>\n",argv[0]);
        exit(0);
      }
    
      strcpy(filename, argv[1]);
      nparts = atoi(argv[2]);
      OpType = OP_MLKP;

      CoarsenTo = KMETIS_CTO;
      MatchType = KMETIS_MTYPE;
      InitPartType = KMETIS_IPTYPE;
      RefineType = KMETIS_RTYPE;

      if (nparts <= 1)
        errexit("Nparts must be greater than 1!");
      break;
    case ET_METIS_O:
      if (argc != 2) {
        printf("Usage: %s <GraphFile>\n",argv[0]);
        exit(0);
      }
    
      nparts = 0;
      strcpy(filename, argv[1]);
      OpType = OP_MLND;

      CoarsenTo = OMETIS_CTO;
      MatchType = OMETIS_MTYPE;
      InitPartType = OMETIS_IPTYPE;
      RefineType = OMETIS_RTYPE;

      break;
    default:
      errexit("Internal Error!");
  }


  InitTimers();
  starttimer(&TotalTmr);

  readgraph(&graph, filename);

  PrintInfo(filename, graph, nparts);

  InitRandom();

  PermuteGraphRandom(&graph);


  switch (OpType) {
    case OP_RMLB:
      printf("Recursive Partitioning... ------------------------------------------\n");
      part = ismalloc(graph.nvtxs, 0, "main: part");
      kpwgts = ismalloc(nparts, 0, "main: part");
      cuts = ismalloc(2*nparts, -1, "main: cuts");

      starttimer(&MlevelTmr);
      totalcut = MultiLevelPart(&graph, nparts, CoarsenTo, MatchType, InitPartType, RefineType, dbglvl, IsWeighted, part, cuts, kpwgts);
      stoptimer(&MlevelTmr);

      if (dbglvl&DBG_OUTPUT)
        WritePartition(filename, part, graph.nvtxs, nparts); 

      if (dbglvl&DBG_PARTSIZES) {
        printf("  Partition Sizes:\n    ");
        for (i=1; i<=nparts; i++) {
          printf("%6d ", kpwgts[i-1]);
          if (i%9 == 0)
            printf("\n    ");
        }
        printf("\n    ");
      }

      if (ispow2(nparts)) {
        PrintPartResults("  Edge-Cuts", nparts, cuts);
        printf("  Balance: %4.3f\n", ComputePartBalance(&graph, nparts, kpwgts));
      }
      else {
        printf("  %d-way Edge-Cut: %7d, Balance: %4.3f\n", nparts, totalcut, ComputePartBalance(&graph, nparts, kpwgts));
      }

      GKfree(part, kpwgts, cuts, -1);
      break;
    case OP_MLKP:
      printf("K-way Partitioning... ----------------------------------------------\n");
      part = ismalloc(graph.nvtxs, 0, "main: part");
      kpwgts = ismalloc(nparts, 0, "main: part");

      starttimer(&MlevelTmr);
      totalcut = KWayPart(&graph, nparts, CoarsenTo, MatchType, RefineType, dbglvl, IsWeighted, part, kpwgts);
      stoptimer(&MlevelTmr);

      if (dbglvl&DBG_OUTPUT)
        WritePartition(filename, part, graph.nvtxs, nparts); 

      if (dbglvl&DBG_PARTSIZES) {
        printf("  Partition Sizes:\n    ");
        for (i=1; i<=nparts; i++) {
          printf("%6d ", kpwgts[i-1]);
          if (i%9 == 0)
            printf("\n    ");
        }
        printf("\n    ");
      }
      printf("  %d-way Edge-Cut: %7d, Balance: %4.3f\n", nparts, totalcut, ComputePartBalance(&graph, nparts, kpwgts));

      GKfree(part, kpwgts, -1);
      break;
    case OP_MLND:
      printf("Ordering... --------------------------------------------------------\n");
      perm = imalloc(graph.nvtxs, "main: perm");
      smbnz = imalloc(graph.nvtxs, "main: smbnz");
      stree = (SepNodeType *)GKmalloc(sizeof(SepNodeType)*graph.nvtxs, "main: stree");
      for (i=0; i<graph.nvtxs; i++)
        stree[i].nvtxs = -1;

      starttimer(&MlevelTmr);
      MultiLevelOrder(&graph, CoarsenTo, MatchType, InitPartType, RefineType, dbglvl, perm, stree);
      stoptimer(&MlevelTmr);

      if (dbglvl&DBG_OUTPUT)
        WriteOrder(filename, perm, graph.nvtxs); 

      ComputeFillIn(filename, perm, graph.nvtxs, 3*graph.nedges, smbnz);
      PrintOrderResults(graph.nvtxs, stree, smbnz);

      GKfree(perm, smbnz, stree, -1);
      break;
    default:
      errexit("Unknown");
  }

  stoptimer(&TotalTmr);

  if (dbglvl&DBG_TIME)
    PrintTimers(exectype);

  exit(0);
}  


/*************************************************************************
* This function parses the string used to invoke the rutine to determine
* what operation will perform
**************************************************************************/
int parsecmd(char *cmd)
{
  int i;

  for (i=strlen(cmd)-1; i>=0; i--)
    if (cmd[i] == '/')
      break;

  if (!strcmp(cmd+i+1, "pmetis"))
    return ET_METIS_P;

  if (!strcmp(cmd+i+1, "kmetis"))
    return ET_METIS_P_K;

  if (!strcmp(cmd+i+1, "ometis"))
    return ET_METIS_O;

  return ET_METIS;
}

/*************************************************************************
* This function prints run parameters
**************************************************************************/
void PrintInfo(char *filename, CoarseGraphType graph, int nparts)
{
  printf("********************************************************************\n");
  printf(" METIS 2.0   Copyright 1995, Regents of the University of Minnesota\n");
  printf("\nGraph Information --------------------------------------------------\n");
  printf("   Graph: %s, Size: %5d, %7d, Parts: %3d, Cto: %d\n",
                  filename, graph.nvtxs, graph.nedges, nparts, CoarsenTo);
  printf(" Options: ");
  switch (MatchType) {
    case MATCH_RM:
      printf("RM");
      break;
    case MATCH_HEM: 
      printf("HEM");
      break;
    case MATCH_LEM: 
      printf("LEM");
      break;
    case MATCH_HCM:
      printf("HCM");
      break;
    case MATCH_MHEM:
      printf("HEM*");
      break;
    case MATCH_SRM:
      printf("SRM");
      break;
    case MATCH_SHEM:
      printf("SHEM");
      break;
    case MATCH_SMHEM:
      printf("SHEM*");
      break;
    default:
      printf("??");
  }
  printf(", ");
  switch (RefineType) {
    case REFINE_GR:
      printf("GR");
      break;
    case REFINE_KLR:
      printf("KLR");
      break;
    case REFINE_GKLR:
      printf("GKLR");
      break;
    case REFINE_BGR:
      printf("BGR");
      break;
    case REFINE_BKLR:
      printf("BKLR");
      break;
    case REFINE_BGKLR:
      printf("BGKLR");
      break;
    case REFINE_NONE:
      printf("NR");
      break;
    default:
      printf("??");
  }
  printf(", ");
  switch (InitPartType) {
    case INITPART_GGP:
      printf("GGP");
      break;
    case INITPART_GGGP:
      printf("GGGP");
      break;
    case INITPART_EIG:
      printf("Eig");
      break;
    case INITPART_GGPKL:
      printf("GGPKL");
      break;
    default:
      printf("??");
  }
  printf(", ");
  switch (OpType) {
    case OP_RMLB:
      printf("Rec-Partition");
      break;
    case OP_MLKP:
      printf("K-Partition");
      break;
    case OP_MLND:
      printf("Order");
      break;
    default:
      errexit("??");
  }
  printf("\n\n");
}


/*************************************************************************
* This function prints run parameters
**************************************************************************/
void PrintOptions(void)
{
  printf(" MType: RM: %d, HEM: %d, LEM: %d, HCM: %d, HEM*: %d, SRM: %d, SHEM: %d, SHEM*: %d\n",
                MATCH_RM, MATCH_HEM, MATCH_LEM, MATCH_HCM, MATCH_MHEM, MATCH_SRM, MATCH_SHEM, MATCH_SMHEM);
  printf(" RType: GR: %d, KLR: %d, GKLR: %d, BGR: %d, BKLR: %d, BGKLR: %d, NR: %d\n",
                REFINE_GR, REFINE_KLR, REFINE_GKLR, REFINE_BGR, REFINE_BKLR, REFINE_BGKLR, REFINE_NONE);
  printf("IPType: GGP: %d, GGGP: %d, EIG: %d, GGPKL: %d\n",
                INITPART_GGP, INITPART_GGGP, INITPART_EIG, INITPART_GGPKL);
  printf("OPtype: Rec-Partition: %d, Kway-Partition: %d, Order: %d\n",
                OP_RMLB, OP_MLKP, OP_MLND);
}


/*************************************************************************
* This function initializes the timers
**************************************************************************/
void InitTimers(void)
{
  cleartimer(&TotalTmr);
  cleartimer(&CoarsenTmr);
  cleartimer(&GreedyTmr);
  cleartimer(&GreedyInitTmr);
  cleartimer(&GreedyIterTmr);
  cleartimer(&GreedyWrapUpTmr);
  cleartimer(&MlevelTmr);
  cleartimer(&InitPartTmr);
  cleartimer(&ProjectTmr);
  cleartimer(&SplitTmr);
  cleartimer(&IOTmr);
  cleartimer(&BalanceTmr);
  cleartimer(&UncrsTmr);
}

/*************************************************************************
* This function prints the various timers
**************************************************************************/
void PrintTimers(int exectype)
{
  printf("\nTiming Information -------------------------------------------------");
  printf("\n Multilevel: \t\t %7.3f", gettimer(&MlevelTmr));
  printf("\n     Coarsening: \t\t %7.3f", gettimer(&CoarsenTmr));
  printf("\n     Initial Partition: \t %7.3f", gettimer(&InitPartTmr));
  printf("\n     Uncoarsening: \t\t %7.3f", gettimer(&UncrsTmr));
  if (exectype == ET_METIS) {
    printf("\n        Refinement: \t\t\t %7.3f", gettimer(&GreedyTmr));
    printf("\n            Initialize: \t\t\t %7.3f", gettimer(&GreedyInitTmr));
    printf("\n            Iterate: \t\t\t\t %7.3f", gettimer(&GreedyIterTmr));
    printf("\n            Wrap Up: \t\t\t\t %7.3f", gettimer(&GreedyWrapUpTmr));
    printf("\n        Projection: \t\t\t %7.3f", gettimer(&ProjectTmr));
    if (OpType != OP_MLND)
      printf("\n        Balancing:  \t\t\t %7.3f", gettimer(&BalanceTmr));
  }
  if (OpType != OP_MLKP)
    printf("\n     Splitting: \t\t %7.3f", gettimer(&SplitTmr));
  printf("\n I/O: \t\t\t %7.3f", gettimer(&IOTmr));
  printf("\n Total: \t\t %7.3f", gettimer(&TotalTmr));
  printf("\n********************************************************************\n");
}


/*************************************************************************
* This function checks if n is a power of two
**************************************************************************/
int ispow2(int n)
{
  int i=1;

  while (i < n)
    i = 2*i;

  if (i == n)
    return 1;
  else
    return 0;
}

