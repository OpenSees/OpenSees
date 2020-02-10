#!/usr/bin/env python

from __future__ import print_function
Sorted = ["PetscBarrier",

         "ThreadCommRunKer",
         "ThreadCommBarrie",

         "TSSetUp",
         "TSStep",
         "TSFunctionEval",
         "TSJacobianEval",

         "SNESSolve",
         "SNESFunctionEval",
         "SNESJacobianEval",
         "MatFDColorCreate",
         "MatFDColorSetUp",
         "MatFDColorApply",
         "MatFDColorFunc",
         "MatAssemblyBegin",
         "MatAssemblyEnd",

         "SNESNGSEval",
         "SNESNGSFuncEval",
         "SNESLineSearch",
         "SNESNPCSolve",

         "KSPSetUp",
         "PCSetUp",
         "MatLUFactor",
         "MatLUFactorSym",
         "MatLUFactorNum",
         "MatCholeskyFctr",
         "MatCholFctrSym",
         "MatCholFctrNum",
         "MatILUFactor",
         "MatILUFactorSym",
         "MatICCFactorSym",

         "KSPSolve",
         "MatMult",
         "MatMult MF",
         "KSPGMRESOrthog",
         "VecMDot",
         "VecMAXPY",
         "VecNormalize",

         "PCSetUpOnBlocks",
         "PCApplyOnBlocks",
         "PCApply",
         "MatSolve",
         "MatSOR",

         "PCApplyCoarse",
         "PCApplyMultiple",
         "PCApplySymmLeft",
         "PCApplySymmRight",
         "",
         "PCModifySubMatri",

         "PCGAMGGraph_AGG",
         "PCGAMGGraph_GEO",
         "PCGAMGCoarse_AGG",
         "PCGAMGCoarse_GEO",
         "PCGAMGProl_AGG",
         "PCGAMGProl_GEO",
         "PCGAMGPOpt_AGG",


         "DMConvert",
         "DMGlobalToLocal",
         "DMLocalToGlobal",
         "DMDALocalADFunc",
         "DMPlexInterpolate",
         "DMPlexPartition",
         "DMPlexDistribute",
         "DMPlexDistCones",
         "DMPlexDistLabels",
         "DMPlexDistSF",
         "DMPlexDistOvrlp",
         "DMPlexDistField",
         "DMPlexDistData",
         "DMPlexStratify",
         "DMPlexPreallocate",
         "DMPlexResidualFEM",
         "DMPlexJacobianFEM",


         "MatMults",
         "MatMultConstr",
         "MatMultAdd",
         "MatMultTranspose",
         "MatMultTrConstr",
         "MatMultTrAdd",
         "MatSolves",
         "MatSolveAdd",
         "MatSolveTranspos",
         "MatSolveTrAdd",
         "MatForwardSolve",
         "MatBackwardSolve",
         "MatCopy",
         "MatConvert",
         "MatScale",
         "MatResidual",
         "MatSetValues",
         "MatGetValues",
         "MatGetRow",
         "MatGetRowIJ",
         "MatGetSubMatrice",
         "MatGetOrdering",
         "MatIncreaseOvrlp",
         "MatPartitioning",
         "MatCoarsen",
         "MatZeroEntries",
         "MatLoad",
         "MatView",
         "MatAXPY",
 
         "MatTranspose",
         "MatMatMult",
         "MatMatSolve",
         "MatMatMultSym",
         "MatMatMultNum",
         "MatMatMatMult",
         "MatMatMatMultSym",
         "MatMatMatMultNum",
         "MatPtAP",
         "MatPtAPSymbolic",
         "MatPtAPNumeric",
         "MatRARt",
         "MatRARtSym",
         "MatRARtNum",
         "MatMatTransMult",
         "MatMatTrnMultSym",
         "MatMatTrnMultNum",
         "MatTrnMatMult",
         "MatTrnMatMultSym",
         "MatTrnMatMultNum",
         "MatTrnColorCreate",
         "MatGetRedundant",
         "MatGetSeqNZStrct",
         "MatGetMultiProcBlock",
         "MatMPISumSeqNumeric",
         "MatMPISumSeqSymbolic",
         "MatMPISumSeq",
         "MatMPIConcateSeq",
         "MatGetLocalMat",
         "MatGetLocalMatCondensed",
         "MatGetBrowsOfAcols",
         "MatGetBrAoCol",
         "MatApplyPAPt_Symbolic",
         "MatApplyPAPt_Numeric",
         "MatApplyPAPt",
         "MatGetSymTrans",
         "MatGetSymTransR",
         "MatTranspose_SeqAIJ_FAST",
         "MatCUSPCopyTo",
         "MatCUSPARSECopyTo",
         "MatViennaCLCopyTo",
         "MatSetValBatch",
         "MatSetValBatch1",
         "MatSetValBatch2",
         "MatSetValBatch3",
         "MatSetValBatch4",
         "MatColoringApply",
         "MatColoringComm",
         "MatColoringLocal",
         "MatColoringIS",
         "MatColoringSetUp",

         "VecView",
         "VecMax",
         "VecMin",
         "VecDot",
         "VecDotNorm2",
         "VecTDot",
         "VecMTDot",
         "VecNorm",
         "VecScale",
         "VecCopy",
         "VecSet",
         "VecAXPY",
         "VecAYPX",
         "VecAXPBYCZ",
         "VecWAXPY",
         "VecSwap",
         "VecOps",
         "VecAssemblyBegin",
         "VecAssemblyEnd",
         "VecPointwiseMult",
         "VecSetValues",
         "VecLoad",
         "VecScatterBegin",
         "VecScatterEnd",
         "VecSetRandom",
         "VecReduceArith",
         "VecReduceComm",
         "VecReduceBegin",
         "VecReduceEnd",


         "PetscSFSetGraph",
         "PetscSFBcastBegin",
         "PetscSFBcastEnd",
         "PetscSFReduceBegin",
         "PetscSFReduceEnd",
         "PetscSFFetchOpBegin",
         "PetscSFFetchOpEnd"]


def ComputeTotals(localTimes,localFlops,localMessages,localMessageLens,localReductions):
  '''Computes time spent by all processes'''
  time  = 0
  flops = 0
  for t in localTimes:
    time += localTimes[t]
  flops = 0
  for t in localFlops:
    flops += localFlops[t]
  messages = 0
  for t in localMessages:
    messages += localMessages[t]
  messagelens = 0
  for t in localMessageLens:
    messagelens += localMessageLens[t]
  reductions = 0
  for t in localReductions:
    reductions += localReductions[t]
  return time,flops,messages,messagelens,reductions

def ComputeSums(Stages):
  ''' Computes the sum over all processes for each event; removes events with zero count and summaries'''
  sumStages = {}
  for stages in Stages:
    sumStages[stages] = {}
    for events in Stages[stages]:
      sumStages[stages][events] = {}
      for t in Stages[stages][events][0]:
        sumStages[stages][events][t] = 0
        for s in Stages[stages][events]:
          sumStages[stages][events][t] += Stages[stages][events][s][t]

  for stages in Stages:
    for events in Stages[stages]:
      if events == "summary" or not sumStages[stages][events]["count"] or events.startswith('Thread') or events == 'TSFunctionEval' or events == 'TSJacobianEval':
        sumStages[stages].pop(events,None)
  return sumStages

def ObjectsCompare(a,b):
  return Sorted.index(a) - Sorted.index(b)

def PrintPercentTable(localTimes,localFlops,localMessages,localMessageLens,localReductions,Stages,Latex = False):
  ''' Prints a simple table that displays the percent of time, flops, etc for each event in each stage'''
  if Latex:
    print("\documentclass{article}")
    print("\\begin{document}")
    print("\\begin{table}[!htbp]")
    print("\centering")
    if len(localTimes) > 1:
      print("\\begin{tabular}{lcccccc}")
      print(" &  & \multicolumn{4}{c}{--------------- Percent of -------------} &  \\\\")
      print("Event & Count & Time & Flops & Messages & Reductions & Flop rate \\\\")
      print("\hline")
  else:
    if len(localTimes) > 1:
       print("                                       ---------  Percent of  ------")
       print("Event                       Count     Time  Flops Messages Reductions     Flop rate")
       print("============================================================================")
    else:
       print("                                      Percent of")
       print("Event                       Count     Time  Flops      Flop rate")
       print("=========================================================")

  time,flops,numMessages,numMessageLen,numReductions  = ComputeTotals(localTimes,localFlops,localMessages,localMessageLens,localReductions)
  if not numMessages: numMessages = 1
  if not numReductions: numReductions = 1
  sumStages = ComputeSums(Stages)
  for stage in sumStages:
    if len(sumStages) > 1: 
      if Latex:
        pass
      else:
        print("Stage: "+stage)
    L = sumStages[stage].keys()
    L.sort(cmp=ObjectsCompare)
    seperatoradded = False
    for i in xrange(len(L)):
      event  = L[i]
      space  = 0
      if event.startswith("SNESSolve"):        space = 1
      elif event.startswith("SNES"):           space = 2
      elif event.startswith("KSPGMRESOrthog"): space = 3
      elif event.startswith("KSP"):            space = 2
      elif event.startswith("PC"):             space = 3
      elif event.startswith("MatMult"):        space = 3
      elif event.startswith("Mat"):            space = 4
      elif event.startswith("Vec"):            space = 5
      if Sorted.index(event) > Sorted.index("PCGAMGPOpt_AGG"):
        space = 2
        if not seperatoradded: 
          if Latex:
            print("--Overlapping events---\\\\")
          else:
            print("--Overlapping events---")
          seperatoradded = True

      if len(localTimes) > 1:
        values = [100*sumStages[stage][event]["time"]/time,100*sumStages[stage][event]["flops"]/flops,100*sumStages[stage][event]["numMessages"]/numMessages,100*sumStages[stage][event]["numReductions"]/numReductions]
        if max(values) > .5:
          if Latex:
            print('\\hspace{%1dem}' % space,event,"&",'%6.0f' % sumStages[stage][event]["count"],"&",'%5.0f' % values[0],"&",'%5.0f' % values[1],"&",'%5.0f' % values[2],"&",'%5.0f' % values[3],"&",'%8.0f' % ((sumStages[stage][event]["flops"]/sumStages[stage][event]["time"])/1000000.0),"\\\\")
          else:
            print("            "[0:space],event.ljust(26-space),'%6.0f' % sumStages[stage][event]["count"],"   ",'%5.0f' % values[0],'%5.0f' % values[1],'%5.0f' % values[2],'%5.0f' % values[3],"        ",'%8.0f' % ((sumStages[stage][event]["flops"]/sumStages[stage][event]["time"])/1000000.0))
      else:
        values = [100*sumStages[stage][event]["time"]/time,100*sumStages[stage][event]["flops"]/flops]
        if max(values) > .5:
          if Latex:
            print('\\hspace{%1dem}' % space,event,"&",'%6.0f' % sumStages[stage][event]["count"],"&",'%5.0f' % values[0],"&",'%5.0f' % values[1],"&",'%8.0f' % ((sumStages[stage][event]["flops"]/sumStages[stage][event]["time"])/1000000.0),"\\\\")
          else:
            print("  ",event.ljust(24),'%6.0f' % sumStages[stage][event]["count"],"   ",'%5.0f' % values[0],'%5.0f' % values[1],"   ",'%8.0f' % ((sumStages[stage][event]["flops"]/sumStages[stage][event]["time"])/1000000.0))

  if Latex:
    print("\end{tabular}")
    print("\end{table}")
    print("\end{document}")

if __name__ == '__main__':
  import sys
  import os
  sys.path.append(os.getcwd())
  datafile = sys.argv[1]
  if datafile.endswith('.py'): datafile = datafile[0:-3]
  exec('import '+datafile+' as data')
  latex = False
  if len(sys.argv) > 2: latex = True

  PrintPercentTable(data.LocalTimes,data.LocalFlops,data.LocalMessages,data.LocalMessageLens,data.LocalReductions,data.Stages,Latex = latex)
