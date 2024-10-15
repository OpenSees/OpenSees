//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
#include <tcl.h>
#include <vector>
#include <OPS_Globals.h>
// #include <mpi.h>
#include <Channel.h>
#include <MachineBroker.h>

// #  include <DistributedDisplacementControl.h>
// #  include <ShedHeaviest.h>
// #  include <MPIDiagonalSOE.h>
// #  include <MPIDiagonalSolver.h>
#include <ShadowSubdomain.h>
#include <Metis.h>
#include <FEM_ObjectBroker.h>
#include <DomainPartitioner.h>
#include <domain/domain/partitioned/PartitionedDomain.h>
#include <GraphPartitioner.h>
#include <Subdomain.h>
#include <SubdomainIter.h>
#include <MachineBroker.h>
#include <StaticDomainDecompositionAnalysis.h>
#include <TransientDomainDecompositionAnalysis.h>

// #  define MPIPP_H
// #  include <DistributedSuperLU.h>
// #  include <DistributedProfileSPDLinSOE.h>

 struct PartitionRuntime {
   MachineBroker       *machine            = nullptr;
   FEM_ObjectBroker    *broker             = nullptr;
   DomainPartitioner   *DOMAIN_partitioner = nullptr;
   GraphPartitioner    *GRAPH_partitioner  = nullptr;
// LoadBalancer        *balancer           = nullptr;
   Channel             **channels          = nullptr;  
   int  num_subdomains    = 0;
   bool partitioned       = false;
   bool using_main_domain = false;
   bool setMPIDSOEFlag    = false;
   int  main_partition    = 0;
   PartitionedDomain     theDomain;
 };


static int partitionModel(PartitionRuntime& part, int eleTag);
static Tcl_CmdProc opsPartition;
static Tcl_CmdProc wipePP;
extern Tcl_CmdProc TclCommand_specifyModel;

void 
Init_PartitionRuntime(Tcl_Interp* interp, MachineBroker* theMachineBroker, FEM_ObjectBroker* theBroker)
{
  PartitionRuntime *part = new PartitionRuntime{theMachineBroker, theBroker};

  //
  // set some global parameters
  //
  if (theMachineBroker->getPID() == 0) {
  
    // always use p0 even if ODD number of partitions
    part->num_subdomains    = theMachineBroker->getNP();
    part->using_main_domain = true;
    part->main_partition    = 1;

  } else {
    part->using_main_domain = false;
    part->num_subdomains = 0;
    part->partitioned = false;

  }
  
  
  Tcl_CreateCommand(interp, "partition", &opsPartition, (ClientData)part, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "wipePP",    &wipePP,       (ClientData)part, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "model",     &TclCommand_specifyModel,  (ClientData)&part->theDomain, (Tcl_CmdDeleteProc *)NULL);
}



int
opsPartition(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  PartitionRuntime& part = *static_cast<PartitionRuntime*>(clientData);

  int eleTag;
  if (argc == 2) {
    if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
      ;
    }
  }
  partitionModel(part, eleTag);
  return TCL_OK;
}

static int
partitionModel(PartitionRuntime& part, int eleTag)
{
  if (part.partitioned == true)
    return 0;

  int result = 0;

  if (part.channels != nullptr)
    delete[] part.channels;

  part.channels = new Channel *[part.num_subdomains];

  // create some subdomains
  for (int i = 1; i <= part.num_subdomains; i++) {
    if (i != part.main_partition) {
      ShadowSubdomain *theSubdomain =
          new ShadowSubdomain(i, *part.machine, *part.broker);
      part.theDomain.addSubdomain(theSubdomain);
      part.channels[i - 1] = theSubdomain->getChannelPtr();
    }
  }

  // create a partitioner & partition the domain
  if (part.DOMAIN_partitioner == nullptr) {
    //      part.balancer = new ShedHeaviest();
    // OPS_DOMAIN_partitioner = new DomainPartitioner(*OPS_GRAPH_partitioner, *part.balancer);
    part.GRAPH_partitioner = new Metis;
    part.DOMAIN_partitioner = new DomainPartitioner(*part.GRAPH_partitioner);
    part.theDomain.setPartitioner(part.DOMAIN_partitioner);
  }

  result = part.theDomain.partition(part.num_subdomains, part.using_main_domain,
                                    part.main_partition, eleTag);

  if (result < 0)
    return result;

  part.partitioned = true;

  DomainDecompositionAnalysis *theSubAnalysis;
  SubdomainIter &theSubdomains = part.theDomain.getSubdomains();
  Subdomain *theSub = nullptr;

  void* the_static_analysis = nullptr;
#if 0
  // create the appropriate domain decomposition analysis
  while ((theSub = theSubdomains()) != nullptr) {
    if (the_static_analysis != nullptr) {
      theSubAnalysis = new StaticDomainDecompositionAnalysis(
          *theSub, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
          *theSOE, *theStaticIntegrator, theTest, false);

    } else {
      theSubAnalysis = new TransientDomainDecompositionAnalysis(
          *theSub, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
          *theSOE, *theTransientIntegrator, theTest, false);
    }
    theSub->setDomainDecompAnalysis(*theSubAnalysis);
  }
#endif
  return result;
}


static int 
wipePP(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  PartitionRuntime& part = *static_cast<PartitionRuntime*>(clientData);

  if (part.partitioned == true && part.num_subdomains > 1) {
    SubdomainIter &theSubdomains = part.theDomain.getSubdomains();
    Subdomain *theSub =nullptr;
    
    // create the appropriate domain decomposition analysis
    while ((theSub = theSubdomains()) != nullptr)
      theSub->wipeAnalysis();
  }
  return TCL_OK;  
}

