CONTENTS:
As to the organisation: For every class i provide there will eventually
be 3 files, a class.h, a class.C and a latex class.tex file. all files
are in the following subdirectories:

doc: - contains some files concerning documentation, model.tex a latex
file and model.ps the postscript (every class does not yet have an upto
date class file .. prof. fenves is deciding on the format the
documentation will have .. the .ps file is around 200 pages if you want
to look at it and give your own comments) there is a latex2html
version of this file at:
http://eig.stanford.edu/~g3

element: - contains the element classes, all in subdirectories, i.e. the
Truss.h and Truss.C files are in a subdirectory truss.

domain: - contains the classes for the domain. all in subdirectores:
        domain - the Domain class and abstract Iters for the objects of
                the Domain, two subdirectories single- contains iters  
                for Domain, partitioned contains the PartitionedDomain
                class and the iters for this class.
        node - the Node and NodalLoad classes.
        constraints - the SP_Constraint and MP_Constraint classes.
        subdomain - the Subdomain, ActorSubdomain and ShadowSubdomain
                classes.
        pattern - the abstact classes LoadPattern and TimeSeries and some
		concrete classes.
        partitioner - the DomainPartitioner class.
        loadBalancer - the LoadBalancer class and some concrete
                subclasses.

analysis: - contains the classes for the analysis. again the classes in
        subdirectories:
        analysis - the Analysis, StaticAnalysis, and some others
        algorithm - SolutionAlgorithm class and 2 subdirectories:
                equiSolnAlgo - Linear, NewtonRaphson and ModifiedNewton
                domainDecompAlgo - DomainDecompAlgo
        dof_grp - DOF_Group
        fe_ele - FE_Element class and subdirectory penalty containing
                PenaltySP_FE and PenaltyMP_FE
        penalty - ConstraintHandler, PenaltyConstraintHandler and
                PlainHandler
        integrator - Integrator, IncrementalIntegrator, LoadControl,
                ArcLength, StaticIntegrator, TransientIntegrator and
                Newmark
        model - AnalysisModel and its iters.
        numberer - DOF_Numberer
                
graph: - contains the graph classes, again in subdirectories
        graph - Graph, Vertex, and others
        partitioner - GraphPartitioner and Metis, also in metis-2.0
                contains the code downloaded to build metis.
        numberer - GraphNumberer, RCM and some others

system_of_eqn: - contains the SystemofEqn class and linearSOE subclass:
        linearSOE: contains LinearSOE, LinearSolver and DomainSolver
        classes and a bunch of subdirectories for each soe and solver:
        
        fullGEN - for a full general solver
        bandGEN - for a banded general solver
        profileSPD - for my profile solver and a solver using skypack,
                developed by O.Marques, now working with Jim Demmel
                (the code for skypack in a subdirectory)
        petsc - for the petsc solver
        sparseGEN - for superLU and thraeded superLU.
        symSparse - for Kincho Law's symmetric sparse solver

        slowMatrix - to be removed.
        eleByEle - to be filled in at some stage.

tagged: - as not all c++ compilers do templates yet i use TagggedObject
        and the files in it's subdirectory storage (TaggedObjectStorage,
        ArrayOfTaggedObjects and their iters) for my containers to store
        the objects in Domain and AnalysisModel.

modelbuilder: - contains ModelBuilder class and some others. the
        subdirectory triangle contains TrianglePlaneStress class and the
        subdirectory Triangle contains J. Shewchucks triangle code.
	
	subdircetory tcl containts the TclModelBuilder.
        
tcl: - contains the stuff for the interpreter (tkAppInit.C, tkConfig.h
        tkMain.C, commands.C)

utility: - contains the Timer class.
 
material: - Material, UniaxialMaterial and some concrete classes need by
        my Truss ele and Filip's element.
        
actor: the classes for my parallel stuff, again in subdirectories:
        shadow - the Shadow class.  
        actor - the Actor class
        channel - the Channel, TCP_Socket and MPI_Channel classes
        machineBroker - the MachineBroker and some classes for the
                machiines over here.
        objectBroker - the FEM_ObjectBroker class.
        message - the Message class
        address - the Address class.

matrix: contains the Matrix, Vector and ID classes.


in this directory there is the Makefile, Makefile.incl for the current
include paths needed for the files in the framework, classTags.h
(contains class tags), and bool.h (some compilers define bool, some don't)
        
there is a bit of fluff in what i have given you and it needs to be
cleaned out. NEXT TIME.



NOTES FOR MYSELF:
notes for alpha cluster:
  1) change 'ssh' in ~/remote/remote.c to 'rsh'
  2) remove the threaded profile solver

notes for holden:
  1) create a HoldenMachineBroker which starts processes on the
     alpha cluster.

notes for millenium machines:
  1) not yet working for CC compiler, need to bcopy() struct addr to
     addr_in
  2) the parallel domain decomposition is using ssh to start the
     remote processes, has to be set up so don't need to supply password -
     this can be done from mill.cs - but not from any other millenium machine!
     they must have some problem with ssh. talk to eric again.
     rsh would be nice! see if eric will allow rsh until fix ssh?



