#include <CulaSparseSolverS5.h>

#define culaSparse cula;

CulaSparseSolverS5::CulaSparseSolverS5(void)
:SparseGenRowLinSolver(SOLVER_TAGS_CulaSparseS5)
{  
  //config.debug = 1;
  
  single=0;
  status = culaSparseCreate(&handle);
  
  if ( status != culaSparseNoError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  
  status= culaSparsePreinitializeCuda(handle);
  
  if ( status != culaSparseNoError )
    {
      
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  
  
  status = culaSparseCreatePlan(handle, &plan);
  if ( status != culaSparseNoError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  
  status = culaSparseSetCudaPlatform(handle, plan, 0);
  if ( status != culaSparseNoError)
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  status = culaSparseConfigInit(handle,&config);
  
  
  config.relativeTolerance=1.0E-6;
  //config.absoluteTolerance=10.0;
  //config.useInitialResultVector=1;
  config.maxIterations=100000;
  
  
  if ( status != culaSparseNoError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  
  status = culaSparseCsrOptionsInit(handle,&CsrOpt);
  if ( status != culaSparseNoError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  CsrOpt.indexing=0;
  
  culaSparseSetBicgSolver(handle, plan, 0);
  
  culaSparseFainvOptions fainv;
  culaSparseFainvOptionsInit(handle,&fainv);
  fainv.dropTolerance=0.0;
  
  
  culaSparseSetFainvPreconditioner(handle,plan,&fainv);
  
  
}


CulaSparseSolverS5::CulaSparseSolverS5(double relTol,
				       int maxInteration,
				       int preCond,
				       int solver,
				       int single,
				       int host)
 :SparseGenRowLinSolver(SOLVER_TAGS_CulaSparseS5)
{  
  //config.debug = 1;
  
  
  this->relTol=relTol;
  this->maxInteration=maxInteration;
  this->preCond=preCond;
  this->solver=solver;
  this->single=single;
  this->host=host;
  
  
  status = culaSparseCreate(&handle);
  
  if ( status != culaSparseNoError )
    {
      
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  
  if (host==0)
    {
      status= culaSparsePreinitializeCuda(handle);
      
      if ( status != culaSparseNoError )
	{
	  
	  culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
	  opserr<<buf<<"\n";
	  return;
	}
      
    }
  
  
  
  status = culaSparseCreatePlan(handle, &plan);
  if ( status != culaSparseNoError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  
  
  status = culaSparseConfigInit(handle,&config);
  
  
  config.relativeTolerance=relTol;
  config.maxIterations=maxInteration;
  
  
  if ( status != culaSparseNoError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  
  
  
  status = culaSparseCsrOptionsInit(handle,&CsrOpt);
  if ( status != culaSparseNoError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return;
    }
  CsrOpt.indexing=0;
  
  
  
  if (host==1)
    {
      culaSparseHostOptions hstOpt;
      status=culaSparseHostOptionsInit(handle,&hstOpt);
      if ( status != culaSparseNoError )
	{
	  culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
	  opserr<<buf<<"\n";
	  return;
	}
      hstOpt.debug=0;
      
      status=culaSparseSetHostPlatform(handle,plan,&hstOpt);
      if ( status != culaSparseNoError )
	{
	  culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
	  opserr<<buf<<"\n";
	  return;
	}
      
      
    }else if (host==0)
    {
      culaSparseCudaOptions cudaOpt;
      
      status=culaSparseCudaOptionsInit(handle,&cudaOpt);
      if ( status != culaSparseNoError )
	{
	  culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
	  opserr<<buf<<"\n";
	  return;
	}
      cudaOpt.useHybridFormat=0;
      //cudaOpt.reordering=culaSparseSymamdReordering;
      
      status=culaSparseSetCudaPlatform(handle,plan,&cudaOpt);
      if ( status != culaSparseNoError )
	{
	  culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
	  opserr<<buf<<"\n";
	  return;
	}
      
      
    }else
    {
      status=culaSparseArgumentError;
      return;
    }
  
  
  
  switch (preCond)
    {
    case 0: status=culaSparseSetNoPreconditioner(handle,plan,0);
      break;
    case 1: status=culaSparseSetJacobiPreconditioner(handle,plan,0);
      break;
    case 2: status=culaSparseSetBlockJacobiPreconditioner(handle,plan,0);
      break;
    case 3: status=culaSparseSetIlu0Preconditioner(handle,plan,0);
      break;
    case 4: status=culaSparseSetAinvPreconditioner(handle,plan,0);
      break;
    case 5: culaSparseFainvOptions fainv;
      status=culaSparseFainvOptionsInit(handle,&fainv);
      fainv.dropTolerance=0.0;
      status=culaSparseSetFainvPreconditioner(handle,plan,&fainv);
      break;
    default:
      status=culaSparseArgumentError;
      return;
    }
  
  switch (solver)
    {
    case 0: status=culaSparseSetCgSolver(handle,plan,0);
      break;
    case 1: status=culaSparseSetBicgSolver(handle,plan,0);
      break;
    case 2: status=culaSparseSetBicgstabSolver(handle,plan,0);
      break;
    case 3: status=culaSparseSetBicgstablSolver(handle,plan,0);
      break;
    case 4: status=culaSparseSetGmresSolver(handle,plan,0);
      break;
    case 5: status=culaSparseSetMinresSolver(handle,plan,0);
      break;
    default:
      status=culaSparseArgumentError;
      return;
    }
}



CulaSparseSolverS5::~CulaSparseSolverS5(void)
{
  culaSparseDestroyPlan(plan);
  culaSparseDestroy(handle);
  if (single==1)
    {
      
      delete Xsingle;
      delete Bsingle;
      delete Asingle;
    }
}


int
CulaSparseSolverS5::setSize()
{
  int n=theSOE->size;
  int nnz=theSOE->nnz;
  
  if (single==1)
    {
      
      Xsingle=new float[n];
      Bsingle=new float[n];
      Asingle=new float[nnz];
    }
  
  
  
  return 0;
}


//int
//CulaSparseSolverS5::setLinearSOE(SparseGenColLinSOE &theLinearSOE)
int
CulaSparseSolverS5::setLinearSOE(SparseGenRowLinSOE &theLinearSOE)
{
  theSOE = &theLinearSOE;
  return 0;
}


int
CulaSparseSolverS5::sendSelf(int cTAg, Channel &theChannel)
{
  // doing nothing
  return 0;
}


int
CulaSparseSolverS5::recvSelf(int cTag,
			   Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // nothing to do
  return 0;
}

int CulaSparseSolverS5::solve(void)
{
  
  if (theSOE == 0) {
    opserr << "WARNING CulaSparse::solve(void)- ";
    opserr << " No LinearSOE object has been set\n";
    return -1;
  }
  
  if ( status != culaSparseNoError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return -1;
    }
  
  n = theSOE->size;
  nnz = theSOE->nnz;
  
  double *Xptr = theSOE->X;
  double *Bptr = theSOE->B;
  double *Aptr = theSOE->A;
  
  
  int *rowPtr=theSOE->rowStartA;
  int *colInd=theSOE->colA;
  
  
  if (single==0)
    {
      culaSparseSetDcsrData(handle, plan,&CsrOpt,n, nnz, Aptr, rowPtr, colInd, Xptr, Bptr);
    }else if(single ==1)
    {
      
      for (int i=0;i<n;++i)
	{
	  Xsingle[i]=Xptr[i];
	  Bsingle[i]=Bptr[i];
	  Asingle[i]=Aptr[i];
	}
      for (int i=n;i<nnz;++i)
	{
	  Asingle[i]=Aptr[i];
	}
      culaSparseSetScsrData(handle, plan,&CsrOpt,n, nnz, Asingle, rowPtr, colInd, Xsingle, Bsingle);
      
    }else
    
    {
      opserr<<"CulaSparseDataError!\n";
      return -1;
    }
  
  // execute plan
  status = culaSparseExecutePlan(handle, plan, &config, &result);
  
  
  // see if solver failed for a non-data related reason
  if ( status != culaSparseNoError && status != culaSparseDataFormatError )
    {
      culaSparseGetLastStatusString(handle, buf, sizeof(buf) );
      opserr<<buf<<"\n";
      return -1;
    }
  // print result string
  /*
  culaSparseGetResultString(handle,&result, buf, sizeof(buf) );
  opserr<<buf<<"\n";
  */  
  
  if (single==1)
    {
      for (int i=0;i<n;++i)
	{
	  Xptr[i]=Xsingle[i];		
	}
      
    }

  return 0;
}


