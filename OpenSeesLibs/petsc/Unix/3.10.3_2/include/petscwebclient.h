/*
  This is the include file for source code that accesses HTTPS
*/
#if !defined(__PETSCWEBCLIENT_H)
#define __PETSCWEBCLIENT_H

/*  complex number I conflicts with SSL include files */
#define PETSC_SKIP_COMPLEX
#include <petscsys.h>

#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>

#include <openssl/ssl.h>
#include <openssl/err.h>

PETSC_EXTERN PetscErrorCode PetscSSLInitializeContext(SSL_CTX **);
PETSC_EXTERN PetscErrorCode PetscSSLDestroyContext(SSL_CTX *);
PETSC_EXTERN PetscErrorCode PetscHTTPSRequest(const char[],const char[],const char[],const char[],const char[],SSL*,char[],size_t);
PETSC_EXTERN PetscErrorCode PetscHTTPSConnect(const char[],int,SSL_CTX *,int *,SSL **);

PETSC_EXTERN PetscErrorCode PetscHTTPRequest(const char[],const char[],const char[],const char[],const char[],int,char[],size_t);
#endif
