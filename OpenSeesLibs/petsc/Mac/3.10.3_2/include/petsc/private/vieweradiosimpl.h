
#ifndef __VIEWERADIOSIMPL_H
#define __VIEWERADIOSIMPL_H

typedef struct {
  char          *filename;
  PetscFileMode btype;
  PetscInt      timestep;
  int64_t       adios_handle;
  ADIOS_FILE    *adios_fp;
} PetscViewer_ADIOS;

#endif
