#ifndef OPS_PACKAGES_H
#define OPS_PACKAGES_H
/*                                                                        
** $Revision: 1.1 $
** $Date: 2008-12-01 23:33:16 $
** $Source: /usr/local/cvs/OpenSees/SRC/api/packages.h,v $
                                                                        
** Written: fmk 
*/


int getLibraryFunction(const char *libName, const char *functName, void **libHandle, void **funcHandle);
#endif // OPS_PACKAGES_H
