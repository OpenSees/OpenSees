//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// Written: fmk
//
#ifndef OPS_PACKAGES_H
#define OPS_PACKAGES_H
int getLibraryFunction(const char *libName, const char *functName,
                       void **libHandle, void **funcHandle);
#endif // OPS_PACKAGES_H
