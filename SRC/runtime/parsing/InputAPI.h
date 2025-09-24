#ifndef G3PARSE_H
#define G3PARSE_H
#ifndef G3_RUNTIME_H
#include <runtime/runtime/G3_Runtime.h>
#define G3_RUNTIME_H
#endif


#define G3_Char TCL_Char

enum SuccessFlag {
  G3_OK    = TCL_OK, 
  G3_ERROR = TCL_ERROR
};

typedef enum SuccessFlag SuccessFlag;

#endif // G3PARSE_H
