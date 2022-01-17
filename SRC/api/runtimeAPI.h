#ifndef OPS_RUNTIME_API_H
#define OPS_RUNTIME_API_H

class G3_Runtime;
#ifdef OPS_USE_RUNTIME
#  define OPS_ADD_RUNTIME_VPV(func) (func)([[maybe_unused]] G3_Runtime *rt)
#  include <g3_api.h>
#else
#  define OPS_ADD_RUNTIME_VPV(func) (func)()
#endif

#endif // OPS_RUNTIME_API_H
