#ifndef __ARPACKDEF_H__
#define __ARPACKDEF_H__

#define INTERFACE64 0

#if INTERFACE64
#define c_int  c_int64_t
#define a_int    int64_t
#define a_uint  uint64_t
#else
#define a_int            int
#define a_uint  unsigned int
#endif

#endif
