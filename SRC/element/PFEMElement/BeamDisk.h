#ifndef BeamDisk_h
#define BeamDisk_h

#include <vector>

#include "Flume.h"

class BeamDisk : public Flume {
   public:
    // data: thk, radius, size, dir
    BeamDisk(int tag, const std::vector<double>& center,
             const std::vector<double>& data);
    ~BeamDisk();
    int mesh();
};

#endif
