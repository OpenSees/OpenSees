#ifndef BeamBrick_h
#define BeamBrick_h

#include <map>

#include "Flume.h"

class BeamBrick : public Flume {
   public:
    BeamBrick(int tag, const std::vector<double>& crds,
              const std::vector<double>& dimensions);
    ~BeamBrick();

    virtual int create_line(Node* nd1, Node* nd2, int& nodeTag,
                            int dir);
    virtual int create_face(Node* nd1, Node* nd2, int& nodeTag,
                            int dir1, int dir2);

    int mesh();

   private:
    std::vector<int> elenodes;
    std::map<std::pair<int, int>, std::vector<int>>
        linenodes;  // key is (nodeTag, dir), value is line node tags
};

#endif
