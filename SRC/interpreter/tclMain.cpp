
#include "TclInterpreter.h"

int main(int argc, char **argv) {
  TclInterpreter *theInterpreter = new TclInterpreter(argc, argv);
  int res = theInterpreter->run();
  
  delete theInterpreter;
  return res;
}
