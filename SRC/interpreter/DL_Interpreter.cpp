
#include "DL_Interpreter.h"

DL_Interpreter  *ops_TheActiveInterpreter = 0;

DL_Interpreter::DL_Interpreter() {

}

DL_Interpreter::~DL_Interpreter() {
  // does nothing
}

int 
DL_Interpreter::addCommand(const char *, Command &) {
  return -1;
}

int 
DL_Interpreter::removeCommand(const char *) {
  return -1;
}

int 
DL_Interpreter::getNumRemainingInputArgs(void) {
  return -1;
}

int 
DL_Interpreter::getInt(int *, int numArgs) {
  return -1;
}

int 
DL_Interpreter::getDouble(double *, int numArgs) {
  return -1;
}

int 
DL_Interpreter::getString(char *cArray, int size) {
  return -1;
}

int 
DL_Interpreter::getStingCopy(char **stringPtr) {
  return -1;
}
