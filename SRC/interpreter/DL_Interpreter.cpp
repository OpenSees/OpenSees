
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

const char*
DL_Interpreter::getString() {
  return 0;
}

int 
DL_Interpreter::getStingCopy(char **stringPtr) {
  return -1;
}

void
DL_Interpreter::resetInput(int cArg)
{
}

int
DL_Interpreter::setInt(int *, int numArgs) {
    return -1;
}

int
DL_Interpreter::setDouble(double *, int numArgs) {
    return -1;
}

int
DL_Interpreter::setString(const char*)
{
    return -1;
}
