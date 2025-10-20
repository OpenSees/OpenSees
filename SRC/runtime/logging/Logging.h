#pragma once
#include <handler/OPS_Stream.h>
#include <logging/AnsiColors.h>

class OPS_Stream;
namespace OpenSees {
   extern OPS_Stream *opserrPtr;
   extern OPS_Stream *opslogPtr;
   extern OPS_Stream *opswrnPtr;
   extern OPS_Stream *opsdbgPtr;
}


#define opsdbg (*OpenSees::opsdbgPtr)
#define opswrn ((*OpenSees::opswrnPtr) << G3_WARN_PROMPT)
#define opslog (*OpenSees::opslogPtr)
#ifndef opserr
#  define opserr (*OpenSees::opserrPtr)
#  define endln "\n"
#endif
#define LOG_TEST ":: "
#define LOG_ITERATE BLU "   ITERATE" COLOR_RESET " :: "
#define LOG_FAILURE RED "   FAILURE" COLOR_RESET " :: "
#define LOG_SUCCESS GRN "   SUCCESS" COLOR_RESET " :: "
#define LOG_CONTINUE "\n              "



class G3_Runtime;

extern const char *G3_WARN_PROMPT;
extern const char *G3_DEBUG_PROMPT;

namespace OpenSees {
// NOTE: These are defined in logging.cpp

  // maybe change "Prompt" to "Signal"
  extern const char * SignalMessageEnd;
  extern const char * PromptParseError;
  extern const char * PromptValueError;
  extern const char * PromptModelError;
  extern const char * SignalWarning;
  extern const char * PromptAnalysisFailure;
  extern const char * PromptAnalysisSuccess;
  extern const char * PromptAnalysisIterate;
};


enum G3_Stream {
  G3_StdOut, G3_StdIn, G3_StdErr, G3_Null
};

enum G3_StreamLevel {
  G3_LevelError, 
  G3_LevelDebug, 
  G3_LevelLog, 
  G3_LevelWarn
};

int G3_SetStreamColor(G3_Runtime* rt, int strm, int flag);
int G3_SetStreamLevel(int stream, bool on);

