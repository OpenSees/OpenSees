//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file provides basic file format handling commands,
// such as naive XML processing and binary conversion.
//
#include <tcl.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <OPS_Globals.h>

extern int binaryToText(const char *inputFile, const char *outputFile);
extern int textToBinary(const char *inputFile, const char *outputFile);

int
convertBinaryToText(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv)
{
  if (argc < 3) {
    opserr << "ERROR incorrect # args - convertBinaryToText inputFile "
              "outputFile\n";
    return -1;
  }

  const char *inputFile = argv[1];
  const char *outputFile = argv[2];

  return binaryToText(inputFile, outputFile);
}

int
convertTextToBinary(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv)
{
  if (argc < 3) {
    opserr << "ERROR incorrect # args - convertTextToBinary inputFile "
              "outputFile\n";
    return -1;
  }

  const char *inputFile = argv[1];
  const char *outputFile = argv[2];

  return textToBinary(inputFile, outputFile);
}

int
stripOpenSeesXML(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{

  if (argc < 3) {
    opserr << "ERROR incorrect # args - stripXML input.xml output.dat "
              "<output.xml>\n";
    return -1;
  }

  const char *inputFile = argv[1];
  const char *outputDataFile = argv[2];
  const char *outputDescriptiveFile = 0;

  if (argc == 4)
    outputDescriptiveFile = argv[3];

  // open files
  std::ifstream theInputFile;
  theInputFile.open(inputFile, std::ios::in);
  if (theInputFile.bad()) {
    opserr << "stripXML - error opening input file: " << inputFile << endln;
    return -1;
  }

  std::ofstream theOutputDataFile;
  theOutputDataFile.open(outputDataFile, std::ios::out);
  if (theOutputDataFile.bad()) {
    opserr << "stripXML - error opening input file: " << outputDataFile
           << endln;
    return -1;
  }

  std::ofstream theOutputDescriptiveFile;
  if (outputDescriptiveFile != 0) {
    theOutputDescriptiveFile.open(outputDescriptiveFile, std::ios::out);
    if (theOutputDescriptiveFile.bad()) {
      opserr << "stripXML - error opening input file: " << outputDescriptiveFile
             << endln;
      return -1;
    }
  }

  std::string line;
  bool spitData = false;
  while (!theInputFile.eof()) {
    getline(theInputFile, line);
    const char *inputLine = line.c_str();

    if (spitData == true) {
      if (strstr(inputLine, "</Data>") != 0)
        spitData = false;
      else {
        ; // theOutputDataFile << line << endln;
      }
    } else {
      const char *inputLine = line.c_str();
      if (strstr(inputLine, "<Data>") != 0)
        spitData = true;
      else if (outputDescriptiveFile != nullptr) {
        ; // theOutputDescriptiveFile << line << endln;
      }
    }
  }

  theInputFile.close();
  theOutputDataFile.close();

  if (outputDescriptiveFile != nullptr)
    theOutputDescriptiveFile.close();

  return 0;
}
