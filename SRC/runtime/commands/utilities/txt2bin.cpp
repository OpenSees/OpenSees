#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

using std::cerr;
using std::ios;
using std::setiosflags;
using std::ifstream;
using std::ofstream;
using std::string;
using std::getline;

int 
binaryToText(const char *inputFilename, const char *outputFilename, int prec = 0)
{
  //
  // open the files
  //
  ifstream input(inputFilename, ios::in | ios::binary);
  ofstream output(outputFilename, ios::out);
  if (input.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << inputFilename << std::endl;
    output.close();
    return -1;
  }
  if (output.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << outputFilename << std::endl;
    output.close();
    return -1;
  }

  //
  // until done
  //   read input consisting doubles till \n and write to output file
  //
  double data;
  char *c = (char *)&data;
  int numNumbers = 0;
  /* ORIGINAL */
  while ( !input.eof()) {
    input.read(c, 1);
    if (*c == '\r') {
      std::cerr << "\nc: " << c << "\n";
      input.read(c, 1);
    }
    if (*c != '\n') {
      input.read(&c[1],7);
      output << data << " ";
      numNumbers++;
    } else {
      std::cerr << c;
      if (numNumbers != 0) {
        std::cerr << "\nc: " << c << "\n";
        output << "\n";
      }
      numNumbers = 0;
    }
  }
  output << "\n";
  // 
  // close the files
  //

  input.close();
  output.close();

  return 0;
}

int 
binaryToText2(const char *inputFilename, const char *outputFilename, int prec = 0)
{
  //
  // open the files
  //
  ifstream input(inputFilename, ios::in | ios::binary);
  ofstream output(outputFilename, ios::out);
  if (input.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << inputFilename << std::endl;
    output.close();
    return -1;
  }
  if (output.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << outputFilename << std::endl;
    output.close();
    return -1;
  }

  //
  // until done
  //   read input consisting doubles till \n and write to output file
  //

  double data;
  char *c = (char *)&data;
  int numNumbers = 0;
  /* ORIGINAL
  while ( !input.eof()) {
    input.read(c, 1);
    if (*c == '\r') {
      std::cerr << c;
      input.read(c, 1);
    }
    if (*c != '\n') {
      input.read(&c[1],7);
      output << data << " ";
      numNumbers++;
    } else {
      std::cerr << c;
      if (numNumbers != 0)
        output << "\n";
      numNumbers = 0;
    }
  }
  REPLACED WITH: */

  int nLF = 0;
  int intervalLF = 0;
  double aveIntervalLF = 0;
  int numCol;

  int dataLen = 0;
  int numRow;

  while (!input.eof()) {
    input.read(c, 1);
    dataLen++;

    if ((*c == '\n') && (intervalLF%8 == 0)) {    
      aveIntervalLF = (aveIntervalLF * nLF + intervalLF)/(nLF + 1);
      nLF++;
      intervalLF = 0;
    } else {
      intervalLF++;      
    }

  }
  numCol = (int)(aveIntervalLF/8.0 + 0.5);
  numRow = (int)(dataLen/(8.0*numCol+1.0) + 0.5);
  std::cerr << "\nnr: " << numRow << "\nnc: " << numCol << "\n";

  // rewind
  input.clear();
  input.seekg(0,std::ios::beg); 

  // read, output
  for (int ii=0; ii<numRow; ii++) {
    for (int jj=0; jj<numCol; jj++) {

      input.read(&c[0],8);
      output << data ;
      if (jj<(numCol-1)) output << " ";

    }
    input.read(c, 1);//LF
    output << "\n";

  }
  /* END REPLACEMENT BLOCK */

  // 
  // close the files
  //

  input.close();
  output.close();

  return 0;
}

int 
textToBinary(const char *inputFilename, const char *outputFilename)
{
  //
  // open the files
  //

  ifstream input(inputFilename, ios::in);
  ofstream output(outputFilename, ios::out | ios::binary);
  if (input.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << inputFilename << std::endl;
    output.close();
    return -1;
  }
  if (output.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << outputFilename << std::endl;
    output.close();
    return -1;
  }

  //
  // until done
  //   read input consisting doubles till \n and write to output file
  //

  char data[100];
  char *dataNext;
  double d;

  while (!input.eof()) {
    string inputLine;
    getline(input, inputLine);
    const char *c = inputLine.data();
    const char *cNext = c;

    int loc = 0;
    int endLoc = int(inputLine.length());
    int numNumbers = 0;

    while (loc < endLoc) {
      int dataCount = 0;
      while ((loc < endLoc) && 
             (*cNext != ' ') && 
             (*cNext != '\n')) {
        data[dataCount++] = cNext[0];
        cNext++;
        loc++;
      }
      
      if (dataCount != 0) {
        data[dataCount] = '\n';
        d = strtod(&data[0], &dataNext);
        output.write((char *)&d, 8);
        numNumbers++;
      }
      
      cNext++;
      loc++;
    }
    
    if (numNumbers != 0)
      output << '\n';
  }
  
  // 
  // close the files
  //

  input.close();
  output.close();

  return 0;
}


#include "string.h"
const char *file_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

int main(int argc, char *argv[])
{
  int argi = 0;
  std::cerr << "argc = " << argc << "\n\t" << argv[0];
  if (strcmp(file_ext(argv[++argi]), "txt")==0) {
    textToBinary(argv[argi], argv[argi+1]);
  } else {
    if (argc == 3)
      binaryToText2(argv[argi], argv[argi+1]);
    else if (strcmp(argv[argi+2], "-b2")==0)
      binaryToText2(argv[argi], argv[argi+1]);
    else if (strcmp(argv[argi+2], "-b1")==0)
      binaryToText(argv[argi], argv[argi+1]);
  }
  
}

