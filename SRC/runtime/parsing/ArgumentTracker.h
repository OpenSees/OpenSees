//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//                   This file is not yet licensed.
//
//===----------------------------------------------------------------------===//
//
// Description: ArgumentTracker is a template class designed to aid parsing
// command arguments.
//
// Written: cmp
//
#pragma once
#include <set>
#include <array>
#include <vector>

template <class Enum>
class ArgumentTracker {
public:
  template<int n>
  ArgumentTracker(std::array<Enum, n>&p) {
    for (Enum i : p )
      positions.push_back(i);
  }

  ArgumentTracker() {
    int i = 0;
    while (static_cast<Enum>(i) != Enum::End) {
      positions.push_back(static_cast<Enum>(i));
      i++;
    }
  }

  Enum current() const {
    for (Enum e : positions)
      if (consumed.find(e) == consumed.end())
        return e;

    return Enum::End;
  }

  bool contains(Enum arg) const {
    return consumed.find(arg) == consumed.end();
  }

  void consume(Enum argument) {
    consumed.insert(argument);
  }

  void increment() {
    consumed.insert(current());
  }

  template<class T>
  void print(T& strm) const {
    for (auto i : consumed) {
      strm << static_cast<int>(i) << "\n";
    }
  }


private:
  std::set<Enum>    consumed;
  std::vector<Enum> positions;
};

