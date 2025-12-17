// CatenaryCableMaxsag.cpp 的最小头文件集：

#include <OPS_Globals.h>      // 必须：opserr
#include <Domain.h>           // 必须：Domain 类
#include <Node.h>             // 必须：Node 类、getCrds()
#include <CatenaryCable.h>    // 必须：CatenaryCable 类
#include <Vector.h>           // 可能需要：Vector 类

#include <string>             // std::string
#include <algorithm>          // std::transform
#include <cctype>             // ::tolower
#include <math.h>             // sqrt 等#pragma once
// CatenaryCableMaxsag.h
#ifndef CATENARY_CABLE_MAXSAG_H
#define CATENARY_CABLE_MAXSAG_H

int processCatenaryWithMaxsag(
    Domain* theDomain,  // 添加这个参数，直接传入 Domain
    int eleTag, int nodeI, int nodeJ,
    double weight, double E, double A,
    double* pL0,
    double alpha, double tempChange, double rho, double errorTol);

#endif