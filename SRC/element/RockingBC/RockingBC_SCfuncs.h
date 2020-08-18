/* ****************************************************************** **
**    RockingBC Shear calc functions  			    				  **
** ****************************************************************** */

// Written by: Evangelos Avgenakis

#ifndef RockingBC_SCfuncs_h
#define RockingBC_SCfuncs_h

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

#include <Matrix.h>
#include <Vector.h>

double inline SC_A(double x) { return 2.436222252877402 - 2.3818059387327604*x + 0.7078998718614156*x*x; };
double inline SC_DA(double x) { return -2.3818059387327604 + 2 * 0.7078998718614156*x; };

double inline SC_B(double x) { return 0.6982001887951753 - 1.098308073905204*x + 1.9266756798514126*x*x + -1.1270666845181774*x*x*x + 0.688867046041808*x*x*x*x; };
double inline SC_DB(double x) { return -1.098308073905204 + 2 * 1.9266756798514126*x + 3 * (-1.1270666845181774)*x*x + 4 * 0.688867046041808*x*x*x; };

double inline SC_C(double x) { return 0.8134604447686402*pow(1. - x, 3.770057533864266) + 1.; };
double inline SC_DC(double x) { return -0.8134604447686402*3.770057533864266*pow(1. - x, 2.770057533864266); };

double inline SC_D(double x) { return (1. - x)*(2.340417693163326 - 1.9592356132890616*x + 0.8914260492531663*x*x); };
double inline SC_DD(double x) { return -2.340417693163326 - (x - 1.)*(-1.9592356132890616 + 2 * 0.8914260492531663*x) + 1.9592356132890616*x - 0.8914260492531663*x*x; };

double inline SC_E(double x) { return 1.4043226196463283 + 0.1302424508017461*pow(1. - x, 3.6564163357661053) - 0.0549296131209048*x; };
double inline SC_DE(double x) { return -0.0549296131209048 - 0.1302424508017461*3.6564163357661053*pow(1. - x, 2.6564163357661053); };

double inline SC_F(double x) { return 0.4343458286281541*x + 3.107476490749382*x*x + (-6.967836976078876)*x*x*x + 6.501720103798543*x*x*x*x + (-2.284276614857206)*x*x*x*x*x; };
double inline SC_DF(double x) { return 0.4343458286281541 + 2 * 3.107476490749382*x + 3 * (-6.967836976078876)*x*x + 4 * 6.501720103798543*x*x*x + 5 * (-2.284276614857206)*x*x*x*x; };

void Dt_calc(const Vector& P, double& d, Vector& dddP);
void Rt_calc(const Vector& P, double& th, Vector& dthdP);

void se_shear_1der(const Vector& Youter, Vector& Ut, Matrix& dUt_dYouter);

#endif