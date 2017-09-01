#ifndef NURBS_H
#define NURBS_H


#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>


double g(double u,double i,std::vector<double> knots,int n) ;

double f(double u,double i,std::vector<double> knots, int n);

double N(double u,int i, int n, const std::vector<double> knots);

std::pair<double,double> getPointinNurbsCurve (double u,int n, const std::vector<std::pair<double,double>> controlPoints, const std::vector<double> knots,std::vector<double> weights);

void  Nurbs_Boundaries_getNurbsData (std::string filename);


#endif
