#include "Nurbs.h"

struct NurbsCurve {
	std::vector<std::pair<double,double>> controlPoints;
	std::vector<double> knots;
	std::vector<int> knotMultiplicity;
	bool sense;
};

typedef std::list<NurbsCurve> NurbsBoundary;

typedef std::list<NurbsBoundary> NurbsBoundaries; //limited to one surface with arbitrary number of holes

void  Nurbs_Boundaries_getNurbsData (std::string filename) {
//Recieves parasolid text datafile and outputs the 
}

double g(double u,double i,std::vector<double> knots,int n) {
	std::cout<<"g is: "<< (knots[i+n]-u)/(knots[i+n]-knots[i])<<"\n";
	return (knots[i+n]-u)/(knots[i+n]-knots[i]);

}

double f(double u,double i,std::vector<double> knots, int n) {
	std::cout<<"f is: "<< (u-knots[i])/(knots[i+n]-knots[i])<<"\n";
	return (u-knots[i])/(knots[i+n]-knots[i]);

}

double N(double u,int i,std::pair<double,double> controlPoint, int n, const std::vector<double> knots) {

	if (n==0) {

		if (i<=u && u<=i+1) {
			return 1;
		} else {
			return 0;
		}

	} else {
	return (f(u,i,knots,n)*N(u,i,controlPoint,n-1,knots)+g(u,i,knots,n)*N(u,i-1,controlPoint,n-1,knots));

	}

}



std::pair<double,double> getPointinNurbsCurve (double u,int n, const std::vector<std::pair<double,double>> controlPoints, const std::vector<double> knots, const std::vector<int> KnotMultiplicity,std::vector<double> weights) {

	double denominator = 0;
	std::pair<double,double> nominator;

	nominator.first = 0; nominator.second =0;

	for (int j=0; j<=controlPoints.size()-1; j++) {

			denominator += N(u,j,controlPoints[j], n,knots)*weights[j];
			std::cout<<"j is: "<<j<<" denom is: "<<denominator<<"\n";
			nominator.first += denominator*controlPoints[j].first;
			nominator.second += denominator*controlPoints[j].second;


}

	return nominator;


}
