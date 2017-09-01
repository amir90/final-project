#include "Nurbs.h"

struct NurbsCurve {
	std::vector<std::pair<double,double>> controlPoints;
	std::vector<double> knots;
};

typedef std::list<NurbsCurve> NurbsBoundary;

typedef std::list<NurbsBoundary> NurbsBoundaries; //limited to one surface with arbitrary number of holes

void  Nurbs_Boundaries_getNurbsData (std::string filename) {
//Recieves parasolid text datafile and outputs the 
}

double g(double u,double i,std::vector<double> knots,int n) {
	if (knots[i+n]!=knots[i]) {
	//	std::cout<<"g is: "<< (knots[i+n]-u)/(knots[i+n]-knots[i])<<"\n";
	return (knots[i+n]-u)/(knots[i+n]-knots[i]);
	} else {
	//	std::cout<<"g is: 0"<<"\n";
		return 0;
	}

}

double f(double u,double i,std::vector<double> knots, int n) {
	// u - position on curve to evaluate
	// i - knot number
	// knots - knot vector
	// n - degree
	//degree

	if (knots[i+n]!=knots[i]) {
	//std::cout<<"f is: "<< (u-knots[i])/(knots[i+n]-knots[i])<<"\n";
	return (u-knots[i])/(knots[i+n]-knots[i]);
	} else {
	//	std::cout<<"f is: 0"<<"\n";
		return 0;
	}

}

double N(double u,int i, int n, const std::vector<double> knots) {
	// u - position on curve to evaluate
	// i - current knot
	// n - degree
	if (n==0) {
		if ((knots[i]<=u && u<=knots[i+1])) {
		//	std::cout<<"here: "<<i<<" , "<<knots[i]<<" , "<<u<<" , "<<knots[i+1]<<" returned 1 \n";
			return 1;
		} else {
		//	std::cout<<"here: "<<i<<" , "<<knots[i]<<" , "<<u<<" , "<<knots[i+1]<<" returned 0 \n";
			return 0;
		}

	} else {

			auto test =  (f(u,i,knots,n)*N(u,i,n-1,knots)+g(u,i+1,knots,n)*N(u,i+1,n-1,knots));
		//	std::cout<<"what i am actually returning is: "<<test<<"\n";
			return test;

	}

}


std::pair<double,double> getPointinNurbsCurve (double u,int n, const std::vector<std::pair<double,double>> controlPoints, const std::vector<double> knots,std::vector<double> weights) {

	double denominator = 0;
	std::pair<double,double> nominator;

	nominator.first = 0; nominator.second =0;

	for (int j=0; j<=controlPoints.size()-1; j++) {

			auto curr_value = N(u,j, n,knots)*weights[j];
			denominator += curr_value;
		//	std::cout<<"j is: "<<j<<" denom is: "<<denominator<<" , "<<controlPoints[j].first<<"\n";
			nominator.first += curr_value*controlPoints[j].first;
			nominator.second +=  curr_value*controlPoints[j].second;


}

	nominator.first=nominator.first/denominator;
	nominator.second=nominator.second/denominator;

	return nominator;


}
