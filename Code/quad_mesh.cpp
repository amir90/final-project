#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
// workaround a bug in Boost-1.54
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <ipelet.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>
#include <boost/math/special_functions/binomial.hpp>
#include <math.h>
#include <CGAL/squared_distance_2.h>
typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_2                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Aff_transformation_2<Kernel> 						Transformation;
#define 	CGAL_PI   3.14159265358979323846


int main(int /* argc */, char* argv[])
{

	Transformation rotate(CGAL::ROTATION, std::sin(CGAL_PI/2), std::cos(CGAL_PI/2));

	//step 1: import an parasolid

	//step 2: find the boundary, describe it with a countinous curve

	// currently, not implemented

	//step 3:  create nodes on curve (how to best space out the nodes? - paper by talbert)

	//currently, use stub. rectangle with rectangular hole. specify bezier curve with control points


	struct bezier_curve
	 {
		std::list<Point> controlPoints;
	};


	std::list<bezier_curve> boundary;

	std::list<std::list<bezier_curve>> boundaries; //list which holds all the boundaries

	// ... input boundaries .... currently manual


	bezier_curve b1;

	b1.controlPoints.push_back(Point(-100,-100));
	b1.controlPoints.push_back(Point(100,-100));

	boundary.push_back(b1);

	b1.controlPoints.clear();

	b1.controlPoints.push_back(Point(100,-100));
	b1.controlPoints.push_back(Point(100,100));

	boundary.push_back(b1);

	b1.controlPoints.clear();

	b1.controlPoints.push_back(Point(100,100));
	b1.controlPoints.push_back(Point(-100,100));

	boundary.push_back(b1);

	b1.controlPoints.clear();

	b1.controlPoints.push_back(Point(-100,100));
	b1.controlPoints.push_back(Point(-100,-100));

	boundary.push_back(b1);

	boundaries.push_back(boundary);

	Mesh m;

	// begin meshing - populate boundaries with nodes

int max_nodes=6; //how many nodes per edge

	for (auto i = boundaries.begin(); i!=boundaries.end(); i++) { //for all boundaries
		for (auto j = i->begin(); j!=i->end(); j++) { //for all bezier curve in a boundary
		int n = j->controlPoints.size()-1;
			for (int t=1; t<=max_nodes; t++) {
				double seed=(double)t/max_nodes;
				double x=0;
				double y=0;
				int count=0;
			for (auto k = j->controlPoints.begin(); k!=j->controlPoints.end(); k++) { //for all points in the bezier curves
				double coeff =(double) pow(1-seed,n-count)*(double) pow(seed,count)*(double) boost::math::binomial_coefficient<double>(n,count);
				if (t==max_nodes && k==(--j->controlPoints.end())) {
					coeff = 1;
				}
				count++;
				x+=coeff*k->x();
				y+=coeff*k->y();
			}

			m.add_vertex(Point(x,y));

			}
		}
	}

	auto st = m.vertices_begin();
	auto fi = m.vertices_end();

	// begin "real algorithm"

	std::list<std::list<Mesh::Vertex_index>> sm_boundaries; //list of boundaries!

	std::list<Mesh::Vertex_index> sm_boundary; //populate boundary list with vertex index

	for (auto i = m.vertices_begin(); i!=m.vertices_end(); i++) {
		sm_boundary.push_back(*i); //getting correct vertex handle!
	}

	sm_boundaries.push_back(sm_boundary);

	//test sm_boundary

	//offset nodes according to paper


	//first node

	auto last_v = m.point(*(--sm_boundary.end()));

	auto  curr_v = m.point(*(sm_boundary.begin()));

	auto next_v =m.point( *(++sm_boundary.begin()));

	Kernel::Vector_2 vec1 = Kernel::Vector_2(curr_v,last_v);

	Kernel::Vector_2 vec2 = Kernel::Vector_2(curr_v,next_v);

	auto vec1_norm = vec1/std::sqrt(vec1.squared_length());

	auto vec2_norm = vec2/std::sqrt(vec2.squared_length());

	Kernel::RT D = (std::sqrt(vec1.squared_length())+std::sqrt(vec2.squared_length()))/2;

	Kernel::Vector_2 Theta;

	if (CGAL::collinear(last_v,curr_v,next_v)) {
		//rotate vec2 90 degrees counterclockwise (positive)
		std::cout<<"here\n";
		Theta = rotate(vec2_norm);

	} else {

		Theta =  vec1_norm+vec2_norm; //CGAL::angle(last_v,curr_v,next_v)/2; //calculateAngle(last_v,next_v);

	}
	auto new_v = curr_v + Theta*D;

	auto new_index = m.add_vertex(new_v);

	auto old_index = new_index;

	//loop over other nodes (Except the last node)

	for (auto i=(++sm_boundary.begin()); i!=(--sm_boundary.end()) ;i++) {


		last_v = curr_v;

		curr_v = m.point(*i);

		auto j = i;

		j++;

		next_v = m.point(*(j));

		vec1 = Kernel::Vector_2(curr_v,last_v);

		vec2 = Kernel::Vector_2(curr_v,next_v);

		vec1_norm = vec1/std::sqrt(vec1.squared_length());

		vec2_norm = vec2/std::sqrt(vec2.squared_length());

		D = (std::sqrt(vec1.squared_length())+std::sqrt(vec2.squared_length()))/2;

		if (CGAL::collinear(last_v,curr_v,next_v)) {
			//rotate vec2 90 degrees counterclockwise (positive)
			Theta = rotate(vec2_norm);

		} else {

			Theta =  vec1_norm+vec2_norm; //CGAL::angle(last_v,curr_v,next_v)/2; //calculateAngle(last_v,next_v);

		}


		new_v = curr_v + Theta*D;

		new_index = m.add_vertex(new_v);

		j=i;

		j--;

		m.add_face(*i,*j,old_index,new_index);

		old_index = new_index;

	}

	//TODO: handle last node (like first node)

	last_v = curr_v;

	curr_v = m.point(*(--sm_boundary.end()));

	next_v =m.point( *(sm_boundary.begin()));

	vec1 = Kernel::Vector_2(curr_v,last_v);

	vec2 = Kernel::Vector_2(curr_v,next_v);

	vec1_norm = vec1/std::sqrt(vec1.squared_length());

	vec2_norm = vec2/std::sqrt(vec2.squared_length());

	D = (std::sqrt(vec1.squared_length())+std::sqrt(vec2.squared_length()))/2;

	if (CGAL::collinear(last_v,curr_v,next_v)) {
		//rotate vec2 90 degrees counterclockwise (positive)
		Theta = rotate(vec2_norm);

	} else {

		Theta =  vec1_norm+vec2_norm; //CGAL::angle(last_v,curr_v,next_v)/2; //calculateAngle(last_v,next_v);

	}
	new_v = curr_v + Theta*D;

	new_index = m.add_vertex(new_v);

	auto j=--sm_boundary.end();

	j--;

	m.add_face(*(--sm_boundary.end()),*j,old_index,new_index);

	//TODO: consider a circular list (better implementation)




//testing node placement using ipe
std::ofstream myFile;
myFile.open("Ipe.xml");
myFile << "<page>\n";
for (auto i = m.vertices().begin(); i!=m.vertices().end(); i++) {
//	std::cout <<m.point(*i);
//	std::cout<<"\n";


myFile << "<use name=\"mark/disk(sx)\" " << "pos= \"" << m.point(*i).x() << " " << m.point(*i).y() << "\" size=\"normal\" stroke=\"black\"/>\n";



}
for (auto i = m.edges().begin(); i!=m.edges().end(); i++) {
//	std::cout <<m.point(*i);
//	std::cout<<"\n";

Point p1 = m.point((m.vertex(*i,0)));

Point p2 = m.point((m.vertex(*i,1)));

myFile << "<path stroke = \"black\"> \n"  << p1 <<" m \n" << p2 << " l \n" << "</path> \n";


}
myFile << "</page>\n";
myFile.close();




}


