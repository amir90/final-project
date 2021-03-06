#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
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
typedef CGAL::Simple_cartesian<double>                     Kernel;
typedef Kernel::Point_2                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef Kernel::Segment_2									Segment;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Aff_transformation_2<Kernel> 						Transformation;
typedef std::list<Mesh::Vertex_index>  						Boundary;
typedef std::list<Mesh::Vertex_index>::iterator  BoundaryIterator;
#define PI 3.14159
#define beta 0.226893
#define	Crit1_lower 0.785398 //in radians
#define Crit1_higher 2.35619 //in radians
#define Min_Collision_Value 0 //cos of angle in radians
#define Relaxation_Factor 1.3

Point getOffsetPoint (Point last_v, Point curr_v, Point next_v) {
	Kernel::Vector_2 vec1 = Kernel::Vector_2(curr_v,last_v);

	Kernel::Vector_2 vec2 = Kernel::Vector_2(curr_v,next_v);

	auto vec1_norm = vec1/Kernel::FT(std::sqrt(CGAL::to_double(vec1.squared_length())));

	auto vec2_norm = vec2/Kernel::FT(std::sqrt(CGAL::to_double(vec2.squared_length())));

	Kernel::FT D = Kernel::FT((std::sqrt(CGAL::to_double(vec1.squared_length()))+std::sqrt(CGAL::to_double(vec2.squared_length())))/2);
	Kernel::Vector_2 Theta;

	Kernel::FT Angle;


	if (CGAL::left_turn(last_v,curr_v,next_v)) {
		if ((vec1_norm*vec2_norm)<-1) {
			Angle=PI;
		} else {
			Angle = acos(vec1_norm*vec2_norm);
		}
	} else {
		if ((vec1_norm*vec2_norm)<-1) {
			Angle=PI;
		} else
			Angle = 2*PI-acos(vec1_norm*vec2_norm);
		}

	Transformation rotate2(CGAL::ROTATION, std::sin(Angle/2), std::cos(Angle/2));

	Theta = rotate2(vec2_norm);
	auto new_v = curr_v + Theta*D;

	if (new_v.x()!=new_v.x()) {
		std::cout << "Found NAN value\n";

	}

	return new_v;
}

BoundaryIterator getNext (BoundaryIterator begin,BoundaryIterator end,BoundaryIterator curr, int num) {
	if (num>0) {
		while (num>0) {
			num--;
			curr == std::next(end,-1)? curr=begin:curr=std::next(curr,1);
			}
	} else if (num<0) {
		while (num<0) {
			num++;
			curr == begin? curr=std::next(end,-1):curr=std::next(curr,-1);
		}
	}
	return curr;
}

Mesh getTestMesh() {
	Mesh m;

	//Hexagon
/*
	double h = std::sqrt(6*6-3*3);
	m.add_vertex(Point(3,h));
	m.add_vertex(Point(-3,h));
	m.add_vertex(Point(-6,0));
	m.add_vertex(Point(-3,-h));
	m.add_vertex(Point(3,-h));
	m.add_vertex(Point(6,0));
	*/


	//Convex six-loop with one 180-degree angle
	/*
	m.add_vertex(Point(0,20));
	m.add_vertex(Point(10,10));
	m.add_vertex(Point(10,0));
	m.add_vertex(Point(0,0));
	m.add_vertex(Point(-10,0));
	m.add_vertex(Point(-10,10));
	*/


	//Convex six-loop with two 180-degree angle - case 1
	/*
	m.add_vertex(Point(-20,20));
	m.add_vertex(Point(-20,0));
	m.add_vertex(Point(0,0));
	m.add_vertex(Point(20,0));
	m.add_vertex(Point(20,20));
	m.add_vertex(Point(0,20));
/*
/*
	//Convex six-loop with two 180-degree angle - case 2
	m.add_vertex(Point(-30,20));
	m.add_vertex(Point(-30,0));
	m.add_vertex(Point(20,0));
	m.add_vertex(Point(20,15));
	m.add_vertex(Point(20,20));
	m.add_vertex(Point(0,20));
	*/

/*
	//Convex six-loop with two 180-degree angle - case 3
	m.add_vertex(Point(-30,20));
	m.add_vertex(Point(-15,0));
	m.add_vertex(Point(15,0));
	m.add_vertex(Point(30,20));
	m.add_vertex(Point(15,20));
	m.add_vertex(Point(-15,20));
	*/

/*
	//Convex six-loop with three 180-degree angle - case 1
	m.add_vertex(Point(0,0));
	m.add_vertex(Point(10,10));
	m.add_vertex(Point(20,20));
	m.add_vertex(Point(30,30));
	m.add_vertex(Point(30,0));
	m.add_vertex(Point(10,0));
	*/

	//Convex six-loop with three 180-degree angle - case 2

/*
	m.add_vertex(Point(30,0));
	m.add_vertex(Point(15,15));
	m.add_vertex(Point(0,30));
	m.add_vertex(Point(-15,15));
	m.add_vertex(Point(-30,0));
	m.add_vertex(Point(0,0));
*/

	m.add_vertex(Point)


	/*
	//Convex six-loop with three 180-degree angle - case 4
	m.add_vertex(Point(-30,0));
	m.add_vertex(Point(-15,0));
	m.add_vertex(Point(0,0));
	m.add_vertex(Point(15,0));
	m.add_vertex(Point(30,0));
	m.add_vertex(Point(0,60));
	*/

	//Testing Bad element removal! - Case 1
/*
	auto e1_v1 = m.add_vertex(Point(0,20));
	auto e1_v2 = m.add_vertex(Point(-20,0));
	auto e1_v3 = m.add_vertex(Point(0,-20));
	auto e1_v4 = m.add_vertex(Point(20,0));
	auto e1_v5 = m.add_vertex(Point(0,0));

	m.add_face(e1_v1,e1_v2,e1_v3,e1_v5);
	m.add_face(e1_v3,e1_v4,e1_v1,e1_v5);
*/

	//Testing Bad element removal! - Case 2
/*
	auto e1_v1 = m.add_vertex(Point(0,0));
	auto e1_v2 = m.add_vertex(Point(15,15));
	auto e1_v3 = m.add_vertex(Point(15,30));
	auto e1_v4 = m.add_vertex(Point(0,30));
	m.add_face (e1_v4,e1_v3,e1_v2,e1_v1);

	auto e2_v1 = m.add_vertex(Point(30,30));
	auto e2_v4 = m.add_vertex(Point(30,0));
	m.add_face (e1_v3,e2_v1,e2_v4,e1_v2);

	auto e3_v4 = m.add_vertex(Point(15,-15));
	m.add_face(e1_v1,e1_v2,e2_v4,e3_v4);

	auto e4_v1 = m.add_vertex(Point(0,-30));
	auto e4_v2 = m.add_vertex(Point(15,-30));
	m.add_face(e1_v1,e3_v4,e4_v2,e4_v1);

	auto e5_v1 = m.add_vertex(Point(30,-30));
	m.add_face(e3_v4,e2_v4,e5_v1,e4_v2);

*/

	//Case 3
/*
	auto v1 = m.add_vertex(Point(0,20));
	auto v2 = m.add_vertex(Point(10,10));
	auto v3 = m.add_vertex(Point(10,0));
	auto v4 = m.add_vertex(Point(0,0));
	auto v5 = m.add_vertex(Point(-10,0));
	auto v6 = m.add_vertex(Point(-10,10));
	auto v7=m.add_vertex(Point(0,5));
	m.add_face(v1,v2,v3,v7);
	m.add_face(v3,v4,v5,v7);
	m.add_face(v5,v6,v1,v7);
	*/


	//Case 4
/*
	double h = std::sqrt(6*6-3*3);
	auto v1 =m.add_vertex(Point(3,h));
	auto v2 =m.add_vertex(Point(-3,h));
	auto v3 =m.add_vertex(Point(-6,0));
	auto v4 =m.add_vertex(Point(-3,-h));
	auto v5 =m.add_vertex(Point(3,-h));
	auto v6 =m.add_vertex(Point(6,0));
	auto v7 =m.add_vertex(Point(2,0));
	auto v8 =m.add_vertex(Point(-2,0));
	m.add_face(v1,v2,v8,v7);
	m.add_face(v2,v3,v4,v8);
	m.add_face(v4,v5,v7,v8);
	m.add_face(v5,v6,v1,v7);

*/
	return m;
}


bool isCollinear (Point last_v, Point curr_v, Point next_v) {

	auto vec1 = curr_v-last_v;

	auto vec2 = next_v-curr_v;

	auto vec1_norm = vec1/Kernel::FT(std::sqrt(CGAL::to_double(vec1.squared_length())));

	auto vec2_norm = vec2/Kernel::FT(std::sqrt(CGAL::to_double(vec2.squared_length())));

	if (acos(vec1_norm*vec2_norm)<=beta) {
		return true;
	} else {
		return false;
	}
}

//TODO: test to make sure works for every starting position
//TODO: define *beta* - according to finalProject.pdf to decide if point is colinear
void applySixNodeSplitter (Boundary boundary, Mesh& m) {


	std::cout<<"Applying 6 node splitter\n";

	//total number of collinear vertices
	int numOfCollinearVertex=0;
	//temporary number of consecutive collinear vertices
	int temp=0;
	//maximum number of consecutive collinear vertices
	int maxNumOfConsecutiveCollinearVertex=0;
	//number of non-collinear vertices between current most forward collinear vertex and next collinear vertex
	int numOfNormalVertices = 0;
	// save current numOfNormalVertices
	int SaveNumOfNormalVertices = 0;
	//assume splitter vertex is at beginning of boundary
	auto splitterVertex = (boundary.begin());
	//j is pointer to current vertex
	auto j=boundary.begin();
	bool normalVerticesFlag=false;

	//go *TWICE* around the boundary - makes starting vertex of boundary irrelevent
	for (int i=1; i<=2*boundary.size(); i++) {

		auto next_vertex = getNext(boundary.begin(),boundary.end(),j,1);
		auto last_vertex = getNext(boundary.begin(),boundary.end(),j,-1);

			if (isCollinear(m.point(*last_vertex),m.point(*j),m.point(*next_vertex))) {
			numOfCollinearVertex++; //found another collinear vertex
			temp++;
			if (normalVerticesFlag) {
			SaveNumOfNormalVertices = numOfNormalVertices;
			normalVerticesFlag =false;
			}
			numOfNormalVertices=0;
			} else {
				if (temp>maxNumOfConsecutiveCollinearVertex) {
					maxNumOfConsecutiveCollinearVertex=temp;
					splitterVertex = last_vertex;
					normalVerticesFlag =true; //start keeping count of distance to next collinear vertex
				}
			numOfNormalVertices++;
			temp=0;
			}

		j=next_vertex;
	}

	numOfCollinearVertex=numOfCollinearVertex/2;

	BoundaryIterator si = boundary.begin();
	BoundaryIterator ei = boundary.end();

	std::cout<<"NumberOfNormalVertices: " << SaveNumOfNormalVertices << "\nNumber of consecutive colinear vertices: " << maxNumOfConsecutiveCollinearVertex << "\nNumber of colinear vertices: " <<numOfCollinearVertex <<"\n";

	if (numOfCollinearVertex==0) {
		if (splitterVertex==si) {
		}
		m.add_face(*boundary.begin(),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),1),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),2),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),3));
		m.add_face(*getNext(boundary.begin(),boundary.end(),(boundary.begin()),3),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),4),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),5),*boundary.begin());

	}
	if (numOfCollinearVertex==1) {
		std::cout<<"Applied six-node splitter - Begin \n";
		std::cout<<m.point(*getNext(si,ei,splitterVertex,3)) <<" , "<<m.point(*getNext(si,ei,splitterVertex,2))<<" , "<<m.point(*getNext(si,ei,splitterVertex,1))<<" , "<<m.point(*getNext(si,ei,splitterVertex,0))<<"\n";
		std::cout<<m.point(*getNext(si,ei,splitterVertex,5)) <<" , "<<m.point(*getNext(si,ei,splitterVertex,4))<<" , "<<m.point(*getNext(si,ei,splitterVertex,3))<<" , "<<m.point(*getNext(si,ei,splitterVertex,0))<<"\n";
		//	m.add_face(*boundary.begin(),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),1),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),2),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),3));
		//	m.add_face(*getNext(boundary.begin(),boundary.end(),(boundary.begin()),3),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),4),*getNext(boundary.begin(),boundary.end(),(boundary.begin()),5),*boundary.begin());
		//	m.add_face (*getNext(si,ei,splitterVertex,3),*getNext(si,ei,splitterVertex,2),*getNext(si,ei,splitterVertex,1),*splitterVertex);
		//	m.add_face (*splitterVertex,*getNext(si,ei,splitterVertex,5),*getNext(si,ei,splitterVertex,4),*getNext(si,ei,splitterVertex,3));
		m.add_face (*getNext(si,ei,splitterVertex,0),*getNext(si,ei,splitterVertex,1),*getNext(si,ei,splitterVertex,2),*getNext(si,ei,splitterVertex,3));
		m.add_face (*getNext(si,ei,splitterVertex,3),*getNext(si,ei,splitterVertex,4),*getNext(si,ei,splitterVertex,5),*getNext(si,ei,splitterVertex,0));

		std::cout<<"Applied six-node splitter - End \n";
	}
	if (numOfCollinearVertex==2 && maxNumOfConsecutiveCollinearVertex==1) {
		if (SaveNumOfNormalVertices==1) {
			auto lastCollinearVertex = getNext(si,ei,splitterVertex,-2);
			auto originVertex = getNext(si,ei,lastCollinearVertex,1);
			auto newVertex = m.add_vertex(m.point(*originVertex) + (m.point(*lastCollinearVertex)-m.point(*originVertex))+(m.point(*splitterVertex)-m.point(*originVertex)));
			std::cout <<"\nhere "<<m.point(*splitterVertex)<<" , "<<m.point(*originVertex) << " , "<<m.point(newVertex)<<"\n";
			m.add_face (*splitterVertex,*originVertex,*lastCollinearVertex,newVertex);
			auto v1 = getNext(si,ei,splitterVertex,2);
			m.add_face(*v1,*getNext(si,ei,v1,-1),*splitterVertex,newVertex);
			m.add_face(*v1,newVertex,*lastCollinearVertex,*getNext(si,ei,v1,1));
		} else if (SaveNumOfNormalVertices==3) {
			auto lastCollinearVertex = getNext(si,ei,splitterVertex,-3);
			auto originVertex = getNext(si,ei,lastCollinearVertex,-1);
			auto newVertex = m.add_vertex(m.point(*originVertex) + (m.point(*lastCollinearVertex)-m.point(*originVertex))+(m.point(*splitterVertex)-m.point(*originVertex)));
			m.add_face (*splitterVertex,*originVertex,*lastCollinearVertex,newVertex);
			auto v1 = getNext(si,ei,splitterVertex,-2);
			m.add_face(*v1,*getNext(si,ei,v1,-1),*splitterVertex,newVertex);
			m.add_face(*v1,newVertex,*lastCollinearVertex,*getNext(si,ei,v1,1));
		} else {
	//		std::cout <<"\nhere "<< m.point(*splitterVertex)<<"\n";
			m.add_face (*getNext(si,ei,splitterVertex,3),*getNext(si,ei,splitterVertex,2),*getNext(si,ei,splitterVertex,1),*splitterVertex);
			m.add_face (*splitterVertex,*getNext(si,ei,splitterVertex,5),*getNext(si,ei,splitterVertex,4),*getNext(si,ei,splitterVertex,3));
				}
	}
	if (numOfCollinearVertex==2 && maxNumOfConsecutiveCollinearVertex==2) {
		auto lastCollinearVertex = getNext(si,ei,splitterVertex,-1);
		auto v1 = getNext(si,ei,splitterVertex,2);
		auto newV1 = m.add_vertex(m.point(*v1) + (m.point(*getNext(si,ei,v1,1))-m.point(*v1))+(m.point(*getNext(si,ei,v1,-1))-m.point(*v1)));
		auto v2 = getNext(si,ei,lastCollinearVertex,-2);
		auto newV2 = m.add_vertex(m.point(*v2) + (m.point(*getNext(si,ei,v2,1))-m.point(*v2))+(m.point(*getNext(si,ei,v2,-1))-m.point(*v2)));
		m.add_face (*splitterVertex, *lastCollinearVertex, newV2, newV1);
		m.add_face (*getNext(si,ei,splitterVertex,1),*splitterVertex,newV1,*v1);
		m.add_face(*lastCollinearVertex,*getNext(si,ei,lastCollinearVertex,-1),*v2,newV2);
		m.add_face (newV1,newV2,*v2,*v1);
	}
	if (numOfCollinearVertex==3 && maxNumOfConsecutiveCollinearVertex==3) {
		auto newV1 = m.add_vertex(m.point(*getNext(si,ei,splitterVertex,-2))+(m.point(*getNext(si,ei,splitterVertex,2)) - (m.point(*getNext(si,ei,splitterVertex,-2))))*2/3);
		auto newV2 = m.add_vertex(m.point(*getNext(si,ei,splitterVertex,0))+(m.point(*getNext(si,ei,splitterVertex,2)) - (m.point(*getNext(si,ei,splitterVertex,0))))*2/3);
		auto newV3 = m.add_vertex(m.point(*getNext(si,ei,splitterVertex,-1))+(m.point(*getNext(si,ei,splitterVertex,2)) - (m.point(*getNext(si,ei,splitterVertex,-1))))*1/3);

//		std::cout<<"\nhere , "<<m.point(*getNext(si,ei,splitterVertex,2))<<" , "<<m.point(*getNext(si,ei,splitterVertex,-2))<<" , "<<(m.point(*getNext(si,ei,splitterVertex,2)) - (m.point(*getNext(si,ei,splitterVertex,-2))))*2/3<<" , "<<m.point(newV1) << " , " <<m.point(newV2)<< " ," <<m.point(newV3)<<"\n";

		m.add_face(*splitterVertex,*getNext(si,ei,splitterVertex,-1),newV3,newV2);
		m.add_face(*getNext(si,ei,splitterVertex,-1),*getNext(si,ei,splitterVertex,-2),newV1,newV3);
		m.add_face(*getNext(si,ei,splitterVertex,-2),*getNext(si,ei,splitterVertex,-3),*getNext(si,ei,splitterVertex,-4),newV1);
		m.add_face(newV2,*getNext(si,ei,splitterVertex,2),*getNext(si,ei,splitterVertex,1),*splitterVertex);
		m.add_face(*getNext(si,ei,splitterVertex,2),newV2,newV3,newV1);


	}
	if (numOfCollinearVertex==3 && maxNumOfConsecutiveCollinearVertex==2) {
		m.add_face (*getNext(si,ei,splitterVertex,3),*getNext(si,ei,splitterVertex,2),*getNext(si,ei,splitterVertex,1),*splitterVertex);
								m.add_face (*splitterVertex,*getNext(si,ei,splitterVertex,5),*getNext(si,ei,splitterVertex,4),*getNext(si,ei,splitterVertex,3));

	}
	if (numOfCollinearVertex==3 && maxNumOfConsecutiveCollinearVertex==1) {
			auto NewV = m.add_vertex(CGAL::ORIGIN+((m.point(*getNext(si,ei,splitterVertex,0))-CGAL::ORIGIN)+(m.point(*getNext(si,ei,splitterVertex,2))-CGAL::ORIGIN)+(m.point(*getNext(si,ei,splitterVertex,4))-CGAL::ORIGIN))/3);
					m.add_face(*splitterVertex,*getNext(si,ei,splitterVertex,1),*getNext(si,ei,splitterVertex,2),NewV);
					m.add_face(*getNext(si,ei,splitterVertex,2),*getNext(si,ei,splitterVertex,3),*getNext(si,ei,splitterVertex,4),NewV);
					m.add_face(*getNext(si,ei,splitterVertex,4),*getNext(si,ei,splitterVertex,5),*splitterVertex,NewV);

	}
}

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
	b1.controlPoints.push_back(Point(50, -200));
	b1.controlPoints.push_back(Point(100,-100));

	boundary.push_back(b1);

	b1.controlPoints.clear();
	b1.controlPoints.push_back(Point(100,-100));
	b1.controlPoints.push_back(Point(150, 50));
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

int max_nodes=25; //how many nodes per edge

	for (auto i = boundaries.begin(); i!=boundaries.end(); i++) { //for all boundaries
		for (auto j = i->begin(); j!=i->end(); j++) { //for all bezier curve in a boundary
		int n = j->controlPoints.size()-1;
			for (int t=1; t<=(max_nodes); t++) {
				double seed=(double)t/(max_nodes);
				double x=0;
				double y=0;
				int count=0;
			for (auto k = j->controlPoints.begin(); k!=j->controlPoints.end(); k++) { //for all points in the bezier curves
				double coeff =(double) pow(1-seed,n-count)*(double) pow(seed,count)*(double) boost::math::binomial_coefficient<double>(n,count);
				if (t==(max_nodes) && k==(--j->controlPoints.end())) {
					coeff = 1;
				}
				count++;
				x+=coeff*CGAL::to_double(k->x());
				y+=coeff*CGAL::to_double(k->y());
			}


			m.add_vertex(Point(x,y));

			}
		}
	}

//	auto st = m.vertices_begin();
//	auto fi = m.vertices_end();

//	m = getTestMesh();

	// begin meshing algorithm

	std::list<std::list<Mesh::Vertex_index>> sm_boundaries; //list of boundaries!

	std::list<Mesh::Vertex_index> sm_boundary; //populate boundary list with vertex index

	std::list<Mesh::Vertex_index> new_sm_boundary; //used for managing boundary updates

	for (auto i = m.vertices_begin(); i!=m.vertices_end(); i++) {
		sm_boundary.push_back(*i); //getting correct vertex handle!
	}

	sm_boundaries.push_back(sm_boundary);

while (!sm_boundaries.empty()) { //as long a14s you still have elements to fill in...
//TODO: mark
	//for (int s=1; s<4; s++) {

	//if (s==5) {
	//	std::cout<< "left: " << sm_boundaries.size()<< " size: " << sm_boundaries.begin()->size()<<"\n";
//	}


	std::list<Mesh::Vertex_index> new_sm_boundary; //used for managing boundary updates


		if (sm_boundaries.begin()->size()==6 ){ //down to 6 nodes - mesh using predefined templates
			applySixNodeSplitter(*sm_boundaries.begin(),m);
			sm_boundaries.pop_front();
			continue;
		}

		if (sm_boundaries.begin()->size()==4) { //down to 4 nodes - make element

			auto j = sm_boundaries.begin()->end();
			std::cout << "Creating 4-element Node!\n";
			m.add_face (*std::next(j,-1),*std::next(j,-2),*std::next(j,-3),*std::next(j,-4));
			sm_boundaries.pop_front();
			continue;
		}


		std::cout << "Check if elements need to be split\n";

		bool nullFlag = false; //Are there any legal candidate nodes?

		bool intersectFlag = false; //Does a segment from current to node to candidate node intersect the boundary?

		Kernel::FT D1,D2,D3,D4; //Distances according to paper

		std::list<Mesh::Vertex_index>::iterator v1,v2;

		Kernel::FT minDistance;

		// i is pointer to current node
		for (auto i=sm_boundaries.begin()->begin(); i!=sm_boundaries.begin()->end(); i++) {

			auto next_i = getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),i,1);
			auto last_i = getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),i,-1);

			//normal vector to current node with respect to boundary
			auto vec_i = getOffsetPoint(m.point(*last_i),m.point(*i),m.point(*next_i)) - m.point(*i);
			//normalize vector
			auto norm_vec_i = vec_i/std::sqrt(CGAL::to_double(vec_i.squared_length()));

				//go over all other nodes in the boundary and see if valid candidates. j is pointer to candidate node
				//TODO: can be optimized according using scan to find visible nodes (in talberts paper) - will get complexity down from O(n^3) to O(n^2)
				for (auto j=getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),i,3); j!=getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),i,-3); j=getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),j,1)) {

				//	std::cout<<"PROBLEM! "<< m.point(*j)<<"\n";


					auto next_j = getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),j,1);
					auto last_j = getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),j,-1);



					Segment seg;

					//normal vector for candidate vertex
					auto vec_j = getOffsetPoint(m.point(*last_j),m.point(*j),m.point(*next_j)) - m.point(*j);

					//normalize vector
					auto norm_vec_j = vec_j/std::sqrt(CGAL::to_double(vec_j.squared_length()));

					//make sure elements are on collision path - Criterion: angle between normal vectors between 90 and 180 degrees (not given in paper)
					if ((norm_vec_j*norm_vec_i>Min_Collision_Value)) {
						continue;
					}

							seg = Segment(m.point(*i),m.point(*j));

							//check if segment intersects intersects boundary
							for (auto k=sm_boundaries.begin()->begin(); k!=sm_boundaries.begin()->end(); k++) { //check if no intersections with other segments in boundary (boundary may be concave)

								intersectFlag = false;

								auto next_k = getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),k,1);

								// if segment does includes the candidate nodes - continue to next segment
								if ((k==i) || (next_k==i) || (k==j) || (next_k==j)) {
									continue;
								}

								Segment tempSeg = Segment(m.point(*k),m.point(*next_k));

								// if there exists intersection - go to next candidate loop
								if (CGAL::do_intersect(seg,tempSeg)) {
									intersectFlag=true;
									break;
								}

							}

					// if previous loop was exited with no intersection - i and j are a legal candidate pair
					if (intersectFlag==false) {
						//calculate distance - update node pair if minimal
						if (nullFlag==false) { // if this is the first legal pair of vertices
							//set node pair

					//		std::cout<<"init minimal distance\n";
							v1 = i;
							v2 = j;
							//set current minimal distance
							minDistance= std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*j)));
							nullFlag=true;
							//calculate D1,D2,D3,D4 according to paper
							D1 = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*last_i)));
							D2 = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*next_i)));
							D3 = std::sqrt(CGAL::squared_distance(m.point(*j),m.point(*last_j)));
							D4 = std::sqrt(CGAL::squared_distance(m.point(*j),m.point(*next_j)));
						} else {

					//		std::cout<<"update minimal distance\n";
							//if needed update minimal distance, vertex pair, D1,D2,D3,D4
							if (minDistance>std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*j)))) {
								minDistance = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*j)));
								v1=i;
								v2=j;
								D1 = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*last_i)));
								D2 = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*next_i)));
								D3 = std::sqrt(CGAL::squared_distance(m.point(*j),m.point(*last_j)));
								D4 = std::sqrt(CGAL::squared_distance(m.point(*j),m.point(*next_j)));

							}

						}

					}
				}
		}

		if ((nullFlag!=false) && ((Relaxation_Factor*((D1+D2)/2+(D3+D4)/2))>=minDistance)) {

		//check if minDistance acquired, justifies splitting the boundary, as described in paper

			std::cout<<"splitting boundary\n";
			std::list<Mesh::Vertex_index> new_sm_boundary_1, new_sm_boundary_2;
			//start from i, move counterclockwise until encounter with j
			auto k=v1;
			auto midSplitNode = m.point(*v1)+Kernel::Vector_2(m.point(*v2)-m.point(*v1))/2;
			int counter=0;
			auto midSplitNode_index = *v1;  //initialize variable with auto
			new_sm_boundary_1.push_back(*k);
			counter++;
			while (k!=v2) {
			//	std::cout<<"in loop 1: "<<m.point(*v1)<<" "<< m.point(*v2)<< "\n";
				k++;
				counter++;
				if (k==sm_boundaries.begin()->end()) {
					k=sm_boundaries.begin()->begin();
				}
				new_sm_boundary_1.push_back(*k);

			} //ends with k=v2
			if (counter % 2 ==1) {
				midSplitNode_index = m.add_vertex(midSplitNode);
				new_sm_boundary_1.push_back(midSplitNode_index);
			}
			new_sm_boundary_2.push_back(*k);
			while (k!=v1) {
				k++;
				if (k==sm_boundaries.begin()->end()) {
					k=sm_boundaries.begin()->begin();
			}
			new_sm_boundary_2.push_back(*k);
			} //ends with k=v1

			if (counter %2 ==1) {
				new_sm_boundary_2.push_back(midSplitNode_index);
			}

			sm_boundary = new_sm_boundary_1;
			sm_boundaries.push_back(sm_boundary);
			sm_boundary = new_sm_boundary_2;
			sm_boundaries.push_back(sm_boundary);
			sm_boundaries.pop_front();

			continue;

		}


	// if boundary does not need to be split - offset elements

	std::cout<<"offsetting element row\n";

	//Offseting elements - handle first node

	auto last_v = m.point(*std::next(sm_boundaries.begin()->end(),-1));

	auto  curr_v = m.point(*(sm_boundaries.begin()->begin()));

	auto next_v =m.point( *(std::next(sm_boundaries.begin()->begin(),1)));

	auto new_v = getOffsetPoint(last_v,curr_v,next_v);

	auto new_index = m.add_vertex(new_v);

	new_sm_boundary.push_back(new_index);

	auto old_index = new_index; //initialized with auto

	auto mid_node1_index = new_index; //initialized with auto

	auto firstOffsetNodeIndex= new_index; //initialized with auto

	bool midVertexFlag=false;
	bool nodeEliminatedFlag=false;

	//loop over other nodes

	for (auto i=(++(sm_boundaries.begin()->begin())); i!=sm_boundaries.begin()->end() ;i++) {

		last_v = m.point(*getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),i,-1));

		curr_v = m.point(*i);

		next_v = m.point(*getNext(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end(),i,1));

		if ((last_v==curr_v) || (curr_v==next_v)) {
			std::cout<<"Error here";
		}
		new_v = getOffsetPoint(last_v,curr_v,next_v);

		auto Di = CGAL::squared_distance(new_v,m.point(*i));

		auto Dij = CGAL::squared_distance(new_v,m.point(old_index));

		if (nodeEliminatedFlag) {
			m.add_face(old_index,*i,*std::next(i,-1),*std::next(i,-2));
			nodeEliminatedFlag=false;
			continue;
		}
		if ((Dij>Kernel::FT(1.453*1.453)*Di) && (1!=1)){ //check if middle node needs to be added
			std::cout << "need to add middle node!\n";
			if (!midVertexFlag) {
				auto mid_node1 = m.point(old_index)+(new_v-m.point(old_index))/2;
				mid_node1_index = m.add_vertex(mid_node1);
				new_sm_boundary.push_back(mid_node1_index);
				//add new face
				m.add_face(old_index,mid_node1_index,*i,*(std::next(i,-1)));
				midVertexFlag=true;
				old_index =  mid_node1_index;
			} else { //midVertexFlag=true
				auto mid_node2 =  m.point(old_index)+(new_v-m.point(old_index))/2;
				auto  mid_node2_index = m.add_vertex(mid_node2);
				new_sm_boundary.push_back(mid_node2_index);
				m.add_face(mid_node2_index,*i,*(std::next(i,-1)),mid_node1_index);
				old_index = mid_node2_index;
				midVertexFlag=false;
			}
		}else if ((Dij<Kernel::FT(0.517*0.517)* Di) && i!=std::next(sm_boundaries.begin()->end(),-1)) { //check if nodes need to be eliminated
			nodeEliminatedFlag =true;
			midVertexFlag=false;
		} else { //node is fine*/

		new_index = m.add_vertex(new_v);
		new_sm_boundary.push_back(new_index);
		m.add_face(*i,*std::next(i,-1),old_index,new_index);
		old_index = new_index;
		midVertexFlag=false;

		}

	}

	//stiching to close offset element row

		if (midVertexFlag==true) {
			std::cout <<"exited with mid vertex flag \n";
		}
      m.add_face(*(sm_boundaries.begin()->begin()),*(std::next(sm_boundaries.begin()->end(),-1)),old_index,firstOffsetNodeIndex);

      std::cout<<m.point(*(sm_boundaries.begin()->begin())) << " " <<m.point(*(std::next(sm_boundaries.begin()->end(),-1)))<< " "<<m.point(old_index) <<" " <<m.point(firstOffsetNodeIndex);


      //perform local smoothing

	std::cout<<"performing local smoothing\n";
	for (auto i=new_sm_boundary.begin(); i!=new_sm_boundary.end(); i++) {

		if (m.degree(*i)==3)
		{
	    CGAL::Vertex_around_target_circulator<Mesh> vbegin(m.halfedge(*i),m), done(vbegin);

	    auto vbegin1=vbegin;

	    auto next_vertex = i;
	    auto last_vertex = i;
	    auto prev_v=*i;
	    if (i==new_sm_boundary.begin()) {
	    	 last_vertex = std::next(new_sm_boundary.end(),-1);
	    } else {
	    	 last_vertex = std::next(i,-1);
	    }

	    if (i==--new_sm_boundary.end()) {
	    	 next_vertex = new_sm_boundary.begin();
	    } else {
	    	 next_vertex = std::next(i,1);
	    }

	    do {
	      if ((*vbegin!=*last_vertex) && (*vbegin!=*next_vertex) && (*vbegin!=*i)) {
	    	  prev_v = *vbegin;
	     	  break;
	      }
	      vbegin++;
	    } while(vbegin != done);

	    auto V1 = m.point(*last_vertex)-m.point(prev_v);
	    auto V2 = m.point(*next_vertex)-m.point(prev_v);

	    auto V3 = (V1+V2)/2;

	   m.point(*i) = m.point(prev_v) + V3;


		}



}
	sm_boundary = new_sm_boundary;
	sm_boundaries.push_back(sm_boundary);

	sm_boundaries.pop_front();

}


//perform global smoothing on mesh


for (auto t=1; t<=6 ; t++){

std::cout<<"performing global smoothing\n";

for (auto i=m.vertices_begin(); i!=m.vertices_end(); i++) {
	if (m.is_border(*i,true)) {continue;}
	CGAL::Vertex_around_target_circulator<Mesh> vbegin(m.halfedge(*i),m), done(vbegin);
    Point newPoint = CGAL::ORIGIN;
    int N=0;
    do{
    	newPoint += m.point(*vbegin)-CGAL::ORIGIN;
    	vbegin++;N++;

    } while (vbegin!=done);

   m.point(*i) = Point(newPoint.x()/N,newPoint.y()/N);
}
}

/*

//perform improvements to unsatisfactory elements
//go over all elements

bool badElementFlag = false;
auto f_it = m.faces_begin();
	while (f_it!=m.faces_end()) {
		//Case 1: check if element shares two of its nodes with neighbor
		auto e = m.halfedge(*f_it);
		auto s1 = e;
		auto s2 = m.next(s1);
		bool Case1=false;
		bool firstFlag = true;
		while (s1!=e || firstFlag==true ) {
			firstFlag=false;
		//	std::cout << (m.face(m.opposite(s2)) == m.face(m.opposite(s1))) << ", ";
			if (m.face(m.opposite(s2)) == m.face(m.opposite(s1)) && !m.is_border(m.opposite(s2))) {
				CGAL::Euler::remove_center_vertex(s1,m);
				Case1=true;
	//			m.remove_vertex(m.vertex(m.edge(s1),0));
	//			m.remove_edge(m.edge(s1));
	//			m.remove_edge(m.edge(s2));
	//			std::cout <<"Vertex to delete: " <<m.point(m.vertex(m.edge(s1),0)) << "\n";
				break;

			}
			auto temp = m.next(s1);
			s1 = temp;
			s2 = m.next(s1);
		}

		if (Case1) {
			m.collect_garbage();
			f_it=m.faces_begin();
			continue;
		}*/
/*
		//Case 2: check antipodal pair of two nodes on element.

		bool Case2=false;

		auto vbegin = ((m.vertices_around_face(m.halfedge(*f_it)))).begin();

		 v1=*vbegin;
				vbegin++;
		 v2 = *vbegin;
				vbegin++;
		 v3 = *vbegin;
				vbegin++;
		 v4 = *vbegin;

		if (m.degree(v1)==3 && m.degree(v3)==3 && !m.is_border(v1,true) && !m.is_border(v3,true)) {

			std::cout<<"Delete case 2: Points to delete "<< m.point(v2) << " , " << m.point(v4) <<"\n";
			auto f1 = m.face(m.opposite(m.halfedge(v1,v2)));
			auto f2 = m.face(m.opposite(m.halfedge(v2,v3)));
			auto f3 = m.face(m.opposite(m.halfedge(v3,v4)));
			auto f4 = m.face(m.opposite(m.halfedge(v4,v1)));

			m.remove_face(f1); m.remove_face(f2); m.remove_face(f3); m.remove_face(f4);
			m.remove_face(*f_it);

			auto newVertex = m.point(v1) + (m.point(v3)-m.point(v1))/2;
			m.remove_edge(m.edge(m.halfedge(v1,v2)));
			m.remove_edge(m.edge(m.halfedge(v2,v3)));
			std::cout<<newVertex<<" , "<<m.point(v1) <<"\n";
			m.point(v2)=newVertex;
			m.set_target(m.halfedge(v3),v1);
			std::cout<<newVertex<<" , "<<m.point(v1) <<"\n";
			m.remove_vertex(v3);

			Case2=true;

		}

		if (m.degree(v2)==3 && m.degree(v4)==3 && !m.is_border(v2,true) && !m.is_border(v4,true)) {

			std::cout<<"Delete case 2: Points to delete "<< m.point(v2) << " , " << m.point(v4) <<"\n";
			auto f1 = m.face(m.opposite(m.halfedge(v1,v2)));
			auto f2 = m.face(m.opposite(m.halfedge(v2,v3)));
			auto f3 = m.face(m.opposite(m.halfedge(v3,v4)));
			auto f4 = m.face(m.opposite(m.halfedge(v4,v1)));

			m.remove_face(f1); m.remove_face(f2); m.remove_face(f3); m.remove_face(f4);
			m.remove_face(*f_it);

			auto newVertex = m.point(v2) + (m.point(v4)-m.point(v2))/2;
			m.remove_edge(m.edge(m.halfedge(v3,v4)));
			m.remove_edge(m.edge(m.halfedge(v4,v1)));
			std::cout<<newVertex<<" , "<<m.point(v2) <<"\n";
			m.point(v2)=newVertex;
			m.set_target(m.halfedge(v4),v2);
			std::cout<<newVertex<<" , "<<m.point(v2) <<"\n";
			m.remove_vertex(v4);

			Case2=true;
		}

		if (Case2) {

		}


		//Case 3 & 4: for every vertex check if degree is 3 *and* it is not a border vertex

		auto start = ((m.vertices_around_face(m.halfedge(*f_it)))).begin();

		 e = m.halfedge(*f_it); //just to initialize

		  v1=*start;
		  v2=*(++start);
		  v3=*(++start);
		  v4=*(++start);
		auto v_selected = m.null_vertex();
		bool Case3=false;

		std::cout << "checking borders: " << m.is_border(v1,true) <<" "<<m.is_border(v2,true)<<" "<<m.is_border(v3,true)<<" "<<m.is_border(v4,true)<<"\n";

		if (m.degree(v1)==3 && !m.is_border(v1,true)) {
			std::cout<<"selected vertex is v1: "<<m.point(v1)<<"\n";
			v_selected=v1;
			Case3=true;
		} else if (m.degree(v2)==3  && !m.is_border(v2,true)) {
			std::cout<<"selected vertex is v2: "<<m.point(v2)<<"\n";
			v_selected=v2;
			Case3=true;
		} else if (m.degree(v3)==3  && !m.is_border(v3,true)) {
			std::cout<<"selected vertex is v3: "<<m.point(v3)<<"\n";
			v_selected=v3;
			Case3=true;
		} else if (m.degree(v4)==3  && !m.is_border(v4,true)) {
			std::cout<<"selected vertex is v4: "<<m.point(v4)<<"\n";
			v_selected=v4;
			Case3=true;
		}

		if (Case3) {
			e = CGAL::Euler::remove_center_vertex(m.halfedge(v_selected),m);
			auto it = m.vertices_around_face(e);
			Boundary tempBoundary;
			tempBoundary.assign(it.begin(),it.end());

			//test nodes in boundary
			for (auto k = tempBoundary.begin(); k!=tempBoundary.end(); k++) {
				std::cout<<m.point(*k) <<" , ";
			}
			std::cout<<"\n";

			CGAL::Euler::make_hole(e,m);

			applySixNodeSplitter (tempBoundary,m);

		//	m.add_face(*tempBoundary.begin(),*getNext(tempBoundary.begin(),tempBoundary.end(),(tempBoundary.begin()),1),*getNext(tempBoundary.begin(),tempBoundary.end(),(tempBoundary.begin()),2),*getNext(tempBoundary.begin(),tempBoundary.end(),(tempBoundary.begin()),3));
		//	m.add_face(*getNext(tempBoundary.begin(),tempBoundary.end(),(tempBoundary.begin()),3),*getNext(tempBoundary.begin(),tempBoundary.end(),(tempBoundary.begin()),4),*getNext(tempBoundary.begin(),tempBoundary.end(),(tempBoundary.begin()),5),*tempBoundary.begin());

			if (m.has_garbage()) {
			m.collect_garbage();
			std::cout <<"Number of faces: "<<m.num_faces()<< " Number of edges: " <<m.num_edges()<<"\n";
			//check which face remains
			for (auto i=m.faces_begin(); i!=m.faces_end(); i++) {
				for (auto k=m.vertices_around_face(m.halfedge(*i)).begin(); k!=m.vertices_around_face(m.halfedge(*i)).end(); k++) {
					std::cout<<m.point(*k) << " , ";
				}
				std::cout<<"\n";
			}
			}


			badElementFlag=true;
			f_it = m.faces_begin();
			continue;

		}

*/

//			f_it++;






/*
		//case 4: for every edge of face, check if source and target both have degree 3


		auto ebegin = m.halfedges_around_face(m.halfedge(*i));
		for (auto i=ebegin.begin(); i!=ebegin.end(); i++) {
			if (m.degree((m.source(*i)))==3 && m.degree(m.target(*i))==3) {
				//go to previous edge
				//go to opposite edge
				//add target vertex -1
				//go to next edge
				//add target vertex -2
				//go to next edge
				//add target vertex -3
				//go to next edge
				//go to opposite edge
				//go to next edge
				//add target vertex - 4
				//go to next edge
				//got to opposite edge
				//go to next edge
				//add target vertex -5
				//got to next edge
				//add target vertex - 6

				break;i
			}


		}

*/

//	}




//std::cout <<"Number of faces: "<<m.num_faces()<< " Number of edges: " <<m.num_edges()<<"\n";




for (auto t=1; t<=6; t++){

std::cout<<"performing global smoothing\n";

for (auto i=m.vertices_begin(); i!=m.vertices_end(); i++) {
	if (m.is_border(*i,true)) {continue;}
	CGAL::Vertex_around_target_circulator<Mesh> vbegin(m.halfedge(*i),m), done(vbegin);
    Point newPoint = CGAL::ORIGIN;
    int N=0;
    do{
    	newPoint += m.point(*vbegin)-CGAL::ORIGIN;
    	vbegin++;N++;

    } while (vbegin!=done);

   m.point(*i) = Point(newPoint.x()/N,newPoint.y()/N);
}
}

std::cout<<"got here \n";

//testing node placement using ipe
std::ofstream myFile;
myFile.open("Ipe.xml");
myFile << "<page>\n";
for (auto i = m.vertices().begin(); i!=m.vertices().end(); i++) {
//	std::cout <<m.point(*i);
//	std::cout<<"\n";


myFile << "<use name=\"mark/disk(sx)\" " << "pos= \"" << m.point(*i).x() << " " << m.point(*i).y() << "\" size=\"normal\" stroke=\"black\"/>\n";



}

std::cout<<"got here \n";

for (auto i = m.edges_begin(); i!=m.edges_end(); i++) {
//	std::cout <<m.point(*i);
//	std::cout<<"\n";

Point p1 = m.point((m.vertex(*i,0)));

Point p2 = m.point((m.vertex(*i,1)));

myFile << "<path stroke = \"black\"> \n"  << p1 <<" m \n" << p2 << " l \n" << "</path> \n";


}
std::cout<<"got here \n";
myFile << "</page>\n";
myFile.close();
}








