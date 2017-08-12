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
#define 	CGAL_PI   3.14159265358979323846


Point getOffsetPoint (Point last_v, Point curr_v, Point next_v) {
	Kernel::Vector_2 vec1 = Kernel::Vector_2(curr_v,last_v);

	Kernel::Vector_2 vec2 = Kernel::Vector_2(curr_v,next_v);

	auto vec1_norm = vec1/Kernel::FT(std::sqrt(CGAL::to_double(vec1.squared_length())));

	auto vec2_norm = vec2/Kernel::FT(std::sqrt(CGAL::to_double(vec2.squared_length())));

	Kernel::FT D = Kernel::FT((std::sqrt(CGAL::to_double(vec1.squared_length()))+std::sqrt(CGAL::to_double(vec2.squared_length())))/2);
	Kernel::Vector_2 Theta;

	 auto Angle = acos(vec1_norm*vec2_norm);
	// std::cout << Angle << "\n";
	Transformation rotate2(CGAL::ROTATION, std::sin(Angle/2), std::cos(Angle/2));

	Theta = rotate2(vec2_norm);
	auto new_v = curr_v + Theta*D;

	return new_v;
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

int max_nodes=5; //how many nodes per edge

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

	auto st = m.vertices_begin();
	auto fi = m.vertices_end();

	// begin meshing algorithm

	std::list<std::list<Mesh::Vertex_index>> sm_boundaries; //list of boundaries!

	std::list<Mesh::Vertex_index> sm_boundary; //populate boundary list with vertex index

	std::list<Mesh::Vertex_index> new_sm_boundary; //used for managing boundary updates

	for (auto i = m.vertices_begin(); i!=m.vertices_end(); i++) {
		sm_boundary.push_back(*i); //getting correct vertex handle!
	}

	sm_boundaries.push_back(sm_boundary);

while (!sm_boundaries.empty()) { //as long as you still have elements to fill in...
//for (int s=1; s<100; s++) {

	std::cout<< "Size: " << sm_boundaries.begin()->size() <<"\n" ;

	//std::cout << s << " \n";

	for (auto k=sm_boundaries.begin()->begin(); k!=sm_boundaries.begin()->end(); k++) {
	//std::cout<<m.point(*k)<<" \n ";
	}

	std::list<Mesh::Vertex_index> new_sm_boundary; //used for managing boundary updates


		if (sm_boundaries.begin()->size()==6 ){ //mesh using predefined templates
			if (sm_boundaries.begin()->size()!=4 && sm_boundaries.begin()->size()!=6) {
		//		std::cout << "size=5";
			}

			//implement six node splitters
			// need to keep track of the number of 180 degree angles, and number of consecutive(!) 180 degree angles
			int numOfCollinearVertex=0;
			int temp=0;
			int maxNumOfConsecutiveCollinearVertex=0;
			int numOfNormalVertices = 0;
			auto splitterVertex = sm_boundaries.begin()->begin();
			bool collinearFlag=false;

			for (auto j=sm_boundaries.begin()->begin(); j!=sm_boundaries.begin()->end(); j++) {
				auto next_vertex = std::next(j,1);
				auto last_vertex = std::next(j,-1);
				if (CGAL::collinear(m.point(*last_vertex),m.point(*j),m.point(*next_vertex))) {

					numOfCollinearVertex++;
					temp++;
					splitterVertex = *j;
					numOfNormalVertices=0;
					if (j==--sm_boundaries.begin()->end())  { //If starting at middle of collinear vertex sequence
						temp ++;
						splitterVertex = sm_boundaries.begin()->begin();
					}

				} else {
					if (temp>maxNumOfConsecutiveCollinearVertex) {
						maxNumOfConsecutiveCollinearVertex=temp;
						numOfNormalVertices++;
					}
					temp=0;
				}
			}/*
			if (sm_boundaries.begin()->size()==6) {
				auto j=sm_boundaries.begin()->begin();
				//m.add_face (*j,*(std::next(j,1)),*(std::next(j,2)),*(std::next(j,3)));
				m.add_face (*(std::next(j,3)),*(std::next(j,2)),*(std::next(j,1)),*j);
		//	m.add_face (*(std::next(j,3)),*(std::next(j,4)),*(std::next(j,5)),*j);
				m.add_face (*j,*(std::next(j,5)),*(std::next(j,4)),*(std::next(j,3)));
				std::cout << "here : " << m.point(*j) << " "<< m.point(*(std::next(j,1))) << " " << m.point(*std::next(j,2))  << " " <<  m.point(*std::next(j,3))  << " "<<  m.point(*std::next(j,4)) << " " <<  m.point(*std::next(j,5)) << "\n";
			} */

			auto j=sm_boundaries.begin()->begin();

			if (numOfCollinearVertex==0) {
				m.add_face (*(std::next(splitterVertex,3)),*(std::next(splitterVertex,2)),*(std::next(splitterVertex,1)),*splitterVertex);
				m.add_face (*splitterVertex,*(std::next(splitterVertex,5)),*(std::next(splitterVertex,4)),*(std::next(splitterVertex,3)));
			}
			if (numOfCollinearVertex==1) {
				m.add_face (*getNext(splitterVertex,3),*getNext(splitterVertex,2),*getNext(splitterVertex,1),*splitterVertex);
				m.add_face (*splitterVertex,*getNext(splitterVertex,5),*getNext(splitterVertex,4),*getNext(splitterVertex,3));
			}
			if (numOfCollinearVertex==2 && maxNumOfConsecutiveCollinearVertex==1) {
				if (numOfNormalVertex==1) {
					m.add_face (*(std::next(splitterVertex,3)),*(std::next(splitterVertex,2)),*(std::next(splitterVertex,1)),*splitterVertex);
					m.add_face (*splitterVertex,*(std::next(splitterVertex,5)),*(std::next(splitterVertex,4)),*(std::next(splitterVertex,3)));
				} else if (numOfnormalVertex==3) {

				} else {
					m.add_face (*(std::next(splitterVertex,3)),*(std::next(splitterVertex,2)),*(std::next(splitterVertex,1)),*splitterVertex);
					m.add_face (*splitterVertex,*(std::next(splitterVertex,5)),*(std::next(splitterVertex,4)),*(std::next(splitterVertex,3)));
				}
			}
			if (numOfCollinearVertex==2 && maxNumOfConsecutiveCollinearVertex==2) {

			}
			if (numOfCollinearVertex==3 && maxNumOfConsecutiveCollinearVertex==3) {

			}
			if (numOfCollinearVertex==3 && maxNumOfConsecutiveCollinearVertex==2) {

			}
			if (numOfCollinearVertex==3 && maxNumOfConsecutiveCollinearVertex==1) {

			}

			sm_boundaries.pop_front();
			continue;
		}

		//TODO: for each node on the new boundary, compare with every other node on the boundary to find split lines


	//TODO:: add condition that valid pairs are only those with a [0,-1] scalar product between their normals.
		std::cout << "Check if elements need to be split\n";

		bool nullFlag = false;

		bool intersectFlag;

		Kernel::FT D1,D2,D3,D4;

		std::list<Mesh::Vertex_index>::iterator v1,v2;

		Kernel::FT minDistance;

		for (auto i=sm_boundaries.begin()->begin(); i!=sm_boundaries.begin()->end(); i++) {

			auto next_i=i;
			auto last_i=i;

			if (i==--sm_boundaries.begin()->end()) {
				next_i = sm_boundaries.begin()->begin();
			} else {
				 next_i = std::next(i,1);
			}

			if (i==sm_boundaries.begin()->begin()) {
				last_i = std::next(sm_boundaries.begin()->end(),-1);
			} else {
			    last_i = std::next(i,-1);
			}

			auto vec_i = getOffsetPoint(m.point(*last_i),m.point(*i),m.point(*next_i)) - m.point(*i);

			auto norm_vec_i = vec_i/std::sqrt(CGAL::to_double(vec_i.squared_length()));

				for (auto j=sm_boundaries.begin()->begin(); j!=sm_boundaries.begin()->end(); j++) {

					auto next_j=j;
					auto last_j=j;
						if (j==--sm_boundaries.begin()->end()) {
							next_j = sm_boundaries.begin()->begin();
						} else {
							next_j = std::next(j,1);
						}

						if (j==sm_boundaries.begin()->begin()) {
							last_j = std::next(sm_boundaries.begin()->end(),-1);
						} else {
							last_j = std::next(j,-1);
						}

					Segment seg;


					auto vec_j = getOffsetPoint(m.point(*last_j),m.point(*j),m.point(*next_j)) - m.point(*j);

					auto norm_vec_j = vec_j/std::sqrt(CGAL::to_double(vec_j.squared_length()));

					if ((norm_vec_j*norm_vec_i>0)) {
						break;
					}

						if ((*j!=*i) && (*next_j!=*i) && (*last_j!=*i)) {

							seg = Segment(m.point(*i),m.point(*j));


							for (auto k=sm_boundaries.begin()->begin(); k!=sm_boundaries.begin()->end(); k++) { //check if no intersections with other segments in boundary (boundary may be concave)

								 intersectFlag = false;

								 auto next_k=k;

								if (k==--sm_boundaries.begin()->end()) {
								next_k = sm_boundaries.begin()->begin();
								} else {
								 next_k = std::next(k,1);
								}

								if ((k!=i) && (next_k!=i) && (k!=j) && (next_k!=j)) {

									Segment tempSeg = Segment(m.point(*k),m.point(*next_k));


									if (CGAL::do_intersect(seg,tempSeg)) {
										intersectFlag=true;
										break;
									}
								}
							}
					if (intersectFlag==false) {
							//calculate distance - update node pair if minimal
						if (nullFlag==false) {
							v1 = i;
							v2 = j;
							minDistance= std::sqrt(CGAL::squared_distance(m.point(*v1),m.point(*v2)));
							nullFlag=true;
							D1 = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*last_i)));
							D2 = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*next_i)));
							D3 = std::sqrt(CGAL::squared_distance(m.point(*j),m.point(*last_j)));
							D4 = std::sqrt(CGAL::squared_distance(m.point(*j),m.point(*next_j)));
						} else {
							if (minDistance>std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*j)))) {
								minDistance = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*j)));
								v1=i;
								D1 = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*last_i)));
								D2 = std::sqrt(CGAL::squared_distance(m.point(*i),m.point(*next_i)));
								D3 = std::sqrt(CGAL::squared_distance(m.point(*j),m.point(*last_j)));
								D4 = std::sqrt(CGAL::squared_distance(m.point(*j),m.point(*next_j)));
								v2=j;
							}

						}

					}

				}

				}
		}

		//std::cout << nullFlag << " "  << (D1+D2)/2+(D3+D4)/2 << " " << minDistance << "\n";
	//	std::cout << m.point(*v1) << " "  << m.point(*v2)  << "\n";
		if ((nullFlag!=false) && ((1.3*(D1+D2)/2+(D3+D4)/2)>=minDistance)) {

		//check if minDistance acquired, justifies splitting the boundary, as described in paper

			std::cout<<"splitting boundary\n";

		std::list<Mesh::Vertex_index> new_sm_boundary_1, new_sm_boundary_2;
			//start from i, move counterclockwise until you see j

		auto k=v1;
		auto midSplitNode = m.point(*v1)+Kernel::Vector_2(m.point(*v2)-m.point(*v1))/2;
		auto midSplitNode_index = m.add_vertex(midSplitNode);

		new_sm_boundary_1.push_back(*k);
		while (k!=v2) {
			k++;
			if (k==sm_boundaries.begin()->end()) {
					k=sm_boundaries.begin()->begin();
			}
			new_sm_boundary_1.push_back(*k);

		} //ends with k=v2
		new_sm_boundary_1.push_back(midSplitNode_index);

		new_sm_boundary_2.push_back(*k);
		while (k!=v1) {
			k++;
			if (k==sm_boundaries.begin()->end()) {
					k=sm_boundaries.begin()->begin();
			}
			new_sm_boundary_2.push_back(*k);
		} //ends with k=v1
		new_sm_boundary_2.push_back(midSplitNode_index);

		sm_boundary = new_sm_boundary_1;
		sm_boundaries.push_back(sm_boundary);
		sm_boundary = new_sm_boundary_2;
		sm_boundaries.push_back(sm_boundary);
		sm_boundaries.pop_front();
		continue;

		}




	std::cout<<"offsetting element row\n";

	//CGAL::Circulator_from_iterator<std::list<Mesh::Vertex_index>::iterator> circ(sm_boundaries.begin()->begin(),sm_boundaries.begin()->end());

	//first node

	auto last_v = m.point(*std::next(sm_boundaries.begin()->end(),-1));

	auto  curr_v = m.point(*(sm_boundaries.begin()->begin()));

	auto next_v =m.point( *(std::next(sm_boundaries.begin()->begin(),1)));

	Kernel::Vector_2 vec1 = Kernel::Vector_2(curr_v,last_v);

	Kernel::Vector_2 vec2 = Kernel::Vector_2(curr_v,next_v);

	auto vec1_norm = vec1/Kernel::FT(std::sqrt(CGAL::to_double(vec1.squared_length())));

	auto vec2_norm = vec2/Kernel::FT(std::sqrt(CGAL::to_double(vec2.squared_length())));

	Kernel::FT D = Kernel::FT((std::sqrt(CGAL::to_double(vec1.squared_length()))+std::sqrt(CGAL::to_double(vec2.squared_length())))/2);
	Kernel::Vector_2 Theta;

	 auto Angle = acos(vec1_norm*vec2_norm);
	// std::cout << Angle << "\n";
	Transformation rotate2(CGAL::ROTATION, std::sin(Angle/2), std::cos(Angle/2));

	Theta = rotate2(vec2_norm);
	auto new_v = curr_v + Theta*D;

	auto new_index = m.add_vertex(new_v);

	new_sm_boundary.push_back(new_index);

	auto old_index = new_index;

	auto firstOffsetNodeIndex= new_index;

	//loop over other nodes

	for (auto i=(++(sm_boundaries.begin()->begin())); i!=std::next(sm_boundaries.begin()->end(),-2) ;i++) {


		last_v = m.point(*std::next(i,-1));

		curr_v = m.point(*i);

		auto j = std::next(i,1);

		if (i==--(sm_boundaries.begin()->end())) {

			j=sm_boundaries.begin()->begin();

		}

		next_v = m.point(*j);

		vec1 = Kernel::Vector_2(curr_v,last_v);

		vec2 = Kernel::Vector_2(curr_v,next_v);

		 vec1_norm = vec1/Kernel::FT(std::sqrt(CGAL::to_double(vec1.squared_length())));

		 vec2_norm = vec2/Kernel::FT(std::sqrt(CGAL::to_double(vec2.squared_length())));


		  D = Kernel::FT((std::sqrt(CGAL::to_double(vec1.squared_length()))+std::sqrt(CGAL::to_double(vec2.squared_length())))/2);

		   Angle = acos(vec1_norm*vec2_norm);

		  Transformation rotate2(CGAL::ROTATION, std::sin(Angle/2), std::cos(Angle/2));

		 Theta = rotate2(vec2_norm);


		new_v = curr_v + Theta*D;


		auto Dij = CGAL::squared_distance(new_v,m.point(old_index));

		auto Di = CGAL::squared_distance(new_v,m.point(*i));

		if ((Dij<Kernel::FT(0.517*0.517)* Di)) { //check if nodes need to be eliminated
			m.add_face(*std::next(i,-1),old_index,*j,*i);
			// go over next node
			if (i!=std::next((sm_boundaries.begin()->end()),-3)) {
			i++;
			}
	/*	} else if (Dij>Kernel::FT(1.453*1.453)* Di){ //check if middle node needs to be added
		//still needs some work..unclear about implementation (can create triangular element...)
			auto mid_node = m.point(old_index)+(m.point(old_index)-new_v)/2;
			auto  mid_node_index = m.add_vertex(mid_node);
			m.add_face(*std::next(i,-1),old_index,mid_node_index,*i);
			i--;
*/
		} else { //node is fine*/

			if (i==--sm_boundaries.begin()->end()) {
		//		std::cout << "Test1 \n" << last_v << " "<< curr_v << " " <<next_v << " \n";
		//		std::cout << "Test2 \n" << vec1_norm << " "<< vec2_norm << " " << " \n";
		//		std::cout << "Test3 \n" << std::acos(vec1_norm*vec2_norm)  << " \n";
		//		std::cout << "Test4 \n" << Angle  << " \n";
			}

		new_index = m.add_vertex(new_v);

		new_sm_boundary.push_back(new_index);


		m.add_face(*i,*std::next(i,-1),old_index,new_index);

		old_index = new_index;

		}

	}
		auto thirdToLastv = std::next(sm_boundaries.begin()->end(),-3);
	 auto secondToLastv = std::next(sm_boundaries.begin()->end(),-2);
	 auto Lastv = std::next(sm_boundaries.begin()->end(),-1);
	 auto Firstv = sm_boundaries.begin()->begin();
//	 std::cout << "Test4 " << m.point(*Lastv)  <<" " << m.point(*Firstv) <<" " <<m.point(*secondToLastv) <<" \n";
      m.add_face(*secondToLastv, *thirdToLastv, old_index, firstOffsetNodeIndex);
      m.add_face(*secondToLastv,firstOffsetNodeIndex,*Firstv,*Lastv);

	//TODO: perform local smoothing

	//for each new node in the new boundary...

	//assumption that boundary nodes are connected to exactly 3 other nodes (not always true..)
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
	   // 	  std::cout << m.point(*vbegin) << "\n";
	    	  break;
	      }
	      vbegin++;
	    } while(vbegin != done);

	    auto V1 = m.point(*last_vertex)-m.point(prev_v);
	    auto V2 = m.point(*next_vertex)-m.point(prev_v);

	    auto V3 = (V1+V2)/2;

	   m.point(*i) = m.point(prev_v) + V3;
	   if (m.point(*i)!=m.point(*i)) {
	  // std::cout << m.point(*last_vertex) << " "<<m.point(prev_v) << " " << m.point(*next_vertex) << "\n";
		   do {
		  	  //  	  std::cout << m.point(*vbegin) <<" ";

		  	      vbegin1++;
		  	    } while(vbegin1 != done);
		//   std::cout << "\n";
	   }
		}

	  }

	sm_boundary = new_sm_boundary;
	sm_boundaries.push_back(sm_boundary);

	sm_boundaries.pop_front();




}


//perform global smoothing on mesh

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



//perform improvements to unsatisfactory elements
//go over all elements
	for (auto i=m.faces_begin(); i!=m.faces_end(); i++) {
		//Case 1: check if element shares two of its nodes with neighbor
		//Case 2: check for two nodes on element which are on different neigbors!
	}




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








