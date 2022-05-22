#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <queue>
#include <vector>
#include <iterator>
#include <tuple>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/grid_graph.hpp>
#include <boost/array.hpp>
#include <cmath>
#include<fstream>
#include <time.h>


using namespace boost;




struct NodeInfo // Node Infos
{
	double distance;
	double h_t;
	int parent;
	int x;
	int y;

	NodeInfo(double d = 0, double h = 0, int p = 0, int x_cord = 0, int y_cord = 0)
		: distance(d), h_t(h), parent(p), x(x_cord), y(y_cord) {}

};

struct EdgeInfo // Edge Infos
{
	double weight;
	EdgeInfo(double w = 0) : weight(w) {}
};




//============= Boost Graph  =============//

typedef adjacency_list<vecS, vecS, bidirectionalS, NodeInfo, EdgeInfo> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::vertex_iterator vertex_it;
typedef graph_traits<Graph>::edge_iterator edge_it;
typedef graph_traits<Graph>::out_edge_iterator out;
typedef graph_traits<Graph>::adjacency_iterator adj_it;
typedef graph_traits<Graph>::in_edge_iterator in;
typedef graph_traits<Graph>::vertices_size_type VertexIndex;



typedef property_map<Graph, double EdgeInfo::*>::type WeightMap;
typedef property_map<Graph, double NodeInfo::*>::type DistanceMap;
typedef property_map<Graph, int NodeInfo::*>::type ParentMap;
typedef property_map<Graph, double NodeInfo::*>::type H_Map;



// Functor for priority queue 
struct CompareNodes {
	static property_map<Graph, double NodeInfo::*>::type v_dist;
	bool operator() (Vertex& v, Vertex& u) {
		return v_dist[v] > v_dist[u];
	}
};




property_map<Graph, double NodeInfo::*>::type CompareNodes::v_dist;



// Dijkstra declaration rev_dijkstra=true => performs reverse dijkstra
int Dijkstra_SP(Graph& G, Vertex s, Vertex t, WeightMap& w, DistanceMap& d, ParentMap& p, bool rev_dijkstra);


// For copying the contents of maps
void copy_map(Graph& G, DistanceMap& d, std::map<Vertex, int>& cp)
{
	vertex_it vi, vend;

	for (tie(vi, vend) = vertices(G); vi != vend; vi++)
	{
		cp[*vi] = d[*vi];
	}

}

// Heuristic calculation for random Graphs
void h_calc(Graph& G, WeightMap& w, H_Map& h, DistanceMap d, Vertex t)
{

	Vertex L1, L2;

	vertex_it vi, vi_end;

	edge_it ei;

	Vertex src, trg;

	mt19937 rng;
	rng.seed(uint32_t(time(0)));
	uniform_int<> u(0, num_vertices(G));
	variate_generator<mt19937&, uniform_int<> > rng1(rng, u);


	std::map<Vertex, int> L1_Map;

	std::map<Vertex, int> L2_Map;

	std::map<Vertex, int> L1_Rev;

	std::map<Vertex, int> L2_Rev;


	ParentMap p = get(&NodeInfo::parent, G);

	L1 = random_vertex(G, rng1);

	Dijkstra_SP(G, L1, -1, w, d, p, false);

	copy_map(G, d, L1_Map);

	Dijkstra_SP(G, L1, -1, w, d, p, true);

	copy_map(G, d, L1_Rev);

	L2 = random_vertex(G, rng1);


	while (L1 == L2)
	{
		L2 = random_vertex(G, rng1);
	}

	Dijkstra_SP(G, L2, -1, w, d, p, false);

	copy_map(G, d, L2_Map);

	Dijkstra_SP(G, L2, -1, w, d, p, true);

	copy_map(G, d, L2_Rev);


	for (vi = vertices(G).first; vi != vertices(G).second; vi++) {

		h[*vi] = std::max(L1_Map[t] - L1_Map[*vi], L2_Rev[*vi] - L2_Rev[t]);

	}


	// Alteration of weights for each edge based on Heuristic
	for (ei = edges(G).first; ei != edges(G).second; ++ei) {

		src = source(*ei, G);

		trg = target(*ei, G);

		w[*ei] = w[*ei] + h[trg] - h[src];

	}




}

// Heuristic calculation for 2D Graphs
void h_euclidean(Graph& G, H_Map& h, Vertex t, WeightMap& w)
{
	Vertex src, trg;

	vertex_it vi, vend;

	edge_it ei;

	double dist_x, dist_y;


	for (vi = vertices(G).first; vi != vertices(G).second; vi++) {

		dist_x = G[*vi].x - G[t].x;

		dist_y = G[*vi].y - G[t].y;

		h[*vi] = std::sqrt(std::pow(dist_x, 2) + std::pow(dist_y, 2));

		std::cout << "h[" << *vi << "]" << " =  " << h[*vi] << std::endl;

	}


	// Alteration of weights for each edge based on Heuristic
	for (ei = edges(G).first; ei != edges(G).second; ++ei) {

		src = source(*ei, G);

		trg = target(*ei, G);

		w[*ei] = w[*ei] + h[trg] - h[src];

	}



}


// Function for Graph 2D with coordinates
void Graph_2D(Graph& G, int r, int c, Vertex& s, Vertex& t)
{

	std::vector<std::vector<int>> vertices(r, std::vector<int>(c));
	Vertex tmp;


	// Create the vertices and add them to a 2D vector for random access
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			vertices[i][j] = add_vertex(G);
			tmp = vertices[i][j];

			G[tmp].x = i;
			G[tmp].y = j;

		}
	}

	// Add edges along the columns  v1-v2-v3-v4 ...
	//                              vk-vk+1-vk+2-vk+3 ...
	// initialize them with random values in [1,100]
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c - 1; j++)
		{
			Edge e = add_edge(vertices[i][j], vertices[i][j + 1], G).first;
			s = source(e, G);
			t = target(e, G);
			Edge e_rev = add_edge(t, s, G).first;
			G[e].weight = rand() % 2 + 1;
			G[e_rev].weight = rand() % 2 + 1;
		}

	}


	// Add edges along the rows  v1  v4
	//                           |   |
	//							 v2  v5
	//                           |   |
	//                           v3  v6 
	// initialize with random values in [1,100]
	for (int j = 0; j < c; j++)
	{
		for (int i = 0; i < r - 1; i++)
		{
			Edge e = add_edge(vertices[i][j], vertices[i + 1][j], G).first;
			s = source(e, G);
			t = target(e, G);
			Edge e_rev = add_edge(t, s, G).first;
			G[e].weight = rand() % 2 + 1;
			G[e_rev].weight = rand() % 2 + 1;
		}

	}




	s = vertices[rand() % (r - 1) + 0][0];

	t = vertices[rand() % (r - 1)][c - 1];


}

// A* Algorithm 
int A_Star(Graph& G, Vertex s, Vertex t, WeightMap& w, ParentMap& p, DistanceMap& d, H_Map& h)
{
	std::priority_queue<Vertex, std::vector<Vertex>, CompareNodes> q;

	int cnt_nodes = 0; // for counting the nodes till reaches target

	edge_it ei;

	vertex_it it, end;

	Vertex u, tmp, src, trg;

	Edge e;

	out out_i;


	d[s] = 0;
	h[t] = 0;


	if (s == t)
	{
		d[t] = 0;
		return 0;
	}



	q.push(s);

	while (!q.empty())
	{

		u = q.top();


		if (u == t)
		{

			d[t] += h[s];
			return cnt_nodes;
		}


		q.pop();

		for (out_i = out_edges(u, G).first; out_i != out_edges(u, G).second; out_i++) {

			cnt_nodes++;

			tmp = target(*out_i, G);

			if (d[tmp] > d[u] + w[*out_i])
			{
				d[tmp] = d[u] + w[*out_i];
				p[tmp] = u;
				q.push(tmp);
			}




		}

	}





}






// Dijkstra Algorithm
int Dijkstra_SP(Graph& G, Vertex s, Vertex t, WeightMap& w, DistanceMap& d, ParentMap& p, bool rev_dijkstra)
{
	std::priority_queue<Vertex, std::vector<Vertex>, CompareNodes> Q;
	vertex_it it, end;
	Vertex u, tmp;
	int alt = 0;
	Edge e;
	int cnt_nodes = 0;

	Q.push(s);

	out out_i;

	in in_i;


	if (s == t)
	{
		d[t] = 0;
		return 0;
	}

	while (!Q.empty())
	{
		u = Q.top();

		if (rev_dijkstra == false)
		{

			if (u == t)
			{
				return cnt_nodes;

			}
		}

		Q.pop();


		if (rev_dijkstra == false)
		{

			for (out_i = out_edges(u, G).first; out_i != out_edges(u, G).second; out_i++) {

				cnt_nodes++;

				tmp = target(*out_i, G);

				if (d[tmp] > d[u] + w[*out_i])
				{
					d[tmp] = d[u] + w[*out_i];
					p[tmp] = u;
					Q.push(tmp);
				}


			}
		}

		else
		{
			for (in_i = in_edges(u, G).first; in_i != in_edges(u, G).second; in_i++)
			{
				tmp = source(*in_i, G);

				if (d[tmp] > d[u] + w[*in_i])
				{
					d[tmp] = d[u] + w[*in_i];
					p[tmp] = u;
					Q.push(tmp);
				}
			}
		}

	}

}




int main()
{

	std::ofstream file; // File for writing Time Measures and Number of Nodes until it finalized the d(s,t) 

	int cnt = 0;

	int r, c;

	int a_dist, d_dist; // For checking if the two algorithms agree based on the distance

	double elapsed; // elapsed time for measurements the algorithms

	Graph G;

	Vertex s, t;

	vertex_it vi, vi_end;

	edge_it ei, end;

	WeightMap w_star, w_dijk;

	H_Map h;

	ParentMap p;

	mt19937 rng;
	rng.seed(uint32_t(time(0)));
	uniform_int<> u(0, 1000);
	variate_generator<mt19937&, uniform_int<> > rng1(rng, u);

	int vert, edg;

	/*
	vert = 1000;
	edg = 2000;

	generate_random_graph(G, vert, edg, rng, false, false);
	*/
	
	
	r = 80;
	c = 1000;

	Graph_2D(G, r, c, s, t);
	

	//file.open("Random_Graph.txt", std::ios::app);

	file.open("Test.txt", std::ios::app);


	try
	{
		if (!file)
		{
			throw std::runtime_error("Could not open file");
		}
	}
	catch (std::exception& ex)
	{
		std::cout << ex.what() << std::endl;
	}
	
	/*
	for (ei = edges(G).first; ei != edges(G).second; ++ei) {

		G[*ei].weight = rand() % 2 + 1;
	}

	*/

	s = random_vertex(G, rng);

	t = random_vertex(G, rng);


	CompareNodes::v_dist = get(&NodeInfo::distance, G);


	for (tie(vi, vi_end) = vertices(G); vi != vi_end; vi++) {

		CompareNodes::v_dist[*vi] = INT_MAX;

	}

	CompareNodes::v_dist[s] = 0;



	w_star = get(&EdgeInfo::weight, G);
	w_dijk = get(&EdgeInfo::weight, G);




	h = get(&NodeInfo::h_t, G);

	p = get(&NodeInfo::parent, G);




	clock_t start_t = clock();

	cnt = Dijkstra_SP(G, s, t, w_dijk, CompareNodes::v_dist, p, false);

	clock_t end_t = clock();


	elapsed = double(end_t - start_t) / CLOCKS_PER_SEC;


	file << "Time of Dijkstra for Grid Graph(" << r << "," << c << ")and weights in[1, 2] = "
		<< elapsed << std::endl;

	d_dist = CompareNodes::v_dist[t];


	if (d_dist == INT_MAX)
	{
		file << "No path from " << s << " to " << t << " exists in Graph" << std::endl;
	}

	else {

		file << "Distance in Dijkstra from source " << s << " to " << t << " is " << d_dist << std::endl;

		file << "Nodes Dijkstra accessed till reach target are " << cnt << std::endl;
	}

	for (tie(vi, vi_end) = vertices(G); vi != vi_end; vi++) {

		CompareNodes::v_dist[*vi] = INT_MAX;

	}

	CompareNodes::v_dist[s] = 0;

	h_euclidean(G,h, t,w_star);

	start_t = clock();

	cnt = A_Star(G, s, t, w_star, p, CompareNodes::v_dist, h);

	end_t = clock();


	elapsed = double(end_t - start_t) / CLOCKS_PER_SEC;

	file << std::endl;

	file << "Time of A* for	Grid Graph(" << r << "," << c << ")and weights in[1, 2] = "
		<< elapsed << std::endl;


	a_dist = CompareNodes::v_dist[t];

	if (d_dist == INT_MAX)
	{
		file << "No path from " << s << " to " << t << " exists in Graph" << std::endl;
	}
	else {

		file << "Distance in A* from source " << s << " to " << t << " is " << a_dist << std::endl;

		file << "Nodes A* accessed till reach target are " << cnt << std::endl;
	}


	if (a_dist == d_dist)
	{
		file << "Algos agreed" << std::endl;
	}

	file.close();

	return 0;


}