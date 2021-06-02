//#include "adjList.h"
#include "graph.h"
const Vertex Graph::InvalidVertex = "_CS225INVALIDVERTEX";
const int Graph::InvalidWeight = INT_MIN;
const string Graph::InvalidLabel = "_CS225INVALIDLABEL";
const Edge Graph::InvalidEdge = Edge(Graph::InvalidVertex, Graph::InvalidVertex, Graph::InvalidWeight, Graph::InvalidLabel);
Graph::Graph() {
    
}
vector<Vertex> Graph::getAdjacent(Vertex source) const 
{
    auto lookup = adjacency_list.find(source);

    if(lookup == adjacency_list.end())
        return vector<Vertex>();

    else
    {
        vector<Vertex> vertex_list;
        unordered_map <Vertex, Edge> & map = adjacency_list[source];
        for (auto it = map.begin(); it != map.end(); it++)
        {
            vertex_list.push_back(it->first);
        }
        return vertex_list;
    }
}


Vertex Graph::getStartingVertex() const
{
    return adjacency_list.begin()->first;
}

vector<Vertex> Graph::getVertices() const
{
    vector<Vertex> ret;

    for(auto it = adjacency_list.begin(); it != adjacency_list.end(); it++)
    {
        ret.push_back(it->first);
    }

    return ret;
}

Edge Graph::getEdge(Vertex source , Vertex destination) const
{
    if(assertEdgeExists(source, destination, __func__) == false)
        return Edge();
    Edge ret = adjacency_list[source][destination];
    return ret;
}
vector<Edge> Graph::getEdges() const
{
    if (adjacency_list.empty())
        return vector<Edge>();

    vector<Edge> ret;
    set<pair<Vertex, Vertex>> seen;

    for (auto it = adjacency_list.begin(); it != adjacency_list.end(); it++)
    {
        Vertex source = it->first;
        for (auto its = adjacency_list[source].begin(); its != adjacency_list[source].end(); its++)
        {
            Vertex destination = its->first;
            if(seen.find(make_pair(source, destination)) == seen.end())
            {
                //this pair is never added to seen
                ret.push_back(its->second);
                seen.insert(make_pair(source,destination));
                //if(!directed)
                //{
                    seen.insert(make_pair(destination, source));
                //}
            }
        }
    }

    return ret;
}

bool Graph::vertexExists(Vertex v) const
{
    return assertVertexExists(v, "");
}

bool Graph::edgeExists(Vertex source, Vertex destination) const
{
    return assertEdgeExists(source, destination, "");
}

Edge Graph::setEdgeLabel(Vertex source, Vertex destination, string label)
{
    if (assertEdgeExists(source, destination, __func__) == false)
        return InvalidEdge;
    Edge e = adjacency_list[source][destination];
    Edge new_edge(source, destination, e.getWeight(), label);
    adjacency_list[source][destination] = new_edge;

    /*if(!directed)
    {
        Edge new_edge_reverse(destination,source, e.getWeight(), label);
        adjacency_list[destination][source] = new_edge_reverse;
    }*/
    return new_edge;
}


string Graph::getEdgeLabel(Vertex source, Vertex destination) const
{
    if(assertEdgeExists(source, destination, __func__) == false)
        return InvalidLabel;
    return adjacency_list[source][destination].getLabel();
}

int Graph::getEdgeWeight(Vertex source, Vertex destination) const
{
    if (!weighted)
        error("can't get edge weights on non-weighted graphs!");

    if(assertEdgeExists(source, destination, __func__) == false)
        return InvalidWeight;
    return adjacency_list[source][destination].getWeight();
}

void Graph::insertVertex(Vertex v)
{
    // will overwrite if old stuff was there
    removeVertex(v);
    // make it empty again
    adjacency_list[v] = unordered_map<Vertex, Edge>();
}


Vertex Graph::removeVertex(Vertex v)
{

    if (adjacency_list.find(v) != adjacency_list.end())
    {
        /*if(!directed){
            for (auto it = adjacency_list[v].begin(); it != adjacency_list[v].end(); it++)
            {
                Vertex u = it->first;
                adjacency_list[u].erase(v); 
            }
            adjacency_list.erase(v);
            return v;
        }*/
        
        adjacency_list.erase(v);
        for(auto it2 = adjacency_list.begin(); it2 != adjacency_list.end(); it2++)
        {
            Vertex u = it2->first;
            if (it2->second.find(v)!=it2->second.end())
            {
                it2->second.erase(v);
            }
        }
        return v;
    }

    return InvalidVertex;
}

bool Graph::insertEdge(Vertex source, Vertex destination)
{
    if(adjacency_list.find(source)!= adjacency_list.end() 
    && adjacency_list[source].find(destination)!= adjacency_list[source].end())
    {
        //edge already exit
        return false;
    }

    if(adjacency_list.find(source)==adjacency_list.end())
    {
        adjacency_list[source] = unordered_map<Vertex, Edge>();
    }
        //source vertex exists
    adjacency_list[source][destination] = Edge(source, destination);
    //if(!directed) {
        if(adjacency_list.find(destination)== adjacency_list.end())
        {
            adjacency_list[destination] = unordered_map<Vertex, Edge>();
        }
        adjacency_list[destination][source] = Edge(source, destination);
    //}
    
    return true;
}

Edge Graph::removeEdge(Vertex source, Vertex destination)
{
    if(assertEdgeExists(source, destination, __func__) == false)
        return InvalidEdge;
    Edge e = adjacency_list[source][destination];
    adjacency_list[source].erase(destination);
    // if undirected, remove the corresponding edge
    //if(!directed){
        adjacency_list[destination].erase(source);
   // }
    return e;
}


Edge Graph::setEdgeWeight(Vertex source, Vertex destination, int weight)
{
    if (assertEdgeExists(source, destination, __func__) == false)
        return InvalidEdge;
    Edge e = adjacency_list[source][destination];
    //std::cout << "setting weight: " << weight << std::endl;
    Edge new_edge(source, destination, weight, e.getLabel());
    adjacency_list[source][destination] = new_edge;

   // if(!directed)
     //   {
            Edge new_edge_reverse(destination,source, weight, e.getLabel());
            adjacency_list[destination][source] = new_edge_reverse;
       // }

    return new_edge;
}

bool Graph::assertVertexExists(Vertex v, string functionName) const
{
    if (adjacency_list.find(v) == adjacency_list.end())
    {
        if (functionName != "")
            error(functionName + " called on nonexistent vertices");
        return false;
    }
    return true;
}

bool Graph::assertEdgeExists(Vertex source, Vertex destination, string functionName) const
{
    if(assertVertexExists(source,functionName) == false)
        return false;
    if(adjacency_list[source].find(destination)== adjacency_list[source].end())
    {
        if (functionName != "")
            error(functionName + " called on nonexistent edge " + source + " -> " + destination);
        return false;
    }

    if(!directed)
    {
        if (assertVertexExists(destination,functionName) == false)
            return false;
        if(adjacency_list[destination].find(source)== adjacency_list[destination].end())
        {
            if (functionName != "")
                error(functionName + " called on nonexistent edge " + destination + " -> " + source);
            return false;
        }
    }
    return true;
}

void Graph::error(string message) const
{
    cerr << "\033[1;31m[Graph Error]\033[0m " + message << endl;
}
float Graph::toRadians(float point) {
    return (M_PI / 180) * point;
}

int Graph::distance(tuple<float, float> v1, tuple<float, float> v2) {
    float v1Lat = toRadians(get<0>(v1));
    float v1Lang = toRadians(get<1>(v1));
    float v2Lat = toRadians(get<0>(v2));
    float v2Lang = toRadians(get<1>(v2));
    float dLang = v2Lang - v1Lang;
    float dLat = v2Lat- v1Lat;
    float ans = pow(sin(dLat / 2), 2) + cos(v1Lat) * cos(v2Lat) *  pow(sin(dLang / 2), 2); 
    ans = 2 * asin(sqrt(ans)); 
    ans = ans * 3956;
    return (int) ans;
   
}
void Graph::BFS(Vertex source) {
    std::ofstream outfile;
    outfile.open("BFS.txt");
    int numVerticies = getVertices().size();
    bool *visited = new bool[numVerticies];
	
	
    for (int i = 0; i < numVerticies; i++) {
        visited[i] = false;
    }
    list<Vertex> queue;
    visited[stoi(source)] = true;
    queue.push_back(source);
	//int j = 0;
	while (!queue.empty()) {
        source = queue.front();
        cout<< source << " ";
        outfile << source << " ";
        queue.pop_front();
		unordered_map<Vertex, Edge> alist = adjacency_list[source];
		//might need to use an iterator (even if it is just to see if it works [0 -> 4 -> 8])
		for(const auto& pair : alist ) {
			//if (j < 3) //only print the first three things so it is not cluttered
			  //std::cout << "Vertex:[" << pair.first << "] Edge:[" << pair.second.source << "," << pair.second.dest << "]\n";
			
			
			int vert = stoi(pair.second.dest);
			if (!visited[vert]) {
                visited[vert] = true;
                queue.push_back(pair.second.dest);
            }
        }
		//j++;
		
    }
    outfile.close();
}
void Graph::dijkstra(Vertex source) {
    std::ofstream outfile;
    outfile.open("dijkstra.txt");

    set<pair<int, Vertex>> processedVertices;
    vector<int> dist(getVertices().size(), INF);
	
	vector<int> parent(getVertices().size(), -1);
	
    processedVertices.insert(make_pair(0, source));
    dist[stoi(source)] = 0;
	
    while (processedVertices.size() > 0) {
	//for (unsigned j = 0; j < processedVertices.size(); j++) {
        pair<int, Vertex> tmp = *(processedVertices.begin());
        processedVertices.erase(processedVertices.begin());
        int u = stoi(tmp.second); //vertex label stored in second of pair, within a certain set, gets second value of a pair which is the vertex
        //list<pair<int, Vertex>>::iterator i;
        outfile<< "\n" << u << "\n";
        std::cout << "\n" << u << "\n";

        unordered_map<Vertex, Edge> alist = adjacency_list[to_string(u)];
        for(const auto& pair : alist ) {
            int v = stoi(pair.first);
            int weight = (pair.second).getWeight();
            if (dist[v] > dist[u] + weight) {
                if (dist[v] != INF) {
                    processedVertices.erase(processedVertices.find(make_pair(dist[v], to_string(v))));
                }
                dist[v] = dist[u] + weight;
				parent[v] = u;
				outfile<< "\t" << u << " ~ " << parent[v] << "\n";
                std::cout << "\t" << u << " ~ " << parent[v] << "\n";

                processedVertices.insert(make_pair(dist[v], to_string(v)));
            }
        }
        path.clear();
        outfile << "Vertex\tDistance from Source\tParent\n";
        printf("Vertex   Distance from Source   Parent\n"); 
        for (unsigned i = 0; i < getVertices().size(); ++i) {
            path.insert(make_pair(to_string(i), dist[i]));
            //printf("%d \t\t %d\n", i, dist[i]); 
            outfile << i << "\t\t" << dist[i] << "\t\t" << parent[i] << std::endl;
            printf("%d \t\t %d \t\t %d\n", i, dist[i], parent[i]);

        }
    }
	parentVerts = parent;
    outfile.close();
}

void Graph::printPath(Vertex src, Vertex dst,ofstream& outfile_) {
		vector<int> path;
		int v = stoi(dst);
		while (v != stoi(src)) {
			path.push_back(v);
			v = parentVerts[v];
		}
		
		std::reverse(path.begin(), path.end());
		
		std::cout<< "Path from source : " << src <<  " to destination : "  << dst << "\n";
        outfile_ << "Path from source : " << src << " to destination : " << dst << "\n";
		std::cout<<src << " -> ";
        outfile_ << src << " -> ";
		for (unsigned int i = 0; i < path.size(); i++) {
			std::cout<< path[i] << " -> "; 
            outfile_<< path[i] << " -> ";
		}
		std::cout<<"done \n";
		
}

void Graph::landmark(Vertex source, Vertex destination, Vertex landmark) {
    std::ofstream outfile;
    outfile.open("landmark.txt");
    dijkstra(source);
	printPath(source, landmark,outfile);
    int sourceToLandmark = path[landmark];
    cout<<"sourceToLandmark = "<< sourceToLandmark<<endl;
    outfile << "sourceToLandmark = " << sourceToLandmark << endl;
    dijkstra(landmark);
	printPath(landmark, destination,outfile);
    int landmarkToDest = path[destination];
    cout<<"landmarkToDest = "<< landmarkToDest<<endl;
    outfile << "landmarkToDest = " << landmarkToDest << endl;
    int shortestPath = sourceToLandmark + landmarkToDest;
    cout<<"Shortest Path from "<< source<< "  to "<< destination<< " is " << shortestPath<<endl;
    outfile << "Shortest Path from " << source << "  to " << destination << " is " << shortestPath << endl;
    outfile.close();
}

void Graph::txtToGraph(std::string nodesFilename, std::string edgesFilename)
{
    std::ifstream infile(nodesFilename);
    std::string line;
    int i = 0;
    std::ofstream myfile;
    myfile.open("nodes.txt");
    std::map<std::string, std::tuple<float, float>> dict;

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        int userid;
        float lat, lang;
        std::string checkin, locationId;

        if (!(iss >> userid >> checkin >> lat >> lang >> locationId)) { break; }
        if (abs(lat) < 0.0001 || abs(lang) < 0.0001) {
            if (!(iss >> userid >> checkin >> lat >> lang >> locationId)) { break; }
        }
        if (abs(lat) < 0.0001 || abs(lang) < 0.0001) {
            if (!(iss >> userid >> checkin >> lat >> lang >> locationId)) { break; }
        }
        if (userid >= i)
        {
            std::string user = std::to_string(userid);
            insertVertex(user);
            dict.insert(std::make_pair(user, std::make_tuple(lat, lang)));
            myfile << userid << " " << lat << " " << lang << "\n";
            i = userid + 1;
        }

        if (i > 5)
        {
            break;
        }
    }
    myfile.close();

    std::ifstream input(edgesFilename);

    while (std::getline(input, line))
    {
        std::istringstream stream(line);
        int startUser, endUser;
        if (!(stream >> startUser >> endUser)) { break; }
        insertEdge(std::to_string(startUser), std::to_string(endUser));
        int dist = distance(dict[std::to_string(startUser)], dict[std::to_string(endUser)]);
        setEdgeWeight(std::to_string(startUser), std::to_string(endUser), dist);
    }
}
