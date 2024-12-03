// Project Identifier: 3E33912F8BAA7542FC4A1585D2DB6FE0312725B9
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include <getopt.h>

using namespace std;

struct NodeA {
    int x;
    int y;
    enum Area { WILD, SAFE, BOARDER } area;
    // Constructor for convenience
    NodeA(int x, int y, Area area)
        : x(x)
        , y(y)
        , area(area) {}
};

struct Node {
    int x;
    int y;
    Node(int x, int y)
        : x(x)
        , y(y) {}
};

struct PrimData {
    double d;    // lowest distance to MST
    int32_t p;   // index of vertex's parent in MST
    bool k;      // whether or not vertex is already in MST

    PrimData()
        : d(numeric_limits<double>::infinity())
        , p(-1)
        , k(false) {}
};

static inline double calculateDistanceMST(const NodeA& a, const NodeA& b);


string parseCommandLine(int argc, char* argv[]);


string parseCommandLine(int argc, char* argv[]) {
    opterr = false;   // Disable automatic error printing by getopt
    int choice;
    int option_index = 0;
    option long_options[] = {
        { "mode", required_argument, nullptr,  'm'},
        { "help",       no_argument, nullptr,  'h'},
        {nullptr,                 0, nullptr, '\0'}
    };

    string mode;
    bool mode_set = false;
    while ((choice = getopt_long(argc, argv, "m:h", long_options, &option_index)) != -1) {
        switch (choice) {
        case 'm':
            mode = optarg;
            mode_set = true;
            break;
        case 'h':
            cout << "Usage: ./zoo --mode {MST|FASTTSP|OPTTSP} \n";
            cout << "Options:\n";
            cout << " -m, --mode    Select mode: MST, FASTTSP, or OPTTSP\n";
            cout << " -h, --help    Display this help message and exit\n";
            exit(0);
        default:
            cerr << "Unknown command line option\n";
            exit(1);
        }
    }

    // Check if mode was set
    if (!mode_set) {
        cerr << "Error: Mode not set. Use -m or --mode to set the mode.\n";
        exit(1);
    }

    // Additional validation for mode value
    if (mode != "MST" && mode != "FASTTSP" && mode != "OPTTSP") {
        cerr << "Error: Invalid mode. Use one of MST, FASTTSP, or OPTTSP.\n";
        exit(1);
    }

    return mode;
}


static inline double calculateDistanceMST(const NodeA& a, const NodeA& b) {
    if ((a.area == NodeA::WILD && b.area == NodeA::SAFE) || (a.area == NodeA::SAFE && b.area == NodeA::WILD)) {
        return numeric_limits<double>::infinity();
    }

    long long dx = a.x - b.x;
    long long dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}


NodeA::Area categorizeArea(int x, int y) {
    if ((x == 0 && y <= 0) || (x <= 0 && y == 0)) {
        return NodeA::BOARDER;
    }
    if (x < 0 && y < 0) {
        return NodeA::WILD;
    }
    return NodeA::SAFE;
}

vector<NodeA> readNodeMST() {
    int numberNodes = 0;
    cin >> numberNodes;
    vector<NodeA> nodes;
    nodes.reserve(static_cast<size_t>(numberNodes));
    for (int i = 0; i < numberNodes; ++i) {
        int x;
        int y;
        cin >> x >> y;
        nodes.emplace_back(x, y, categorizeArea(x, y));
    }
    return nodes;
}

vector<PrimData> constructMST(vector<NodeA>& nodes) {
    vector<PrimData> prim_table(nodes.size());
    size_t start_vertex = 0;
    prim_table[start_vertex].d = 0;
    // size_t firstUnvisited = 0;

    while (true) {
        // find the vertex with the smallest distance to the MST
        double min_distance = numeric_limits<double>::infinity();
        int min_vertex = -1;
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (!prim_table[i].k && prim_table[i].d < min_distance) {
                min_distance = prim_table[i].d;
                min_vertex = static_cast<int>(i);
            }
        }

        if (min_vertex == -1) {
            break;
        }   // all vertices are in the MST

        // add min_vertex to the MST
        prim_table[static_cast<size_t>(min_vertex)].k = true;

        // update the distance to the MST for each vertex adjacent to min_vertex
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (!prim_table[i].k && static_cast<int>(i) != min_vertex) {
                double distance = calculateDistanceMST(nodes[static_cast<size_t>(min_vertex)], nodes[i]);
                if (distance < prim_table[i].d) {
                    prim_table[i].d = distance;
                    prim_table[i].p = min_vertex;
                }
            }
        }
    }
    // In the following optimization, "k" is useless.


    return prim_table;
}

void printMST(const vector<PrimData>& mst) {
    double totalWeight = 0.0;
    bool mstPossible = true;

    // Check if MST is possible (i.e., every node has a parent except the starting node)
    for (const auto& data : mst) {
        if (data.p == -1 && data.d != 0.0) {   // d != 0.0 to allow starting node
            mstPossible = false;
            break;
        }
    }

    if (!mstPossible) {
        cerr << "Cannot construct MST\n";
        exit(1);
    }

    for (size_t i = 1; i < mst.size(); ++i) {
        totalWeight += mst[i].d;
    }

    cout << fixed << setprecision(2) << totalWeight << "\n";

    for (size_t i = 1; i < mst.size(); ++i) {
        int a = static_cast<int>(i);
        int b = mst[i].p;
        if (a > b) {
            swap(a, b);
        }
        cout << a << " " << b << "\n";
    }
}

struct PrimDataOPT {
    double d;   // lowest distance to MST
    PrimDataOPT()
        : d(numeric_limits<double>::infinity()) {}
};


// Optimal TSP
class Graph {
private:
    vector<Node> nodes;
    vector<size_t> bestTour;
    vector<size_t> currTour;

    double bestTourCost;
    double curCost;
    double upper;

public:
    bool print = false;
    Graph(const vector<Node>& nodes)
        : nodes(nodes)
        , bestTourCost(numeric_limits<double>::infinity())
        , curCost(0)
        , upper(0) {
        bestTour.reserve(nodes.size());
        currTour.reserve(nodes.size());
    }

    double constructMST(vector<Node>& nodes) {
        vector<PrimDataOPT> prim_table(nodes.size());
        size_t start_vertex = 0;
        prim_table[start_vertex].d = 0;
        size_t firstUnvisited = 0;

        // In the following optimization, "k" is useless.
        while (firstUnvisited < nodes.size()) {
            // find the vertex with the smallest distance to the MST
            double min_distance = numeric_limits<double>::infinity();
            // int min_vertex = -1;
            size_t min_vertex = 0;
            for (size_t i = firstUnvisited; i < nodes.size(); ++i) {
                if (prim_table[i].d < min_distance) {
                    min_distance = prim_table[i].d;
                    min_vertex = i;
                }
            }

            // add min_vertex to the MST
            swap(prim_table[firstUnvisited], prim_table[min_vertex]);
            swap(nodes[firstUnvisited], nodes[min_vertex]);
            min_vertex = firstUnvisited;
            ++firstUnvisited;

            // update the distance to the MST for each vertex adjacent to min_vertex
            for (size_t i = firstUnvisited; i < nodes.size(); ++i) {
                double dx = nodes[min_vertex].x - nodes[i].x;
                double dy = nodes[min_vertex].y - nodes[i].y;
                double distance = dx * dx + dy * dy;
                if (distance < prim_table[i].d) {
                    prim_table[i].d = distance;
                }
            }
        }
        double totalWeight = 0.0;
        for (size_t i = 1; i < prim_table.size(); ++i) {
            totalWeight += sqrt(prim_table[i].d);
        }
        return totalWeight;
    }

    void cheapestInsertionTSP(bool print);
    void cheapestInsertion(bool print);

    static inline double calculateDistance(const Node& a, const Node& b) {
        auto dx = static_cast<double>(a.x - b.x);
        auto dy = static_cast<double>(a.y - b.y);
        return sqrt(dx * dx + dy * dy);
    }

    static inline double insertionCost(const Node& i, const Node& j, const Node& k) {
        return calculateDistance(j, i) + calculateDistance(i, k) - calculateDistance(j, k);
    }


    bool promising(const vector<size_t>& path, size_t permLength) {
        // Implement your logic to decide if the path is promising
        if (path.size() - permLength < 5) {
            return true;
        }
        // First, calculate the current cost, it is just curCost

        // Second, calculate cost of the MST connecting the remaining points not in the partial solution.
        vector<Node> remainingNodes;
        vector<bool> visited(nodes.size(), false);
        for (size_t i = 0; i < permLength; ++i) {
            visited[path[i]] = true;
        }
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (!visited[i]) {
                remainingNodes.push_back(nodes[i]);
            }
        }
        double mstCost = constructMST(remainingNodes);

        // Third, calculate the two connection lines
        double connectionCost = 0.0;
        if (remainingNodes.empty()) {
            connectionCost = 0.0;
        } else {
            Node& startNode = nodes[path.front()];
            Node& endNode = nodes[path[permLength - 1]];

            double minDistToStart = numeric_limits<double>::infinity();
            double minDistToEnd = numeric_limits<double>::infinity();

            for (const auto& remainingNode : remainingNodes) {
                double distToStart = calculateDistance(startNode, remainingNode);
                double distToEnd = calculateDistance(endNode, remainingNode);
                minDistToStart = min(minDistToStart, distToStart);
                minDistToEnd = min(minDistToEnd, distToEnd);
            }
            connectionCost = minDistToEnd + minDistToStart;
        }

        // Fourth, compare
        double lowerBound = curCost + mstCost + connectionCost;
        return lowerBound <= bestTourCost;
    }


    void genPerms(vector<size_t>& path, size_t permLength) {
        if (permLength == path.size()) {
            curCost += calculateDistance(nodes[path.back()], nodes[path.front()]);
            if (curCost <= bestTourCost) {
                bestTourCost = curCost;
                bestTour = path;
            }
            curCost -= calculateDistance(nodes[path.back()], nodes[path.front()]);
            return;
        }

        if (!promising(path, permLength)) {
            return;
        }

        for (size_t i = permLength; i < path.size(); ++i) {
            swap(path[permLength], path[i]);
            curCost += calculateDistance(nodes[path[permLength - 1]], nodes[path[permLength]]);
            genPerms(path, permLength + 1);
            curCost -= calculateDistance(nodes[path[permLength - 1]], nodes[path[permLength]]);
            swap(path[permLength], path[i]);
        }
    }

    void solve() {
        print = false;
        bestTourCost = 0.0;
        cheapestInsertion(print);
        bestTourCost = upper;
        genPerms(currTour, 1);
        cout << bestTourCost << "\n";
        reverse(bestTour.begin() + 1, bestTour.end());
        for (auto node : bestTour) {
            cout << node << " ";
        }   // Start from index 1 as the first node is fixed
    }
};

// Cheapest Insertion: TSP Heuristics

void Graph::cheapestInsertion(bool print) {
    vector<size_t> tour = { 0, 0 };
    vector<bool> visited(nodes.size(), false);
    visited[0] = true;

    double minCost = std::numeric_limits<double>::max();
    size_t minVertex = 0;
    size_t minPos = 0;

    for (size_t k = 0; k < nodes.size(); ++k) {
        minCost = std::numeric_limits<double>::max();
        if (!visited[k]) {
            for (size_t j = 0; j < tour.size() - 1; ++j) {
                double cost = calculateDistance(nodes[tour[j]], nodes[k])
                            + calculateDistance(nodes[tour[j + 1]], nodes[k])
                            - calculateDistance(nodes[tour[j]], nodes[tour[j + 1]]);
                if (cost < minCost) {
                    minCost = cost;
                    minVertex = k;
                    minPos = j + 1;
                }
            }
        }
        tour.insert(tour.begin() + static_cast<int>(minPos), minVertex);
        visited[minVertex] = true;
    }

    tour.pop_back();
    double sum = 0;
    for (size_t i = tour.size() - 1; i > 0; --i) {
        sum = sum + calculateDistance(nodes[tour[i]], nodes[tour[i - 1]]);
    }
    upper = sum;

    for (size_t i = 0; i < tour.size() - 1; i++) {
        bestTour.emplace_back(tour[i]);
        currTour.emplace_back(tour[i]);
    }
    reverse(bestTour.begin() + 1, bestTour.end());
    reverse(currTour.begin() + 1, currTour.end());
    if (print) {
        cout << sum << "\n";
        // Print the order of nodes
        for (size_t i = 0; i < tour.size() - 1; i++) {
            cout << tour[i] << " ";
        }
        cout << "\n";
    }
}


int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);
    cout << fixed;
    cout << setprecision(2);
    string mode = parseCommandLine(argc, argv);

    if (mode == "MST") {
        vector<NodeA> nodes = readNodeMST();
        vector<PrimData> mst = constructMST(nodes);
        printMST(mst);
    } else {
        int numberNodes = 0;
        cin >> numberNodes;
        vector<Node> nodes;
        nodes.reserve(static_cast<size_t>(numberNodes));
        for (int i = 0; i < numberNodes; ++i) {
            int x = 0;
            int y = 0;
            cin >> x >> y;
            nodes.emplace_back(x, y);
        }
        Graph graph(nodes);
        if (mode == "FASTTSP") {
            graph.print = true;
            graph.cheapestInsertion(graph.print);
        }
        if (mode == "OPTTSP") {
            graph.solve();
        }
    }

    return 0;
}
