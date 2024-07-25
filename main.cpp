#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <limits>

// Dijkstra implementation code (provided later)
#include "Dijkstra.cpp"

// Definitions for findClosestIntersection and related structures (provided later)
#include "findClosestIntersection.cpp"

// Polygon structure and global variables
std::map<std::string, std::vector<std::vector<std::pair<int, int>>>> polygons;

// Utility function to calculate distance between two points
double calculateDistance(const std::pair<int, int>& a, const std::pair<int, int>& b) {
    return std::sqrt(std::pow(b.first - a.first, 2) + std::pow(b.second - a.second, 2));
}

// Function to determine if a point is convex
bool isConvex(const std::pair<int, int>& prev, const std::pair<int, int>& current, const std::pair<int, int>& next) {
    int cross_product = (next.first - current.first) * (prev.second - current.second) - (prev.first - current.first) * (next.second - current.second);
    return cross_product > 0;
}

// Function to parse polygons from a file
void parse_polygons(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    std::string current_key;
    std::vector<std::pair<int, int>> current_polygon;
    std::istringstream iss;
    bool in_polygon = false;

    while (std::getline(file, line)) {
        if (line[0] == '\'') {
            if (!current_key.empty()) {
                polygons[current_key].push_back(current_polygon);
                current_polygon.clear();
            }
            current_key = line.substr(1, line.find('\'', 2) - 1);
        } else if (line[0] == 'p') {
            if (in_polygon) {
                polygons[current_key].push_back(current_polygon);
                current_polygon.clear();
            }
            in_polygon = true;
        } else if (in_polygon) {
            iss.clear();
            iss.str(line);
            int x, y;
            if (iss >> x >> y) {
                current_polygon.push_back({x, y});
            }
        }
    }

    if (!current_key.empty() && !current_polygon.empty()) {
        polygons[current_key].push_back(current_polygon);
    }
}

void createEdgesAndCalculateDistances(const std::vector<std::pair<int, int>>& polygon, Graph& graph) {
    int n = polygon.size();
    for (int i = 0; i < n; ++i) {
        int prev = (i - 1 + n) % n;
        int next = (i + 1) % n;
        bool is_convex = isConvex(polygon[prev], polygon[i], polygon[next]);
        double edge_cost = calculateDistance(polygon[i], polygon[next]);
        
        if (is_convex) {
            graph.addEdge(i, next, edge_cost);
        } else {
            graph.addEdge(i, next, edge_cost);
            // Find and add extra edges
            Ray ray1 = {polygon[prev], {polygon[i].first - polygon[prev].first, polygon[i].second - polygon[prev].second}};
            Ray ray2 = {polygon[next], {polygon[i].first - polygon[next].first, polygon[i].second - polygon[next].second}};
            
            // Closest intersection points and edges
            Point intersection1;
            Edge closestEdge1;
            std::tie(intersection1, closestEdge1) = findClosestIntersection(ray1, polygon, i);
            
            Point intersection2;
            Edge closestEdge2;
            std::tie(intersection2, closestEdge2) = findClosestIntersection(ray2, polygon, i);
            
            // Adding extra edges for concave points
            graph.addEdge(i, (i + 1) % polygon.size(), edge_cost);

            // Add edge to the first intersection point
            double intersection_cost1 = calculateDistance(polygon[i], {intersection1.x, intersection1.y});
            graph.addEdge(i, n, intersection_cost1);
            
            // Add edges from the intersection point to the closest edge's start and end
            double closestEdgeStart_cost1 = calculateDistance({intersection1.x, intersection1.y}, {closestEdge1.start.x, closestEdge1.start.y});
            double closestEdgeEnd_cost1 = calculateDistance({intersection1.x, intersection1.y}, {closestEdge1.end.x, closestEdge1.end.y});
            graph.addEdge(n, closestEdge1.start_index, closestEdgeStart_cost1);
            graph.addEdge(n, closestEdge1.end_index, closestEdgeEnd_cost1);
            
            // Repeat for the second intersection point
            double intersection_cost2 = calculateDistance(polygon[i], {intersection2.x, intersection2.y});
            graph.addEdge(i, n + 1, intersection_cost2);
            
            double closestEdgeStart_cost2 = calculateDistance({intersection2.x, intersection2.y}, {closestEdge2.start.x, closestEdge2.start.y});
            double closestEdgeEnd_cost2 = calculateDistance({intersection2.x, intersection2.y}, {closestEdge2.end.x, closestEdge2.end.y});
            graph.addEdge(n + 1, closestEdge2.start_index, closestEdgeStart_cost2);
            graph.addEdge(n + 1, closestEdge2.end_index, closestEdgeEnd_cost2);
        }
    }
}

void export_geojson(const std::string& path) {
    std::ofstream out(path);
    out << "{\n";
    out << "  \"type\": \"FeatureCollection\",\n";
    out << "  \"features\": [\n";

    bool first = true;
    for (const auto& pair : polygons) {
        if (!first) {
            out << ",\n";
        }
        out << "    {\n";
        out << "      \"type\": \"Feature\",\n";
        out << "      \"properties\": {\n";
        out << "        \"name\": \"" << pair.first << "\"\n";
        out << "      },\n";
        out << "      \"geometry\": {\n";
        out << "        \"type\": \"MultiPolygon\",\n";
        out << "        \"coordinates\": [\n";

        bool first_polygon = true;
        for (const auto& polygon : pair.second) {
            if (!first_polygon) {
                out << ",\n";
            }
            out << "          [\n";
            bool first_point = true;
            for (const auto& point : polygon) {
                if (!first_point) {
                    out << ",\n";
                }
                out << "            [" << point.first << ", " << point.second << "]";
                first_point = false;
            }
            out << "          ]\n";
            first_polygon = false;
        }

        out << "        ]\n";
        out << "      }\n";
        out << "    }";
        first = false;
    }

    out << "\n  ]\n";
    out << "}";
}

int main() {
    float start = clock();
    std::string file_path = "CORE_TILE_PERI.drc.results"; 
    parse_polygons(file_path);

    std::string output_path = "output.geojson";
    //export_geojson(output_path);

    std::cout << "GeoJSON has been exported to " << output_path << std::endl;
    std::cout << "M1_polygons_size " << polygons["M1_polygons"].size() << std::endl;

    // Initialize graph
    Graph graph;
    graph.initialize(polygon.size());

    // Create edges and calculate distances
    for (const auto& pair : polygons) {
        std::cout << "Key: " << pair.first << std::endl;
        std::cout << "Key_Length " << pair.second.size() << std::endl;
        for (const auto& polygon : pair.second) {
            createEdgesAndCalculateDistances(polygon, graph);
        }
    }

    std::cout << "Time " << clock() - start << std::endl;
    return 0;
}
