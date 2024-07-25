#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <limits>
#include "Graph.h"  // Include the Graph header

// Dijkstra implementation is now part of Graph class

// Definitions for findClosestIntersection and related structures
#include "findClosestIntersection.cpp"
#include "Graph.cpp"

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
            std::cout << "is_convex " << std::endl;
            graph.addEdge(i, next, edge_cost);
            
        } else {
            std::cout << "is_concave " << std::endl;
            graph.addEdge(i, next, edge_cost);
            // Find and add extra edges
            Ray ray1 = {{(double)polygon[prev].first, (double)polygon[prev].second}, {(double)(polygon[i].first - polygon[prev].first), (double)(polygon[i].second - polygon[prev].second)}};
            Ray ray2 = {{(double)polygon[next].first, (double)polygon[next].second}, {(double)(polygon[i].first - polygon[next].first), (double)(polygon[i].second - polygon[next].second)}};
            
            // Closest intersection points and edges
            Point intersection1;
            Edge closestEdge1;
            std::cout << "findClosestIntersection1 start " << std::endl;
            std::tie(intersection1, closestEdge1) = findClosestIntersection(ray1, polygon, i);
            std::cout << "findClosestIntersection1 finished " << intersection1.x << closestEdge1.start.x << closestEdge1.end.x << std::endl;
            Point intersection2;
            Edge closestEdge2;
            std::tie(intersection2, closestEdge2) = findClosestIntersection(ray2, polygon, i);
            std::cout << "Adding extra edges for concave points" << std::endl;
            // Adding extra edges for concave points
            graph.addEdge(i, (i + 1) % polygon.size(), edge_cost);

            std::cout << " Add edge to the first intersection point" << std::endl;
            // Add edge to the first intersection point
            double intersection_cost1 = calculateDistance(polygon[i], {intersection1.x, intersection1.y});
            std::cout << " calculateDistance finished" << intersection_cost1 << std::endl;
            std::cout << " i n " << i << " " << n << std::endl;
            graph.addEdge(i, n, intersection_cost1);
            
            std::cout << " Add edges from the intersection point to the closest edge's start and end" << std::endl;
            // Add edges from the intersection point to the closest edge's start and end
            double closestEdgeStart_cost1 = calculateDistance({intersection1.x, intersection1.y}, {closestEdge1.start.x, closestEdge1.start.y});
            double closestEdgeEnd_cost1 = calculateDistance({intersection1.x, intersection1.y}, {closestEdge1.end.x, closestEdge1.end.y});
            graph.addEdge(n, closestEdge1.start_index, closestEdgeStart_cost1);
            graph.addEdge(n, closestEdge1.end_index, closestEdgeEnd_cost1);
            
            std::cout << " Repeat for the second intersection point" << std::endl;
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
            out << "            [\n";

            bool first_point = true;
            for (const auto& point : polygon) {
                if (!first_point) {
                    out << ",\n";
                }
                out << "              [" << point.first << ", " << point.second << "]";
                first_point = false;
            }

            out << "\n            ]\n";
            out << "          ]";
            first_polygon = false;
        }

        out << "\n        ]\n";
        out << "      }\n";
        out << "    }";
        first = false;
    }

    out << "\n  ]\n";
    out << "}\n";
    out.close();
}

int main() {
    // Read the input file and parse polygons
    parse_polygons("CORE_TILE_PERI.drc.results");

    // Create graph and calculate distances for each polygon
    for (const auto& pair : polygons) {
        std::cout << "Key: " << pair.first << std::endl;
        std::cout << "Key_Length " << pair.second.size() << std::endl;
        for (const auto& polygon : pair.second) {
            for (const auto& point : polygon) {
        // const auto& polygon = pair.second[0]; // Assuming one polygon per key
                // Graph graph(polygon.size());
                Graph graph(1000);
                std::cout << "polygon.size() " << polygon.size() << std::endl;
                createEdgesAndCalculateDistances(polygon, graph);
                std::cout << "polygon.size() " << polygon.size() << std::endl;

                // Run Dijkstra's algorithm on the graph
                std::vector<double> distMin;
                graph.dijkstra(0, distMin);
                std::cout << "polygon.size() " << polygon.size() << std::endl;

                // Print or use the results as needed
                for (double dist : distMin) {
                    std::cout << dist << " ";
                }
                std::cout << "polygon.size() " << polygon.size() << std::endl;

                std::cout << std::endl;
            }
        }
    }

    // Export polygons to GeoJSON
    export_geojson("output.geojson");

    return 0;
}
