#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

struct Point {
    double x, y;
};

struct Ray {
    Point origin;
    Point direction;
};

struct Edge {
    Point start;
    Point end;
    int start_index;
    int end_index;
};

bool intersectRayWithEdge(const Ray& ray, const Edge& edge, Point& intersection) {
    double t = std::numeric_limits<double>::infinity();

    if (edge.start.y == edge.end.y) { // Horizontal segment
        if (ray.direction.y != 0) {
            t = (edge.start.y - ray.origin.y) / ray.direction.y;
            double x = ray.origin.x + t * ray.direction.x;
            if (t > 0 && std::min(edge.start.x, edge.end.x) <= x && x <= std::max(edge.start.x, edge.end.x)) {
                intersection = {x, edge.start.y};
                return true;
            }
        }
    } else { // Vertical segment
        if (ray.direction.x != 0) {
            t = (edge.start.x - ray.origin.x) / ray.direction.x;
            double y = ray.origin.y + t * ray.direction.y;
            if (t > 0 && std::min(edge.start.y, edge.end.y) <= y && y <= std::max(edge.start.y, edge.end.y)) {
                intersection = {edge.start.x, y};
                return true;
            }
        }
    }

    return false;
}

std::pair<Point, Edge> findClosestIntersection(const Ray& ray, const std::vector<std::pair<int, int>>& polygon, int start_index) {
    Point closestIntersection = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    Edge closestEdge = {{0, 0}, {0, 0}, 0, 0};
    double min_distance = std::numeric_limits<double>::infinity();

    for (int i = 0; i < polygon.size(); ++i) {
        int j = (i + 1) % polygon.size();
        Edge edge = {{static_cast<double>(polygon[i].first), static_cast<double>(polygon[i].second)}, 
                     {static_cast<double>(polygon[j].first), static_cast<double>(polygon[j].second)}, i, j};
        Point intersection;
        if (intersectRayWithEdge(ray, edge, intersection)) {
            double distance = std::sqrt(std::pow(intersection.x - ray.origin.x, 2) + std::pow(intersection.y - ray.origin.y, 2));
            if (distance < min_distance) {
                min_distance = distance;
                closestIntersection = intersection;
                closestEdge = edge;
            }
        }
    }

    return {closestIntersection, closestEdge};
}
