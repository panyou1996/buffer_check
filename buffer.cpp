#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <queue>
#include <map>
#include <functional>
#include <algorithm>
#include <cassert>

#define EPSILON 0.00001

class Point {
public:
    double x, y;
    Point(double x=0, double y=0) : x(x), y(y) {}

    bool operator==(const Point& other) const {
        return std::fabs(x - other.x) < EPSILON && std::fabs(y - other.y) < EPSILON;
    }

    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }

    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y);
    }

    Point operator*(double scalar) const {
        return Point(x * scalar, y * scalar);
    }

    double distance(const Point& other) const {
        return std::hypot(x - other.x, y - other.y);
    }
};

class Vector {
public:
    double x, y;
    Vector(double x=0, double y=0) : x(x), y(y) {}

    Vector operator-() const {
        return Vector(-x, -y);
    }

    Vector operator+(const Vector& other) const {
        return Vector(x + other.x, y + other.y);
    }

    double dot(const Vector& other) const {
        return x * other.x + y * other.y;
    }

    Vector normalized() const {
        double length = std::hypot(x, y);
        return Vector(x / length, y / length);
    }
};

double cross(const Vector& a, const Vector& b) {
    return a.x * b.y - b.x * a.y;
}

class Line {
public:
    Point p;
    Vector v;
    Line(const Point& p, const Vector& v) : p(p), v(v) {}

    std::optional<Point> intersect(const Line& other) const {
        double det = cross(v, other.v);
        if (std::fabs(det) < EPSILON)
            return std::nullopt;

        double t = cross(other.p - p, other.v) / det;
        return p + v * t;
    }
};

class LineSegment {
public:
    Point p1, p2;
    LineSegment(const Point& p1, const Point& p2) : p1(p1), p2(p2) {}

    Vector v() const {
        return Vector(p2.x - p1.x, p2.y - p1.y);
    }

    double distance(const Point& p) const {
        return Line(p1, v()).intersect(Line(p, v())).distance(p);
    }
};

bool approximately_equals(double a, double b) {
    return a == b || (std::fabs(a - b) <= std::max(std::fabs(a), std::fabs(b)) * 0.001);
}

bool approximately_same(const Point& a, const Point& b) {
    return approximately_equals(a.x, b.x) && approximately_equals(a.y, b.y);
}

std::vector<Point> normalize_contour(const std::vector<Point>& contour) {
    std::vector<Point> normalized_contour;
    for (size_t i = 0; i < contour.size(); ++i) {
        Point prev = contour[(i + contour.size() - 1) % contour.size()];
        Point point = contour[i];
        Point next = contour[(i + 1) % contour.size()];
        if (!(point == next || (point - prev).normalized() == (next - point).normalized())) {
            normalized_contour.push_back(point);
        }
    }
    return normalized_contour;
}

struct SplitEvent;
struct EdgeEvent;
struct OriginalEdge;
struct Subtree;
struct LAVertex;
struct SLAV;
struct LAV;
struct EventQueue;

using Event = std::variant<SplitEvent, EdgeEvent>;

struct SplitEvent {
    double distance;
    Point intersection_point;
    LAVertex* vertex;
    LineSegment opposite_edge;

    bool operator<(const SplitEvent& other) const {
        return distance > other.distance;
    }
};

struct EdgeEvent {
    double distance;
    Point intersection_point;
    LAVertex* vertex_a;
    LAVertex* vertex_b;

    bool operator<(const EdgeEvent& other) const {
        return distance > other.distance;
    }
};

struct OriginalEdge {
    LineSegment edge;
    Line bisector_left, bisector_right;
};

struct Subtree {
    Point source;
    double height;
    std::vector<Point> sinks;
};

struct LAVertex {
    Point point;
    LineSegment edge_left, edge_right;
    LAVertex* prev;
    LAVertex* next;
    LAV* lav;
    bool is_valid;
    bool is_reflex;
    Line bisector;

    LAVertex(const Point& point, const LineSegment& edge_left, const LineSegment& edge_right,
             const std::optional<std::pair<Vector, Vector>>& direction_vectors = std::nullopt)
        : point(point), edge_left(edge_left), edge_right(edge_right), prev(nullptr), next(nullptr),
          lav(nullptr), is_valid(true) {
        Vector creator_vectors[2] = { edge_left.v().normalized() * -1, edge_right.v().normalized() };
        auto vectors = direction_vectors ? *direction_vectors : std::make_pair(creator_vectors[0], creator_vectors[1]);
        is_reflex = cross(vectors.first, vectors.second) < 0;
        bisector = Line(point, (vectors.first + vectors.second) * (is_reflex ? -1 : 1));
    }

    std::optional<Event> next_event(const std::vector<OriginalEdge>& original_edges);

    void invalidate() {
        if (lav) {
            lav->invalidate(this);
        } else {
            is_valid = false;
        }
    }
};

struct LAV {
    LAVertex* head;
    SLAV* slav;
    int length;

    LAV(SLAV* slav) : head(nullptr), slav(slav), length(0) {}

    static LAV from_polygon(const std::vector<Point>& polygon, SLAV* slav) {
        LAV lav(slav);
        for (size_t i = 0; i < polygon.size(); ++i) {
            Point prev = polygon[(i + polygon.size() - 1) % polygon.size()];
            Point point = polygon[i];
            Point next = polygon[(i + 1) % polygon.size()];
            lav.length++;
            auto vertex = new LAVertex(point, LineSegment(prev, point), LineSegment(point, next));
            vertex->lav = &lav;
            if (lav.head == nullptr) {
                lav.head = vertex;
                vertex->prev = vertex->next = vertex;
            } else {
                vertex->next = lav.head;
                vertex->prev = lav.head->prev;
                vertex->prev->next = vertex;
                lav.head->prev = vertex;
            }
        }
        return lav;
    }

    static LAV from_chain(LAVertex* head, SLAV* slav) {
        LAV lav(slav);
        lav.head = head;
        for (auto* vertex = head; ; vertex = vertex->next) {
            lav.length++;
            vertex->lav = &lav;
            if (vertex->next == head) break;
        }
        return lav;
    }

    void invalidate(LAVertex* vertex) {
        assert(vertex->lav == this);
        vertex->is_valid = false;
        if (head == vertex) {
            head = head->next;
        }
        vertex->lav = nullptr;
    }

    LAVertex* unify(LAVertex* vertex_a, LAVertex* vertex_b, const Point& point) {
        auto replacement = new LAVertex(point, vertex_a->edge_left, vertex_b->edge_right,
                                        std::make_optional(std::make_pair(vertex_b->bisector.v, vertex_a->bisector.v)));
        replacement->lav = this;
        if (head == vertex_a || head == vertex_b) {
            head = replacement;
        }
        vertex_a->prev->next = replacement;
        vertex_b->next->prev = replacement;
        replacement->prev = vertex_a->prev;
        replacement->next = vertex_b->next;
        vertex_a->invalidate();
        vertex_b->invalidate();
        length--;
        return replacement;
    }
};

struct SLAV {
    std::vector<LAV> lavs;
    std::vector<OriginalEdge> original_edges;

    SLAV(const std::vector<Point>& polygon, const std::vector<std::vector<Point>>& holes) {
        lavs.push_back(LAV::from_polygon(normalize_contour(polygon), this));
        for (const auto& hole : holes) {
            lavs.push_back(LAV::from_polygon(normalize_contour(hole), this));
        }
        for (auto& lav : lavs) {
            for (auto* vertex = lav.head; ; vertex = vertex->next) {
                original_edges.push_back({ LineSegment(vertex->prev->point, vertex->point),
                                           Line(vertex->prev->point, vertex->prev->bisector.v),
                                           Line(vertex->point, vertex->bisector.v) });
                if (vertex->next == lav.head) break;
            }
        }
    }

    bool empty() const {
        return lavs.empty();
    }

    std::pair<Subtree, std::vector<Event>> handle_edge_event(const EdgeEvent& event) {
        std::vector<Point> sinks;
sinks.push_back(event.intersection_point);

        std::vector<Event> new_events;

        auto& lav = *event.vertex_a->lav;
        auto* unified_vertex = lav.unify(event.vertex_a, event.vertex_b, event.intersection_point);

        auto handle_split_event = [&](LAVertex* vertex) {
            auto split_event = vertex->next_event(original_edges);
            if (split_event) {
                new_events.push_back(*split_event);
            }
        };

        handle_split_event(unified_vertex->prev);
        handle_split_event(unified_vertex);
        handle_split_event(unified_vertex->next);

        Subtree subtree{event.intersection_point, event.distance, sinks};

        return {subtree, new_events};
    }

    std::pair<Subtree, std::vector<Event>> handle_split_event(const SplitEvent& event) {
        std::vector<Point> sinks;
        sinks.push_back(event.intersection_point);

        std::vector<Event> new_events;

        auto& lav = *event.vertex->lav;
        auto* vertex = event.vertex;

        auto* chain_head = new LAVertex(event.intersection_point, vertex->edge_left,
                                        LineSegment(event.intersection_point, vertex->next->point),
                                        std::make_optional(std::make_pair(event.vertex->bisector.v, vertex->next->bisector.v)));
        auto* chain_tail = new LAVertex(event.intersection_point, LineSegment(event.intersection_point, vertex->prev->point),
                                        vertex->edge_right, std::make_optional(std::make_pair(vertex->prev->bisector.v, vertex->bisector.v)));

        chain_head->next = chain_tail;
        chain_tail->prev = chain_head;
        chain_head->prev = vertex->prev;
        chain_tail->next = vertex->next;
        vertex->prev->next = chain_head;
        vertex->next->prev = chain_tail;

        lav.invalidate(vertex);

        auto new_lav_1 = LAV::from_chain(chain_head, this);
        auto new_lav_2 = LAV::from_chain(chain_tail, this);

        lavs.push_back(new_lav_1);
        lavs.push_back(new_lav_2);

        auto handle_split_event = [&](LAVertex* vertex) {
            auto split_event = vertex->next_event(original_edges);
            if (split_event) {
                new_events.push_back(*split_event);
            }
        };

        handle_split_event(chain_head->prev);
        handle_split_event(chain_head);
        handle_split_event(chain_tail);
        handle_split_event(chain_tail->next);

        Subtree subtree{event.intersection_point, event.distance, sinks};

        return {subtree, new_events};
    }

    std::vector<Subtree> process_next_events(std::priority_queue<Event>& event_queue) {
        std::vector<Subtree> new_subtrees;

        while (!event_queue.empty()) {
            auto event = event_queue.top();
            event_queue.pop();

            if (std::holds_alternative<EdgeEvent>(event)) {
                auto edge_event = std::get<EdgeEvent>(event);
                if (!edge_event.vertex_a->is_valid || !edge_event.vertex_b->is_valid) continue;

                auto [subtree, new_events] = handle_edge_event(edge_event);
                new_subtrees.push_back(subtree);

                for (const auto& new_event : new_events) {
                    event_queue.push(new_event);
                }
            } else if (std::holds_alternative<SplitEvent>(event)) {
                auto split_event = std::get<SplitEvent>(event);
                if (!split_event.vertex->is_valid) continue;

                auto [subtree, new_events] = handle_split_event(split_event);
                new_subtrees.push_back(subtree);

                for (const auto& new_event : new_events) {
                    event_queue.push(new_event);
                }
            }
        }

        return new_subtrees;
    }
};

std::optional<Event> LAVertex::next_event(const std::vector<OriginalEdge>& original_edges) {
    std::optional<EdgeEvent> edge_event;
    std::optional<SplitEvent> split_event;

    for (const auto& original_edge : original_edges) {
        auto intersection = bisector.intersect(original_edge.bisector_left);
        if (intersection && *intersection != point) {
            double dist = point.distance(*intersection);
            if (!split_event || dist < split_event->distance) {
                split_event = SplitEvent{dist, *intersection, this, original_edge.edge};
            }
        }

        intersection = bisector.intersect(original_edge.bisector_right);
        if (intersection && *intersection != point) {
            double dist = point.distance(*intersection);
            if (!split_event || dist < split_event->distance) {
                split_event = SplitEvent{dist, *intersection, this, original_edge.edge};
            }
        }
    }

    auto next_intersection = bisector.intersect(next->bisector);
    if (next_intersection) {
        double dist = point.distance(*next_intersection);
        if (next_intersection != point && next_intersection != next->point) {
            edge_event = EdgeEvent{dist, *next_intersection, this, next};
        }
    }

    if (split_event && (!edge_event || split_event->distance < edge_event->distance)) {
        return split_event;
    }
    if (edge_event) {
        return edge_event;
    }
    return std::nullopt;
}

std::vector<Subtree> compute_voronoi_subtrees(const std::vector<Point>& polygon, const std::vector<std::vector<Point>>& holes) {
    SLAV slav(polygon, holes);
    std::priority_queue<Event> event_queue;

    for (auto& lav : slav.lavs) {
        for (auto* vertex = lav.head; ; vertex = vertex->next) {
            auto next_event = vertex->next_event(slav.original_edges);
            if (next_event) {
                event_queue.push(*next_event);
            }
            if (vertex->next == lav.head) break;
        }
    }

    std::vector<Subtree> subtrees;

    while (!slav.empty()) {
        auto new_subtrees = slav.process_next_events(event_queue);
        subtrees.insert(subtrees.end(), new_subtrees.begin(), new_subtrees.end());
        slav.lavs.erase(std::remove_if(slav.lavs.begin(), slav.lavs.end(),
                                       [](const LAV& lav) { return lav.length == 0; }), slav.lavs.end());
    }

    return subtrees;
}

int main() {
    std::vector<Point> polygon = { Point(0, 0), Point(4, 0), Point(4, 4), Point(0, 4) };
    std::vector<std::vector<Point>> holes = { { Point(1, 1), Point(2, 1), Point(2, 2), Point(1, 2) } };

    auto subtrees = compute_voronoi_subtrees(polygon, holes);

    for (const auto& subtree : subtrees) {
        std::cout << "Subtree rooted at (" << subtree.source.x << ", " << subtree.source.y << ") with height "
                  << subtree.height << " and sinks:";
        for (const auto& sink : subtree.sinks) {
            std::cout << " (" << sink.x << ", " << sink.y << ")";
        }
        std::cout << std::endl;
    }

    return 0;
}
