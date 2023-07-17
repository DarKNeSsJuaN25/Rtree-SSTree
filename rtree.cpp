#include <SDL2/SDL.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <chrono>

const int M = 50;
const int m = 2; // m <= M/2
const double EPS = 1e-9;
const int NUM_DIMENSIONS = 10;

enum SplitType {
    LINEAR_SPLIT,
    QUADRATIC_SPLIT,
    BROWNIE_SPLIT
};

struct Point {
    double coords[NUM_DIMENSIONS];
    Point() {}
    Point(double x, double y) {
        coords[0] = x;
        coords[1] = y;
    }
};

struct MBB {
    Point lower, upper;
    MBB(){
        lower = Point(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        upper = Point(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    }
    MBB(const Point& p1, const Point& p2) {
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            lower.coords[d] = std::min(p1.coords[d], p2.coords[d]);
            upper.coords[d] = std::max(p1.coords[d], p2.coords[d]);
        }
    }

    double expansion_needed(const Point& p) const {
        double expansion = 0;
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            expansion += std::max(0.0, p.coords[d] - upper.coords[d]);
            expansion += std::max(0.0, lower.coords[d] - p.coords[d]);
        }
        return expansion;
    }

    void expand(const Point& p) {
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            lower.coords[d] = std::min(lower.coords[d], p.coords[d]);
            upper.coords[d] = std::max(upper.coords[d], p.coords[d]);
        }
    }

    double area() const {
        double area = 1;
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            area *= upper.coords[d] - lower.coords[d];
        }
        return area;
    }
};

static double point_distance_squared(const Point& p1, const Point& p2) {
    double distance_squared = 0;
    for (int d = 0; d < NUM_DIMENSIONS; d++) {
        double diff = p1.coords[d] - p2.coords[d];
        distance_squared += diff * diff;
    }
    return distance_squared;
}

static double point_distance(const Point& p1, const Point& p2) {
    return std::sqrt(point_distance_squared(p1, p2));
}

class RTreeNode {
public:
    RTreeNode(bool is_leaf = true, std::weak_ptr<RTreeNode> parent = std::weak_ptr<RTreeNode>());

    bool is_leaf() const { return is_leaf_; }
    const MBB& mbb() const { return mbb_; }
    const std::vector<Point>& points() const { return points_; }
    const std::vector<std::shared_ptr<RTreeNode>>& children() const { return children_; }
    std::weak_ptr<RTreeNode> parent() const { return parent_; }

    void insert(const Point& p);
    std::vector<Point> best_first_search(const Point& target) const;

    static SplitType split_type;

private:
    bool is_leaf_;
    MBB mbb_;
    std::vector<Point> points_;
    std::vector<std::shared_ptr<RTreeNode>> children_;
    std::weak_ptr<RTreeNode> parent_;

    void split();
    void linear_split();
    void quadratic_split();
    void brownie_split();
    std::shared_ptr<RTreeNode> choose_subtree(const Point& p) const;
    double compute_overlap(const MBB& mbb1, const MBB& mbb2) const;
};

RTreeNode::RTreeNode(bool is_leaf, std::weak_ptr<RTreeNode> parent) : is_leaf_(is_leaf), parent_(parent) {
    if (is_leaf) {
        mbb_ = MBB();
    }
}

void RTreeNode::insert(const Point& p) {
    if (is_leaf_) {
        points_.push_back(p);
        mbb_.expand(p);

        // Check if the node overflows
        if (points_.size() > M) {
            split();
        }
    } else {
        // If the current node is not a leaf, choose the best subtree to insert the point
        std::shared_ptr<RTreeNode> subtree = choose_subtree(p);
        subtree->insert(p);
        mbb_.expand(p);
    }
}

std::shared_ptr<RTreeNode> RTreeNode::choose_subtree(const Point& p) const {
    std::shared_ptr<RTreeNode> best_child;
    double min_expansion = std::numeric_limits<double>::max();

    for (const auto& child : children_) {
        double child_expansion = child->mbb().expansion_needed(p);
        if (child_expansion < min_expansion) {
            min_expansion = child_expansion;
            best_child = child;
        } else if (child_expansion == min_expansion) {
            // Choose the child whose MBB has the smallest area
            double child_area = child->mbb().area();
            double best_child_area = best_child->mbb().area();
            if (child_area < best_child_area) {
                best_child = child;
            }
        }
    }

    return best_child;
}

void RTreeNode::split() {
    switch (split_type) {
        case LINEAR_SPLIT:
            linear_split();
            break;
        case QUADRATIC_SPLIT:
            quadratic_split();
            break;
        case BROWNIE_SPLIT:
            brownie_split();
            break;
        default:
            throw std::runtime_error("Unknown split type");
    }
}

SplitType RTreeNode::split_type;
// Find the dimension with the maximum spread
int find_dividing_dimension(const std::vector<Point>& points) {
    double max_spread = -1;
    int dividing_dimension = -1;
    for (int d = 0; d < NUM_DIMENSIONS; d++) {
        double min_coord = std::numeric_limits<double>::max();
        double max_coord = std::numeric_limits<double>::lowest();
        for (const auto& p : points) {
            min_coord = std::min(min_coord, p.coords[d]);
            max_coord = std::max(max_coord, p.coords[d]);
        }
        double spread = max_coord - min_coord;
        if (spread > max_spread) {
            max_spread = spread;
            dividing_dimension = d;
        }
    }
    return dividing_dimension;
}

// Calculate the centroid of a set of points in a given dimension
double calculate_centroid(const std::vector<Point>& points, int dimension) {
    double sum = 0;
    for (const auto& p : points) {
        sum += p.coords[dimension];
    }
    return sum / points.size();
}

void RTreeNode::linear_split() {
    // Find the best dimension to split the points
    int dividing_dimension = find_dividing_dimension(points_);

    // Sort the points along the dividing dimension
    std::sort(points_.begin(), points_.end(), [&](const Point& p1, const Point& p2) {
        return p1.coords[dividing_dimension] < p2.coords[dividing_dimension];
    });

    // Split the points into two groups
    int num_groups = 2;
    std::vector<std::vector<Point>> point_groups(num_groups);
    int group_size = std::ceil(points_.size() * 1.0 / num_groups);
    for (int i = 0; i < num_groups; i++) {
        int start_idx = i * group_size;
        int end_idx = std::min((i + 1) * group_size, static_cast<int>(points_.size()));
        for (int j = start_idx; j < end_idx; j++) {
            point_groups[i].push_back(points_[j]);
        }
    }

    // Create two new child nodes and add the points to them
    for (int i = 0; i < num_groups; i++) {
        auto child = std::make_shared<RTreeNode>(true);
        for (const auto& p : point_groups[i]) {
            child->insert(p);
        }
        child->mbb_ = MBB(point_groups[i][0], point_groups[i][point_groups[i].size() - 1]);
        for (int j = 0; j < point_groups[i].size(); j++) {
            child->mbb_.expand(point_groups[i][j]);
        }
        children_.push_back(child);
    }

    // Clear the points vector and set is_leaf_ to false
    points_.clear();
    is_leaf_ = false;
}
void RTreeNode::quadratic_split(){}
void RTreeNode::brownie_split(){}

std::vector<Point> RTreeNode::best_first_search(const Point& target) const {
    struct NodeDistance {
        const RTreeNode* node;
        double distance;

        bool operator<(const NodeDistance& other) const {
            // Compare based on distance in descending order
            return distance < other.distance;
        }
    };

    std::priority_queue<NodeDistance> pq;
    std::vector<Point> result;

    pq.push({this, point_distance(target, mbb_.lower)});

    while (!pq.empty()) {
        const RTreeNode* current_node = pq.top().node;
        pq.pop();

        if (current_node->is_leaf()) {
            for (const auto& p : current_node->points()) {
                result.push_back(p);
            }
        } else {
            for (const auto& child : current_node->children()) {
                pq.push({child.get(), point_distance(target, child->mbb().lower)});
            }
        }
    }

    return result;
}

class RTree {
public:
    RTree(SplitType split_type = LINEAR_SPLIT);
    const RTreeNode* get_root() const { return root_.get(); }
    void draw(const RTreeNode* node, SDL_Renderer* renderer) const;
    void print_ascii() const;
    void insert(const Point& p);
    std::vector<Point> best_first_search(const Point& target) const;

private:
    std::shared_ptr<RTreeNode> root_;
    void print_ascii_node(const RTreeNode* node, int depth = 0) const;
};

RTree::RTree(SplitType split_type) {
    RTreeNode::split_type = split_type;
    root_ = std::make_shared<RTreeNode>(true);
}

void RTree::insert(const Point& p) {
    root_->insert(p);
}

std::vector<Point> RTree::best_first_search(const Point& target) const {
    if (root_) {
        return root_->best_first_search(target);
    } else {
        return {};
    }
}

void RTree::draw(const RTreeNode* node, SDL_Renderer* renderer) const {
    if (node->is_leaf()) {
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    } else {
        SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
    }

    SDL_Rect rect;
    rect.x = static_cast<int>(node->mbb().lower.coords[0]);
    rect.y = static_cast<int>(node->mbb().lower.coords[1]);
    rect.w = static_cast<int>(node->mbb().upper.coords[0] - node->mbb().lower.coords[0]);
    rect.h = static_cast<int>(node->mbb().upper.coords[1] - node->mbb().lower.coords[1]);

    SDL_RenderDrawRect(renderer, &rect);

    // Draw the points in the leaf nodes
    if (node->is_leaf()) {
        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);

        for (const Point& p : node->points()) {
            // Draw a small rectangle (2x2 pixels) for each point
            SDL_Rect point_rect;
            point_rect.x = static_cast<int>(p.coords[0]) - 1;
            point_rect.y = static_cast<int>(p.coords[1]) - 1;
            point_rect.w = 2;
            point_rect.h = 2;
            SDL_RenderFillRect(renderer, &point_rect);
        }
    }
    for (const auto& child : node->children()) {
        draw(child.get(), renderer);
    }
}

void RTree::print_ascii() const {
    print_ascii_node(root_.get());
}

void RTree::print_ascii_node(const RTreeNode* node, int depth) const {
    if (node == nullptr) {
        return;
    }

    std::string indentation(depth * 2, ' ');

    std::cout << indentation << "MBB: ("
              << node->mbb().lower.coords[0] << ", " << node->mbb().lower.coords[1] << "), ("
              << node->mbb().upper.coords[0] << ", " << node->mbb().upper.coords[1] << ")\n";

    if (node->is_leaf()) {
        std::cout << indentation << "Points:\n";
        for (const Point& p : node->points()) {
            std::cout << indentation << "  (" << p.coords[0] << ", " << p.coords[1] << ")\n";
        }
    } else {
        std::cout << indentation << "Children:\n";
        for (const auto& child : node->children()) {
            print_ascii_node(child.get(), depth + 1);
        }
    }
}

// Genera un punto aleatorio en el rango [0, 100] en cada dimensión
Point generateRandomPoint() {
    Point p;
    for (int i = 0; i < NUM_DIMENSIONS; i++) {
        p.coords[i] = std::rand() % 101;  // Números aleatorios entre 0 y 100
    }
    return p;
}

int main() {
    std::srand(std::time(nullptr));

    RTree rtree(LINEAR_SPLIT);
    std::vector<Point> dataset;

    // Genera un millón de datos aleatorios y los inserta en el RTree
    for (int i = 0; i < 1000000; i++) {
        std::cout << i << std::endl;
        Point p = generateRandomPoint();
        rtree.insert(p);
        dataset.push_back(p);
    }

    int num_dimensions = NUM_DIMENSIONS;

    std::cout << "Tiempo de búsqueda vs número de dimensiones\n";
    std::cout << "5 vecinos más cercanos\n";
    std::cout << "Un millón de datos aleatorios\n";
    std::cout << "M = " << M << "\n";

    // Realiza búsquedas de los 5 vecinos más cercanos y mide el tiempo
    for (int dim = 1; dim <= num_dimensions; dim++) {
        auto start_time = std::chrono::high_resolution_clock::now();

        Point target = generateRandomPoint();

        std::vector<Point> neighbors = rtree.best_first_search(target);

        std::sort(neighbors.begin(), neighbors.end(), [&](const Point& p1, const Point& p2) {
            return point_distance(target, p1) < point_distance(target, p2);
        });

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "Dimensiones: " << dim << " - Tiempo: " << duration.count() << " ms\n";
    }

    return 0;
}
