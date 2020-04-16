#include "Eigen/Dense"
#include <vector>
#include <array>

namespace Triangulator {
    struct Vertex;
    struct Segment;

    struct Vertex {
        Vertex() { };
        Vertex(Eigen::Vector2f init_coordinates) { coordinates = init_coordinates; }
        Eigen::Vector2f coordinates; // x and y coordinates of the Vertex
        bool is_helper = false; // Used to mark whether a vertex is used in the super triangle
        bool orphaned = true; // Used to mark whether a vertex has triangle
        int maps_to = -1; // Used to map the vertices to the input vertices
    };

    struct Edge {
        Edge() { };
        Edge(Vertex &init_orig, Vertex &init_dest) {
            orig = &init_orig;
            dest = &init_dest;
        }
        Vertex *orig;
        Vertex *dest;
        
        bool is_helper = false;
    };

    struct Circle {
        Circle() { };
        Circle(Eigen::Vector2f init_center, float init_radius) {
            center = init_center;
            radius = init_radius;
        }

        Eigen::Vector2f center;
        float radius;
    };

    /* Stores vertices, connected triangles and circle of a single triangle */
    struct Triangle {
        std::array<Vertex *, 3> v{nullptr, nullptr, nullptr};
        std::array<Triangle *, 3> t{nullptr, nullptr, nullptr};
        Circle c;
        bool is_bad = false;
        bool marked_for_deletion = false;
    };

    class Delaunay2D {
    private:
        void SplitSeg(int s);
        void SplitTri(Triangle &t, Vertex &p);
        int GetFirstEncroached(float minseglen);
        Triangle *GetFirstBadTriangle(float alpha);
        bool IsTriangleLowQuality(const Triangle &t, float alpha);
        Circle CalculateCircle(const Triangle &t);
        bool inCircle(const Circle &c, const Vertex &v);
        void DelaunayAddPoint(int i);
        void MarkTriangleDeletion(Triangle &t, Triangle *prev_t);
        void DelaunayTriangulation();
    public:
        ~Delaunay2D();
        std::vector<Vertex *> vertices;
        std::vector<Triangle *> triangles;
        std::vector<Edge> segments;
        void ToFile(std::string oname);
        void RefineRupperts(float alpha);
    };
}