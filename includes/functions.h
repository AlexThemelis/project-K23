#include <iostream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <vector>
#include <list>
#include <iterator>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Object.h>

#include <dirent.h>
#include <queue>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Point_2<K> Point_2;
typedef CGAL::Segment_2<K> Segment_2;
typedef std::vector<Point_2>::iterator veciterator;
typedef std::istream_iterator<Point_2> point2_iterator;

typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Kd_tree<Traits> Tree_2;
typedef CGAL::Fuzzy_iso_box<Traits> Fuzzy_iso_box;

#define NUM_OF_FILES 3
#define NUM_OF_ALG 3

enum algorithm { incremental, convex_hull, onion};
enum main_alg {local_search, simulated_annealing, ant_colony};
enum visible { random_vis, min_vis, max_vis};
enum anneal { local, global, subdivision};

Polygon_2 convex_hull_alg(std::vector<Point_2> , visible);
Polygon_2 incremental_alg(std::vector<Point_2>, int, visible);
Polygon_2 incremental_alg_SUBDIVISIAL(Polygon_2,std::vector<Point_2>, int, visible);
Polygon_2 localSearch(Polygon_2, int, double, bool);
Polygon_2 simulatedAnnealing(std::vector<Point_2>, Polygon_2, double long, int, bool, anneal, algorithm);

struct minTotal
{
    double min_distance;
    Segment_2 min_ch_edge;
    Point_2 min_point;
    double min_area;
};

struct localSearchList{
    Segment_2 edge;
    std::vector<Point_2> points;
    double area;
};

struct subdivision_pointsOfInterest{
    Point_2 q_left;
    Point_2 q_right;
    Point_2 r;
    Point_2 p;
};


struct results{
    int size;
    int algorithm; // 1 for 1st alg, 2 for 2nd etc
    double min_score;
    double max_score;
};

int size_to_index(int);
int local_search_threshold_init(int);
int simulated_annealing_L_sim_init(int);


std::vector<Point_2> readPoints(char *);
unsigned long long readArea(char *);

bool compareDescending_x(Point_2, Point_2);

bool compareAscending_y(Point_2 p1, Point_2 p2);

bool compareDescending_y(Point_2 p1, Point_2 p2);
