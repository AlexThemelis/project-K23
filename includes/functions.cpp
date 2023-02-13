#include "functions.h"

using namespace std;

// Functions for sorting by ascending/descending based on x/y
bool compareDescending_x(Point_2 p1, Point_2 p2)
{
    return (p1.x() > p2.x());
}

bool compareAscending_y(Point_2 p1, Point_2 p2)
{
    return (p1.y() < p2.y());
}

bool compareDescending_y(Point_2 p1, Point_2 p2)
{
    return (p1.y() > p2.y());
}

// Function for reading the points by the input file
vector<Point_2> readPoints(char *name)
{
    ifstream input(name);
    vector<Point_2> vect;
    char *x_str = NULL;
    char *y_str = NULL;
    int x;
    int y;
    char *token = NULL;
    int count = 0;
    char buffer[100];
    for (string line; getline(input, line);)
    {
        count++;
        strcpy(buffer, line.c_str()); // String to char * in order to be an argument to strtok
        if (count <= 2)
            continue;                  // First 2 lines do not contain pointss
        token = strtok(buffer, " \t"); // We get the points by using strtok searching tabs ('/t')
        x_str = strtok(NULL, " \t");
        x = atoi(x_str);
        y_str = strtok(NULL, " \t");
        y = atoi(y_str);
        Point_2 p(x, y);
        vect.push_back(p);
    }
    input.close();
    return vect;
}

// Function for reading the area from the input file
unsigned long long readArea(char *name) // We use unsigned long long because for very big input files, the area becomes much bigger from int etc
{
    ifstream input(name);
    char *token = NULL;
    int count = 0;
    char buffer[100];
    unsigned long long area_ch;
    for (string line; getline(input, line);)
    {
        count++;
        if (count != 2)
            continue;
        strcpy(buffer, line.c_str()); // String to char * in order to be an argument to strtok  // First 2 lines do not contain pointss
        token = strtok(buffer, "\"");
        strtok(NULL, "\"");
        strtok(NULL, " \"");
        strtok(NULL, "\"");
        strtok(NULL, "\"");
        strtok(NULL, "\"");
        area_ch = atoi(strtok(NULL, "\""));
    }
    input.close();
    return area_ch;
}

int size_to_index(int size){
    switch(size){
        case 10:
            return 0;
        case 20:
            return 1;
        case 30:
            return 2;
        case 40:
            return 3;
        case 50:
            return 4;
        case 60:
            return 5;
        case 70:
            return 6;
        case 80:
            return 7;
        case 90:
            return 8;
        case 100:
            return 9;
        case 200:
            return 10;
        case 400:
            return 11;
        case 800:
            return 12;
        case 1000:
            return 13;
        case 2000:
            return 14;
        case 5000:
            return 15;
        case 10000:
            return 16;
        case 100000:
            return 17;
    }
    return -1;
}

int local_search_threshold_init(int size){
    switch(size){
        case 10:
            return 10000;
        case 20:
            return 10000;
        case 30:
            return 10000;
        case 40:
            return 10000;
        case 50:
            return 10000;
        case 60:
            return 10000;
        case 70:
            return 10000;
        case 80:
            return 10000;
        case 90:
            return 10000;
        case 100:
            return 10000;
        case 200:
            return 100000;
        case 400:
            return 300000;
        case 800:
            return 500000;
        case 1000:
            return 700000;
        case 2000:
            return 900000;
        case 5000:
            return 2000000;
        case 10000:
            return 5000000;
        case 100000:
            return 100000000;
    }
    return -1;
}

int simulated_annealing_L_sim_init(int size){
    switch(size){
        case 10:
            return 5000;
        case 20:
            return 5000;
        case 30:
            return 5000;
        case 40:
            return 6000;
        case 50:
            return 7000;
        case 60:
            return 7000;
        case 70:
            return 8000;
        case 80:
            return 8000;
        case 90:
            return 9000;
        case 100:
            return 9000;
        case 200:
            return 10000;
        case 400:
            return 10000;
        case 800:
            return 12000;
        case 1000:
            return 20000;
        case 2000:
            return 25000;
        case 5000:
            return 30000;
        case 10000:
            return 50000;
        case 100000:
            return 100000;
    }
    return -1;
}    