#include "includes/functions.h"

// For the new file
#include <iostream>
#include <fstream>

// For timing
#include <chrono>

using namespace std;

int main(int argc, char **argv)
{

    visible flag_visible;
    main_alg flag_main_algorithm;
    algorithm flag_algorithm;
    anneal flag_anneal;
    bool flag_min;
    double long area_ch;
    double threshold;
    int initialization;
    int L_loc, L_sim;

    // Reading the points from the file and command line argumemt checks
    vector<Point_2> vect;
    char *path;
    char *outputFile;
    bool arg_i = false, arg_o = false, arg_proc = false;
    if (argc < 5 || argc == 6 || argc > 12)
    {
        cout << "Wrong arguments!" << endl;
        return -1;
    }
    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-i"))
        {
            if (argv[i + 1] == nullptr)
            {
                cout << "No point set input file given!" << endl;
                return -1;
            }
            path = argv[i + 1];
            arg_i = true;
        }
        else if (!strcmp(argv[i], "-o"))
        {
            if (argv[i + 1] == nullptr)
            {
                cout << "No output file given!" << endl;
                return -1;
            }
            outputFile = argv[i + 1];
            arg_o = true;
        }
        else if (!strcmp(argv[i], "-preprocess"))
        {
            arg_proc = true;
        }
    }

    if (!arg_i || !arg_o)
    {
        cout << "Wrong arguments !" << endl;
        return -1;
    }

    if (!arg_proc)
    {
        initialization = 1;
        flag_visible = random_vis;
    }

    Polygon_2 polygonChain;

    results temp_res;
    double temp_score;
    results res[NUM_OF_ALG][NUM_OF_FILES][18]; // 18 from all the possible sizes
    int fil_count[NUM_OF_ALG][18];
    for (int i = 0; i < NUM_OF_ALG; i++)
        for (int j = 0; j < NUM_OF_FILES; j++)
            fil_count[i][j] = 0;

    // Create file
    ofstream outfile(outputFile);

    int count = 0;
    char cur[4095]; // xrhsimopoieitai gia na swzoume to path
    memset(cur, 0, 4095);
    struct dirent **fil;
    int entr;
    entr = scandir(path, &fil, NULL, alphasort); // arithmos twn entries
    // scandir thelw efoson readdir diavazei arxeia me entelws tyxaia seira
    struct stat s;
    if (entr)
    { // an einai ta entries perissotera apo 0
        while (count < entr)
        {
            if (strcmp(fil[count]->d_name, ".") && strcmp(fil[count]->d_name, ".."))
            {
                // anaparistoun to paron kai to proghoumeno dir, ara den theloume na metaferoume pros apofygh lathwn
                strcpy(cur, path);
                strcat(cur, "/");
                strcat(cur, fil[count]->d_name);
                if (stat(cur, &s) == 0)
                {
                    cout << "path is " << cur << endl;
                    vect = readPoints(cur);
                    area_ch = readArea(cur);

                    temp_res.size = vect.size();

                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    for (int alg = 1; alg <= NUM_OF_ALG; alg++)
                    {
                        switch (alg)
                        {
                        case 1:
                            flag_algorithm = incremental;
                            flag_main_algorithm = local_search;
                            break;
                        case 2:
                            flag_algorithm = incremental;
                            flag_main_algorithm = simulated_annealing;
                            flag_anneal = local;
                            break;
                        case 3:
                            flag_algorithm = convex_hull; // change to incremental after 1000
                            flag_main_algorithm = simulated_annealing;
                            flag_anneal = local;
                            break;
                        }

                        temp_res.algorithm = alg;

                        for (int i = 0; i < 2; i++)
                        {
                            if (i == 0)
                                flag_min = true;
                            else
                                flag_min = false;
                            // Clock starts now
                            auto start = chrono::high_resolution_clock::now();

                            // cout << "before 1st algorithm" << endl;
                            if (flag_algorithm == incremental)
                            {
                                cout << "incremental ";
                                polygonChain = incremental_alg(vect, initialization, flag_visible);
                            }
                            else
                            { // Convex hull
                                cout << "convexhull ";
                                if (alg == 4)
                                {
                                    if (temp_res.size == 1000)
                                        flag_algorithm = incremental; // Change to incremental

                                    if (flag_min)
                                        flag_visible = min_vis;
                                    else
                                        flag_visible = max_vis;
                                }
                                polygonChain = convex_hull_alg(vect, flag_visible);
                            }

                            Polygon_2 optimizedPolygon;

                            if (flag_main_algorithm == local_search)
                            {
                                cout << "local search " << endl;
                                threshold = local_search_threshold_init(temp_res.size);
                                if (temp_res.size < 100)
                                {
                                    L_loc = 5;
                                }
                                else if (temp_res.size == 200)
                                {
                                    L_loc = 2;
                                }
                                else if (temp_res.size > 400)
                                {
                                    L_loc = 1;
                                }
                                optimizedPolygon = localSearch(polygonChain, L_loc, threshold, flag_min);
                            }
                            else
                            { // Simulated annealing
                                cout << "simulated annealing " << endl;
                                L_sim = simulated_annealing_L_sim_init(temp_res.size);
                                optimizedPolygon = simulatedAnnealing(vect, polygonChain, area_ch, L_sim, flag_min, flag_anneal, flag_algorithm);
                            }

                            // Clock stops
                            auto stop = chrono::high_resolution_clock::now();
                            auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);

                            if (duration.count() > 500 * vect.size())
                            {
                                cout << "Cut-off failed!" << endl;
                                if (flag_min)
                                    temp_score = 1;
                                else
                                    temp_score = 0;
                            }
                            else
                            {
                                cout << "Cut-off successful!" << endl;
                                temp_score = ((double)abs(optimizedPolygon.area())) / ((double)area_ch);
                            }

                            cout << "Time should be " << 500 * vect.size() << " and it is " << duration.count() << endl;

                            if (flag_min)
                                temp_res.min_score = temp_score;
                            else
                                temp_res.max_score = temp_score;
                        }

                        if (fil_count[alg - 1][size_to_index(temp_res.size)] >= NUM_OF_FILES)
                        {
                            cout << "Given more files for same input size than expected" << endl;
                            continue;
                        }
                        res[alg - 1][fil_count[alg - 1][size_to_index(temp_res.size)]][size_to_index(temp_res.size)] = temp_res;
                        fil_count[alg - 1][size_to_index(temp_res.size)]++;
                    }
                }

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
            count++;
        }
    }

    for (int i = 0; i < entr; i++)
    { // apodesmeysh mnhmhs efoson to scandir kanei kai malloc
        free(fil[i]);
    }
    free(fil);

    for (int alg = 1; alg <= NUM_OF_ALG; alg++)
    {
        outfile << " \t \t ||  \t \t \t Algorithm " << alg << " \t \t \t ||";
    }
    outfile << endl;
    outfile << "Size ";
    for (int i = 0; i < NUM_OF_ALG; i++)
        outfile << "\t || \t min_score \t || \t max_score \t || \t min_bound \t || \t max_bound \t ||";
    outfile << endl;

    double score_min, score_max, bound_min, bound_max;
    for (int siz = 0; siz < 18; siz++)
    {
        outfile << res[0][0][siz].size << " \t ||";
        for (int alg = 0; alg < NUM_OF_ALG; alg++)
        {
            score_min = 0;
            score_max = 0;
            for (int fil = 0; fil < NUM_OF_FILES; fil++)
            {
                if (fil == 0)
                {
                    bound_min = res[alg][fil][siz].min_score; // initialization
                    bound_max = res[alg][fil][siz].max_score;
                }
                else
                {
                    if (bound_min < res[alg][fil][siz].min_score)
                        bound_min = res[alg][fil][siz].min_score; // We are finding the maximum of the mins
                    if (bound_max > res[alg][fil][siz].max_score)
                        bound_max = res[alg][fil][siz].max_score; // We are finding the minimum of the maxs
                }
                score_min += res[alg][fil][siz].min_score;
                score_max += res[alg][fil][siz].max_score;
            }
            outfile << "\t " << score_min << "\t || \t" << score_max << " \t || \t " << bound_min << "\t || \t" << bound_max << "\t || \t";
        }
        outfile << endl;
    }

    outfile.close();

    return 0;
}
