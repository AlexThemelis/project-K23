#include "functions.h"

using namespace std;

Polygon_2 simulatedAnnealing(vector<Point_2> vect, Polygon_2 polygonChain, double long convex_hull_area, int L, bool minimazation, anneal ann, algorithm alg){
  
    int n_vertices = polygonChain.size();

    double long initial_area = abs(polygonChain.area());
    vector<localSearchList> list;

    srand(time(NULL));
    double random_threshold =((double) rand() / (RAND_MAX));
    double long old_energy;
    double long new_energy;
    double long DE;

    if(minimazation){
        old_energy = n_vertices*(initial_area/convex_hull_area);
    }else{
        old_energy = n_vertices*(1 - initial_area/convex_hull_area);
    }

    Polygon_2 polygonChainCopy = polygonChain;
    Point_2 p,q,r,s;
    int count = 1;
    int bool_count;
    int poly_size = polygonChain.size();
    bool done = false;
    double T = 1;
    bool unchanged = false;
    Tree_2 kdTree;
    bool bool_intersect = false;
    Segment_2 segm_int;
    Point_2 point_int;
    std::list<Point_2> result;
    double minx,miny,maxx,maxy;
    int k,m;
    Polygon_2 tempCopy,tempCopy2;
    Point_2 tempMinPoint;
    Point_2 sourceMinPoint;
    Polygon_2 mergedPolygonChain;
    if(ann == subdivision){
        m = 35;
        k = ceil((n_vertices-1)/(m-1));

        sort(vect.begin(), vect.end()); // Ascending by x, Actually not needed because vect is already sorted from input but just to be sure

        vector<Polygon_2> chains;
        Polygon_2 tempPol;
        chains.reserve(k); // make at least k size
        auto poi = vect.begin();
        auto prev_poi = vect.begin();
        int count = 1;
        vector<Point_2> tempPoints;
        vector<subdivision_pointsOfInterest> Points_Int;
        subdivision_pointsOfInterest tempPointsInt;


        for(int i = 0; i < k; i++){
            while(count<=m && poi != vect.end()){
                if(count == 2){
                    tempPointsInt.q_left = *prev_poi;
                    tempPointsInt.r = *poi;
                } 
                tempPol.push_back(*poi);
                if(count != m){
                    prev_poi = poi;
                    poi++;
                }
                else{ 
                    if((poi+1) != vect.end()){
                        if(!(prev_poi->y() < poi->y() && poi->y() > (poi+1)->y())){
                            prev_poi = poi;
                            poi++; //Increment poi
                            continue; //Dont increment count, The 2 conditions did not hold so we put another one inside
                        }
                    }
                    // In order to have point q, same in both spatial subdivisions
                }
                count++;
            }
            if((poi+1) != vect.end()){
                tempPointsInt.p = *prev_poi;
                tempPointsInt.q_right = *poi;
            }
            Points_Int.push_back(tempPointsInt);
            count = 1;
            chains.push_back(tempPol);
            tempPol.clear();
        }

        auto it = chains.rbegin();
        if(it->size() == 0) chains.pop_back(); // Sometimes a blank last was coming to play
        

        count = 0;
        for(auto it = chains.rbegin() ; it != chains.rend(); ++it){
            if(it->size() > m/4) break;

             // If here then last subdivision is too small, so we append it to the end of the 2nd last   
            if(count == 0){
                k--;
                count++; 
                for(auto i = it->begin(); i != it->end(); ++i){
                    tempPoints.push_back(*i);
                }
            }else if(count == 1){
                    for(auto i = tempPoints.begin(); i != tempPoints.end(); ++i){
                        it->push_back(*i);
                    }
                    count++;
                }else{
                    chains.pop_back();
                    Points_Int.pop_back();
                    break;
                }
        }
        tempPoints.clear();

        tempPoints.clear();
        for(auto it = chains.begin(); it != chains.end();++it){
            tempPoints.clear();
            for(auto i = it->begin(); i != it->end(); ++i){
                tempPoints.push_back(*i);
            }
            if(alg == incremental){
                *it = incremental_alg(tempPoints, 1, random_vis); // "1a" and "random"
            }else{
                *it = convex_hull_alg(tempPoints,random_vis);
            }
        }

        bool exists = false;
        auto iter_poi = Points_Int.begin();
        for(auto it = chains.begin(); it != chains.end(); ++it){
            for(auto i = it->edges_begin(); i != it->edges_end(); ++i){
                if((i->source() == iter_poi->q_left && i->target() == iter_poi->r) || (i->source() == iter_poi->r && i->target() == iter_poi->q_left)){
                    exists = true;
                    break;
                }
            }
            if(!exists){
                // Manually fix it!

                tempCopy = *it;
                tempPoints.clear();
                for(auto i1 = tempCopy.begin(); i1 != tempCopy.end(); ++i1){
                    if(*i1 == iter_poi->r || *i1 == iter_poi->q_left) continue;
                    tempPoints.push_back(*i1);
                }
                if(alg == incremental){
                    tempCopy2 = incremental_alg(tempPoints, 1, random_vis); // "1a" and "random"
                }else{
                    tempCopy2 = convex_hull_alg(tempPoints,random_vis);
                }
                tempCopy = tempCopy2;
                
                tempMinPoint = *min_element(tempCopy.begin(), tempCopy.end(), [](const auto& a, const auto& b) { return a.x() < b.x(); });

                vector<Point_2> vector_qr;
                vector_qr.push_back(iter_poi->r);
                vector_qr.push_back(iter_poi->q_left);

                *it = incremental_alg_SUBDIVISIAL(tempCopy,vector_qr,1,random_vis);

            }
            exists = false;
            iter_poi++;
        }

        // Now we apply global steps for every chain , We have to be carefull not to change one of the 4 importants
        iter_poi = Points_Int.begin();
        count = 0;
        Polygon_2 convexHull,chp;
        vector<Point_2> resultConv;
        for(auto iter_Chain = chains.begin(); iter_Chain != chains.end(); ++iter_Chain){
            initial_area = iter_Chain->area();
            n_vertices = iter_Chain->size();
            
            convexHull = *iter_Chain;

            // Calculating convex hull
            const Polygon_2::Vertices &range = convexHull.vertices();
            CGAL::convex_hull_2(range.begin(), range.end(), back_inserter(resultConv));

            for(auto it = resultConv.begin(); it!=resultConv.end();++it){
                chp.push_back(*it);
            }
            convex_hull_area = chp.area();
            resultConv.clear();

            if(minimazation){
                old_energy = n_vertices*(initial_area/convex_hull_area);
            }else{
                old_energy = n_vertices*(1 - initial_area/convex_hull_area);
            }
            polygonChainCopy = *iter_Chain;
            while(T>=0){
                unchanged = true;
                done = false;
                for(auto it = iter_Chain->edges_begin(); it != iter_Chain->edges_end(); it++){
                    for(auto it1 = iter_Chain->begin(); it1 != iter_Chain->end(); it1++){
                        if(count == poly_size) break; // if count == poly_size then the target will be the same like the 1st iteration
                        q = it->target();
                        s = *it1;
                        p = it->source();
                        r = (*(it+1)).target();
                        if(q == s || p == s || r == s) continue; // p and q shouldn't be the same or connected with an edge

                        // Make sure that we don't change the 4 important points 
                        if(q == iter_poi->q_left || q == iter_poi->q_right || q == iter_poi->p || q == iter_poi->r) continue;
                        if(s == iter_poi->q_left || s == iter_poi->q_right || s == iter_poi->p || s == iter_poi->r) continue;
                        if(p == iter_poi->q_left || p == iter_poi->q_right || p == iter_poi->p || p == iter_poi->r) continue;
                        if(r == iter_poi->q_left || r == iter_poi->q_right || r == iter_poi->p || r == iter_poi->r) continue;

                        

                        for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){ // remove q
                            if(*i == q){
                                polygonChainCopy.erase(i);
                                break;
                            }
                        }

                        for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){ // insert it after s
                            if(*i == s){
                                polygonChainCopy.insert(i+1,q);
                            }
                        }

                         T -= 1.0/((double) L); // No matter what happens we decrease T

                        if(polygonChainCopy.is_simple()){
                            if(minimazation){
                                new_energy = n_vertices*(abs(polygonChainCopy.area())/convex_hull_area);
                                DE = new_energy - old_energy;
                                if(DE < 0 || exp(-DE/T) > random_threshold){
                                    done = true;
                                    old_energy = new_energy;
                                }
                            }else{
                                new_energy = n_vertices*(1 - abs(polygonChainCopy.area())/convex_hull_area);
                                DE = new_energy - old_energy;
                                if(DE > 0 ||  exp(-DE/T) > random_threshold){
                                    done = true;
                                    old_energy = new_energy;
                                }
                            }
                                if(done){    
                                    break;
                                }
                        }

                        polygonChainCopy = *iter_Chain;

                    }
                    count++;
                    if(done){
                        done = false;
                        unchanged = false; // We changed something 1 time
                        break;
                    }
                }
                if(unchanged == true) break;
            }

            T = 1;
            count = 0;
            iter_poi++;
            *iter_Chain = polygonChainCopy;
        }

        iter_poi = Points_Int.begin();
        for(auto it = chains.begin(); it != chains.end(); ++it){
            iter_poi++;
        }

        iter_poi = Points_Int.begin();
        // Merge all the subdivisions into 1
        for(auto it = chains.begin(); it!= chains.end(); ++it){
            for(auto i = it->edges_begin(); i != it->edges_end(); ++i){
                if(i->source() == iter_poi->p && i->target() == iter_poi->q_right && count){ // Remove q in order to connect p with r
                    for(auto i1 = it->begin(); i1 != it->end();++i1){
                        if(*i1 == i->target()){
                            mergedPolygonChain.erase(i1);
                            break;
                        }
                    }
                    iter_poi++;
                }else{
                    mergedPolygonChain.push_back(i->source());
                }
            }
            iter_poi++;
        }

        polygonChain = mergedPolygonChain;
        polygonChainCopy = polygonChain;



        // Local search in the merged Polygon Chain!!
        T = 1;
        while(T>=0){
            for(auto it = polygonChain.edges_begin(); it != polygonChain.edges_end(); ++it){
                q = it->source();
                r = it->target();

                for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){
                    if(*i == q){
                        polygonChainCopy.erase(i);
                        break;
                    }
                }

                for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){
                    if(*i == r){
                        polygonChainCopy.insert(i+1,q);
                        break;
                    }
                }



                bool_count = 0;
                for(auto i = polygonChainCopy.edges_begin(); i != polygonChainCopy.edges_end();++i){
                    if(i->target() == r){
                        p = i->source();
                        bool_count++;
                    }
                    if(i->source() == q){
                        s = i->target();
                        bool_count++;
                    }
                    if(bool_count == 2) break;
                }



                Segment_2 pr(p,r);
                Segment_2 qs(q,s);
                
                CGAL::Object obj = intersection(pr,qs);
                
                // If pr, qs intersect
                if(CGAL::assign(segm_int,obj) || CGAL::assign(point_int,obj)){
                    polygonChainCopy = polygonChain;
                    continue;
                }

                for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){
                    kdTree.insert(*i);
                }

                kdTree.build();

                minx = p.x();
                miny = p.y();
                if(q.x()<minx) minx = q.x(); if(q.y()<miny) miny = q.y();
                if(r.x()<minx) minx = r.x(); if(r.y()<miny) miny = r.y();
                if(s.x()<minx) minx = s.x(); if(s.y()<miny) miny = s.y();

                maxx = p.x();
                maxy = p.y();
                if(q.x()>maxx) maxx = q.x(); if(q.y()>maxy) maxy = q.y();
                if(r.x()>maxx) maxx = r.x(); if(r.y()>maxy) maxy = r.y();
                if(s.x()>maxx) maxx = s.x(); if(s.y()>maxy) maxy = s.y();
                //

                Fuzzy_iso_box rect_box(Point_2(minx,miny),Point_2(maxx,maxy));

                kdTree.search(back_inserter(result), rect_box);

                // We remove duplicates from list
                // https://stackoverflow.com/questions/4877504/how-can-i-remove-duplicate-values-from-a-list-in-c
                set<Point_2> found;
                for (std::list<Point_2>::iterator x = result.begin(); x != result.end();) {
                    if (!found.insert(*x).second) {
                        x = result.erase(x);
                    }
                    else {
                        ++x;
                    }
                }
                //
                
                bool_intersect = false;
                for(auto i1 = result.begin(); i1 != result.end(); ++i1){
                    for(auto i = polygonChain.edges_begin(); i != polygonChain.edges_end(); i++){
                        if(i->source() == *i1 || i->target() == *i1){
                            Segment_2 test_int(i->source(),i->target());
                            CGAL::Object obj1 = intersection(pr,test_int);
                            CGAL::Object obj2 = intersection(qs,test_int);


                            if(CGAL::assign(segm_int,obj1) || CGAL::assign(segm_int,obj2)){
                                bool_intersect = true;
                                break;
                            }


                            if(const Point_2 * point = CGAL::object_cast<Point_2>(&obj1)){
                                if(*point != i->source() && *point != i->target()){
                                    bool_intersect = true;
                                    break;
                                }
                            }

                            if(const Point_2 * point = CGAL::object_cast<Point_2>(&obj2)){
                                if(*point != i->source() && *point != i->target()){
                                    bool_intersect = true;
                                    break;
                                }
                            }
                        }

                    }
                    if(bool_intersect) break;
                }

                if(!bool_intersect){
                    if(minimazation){
                        new_energy = n_vertices*(abs(polygonChainCopy.area())/convex_hull_area);
                        DE = new_energy - old_energy;
                        if(DE < 0 || exp(-DE/T) > random_threshold){
                            done = true;
                            old_energy = new_energy;
                        }
                    }else{
                        new_energy = n_vertices*(1 - abs(polygonChainCopy.area())/convex_hull_area);
                        DE = new_energy - old_energy;
                        if(DE > 0 ||  exp(-DE/T) > random_threshold){
                            done = true;
                            old_energy = new_energy;
                        }
                    }
                }

                T -= 1.0/((double) L);
                
                result.clear();
                kdTree.clear();
                kdTree.invalidate_build();

                if(done){
                    done = false;
                    polygonChain = polygonChainCopy;
                    break;
                }
                polygonChainCopy = polygonChain;
                
            }
        }
    }else{
        if(ann == local){ //Local Step
            while(T>=0){
                for(auto it = polygonChain.edges_begin(); it != polygonChain.edges_end(); ++it){
                    q = it->source();
                    r = it->target();

                    for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){
                        if(*i == q){
                            polygonChainCopy.erase(i);
                            break;
                        }
                    }

                    for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){
                        if(*i == r){
                            polygonChainCopy.insert(i+1,q);
                            break;
                        }
                    }



                    bool_count = 0;
                    for(auto i = polygonChainCopy.edges_begin(); i != polygonChainCopy.edges_end();++i){
                        if(i->target() == r){
                            p = i->source();
                            bool_count++;
                        }
                        if(i->source() == q){
                            s = i->target();
                            bool_count++;
                        }
                        if(bool_count == 2) break;
                    }



                    Segment_2 pr(p,r);
                    Segment_2 qs(q,s);
                    
                    CGAL::Object obj = intersection(pr,qs);
                    
                    // If pr, qs intersect
                    if(CGAL::assign(segm_int,obj) || CGAL::assign(point_int,obj)){
                        polygonChainCopy = polygonChain;
                        continue;
                    }

                    for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){
                        kdTree.insert(*i);
                    }

                    kdTree.build();

                    minx = p.x();
                    miny = p.y();
                    if(q.x()<minx) minx = q.x(); if(q.y()<miny) miny = q.y();
                    if(r.x()<minx) minx = r.x(); if(r.y()<miny) miny = r.y();
                    if(s.x()<minx) minx = s.x(); if(s.y()<miny) miny = s.y();

                    maxx = p.x();
                    maxy = p.y();
                    if(q.x()>maxx) maxx = q.x(); if(q.y()>maxy) maxy = q.y();
                    if(r.x()>maxx) maxx = r.x(); if(r.y()>maxy) maxy = r.y();
                    if(s.x()>maxx) maxx = s.x(); if(s.y()>maxy) maxy = s.y();
                    //

                    Fuzzy_iso_box rect_box(Point_2(minx,miny),Point_2(maxx,maxy));

                    kdTree.search(back_inserter(result), rect_box);

                    // We remove duplicates from list
                    // https://stackoverflow.com/questions/4877504/how-can-i-remove-duplicate-values-from-a-list-in-c
                    set<Point_2> found;
                    for (std::list<Point_2>::iterator x = result.begin(); x != result.end();) {
                        if (!found.insert(*x).second) {
                            x = result.erase(x);
                        }
                        else {
                            ++x;
                        }
                    }
                    //
                    
                    bool_intersect = false;
                    for(auto i1 = result.begin(); i1 != result.end(); ++i1){
                        for(auto i = polygonChain.edges_begin(); i != polygonChain.edges_end(); i++){
                            if(i->source() == *i1 || i->target() == *i1){
                                Segment_2 test_int(i->source(),i->target());
                                CGAL::Object obj1 = intersection(pr,test_int);
                                CGAL::Object obj2 = intersection(qs,test_int);


                                if(CGAL::assign(segm_int,obj1) || CGAL::assign(segm_int,obj2)){
                                    bool_intersect = true;
                                    break;
                                }


                                if(const Point_2 * point = CGAL::object_cast<Point_2>(&obj1)){
                                    if(*point != i->source() && *point != i->target()){
                                        bool_intersect = true;
                                        break;
                                    }
                                }

                                if(const Point_2 * point = CGAL::object_cast<Point_2>(&obj2)){
                                    if(*point != i->source() && *point != i->target()){
                                        bool_intersect = true;
                                        break;
                                    }
                                }
                            }

                        }
                        if(bool_intersect) break;
                    }

                    if(!bool_intersect){
                        if(minimazation){
                            new_energy = n_vertices*(abs(polygonChainCopy.area())/convex_hull_area);
                            DE = new_energy - old_energy;
                            if(DE < 0 || exp(-DE/T) > random_threshold){
                                done = true;
                                old_energy = new_energy;
                            }
                        }else{
                            new_energy = n_vertices*(1 - abs(polygonChainCopy.area())/convex_hull_area);
                            DE = new_energy - old_energy;
                            if(DE > 0 ||  exp(-DE/T) > random_threshold){
                                done = true;
                                old_energy = new_energy;
                            }
                        }
                    }

                    T -= 1.0/((double) L);
                    
                    result.clear();
                    kdTree.clear();
                    kdTree.invalidate_build();

                    if(done){
                        done = false;
                        polygonChain = polygonChainCopy;
                        break;
                    }
                    polygonChainCopy = polygonChain;
                    
                }
            }
        }else{ //Global Step
            while(T>=0){
                unchanged = true;
                for(auto it = polygonChain.edges_begin(); it != polygonChain.edges_end(); it++){
                    for(auto it1 = polygonChain.begin(); it1 != polygonChain.end(); it1++){
                        if(count == poly_size) break; // if count == poly_size then the target will be the same like the 1st iteration
                        q = it->target();
                        s = *it1;
                        p = it->source();
                        r = (*(it+1)).target();
                        if(q == s || p == s || r == s) continue; // p and q shouldn't be the same or connected with an edge

                        
                        for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){ // remove q
                            if(*i == q){
                                polygonChainCopy.erase(i);
                                break;
                            }
                        }

                        for(auto i = polygonChainCopy.begin(); i != polygonChainCopy.end();++i){ // insert it after s
                            if(*i == s){
                                polygonChainCopy.insert(i+1,q);
                            }
                        }

                        if(polygonChainCopy.is_simple()){
                            if(minimazation){
                                new_energy = n_vertices*(abs(polygonChainCopy.area())/convex_hull_area);
                                DE = new_energy - old_energy;
                                if(DE < 0 || exp(-DE/T) > random_threshold){
                                    done = true;
                                    old_energy = new_energy;
                                }
                            }else{
                                new_energy = n_vertices*(1 - abs(polygonChainCopy.area())/convex_hull_area);
                                DE = new_energy - old_energy;
                                if(DE > 0 ||  exp(-DE/T) > random_threshold){
                                    done = true;
                                    old_energy = new_energy;
                                }
                            }
                            if(done){    
                                break;
                            }
                        }

                        T -= 1.0/((double) L);
                        polygonChainCopy = polygonChain;

                    }
                    count++;
                    if(done){
                        done = false;
                        unchanged = false; // We changed something 1 time
                        polygonChain = polygonChainCopy;
                        break;
                    }
                }
            }
        }
    }
    return polygonChain;
}