#include "functions.h"

using namespace std;

Polygon_2 localSearch(Polygon_2 polygonChain, int L_defined, double threshold,bool flag_min){
    
    vector<localSearchList> list;

    double initial_area = abs(polygonChain.area());

    int L = L_defined;
    int count_L = 0;
    int offset = 0;

    vector<Point_2> randomPoints;
    Polygon_2 polygonChainCopy = polygonChain;
    Polygon_2 temporaryPolygon;
    Point_2 currPoint;
    Point_2 tempPoint;
    localSearchList currElement, max_list;
    double difference = initial_area;
    while(abs(difference) >= threshold){ // Apply local search until the threshold is achieved
        L = L_defined;
        while(L>0){
            polygonChainCopy = polygonChain;
            for(auto it = polygonChain.edges_begin(); it != polygonChain.edges_end();++it){
                    for(offset = 0;offset<polygonChainCopy.size()-1;offset++){ 
                        randomPoints.clear();
                        for(auto first = polygonChainCopy.begin() + offset;first != polygonChainCopy.end();++first){ // Offset helps us in order to take (if L = 3), 1st,2nd,3rd point, then 2nd,3rd,4th point etc
                            randomPoints.push_back(*first);
                            count_L++;
                            if(count_L == L) break; // Take L points each time
                        }
                        count_L = 0;
                        for(auto iter = randomPoints.begin(); iter != randomPoints.end(); ++iter){
                            for(auto del = polygonChainCopy.begin(); del != polygonChainCopy.end();++del){
                                if(*del == *iter){
                                    polygonChainCopy.erase(del); // We delete the selected points from the chain
                                    break;
                                }
                            }
                        }
                        
                        for(auto i3 = polygonChainCopy.begin();i3 !=polygonChainCopy.end();++i3){
                            if(i3->x() == (it->target()).x() && i3->y() == (it->target()).y()){
                                temporaryPolygon = polygonChainCopy;
                                for(auto iter = randomPoints.rbegin(); iter != randomPoints.rend() ; ++iter){
                                    polygonChainCopy.insert(i3,*iter); // Insert them on between the i3 edge
                                }

                                if(flag_min == true){ // minimazation
                                    if(polygonChainCopy.is_simple() && abs(polygonChainCopy.area()) < initial_area){ // if simple and negative win
                                        currElement.edge = *it;
                                        currElement.points = randomPoints;
                                        currElement.area = abs(polygonChainCopy.area());
                                        list.push_back(currElement);
                                    }
                                }else{ // maximazation
                                    if(polygonChainCopy.is_simple() && abs(polygonChainCopy.area()) > initial_area){ // if simple and positive win
                                        currElement.edge = *it;
                                        currElement.points = randomPoints;
                                        currElement.area = abs(polygonChainCopy.area());
                                        list.push_back(currElement);
                                    }             
                                }
                                polygonChainCopy = temporaryPolygon;
                                break;
                            }
                        }
                        // Bring it to normal position
                        polygonChainCopy = polygonChain;
                    }
                }
            L--; // Decrease L each time in order to take in each iteration L, L-1, L-2, ..., 1 point
        }

        max_list.area = initial_area;
        if(flag_min == true){
            for(auto it = list.begin(); it != list.end(); ++it){ // Get the max_list with the highest area

                if(max_list.area > it->area){ // SOS! Here max_list actually means min_list!!!
                    max_list = *it;
                }
            }
        }else{
            for(auto it = list.begin(); it != list.end(); ++it){ // Get the max_list with the highest area

                if(max_list.area < it->area){
                    max_list = *it;
                }
            }
        }

        for(auto iter = max_list.points.begin(); iter != max_list.points.end(); ++iter){
            for(auto del = polygonChain.begin(); del != polygonChain.end();++del){
                if(*del == *iter){
                    polygonChain.erase(del);
                    break;
                }
            }
        }

        for(auto i3 = polygonChain.begin();i3 !=polygonChain.end();++i3){
            if(i3->x() == ((max_list.edge).target()).x() && i3->y() == ((max_list.edge).target()).y()){
                for(auto iter = (max_list.points).rbegin(); iter != (max_list.points).rend() ; ++iter){
                    polygonChain.insert(i3,*iter);
                }
                break;
            }
        }

        difference = max_list.area - initial_area;
        initial_area = max_list.area;
    }
    return polygonChain;
}
