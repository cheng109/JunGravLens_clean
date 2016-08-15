#include <vector> 
#include <climits> 
#include <iostream> 
#include <utility> 
#include <cstdlib> 
#include <ctime> 
#include "Kmeans.h"
using namespace std; 

Kmeans::Kmeans(int _k, vector<double>& _xPos, vector<double>& _yPos, vector<int>& _group, vector<double>& _w): k(_k), n(_xPos.size()), 
						centers(_k), xPos(_xPos), yPos(_yPos), group(_group), w(_w) {
		// initite centers; 
		for(int i=0; i<k ; ++i)
			centers[i] = {_xPos[i], _yPos[i]};  		
	}

void Kmeans::initiate() {  
		double xSmallest = INT_MAX; 
		double xLargest = INT_MIN; 
		for(int i=0; i<n; ++i) {
			if(xPos[i]<xSmallest) {
				xSmallest = xPos[i]; 
				centers[0] = {xPos[i], yPos[i]}; 
			}
			if(xPos[i] >xLargest) {
				xLargest = xPos[i]; 
				centers[k-1] = {xPos[i], yPos[i]} ; 
			}
		}
		if(k==3)  
			centers[1] = {0.5*(centers[0].first + centers[k-1].first),  0.5*(centers[0].second + centers[k-1].second)} ;   

	}

void Kmeans::updateGroup() {
		for(int i=0; i<n; ++i) {
			int minDist = INT_MAX; 
			for(int j=0; j<k; ++j) {
				int xDist = xPos[i] -centers[j].first; 
				int yDist = yPos[i] -centers[j].second; 
				if(xDist*xDist + yDist*yDist < minDist) {
					group[i] = j; 
					minDist = xDist*xDist + yDist*yDist; 

				}
			}
		}
	}


struct Triplet{
	double xPos; 
	double yPos; 
	double w; 
	Triplet(double _xPos, double _yPos, double _w): xPos(_xPos), yPos(_yPos), w(_w) {
	}
}; 


void Kmeans::updateCenters(){
		//centers.resize(k); 
		vector<vector<Triplet>> temp (k); 
		vector<double> weightByGroup(k, 0); 
		for(int i=0; i<n; ++i) {
			temp[group[i]].push_back(Triplet(xPos[i], yPos[i], w[i])) ;
			weightByGroup[group[i]] += w[i]; 
		}
		double totalW = 0; 
		for(int i=0; i<k; ++i) {
			// update the centroid for each group; 
			centers[i] = {0, 0}; 
			for(int j=0 ;j<temp[i].size(); ++j) {
				centers[i].first  +=(temp[i][j].w * temp[i][j].xPos); 
				centers[i].second +=(temp[i][j].w * temp[i][j].yPos); 
			}
			// summation of all the weight; 
			centers[i].first  /= weightByGroup[i]; 
			centers[i].second /= weightByGroup[i];
		}

	}

vector<int> Kmeans::countGroup() {
		vector<int> ret(k, 0); 
		for(int i=0; i<n; ++i) {
			ret[group[i]] ++; 
		}
		return ret; 
	}

// int main() {

// 	srand((int)time(0));
// 	int n= 10000; 
// 	int k = 4; 
// 	vector<double> x(n); 
// 	vector<double> y(n); 
// 	for(int i=0 ; i<n; ++i) {
// 		x[i] = rand()%100 ;
// 		y[i] = rand()%100; 
// 		if(i<50) {
// 			x[i] += 0; 
// 			y[i] += 0; 
// 		}
// 	} 

// 	Kmeans kmeans(k, x, y); 
// 	kmeans.initiate(); 


// 	double diff = 0;
// 	vector<pair<double, double>> pre_centers ; 
// 	for(int i=0; i<20; ++i) {

// 		cout << i << ":\t" << kmeans.centers[0].first << "\t" << kmeans.centers[0].second 
// 				 << "\t" << kmeans.centers[1].first << "\t" << kmeans.centers[1].second 
// 				 << "\t" << diff << endl;

// 		kmeans.updateGroup(); 

// 		// print out the group info; 
// 		vector<int> count = kmeans.countGroup(); 
// 		// for(int i=0; i<k; ++i) 
// 		// 	cout <<  count[i]  << "\t"; 
// 		// cout << endl; 
// 		pre_centers = kmeans.centers; 
// 		kmeans.updateCenters(); 
// 		diff = 0; 
// 		for(int i=0; i<k; ++i) {
// 			int xdiff = pre_centers[i].first - kmeans.centers[i].first; 
// 			int ydiff = pre_centers[i].second -kmeans.centers[i].second; 
// 			diff += (xdiff * xdiff + ydiff * ydiff );  
// 		}

// 	}

// 	return 0; 
// }









