/*
 * Power Crust software, by Nina Amenta, Sunghee Choi and Ravi Krishna Kolluri.
 * Copyright (c) 2000 by the University of Texas
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee under the GNU Public License is hereby granted, 
 * provided that this entire notice  is included in all copies of any software 
 * which is or includes a copy or modification of this software and in all copies 
 * of the supporting documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#include <fstream>
#include <vector>
#include <algorithm>
#include <set>

#include <stdlib.h>

#define SQ(A) ((A)*(A))

// class node for the vertex of the graph

using namespace std;

enum polelabel {UNLABELED,EXT,INT};
enum polestatus {FIXED,PRESENT,REMOVED};
  
class vertex {

 public: 
  
    double x,y,z; // the coordinates of the pole
  
    int index ; // index of the pole
    polelabel label;
    polestatus  status; // the status in the simplification
    double radius; // the radius of the ball
    double distance; // the max distance between samples 

    vector<int> neighbors; // the neigbors in the graph

    // constructor

    vertex(double px,double py,double pz,double pr,polelabel pl,double d) {

        x=px;
        y=py;
        z=pz;

        radius=pr;
        label=pl;
        distance=d;
        status=PRESENT;
    }

};


class powershape {

  
    void addEdge(int,int);
    double distance(int,int);
  
 public:
  
    vector<vertex> verts; // the vertices of the graph
    double noiseThreshold;
    double redThreshold;
  

    void outputPlot(ofstream&);
    void readInput(ifstream& );
    void noiseRemove();
    void redRemove();

    //  int compare(const void *n1,const void *n2);
    //  int compare(int,int);

    void outputPoles(ofstream& );

    powershape(double nt,double rt) {
        noiseThreshold=nt;
        redThreshold=rt;
    }

    powershape() {
    }

};








