#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
#include <stack>

using namespace std;

class face { // face of the input polyhedron


 public:
    unsigned int  orientation;
    bool visited;
    vector<int> nodes;
    //  vector<int> dir; // the directions of the edges

    face() {
        visited=false;
    }
};


class edge {

 public:



    int face1;
    int face2; // the two faces incident on the edge
    int dir1;
    int dir2; // the directions the edge is in the two faces
    int index; // index of the bigger vertex

    edge(int in) {
        index=in;
    }


};






class vertex {

 public:
    double x;
    double y;
    double z;

    vertex(double px,double py,double pz) {

        x=px;
        y=py;
        z=pz;
    }

};
