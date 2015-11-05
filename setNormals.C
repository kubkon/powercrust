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
#include <iostream>
#include <cstring>
#include <vector>
#include <set>
#include <stack>
#include "ndefs.h"
#define CORRECT 2
#define WRONG 3

using namespace std;

bool reverseFaces=false;
char inFileName[32]="input",outFileName[32]="output";

class neighbors {

 public:
    vector<edge> neigh;
};


class powercrust {  // class that contains the whole powercrust.


    void dealEdge(int,int,int,stack<int>&);

 public:
    vector<vertex> verts;
    vector<neighbors> edges;
    vector<face> faces;
    void readInput(ifstream&);
    void addFace(vector<int>&);
    void addEdge(int,int,int);
    void correctOrientation();
    void reverseAll();
    void writeOutput(ofstream&);


    int numedges;
    int numfaces;
    powercrust() {
        numedges=0;
        numfaces=0;
    }

};


void powercrust::addFace(vector<int>& indices){

    face currFace;

    numfaces++;
    for(int i=0;i<indices.size()-1;i++) {
        currFace.nodes.push_back(indices[i]);
        addEdge(indices[i],indices[i+1],numfaces-1);
    }
    currFace.nodes.push_back(indices[indices.size()-1]);

    if(indices.size()>2)
        addEdge(indices[indices.size()-1],indices[0],numfaces-1);

    faces.push_back(currFace);

}

void powercrust::addEdge(int index1,int index2,int faceindex) {

    int sindex,lindex,dir;
    if(index1<index2) {
        sindex=index1;
        lindex=index2;
        dir=0;
    }
    else {
        sindex=index2;
        lindex=index1;
        dir=1;
    }

    // check if this edge is already present
    int eindex=0;
    for(;eindex<edges[sindex].neigh.size();eindex++) {
        if(edges[sindex].neigh[eindex].index==lindex)
            break;
    }

    if(eindex<edges[sindex].neigh.size()) { // if a match was found
        edges[sindex].neigh[eindex].face2=faceindex;
        edges[sindex].neigh[eindex].dir2=dir;
    }
    else {
        edge currEdge(lindex);
        currEdge.face1=faceindex;
        currEdge.dir1=dir;
        edges[sindex].neigh.push_back(currEdge);
    }



}



void powercrust::readInput(ifstream& inFile) {

    // readin the powercrust

    char temp[12];
    inFile>>temp;
    int nverts,nfaces,nedges;
    inFile>>nverts;
    inFile>>nfaces;
    inFile>>nedges;
    double tx,ty,tz;
    vector<int> indices;
    int nindices,tempint;

    for(int i=0;i<nverts;i++) {
        inFile>>tx;
        inFile>>ty;
        inFile>>tz;
        //  cout<<tx<<"\t"<<ty<<"\t"<<tz<<"\n";

        verts.push_back(vertex(tx,ty,tz));
        edges.push_back(neighbors());

    }

    // now start reading in the faces

    for(int i=0;i<nfaces;i++) {
        inFile>>nindices;
        //cout<<i<<"\t"<<numedges<<"\n";
        indices.clear();
        for(int j=0;j<nindices;j++){
            inFile>>tempint;
            indices.push_back(tempint);
        }
        addFace(indices);
    }

}


void powercrust::correctOrientation() {

    int Or[2]; // number of faces with each orientation
    stack<int> remFaces;
    stack<int> allFaces;

    //  remFaces.push(0);
    int currFace,curror;
    int currsize,corr;
    int e1,e2;
    Or[0]=0;
    Or[1]=0;

    for(int i=0;i<faces.size();i++) {

        if(faces[i].visited==false) {
            Or[0]=0;
            Or[1]=0;
            faces[i].orientation=0;
            faces[i].visited=true;
            cout<<"starting a new component ..\n";
            remFaces.push(i);

            while(!remFaces.empty()) {
                currFace=remFaces.top();
                allFaces.push(currFace);
                remFaces.pop();
                curror=faces[currFace].orientation;
                Or[curror]++;
                currsize=faces[currFace].nodes.size();
                for(int i=0;i<currsize;i++) {
                    dealEdge(faces[currFace].nodes[i],
                             faces[currFace].nodes[(i+1)%currsize],currFace,remFaces);
                }

            }
            // done with the whole component , now correct the orientation
            cout<<"The number of faces with the two orientations are "<<Or[0]<<"\t"<<Or[1]<<"\n";
            if(Or[0]>Or[1])
                corr=0;
            else
                corr=1;
            while(!allFaces.empty()) {
                currFace=allFaces.top();
                allFaces.pop();
                if(faces[currFace].orientation==corr)
                    faces[currFace].orientation=CORRECT;
                else
                    faces[currFace].orientation=WRONG;
            }

        }
    }

}



void powercrust::dealEdge(int index1,int index2,int faceindex,stack<int>& remFaces) {

    int sindex=(index1<index2)?index1:index2;
    int lindex=(index1<index2)?index2:index1;

    // cout<<"the sindex is \t"<<sindex<<"\n";
    //cout<<"the lindex is \t"<<lindex<<"\n";
    int eindex=0;
    for(;eindex<edges[sindex].neigh.size();eindex++) {
        if(edges[sindex].neigh[eindex].index==lindex)
            break;
    }

    // the required edge is given by eindex
    if(eindex==edges[sindex].neigh.size())
        cout<<"Something is wrong ...\n";

    int newface;
    if(faceindex== edges[sindex].neigh[eindex].face1)
        newface= edges[sindex].neigh[eindex].face2;
    else
        newface= edges[sindex].neigh[eindex].face1;
    if(!faces[newface].visited) {
        faces[newface].visited=true;
        if( edges[sindex].neigh[eindex].dir1!= edges[sindex].neigh[eindex].dir2)
            faces[newface].orientation=faces[faceindex].orientation;
        else
            faces[newface].orientation=1-faces[faceindex].orientation;
        remFaces.push(newface);
    }

}


void powercrust::reverseAll() {

    for(int i=0;i<faces.size();i++)
        faces[i].orientation=WRONG;

}

void powercrust::writeOutput(ofstream& outFile) {

    cout<<"Writing output..\n";
    outFile<<"OFF\n";
    outFile<<verts.size()<<"\t"<<faces.size()<<"\t0\n";
    outFile.precision(12);
    for(int i=0;i<verts.size();i++)
        outFile<<verts[i].x<<"\t"<<verts[i].y<<"\t"<<verts[i].z<<"\n";

    for(int i=0;i<faces.size();i++) {
        if(faces[i].orientation==CORRECT) {
            outFile<<faces[i].nodes.size()<<"\t";
            for(int j=0;j<faces[i].nodes.size();j++)
                outFile<<faces[i].nodes[j]<<"\t";
            outFile<<"\n";
        }
        else {
            if(!reverseFaces)
                cout<<"Reversing..\n";

            outFile<<faces[i].nodes.size()<<"\t";
            for(int j=faces[i].nodes.size()-1;j>=0;j--)
                outFile<<faces[i].nodes[j]<<"\t";
            outFile<<"\n";
        }

    }

}




int parseCommand(int argc,char **argv) {

    for(int i=1;i<argc;) {
        //  cout<<"Argv is\t"<<argv[i]<<"\n";

        if(argv[i][0]!='-')
            return 0;
        switch(argv[i][1]) {
        case 'o':strcpy(outFileName,argv[i+1]);
            i=i+2;
            break;
        case 'i':strcpy(inFileName,argv[i+1]);
            i=i+2;
            break;

        case 'r':reverseFaces=true;
            i=i+1;
            break;


        default:return 0;
        }

    }
    return 1;

}


int main(int argc,char** argv) {

    powercrust crust;
    ifstream inFile;
    ofstream outFile;
    if(!parseCommand(argc,argv)) {
        cout<<"Error in arguments..\n";
        exit(0);
    }

    inFile.open(inFileName);
    outFile.open(outFileName);

    if(!inFile){
        cout<<"Cannot open input File\n";
        exit(0);
    }
    if(!outFile){
        cout<<"Cannot open output File\n";
        exit(0);
    }



    crust.readInput(inFile);
    cout<<"Done with the reading .\n";

    if(crust.faces.size()==0)
        exit(0);

    if(!reverseFaces)
        crust.correctOrientation();
    else
        crust.reverseAll();

    crust.writeOutput(outFile);


}
