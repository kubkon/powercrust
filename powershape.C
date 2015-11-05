#include <iostream>
#include <cstring>
#include <math.h>
#include "sdefs.h"


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

powershape pshape; // a global variable
char inFileName[32]="input",outFileName[32]="output";
bool plotGraph=true;
double nt=0,rt=0;
bool outputUnlabeled=false;


void powershape::addEdge(int i,int j) {
    if(find(verts[i].neighbors.begin(),verts[i].neighbors.end(),j)==
       verts[i].neighbors.end()){
        verts[i].neighbors.push_back(j);
        verts[j].neighbors.push_back(i);

    }
}


void powershape::readInput(ifstream& inFile) {

    int nverts,nfaces;

    inFile>>nverts;
    inFile>>nfaces;
    double x,y,z,r,d;
    int stat;
    polelabel pl;
    vector<int> findices;
    int nfverts;
    int temp;

    for(int i=0;i<nverts;i++) {
        inFile>>x;
        inFile>>y;
        inFile>>z;
        inFile>>r;
        inFile>>stat;
        inFile>>d;

        if(stat==0)
            pl=UNLABELED;
        else
            if(stat==1)
                pl=EXT;
            else
                pl=INT;


        verts.push_back(vertex(x,y,z,r,pl,d));
    }

    // now for the faces
    for(int i=0;i<nfaces;i++) {
        findices.clear();
        inFile>>nfverts;
        for(int j=0;j<nfverts;j++) {
            inFile>>temp;
            findices.push_back(temp);
        }

        // actually add edges to the graph

        for(int j=0;j<findices.size()-1;j++)
            addEdge(findices[j],findices[j+1]);

        if(findices.size()>2)
            addEdge(findices[0],findices[findices.size()-1]);

    }
}



void powershape::noiseRemove() {

    int numPoints=0;

    for(int i=0;i<verts.size();i++)
        if(verts[i].distance<= noiseThreshold){
            verts[i].status=REMOVED;
            numPoints++;
        }
    cout<<"Number of points removed by the noise criterion\t"<<numPoints<<"\n";
}


int compare(const void *n1,const void *n2){
    //int compare(int index1,int index2) {

    int index1=*((int *)n1);
    int index2=*((int *)n2);

    if(index1 >= pshape.verts.size() || index2 >= pshape.verts.size()){
        cout<<"Error!..\n"<<index1<<"\t"<<index2<<"\n";
        return 0;
    }

    if(pshape.verts[index1].radius > pshape.verts[index2].radius) return -1;
    if(pshape.verts[index1].radius <  pshape.verts[index2].radius) return 1;

    return 0;
}


int compareDistance(const void *n1,const void *n2){
    //int compare(int index1,int index2) {

    int index1=*((int *)n1);
    int index2=*((int *)n2);

    if(index1 >= pshape.verts.size() || index2 >= pshape.verts.size()){
        cout<<"Error!..\n"<<index1<<"\t"<<index2<<"\n";
        return 0;
    }

    if(pshape.verts[index1].distance > pshape.verts[index2].distance) return 1;
    if(pshape.verts[index1].distance <  pshape.verts[index2].distance) return -1;

    return 0;
}


void powershape::outputPlot(ofstream& outFile) {

    int *pqueue=(int*)malloc(verts.size()*sizeof(int));

    for(int i=0;i<verts.size();i++) pqueue[i]=i;

    qsort(pqueue,verts.size(),sizeof(int),compareDistance);


    for(int i=0;i<verts.size();i++) {
        outFile<<i<<"\t"<<verts[pqueue[i]].distance<<"\n";
    }

}

double powershape::distance(int i,int j) {

    return (sqrt(SQ(verts[i].x-verts[j].x) +
                 SQ(verts[i].y-verts[j].y) +
                 SQ(verts[i].z-verts[j].z)));

}

void powershape::redRemove() {

    int *pqueue;
    int currnode,currneighbor;
    int numPoints=0;

    pqueue=(int *)malloc(verts.size()*sizeof(int));



    for(int i=0;i<verts.size();i++)  pqueue[i]=i;

    qsort(pqueue,verts.size(),sizeof(int),compare);



    for(int i=0;i<verts.size();i++) {



        currnode=pqueue[i];

        if(verts[currnode].status!=PRESENT)
            continue;

        verts[currnode].status=FIXED;

        for(int j=0;j<verts[currnode].neighbors.size();j++) {
            currneighbor=verts[currnode].neighbors[j];
            if(verts[currneighbor].status==REMOVED || verts[currneighbor].status==FIXED)
                continue;

            if((distance(currnode,currneighbor)-
                (sqrt(verts[currnode].radius)-sqrt(verts[currneighbor].radius)))<
               redThreshold){
                verts[currneighbor].status=REMOVED;

                numPoints++;

                // add the neighbors to the original graph
                for(int k=0;k<verts[currneighbor].neighbors.size();k++)
                    addEdge(verts[currneighbor].neighbors[k],currnode);
            }

        }

        verts[currnode].neighbors.clear();

    }
    cout<<"Number of points removed by the redundancy criterion\t"<<numPoints<<"\n";
}



void powershape::outputPoles(ofstream& outFile) {

    outFile.precision(12);
    int numverts=0;
    for(int i=0;i<verts.size();i++)
        if(verts[i].status!=REMOVED &&
           (verts[i].label!=UNLABELED || outputUnlabeled))
            numverts++;

    outFile<<numverts<<"\n";




    for(int i=0;i<verts.size();i++)
        if(verts[i].status!=REMOVED &&
           (verts[i].label!=UNLABELED || outputUnlabeled))
            outFile<<verts[i].x<<"\t"
                   <<verts[i].y<<"\t"
                   <<verts[i].z<<"\t"
                   <<verts[i].radius<<"\t"
                   <<verts[i].label<<"\t"
                   <<verts[i].distance<<"\n";


}

/* parse the command line */

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
        case 'n':nt=atof(argv[i+1]);
            i=i+2;
            break;
        case 'r':rt=atof(argv[i+1]);
            i=i+2;
            break;
        case 'p':plotGraph=atoi(argv[i+1]);
            i=i+2;
            break;
        case 'u':outputUnlabeled = atoi(argv[i+1]);
            i+=2;
            break;

        default:return 0;
        }

    }
    return 1;

}







int main (int argc,char** argv) {



    if(!parseCommand(argc,argv)) {
        cout<<"Error in arguments..\n";
        exit(0);
    }

    ifstream inFile(inFileName);
    ofstream outFile(outFileName);

    pshape.readInput(inFile);

    cout<<"Number of vertices in the input File \t"<< pshape.verts.size()<<"\n";
    cout<<"The input file is\t"<<inFileName<<"\n";
    cout<<"The output file is\t"<<outFileName<<"\n";
    cout<<"The noise threshold is\t "<<nt<<"\n";
    cout<<"The redundancy threshold  is\t"<<rt<<"\n";




    ofstream plotLog("distancelog");


    pshape.noiseThreshold=nt;
    pshape.redThreshold=rt;


    if(plotGraph)
        pshape.outputPlot(plotLog);

    pshape.noiseRemove();
    pshape.redRemove();

    pshape.outputPoles(outFile);

}
