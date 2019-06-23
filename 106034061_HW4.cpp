//
//  106034061_HW4.cpp
//  Program
//
//  Created by 曾靖渝 on 2019/4/29.
//  Copyright © 2019年 曾靖渝. All rights reserved.
//
//  Reference:https://www.techiedelight.com/travelling-salesman-problem-using-branch-and-bound/
//
#include <algorithm>
#include <cstring>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <utility>
#include <climits>
#define BLOCK -1
#define INF INT_MAX
using namespace std;
//////////////////////////////////////
// State Space Tree nodes
struct Node
{
    // stores edges of state space tree helps in tracing path when answer is found
    vector<pair<int, int> > path;
    // stores the reduced matrix
    int reducedMatrix[16][16];
    // stores the lower bound
    int cost;
    //stores current city number
    int vertex;
    // stores number of cities visited so far
    int level;
};

//int N;
//Generate Test case
int** generate_test(int N){
    int **adj;
    int seed;
    seed = (int)time(NULL);
    srand(seed);
    adj=(int **)malloc(sizeof(int*)*(N));
    for(int i=0;i<N;i++)
        adj[i]=(int *)malloc(sizeof(int)*(N));
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            if(i!=j)
                adj[i][j]=rand()%100;
            else
                adj[i][j]=BLOCK;
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++){
            if(j!=N-1)
                cout<<adj[i][j]<<" ";
            else
                cout<<adj[i][j]<<endl;
        }
    
    return adj;
}
//using Brutal Force
class TSP_BF{
public:
    TSP_BF(int N, int **adj){
        this->final_res=INT_MAX;
        this->N=N;
        this->adj=adj;
        final_path=new int(N+1);
        memset(final_path, 0, sizeof(final_path));
        visited =new bool(N);
        memset(visited, false, sizeof(visited));
    }
    // implementation of traveling Salesman Problem
    int TSP(int s)
    {
        // store all vertex apart from source vertex
        vector<int> vertex;
        for (int i = 0; i < N; i++)
            if (i != s)
                vertex.push_back(i);
        // store minimum weight Hamiltonian Cycle.
        int min_path = INT_MAX;
        do {
            // store current Path weight(cost)
            int current_pathweight = 0;
            // compute current path weight
            int k = s;
            for (int i = 0; i < vertex.size(); i++) {
                current_pathweight += adj[k][vertex[i]];
                k = vertex[i];
            }
            current_pathweight += adj[k][s];
            // update minimum
            min_path = min(min_path, current_pathweight);
            if(min_path==current_pathweight){
                for(int i=0;i<N;i++){
                    final_path[i+1]=vertex[i];
                }
            }
        } while (next_permutation(vertex.begin(), vertex.end()));
        final_res=min_path;
        return min_path;
    }
    void TSP_BF_output(){
        clock_t  start, end;
        start=clock();
        TSP(0);
        end=clock();
        time=(double)end-(double)start;
        cout<<"Solution : ";
        for (int i=0; i<=N; i++)
            if(i!=N)
                cout<<this->final_path[i]<<" - ";
            else
                cout<<this->final_path[0]<<"\n";
        cout<<"Cost     : "<<this->final_res<<endl;
        cout<<"Time     : "<<time/CLOCKS_PER_SEC<<" s\n";
    }
    void TSP_BF_foutput(){
        ofstream fout1;
        fout1.open("BF.txt");
        if(!fout1)cout<<"Output1 Error..."<<endl;
        clock_t start, end;
        start=clock();
        TSP(0);
        end=clock();
        time=(double)end-(double)start;
        fout1<<"Solution : ";
        for (int i=0; i<=N; i++)
            if(i!=N)
                fout1<<this->final_path[i]<<" - ";
            else
                fout1<<this->final_path[0]<<"\n";
        fout1<<"Cost     : "<<this->final_res<<endl;
        fout1<<"Time     : "<<time/CLOCKS_PER_SEC<<" s\n";
        fout1.close();
    }
private:
    int N;
    int *final_path;
    bool *visited;
    int final_res;
    int **adj;
    double time;
};

// using Branch and Bound.
class TSP_BB{
public:
    TSP_BB(int N, int** iadj){
        this->N=N;
        for(int i=0;i<N;i++)
            for(int j=0;j<N;j++)
                this->adj[i][j]=iadj[i][j];
        for(int i=0;i<N;i++)
            final_path[i]=0;
    }
    // Function to allocate a new node (i, j) corresponds to visiting city j from city i
    Node* newNode(int parentMatrix[16][16], vector<pair<int, int> > const &path,int level, int i, int j){
        Node* node = new Node;
        // stores ancestors edges of state space tree
        node->path = path;
        // skip for root node
        if (level != 0)
            // add current edge to path
            node->path.push_back(make_pair(i, j));
        // copy data from parent node to current node
        memcpy(node->reducedMatrix, parentMatrix, sizeof node->reducedMatrix);
        // Change all entries of row i and column j to infinity skip for root node
        for (int k = 0; level != 0 && k < N; k++){
            // set outgoing edges for city i to infinity
            node->reducedMatrix[i][k] = INF;
            // set incoming edges to city j to infinity
            node->reducedMatrix[k][j] = INF;
        }
        // Set (j, 0) to infinity ,here start node is 0
        node->reducedMatrix[j][0] = INF;
        // set number of cities visited so far
        node->level = level;
        // assign current city number
        node->vertex = j;
        return node;
    }
    
    // Function to reduce each row in such a way that there must be at least one zero in each row
    void rowReduction(int reducedMatrix[16][16], int row[16])
    {
        // initialize row array to INF
        fill_n(row, N, INF);
        // row[i] contains minimum in row i
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if (reducedMatrix[i][j] < row[i])
                    row[i] = reducedMatrix[i][j];
        // reduce the minimum value from each element in each row
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if (reducedMatrix[i][j] != INF && row[i] != INF)
                    reducedMatrix[i][j] -= row[i];
    }
    // Function to reduce each column in such a way that
    // there must be at least one zero in each column
    void columnReduction(int reducedMatrix[16][16], int col[16])
    {
        // initialize col array to INF
        fill_n(col, N, INF);
        // col[j] contains minimum in col j
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if (reducedMatrix[i][j] < col[j])
                    col[j] = reducedMatrix[i][j];
        // reduce the minimum value from each element in each column
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                if (reducedMatrix[i][j] != INF && col[j] != INF)
                    reducedMatrix[i][j] -= col[j];
    }
    // Function to get the lower bound on
    // on the path starting at current min node
    int calculateCost(int reducedMatrix[16][16])
    {
        // initialize cost to 0
        int cost = 0;
        // Row Reduction
        int row[16];
        rowReduction(reducedMatrix, row);
        //        for(int i=0;i<N;i++)
        //            cout<<row[i]<<" ";
        //        cout<<endl;
        // Column Reduction
        int col[16];
        //        for(int i=0;i<N;i++)
        //            cout<<col[i]<<" ";
        //        cout<<endl;
        columnReduction(reducedMatrix, col);
        //        for(int i=0;i<N;i++)
        //            cout<<col[i]<<" ";
        //        cout<<endl;
        // the total expected cost is the sum of all reductions
        for (int i = 0; i < N; i++){
            cost += (row[i] != INT_MAX) ? row[i] : 0;
            cost += (col[i] != INT_MAX) ? col[i] : 0;
        }
        return cost;
    }
    void copy_path(vector<pair<int, int> > const &list){
        for (int i = 0; i < list.size(); i++){
            final_path[i]=list[i].first;
        }
        final_path[N]=list[0].first;
    }
    // Comparison object to be used to order the heap
    struct comp {
        bool operator()(const Node* lhs, const Node* rhs) const{
            return lhs->cost > rhs->cost;}
    };
    // Function to solve Traveling Salesman Problem using Branch and Bound
    int TSP(int costMatrix[16][16])
    {
        // Create a priority queue to store live nodes of search tree;
        priority_queue<Node*, std::vector<Node*>, comp> pq;
        vector<pair<int, int> > v;
        // create a root node and calculate its cost
        // The TSP starts from first city i.e. node 0
        Node* root = newNode(costMatrix, v, 0, -1, 0);
        // get the lower bound of the path starting at node 0
        root->cost = calculateCost(root->reducedMatrix);
        // Add root to list of live nodes;
        pq.push(root);
        // Finds a live node with least cost, add its children to list of live nodes and finally deletes it from the list
        while (!pq.empty()){
            // Find a live node with least estimated cost
            Node* min = pq.top();
            // The found node is deleted from the list of live nodes
            pq.pop();
            // i stores current city number
            int i = min->vertex;
            // if all cities are visited
            if (min->level == N - 1)
            {
                // return to starting city
                min->path.push_back(make_pair(i, 0));
                copy_path(min->path);
                // return optimal cost
                return min->cost;
            }
            // do for each child of min (i, j) forms an edge in space tree
            for (int j = 0; j < N; j++)
            {
                if (min->reducedMatrix[i][j] != INF)
                {
                    // create a child node and calculate its cost
                    Node* child = newNode(min->reducedMatrix, min->path, min->level + 1, i, j);
                    /* Cost of the child =cost of parent node+cost of the edge(i, j)+lower bound of the path starting at node j */
                    child->cost = min->cost + min->reducedMatrix[i][j]+calculateCost(child->reducedMatrix);
                    // Add child to list of live nodes
                    pq.push(child);
                }
            }
            delete min;
        }
        return -1;
    }
    void TSP_BB_output(){
        clock_t  start, end;
        start=clock();
        final_res=TSP(this->adj);
        end=clock();
        time=(double)end-(double)start;
        cout<<"Solution : ";
        for (int i=0; i<=N; i++)
            if(i!=N)
                cout<<this->final_path[i]<<" - ";
            else
                cout<<this->final_path[0]<<"\n";
        cout<<"Cost     : "<<this->final_res<<endl;
        cout<<"Time     : "<<time/CLOCKS_PER_SEC<<" s\n";
    }
    void TSP_BB_foutput(){
        ofstream fout2;
        fout2.open("BB.txt");
        if(!fout2)cout<<"Output2 Error..."<<endl;
        clock_t start, end;
        start=clock();
        final_res=TSP(this->adj);
        end=clock();
        time=(double)end-(double)start;
        fout2<<"Solution : ";
        for (int i=0; i<=N; i++)
            if(i!=N)
                fout2<<this->final_path[i]<<" - ";
            else
                fout2<<this->final_path[0]<<"\n";
        fout2<<"Cost     : "<<this->final_res<<endl;
        fout2<<"Time     : "<<time/CLOCKS_PER_SEC<<" s\n";
        fout2.close();
    }
    
private:
    int final_path[16];
    int adj[16][16];
    int N;
    double time;
    int final_res;
};
//////////////////////////////////////
// Driver code
int main()
{
    int N;
    int **adj;
    fstream  fout1,fout2;
    //    fstream fin;
    ifstream fin;
    fin.open("input.txt");
    fout1.open("BF.txt",ios::out);
    fout2.open("BB.txt",ios::out);
    
    if(!fin)cout<<"File Input Error"<<endl;
    if(!fout1)cout<<"File Output1 Error"<<endl;
    if(!fout2)cout<<"File Output2 Error"<<endl;
    
    cout<<"Start..."<<endl;
    //    cin>>N;
    //    adj=(int **)malloc(sizeof(int*)*N);
    //    for(int i=0;i<N;i++)
    //        adj[i]=(int *)malloc(sizeof(int)*N);
    //    for(int i=0;i<N;i++)
    //        for(int j=0;j<N;j++)
    //        {
    //            int cost;
    //            cin>>cost;
    //            adj[i][j]=cost;
    //        }
    fin>>N;
    adj=(int **)malloc(sizeof(int*)*N);
    for(int i=0;i<N;i++)
        adj[i]=(int *)malloc(sizeof(int)*N);
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
        {
            int cost;
            fin>>cost;
            adj[i][j]=cost;
        }
    //        adj=generate_test(N);
    TSP_BB tsp_BB(N, adj);
    TSP_BF tsp_BF(N, adj);
    tsp_BB.TSP_BB_foutput();
    tsp_BF.TSP_BF_foutput();
    return 0;
}
/*
 7
 -1 3 93 13 33 9 57
 4 -1 77 42 21 16 34
 45 17 0 36 16 28 25
 39 90 80 0 56 7 91
 28 46 88 33 0 25 57
 3 88 18 46 92 0 7
 44 26 33 27 84 39 0
 
 4
 -1 10 15 20
 10 -1 35 25
 15 35 -1 30
 20 25 30 -1
 
 4
 -1 4 1 3
 4 -1 2 1
 1 2 -1 5
 3 1 5 -1
 */
