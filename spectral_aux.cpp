#include <Rcpp.h>
#include <limits>
#include <queue>
#include <vector>
#include <algorithm>
using namespace Rcpp;

double dist_squared(const NumericMatrix::Row& x, const NumericMatrix::Row& y);

// [[Rcpp::export]]
IntegerMatrix Mnn(NumericMatrix x, int m) {
  int nrow = x.nrow();
  double big_double = std::numeric_limits<double>::max(); 
  IntegerMatrix neigh_index = IntegerMatrix(nrow, m);
  NumericVector neigh_dist = NumericVector(m,big_double);
  double d;
  for(int i=0; i<nrow; i++){
    std::fill(neigh_dist.begin(),neigh_dist.end(),big_double);
    for(int j=0; j<nrow; j++){
      if(i!=j){
        d = dist_squared(x(i, _),x(j,_));
        int k = m-1;
        //k is the index for new neighbour index (if it falls into first m values)
        while(k > -1 && d < neigh_dist[k]){
          --k;
        }
        ++k;
        if(k < m)
        {
          for(int l=m-1; l>k; l--){
            neigh_dist[l] = neigh_dist[l-1];
            neigh_index(i,l) = neigh_index(i,l-1);
          }
          neigh_dist[k] = d;
          neigh_index(i,k) = j;
        }
      }
    }
  }
  return neigh_index;
}

// [[Rcpp::export]]
IntegerMatrix Mnn_graph(IntegerMatrix s){
  int snrow = s.nrow(), sncol = s.ncol();
  //create adjacency matrix
  IntegerMatrix g = IntegerMatrix(snrow,snrow);
  for(int i=0; i<snrow; i++){
    for(int j=0; j<sncol; j++){
      g(i,s(i,j)) = 1;
      g(s(i,j),i) = 1;
    }
  }
  
  //check if connected, BFS
  std::queue<int> open = std::queue<int>();
  std::vector<int> compind = std::vector<int>(snrow,0);
  int ind=0;
  int cmpnent=0;
  int vcount=0;
  while(vcount<snrow){
    if(open.empty()){
      cmpnent++;
      while(compind[ind]!=0){
        ind++;
      }
      compind[ind]=cmpnent;
      open.push(ind);
      vcount++;
    }
    int& current = open.front();
    for(int i=0; i<snrow; i++){
      if(g(current,i)==1){
        if(compind[i]==0){
          compind[i]=cmpnent;
          open.push(i);
          vcount++;
        }
      }
    }
    open.pop();
  }
  
  //if g not connected then add edges
  if(cmpnent>1){
    std::vector<int> reps = std::vector<int>(cmpnent);
    for(int i=0; i<snrow; i++){
      reps[compind[i]-1]=i;
    }
    for(int i=0; i<cmpnent-1; i++){
      g(reps[i],reps[i+1])=1;
      g(reps[i+1],reps[i])=1;
    }
  }
  return g;
}

double dist_squared(const NumericMatrix::Row& x, const NumericMatrix::Row& y){
  double d = 0;
  int len = x.size();
  for(int i = 0; i < len; i++){
    d += (x[i] - y[i])*(x[i] - y[i]);
  }
  return d;
}