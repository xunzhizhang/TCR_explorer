#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>            
using namespace Rcpp;

int  BLOSUM62[128][128];     
// [[Rcpp::export]]
int loadBLOSUM62(){
  for (int i = 0; i < 128; ++i) {
    for (int j = 0; j < 128; ++j) {
      BLOSUM62[i][j] = 0;
    }
  }
  std::ifstream file("blosum62.3col.txt");
  std::string   line;
  
  while(std::getline(file, line)){
    std::stringstream   linestream(line);
    std::string         a1;
    std::string         a2;
    int                 val;
    std::getline(linestream, a1, '\t');
    std::getline(linestream, a2,'\t');
    linestream >> val;
    BLOSUM62[(unsigned char) a1[0]][(unsigned char) a2[0]] = val;  
  }
  return 0;
}

// [[Rcpp::export]]
int di(const std::string& v1, const std::string& v2){
  size_t len = v1.size();
  if (len > v2.size()) {
    len = v2.size();
  }
  int d = 0;
  for (size_t i = 0 ; i != len; ++i) {
    d += BLOSUM62[(unsigned char)v1[i]][(unsigned char)v2[i]];
  }
  return d;
}
// [[Rcpp::export]]
NumericMatrix score(CharacterVector s, size_t len){
  loadBLOSUM62();
  std::vector<std::string> seq = Rcpp::as <std::vector<std::string> > (s);
  NumericMatrix smatrix(len, len);
  for (size_t i = 0; i < len; i++){
    for (size_t j = 0; j <= i; j++){
      smatrix(i, j) = di(seq[i], seq[j]);
    }
  }
  return smatrix;
}


// [[Rcpp::export]]
NumericMatrix loc_m(CharacterVector s, size_t len){
  std::vector<std::string> seq = Rcpp::as <std::vector<std::string> > (s);
  NumericMatrix dmatrix = score (s, len);
  for (size_t i = 0; i < len; i++){
    for (size_t j = 0; j < i; j++){
      dmatrix(i, j) = 1.0 - 2.0 * di(seq[i], seq[j]) / (dmatrix(i,i) + dmatrix(j, j));
      dmatrix(j, i) = dmatrix(i, j);
    }
  }
  for (size_t i = 0; i < len; i++){
    dmatrix(i, i) = 0.0;
  }
  
  return dmatrix;
}

// [[Rcpp::export]]
IntegerVector colv(CharacterVector s, size_t len){
  std::vector<std::string> seq = Rcpp::as <std::vector<std::string> > (s);
  NumericMatrix smatrix = score (s, len);
  const int THREAHOLD = 20;         
  int totalCluster = 1;  
  IntegerVector cluster(len, -1);
  cluster[0] = 1;
  int maxScore = -1000;
  int maxIdx = -1000;
  const int NUM_THREADS = 16;
#ifdef _OPENMP
  omp_set_num_threads(NUM_THREADS);
#endif
  int j;
  IntegerVector d(len);
  for (int i = 1; i < len; ++i){
    maxScore = -1000;
#pragma omp parallel for
    for (j = 0; j <= i - 1; ++j){
      if (smatrix(i,j) > maxScore){
        maxScore = smatrix(i,j);
        maxIdx = j;
      }
    }
    if (maxScore > THREAHOLD){         
      cluster[i] = cluster[maxIdx];
    } else {
      totalCluster ++;
      cluster[i] = totalCluster;
    }
    
  }
  return cluster;
}
