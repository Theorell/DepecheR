
/* 
 * File:   Clusterer.cpp
 * Author: theorell
 * 
 * Created on May 3, 2017, 8:40 AM
 */

#include "Clusterer.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/src/Core/util/Constants.h>
#include <cmath>
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;



//help functions for the interfacing
RowMatrixXd numeric_to_eigen(NumericMatrix X){
  const unsigned int rows = X.nrow();
  const unsigned int cols = X.ncol();
  RowMatrixXd data = RowMatrixXd::Zero(rows,cols);
  for(unsigned int i = 0; i<rows; i++){
    for(unsigned int j = 0; j<cols; j++){
      data(i,j)=X(i,j);
    }
  }
  return data;   
}

Eigen::VectorXd numericVector_to_eigen(NumericVector X){
  const unsigned int rows = X.size();
  Eigen::VectorXd data = Eigen::VectorXd::Zero(rows);
  for(unsigned int i = 0; i<rows; i++){
    data(i)=X(i);
  }
  return data;   
}

Eigen::Matrix<unsigned int, -1,1>  integer_to_eigen(IntegerVector X){
  const unsigned int rows = X.size();
  Eigen::Matrix<unsigned int, -1,1>  data = Eigen::Matrix<unsigned int, -1,1>::Zero(rows);
  for(unsigned int i = 0; i<rows; i++){
    data(i)=X(i);
  }
  return data;   
}

NumericMatrix eigen_to_numeric(RowMatrixXd X){
  const unsigned int rows = X.rows();
  const unsigned int cols = X.cols();
  NumericMatrix data(rows,cols);
  for(unsigned int i = 0; i<rows; i++){
    for(unsigned int j = 0; j<cols; j++){
      data(i,j) = X(i,j);
    }
  }
  return data;
}
//Interfacing with R!
// [[Rcpp::export]]
List sparse_k_means(NumericMatrix X, const unsigned int k, const double reg, const bool no_zero,const unsigned long seed_off_set = 0) {
  //use eigen internally
  const unsigned int rows = X.nrow();
  //const unsigned int cols = X.ncol();
  RowMatrixXd data = numeric_to_eigen(X);
  
  Clusterer c = Clusterer();
  c.reseed(seed_off_set);
  const Return_values partitioning =c.find_centers(data, k, reg, no_zero);
  IntegerVector ints(rows);
  IntegerVector ints_no_zero(rows);
  
  for(unsigned int i = 0; i<rows; i++){
    ints(i) = partitioning.indexes(i);
    ints_no_zero(i) = partitioning.indexes_no_zero(i);
  }
  NumericMatrix centers = eigen_to_numeric(partitioning.centers);
  NumericMatrix centers_no_zero = eigen_to_numeric(partitioning.centers_no_zero);
  
  List ret;
  ret["i"]=ints;
  ret["o"]=ints_no_zero;
  ret["c"]=centers;
  ret["v"]=centers_no_zero;
  ret["n"]=partitioning.norm;
  ret["m"]=partitioning.norm_no_zero;
  //std::cout<<"the norm is: "<< partitioning.norm<<std::endl;
  return ret;
}



// [[Rcpp::export]]
List grid_search(NumericMatrix X, IntegerVector k, NumericVector reg, const unsigned int iterations, const unsigned int bootstrapSamples, const unsigned long seed_off_set = 0) {
  //const unsigned int rows = X.nrow();
  
  RowMatrixXd data = numeric_to_eigen(X);
  Eigen::Matrix<unsigned int, -1,1>  k_vector = integer_to_eigen(k);
  Eigen::VectorXd reg_vector = numericVector_to_eigen(reg);
  Clusterer c = Clusterer();
  c.reseed(seed_off_set);
  Optimization_values vals = c.optimize_param(data,k_vector,reg_vector,iterations, bootstrapSamples);
  NumericMatrix distances = eigen_to_numeric(vals.distances);
  NumericMatrix distances_no_zero = eigen_to_numeric(vals.distances_no_zero);
  NumericMatrix found_cluster = eigen_to_numeric(vals.found_cluster);
  NumericMatrix found_cluster_no_zero = eigen_to_numeric(vals.found_cluster_no_zero);
  std::vector<std::vector<NumericMatrix> > numCenters(vals.centers.size());
  for(unsigned int j = 0; j<vals.centers.size();j++){
    //numCenters.push_back();
    for(unsigned int i = 0; i<4;i++){
      numCenters[j].push_back(eigen_to_numeric(vals.centers[j][i]));
    }
  }
  
  List ret;
  
  
  
  ret["d"] = distances;
  ret["z"] = distances_no_zero;
  ret["n"] = found_cluster;
  ret["m"] = found_cluster_no_zero;
  ret["c"] = numCenters;
  return ret;
  
}

// [[Rcpp::export]]
List allocate_points(NumericMatrix X, NumericMatrix mu,const bool no_zero ) {
  
  Clusterer c = Clusterer();
  RowMatrixXd data = numeric_to_eigen(X);//c.m_rescale(numeric_to_eigen(X));
  RowMatrixXd eigen_mu = numeric_to_eigen(mu);
  Eigen::VectorXi indices = c.allocate_clusters(data, eigen_mu, no_zero);
  unsigned int rows = X.rows();
  IntegerVector inds(rows);
  for(unsigned int i = 0; i<rows; i++){
    inds(i) = indices(i);
  }
  List ret;
  ret["i"]=inds;
  return ret;
}

// [[Rcpp::export]]
const double rand_index(IntegerVector inds1,IntegerVector inds2, unsigned int k){
  Clusterer c = Clusterer();
  const unsigned int rows = inds1.size();
  Eigen::VectorXi inds1_eigen = Eigen::VectorXi::Zero(rows);
  Eigen::VectorXi inds2_eigen = Eigen::VectorXi::Zero(rows);
  for(unsigned int i = 0; i<rows; i++){
    inds1_eigen(i)= inds1(i);
    inds2_eigen(i)= inds2(i);
  }
  return c.cluster_distance(inds1_eigen,inds2_eigen,k);
}
