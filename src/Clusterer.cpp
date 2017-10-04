/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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
#include <eigen3/Eigen/src/Core/util/Constants.h>
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
    const unsigned int rows = X.nrow();
    
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
    List ret;
    IntegerVector inds1(rows);
    IntegerVector inds2(rows);
    IntegerVector inds1_no_zero(rows);
    IntegerVector inds2_no_zero(rows);
    
    for(unsigned int i = 0; i<rows; i++){
        inds1 = vals.indexes1(i);
        inds2 = vals.indexes2(i);
        inds1_no_zero = vals.indexes1_no_zero(i);
        inds2_no_zero = vals.indexes2_no_zero(i);
    }
    
    
    ret["d"] = distances;
    ret["z"] = distances_no_zero;
    ret["n"] = found_cluster;
    ret["m"] = found_cluster_no_zero;
    ret["e"] = inds1;
    ret["r"] = inds2;
    ret["t"] = inds1_no_zero;
    ret["y"] = inds2_no_zero;
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
Clusterer::Clusterer() {
    srand((unsigned int) time(0));
}

Clusterer::Clusterer(const Clusterer& orig) {
}

//This function allocates each point to its closest cluster


const Eigen::VectorXi Clusterer::allocate_clusters(const RowMatrixXd& X, const RowMatrixXd& mu, const bool no_zero ){
    const unsigned int mu_rows = mu.rows();
    const unsigned int X_rows = X.rows();
    
    //make a vector to iterate over
    std::vector<unsigned int> active_indices;
    for(unsigned int i = 0; i< mu_rows; i++){
        if(no_zero){
            if(!mu.row(i).isZero(0)){
                //std::cout<<"muRow: "<< mu.row(i)<< " and I:"<<i<<std::endl;
                active_indices.push_back(i);
            }
        }else{
            active_indices.push_back(i);
        }
    }
    if(active_indices.size()<2){
        //std::cout<<"Less than 2 clusters produced, resulting in trivial clustering"<<std::endl;
        return Eigen::VectorXi::Zero(X_rows);
    }
    //make a matrix containing all distances
    unsigned int non_zero_size = active_indices.size();
    RowMatrixXd distances = RowMatrixXd::Zero(X_rows,non_zero_size);
    for(unsigned int i = 0; i <non_zero_size; i++ ){
        distances.col(i)= ((-X).rowwise()+mu.row(active_indices[i])).rowwise().norm();
    }
    //now find the min indices
    Eigen::VectorXi mu_ind =  Eigen::VectorXi::Zero(X_rows);
    Eigen::VectorXd::Index min_index;
    for(unsigned int j = 0; j<X_rows; j++){
        (distances.row(j)).minCoeff(&min_index);
        mu_ind(j)= active_indices[min_index];
                
    //std::cout<<"Chosen indice: "<< mu_ind(j)<<"from vector: "<<distances.row(j)<<std::endl;
    }
    
//    //allocate all points to the first active cluster.
//    Eigen::VectorXi mu_ind =  Eigen::VectorXi::Ones(X_rows)*active_indices[0];
//    
//    //dists starts as the distances to the first cluster
//    Eigen::VectorXd dists =  ((-X).rowwise()+mu.row(active_indices[0])).rowwise().norm();
//    for(unsigned int i = 1; i <non_zero_size; i++ ){
//        //temp_dists is the distance to the next cluster
//        Eigen::VectorXd temp_dists =  ((-X).rowwise()+mu.row(active_indices[i])).rowwise().norm();
//        for(unsigned int j = 0; j<X_rows; j++){
//            if(temp_dists(j)<dists(j)){
//                mu_ind(j)=active_indices[i];
//                dists(j)=temp_dists(j);
//            }
//        }
//    }
    //;
    return mu_ind;
}

const RowMatrixXd Clusterer::reevaluate_centers(const RowMatrixXd& X, const Eigen::VectorXi inds, const unsigned int k, const double reg ){
    
    const unsigned int X_rows = X.rows();
    const unsigned int X_cols = X.cols();
    RowMatrixXd mu_p = RowMatrixXd::Zero(k,X_cols);
    RowMatrixXd mu_n = RowMatrixXd::Zero(k,X_cols);
    // calculate the center of each cluster
    RowMatrixXd centers = RowMatrixXd::Zero(k,X_cols);
    Eigen::VectorXi nums = Eigen::VectorXi::Zero(k);

    for(unsigned int i = 0; i<X_rows; i++){
        centers.row(inds(i))+=X.row(i);
        nums(inds(i))+=1;        
    }
    for(unsigned int i = 0; i< k; i++){
        mu_p.row(i)= ((centers.row(i).array()-reg/2)/nums(i)).matrix();
        mu_n.row(i)= ((centers.row(i).array()+reg/2)/nums(i)).matrix();
    }
    RowMatrixXd mu_new = mu_n.cwiseMin(mu_p.cwiseMax(0));
    return mu_new;
    //now do the three essential cases
    //#pragma omp parallel for shared(mu_p,mu_n)
//    for(unsigned int i = 0; i< k; i++){
//        for(unsigned int j = 0; j< X_cols; j++){
//            if(mu_p(i,j)<0){
//                mu_p(i,j)=mu_n(i,j);
//                if(mu_p(i,j)>0){
//                    mu_p(i,j)=0;
//                }
//            }
//        }
//    }
//    return mu_p;
}
unsigned int Clusterer::element_from_vector(Eigen::VectorXd elements){
    //create a map to put all values in

    unsigned int e_size = elements.size();
    std::map<double, unsigned int> map;
    double cum = 0;
    for(unsigned int i=0; i<e_size; i++){
        if(elements(i)>0){
            cum+=elements(i);
            map[cum]= i;
        }
    }
    //get an index
    double x = ((double) rand() / (RAND_MAX))*cum;
    unsigned int ret =0;
    if(x<cum){
        std::map<double, unsigned int>::iterator iter = map.upper_bound(x);
        ret = iter->second;
    } else{
        ret = map[x];
    }
    
    
    
    return ret;
}

RowMatrixXd Clusterer::initialize_mu(const RowMatrixXd& X, const unsigned int k){
    const unsigned int X_cols = X.cols();
    const unsigned int X_rows = X.rows();
    RowMatrixXd mu =  RowMatrixXd::Zero(k,X_cols);
    
    //pick a random starting point
    int randi = std::rand() % X_rows;
    mu.row(0)=X.row(randi);
    //compute distances
    Eigen::VectorXd dists =  ((-X).rowwise()+mu.row(0)).rowwise().squaredNorm();
    //std::cout<<"element from vector 1: dists"<<dists <<"X"<<X<<"mu"<<mu.row(0)<<std::endl;
    unsigned int element_index = element_from_vector(dists);
    mu.row(1)=X.row(element_index);
    for(unsigned int i = 2; i<k; i++){     
        //get the new minimal distances
        Eigen::VectorXd temp_dists =  ((-X).rowwise()+mu.row(i-1)).rowwise().squaredNorm();
        for(unsigned int j = 0; j<X_rows; j++){
            if(temp_dists(j)<dists(j)){
                dists(j)=temp_dists(j);
            }
        }
        //std::cout<<"element from vector "<<i<<": dists"<<dists <<"X"<<X<<"mu"<<mu.row(0)<<std::endl;
        element_index = element_from_vector(dists);
        mu.row(i)=X.row(element_index);
    }
    return mu;
}

RowMatrixXd Clusterer::m_rescale(const RowMatrixXd& Xin){
    //subtract the median
    //scale the matrix to be central
    RowMatrixXd X = Xin;
    //Eigen::VectorXd median = Eigen::VectorXd::Zero(x_rows);
    unsigned int x_cols = X.cols();
    unsigned int upper = (x_cols-1)*0.99;
    unsigned int lower = (x_cols-1)*0.01;
    for(unsigned int i = 0; i<x_cols; i++){
        Eigen::VectorXd tempCol = X.col(i);
        std::sort(tempCol.data(),tempCol.data()+tempCol.size());
        double scaling = 1;
        if(tempCol(upper)-tempCol(lower)!=0.0){
            scaling = tempCol(upper)-tempCol(lower);
        }
        X.col(i)=(X.col(i).array()/scaling).matrix();
    } 
    X = X.rowwise()-X.colwise().mean();
    //get the mean absolute deviation
//    Eigen::RowVectorXd abs_dev = X.cwiseAbs().colwise().sum();
//    //Eigen::RowVectorXd mean = Xin.colwise().mean();
//    //Eigen::RowVectorXd st_dev = ((Xin.rowwise() - mean).array().square().colwise().sum() / (Xin.rows() - 1)).sqrt();
//    //Check if there are degenerate dimensions
//    for(unsigned int i = 0; i <abs_dev.size(); i++){
//        if(abs_dev(i)==0.0){
//            abs_dev(i) = 1;
//        }
//    }
//    // get rid of the mean abs_dev
//    X = X.array().rowwise()/abs_dev.array();
    
//    if((st_dev.array()==0.0).any()){
//        std::cout<<"Input matrix contains degenrate dimensions with no variance: Aborting"<<std::endl;
//        RowMatrixXd ret;
//        return ret;
//    }
    //RowMatrixXd X = (Xin.rowwise() - mean).array().rowwise() / st_dev.array();
    
    //subtract the median
    //unsigned int x_cols = X.cols();

    return X;
}
double Clusterer::cluster_norm(const RowMatrixXd& X, const RowMatrixXd& centers, const Eigen::VectorXi ind, const double reg){
    const unsigned int rows = X.rows();
    double ssq = 0;
    for(unsigned int i = 0; i<rows; i++){
        ssq+= (X.row(i)-centers.row(ind(i))).squaredNorm();
    }
    for(unsigned int i = 0; i<centers.cols(); i++){
        for(unsigned int j = 0; j<centers.rows(); j++){
            ssq+=centers(j,i)*reg;
        }
    }
    //std::cout<<"The ssq is: "<<ssq<<std::endl; 
    return ssq;
}

const Return_values Clusterer::find_centers(const RowMatrixXd& Xin, const unsigned int k, const double reg, const bool no_zero) {
    //parse the block of memory to an eigen matrix
    const unsigned int rows = Xin.rows();
    //const unsigned int cols = X.cols();
    //scale the matrix to be central
    RowMatrixXd X = Xin;//m_rescale(Xin);
    
    //initialize mu
    
    RowMatrixXd mu = initialize_mu(X, k);
    
    
    Eigen::VectorXi mu_ind =  Eigen::VectorXi::Ones(rows);
    Eigen::VectorXi mu_ind_old =  Eigen::VectorXi::Zero(rows);
    unsigned int count_limit = 1000;
    for(unsigned int i = 0; i<count_limit; i++){
        //std::cout << "Iteration: " <<i<< std::endl;
        mu_ind=allocate_clusters(X,mu);
        if((mu_ind_old-mu_ind).isZero(0)){
            break;
        }
        mu_ind_old=mu_ind;
        mu = reevaluate_centers(X,mu_ind,k,reg);
        
        
    }
    Return_values ret;
    ret.indexes = mu_ind;
    ret.centers = mu;
    ret.norm = cluster_norm(X,mu,mu_ind, reg);
    //if no zero is true, continue to deallocate all points from the zero clusters.
    if(no_zero){
        //std::cout << "Removing zero clusters"<< std::endl;
        for(unsigned int i = 0; i<count_limit; i++){
            //std::cout << "Iteration: " <<i<< std::endl;
            mu_ind=allocate_clusters(X,mu, no_zero);
            if((mu_ind_old-mu_ind).isZero(0)){
                break;
            }
            mu_ind_old=mu_ind;
            mu = reevaluate_centers(X,mu_ind,k,reg);
        }
    }
    ret.indexes_no_zero = mu_ind;
    ret.centers_no_zero = mu;
    ret.norm_no_zero = cluster_norm(X,mu,mu_ind,reg);
    
    //std::tuple<Eigen::VectorXi,RowMatrixXd> tuple;
    
    return ret;
}

const double Clusterer::cluster_distance(const Eigen::VectorXi c1, const Eigen::VectorXi c2, const unsigned int k){
    //first get the cluster distribution
    const unsigned int intSize = c1.size();
    double size = (double) intSize;
    Eigen::VectorXd population1 = Eigen::VectorXd::Zero(k);
    Eigen::VectorXd population2 = Eigen::VectorXd::Zero(k);
    //unsigned int non_zero = 0;
    
    for(unsigned int i = 0; i < intSize; i ++){
        population1(c1(i))+=1.0/size;
        population2(c2(i))+=1.0/size;
    }
    //now calculate the false positives
    //std::cout<<"The populations are: "<< population1 << "and "<< population2 <<std::endl;
    const double false1 = (population1.array()*(population1.array()-1.0/size)*(size/(size-1))).sum()*(population2.array()*(population2.array()-1.0/size)*(size/(size-1))).sum();
    const double false2 = (population1.array()*(1-(population1.array()-1.0/size)*(size/(size-1)))).sum()*(population2.array()*(1-(population2.array()-1.0/size)*(size/(size-1)))).sum();
    
    double stability = 0;
    //always test 10000 indices!
    const unsigned int indice_num = 10000;
    
    for(unsigned int num = 0; num<indice_num; num++){
        unsigned int i = 0;
        unsigned int j = 0;
        while(i==j){
            Eigen::Matrix<unsigned int, 2,1> rand_inds = Eigen::Matrix<unsigned int, 2,1>::Random();
            i = rand_inds(0)%intSize;
            j = rand_inds(1)%intSize;
        }
        if(c1(i)==c1(j)&& c2(i)==c2(j)){
            stability+= 1 - false1;
        } else if((c1(i)!=c1(j)&& c2(i)!=c2(j))){
            stability+=1 - false2;
        }
            
    }
    
    return 1-(stability/indice_num);
}

//generates bootstrapped data sets by resampling an old data set
const RowMatrixXd Clusterer::bootstrap_data(const RowMatrixXd& X, const unsigned int bootstrapSamples){
    const unsigned int bootstrapSize = bootstrapSamples;
    const unsigned int brows = bootstrapSize;
    const unsigned int xrows = X.rows();
    RowMatrixXd b_sample = RowMatrixXd::Zero(brows,X.cols());
    
    const Eigen::Matrix<unsigned int, -1,1> rand_inds = Eigen::Matrix<unsigned int, -1,1>::Random(brows);
    for(unsigned int i = 0; i<brows; i++){
        unsigned int rand_num = rand_inds(i)%xrows;

        b_sample.row(i)=X.row(rand_num);
    }
    return b_sample;
}

unsigned int Clusterer::n_used_clusters(const unsigned int k, const Eigen::VectorXi inds){
    Eigen::VectorXi population = Eigen::VectorXi::Zero(k);
    //unsigned int non_zero = 0;
    unsigned int len = inds.size();
    
    for(unsigned int i = 0; i < len; i ++){
        population(inds(i))=1;
    }
    return population.sum();
}

const Optimization_values Clusterer::optimize_param(const RowMatrixXd& Xin, const Eigen::Matrix<unsigned int, -1,1> k, const Eigen::VectorXd reg, const unsigned int iterations, const unsigned int bootstrapSamples){
    //create the matrix holding the distances, k in the rows, reg in the columns
    RowMatrixXd X = Xin;//m_rescale(Xin);
    const unsigned int k_size = k.size();
    const unsigned int reg_size = reg.size();
    RowMatrixXd distances = RowMatrixXd::Zero(k_size,reg_size);
    RowMatrixXd distances_no_zero = RowMatrixXd::Zero(k_size,reg_size);
    RowMatrixXd found_cluster = RowMatrixXd::Zero(k_size,reg_size);
    RowMatrixXd found_cluster_no_zero = RowMatrixXd::Zero(k_size,reg_size);
    Eigen::VectorXi ind1;
    Eigen::VectorXi ind2;
    //allocate with no_zero
    Eigen::VectorXi ind1_no_zero;
    Eigen::VectorXi ind2_no_zero;
    for(unsigned int i = 0; i<iterations; i++){
        //std::cout<<"Overall iteration: "<< i << std::endl;
        for(unsigned int j = 0; j<k_size; j++){
            for(unsigned int l = 0; l<reg_size; l++){
                //create three bootstrap samples
                const RowMatrixXd b1 = bootstrap_data(X,bootstrapSamples);
                const RowMatrixXd b2 = bootstrap_data(X,bootstrapSamples);
                //const RowMatrixXd b3 = bootstrap_data(X,bootstrapSamples);
                //create two partitionings
                //std::cout<<"Finding centers: "<< std::endl;
                //try to remove zeros by default
                bool remove_zeros = true;
                const Return_values ret1 = find_centers(b1,k(j), reg(l),remove_zeros);
                const Return_values ret2 = find_centers(b2,k(j), reg(l),remove_zeros);
                //allocate X using the new mus
                ind1 = allocate_clusters(X, ret1.centers );
                ind2 = allocate_clusters(X, ret2.centers );
                //allocate with no_zero
                ind1_no_zero = allocate_clusters(X, ret1.centers_no_zero, true);
                ind2_no_zero = allocate_clusters(X, ret2.centers_no_zero, true);
                //std::cout<<"calculating distances: "<< std::endl;
                distances(j,l)+=cluster_distance(ind1,ind2,k(j));
                distances_no_zero(j,l)+=cluster_distance(ind1_no_zero,ind2_no_zero,k(j));
                //std::cout<<"countng non-zero: "<< std::endl;
                found_cluster(j,l)+=n_used_clusters(k(j),ret1.indexes);
                found_cluster(j,l)+=n_used_clusters(k(j),ret2.indexes);
                found_cluster_no_zero(j,l)+=n_used_clusters(k(j),ret1.indexes_no_zero);
                found_cluster_no_zero(j,l)+=n_used_clusters(k(j),ret2.indexes_no_zero);
            }
        }
    }
    Optimization_values ret;
    ret.distances= distances/iterations;
    ret.distances_no_zero= distances_no_zero/iterations;
    ret.found_cluster = found_cluster/(iterations*2);
    ret.found_cluster_no_zero = found_cluster_no_zero/(iterations*2);
    ret.indexes1 = ind1;
    ret.indexes2 = ind2;
    ret.indexes1_no_zero = ind1_no_zero;
    ret.indexes2_no_zero = ind2_no_zero;
    return ret;
    
}

Clusterer::~Clusterer() {
}

