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
    //Rcout<<"the norm is: "<< partitioning.norm<<std::endl;
  return ret;
}

// [[Rcpp::export]]
List grid_search(NumericMatrix X, IntegerVector k, NumericVector reg, const unsigned int iterations, const unsigned int bootstrapSamples, const unsigned long seed_off_set = 0) {
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

    ret["d"] = distances;
    ret["z"] = distances_no_zero;
    ret["n"] = found_cluster;
    ret["m"] = found_cluster_no_zero;
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
const Eigen::VectorXi Clusterer::allocate_to_model_based(const RowMatrixXd& X, const RowMatrixXd& mu, const RowMatrixXd& sigma, const Eigen::VectorXd& pi){
      //calculate the fkj matrix
  Eigen::MatrixXd wfkj = Eigen::MatrixXd::Zero(this->k_, this->j_);
  for(unsigned int i = 0; i<this->k_; i++ ){
    for(unsigned int l = 0; l< this->j_; l++ ){
      Eigen::ArrayXd diff = X.row(l).array()-mu.row(i).array();
      Eigen::ArrayXd inv = 1.0/sigma.row(i).array();
      wfkj(i,l)=pi(i)*inv.prod()*std::exp(-0.5*(diff*diff*inv*inv).sum());
    }
  }
    //ignore one zero cluster, kill the remaining
  bool zero_found = false;
  for(unsigned int i=0;i<this->k_; i++){
      if(mu.row(i).isZero(0)){
          if(zero_found){
              wfkj.row(i).setZero();
          } else {
              zero_found=true;
          }
        
      }
  }
  Eigen::MatrixXd tau = Eigen::MatrixXd::Zero(this->k_, this->j_);
  for(unsigned int l = 0; l< this->j_; l++ ){
      double sum = wfkj.col(l).array().sum();
      if(sum==0){
          tau.col(l).setZero();
      }else{
          tau.col(l)=wfkj.col(l)/sum;
      }    
  }
    
  if(!((tau.array() == tau.array())).all()){
    std::cout<< "Nan values in tau_ found"<<std::endl;
  }
  
  return index_from_tau(tau);
}

const Eigen::VectorXi Clusterer::allocate_clusters(const RowMatrixXd& X, const RowMatrixXd& mu, const bool no_zero ){
    const unsigned int mu_rows = mu.rows();
    const unsigned int X_rows = X.rows();

    //make a vector to iterate over
    std::vector<unsigned int> active_indices;
    for(unsigned int i = 0; i< mu_rows; i++){
        if(no_zero){
            if(!mu.row(i).isZero(0)){
                //Rcout<<"muRow: "<< mu.row(i)<< " and I:"<<i<<std::endl;
                active_indices.push_back(i);
            }
        }else{
            active_indices.push_back(i);
        }
    }
    
    
    if(active_indices.size()==0){
        active_indices.push_back(0);
    }
    //allocate all points to the first active cluster.
    Eigen::VectorXi mu_ind =  Eigen::VectorXi::Ones(X_rows)*active_indices[0];
    unsigned int non_zero_size = active_indices.size();
    //dists starts as the distances to the first cluster
    Eigen::VectorXd dists =  ((-X).rowwise()+mu.row(active_indices[0])).rowwise().norm();
    for(unsigned int i = 1; i <non_zero_size; i++ ){
        //temp_dists is the distance to the next cluster
        Eigen::VectorXd temp_dists =  ((-X).rowwise()+mu.row(active_indices[i])).rowwise().norm();
        for(unsigned int j = 0; j<X_rows; j++){
            if(temp_dists(j)<dists(j)){
                mu_ind(j)=active_indices[i];
                dists(j)=temp_dists(j);
            }
        }
    }
    //Rcout<<"mu ind: "<< mu_ind<<std::endl;
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
    //now do the three essential cases
    //#pragma omp parallel for shared(mu_p,mu_n)
    for(unsigned int i = 0; i< k; i++){
        for(unsigned int j = 0; j< X_cols; j++){
            if(mu_p(i,j)<0){
                mu_p(i,j)=mu_n(i,j);
                if(mu_p(i,j)>0){
                    mu_p(i,j)=0;
                }
            }
        }
    }
    return mu_p;
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
    //Rcout<<"element from vector 1: dists"<<dists <<"X"<<X<<"mu"<<mu.row(0)<<std::endl;
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
        //Rcout<<"element from vector "<<i<<": dists"<<dists <<"X"<<X<<"mu"<<mu.row(0)<<std::endl;
        element_index = element_from_vector(dists);
        mu.row(i)=X.row(element_index);
    }
    return mu;
}

RowMatrixXd Clusterer::m_rescale(const RowMatrixXd& Xin){

    RowMatrixXd X = Xin;

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


    return X;
}
double Clusterer::cluster_norm(const RowMatrixXd& X, const double reg){
  double norm = 0;
  for(unsigned int i = 0; i<this->k_; i++ ){
    for(unsigned int l = 0; l< this->j_; l++ ){
      Eigen::ArrayXd diff = X.row(l).array()-this->mu_.row(i).array();
      Eigen::ArrayXd inv = 1.0/this->sigma_.row(i).array();
      if(this->pi_(i)>0){
        norm+=this->tau_(i,l)*(std::log(this->pi_(i))+std::log(inv.prod())-0.5*(diff*diff*inv*inv).sum());
      }
    }
  }
  norm-=this->mu_.lpNorm<1>()*reg;

    return norm;
}

void Clusterer::initialize_members(const RowMatrixXd& X, const RowMatrixXd& mu){
  const unsigned int mu_rows = mu.rows();
  this->k_ = mu_rows;
  this->j_ = X.rows();
  this->p_ = X.cols();
  this->tau_=Eigen::MatrixXd::Zero(this->k_,this->j_);
  this->mu_= mu;
  this->sigma_=RowMatrixXd::Ones(this->k_,this->p_);
  this->pi_ = Eigen::VectorXd::Ones(this->k_)/this->k_;
}

void Clusterer::e_step(const RowMatrixXd& X){
  //calculate the fkj matrix
  Eigen::MatrixXd wfkj = Eigen::MatrixXd::Zero(this->k_, this->j_);
  //for each cluster, get the weigthed probability
  //std::cout<<"Some random outputs, mu , pi sigma: "<<mu_<<pi_<<sigma_<<std::endl;
  //std::cout<<"k and j: "<<k_<<j_<<std::endl;
  for(unsigned int i = 0; i<this->k_; i++ ){
    for(unsigned int l = 0; l< this->j_; l++ ){
      Eigen::ArrayXd diff = X.row(l).array()-this->mu_.row(i).array();
      Eigen::ArrayXd inv = 1.0/this->sigma_.row(i).array();
      //std::cout<<"diff and inv: "<<diff<<inv<<std::endl;
      //std::cout<<"assignment: "<<this->pi_(i)*inv.prod()*std::exp(-0.5*(diff*diff*inv).sum())<<std::endl;
      wfkj(i,l)=this->pi_(i)*inv.prod()*std::exp(-0.5*(diff*diff*inv*inv).sum());
    }
  }
    //ignore one zero cluster, kill the remaining
  bool zero_found = false;
  for(unsigned int i=0;i<this->k_; i++){
      if(this->mu_.row(i).isZero(0)){
          if(zero_found){
              wfkj.row(i).setZero();
          } else {
              zero_found=true;
          }
        
      }
  }
  for(unsigned int l = 0; l< this->j_; l++ ){
      double sum = wfkj.col(l).array().sum();
      if(sum==0){
          this->tau_.col(l).setZero();
      }else{
          this->tau_.col(l)=wfkj.col(l)/sum;
      }
        
  }

  if(!((tau_.array() == tau_.array())).all()){
    std::cout<< "Nan values in tau_ found"<<std::endl;
  }
}

int Clusterer::m_step(const RowMatrixXd& X, const double pen_term){
  //update pi
  RowMatrixXd old_pi=this->pi_;
  for(unsigned int i = 0; i<this->k_; i++){
    this->pi_(i)=this->tau_.row(i).array().sum()/this->j_;
  }
  
  //Axels analytic solution!
  //loop through k,p to get some important statistics
  RowMatrixXd no_x = RowMatrixXd::Zero(this->k_, this->p_);
  RowMatrixXd lin_x = RowMatrixXd::Zero(this->k_, this->p_);
  RowMatrixXd square_x = RowMatrixXd::Zero(this->k_, this->p_);
  for(unsigned int m = 0; m<this->p_; m++){
        for(unsigned int i = 0; i<this->k_; i++ ){
      
          no_x(i,m)=this->tau_.row(i).array().sum();
          lin_x(i,m)=(this->tau_.row(i).array().transpose()*X.col(m).array()).sum();
          square_x(i,m)=(this->tau_.row(i).array().transpose()*X.col(m).array().square()).sum();        
      }     
  }
  //std::cout<<"X: "<<X<<"tau: "<<this->tau_<< "noX: "<<noX<<"linX: "<<linX<<"squareX: "<<squareX<<std::endl;
  RowMatrixXd old_mu=this->mu_;
  RowMatrixXd old_sigma=this->sigma_;
  //now calculate the mus
  for(unsigned int m = 0; m<this->p_; m++){
    for(unsigned int i = 0; i<this->k_; i++ ){
        double mu = 0;
        double sigma = 0.001;
        if(no_x(i,m)>std::exp(-20)){
            double a = lin_x(i,m)/pen_term;
             

            if(lin_x(i,m)>0){
                double q_pos = square_x(i,m)/no_x(i,m)-a;
                double b_pos = no_x(i,m)/pen_term;
                double p_pos = b_pos-2*lin_x(i,m)/no_x(i,m);
                mu = -p_pos/2+std::sqrt(std::pow(p_pos,2)/4-q_pos);
                if(mu<0 || mu!=mu){
                    mu=0;
                    sigma = std::max(std::sqrt(square_x(i,m)/no_x(i,m)),0.001);
                }else{
                    if(a-mu*b_pos>0){
                        sigma = std::max(std::sqrt(a-mu*b_pos),0.001);  
                    }
                }


            } else if (lin_x(i,m)<0){
                double q_neg = square_x(i,m)/no_x(i,m)+a;
                double b_neg = -no_x(i,m)/pen_term;
                double p_neg = b_neg-2*lin_x(i,m)/no_x(i,m);
                mu = -p_neg/2-std::sqrt(std::pow(p_neg,2)/4-q_neg);
                if(mu>0 || mu!=mu){
                    mu=0;
                    sigma = std::max(std::sqrt(square_x(i,m)/no_x(i,m)),0.001);
                }else{
                    if(-a-mu*b_neg>0){
                        sigma = std::max(std::sqrt(-a-mu*b_neg),0.001); 
                    }
                }

            }
        }
        this->mu_(i,m)=mu;
        this->sigma_(i,m)=sigma;
        
    }
  }
  //check for convergence
  double diffMu = (this->mu_-old_mu).squaredNorm()/(this->mu_.squaredNorm()+0.001);
  double diffSigma = (this->sigma_-old_sigma).squaredNorm()/(this->sigma_.squaredNorm()+0.001);
  double diffPi = (this->pi_-old_pi).squaredNorm()/(this->pi_.squaredNorm()+0.001);

  std::cout<< "The diff Values for Mu, Sigma and Pi are: "<< diffMu <<", " <<diffSigma<<", " << diffPi<< std::endl;

  if (diffMu<std::exp(-10)){// && diffPi<0.001 && diffSigma<0.001){
    return 1;
  } else {
    return 0;
  }
}

const Eigen::VectorXi Clusterer::index_from_tau(Eigen::MatrixXd& tau){
  Eigen::VectorXi mu_ind = Eigen::VectorXi::Zero(this->j_);
  Eigen::MatrixXf::Index max_index;
  //std::cout<<"Tau is: "<<tau<<std::endl;
  for(unsigned int l=0; l<this->j_;l++){
      double max = 0;
     for(unsigned int i=0; i<this->k_;i++){ 
         if(tau(i,l)>max){
             mu_ind(l)=i;
             max = tau(i,l);
         }
        
     }
  }
  //std::cout<<"Mu ind is: "<<mu_ind<<std::endl;
  return mu_ind;
}

const Return_values Clusterer::find_centers(const RowMatrixXd& Xin, const unsigned int k, const double reg, const bool no_zero) {
    //parse the block of memory to an eigen matrix
    const unsigned int rows = Xin.rows();
    //const unsigned int cols = X.cols();
    //scale the matrix to be central
    RowMatrixXd X = Xin;//
    //RowMatrixXd X = m_rescale(Xin);

    //initialize mu

    RowMatrixXd mu = initialize_mu(X, k);


    Eigen::VectorXi mu_ind =  Eigen::VectorXi::Ones(rows);
    Eigen::VectorXi mu_ind_old =  Eigen::VectorXi::Zero(rows);
    unsigned int count_limit = 1000;
    //modified to normal k-means
    for(unsigned int i = 0; i<count_limit; i++){
        //Rcout << "Iteration: " <<i<< std::endl;
        mu_ind=allocate_clusters(X,mu);
        if((mu_ind_old-mu_ind).isZero(0)){
            break;
        }
        mu_ind_old=mu_ind;
        mu = reevaluate_centers(X,mu_ind,k,0);


    }
    //extend to model based clustering
    //initialize pi tao mu and sigma
    this->initialize_members(X,mu);
    int conv =0;
    int count=0;
    std::cout<<"Starting PMC"<<std::endl;
    //double old_norm = 0;
    //double norm = 0;
    while(conv==0 && count < 100){
      this->e_step(X);
      conv=this->m_step(X,reg);
      //old_norm=norm;
      //norm= this->cluster_norm(X,reg);
      //std::cout<<"The norm diff is: "<< norm-old_norm << std::endl;
      count++;
    }
    //index from tau
    mu_ind=this->index_from_tau(this->tau_);
    Return_values ret;
    ret.indexes = mu_ind;
    ret.centers = this->mu_;
    ret.norm = cluster_norm(X,reg);
    ret.indexes_no_zero = mu_ind;
    ret.centers_no_zero = this->mu_;
    ret.norm_no_zero = ret.norm;

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
    //Rcout<<"The populations are: "<< population1 << "and "<< population2 <<std::endl;
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
            stability+= 1;
        } else if((c1(i)!=c1(j)&& c2(i)!=c2(j))){
            stability+=1;
        }

    }

    return 1-((stability/indice_num)-(false1+false2))/(1-(false1+false2));
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

    for(unsigned int i = 0; i<iterations; i++){
        //Rcout<<"Overall iteration: "<< i << std::endl;
        for(unsigned int j = 0; j<k_size; j++){
            for(unsigned int l = 0; l<reg_size; l++){
                //create three bootstrap samples
                const RowMatrixXd b1 = bootstrap_data(X,bootstrapSamples);
                const RowMatrixXd b2 = bootstrap_data(X,bootstrapSamples);
                const RowMatrixXd b3 = bootstrap_data(X,bootstrapSamples);
                //create two partitionings
                //Rcout<<"Finding centers: "<< std::endl;
                bool remove_zeros = false;
                const Return_values ret1 = find_centers(b1,k(j), reg(l),remove_zeros);
                const Eigen::VectorXi ind1 = allocate_to_model_based(b3, this->mu_, this->sigma_, this->pi_);
                const Return_values ret2 = find_centers(b2,k(j), reg(l),remove_zeros);
                const Eigen::VectorXi ind2 = allocate_to_model_based(b3, this->mu_, this->sigma_, this->pi_);

                distances(j,l)+=cluster_distance(ind1,ind2,k(j));

                found_cluster(j,l)+=n_used_clusters(k(j),ret1.indexes);
                found_cluster(j,l)+=n_used_clusters(k(j),ret2.indexes);

            }
        }
    }
    Optimization_values ret;
    ret.distances= distances/iterations;
    ret.distances_no_zero= distances/iterations;
    ret.found_cluster = found_cluster/(iterations*2);
    ret.found_cluster_no_zero = found_cluster/(iterations*2);
    return ret;

}

Clusterer::~Clusterer() {
}

