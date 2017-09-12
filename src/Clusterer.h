/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   Clusterer.h
 * Author: theorell
 *
 * Created on May 3, 2017, 8:40 AM
 */

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Rcpp.h>
#ifndef CLUSTERER_H
#define CLUSTERER_H
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;

//This struct is returned from the method find_centers.
struct Return_values {
    Eigen::VectorXi indexes;
    Eigen::VectorXi indexes_no_zero;
    RowMatrixXd centers;
    RowMatrixXd centers_no_zero;
    double norm;
    double norm_no_zero;
};

//this struct is returned from the mehod optimize parameters.
struct Optimization_values {
    RowMatrixXd found_cluster;
    RowMatrixXd found_cluster_no_zero;
    RowMatrixXd distances;
    RowMatrixXd distances_no_zero;
};

Rcpp::NumericMatrix eigen_to_numeric(RowMatrixXd X);
Rcpp::List allocate_points(Rcpp::NumericMatrix X, Rcpp::NumericMatrix mu,const bool no_zero);

class Clusterer {
public:
    Clusterer();
    Clusterer(const Clusterer& orig);
    const Return_values find_centers(const RowMatrixXd& X, const unsigned int k, const double reg, const bool no_zero=false);
    const Optimization_values optimize_param(const RowMatrixXd& X, const Eigen::Matrix<unsigned int, -1,1> k, const Eigen::VectorXd reg, const unsigned int iterations, const unsigned int bootstrapSamples);
    const Eigen::VectorXi allocate_clusters(const RowMatrixXd& X, const RowMatrixXd& mu, const bool no_zero=false);
    void reseed(const unsigned long seed_off_set){
        srand((unsigned int) time(0)+seed_off_set);;
    }
    virtual ~Clusterer();
    const Eigen::VectorXi allocate_to_model_based(const RowMatrixXd& X, const RowMatrixXd& mu, const RowMatrixXd& sigma, const Eigen::VectorXd& pi);
    friend class Tester;
    friend Rcpp::List allocate_points(Rcpp::NumericMatrix X, Rcpp::NumericMatrix mu,const bool no_zero );


private:
    int m_step(const RowMatrixXd& X, const double regVec);
    void e_step(const RowMatrixXd& X);
    void initialize_members(const RowMatrixXd& X, const RowMatrixXd& mu);
    double cluster_norm(const RowMatrixXd& X, const double reg);
    RowMatrixXd m_rescale(const RowMatrixXd& Xin);
    const RowMatrixXd  reevaluate_centers(const RowMatrixXd& X, const Eigen::VectorXi inds, const unsigned int k, const double reg );
    const Eigen::VectorXi index_from_tau(Eigen::MatrixXd& tau);
    RowMatrixXd initialize_mu(const RowMatrixXd& X, const unsigned int k);
    unsigned int element_from_vector(Eigen::VectorXd elements);
    const double cluster_distance(const Eigen::VectorXi c1, const Eigen::VectorXi c2, const unsigned int k);
    const RowMatrixXd bootstrap_data(const RowMatrixXd& X, const unsigned int bootstrapSamples);
    unsigned int n_used_clusters(const unsigned int k, const Eigen::VectorXi inds);

    Eigen::MatrixXd tau_;
    RowMatrixXd mu_;
    RowMatrixXd sigma_;
    Eigen::VectorXd pi_;
    unsigned int p_;
    unsigned int k_;
    unsigned int j_;


};

#endif /* CLUSTERER_H */

