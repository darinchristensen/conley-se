#include <iostream>
#include <math.h>
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// ------------------------------------------


double to_radians_cpp(double degrees){
    return(degrees * 3.141593 / 180);
}


// ------------------------------------------

// Haversine Formula:
double haversine_cpp(double lat1, double long1,
                     double lat2, double long2,
                     std::string unit="km"){
    int radius = 6378;
    double delta_phi = to_radians_cpp(lat2 - lat1);
    double delta_lambda = to_radians_cpp(long2 - long1);
    double phi1 = to_radians_cpp(lat1);
    double phi2 = to_radians_cpp(lat2);
    double term1 = pow(sin(delta_phi / 2), 2);
    double term2 = cos(phi1) * cos(phi2) * pow(sin(delta_lambda/2), 2);
    double the_terms = term1 + term2;
    double delta_sigma = 2 * atan2(sqrt(the_terms), sqrt(1-the_terms));
    double distance = radius * delta_sigma;

    /* if it is anything *but* km it is miles */
    if(unit != "km"){
        return(distance*0.621371);
    }

    return(distance);
}


// ------------------------------------------


// Distance Formula used by SH:
double sh_cpp(double lat1, double long1,
             double lat2, double long2) {
    double distance = pow(pow(111 * (lat1 - lat2), 2) + pow(cos(lat1 * 3.1415927 / 180) * 111 * (long1 - long2), 2), .5);
    return(distance);
}


// ------------------------------------------

// [[Rcpp::export]]
arma::mat DistMat(arma::mat M, double cutoff,
   std::string kernel="bartlett",
   std::string dist_fn="Haversine"){

    long long int nrow = M.n_rows;
    arma::mat dmat(nrow, nrow, fill::zeros);

    for( long long int i = 0; i < nrow; i++ ){
        dmat(i, i) = 1;

        for( long long int j = i+1; j < nrow; j++ ){
            // Distance Function:
            double d;
            if(dist_fn != "Haversine") {
                d = sh_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
            } else {
                d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
            }

            // Kernel:
            int v = d <= cutoff;
            // if(v == 0) continue;

            if(kernel != "bartlett") {
                dmat(i,j) = dmat(j,i) = v;
            } else {
                dmat(i,j) = dmat(j,i) = (1 - d / cutoff) * v;
            }
        }
    }
    // dmat += eye(nrow, nrow);
    return dmat;
}


// ------------------------------------------


// [[Rcpp::export]]
arma::mat XeeXhC(arma::mat M, double cutoff,
   arma::mat X, arma::vec e, int n1, int k,
   std::string kernel="bartlett",
   std::string dist_fn="Haversine"){

    long long int nrow = M.n_rows;
    arma::mat dmat(nrow, nrow, fill::zeros);

    for( long long int i = 0; i < nrow; i++ ){
        dmat(i, i) = 1;

        for( long long int j = i+1; j < nrow; j++ ){
            // Distance Function:
            double d;
            if(dist_fn != "Haversine") {
                d = sh_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
            } else {
                d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
            }

            // Kernel:
            int v = d <= cutoff;
            // if(v == 0) continue;

            if(kernel != "bartlett") {
                dmat(i,j) = dmat(j,i) = v;
            } else {
                dmat(i,j) = dmat(j,i) = (1 - d / cutoff) * v;
            }
        }
    }
    // dmat += eye(nrow, nrow);
    // return dmat;

    arma::mat XeeXh(k, k, fill::zeros);
    for( long long int i = 0; i < nrow; i++ ){
        arma::mat e_mat(1, n1, fill::zeros);
        e_mat.fill(e[i]);

        arma::mat k_mat(k, 1, fill::ones);

        arma::mat d_row(1, n1, fill::ones);
        d_row %= dmat.row(i); d_row %= e.t();

        arma::mat X_row(k, 1, fill::ones);
        X_row %= X.row(i).t();

        XeeXh += (X_row * e_mat % (k_mat * d_row)) * X;
    }
    return XeeXh;
}


// ------------------------------------------

// [[Rcpp::export]]
arma::mat Bal_XeeXhC(arma::mat dmat,
   arma::mat X, arma::vec e, int n1, int k){
    
    long long int nrow = dmat.n_rows;
    
    arma::mat XeeXh(k, k, fill::zeros);
    for( long long int i = 0; i < nrow; i++ ){
        arma::mat e_mat(1, n1, fill::zeros);
        e_mat.fill(e[i]);

        arma::mat k_mat(k, 1, fill::ones);

        arma::mat d_row(1, n1, fill::ones);
        d_row %= dmat.row(i); d_row %= e.t();

        arma::mat X_row(k, 1, fill::ones);
        X_row %= X.row(i).t();

        XeeXh += (X_row * e_mat % (k_mat * d_row)) * X;
    }
    return XeeXh;
}


// ------------------------------------------

// [[Rcpp::export]]
arma::mat XeeXhC_Lg(arma::mat M, double cutoff,
   arma::mat X, arma::vec e, int n1, int k,
   std::string kernel="bartlett",
   std::string dist_fn="Haversine"){

    long long int nrow = M.n_rows;
    arma::mat XeeXh(k, k, fill::zeros);

    for( long long int i = 0; i < nrow; i++ ){
        arma::vec d_row(nrow);

        for( long long int j = 0; j < nrow; j++ ){
            double d;
            // Distance Function:
            if(dist_fn != "Haversine") {
                d = sh_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
            } else {
                d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
            }

            // Kernel:
            int v = d <= cutoff;
            if(kernel != "bartlett") {
                d_row[j] = v;
            } else {
                d_row[j] = (1 - d / cutoff) * v;
            }
        }

        arma::mat e_mat(1, nrow, fill::zeros);
        e_mat.fill(e[i]);

        arma::mat k_mat(k, 1, fill::ones);

        d_row %= e;

        arma::mat X_row(1, k, fill::ones);
        X_row %= X.row(i);

        XeeXh += (X_row * e_mat % (k_mat * d_row.t())) * X;
    }

    return XeeXh;
}


// ------------------------------------------


// [[Rcpp::export]]
arma::mat TimeDist(arma::vec times, double cutoff,
   arma::mat X, arma::vec e, int n1, int k){

    long long int nrow = times.n_elem;
    arma::mat dmat(nrow, nrow, fill::ones);
    // dmat.each_row() %= times.t();

    for( long long int i = 0; i < nrow; i++ ){
        arma::vec t_diff(nrow);
        t_diff = times;

        t_diff -= times[i];
        t_diff = abs(t_diff);

        NumericVector v1(nrow); NumericVector v2(nrow);
        for( long long int j = 0; j < nrow; j++ ) {
            v1[j] = t_diff[j] <= cutoff;
            v2[j] = t_diff[j] != t_diff[i];
            t_diff[j] = v1[j] * v2[j] * (1 - t_diff[j] / (cutoff + 1));
            // t_diff[j] = v1[j] * v2[j] * (1 - t_diff[j]) / (cutoff + 1);
        }

        dmat.row(i) %= t_diff.t();
    }

    arma::mat XeeXh(k, k, fill::zeros);
    for( long long int i = 0; i < nrow; i++ ){
        arma::mat e_mat(1, n1, fill::zeros);
        e_mat.fill(e[i]);

        arma::mat k_mat(k, 1, fill::ones);

        arma::mat d_row(1, n1, fill::ones);
        d_row %= dmat.row(i); d_row %= e.t();

        arma::mat X_row(k, 1, fill::ones);
        X_row %= X.row(i).t();

        XeeXh += (X_row * e_mat % (k_mat * d_row)) * X;
    }

    return XeeXh;
}
