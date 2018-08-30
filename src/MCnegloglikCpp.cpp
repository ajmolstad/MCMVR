#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
double MCnegloglikCpp(arma::mat beta, arma::mat z, arma::mat y, double tau){
	double n = z.n_rows;
	mat temp = y - z*beta; 
	mat cov = beta.t()*beta; 
	cov.diag() += tau; 
	cov = inv(cov); 
	double keep = 0.0; 
	for(int i=0; i<n; ++i){
		keep += (temp(span(i,i), span::all)*cov*temp(span(i,i),span::all).t()).eval()(0,0); 
	}

	double out = pow(n, -1)*keep; 
	return out; 
}

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
arma::mat MCgradbetaCpp(arma::mat beta, arma::mat x, arma::mat y, arma::mat xtx, arma::mat xty, double tau){
	double n = x.n_rows;
	mat cov = beta.t()*beta; 
	cov.diag() += tau; 
	cov = inv(cov); 
	mat temp = y - x*beta; 
	mat T1 = - beta*cov.t()*temp.t()*temp - xty + xtx*beta;
	mat out = 2*pow(n,-1)*(T1*cov);
	return out; 
}
