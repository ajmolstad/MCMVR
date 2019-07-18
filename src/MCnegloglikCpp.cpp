#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
double MCnegloglikCpp(arma::mat beta, arma::mat x, arma::mat y, double tau, arma::mat d){
	double n = x.n_rows;
	mat temp = y - x*beta; 
	mat cov = (beta.t()%d)*beta; 
	cov.diag() += tau; 
	cov = solve(cov, temp.t()); 
	//double keep = 0.0; 
	//for(int i=0; i<n; ++i){
	//	keep += (temp(span(i,i), span::all)*cov*temp(span(i,i),span::all).t()).eval()(0,0); 
	//}
	double out = pow(n, -1)*accu(temp%cov.t()); 
	return out; 
}

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
arma::mat MCgradbetaCpp(arma::mat beta, arma::mat x, arma::mat y, arma::mat xtx, arma::mat xty, double tau, arma::mat d){
	double n = x.n_rows;
	mat cov = (beta.t()%d)*beta; 
	cov.diag() += tau; 
	cov = inv(cov); 
	mat temp = y - x*beta; 
	mat T1 = - (d.t()%beta)*cov.t()*temp.t()*temp - xty + xtx*beta;
	mat out = 2*pow(n,-1)*(T1*cov);
	return out; 
}
