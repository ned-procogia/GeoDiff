#include <cmath>
#include "GeoDiff.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]
using namespace std;
using namespace Rcpp;
using namespace roptim;

// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
// // [[Rcpp::export]]
// arma::vec dnbinom_mu_vec(arma::vec x, double sz, arma::vec mu, int lg){
//   int N = x.n_elem;
//   arma::vec prob(N);
//   //Rcpp::dnbinom_mu(x, sz, mu, lg)
//   for(int i=0; i<N; i++)
//     prob(i) = R::dnbinom_mu(x(i), sz, mu(i), lg);
// 
//   return(prob);
// }

class NBthDE_paranll : public Functor {
public:
  arma::vec y; // what is y?
  arma::mat X; // what is x?
  arma::vec alpha0; // what is alpha0?
  arma::vec alpha; // what is alpha0?
  arma::mat preci1; // what is alpha0?
  double preci2; // what is alpha0?
  double threshold0; // what this a threashold for?


  double operator()(const arma::vec &x) override {
    int n = X.n_cols;
    arma::vec beta = x(arma::span(0,n-1));
    
    
    double r = x(n);
    //Rcout<< "r:" << r << ";   \n";
    if(std::isnan(r)){
      throw 20;
    }
    double threshold = x(n+1);


    arma::vec tmp0 = arma::exp2(X*beta);
    // armadillo vectors: % element-wise multiplication of two objects
    arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
    arma::vec llk = dnbinom_mu_vec(y, r, tmp1, 1);
    //arma::mat pen
    arma::mat pen10 = beta.t()*preci1*beta;
    //Rcout << "pen10(0,0): " << pen10(0,0) << "\n";
    //Rcout << "sum(llk): " << sum(llk) << "\n";
    double pen1 = pen10(0,0)/2.0;
    //+nmh*(1.0/2.0)*pow((threshold-threshold0),2)*preci2
    if (std::isinf(sum(llk))){
      throw 20;
    }
    return(-arma::sum(llk)+pen1+(1.0/2.0)*pow((threshold-threshold0),2)*preci2);
  }



  void Gradient(const arma::vec &x, arma::vec &gr) override {
    int n = X.n_cols;
    int m = y.n_elem;


    gr = arma::zeros<arma::vec>(n+2);

    arma::vec beta = x(arma::span(0,n-1));

    double r = x(n);
    double threshold = x(n+1);

    // can r = 0?
    // 
    arma::vec tmp0 = arma::exp2(X*beta);

    arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
    arma::vec tmp2 = (y/tmp1-1.0)/(1.0+tmp1/r);

    gr(arma::span(0,n-1)) = (-log(2.0)*(tmp2%alpha%tmp0).t()*X+beta.t()*preci1).t();

    arma::vec pLr = -arma::log(1.0+tmp1/r);

    for(int k = 0; k < m; k++){
      for(int j = 0; j < y(k); j++){
        pLr(k) += 1.0/(j+r);
      }
    }

    pLr += -(y-tmp1)/(r+tmp1);

    gr(n) = -arma::sum(pLr);

    arma::mat tmp3 = tmp2.t()*alpha0;
    gr(n+1) = -tmp3(0,0)+(threshold-threshold0)*preci2;

    }

};

// [[Rcpp::export]]
List NBthDE_paraOptfeat(arma::mat& X, //t(object_mat[features_high, ])
               arma::vec y, // also confusingly called X in the layer above // y = model.matrix(~fov, data = annot[high_ROIs, ])
               arma::vec alpha0, // sizefact_BG
               arma::vec alpha, // sizefact
               arma::mat& preci1, // preci1
               double threshold0, //threshold_mean * probenum[features_high]
               double preci2, // preci2, conveniently
               arma::vec& x0, // startpara
               bool calhes) { // sizescalebythreshold
  NBthDE_paranll f;
  f.X=X;
  f.y=y;
  f.alpha0 = alpha0;
  f.alpha = alpha;
  f.preci1=preci1;
  f.threshold0=threshold0;
  f.preci2=preci2;
  //Rcout <<"x0: "<< x0 << "\n";
  int n = X.n_cols;
  arma::vec x = x0;
  // arma::vec beta = x(arma::span(0,n-1));
  // double r = x(n);
  // double threshold = x(n+1);
  // arma::vec tmp0 = arma::exp2(X*beta);
  // arma::vec tmp1 = alpha0*threshold+alpha%tmp0;
  // if (isinf(arma::accu(tmp1))){
  //   Rcout << "accu(tmp1): " << accu(tmp1) << "\n";
  // }
  arma::vec lower = arma::ones<arma::vec>(n+2) * (-100);
  lower(n) = 0.01;
  lower(n+1) = 0.01;
  arma::vec upper = arma::ones<arma::vec>(n+2) * 100;
  upper(n) = 1000;
  upper(n+1) = 1000000;
  Roptim<NBthDE_paranll> opt("L-BFGS-B");
  opt.set_lower(lower);
  opt.set_upper(upper);
  //opt.control.maxit = maxit;
  
  opt.set_hessian(calhes);
  
  // opt.set_hessian(true);
  opt.control.pgtol=1e-3;
  // arma::zeros<arma::vec>(n+2);
  // x(arma::span(0,n-1))=arma::solve(X, arma::log2(y/alpha + 0.001));
  // x(n) = 1;
  // x(n+1)=threshold0;
  
  opt.minimize(f, x);
  
  // arma::mat hes = opt.hessian();
  // double hes_det = arma::log_det(hes);
  // double hes_det1 = arma::log_det(hes(arma::span(1,n-1), arma::span(1,n-1)));
  return List::create(Named("par") = opt.par(),
                      Named("conv") = opt.convergence(),
                      Named("hes") = opt.hessian());

}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List NBthDE_paraOptall(arma::mat& Y,
                   arma::mat& X,
                   arma::vec& alpha0,
                   arma::vec& alpha,
                   arma::mat& preci1,
                   arma::vec& threshold0,
                   double preci2,
                   arma::vec& x0,
                   bool sizescale,
                   bool calhes){

  int n = X.n_cols;
  int m = Y.n_cols;
  Rcout << "number of rows: "<< m << " \n";
  
  //initialise parameter matrix, hessian list, and conv vector
  arma::mat par(n+2,m);
  arma::mat fake_par(n+2,m);
  List hes(m);

  arma::vec conv(m);
  double mynan = std::numeric_limits<double>::quiet_NaN();
  //set each element of the fake parameter matrix to nan
  for (int k=0; k<m;k++){
    for (int j=0; j<n+2; j++){
      fake_par.col(k)[j]=mynan;
    }
  }
  
  int failcount = 0;
  // arma::vec fake_parameters(n);
  // fake_parameters.fill(0);
  // arma::vec fake_hessian(m);
  // fake_hessian.fill(0);
  // arma::vec fake_conv(m);
  // fake_conv.fill(0);
  if(sizescale){
    for(int i=0; i < m; i++){
      try{
        List result = NBthDE_paraOptfeat(X, Y.col(i),
                                         threshold0(i)*alpha0, threshold0(i)*alpha,
                                         preci1, 1.0, preci2, x0, calhes);
  
        par.col(i) = (as<arma::vec>(result["par"]));
        hes[i] = result["hes"];
        conv(i) = result["conv"];
      }
      catch (...){
        par.col(i) = fake_par.col(i);
        hes[i] = mynan;
        conv(i) = mynan;
        failcount++;
        //Rcout << i << "failed ";
      }
    }

  } else {
    for(int i=0; i < m; i++){
      try{
        List result = NBthDE_paraOptfeat(X, Y.col(i),
                               alpha0, alpha,
                               preci1, threshold0(i), preci2, x0, calhes);
  
        par.col(i) = (as<arma::vec>(result["par"]));
        hes[i] = result["hes"];
        conv(i) = result["conv"];
      }
        catch (...){
          ;
          par.col(i) = fake_par.col(i);
          hes[i] = mynan;
          conv(i) = mynan;
          failcount++;
          //Rcout << i << "failed ";
      }
    }
  }
  Rcout << "failed rows: " << failcount << "\n";
  return List::create(Named("par") = par,
                      Named("conv") = conv,
                      Named("hes") = hes);
}

