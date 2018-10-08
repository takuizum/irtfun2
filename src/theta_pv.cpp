#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//'Sampling plausible values(PVs) based on rejection sampling.
//'
//'@param x response data
//'@param nofrands the number of PVs
//'@param eap_apply a numeric vector of estimated EAP.
//'@param const_apply a numeric vector of const of a posterior distribution.
//'@param map_apply a numeric vector of estimated MAP.
//'@param n the number of subjects.
//'@param maxtheta a max value of theta.
//'@param mintheta a mininum value of theta.
//'@param a an item slope parameter vector.
//'@param b an item location parameter vector.
//'@param c an item asymptote parameter vector.
//'@param D a factor constant.
//'@param mu p hyper parameter of a posterior distribution.
//'@param sigma same as above.
//'@export
// [[Rcpp::export]]

NumericMatrix theta_pv(DataFrame x, const int nofrands, NumericVector eap_apply, NumericVector const_apply, NumericVector map_apply,
                       const int n, double maxtheta, double mintheta,
                       NumericVector a, NumericVector b, NumericVector c, const double D, const double mu, const double sigma){
  // output matrix
  NumericMatrix pv (n,nofrands);

  double m = a.length();
  // convert DataFrame to Numeric Matrix
  NumericMatrix xall (n,m);
  for (int j=0; j<m; j++){ // もとのデータフレームから項目反応パタンだけを抜き出し，行列として保存
    NumericVector temp = x[j];
    for(int i=0; i<n; i++){
      double temp2 = temp[i];
      xall(i,j) = temp2;
    }
  }

  // start rejection sampling
  int NO = 0;
  int YES = 0;
  int times = 0;
  for(int k=0; k<n; k++){
    NumericVector xi = xall(k,_);
    times += 1;
    double eap = eap_apply[k];
    double const_cpp = const_apply[k];
    double map = map_apply[k];
    //Rprintf("exrtact single subject eap is %f,map is %f,const is %f \n",eap,map,const_cpp);

    double zmin = mintheta + eap;
    double zmax = maxtheta + eap;
    double prior_m = R::dnorm(map,mu,sigma,false);
    double z,y,a_j,b_j,c_j,xij,yheight;
    double mpdc = 0;
    // calculate log likelihood in MAP estimator.
    for(int j=0; j<m; j++){
      a_j = a[j];
      if(a_j == 0) continue;
      b_j = b[j];
      c_j = c[j];
      xij = xi[j];
      double mpdc_t;
      if(xij == 1){
        mpdc_t = c_j+(1.0-c_j)/(1.0+exp(-D*a_j*(map-b_j)));
      } else if(xij==0){
        mpdc_t = 1 - (c_j+(1.0-c_j)/(1.0+exp(-D*a_j*(map-b_j))));
      } else {
        mpdc_t = 0;
      }
      mpdc += log(mpdc_t);
      //Rprintf("mpdc is %f\n",mpdc);

    }
    yheight = exp(mpdc)*prior_m/const_cpp * 1.001;

    //Rprintf("calculate y-height is %f, mpdc is %f, const is %f.\n",yheight,exp(mpdc),const_cpp);

    // calculate yheight
    int nofpv = 0;
    //Rprintf("\nPV sampling\n");
    while(nofpv < nofrands ){
      z = R::runif(zmin,zmax);
      y = R::runif(0,yheight);
      // calculate fg
      double prior_z = R::dnorm(z,mu,sigma,false);
      double fg = 0;
      for(int j=0; j<m; j++){
        a_j = a[j];
        if(a_j == 0) continue;
        b_j = b[j];
        c_j = c[j];
        double fg_t;
        if(xij == 1){
          fg_t = c_j+(1.0-c_j)/(1.0+exp(-D*a_j*(z-b_j))) ;
        } else if(xij==0){
          fg_t = 1 - (c_j+(1.0-c_j)/(1.0+exp(-D*a_j*(z-b_j))));
        } else {
          fg_t = 0;
        }
        fg += log(fg_t) ;
      }
      fg = exp(fg)*prior_z/const_cpp;
      //Rprintf("y-height is %f, y is %f, z is %f, fgvalue is %f \r",yheight,y,z,fg);
      if(y <= fg){
        pv(k,nofpv) = z;
        nofpv += 1;
        YES += 1;
      } else {
        NO += 1;
      }
    }
    Rprintf("NOW---%d / %d , ACCEPT---%d, REJECT---%d \r",k+1,n,YES,NO);

    //if(times == counter){
    //  Rprintf("%d / n \r",k);
    //  times = 0;
    //}
  }
  return pv;
}

