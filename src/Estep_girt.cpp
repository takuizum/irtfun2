#include <Rcpp.h>
#include <boost/multi_array.hpp>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(BH)]]

//'Estep in girt.
//'
//'@param x DataFrame.
//'@param a0 discrinimation parameter vector.
//'@param b0 difficulty parameter vector.
//'@param Xq node of theta dist.
//'@param AX weight of theta dist.
//'@param Yr node of phi dist.
//'@param BY weighit of phi dist.
//'@param D factor constant.
//'@export
// [[Rcpp::export]]


List Estep_girt(DataFrame x,
                NumericVector a0,
                NumericVector b0,
                NumericVector Xq,
                NumericVector AX,
                NumericVector Yr,
                NumericVector BY,
                double D
){

  const int nj = x.length();// item n
  const int nn = x.nrows(); // subject n
  const int nq = Xq.length();
  const int nr = Yr.length();

  boost::multi_array <double, 3> knqr (boost::extents[nn][nq][nr]);
  //Dimension d(nn,nq,nr);
  //NumericVector knqr(d);

  // もとのデータフレームから項目反応パタンだけを抜き出し，行列として保存
  NumericMatrix xall(nn, nj);
  for (int j=0; j<nj; j++){
    NumericVector kk = x[j];
    for(int i=0; i<nn; i++){
      double k = kk[i];
      xall(i,j) = k;
    }
  }

  int u;
  double phi, theta, a, b, tt, t, A, e, p;
  for(int r=0; r<nr; r++){
    phi = Yr[r];
    //Rprintf("%d % / 100%\r", 100/nr*(r+1));
    for(int q=0; q<nq; q++){
      theta = Xq[q];
      for(int i=0; i<nn; i++){
        t = 1.0;
        for(int j=0; j<nj; j++){
          a = a0[j];
          b = b0[j];
          A = sqrt(1+phi*phi*a*a);
          e = exp(-D*a/A*(theta-b));
          p = 1/(1+e);
          u = xall(i,j);
          if(u == 1){ // correct
            tt = p;
          } else if (u == 0) { // incorrect
            tt = 1.0-p;
          } else { // NA
            tt = 1;
          }
          t = t*tt;
        }
        knqr[i][q][r] = t*AX[q]*BY[r];
      }
    }
  }


  return List::create(_["knqr"]=knqr);


}

