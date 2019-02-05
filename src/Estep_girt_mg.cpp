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
//'@param group a vector.
//'@param ind a design matrix for group.
//'@param resp a design matrix for person.
//'@param MLL a vector
//'@export
// [[Rcpp::export]]


List Estep_girt_mg(DataFrame x,
                NumericVector a0,
                NumericVector b0,
                NumericVector Xq,
                NumericMatrix AX, // multigroup
                NumericVector Yr,
                NumericMatrix BY, // multigroup
                double D,
                IntegerVector group,
                IntegerMatrix ind, // design matrix
                IntegerMatrix resp, // design matrix
                NumericVector MLL
){

  const int nj = x.length();// item n
  const int nn = x.nrows(); // subject n
  const int nq = Xq.length(); // node of theta
  const int nr = Yr.length(); // node of phi
  const int ng = max(group); // group n

  boost::multi_array <double, 4> knqr (boost::extents[ng][nn][nq][nr]);
  boost::multi_array <double, 4> hqr (boost::extents[ng][nn][nq][nr]);
  boost::multi_array <double, 4> Njqr (boost::extents[ng][nj][nq][nr]);
  boost::multi_array <double, 4> rjqr (boost::extents[ng][nj][nq][nr]);
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
  for(int g=0; g<ng; g++){
    for(int r=0; r<nr; r++){
      phi = Yr[r];
      //Rprintf("%d % / 100%\r", 100/nr*(r+1));
      for(int q=0; q<nq; q++){
        theta = Xq[q];
        for(int i=0; i<nn; i++){
          if(group[i] != g+1) continue; // 集団に属さない受験者の部分はスキップ
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
          knqr[g][i][q][r] = t*AX[q]*BY[r];
        }
      }
    }
  }

  // 各受検者のtheta, phiごとに事後分布の重みを計算する。
  // 規格化
  double f = 0; // 周辺対数尤度代入用
  for(int g=0; g<ng; g++){
    for(int i=0; i<nn; i++){
      if(group[i] != g+1) continue; // 集団に属さない受験者の部分はスキップ
      double Z = 0;
      for(int q=0; q<nq; q++){ // 総和を1にするための分母の計算 // sum
        for(int r=0; r<nr; r++){
          Z += knqr[g][i][q][r];
        }
      }
      f += log(Z);
      for(int q=0; q<nq; q++){ // 総和を1にするための分母の計算 // sum
        for(int r=0; r<nr; r++){
          double l = knqr[g][i][q][r];
          hqr[g][i][q][r] =  l / Z;
        }
      }
    }
  }
  MLL.push_back(f);
  if(traits::is_nan<REALSXP>(f)){
    // 対数尤度の計算に失敗したら，計算を中止する。
    stop("Can't calculate marginal log likelihood.");
  }

  for(int g=0; g<ng; g++){
    for(int j=0; j<nj; j++){ // 各分点の期待度数
      if(ind(g,j) == 0) continue;
      for(int q=0; q<nq; q++){
        for(int r=0; r<nr; r++){
          double k = 0;
          for(int i=0; i<nn; i++){ // 欠測値がある場合，項目ごとに受検者数が異なる。
            if(resp(i,j)==0) continue;
            k += hqr[g][i][q][r];
          }
          Njqr[g][j][q][r]= k;
        }
      }
    }
  }

  for(int g=0; g<ng; g++){
    for(int j=0; j<nj; j++){ // 各分点の正答受検者の期待度数
      if(ind(g,j) == 0) continue;
      for(int q=0; q<nq; q++){
        for(int r=0; r<nr; r++){
          double h = 0;
          for(int i=0; i<nn; i++){ // sum
            //double d = resp(i,j);
            if(resp(i,j) == 0) continue;
            h += xall(i,j)*hqr[g][i][q][r];
          }
          rjqr[g][j][q][r] = h;
        }
      }
    }
  }



  return List::create(_["knqr"]=knqr, _["hqr"]=hqr, _["Njqr"]=Njqr, _["rjqr"]=rjqr, _["MLL"]=MLL);


}

