#include <Rcpp.h>
#include <boost/multi_array.hpp>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(BH)]]

//'Estep in girt.
//'
//'@param xall Item response matrix.
//'@param t0 item parameter vector
//'@param Xm node of theta dist.
//'@param Wm weight of theta dist.
//'@param group a vector.
//'@param ind a design matrix for group.
//'@param resp a design matrix for person.
//'@param D factor constant.
//'@param MLL a vector
//'@export
// [[Rcpp::export]]

List Estep_irt(
  IntegerMatrix xall,
  NumericMatrix t0, // a, b, c
  NumericVector Xm,
  NumericMatrix Wm,
  IntegerVector group,
  IntegerMatrix ind, // design matrix
  IntegerMatrix resp, // design matrix
  double D,
  NumericVector MLL
){

  int nj = xall.ncol();// item n
  int ni = xall.nrow(); // subject n
  int N  = Xm.length(); // node n
  int ng = max(group); // group n

  double a,b,c,t,tt,u,x;
  int g,i,m,j;

  boost::multi_array <double, 3> Lim (boost::extents[ng][ni][N]);
  boost::multi_array <double, 3> Gim (boost::extents[ng][ni][N]);
  boost::multi_array <double, 3> Njm (boost::extents[ng][nj][N]);
  boost::multi_array <double, 3> rjm (boost::extents[ng][nj][N]);

  for(g=0; g<ng; g++){
    for(i=0; i<ni; i++){
      if(group[i] != g+1) continue; // 集団に属さない受験者の部分はスキップ
      for(m=0; m<N; m++){
        t = 1;
        x = Xm[m];
        tt = 0;
        for(j=0; j<nj; j++){
          //if(resp(i,j)==0) continue; // NAの反応パタンのところは計算ループから外れる。
          a = t0(j,0);
          if(a == 0) continue;
          b = t0(j,1);
          c = t0(j,2);
          u = xall(i,j);
          if(u == 1){ // 正答した場合の尤度
            tt = c+(1.0-c)/(1.0+exp(-D*a*(x-b)));
          } else if (u == 0) { // 誤答した場合の尤度
            tt = 1.0-(c+(1.0-c)/(1.0+exp(-D*a*(x-b))));
          } else { // 欠測値の場合は尤度を計算しないので，1
            tt = 1;
          }
          t = t * tt; // sum
        }
        Lim[g][i][m] = t;
      }
    }
  }

  // 各受検者のthetaごとに事後分布の重みを計算する。
  double uu;
  double f = 0; // 周辺対数尤度代入用
  double l,w;
  for(int g=0; g<ng; g++){
    for(int i=0; i<ni; i++){
      if(group[i] != g+1) continue; // 集団に属さない受験者の部分はスキップ
      u = 0;
      for(int m=0; m<N; m++){ // 総和を1にするための分母の計算 // sum
        l = Lim[g][i][m];
        w = Wm(m,g);
        uu = l*w;
        u += uu;
      }
      f += log(u);
      for(int m=0; m<N; m++){
        l = Lim[g][i][m];
        w = Wm(m,g);
        Gim[g][i][m] =  l * w / u;
      }

    }
  }
  MLL.push_back(f);
  if(traits::is_nan<REALSXP>(f)){
    // 対数尤度の計算に失敗したら，計算を中止する。
    stop("Can't calculate marginal log likelihood.");
  }

  //Rcout<<"expected frequency of subjects in each nodes calculation.\n";
  double k,kk;
  for(int g=0; g<ng; g++){
    for(int j=0; j<nj; j++){ // 各分点の期待度数
      if(ind(g,j) == 0) continue;
      for(int m=0; m<N; m++){
        k = 0;
        for(int i=0; i<ni; i++){ // 欠測値がある場合，項目ごとに受検者数が異なる。
          if(resp(i,j)==0) continue;
          //double d = resp(i,j);
          kk = Gim[g][i][m];
          k += kk; //* d;
        }
        Njm[g][j][m]= k;
      }
    }
  }

  double h,gg,hh;
  for(int g=0; g<ng; g++){
    for(int j=0; j<nj; j++){ // 各分点の正答受検者の期待度数
      if(ind(g,j) == 0) continue;
      for(int m=0; m<N; m++){
        h = 0;
        for(int i=0; i<ni; i++){ // sum
          //double d = resp(i,j);
          if(resp(i,j) == 0) continue;
          u = xall(i,j); // そもそもその項目に回答していない場合は，度数に数え上げない
          gg = Gim[g][i][m];
          hh = u*gg;
          h += hh;
        }
        rjm[g][j][m] = h;
      }
    }
  }

  return List::create(_["Njm"]=Njm, _["rjm"]=rjm, _["Lim"]=Lim, _["Gim"]=Gim, _["MLL"]=MLL);

}
