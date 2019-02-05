

#include <Rcpp.h>
#include <boost/multi_array.hpp> // 多次元配列用

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(BH)]]
//'Estimation of population distribution via EM algorithm.
//'
//'This function is new version of \code{MML_EM_dist}
//'@param x Item response matrix.
//'@param para item parameter data.frame estimated by \code{\link{estip}}
//'@param N the number of nodes in integration.
//'@param fc0 a column of first item response.
//'@param ng the number of groups
//'@param gc0 a column of group. the element must be integer and the minimum number must be 1.
//'@param eMLL a convergence criteria(CC) for marginal log likelihood.
//'@param eDIST a CC for population distribution.
//'@param D factor constant.
//'@param fix Don't use. If 1, fix population distribution mean and sigma each EM cycle.
//'@param print How much information you want to display? from 1 to 3. The larger, more information is displayed.
//'@param max maximum value of theta in integration.
//'@param min minimum value of theta in integration.
//'@param maxiter_em maximum iteration time for EM cycle.
//'@param rm_list a vector of item U want to remove for estimation. NOT list.
//'@export
// [[Rcpp::export]]


List estdist(DataFrame x, DataFrame para, const int N = 31, const double eMLL = 1e-6, const double eDIST = 1e-4,
             int fc0 = 2,int ng = 1, int gc0 = 2, const double D = 1.702, const int fix = 1, const int print = 0,
             const double max = 4.0, const double min = -4.0, const int maxiter_em = 200,
             CharacterVector rm_list = CharacterVector::create("NONE")
){
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Update Information
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  struct LocalFunc{ // R function define

    static NumericVector colSums_cpp(NumericMatrix mat){ // R::colSums
      Function f("colSums");
      return f(mat, _["na.rm"] = true);
    }

    static NumericVector rowSums_cpp(NumericMatrix mat){ // R::rowSums
      Function f("rowSums");
      return f(mat, _["na.rm"] = true);
    }

    static NumericVector cor_cpp(NumericVector vec, NumericMatrix mat){ //R::cor
      Function f("cor");
      return f(vec, _["y"] = mat, _["use"] = "pairwise.complete.obs");
    }

    static NumericVector colMeans_cpp(NumericMatrix mat){ // R::colMeans
      Function f("colMeans");
      return f(mat, _["na.rm"] = true);
    }
  };

  const int lc = x.length();// item n
  const int fc = fc0 -1;// Rcpp上では0から要素数がカウントされるので，数値を調整する。

  // 入力データの確認
  const int nj = lc - fc;
  const int ni = x.nrows(); // subject n
  int rm_n = 0;
  if(rm_list[0] != "NONE") rm_n = rm_list.length();  // 削除項目

  const int gc = gc0 -1;

  IntegerVector group (ni,1);
  if(ng != 1) group = x[gc];
  //NumericVector group = x[gc];

  group = group -1; // 集団変数を，C++の要素数に適合するように，-1する。

  NumericMatrix xall (ni, nj);
  CharacterVector Item_temp = x.names();
  CharacterVector Item (nj);
  for (int j=0; j<nj; j++){ // もとのデータフレームから項目反応パタンだけを抜き出し，行列として保存
    int jj = j + fc;
    NumericVector kk = x[jj];
    Item[j] = Item_temp[jj];
    for(int i=0; i<ni; i++){
      double k = kk[i];
      xall(i,j) = k;
    }
  }
  // 削除した項目の列番号（データフレーム）を取得
  IntegerVector rm_id(1,9999);

  if(rm_list[0] != "NONE"){
    Rcout<<"You select Item; "<<rm_list<<"as not use.\n";
    rm_id = IntegerVector(rm_n); // result
    for(int l=0; l<rm_n; l++){
      String r = rm_list[l];
      rm_id[l] = x.findName(r)-(fc-1);
    }
  } else {

  }


  // 欠測値がある箇所を0，そうではない箇所を1とした行列を作成
  NumericMatrix resp (ni,nj);
  for(int i=0; i<ni; i++){
    for(int j=0; j<nj; j++){
      if(NumericVector::is_na(xall(i,j))) resp(i,j) = 0; // NAならば0を
      else resp(i,j) = 1; // NA以外（0か1）ならば1を
    }
  }


  NumericMatrix ind (ng,nj); // 当該集団が，項目に回答したか否か
  for(int g=0; g<ng; g++){
    // 集団が受検した項目を確認
    for(int i=0; i<ni; i++){
      if(group[i] != g) continue;
      for(int j=0; j<nj; j++){
        int q = ind(g,j);
        ind(g,j)= q + resp(i,j);
      }
    }

    for( int j=0; j<nj; j++){
      int q = ind(g,j);
      if(q != 0) ind(g,j) = 1; // 誰か一人でも回答していたら1，集団内が全員無回答なら0
    }
  }


  Rcout << "The number of subject is " << ni << ".\nThe number of item is " << nj-rm_n <<".\nThe number of remove item is "<<rm_n<<".\n";


  // initial values
  NumericMatrix t0(nj,3);
  NumericVector a0 = para["a"];
  NumericVector b0 = para["b"];
  NumericVector c0 = para["c"];
  // NumericVector a0 = para[0];
  // NumericVector b0 = para[1];
  // NumericVector c0 = para[2];
  for(int j=0; j<nj; j++){
    t0(j,0) = a0[j];
    t0(j,1) = b0[j];
    t0(j,2) = c0[j];
  }

  // prior distribution
  double seq = (max - min)/(N-1); // 分点の間隔。これは全母集団で固定。
  NumericVector Xm (N);
  Xm[0] = min;
  for(int m=1; m<N; m++){
    double pre = Xm[m-1];
    Xm[m] = pre + seq;
  }


  // 母集団分布の初期値を決定。
  NumericMatrix Um (N,ng); // 推定母集団分布計算用
  NumericMatrix dist (0);


  // 結果代入用の行列，ベクトルを準備　＝　こうしないと，最終的にスコープできない。
  boost::multi_array <double, 3> Gim (boost::extents[ng][ni][N]);
  boost::multi_array <double, 3> Lim (boost::extents[ng][ni][N]);
  boost::multi_array <double, 3> Nm (boost::extents[ng][nj][N]);
  boost::multi_array <double, 3> rjm (boost::extents[ng][nj][N]);
  NumericMatrix skip_para (nj,3); // 推定に失敗した項目番号の代入用

  /////////////////////////
  // EM step starts here //
  /////////////////////////
  // 推定母集団分布の計算
  // 初期値に一様分布をもちいた場合，Easy Estimationの計算結果とおおむね一致するはず。

  Rprintf("Start calculating estimated population distribution.\n");

  // result
  NumericVector mean_pop(ng);
  NumericVector sd_pop(ng);
  double diff2 = 0;
  NumericVector MLL2 (0);

  // 一様分布を初期値に設定する。
  NumericMatrix unif_dist (0);
  double NN = Um.nrow();
  for(int i=0; i<N; i++){
    for(int j=0; j<ng; j++){
      Um(i,j) = 1/NN;
    }
  }


  // start loop
  int count3 =0;
  int conv3 = 0;
  LogicalVector conv0(ng);

  while(conv3 == 0){
    count3 += 1;

    // calculate Lim
    double a,b,c,t,xx,tt,u;
    int g,i,m,j;
    for(g=0; g<ng; g++){
      for(i=0; i<ni; i++){
        if(group[i] != g) continue; // 集団に属さない受験者の部分はスキップ
        for(m=0; m<N; m++){
          t = 1;
          xx = Xm[m];
          tt = 0;
          for(j=0; j<nj; j++){
            //if(resp(i,j)==0) continue; // NAの反応パタンのところは計算ループから外れる。
            a = t0(j,0);
            if(a==0)continue;
            b = t0(j,1);
            c = t0(j,2);
            u = xall(i,j);
            if(u == 1){ // 正答した場合の尤度
              tt = c+(1.0-c)/(1.0+exp(-D*a*(xx-b)));
            } else if (u == 0) { // 誤答した場合の尤度
              tt = 1.0-(c+(1.0-c)/(1.0+exp(-D*a*(xx-b))));
            } else { // 欠測値の場合は尤度を計算しないので，1
              tt = 1;
            }
            t = t * tt; // sum
          }
          //boost::multi_array<double, 3>::index idx = {{g,i,m}};
          Lim[g][i][m] = t;
        }
      }
    }


    //Rcout<<"marginal log likelihood calculation.\n";
    // Calculation marginal log likelihood

    double f = 0; // 周辺対数尤度代入用
    double ff,fff,l,w;

    for(int g=0; g<ng; g++){
      for(int i=0; i<ni; i++){
        if(group[i] != g) continue; // 集団に属さない受験者の部分はスキップ
        ff = 0; // 受検者iの周辺対数尤度代入用
        for(int m=0; m<N; m++){
          l = Lim[g][i][m];
          w = Um(m,g);
          fff =  l*w ;
          ff = ff + fff;
        }
        f = f + log(ff);
      }
    }

    MLL2.push_back(f);
    diff2 = fabs(MLL2[count3-1] - MLL2[count3]); // 周辺対数尤度の差

    //Rcout << "-2 Marginal Log Likelihood is " << -2 * f <<"\r";
    Rprintf("%d times -2 Marginal Log Likelihood is %f \r",count3,-2*f); // 小数点形式の出力に対応


    // calculate Gim
    double uu;
    for(int g=0; g<ng; g++){
      for(int i=0; i<ni; i++){
        if(group[i] != g) continue; // 集団に属さない受験者の部分はスキップ
        u = 0;
        for(int m=0; m<N; m++){ // 総和を1にするための分母の計算 // sum
          l = Lim[g][i][m];
          w = Um(m,g);
          uu = l*w;
          u = u + uu;
        }
        for(int m=0; m<N; m++){
          l = Lim[g][i][m];
          w = Um(m,g);
          Gim[g][i][m] =  l * w / u;
        }
      }
    }

    // calculate Nm

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
            k = k + kk; //* d;
          }
          Nm[g][j][m] = k;
        }
      }
    }

    // 一様分布の確率密度をリセット
    for(int m=0; m<N; m++){
      for(int g=0; g<ng; g++){
        Um(m,g) = 0;
      }
    }

    // 母集団ごとの平均と標準偏差を計算
    for(int g=0; g<ng; g++){
      // probability
      double mean_g = 0;
      double sd_g = 0;
      //int counti = 0; // 母集団の受検した項目数カウント用
      for(int j=0; j<nj; j++){
        // 集団内のだれも回答をしていなかったら，次の項目のループへスキップ。
        if(ind(g,j) == 0) continue;

        //counti = counti + 1;
        double mean_j = 0;
        double sNjm = 0; // 項目ごとの分点の受験者の期待度数の総和代入用
        for(int m=0; m<N; m++){  // sum
          double Njm = Nm[g][j][m];
          sNjm = sNjm + Njm;
        }

        for(int m=0; m<N; m++){ // sum
          double Njm = Nm[g][j][m];
          double xm = Xm[m];
          double tempM = xm * Njm/sNjm; // E(X) = \Sigma(x_i * prob_i)
          Um(m,g) += Njm/sNjm;
          mean_j += tempM; //
        } // end of m
        mean_g = mean_g + mean_j; // 項目ごとの母集団平均を足し上げたもの

        // sd

        double var_j = 0;
        double sd_j = 0;
        for(int m=0; m<N; m++){
          double Njm = Nm[g][j][m];
          double xm = Xm[m];
          double tempV = (xm -mean_j)*(xm -mean_j)*Njm/sNjm; // V(X) = Sigma(x_i-E(X))^2*prob_i
          var_j = var_j + tempV;
        }// end of m
        sd_j = std::sqrt(var_j); // 項目ごとに計算した標準偏差
        sd_g = sd_g + sd_j; // 項目ごとの母集団標準偏差を足し上げたもの

      } // end of j

      int counti = sum(ind(g,_));
      // 収束判定用
      if(fabs(mean_pop[g] - mean_g/counti)<eDIST && fabs(sd_pop[g] - sd_g/counti)<eDIST ) conv0[g] = 1;

      mean_pop[g] = mean_g/counti; // 全項目における母集団平均の，平均
      sd_pop[g] = sd_g/counti; // 全項目における母集団標準偏差の，平均
      if(print >= 1) Rcout << "mean: "<< mean_pop[g]<<", sd: "<<sd_pop[g]<<"\n";

      for(int m=0; m<N; m++){
        Um(m,g) /= counti;
      }

    } // end of g

    if(diff2 < eMLL || is_true(all(conv0)) ) conv3 = 1;
    //break;
  }

  unif_dist = cbind(Xm,Um);

  List res = List::create(_["population_dist"] = unif_dist, _["population_mean"] = mean_pop, _["population_sd"] = sd_pop,
                          _["MLL"] = MLL2, _["conv"] = conv3, _["count"] = count3,
                          _["skip_para"] = skip_para, _["rm_n"] = rm_n, _["rm_id"] = rm_id
  );

  return res;

}


