#include <Rcpp.h>
#include <boost/multi_array.hpp>

using namespace Rcpp;


// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(BH)]]

//'Estimate item parameter for binary{0,1} response data.
//'
//'1PL,2PL,3PL,Bayes1PL,Bayes2PL and multigroup estimation is avairable now. U must install C++ compiler(Rtools for windows or Xcode for Mac)in your PC or Mac.
//'@param x an item response data which class is data.frame object.
//'@param model0 Character.U can select which one, "1PL","2PL","3PL".
//'@param N the number of nodes in integration.
//'@param bg0 the number of base grade.
//'@param eMLL a convergence criteria(CC) for marginal log likelihood.
//'@param eEM a CC in EM cycle.
//'@param eM a CC in M step.
//'@param emu a CC for population distribution mean.
//'@param esd a CC for population distribution standard deviation.
//'@param fc0 a column of first item response.
//'@param ng the number of groups
//'@param gc0 a column of group. the element must be integer and the minimum number must be 1.
//'@param D a scaling constant.
//'@param fix Don't use. If 1, fix population distribution mean and sigma each EM cycle.
//'@param print How much information you want to display? from 1 to 3. The larger, more information is displayed.
//'@param ic initial value for guessing parameter.
//'@param max maximum value of theta in integration.
//'@param min minimum value of theta in integration.
//'@param mu hyperparameter for theta dist.
//'@param sigma same as above.
//'@param Bayes If 1, marginal Bayesian estimation runs.
//'@param method Optimising method in M step. "Fisher_Scoring" or "Newton_Raphson" can be selected.
//'@param mu_a hyperparameter of log normal dist for slope parameter.
//'@param sigma_a same as above.
//'@param mu_b hyperparameter of normal dist for location parameter.
//'@param sigma_b same as above.
//'@param mu_c hyperparameter for lower asymptote parameter.
//'@param w_c weight of a Beta dist.
//'@param min_a minimum value of slope parameter.
//'@param maxabs_b maximum absolute value of location parameter.
//'@param maxiter_em maximum iteration time for EM cycle.
//'@param maxiter_j maximum iteration time for Newton Raphton in M step.
//'@param maxskip_j Dont use.
//'@param rm_list a vector of item U want to remove for estimation. NOT list.
//'@param thdist Which distribution do you want `normal` or `empirical` for E step.
//'@param EM_dist If 1, calculate esimated population distribution via EM argorithm.
//'@examples
//'res <- estip(x=sim_data_1, fc=2)
//'# check the parameters
//'res$para
//'res$SE
//'@export
// [[Rcpp::export]]


List estip (DataFrame x,
            CharacterVector model0 = CharacterVector::create("2PL") ,
            const int N = 31,
            const int bg0 = 1,
            int fc0 = 2,
            int ng = 1,
            int gc0 = 2,
            const double eMLL = 1e-6,
            const double eEM = 1e-4,
            const double eM = 1e-3,
            const double emu = 1e-3,
            const double esd = 1e-3,
            const double D = 1.702,
            const double ic = 1/5,
            const double max = 6.0,
            const double min = -6.0,
            const double mu = 0,
            const double sigma = 1,
            const int Bayes = 0,
            const String method = "Fisher_Scoring",
            const double mu_a = 0,
            const double sigma_a = 1,
            const double mu_b = 0,
            const double sigma_b = 2,
            const double mu_c = 1/5,
            const double w_c = 5,
            const int fix = 1,
            const int print = 0,
            const double min_a = 0.1,
            const double maxabs_b = 20,
            const int maxiter_em = 200,
            const int maxiter_j = 20,
            const int maxskip_j = 0,
            CharacterVector rm_list = CharacterVector::create("NONE"),
            const String thdist = "normal",
            const int EM_dist = 1
){
  // argument check
  for(int j=0; j<model0.length(); j++){
    Rprintf("Checking the model string vector %i !\r", j+1);
    if(model0[j] != "1PL" && model0[j] != "2PL" && model0[j] != "3PL") stop("Errpr! option string of 'model' is incorrect. You should select '1~3PL'.");
  }
  Rprintf("\n");
  if(thdist != "normal" && thdist != "empirical") stop("Errpr! option string of 'thdist' is incorrect. You should select 'normal' or 'empirical'.");

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
  const int gc = gc0 -1; // the column of group n
  const int bg = bg0 -1; // base grade

  IntegerVector group (ni,1);
  if(ng != 1) group = x[gc];
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
  int rm_n = 0;

  if(rm_list[0] != "NONE"){
    rm_n = rm_list.length();  // 削除項目
    Rcout<<"You select Item; "<<rm_list<<" as not use.\n";
    rm_id = IntegerVector(rm_n); // result
    for(int l=0; l<rm_n; l++){
      String r = rm_list[l];
      rm_id[l] = x.findName(r)-(fc0-1);
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
  CharacterVector model (nj, model0[0]);
  const int nmodel = model0.length();
  if(nmodel == nj){
    Rprintf("You selected model variable yourself.\n");
    for(int j=0; j<nj; j++){
      model[j] = model0[j];
    }
  }
  IntegerVector cat (nj);
  NumericMatrix t0 (nj, 3); // 初期値代入用の行列　いまのところ2PLのみ
  NumericMatrix t0m (nj,3);
  NumericMatrix initial (nj,3);
  //NumericMatrix t1m (nj,cat);

  NumericVector rs = LocalFunc::rowSums_cpp(xall);
  NumericVector r  = LocalFunc::cor_cpp(rs, xall);
  NumericVector p  = LocalFunc::colMeans_cpp(xall);

  double a0, b0, c0;
  for(int j=0; j<nj; j++){

    if(model[j] == "1PL") cat[j] = 1;
    if(model[j] == "2PL") cat[j] = 2;
    if(model[j] == "3PL") cat[j] = 3;

    a0 = D*r[j]/sqrt(1-r[j]*r[j]); // P.BIS
    b0 = -log(p[j]/(1.0-p[j]));
    c0 = ic;

    if(a0 < 0.1) a0 = 0.3; // 初期値が小さすぎる場合に，大きい値に補正する。

    // a parameter
    if(model[j] == "2PL" || model[j] == "3PL"){
      t0(j,0) =  a0;
      t0m(j,0) = a0;
      initial(j,0) = a0;
    } else {
      t0(j,0) = 1.702/D;
      t0m(j,0) = 1.702/D;
      initial(j,0) = 1.702/D;
    }
    // b parameter
    t0(j,1) =  b0;
    t0m(j,1) = b0;
    initial(j,1) = b0;

    if(model[j] == "3PL"){
      // c parameter
      t0(j,2) =  c0;
      t0m(j,2) = c0;
      initial(j,2) = c0;
    }

  }

  double alpha_c = w_c*mu_c + 1; // 識別力の事前分布計算に用いる定数
  double beta_c = w_c*(1-mu_c) + 1;
  // default alpha=2, beta=5

  // 削除項目はすべて計算をスキップするので
  if(rm_list[0] != "NONE"){
    for(int l : rm_id){
      for(int g=0; g<ng; g++){
        ind(g,l) = 0;
      }
      for(int i=0; i<ni; i++){
        resp(i,l) = 0;
      }
      for(int p=0; p<3; p++){
        t0(l,p) = 0;
        t0m(l,p) = 0;
        initial(l,p) = 0;
      }
    }
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
  NumericMatrix Wm(N,ng); // 母集団ごとに事前分布の重みを計算（初期値なので，どの集団もおなじ
  NumericMatrix Um (N,ng); // 推定母集団分布計算用
  NumericVector mean (ng);
  NumericVector sd (ng);
  NumericMatrix dist (N,ng+1);
  dist(_,0) = Xm;
  mean = mean + mu;
  sd = sd + sigma;


  // 結果代入用の行列，ベクトルを準備　＝　こうしないと，最終的にスコープできない。
  boost::multi_array <double, 3> Gim (boost::extents[ng][ni][N]);
  boost::multi_array <double, 3> Lim (boost::extents[ng][ni][N]);
  boost::multi_array <double, 3> Nm (boost::extents[ng][nj][N]);
  boost::multi_array <double, 3> rjm (boost::extents[ng][nj][N]);
  LogicalVector conv (sum(cat)); // 収束判定代入用ベクトル。
  // 周辺対数尤度代入用ベクトル
  NumericVector MLL(1);
  double diff = 0;
  NumericMatrix SE (nj,3); // 推定の標準誤差代入用
  NumericMatrix skip_para (nj,3); // 推定に失敗した項目番号の代入用

  /////////////////////////
  // EM step starts here //
  /////////////////////////
  Rprintf("Start calculating estimated item parameters.\n");
  //convergence criteria
  int conv1 = 0; // for EM step

  //n of iteraiton
  int count1 = 0;
  int count2 = 0;

  while(conv1 != 1){ // EM step

    ///////////////////////
    //for debug.///////////
    //break; //////////////
    ///////////////////////

    // check cancel command Ctrl + c
    Rcpp::checkUserInterrupt();

    count1 = count1 + 1; // count up

    checkUserInterrupt(); // If user input Ctrl + c, stop EM cycle.

    if(print>=1)Rcout <<count1 << " times E-step calculation. \r";


    // 事前分布の重みを，母集団ごとに計算（前回のMステップの計算結果を考慮する。
    if(thdist == "normal" || (thdist == "empirical" && count1 == 1)){
      for(int g=0; g<ng; g++){
        double dn = 0;
        NumericVector dnn (N);
        for(int m=0; m<N; m++){
          dnn[m] = R::dnorm(Xm[m], mean[g], sd[g], false);
          dn = dn + dnn[m];
        }
        for(int m=0; m<N; m++){
          double wn = dnn[m] / dn;
          Wm(m,g) = wn;
        }
        dist(_,g+1) = Wm(_,g);
      }
    } else if(thdist == "empirical" && count1 != 1){
      for(int g=0; g<ng; g++){
        for(int m=0; m<N; m++){
          Wm(m,g) = Um(m,g);
        }
        dist(_,g+1) = Wm(_,g);
      }
    }


    //Rcout<<"likelihood of all subject in each nodes calculation.\n";
    // 受検者iの項目反応パタンに基づく尤度を，求積点ごとに求める。

    double a,b,c,t,tt,x,u;
    int g,i,m,j;

    for(g=0; g<ng; g++){
      for(i=0; i<ni; i++){
        if(group[i] != g) continue; // 集団に属さない受験者の部分はスキップ
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
            // std::powやlog，expの関数を持ってくると急激に遅くなる。それよりかif文の方がいくらか高速。
            //tt =
            //  std::pow(c+(1.0-c)/(1.0+exp(-D*a*(x-b))),u)* // log of correct probabbility
            //  std::pow(1.0-(c+(1.0-c)/(1.0+exp(-D*a*(x-b)))),(1.0-u)); // log of incorrect probability
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
          w = Wm(m,g);
          fff =  l*w ;
          ff += fff;
        }
        f += log(ff);
      }
    }

    MLL.push_back(f);
    diff = fabs(MLL[count1-1] - MLL[count1]); // 周辺対数尤度の差

    //Rcout << "-2 Marginal Log Likelihood is " << -2 * f <<"\r";
    Rprintf("%d times -2 Marginal Log Likelihood is %f \r",count1,-2*f); // 小数点形式の出力に対応

    if(traits::is_nan<REALSXP>(f)){
      // 対数尤度の計算に失敗したら，計算を中止する。
      stop("Can't calculate marginal log likelihood.");
    }

    //Rcout<<"weight of a posterir distribution calculation.\n";
    // 各受検者のthetaごとに事後分布の重みを計算する。
    double uu;
    for(int g=0; g<ng; g++){
      for(int i=0; i<ni; i++){
        if(group[i] != g) continue; // 集団に属さない受験者の部分はスキップ
        u = 0;
        for(int m=0; m<N; m++){ // 総和を1にするための分母の計算 // sum
          l = Lim[g][i][m];
          w = Wm(m,g);
          uu = l*w;
          u += uu;
        }
        for(int m=0; m<N; m++){
          l = Lim[g][i][m];
          w = Wm(m,g);
          Gim[g][i][m] =  l * w / u;
        }
      }
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
          Nm[g][j][m]= k;
        }
      }
    }

    //Rcout<<"expected freqency of correct response number in each nodes calculation.\n";

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

    //Rcout<<rjm(10,20)<<"\n";

    //Rcout<<"expected log complete likelihood calculation.\n";

    //const int e_ell = 0,//'@param e_ell If 1 use a CC of expected log likelihood.


    //double ell=0;
    //double rr,nm,elll;
    //int ell_check = 0;
    //if(e_ell==1){

    //  if(Bayes == 0){
    //    // 期待対数完全データ尤度関数
    //    for(int g=0; g<ng; g++){
    //      for(int j=0; j<nj; j++){
    //        if(ind(g,j) == 0) continue;
    //        a = t0(j,0);
    //        b = t0(j,1);
    //        c = t0(j,2);
    //        for(int m=0; m<N; m++){
    //          rr = rjm[g][j][m];
    //          nm = Nm[g][j][m];
    //          x = Xm[m];
    //          elll = rr*log(c+(1.0-c)/(1.0+exp(-D*a*(x-b)))) +
    //            (nm - rr) * log(1-(c+(1.0-c)/(1.0+exp(-D*a*(x-b)))));
    //          ell += + elll;
    //          if(traits::is_nan<REALSXP>(ell)) ell_check = j+1;
    //        }
    //      }
    //    }
    //  } else {    // Bayes parameter estimation
    //    double da,db,dc;
    //    for(int g=0; g<ng; g++){
    //      for(int j=0; j<nj; j++){
    //        if(ind(g,j) == 0) continue;
    //        a = t0(j,0);
    //        b = t0(j,1);
    //        c = t0(j,2);
    //        for(int m=0; m<N; m++){
    //          rr = rjm[g][j][m];
    //          nm = Nm[g][j][m];
    //          x = Xm[m];
    //          elll = rr*log(c+(1.0-c)/(1.0+exp(-D*a*(x-b)))) +
    //            (nm - rr) * log(1-(c+(1.0-c)/(1.0+exp(-D*a*(x-b)))));
    //          ell += elll;
    //        }
    //        da = -1/a - (log(a)-mu_a)/a*sigma_a*sigma_a;
    //        db = -(b-mu_b)/sigma_b*sigma_b;
    //        dc = (alpha_c-2)/c - (beta_c-2)/(1-c); // たぶんここのパラメタの計算が非数になる原因だと思われる。
    //        ell += (da + db + dc);
    //        if(traits::is_nan<REALSXP>(ell)) ell_check = j+1;
    //      }
    //    }
    //  }


    //  if(print >= 1) Rcout << "\n expected log complete data likelihood is " <<ell;

    //  if(traits::is_nan<REALSXP>(ell) && Bayes != 1){
    //    // 対数尤度の計算に失敗したら，計算を中止する。
    //    stop("Can't calculate expected log complete data likelihood.\nItem %d is not converged. Check P.BIS!", ell_check);
    //  }
    //}
    ////////////////////////
    // M step starts here //
    ////////////////////////
    if(print>=1) Rcout << count1 << " times M-step calculation.\r";

    // 項目ひとつごとに最大化
    for(int j=0; j<nj; j++){
      //Rcout << j << " item calculation.\r";
      LogicalVector rm_j = j==rm_id;
      if(is_true(any(rm_j))){
        //Rcout<<"CHECK!!!"<< rm_j <<"\n";
        continue;
      }

      if(skip_para(j,0) > maxskip_j || skip_para(j,1) > maxskip_j ||  skip_para(j,2) > maxskip_j){
        //Rcout<<"skip item "<<j+1<<"\r";
        continue;
      }

      // convergence criteria
      int conv2 = 0; // for M step

      NumericVector indj = ind(_,j); // item index

      while(conv2 != 1){ // Fisher scoring starts in item j

        count2 = count2 + 1;
        if(count2 >= maxiter_j)conv2 = 1;

        double a = t0m(j,0);
        double b = t0m(j,1);
        double c = t0m(j,2);
        // each element of Fisher scoring matrix and gradient vector

        double daa=0;
        double dbb=0;
        double dcc=0;
        double dab=0;
        double dac=0;
        double dbc=0;
        double da=0;
        double db=0;
        double dc=0;
        double d=0;
        double rr=0;
        double rrr=0;
        double nm=0;
        double nmm=0;
        double P=0;
        double Q=0;

        if(method == "Fisher_Scoring"){

          for(int m=0; m<N; m++){ // 行列計算
            nm = 0;
            for(int g=0; g<ng; g++){
              nmm = Nm[g][j][m];
              d = indj[g]; // 当該集団が受検していれば1，そうでなければ0
              nm = nm + nmm * d;
            }
            rr = 0;
            for(int g=0; g<ng; g++){
              rrr = rjm[g][j][m];
              d = indj[g]; // 当該集団が受検していれば1，そうでなければ0
              rr = rr + rrr * d;
            }
            double xm = Xm[m];
            P = c + (1-c)/(1+exp(-D*a*(xm-b)));
            Q = 1-P;
            // the element of Fisher Information Matrix
            // matrix[0,0]
            double daam = nm * (P-c) * (P-c) * Q * (xm - b) * (xm - b) * D * D
              / ((1-c) * (1-c) * P);
            daa = daa + daam;
            // matrix[1,1]
            double dbbm = nm * (P-c) * (P-c) * Q * a * a * D * D
              / ((1-c) * (1-c) * P);
            dbb = dbb + dbbm;
            // matrix[2,2]
            double dccm = nm * Q
              / (P * (1-c) * (1-c));
            dcc = dcc + dccm;
            //matrix[1,0], [0,1]
            double dabm = -nm * (P-c) * (P-c) * Q * a * (xm - b) * D * D
              / ((1-c) * (1-c) * P);
            dab = dab + dabm;
            // matrix[2,0],[0,2]
            double dacm = nm * (P-c) * Q * (xm - b) * D
              / ((1-c) * (1-c) * P);
            dac = dac + dacm;
            // matrix[2,1],[1,2]
            double dbcm = -nm * (P-c) * Q * a * D
              / ((1-c) * (1-c) * P);
            dbc = dbc + dbcm;

            // the elements of gradient vector
            double dad = D*(xm - b) * (rr - nm * P) * (P-c)
              / ((1-c) * P);
            da = da + dad;
            double dbd = -D*a*(rr - nm * P) * (P-c)
              / ((1-c) * P);
            db = db + dbd;
            double dcd = (rr - nm * P)
              /( (1-c) * P);
            dc = dc + dcd;
          } // 行列計算の下準備終了 // end of m

          if(Bayes==1){
            // the elements of gradient vector
            da += -1/a - (log(a)-mu_a)/(a*sigma_a*sigma_a);
            db += -(b-mu_b)/(sigma_b*sigma_b);
            dc += (alpha_c-2)/c - (beta_c-2)/(1-c);
            // Information matrix
            //daa += 1/(a*a) - (1-log(a)+mu_a)/(a*a*sigma_a*sigma_a);
            //dbb += -1/(sigma_b*sigma_b);
            //dcc += (alpha_c-2)/(c*c) - (beta_c-2)/(1-c)*(1-c);
            //非対角要素は通常の二階偏微分と同じ。
          }
        }else if(method == "Newton_Raphson"){
          for(int m=0; m<N; m++){ // 行列計算
            nm = 0;
            for(int g=0; g<ng; g++){
              nmm = Nm[g][j][m];
              d = indj[g]; // 当該集団が受検していれば1，そうでなければ0
              nm = nm + nmm * d;
            }
            rr = 0;
            for(int g=0; g<ng; g++){
              rrr = rjm[g][j][m];
              d = indj[g]; // 当該集団が受検していれば1，そうでなければ0
              rr = rr + rrr * d;
            }
            double xm = Xm[m];
            P = c + (1-c)/(1+exp(-D*a*(xm-b)));
            Q = 1-P;
            // the element of Hessian Matrix
            // matrix[0,0]
            double daam = (P-c)*Q*(xm - b)*(xm - b)*D*D
            / ((1-c)*(1-c)*P*P);
            daam *= -nm * P*P + rr*c;
            daa += daam;
            // matrix[1,1]
            double dbbm = (P-c)*Q*a*a*D*D
              / ((1-c)*(1-c)*P*P);
            dbbm *= -nm*P*P + rr*c;
            dbb += dbbm;
            // matrix[2,2]
            double dccm = (rr/(nm*P) -1) - rr*Q/(P*P);
            dccm /= (1-c)*(1-c);
            dcc = dcc + dccm;
            //matrix[1,0], [0,1]
            double dabm = 1;
            dabm -= a*D*(xm-b)*Q/(1-c)*(rr/(nm*P)-1);
            dabm += D*a/(1-c)*(xm-b)*Q/P*(rr*c/P-rr);
            dabm *= -D /(1-c)*(P-c);
            dab += dabm;
            // matrix[2,0],[0,2]
            double dacm = -D/((1-c)*(1-c))*(xm-b)*(P-c)*Q/(P*P)*rr;
            dac += dacm;
            // matrix[2,1],[1,2]
            double dbcm = D/((1-c)*(1-c))*a*(P-c)*Q/(P*P)*rr;
            dbc += dbcm;

            // the elements of gradient vector
            double dad = D*(xm - b) * (rr - nm * P) * (P-c)
              / ((1-c) * P);
            da = da + dad;
            double dbd = -D*a*(rr - nm * P) * (P-c)
              / ((1-c) * P);
            db = db + dbd;
            double dcd = (rr - nm * P)
              /( (1-c) * P);
            dc = dc + dcd;
          } // 行列計算の下準備終了 // end of m

          if(Bayes==1){
            // the elements of gradient vector
            da += -1/a - (log(a)-mu_a)/(a*sigma_a*sigma_a);
            db += -(b-mu_b)/(sigma_b*sigma_b);
            dc += (alpha_c-2)/c - (beta_c-2)/(1-c);
            // Information matrix
            daa += 1/(a*a) - (1-log(a)+mu_a)/(a*a*sigma_a*sigma_a);
            dbb += -1/(sigma_b*sigma_b);
            dcc += (alpha_c-2)/(c*c) - (beta_c-2)/(1-c)*(1-c);
            //非対角要素は通常の二階偏微分と同じ。
          }
        }

        double a1 = 1.702/D;
        double b1 = 0;
        double c1 = 0;

        // 行列，ベクトルに代入
        if(model[j] == "3PL"){
          // determinant of Fisher's information matrix
          double det = 1/(daa*dbb*dcc + dab*dbc*dac + dac*dab*dbc - dac*dbb*dac - dab*dab*dcc - daa*dbc*dbc);

          // each element of inverse of Fisher's information matrix
          double inv_daa = (dbb*dcc-dbc*dbc) * det;
          double inv_dbb = (daa*dcc-dac*dac) * det;
          double inv_dcc = (daa*dbb-dab*dab) * det;
          double inv_dab = (dbc*dac-dab*dcc) * det;
          double inv_dac = (dab*dbc-dbb*dac) * det;
          double inv_dbc = (dab*dac-daa*dbc) * det;

          if(method == "Newton_Raphson"){
            a1 = a - inv_daa*da - inv_dab*db - inv_dac*dc;
            b1 = b - inv_dab*da - inv_dbb*db - inv_dbc*dc;
            c1 = c - inv_dac*da - inv_dbc*db - inv_dcc*dc;
          }
          if(method == "Fisher_Scoring") {
            a1 = a + inv_daa*da + inv_dab*db + inv_dac*dc;
            b1 = b + inv_dab*da + inv_dbb*db + inv_dbc*dc;
            c1 = c + inv_dac*da + inv_dbc*db + inv_dcc*dc;
          }

        } else if (model[j] == "2PL"){
          // determinant of Fisher's information matrix
          double det_Itj = 1/(daa * dbb - dab * dab);

          // each element of inverse of Fisher's information matrix
          double inv_daa = dbb * det_Itj;
          double inv_dbb = daa * det_Itj;
          double inv_dab = -dab * det_Itj;

          if(method == "Newton_Raphson"){
            a1 = a - inv_daa * da - inv_dab * db;
            b1 = b - inv_dab * da - inv_dbb * db;
          }
          if(method == "Fisher_Scoring") {
            a1 = a + inv_daa * da + inv_dab * db;
            b1 = b + inv_dab * da + inv_dbb * db;
          }

        } else {
          if(method == "Newton_Raphson"){
            b1 = b - db/dbb;
          }
          if(method == "Fisher_Scoring") {
            b1 = b + db/dbb;
          }

        }

        if(print >= 3) Rprintf("item %d -- abs a is %.7f, abs b is %.7f, abs c is %.7f in count %d\r", j+1, fabs(a-a1), fabs(b-b1), fabs(c-c1), count2);

        if( fabs(a1 - a) < eM && fabs(b1 - b) < eM && fabs(c1 - c) < eM ) conv2 = 1;

        // 次回の更新用に代入

        // discriminating parameter
        // aパラメタが0以下だったら，計算を中止
        if(a1 <= min_a || a1 > 10 || traits::is_nan<REALSXP>(fabs(a1 - a))){
          skip_para(j,0) += 1; // count up skip vector
          //t0m(j,0) = min_a;
          warning("'a' must be positive value. Error of item %i in %i time iteration.", j+1, count1);
          conv2 = 1;; // 次の項目パラメタの更新へスキップ
        } else {
          t0m(j,0) = a1;
        }


        // difficulty parameter

        if(fabs(b1) >= maxabs_b || traits::is_nan<REALSXP>(fabs(b1 - b))){
          skip_para(j,1) += 1; // count up skip vector
          if(b1 < 0){
            warning("'b' must be larger than %d. Error of item %i in %i time iteration.", maxabs_b, j+1, count1);
            //t0m(j,1) = -maxabs_b;
          } else if(b1 >= 0){
            warning("'b' must be smaller than %d. Error of item %i in %i time iteration.", maxabs_b, j+1, count1);
            //t0m(j,1) = maxabs_b;
          }
          conv2 = 1;; // 次の項目パラメタの更新へスキップ
        } else {
          t0m(j,1) = b1;
        }

        // guessing parameter
        if(c1 < 0 || c1 > 1 || traits::is_nan<REALSXP>(fabs(c1 - c))){
          //skip_para(j,2) += 1; // count up skip vector
          warning("'c' must be real value between 0 to 1. Error of item %i in %i time iteration.", j+1, count1);
          Rprintf("\nThe model of item %i is changed from '3PL' to '2PL'." , j+1);
          t0m(j,2) = 0;
          model[j] = "2PL"; // 2PLに変更する。
          conv2 = 1; // 次の項目パラメタの更新へスキップ
        } else {
          t0m(j,2) = c1;
        }

      }// end of Fisher scoring for item j

    }// end of Fisher scoring in a M step.


    // EMステップの収束判定


    int vv = 0;
    for(int j=0; j<nj; j++){
      if(model[j]=="1PL"){
        for(int p=1; p<2; p++){
          double p1 = t0m(j,p);
          double p2 = t0(j,p);
          double p3 = skip_para(j,p);
          if(fabs(p1 - p2) < eEM || p3>maxskip_j){// 項目パラメタの差の絶対値，と吹っ飛んだパラメタの抑制
            conv[vv] = 1;
          } else {
            conv[vv] = 0;
          }
          vv = vv + 1;
        }
      }else if(model[j]=="2PL"){
        for(int p=0; p<2; p++){
          double p1 = t0m(j,p);
          double p2 = t0(j,p);
          double p3 = skip_para(j,p);
          if(fabs(p1 - p2) < eEM || p3>maxskip_j){// 項目パラメタの差の絶対値，と吹っ飛んだパラメタの抑制
            conv[vv] = 1;
          } else {
            conv[vv] = 0;
          }
          vv = vv + 1;
        }
      }else if(model[j]=="3PL"){
        for(int p=0; p<3; p++){
          double p1 = t0m(j,p);
          double p2 = t0(j,p);
          double p3 = skip_para(j,p);
          if(fabs(p1 - p2) < eEM || p3>maxskip_j){// 項目パラメタの差の絶対値，と吹っ飛んだパラメタの抑制
            conv[vv] = 1;
          } else {
            conv[vv] = 0;
          }
          vv = vv + 1;
        }
      }
    }

    // update group distribution
    // 再計算

    // 母集団ごとの平均と標準偏差を計算
    LogicalVector conv4 (ng); // 平均と標準偏差の収束確認用

    if(thdist == "normal"){
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
            mean_j = mean_j + tempM; //
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
        mean[g] = mean_g/counti; // 全項目における母集団平均の，平均
        sd[g] = sd_g/counti; // 全項目における母集団標準偏差の，平均

        //Rcout << "mean: "<< mean[g]<<", sd: "<<sd[g]<<"\n";
      }
    } else if(thdist == "empirical"){ // 受検者の分点ごとの期待度数を用いて，次回のEステップの事前分布を計算する。
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
        if(fabs(mean[g] - mean_g/counti)<emu && fabs(sd[g] - sd_g/counti)<esd ) conv4[g] = 1;

        mean[g] = mean_g/counti; // 全項目における母集団平均の，平均
        sd[g] = sd_g/counti; // 全項目における母集団標準偏差の，平均
        if(print >= 1) Rcout << "mean: "<< mean[g]<<", sd: "<<sd[g]<<"\n";

        for(int m=0; m<N; m++){
          Um(m,g) /= counti;
        }

      } // end of g
    }


    double A = 1;
    double K = 0;

    if(fix == 1){// && thdist != "empirical"){
      A = sigma/sd[bg];
      K = mu - A * mean[bg];
    }

    // 平均0，標準偏差1になるように線形変換
    // 基準集団は１


    // 平均と標準偏差をキャリブレーション
    if(print == 1) Rprintf("\n----------------------------------");
    for(int g=0; g<ng; g++){
      mean[g] = A * mean[g] + K;
      sd[g] = A * sd[g];
      if (print >= 2) Rprintf("\n mean: %.5f, sd: %.5f ", mean[g],sd[g] );
      //Rcout <<"\n mean: "<< mean[g]<<", sd: "<<sd[g];
    }
    if(print == 1) Rprintf("\n----------------------------------");

    // 次回のEステップ用の項目パラメタもキャリブレーション

    for(int j=0; j<nj; j++){
      a = t0m(j,0); // a parameter
      b = t0m(j,1); // b parameter
      c=  t0m(j,2); // c parameter
      LogicalVector rm_j = 0;
      rm_j = j==rm_id;
      if(skip_para(j,0)>maxskip_j || skip_para(j,1)>maxskip_j || skip_para(j,2)>maxskip_j || is_true(any(rm_j))){
        // Mステップの更新をスキップした項目は，そのまま
        t0(j,0) = a;
        t0(j,1) = b;
        t0(j,2) = c;
        t0m(j,0) = a;
        t0m(j,1) = b;
        t0m(j,2) = c;
      } else {
        // きちんと更新をおこなった項目は，キャリブレーション
        if(model[j] != "1PL"){ // 1PLの場合は識別力の尺度調整を行わない
          t0(j,0) = a/A;
          t0m(j,0) = a/A;
        } else {
          t0(j,0) = a;
          t0m(j,0) = a;
        }
        t0(j,1) = A * b + K;
        t0(j,2) = c;
        t0m(j,1) = A * b + K;
        t0m(j,2) = c;
      }
    }

    if(print >= 2){
      Rprintf("\n Item Parameters in %d times iteration", count1);
      for(int j=0; j<nj; j++){
        //Rcout << "\n a: "<<t0m(j,0)<<", b: "<<t0m(j,1)<< ", c: "<<t0m(j,2);
        Rprintf("\n a: %.5f, b: %.5f, c: %.5f",t0m(j,0), t0m(j,1), t0m(j,2));
      }
      Rcout <<"\n MLL difference is "<<diff;
      Rprintf("\n------------------------------------------------------\n");

    }


    ////////////////////////////////////////////
    // 収束判定
    ////////////////////////////////////////////

    //  項目パラメタの差の絶対値，か周辺対数尤度の差の絶対値が一定数よりも下回ったら収束
    // もしくは最大繰り返し回数に達したら，終了。
    if((is_true(all(conv)))|| count1 == maxiter_em || diff < eMLL || is_true(all(conv4))){
      //Rcout << "eval EM conv\n";
      conv1 = 1; // 収束していたら繰り返しを打ち切り

      if(print>=1)Rcout <<"\nEstimating Standard Error";

      ///////////////////////
      // SE of estimation. //
      ///////////////////////


      for(int j=0; j<nj; j++){

        NumericVector indj = ind(_,j); // item index

        double a = t0m(j,0);
        double b = t0m(j,1);
        double c = t0m(j,2);
        // each element of Fisher's information matrix and gradient vector

        double daa,dbb,dcc,dab,dac,dbc;
        double da, db, dc;
        double d, rr, rrr, nm, nmm;
        double P,Q;
        // initialize
        daa=0;
        dbb=0;
        dcc=0;
        dab=0;
        dac=0;
        dbc=0;
        da=0;
        db=0;
        dc=0;
        P=0;
        Q=0;

        for(int m=0; m<N; m++){ // 行列計算

          nm = 0;
          for(int g=0; g<ng; g++){
            nmm = Nm[g][j][m];
            d = indj[g]; // 当該集団が受検していれば1，そうでなければ0
            nm = nm + nmm * d;
          }

          rr = 0;
          for(int g=0; g<ng; g++){
            rrr = rjm[g][j][m];
            d = indj[g]; // 当該集団が受検していれば1，そうでなければ0
            rr = rr + rrr * d;
          }

          double xm = Xm[m];

          P = c + (1-c)/(1+exp(-D*a*(xm-b)));
          Q = 1-P;

          // matrix[0,0]
          double daam = nm * (P-c) * (P-c) * Q * (xm - b) * (xm - b) * D * D
            / ((1-c) * (1-c) * P);
          daa = daa + daam;

          // matrix[1,1]
          double dbbm = nm * (P-c) * (P-c) * Q * a * a * D * D
            / ((1-c) * (1-c) * P);
          dbb = dbb + dbbm;

          // matrix[2,2]
          double dccm = nm * Q
            / (P * (1-c) * (1-c));
          dcc = dcc + dccm;

          //matrix[1,0], [0,1]
          double dabm = -nm * (P-c) * (P-c) * Q * a * (xm - b) * D * D
            / ((1-c) * (1-c) * P);
          dab = dab + dabm;

          // matrix[2,0],[0,2]
          double dacm = nm * (P-c) * Q * (xm - b) * D
            / ((1-c) * (1-c) * P);
          dac = dac + dacm;

          // matrix[2,1],[1,2]
          double dbcm = -nm * (P-c) * Q * a * D
            / ((1-c) * (1-c) * P);
          dbc = dbc + dbcm;

          // the elements of gradient vector
          double dad = D*(xm - b) * (rr - nm * P) * (P-c) // cが消える
            / ((1-c) * P);
          da = da + dad;

          double dbd = -D*a*(rr - nm * P) * (P-c)
            / ((1-c) * P);
          db = db + dbd;

          double dcd = (rr - nm * P)
            / ((1-c) * P);
          dc = dc + dcd;

        } // 行列計算の下準備終了 // end of m


        // 行列，ベクトルに代入
        if(model[j] == "3PL"){
          // determinant of Fisher's information matrix
          double det = 1/(daa*dbb*dcc+dab*dbc*dac+dac*dab*dbc-dac*dbb*dac-dab*dab*dcc-daa*dbc*dbc);
          // each element of inverse of Fisher's information matrix

          double inv_daa = (dbb*dcc-dbc*dbc) * det;
          double inv_dbb = (daa*dcc-dac*dac) * det;
          double inv_dcc = (daa*dbb-dab*dab) * det;

          SE(j,0) = std::sqrt(inv_daa); //
          SE(j,1) = std::sqrt(inv_dbb);
          SE(j,2) = std::sqrt(inv_dcc);

        } else if (model[j] == "2PL"){
          // determinant of Fisher's information matrix
          double det_Itj = 1/(daa * dbb - dab * dab);

          // each element of inverse of Fisher's information matrix
          double inv_daa = dbb * det_Itj;
          double inv_dbb = daa * det_Itj;

          SE(j,0) = std::sqrt(inv_daa); //
          SE(j,1) = std::sqrt(inv_dbb);

        } else {
          SE(j,1) = std::sqrt(1/dbb);
        }

      }// end of estimate SE for item j

      Rcout << "\nEM algorithm has been converged.\nTotal iteration time is " << count1 << "\n";

    }


  } // end of EM step



  // 推定母集団分布の計算
  // 初期値に一様分布をもちいた場合，Easy Estimationの計算結果とおおむね一致するはず。
  if(EM_dist == 1) Rprintf("Start calculating estimated population distribution.\n");

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

  if(EM_dist == 0) conv3 = 1;
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
      if(fabs(mean_pop[g] - mean_g/counti)<emu && fabs(sd_pop[g] - sd_g/counti)<esd ) conv0[g] = 1;

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


  // Item Fit index
  // reference: Mayekawa, S. (1991) Estimation of Item parameter.
  // Xm; node
  //
  boost::multi_array <double, 3> EMfit (boost::extents[ng][nj][N]);

  if( print > 0) Rcout<<"expected frequency of subjects in each nodes calculation.\n";
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
        Nm[g][j][m]= k;
      }
    }
  }

  if( print > 0) Rcout<<"expected freqency of correct response number in each nodes calculation.\n";

  double h,gg,hh,u;
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

  for(int g=0; g<ng; g++){
    for(int j=0; j<nj; j++){
      double a = t0(j,0);
      double b = t0(j,1);
      double c = t0(j,2);
      //
      for(int m=0; m<N; m++){
        double tpjm = c+(1.0-c)/(1.0+exp(-D*a*(Xm[m]-b))); // theoritical item correct response rate
        double rpjm = rjm[g][j][m]/Nm[g][j][m]; // real item passing rate
        double ejm = rpjm - tpjm; // x-t
        EMfit[g][j][m] = ejm;
      }// end of m
    }// end of j(item)
  }// end of g





  // アウトプットファイルの調整

  NumericMatrix ms = cbind(mean,sd);
  DataFrame Para = DataFrame::create(Named("Item")=Item, Named("a") = t0(_,0), Named("b") = t0(_,1), Named("c") = t0(_,2), Named("model") = model);
  DataFrame SE_d = DataFrame::create(Named("Item")=Item, Named("a") = SE(_,0), Named("b") = SE(_,1), Named("c") = SE(_,2), Named("model") = model);

  List res = List::create(_["para"] = Para, _["SE"] = SE_d, _["initial"] = initial,
                          _["theta.dist"] = dist, _["ms"] = ms,_["mean"] = mean, _["sd"] = sd,
                          _["population_dist"] = unif_dist, _["population_mean"] = mean_pop, _["population_sd"] = sd_pop,
                          _["MLL"] = MLL, _["conv"] = conv, _["count1"] = count1, _["count2"] = count2,
                          _["skip_para"] = skip_para, _["rm_n"] = rm_n, _["rm_id"] = rm_id ,
                          _["p_bis"] = r, _["passing_rate"] = p, _["itemfit_EM"] = EMfit
  );

  return res;

}


