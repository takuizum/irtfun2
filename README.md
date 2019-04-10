# irtfun2
------
## This package provides some usefule function fir IRT.
このパッケージではIRTの分析に使える関数を使用できます。一部の関数を使用するためにはRtoolsなどのコンパイラをPCにインストールする必要があります。(最終更新：2019/04/09)

This package contains some useful function for IRT analysis. It is necessary to install c++ compiler(Rtool etc.) in your PC.

### バージョン0.4.0での変更点。esthetaのPVsオプションで，C++による高速なリジェクションサンプリングを可能にしました。argument'sampling_engine="Cpp"'で実行可能です。

### バージョン0.6.6.1での変更点。ベイズ周辺最尤推定法を実行する機能をestipに追加。Bayes=1で実行できます。一般項目反応理論(General IRT)のパラメタ推定関数estGipを追加しました。GIRTについては芝（1991）「項目反応理論」の2章は6章を参照してください。

### バージョン0.6.7.3での変更点。正則化周辺最尤推定法や複数の反復計算手法に対応した項目パラメタ推定関数estip2を実装しました。estGipを多母集団モデルに適応させました。その他，細かいバグを修正しました。

### バージョン0.6.7.5での変更点。esthetaのMAP推定に関するバグを修正。(2019/04/10)

### estip & estip2
"estip" is a function for marginal maximum likelihood (MML) estimation of item parameter. For binary response data only. MML via EM is a standard estimation method in IRT item parameter estimaiton.

この関数はEMアルゴリズムを用いた周辺最尤推定法により項目パラメタを推定するためのものです。
```{r}
devtools::install_github("takuizum/irtfun2", dependency = TRUE)
library(irtfun2)
# estimation item parameter using simulation data in irtfun2 package.
# 2 parameter logistic model
res <- estip2(sim_dat_2, model = "2PL", fc = 2)
res$para # estimated item parametes(data.frame)
res$se # standard error(data.frame)
```

### estheta
"estheta" is a function for estimating MLE, EAP, MAP and PVs.　Rejection sampling method is implemented for PVs sub routine.

この関数では'estip'などを使用して推定した項目パラメタを使って，受験者特性値\thetaを推定することができます。使用できるオプションにはMLE（最尤推定法），EAP（事後分布期待値），MAP（事後分布最頻値），PVs（推算値）があります。推算値の計算には棄却サンプリング法とスライスランプリングを使用しています。

### estGip
この関数は一般項目反応モデルに基づいて項目パラメタと受検者の能力パラメタ，能力パラメタの標準偏差($\phi$)などを推定します。現時点では推定の標準誤差や適合度には対応していませんが，今後対応予定です。

筆者の英語とパッケージ作成の技術力が乏しいせいで，マニュアルを見ただけではどう使って良いかがよく分からない仕様になっています。興味のある方は直接筆者（澁谷，sep10.taku.izum(at)gmail.com）にご連絡ください。(at)を＠に変えてください。
