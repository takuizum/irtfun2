# irtfun2
------
## This package provides some usefule function.
このパッケージではIRTの分析に使える関数を使用できます。一部の関数を使用するためにはRtoolsなどのコンパイラをPCにインストールする必要があります。(最終更新：2018/10/09)

This package contains some useful function for IRT analysis. It is necessary to install c++ compiler(Rtool etc.) in your PC before use "estip" and "estheta".

### バージョン0.4.0での変更点。esthetaのPVsオプションで，C++による高速なリジェクションサンプリングを可能にしました。argument'sampling_engine="Cpp"'で実行可能です。

### estip
"estip" is a function for marginal maximum likelihood estimation of item parameter. For binary response data only. MML via EM is a standard estimation method in IRT item parameter estimaiton.

この関数はEMアルゴリズムを用いた周辺最尤推定法により項目パラメタを推定するためのものです。

### estheta
"estheta" is a function for estimating MLE, EAP, MAP and PVs.　Rejection sampling method is implemented for PVs sub routine.

この関数では'estip'などを使用して推定した項目パラメタを使って，受験者特性値\thetaを推定することができます。使用できるオプションにはMLE（最尤推定法），EAP（事後分布期待値），MAP（事後分布最頻値），PVs（推算値）があります。推算値の計算には棄却サンプリング法を使用しています。
