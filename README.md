# irtfun2
------
## This package provides some usefule function.
このパッケージではIRTの分析に使える関数を使用できます。

### estip
"estip" is a function for marginal maximum likelihood estimation of item parameter. For binary response data only. MML via EM is a standard estimation method in IRT item parameter estimaiton.

この関数はEMアルゴリズムを用いた周辺最尤推定法により項目パラメタを推定するためのものです。

### estheta
"estheta" is a function for estimating MLE, EAP, MAP and PVs.

この関数では'estip'などを使用して推定した項目パラメタを使って，受験者特性値\thetaを推定することができます。使用できるオプションにはMLE（最尤推定法），EAP（事後分布期待値），MAP（事後分布最頻値），PVs（推算値）があります。
