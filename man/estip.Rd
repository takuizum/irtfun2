\name{estip}
\alias{estip}
\title{
  A function estimate item parameter using marginal maximum likelihood estimation method via EM algorithm
}
\usage{
estip(x, model = "2PL", N = 31L, eMLL = 1e-6, eEM = 1e-4, eM = 1e-3,
  emu = 1e-4, fc0 = 2L, ng = 1L, gc0 = 2L, D = 1.702, fix = 1L,
  print = 0L, ic = 1/5, max = 6.0, min = -6.0, mu = 0, sigma = 1,
  Bayes = 0L, mu_a = 0, sigma_a = 1, mu_b = 0, sigma_b = 2,
  mu_c = 4/13, w_c = 13, min_a = 0.1, maxabs_b = 20, maxiter_em = 200L,
  maxiter_j = 20L, maxskip_j = 5L, rm_list = as.character( c("NONE")),
  thdist = "normal", e_ell = 1L, EM_dist = 1L)
}
\arguments{
\item{x}{an item response data which class is data.frame object.}

\item{model}{Character.U can select which one, "1PL","2PL","3PL".}

\item{N}{the number of nodes in integration}

\item{eMLL}{a convergence criteria(CC) forf marginal log likelihood.}

\item{eEM}{a CC in EM cycle.}

\item{eM}{a CC in M step.}

\item{emu}{a CC for population distribution mean.}

\item{fc0}{a column of first item response.}

\item{ng}{the number of groups.}

\item{gc0}{a column of group.}

\item{D}{a scaling constant.}

\item{fix}{if 1, fix population distribution mean and sigma each EM cycle.}

\item{print}{How much information you want to display? from 1 to 3. The larger, more information is displayed.}

\item{ic}{a initial value for guessing parameter.}

\item{max}{a maximum value of theta in integration.}

\item{min}{a minimum value of theta in integration.}

\item{mu}{a hyperparameter for theta dist.}

\item{sigma}{same as above.}

\item{Bayes}{If 1, marginal Bayesian estimation runs.}

\item{mu_a}{a hyperparameter of log normal dist for slope parameter.}

\item{sigma_a}{same as above.}

\item{mu_b}{a hyperparameter of normal dist for location parameter.}

\item{sigma_b}{same as above.}

\item{mu_c}{a hyperparameter for lower asymptote parameter.}

\item{w_c}{same as above.}

\item{min_a}{a minimum value of slope parameter.}

\item{maxabs_b}{a maximum absolute value of location parameter.}

\item{maxiter_em}{a maximul iteration time for EM cycle.}

\item{maxiter_j}{a maximul iteration time for Newton Raphton in M step.}

\item{maxskip_j}{Dont use.}

\item{rm_list}{a vector of item U want to remove for estimation. NOT list.}

\item{thdist}{Dont use.}

\item{e_ell}{If 1, use a CC of expected log likelihood.}

\item{EM_dist}{If 1, calculate esimated population distribution by EM argorithm.}

}
\description{
  U can select 'model' option. 1PL, 2PL, 3PL, Bayesian estimation is available.
}
\value{
\item{para}{a df of estimated item parameter.}

\item{SE}{a df of standard error of estimation.}

\item{initial}{initial values were used in estimation}

\item{theta.dist}{theta distribution calculated in last EM cycle.}

\item{ms}{mean and sigma calculated by above distribution file.}

\item{mean}{same as above.}

\item{sd}{same as above.}

\item{population_dist}{estimated population distribution via EM algorithm. Initial distribution is uniform.}

\item{population_mean}{mean calculated by above distribution file.}

\item{population_sd}{sd calculated by above distribution file.}

\item{MLL}{a log of marginal log likelihood in all EM cycle.}

\item{conv}{logical.}

\item{count1}{???}

\item{count2}{???}

\item{skip_para}{the item name skipped in EM cycle.}

\item{rm_n}{the number of remove items.}

\item{rm_id}{a removed item list.}

\item{p_bis}{point biserial correlation.}

\item{raw_score}{a raw score matrix.}

\item{passing_rate}{passing rates used in initial value.}

}

\author{
Takumi Shibuya, sep10.taku.izum@gmail.com

Maintainer: same as above.
}

