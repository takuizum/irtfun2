# GRM analysis tools

# response probability----
pgrm <- function(theta, a, b, k, D = 1.0){
  K <- length(b) + 1 # n of category(n of b parameter + 1)
  if(k == 1){
    p1 <- 1/(1 + exp(-D * a * (theta - b[k])))
    p0 <- 1 # category 0
  } else if (k == K){
    p1 <- 0 # category K + 1
    p0 <- 1/(1 + exp(-D * a * (theta - b[k -1])))
  } else if (k < K && k > 1) {
    p0 <- 1/(1 + exp(-D * a * (theta - b[k-1])))
    p1 <- 1/(1 + exp(-D * a * (theta - b[k])))
  } else { # NA
    p0 <- 1
    p1 <- 0
  }
  p0-p1
}

# bt <- c(-1.5, 0, 1)
# pgrm(0, a = 2, b = bt, k = 1)
# tibble(theta = c(-4:4)) %>%
#   ggplot(aes(x = theta))+
#   stat_function(fun = pgrm, args = list(a = 2, b = bt, k = 1, D = 1), colour = 2)+
#   stat_function(fun = pgrm, args = list(a = 2, b = bt, k = 2, D = 1), colour = 3)+
#   stat_function(fun = pgrm, args = list(a = 2, b = bt, k = 3, D = 1), colour = 4)+
#   stat_function(fun = pgrm, args = list(a = 2, b = bt, k = 4, D = 1), colour = 6)

# return cumulative probability----
pgrm_cum <- function(theta, para, D = 1.0){
  a <- para[1]
  b <- para[-1]
  K <- length(b) + 1
  p <- numeric(K)
  for(k in 1:K){
    if(k == 1){
      p1 <- 1/(1 + exp(-D * a * (theta - b[k])))
      p0 <- 1 # category 0
    } else if (k == K){
      p1 <- 0 # category K + 1
      p0 <- 1/(1 + exp(-D * a * (theta - b[k -1])))
    } else if (k < K && k > 1){
      p1 <- 1/(1 + exp(-D * a * (theta - b[k])))
      p0 <- 1/(1 + exp(-D * a * (theta - b[k-1])))
    } else {
      p0 <- 0 # それ以外の反応はあり得ないので，確率0
      p1 <- 0
    }
    p[k] <- p0 - p1
  }
  sample(c(1:(length(b) + 1)), 1, prob = p)
}
# # hist
# result <- numeric(10000)
# for(t in 1:10000) result[t] <- pgrm_cum(0, c(2, bt))
# hist(result)


# return response ----
# return a list of item parameters
# 引数の指定がうまくいっていないので改良の余地あり
grm_para_gen <- function(category_vector, b_dist = "norm", a_dist = "lnorm",
                         args_b_dist = list(0, 1), args_a_dist = list(0, 1),
                         item_name = NULL, index = "Item"){
  assign("rand_dist1", get(paste0("r", b_dist), mode = "function"))
  assign("rand_dist2", get(paste0("r", a_dist), mode = "function"))
  if(a_dist == "beta"){
    res <- purrr::map(as.list(category_vector - 1),
                      function(x){c(rand_dist2(1, args_a_dist[[1]], args_a_dist[[2]])*1.5+0.4,
                                    sort(rand_dist1(n = x, args_b_dist[[1]], args_b_dist[[2]])))})
  } else {
    res <- purrr::map(as.list(category_vector - 1),
                      function(x){c(rand_dist2(1, args_a_dist[[1]], args_a_dist[[2]]),
                                    sort(rand_dist1(n = x, args_b_dist[[1]], args_b_dist[[2]])))})
  }
  if(is.null(item_name)){
    item_name <- formatC(c(1:length(category_vector)), width = 3, flag = 0)
    item_name <- as.vector(apply(matrix(index, ncol = 1), 1, paste0, item_name))
  }
  names(res) <- item_name
  return(res)
}

# parameter check
# tibble(theta = c(-4:4)) %>%
#   ggplot(aes(x = theta))+
#   stat_function(fun = pgrm, args = list(a = test[[1]][1], b = test[[1]][-1], k = 1, D = 1), colour = 2)+
#   stat_function(fun = pgrm, args = list(a = test[[1]][1], b = test[[1]][-1], k = 2, D = 1), colour = 3)

# test <- grm_para_gen(c(2,2,3))
# test2 <- pgrm_cum(0, test[[1]])

# return response (double) ----
sub_grm <- function(theta, para, D = 1.0){
  purrr::map(para, pgrm_cum, theta = theta, D = D)
}


# simulation data.frame generation(variable type is 'tibble') ----
sim_grm_gen <- function(theta, para, D = 1.0, min_resp = 1){
  ID <- 1:length(theta) # subject ID
  if(min_resp == 1){
    res <- as.list(theta) %>% purrr::map_df(sub_grm, para = para, D = D) %>%
      purrr::map_df(as.integer) %>%
      tibble::add_column(ID = ID, .before = 1)
  } else if (min_resp == 0){
    res <- as.list(theta) %>% purrr::map_df(sub_grm, para = para, D = D) %>%
      purrr::map_df(magrittr::subtract, e2 = 1) %>%
      purrr::map_df(as.integer) %>%
      tibble::add_column(ID = ID, .before = 1)
  }
}
# # function check
# set.seed(11)
# theta_sample <- rnorm(10000)
# item_cate <- c(rep(5, 5),rep(3, 5))
# item_sample <- grm_para_gen(item_cate, args_b_dist = list(-1, 2))
# dat1 <- sim_grm_gen(theta_sample, item_sample, min_resp = 0) # a little bit slow
# dat1
#
# # category check
# dat1[,-1] %>%
#   map(unique) %>%
#   map(sort)
# dat1[,-1] %>%
#   map(table)


# ?ltm::grm
#
# grm(data.matrix(dat_2[,-1]))
#
# item_sample
#
# # output for essy estimation(Kumagai, 2009)
# dat_2$ID <- formatC(c(1:50000), width = 5, flag = 0)
# write.table(dat_2, file = "grm1.dat", row.names = F, col.names = F, sep = "")


# merging function ----
merge_category2 <- function(x){
  if(x == 1){
    1
  } else if(x > 1 && x <5){
    2
  } else if(x > 4 && x < 8){
    3
  } else if(x > 7 && x < 11){
    4
  } else if(x > 10 && x < 14){
    5
  } else if(x > 13 && x < 17){
    6
  } else {
    7
  }
}

merge_category <- function(x){
  apply(matrix(x), 1, merge_category2)
}

# test <- c(sample(c(1:17), 10000, replace = T))
# system.time(
#   merge_category(test)
# )
# resp %>% purrr::map_df(merge_category)

#----imputation merged category parameter function----
cal_imp_n <- function(k0, k1){
  imp_n <- floor((k1-1) / (k0-1)) -1 # imputetion N to substitute between 2 adjacent categories
  total_category <- (imp_n+1) * (k0-1) + 1
  imp_ind1 <- rep(imp_n, (k0-1))
  # 割り切れない数にあたった場合，低いカテゴリから順にぶち込んでいく
  if(total_category != k1){
    diff_n <- abs(total_category - k1)
    for(k in 1:diff_n){
      imp_ind1[k] <- imp_ind1[k] + 1
    }
  }
  imp_ind2 <- numeric(k1)
  t <- 0 # counter
  imp_ind2[1] <- 1
  for(i in imp_ind1){
    t <- t + 1
    imp_ind2[t] <- 1
    t <- t + i
  }
  imp_ind2[k1] <- 1
  list(ind1 = imp_ind1, ind2 = imp_ind2 %>% as.logical())
}

imp_ind <- cal_imp_n(3, 5)

# calculate imputed category parameter, based on uniform distribution
imp_category_sub <- function(para, imp_ind){
  K <- length(imp_ind$ind2) # total length
  res <- numeric(K+1) # to involve discrimination parameter
  res[1] <- para[1] # a para
  i <- 1
  for(k in 1:K){
    if(imp_ind$ind2[k]){
      i <- i+1
      res[k+1] <- para[i]
    } else if((i+1) <= length(para)) {
      b1 <- para[i]
      b2 <- para[i+1]
      distance <- (b2-b1)/(imp_ind$ind1[i-1] + 1)
      res[k+1] <- res[k+1] + b1
      res[(k+1):(k+imp_ind$ind1[i-1])] <- res[(k+1):(k+imp_ind$ind1[i-1])] + distance
    }
  }
  res
}

# calculate imputed category parameter, based on empirical distribution
imp_category_sub2 <- function(para, imp_ind, distribution){
  K <- length(imp_ind$ind2) # total length
  res <- numeric(K+1) # to involve discrimination parameter
  res[1] <- para[1] # a para
  # increasement of b parameter
  i <- 1
  increasement <- numeric(length(imp_ind$ind1))
  for(j in 1:length(imp_ind$ind1)){
    denominator <- sum(distribution[seq(i, length.out = imp_ind$ind1[j] + 1)])
    numerator <- abs(para[j+1] - para[j+2])
    increasement[j] <- numerator/denominator
  }
  # imputation
  i <- 1
  for(k in 1:K){
    if(imp_ind$ind2[k]){
      i <- i+1 # index for original b parameter
      res[k+1] <- para[i] # insert original b parameter to 'res' vector(res[1] is a parameter)
    } else if((i+1) <= length(para)) { # not to reach K + 1(not k + 1!) in original parameter vector
      b1 <- para[i]
      b2 <- para[i+1]
      inc <- increasement[i-1]
      # このあとは適切なヒストグラムの度数を求めて，それをincrasementに乗じて，b1に足せばOK
      distance <- (b2-b1)/(imp_ind$ind1[i-1] + 1)
      res[k+1] <- res[k+1] + b1
      res[(k+1):(k+imp_ind$ind1[i-1])] <- res[(k+1):(k+imp_ind$ind1[i-1])] + distance
    }
  }
  res
}

test_cate <- c(rep(1,100), rep(2,500), rep(3,1500), rep(4,500), rep(5,200), rep(6, 50), rep(7, 150))
test_table <- table(test_cate)
test_para <- c(1, -1, 0, 1)
imp_ind <- cal_imp_n(3, 6)

# test <- est_para_7[[1]][-6]
# test <-  test$Item001[-2,]
#
# imp_category_sub(test, imp_ind)

# imputation function
imp_category <- function(coef_object, K = 16){
  J <- length(coef_object)
  para <- coef_object[-J] %>% # remove estimated group parameter in list
    map(function(x)x["pars",]) # remove SE row
  # the number of category parameter of pre imputation(merged)
  k0 <- para %>% map(function(x){length(x)-1})
  # the target number of category parameter
  k1 <- para %>% map(function(x){K}) # K in not the number of categories but category difficulty parameters
  imp_ind <- map2(k0, k1, cal_imp_n)
  map2(para, imp_ind, imp_category_sub)
}

# imputation function that is based on subject frequency observed in merged category
imp_category <- function(coef_object, reponse_data, K = 16){
  #----NOW CODING
}

#----MAP theta estimation algorithm----

# Log likelihood
LLgrm_sub <- function(u, theta, para, D){
  a <- para[1]
  b <- para[-1]
  p1 <- log(pgrm(theta = theta, b = b, a = a, k = u, D = D))
}
LLgrm <- function(u, theta, para, D = 1.0, mu = 0, sigma = 1){
  ul <- as.list(u)
  lp <- purrr::map2_df(ul, para, LLgrm_sub, theta = theta, D = D)
  sum(lp) + log(dnorm(theta, mu, sigma))
}
LLgrm_apply <- function(u, para, D = 1.0, mu = 0, sigma = 1,  optimise_range = c(-8,8)){
  opt <- optimise(LLgrm, interval = optimise_range, u = u, para = para, maximum = T, mu = mu, sigma = sigma)
  opt$maximum
}

ext_par <- function(x){x["par",]} # extract only item parameter
# test <- est_para_17[[1]][-6] %>% map(pluck, ext_par) # map & pluck
modif_mirt_para <- function(para, n_of_items){
  tmp <- para[-(1+n_of_items)]
  tmp %>% map(pluck, ext_par)
}

# # optimize
# theta_sample[[1]] %>% head
# map1 <- apply(resp[[1]][,-1], 1, LLgrm_apply, para = test, D = 1.0, sigma = 3) # map: prior is N(0, 1)
#
# plot(theta_sample_list[[1]], map1, pch = 20, cex = 0.1)
# cor(theta_sample[[1]], map1)

# # validation
# # set.seed(old_seed[1])
# theta_sample <- rnorm(5000, 1, 1) # scale metric
# item_cate <- rep(5, 100)
# item_sample <- grm_para_gen(item_cate, a_dist = "beta", args_b_dist = list("mean" = 0.5, "sd" = 2),
#                             args_a_dist = list("shape1" = 2,"shape2" = 2))
# dat1 <- sim_grm_gen(theta_sample, item_sample, min_resp = 1) # a little bit slow
# dat1[,-1] %>% map(table) %>% map(length) %>% unlist %>% min
# mod <- mirt.model('B = 1-100
#                    START = (GROUP, MEAN_1, 1)')
# fit <- mirt(dat1[,-1], mod, itemtype = "graded", SE = T)
# para1 <- coef(fit, IRTpars = T, printSE = T)
# para2 <- modif_mirt_para(para1, 100)
# # map estimate
# map2 <- apply(dat1[,-1], 1, LLgrm_apply, para = para2, D = 1.0, sigma = 3)
# # visualization
# plot(theta_sample, map2, pch = 20, cex = 0.1)
# cor(theta_sample, map2)
#
# # fpd
# fpdgrm <- function(){
#
# }
# fpdgrm_sub <- function(u, theta, para, D){
#   a <- para[1]
#   b <- para[-1]
#   p1 <- pgrm(theta, b, a, u, D = D)
#   p0 <- pgrm(theta, b, a, u-1, D = D)
#
# }
# # spd?
#
# # fi?


# wrapper function
estgrmtheta <- function(dat, para, fc = 2, gc = 0, IDc = 1, D = 1.0, mu = 0, sigma = 1, optimise_range = c(-8,8) ){
  if(IDc == 0){
    id <- c(1:nrow(dat1))
  } else {
    id <- dat[IDc,]
  }
  xall <- dat[,fc:ncol(dat)]
  if(gc == 0){
    group <- rep(1, nrow(xall))
  } else {
    group <- dat[,gc]
  }
  cat("START MAP ESTIMATION!\n")
  map1 <- apply(xall, 1, LLgrm_apply, para = para, D = D, mu = mu, sigma = sigma, optimise_range = optimise_range)
  tibble::tibble(ID = id, GROUP = group, SCORE = rowSums(xall, na.rm = T), MAP = map1)
}
# estgrmtheta(resp[[1]][,-1], test, fc = 1, IDc = 0)

# GRM analysis tools

# response probability----
pgrm <- function(theta, a, b, k, D = 1.0){
  K <- length(b) + 1 # n of category(n of b parameter + 1)
  if(k == 1){
    p1 <- 1/(1 + exp(-D * a * (theta - b[k])))
    p0 <- 1 # category 0
  } else if (k == K){
    p1 <- 0 # category K + 1
    p0 <- 1/(1 + exp(-D * a * (theta - b[k -1])))
  } else if (k < K && k > 1) {
    p0 <- 1/(1 + exp(-D * a * (theta - b[k-1])))
    p1 <- 1/(1 + exp(-D * a * (theta - b[k])))
  } else { # NA
    p0 <- 1
    p1 <- 0
  }
  p0-p1
}

# bt <- c(-1.5, 0, 1)
# pgrm(0, a = 2, b = bt, k = 1)
# tibble(theta = c(-4:4)) %>%
#   ggplot(aes(x = theta))+
#   stat_function(fun = pgrm, args = list(a = 2, b = bt, k = 1, D = 1), colour = 2)+
#   stat_function(fun = pgrm, args = list(a = 2, b = bt, k = 2, D = 1), colour = 3)+
#   stat_function(fun = pgrm, args = list(a = 2, b = bt, k = 3, D = 1), colour = 4)+
#   stat_function(fun = pgrm, args = list(a = 2, b = bt, k = 4, D = 1), colour = 6)

# return cumulative probability----
pgrm_cum <- function(theta, para, D = 1.0){
  a <- para[1]
  b <- para[-1]
  K <- length(b) + 1
  p <- numeric(K)
  for(k in 1:K){
    if(k == 1){
      p1 <- 1/(1 + exp(-D * a * (theta - b[k])))
      p0 <- 1 # category 0
    } else if (k == K){
      p1 <- 0 # category K + 1
      p0 <- 1/(1 + exp(-D * a * (theta - b[k -1])))
    } else if (k < K && k > 1){
      p1 <- 1/(1 + exp(-D * a * (theta - b[k])))
      p0 <- 1/(1 + exp(-D * a * (theta - b[k-1])))
    } else {
      p0 <- 0 # それ以外の反応はあり得ないので，確率0
      p1 <- 0
    }
    p[k] <- p0 - p1
  }
  sample(c(1:(length(b) + 1)), 1, prob = p)
}
# # hist
# result <- numeric(10000)
# for(t in 1:10000) result[t] <- pgrm_cum(0, c(2, bt))
# hist(result)


# return response ----
# return a list of item parameters
# 引数の指定がうまくいっていないので改良の余地あり
grm_para_gen <- function(category_vector, b_dist = "norm", a_dist = "lnorm",
                         args_b_dist = list(0, 1), args_a_dist = list(0, 1),
                         item_name = NULL, index = "Item"){
  assign("rand_dist1", get(paste0("r", b_dist), mode = "function"))
  assign("rand_dist2", get(paste0("r", a_dist), mode = "function"))
  if(a_dist == "beta"){
    res <- purrr::map(as.list(category_vector - 1),
                      function(x){c(rand_dist2(1, args_a_dist[[1]], args_a_dist[[2]])*1.5+0.4,
                                    sort(rand_dist1(n = x, args_b_dist[[1]], args_b_dist[[2]])))})
  } else {
    res <- purrr::map(as.list(category_vector - 1),
                      function(x){c(rand_dist2(1, args_a_dist[[1]], args_a_dist[[2]]),
                                    sort(rand_dist1(n = x, args_b_dist[[1]], args_b_dist[[2]])))})
  }
  if(is.null(item_name)){
    item_name <- formatC(c(1:length(category_vector)), width = 3, flag = 0)
    item_name <- as.vector(apply(matrix(index, ncol = 1), 1, paste0, item_name))
  }
  names(res) <- item_name
  return(res)
}

# parameter check
# tibble(theta = c(-4:4)) %>%
#   ggplot(aes(x = theta))+
#   stat_function(fun = pgrm, args = list(a = test[[1]][1], b = test[[1]][-1], k = 1, D = 1), colour = 2)+
#   stat_function(fun = pgrm, args = list(a = test[[1]][1], b = test[[1]][-1], k = 2, D = 1), colour = 3)

# test <- grm_para_gen(c(2,2,3))
# test2 <- pgrm_cum(0, test[[1]])

# return response (double) ----
sub_grm <- function(theta, para, D = 1.0){
  purrr::map(para, pgrm_cum, theta = theta, D = D)
}


# simulation data.frame generation(variable type is 'tibble') ----
sim_grm_gen <- function(theta, para, D = 1.0, min_resp = 1){
  ID <- 1:length(theta) # subject ID
  if(min_resp == 1){
    res <- as.list(theta) %>% purrr::map_df(sub_grm, para = para, D = D) %>%
      purrr::map_df(as.integer) %>%
      tibble::add_column(ID = ID, .before = 1)
  } else if (min_resp == 0){
    res <- as.list(theta) %>% purrr::map_df(sub_grm, para = para, D = D) %>%
      purrr::map_df(magrittr::subtract, e2 = 1) %>%
      purrr::map_df(as.integer) %>%
      tibble::add_column(ID = ID, .before = 1)
  }
}
# # function check
# set.seed(11)
# theta_sample <- rnorm(10000)
# item_cate <- c(rep(5, 5),rep(3, 5))
# item_sample <- grm_para_gen(item_cate, args_b_dist = list(-1, 2))
# dat1 <- sim_grm_gen(theta_sample, item_sample, min_resp = 0) # a little bit slow
# dat1
#
# # category check
# dat1[,-1] %>%
#   map(unique) %>%
#   map(sort)
# dat1[,-1] %>%
#   map(table)


# ?ltm::grm
#
# grm(data.matrix(dat_2[,-1]))
#
# item_sample
#
# # output for essy estimation(Kumagai, 2009)
# dat_2$ID <- formatC(c(1:50000), width = 5, flag = 0)
# write.table(dat_2, file = "grm1.dat", row.names = F, col.names = F, sep = "")


# merging function ----
merge_category2 <- function(x){
  if(x == 1){
    1
  } else if(x > 1 && x <5){
    2
  } else if(x > 4 && x < 8){
    3
  } else if(x > 7 && x < 11){
    4
  } else if(x > 10 && x < 14){
    5
  } else if(x > 13 && x < 17){
    6
  } else {
    7
  }
}

merge_category <- function(x){
  apply(matrix(x), 1, merge_category2)
}

# test <- c(sample(c(1:17), 10000, replace = T))
# system.time(
#   merge_category(test)
# )
# resp %>% purrr::map_df(merge_category)

#----imputation merged category parameter function----
cal_imp_n <- function(k0, k1){
  imp_n <- floor((k1-1) / (k0-1)) -1 # imputetion N to substitute between 2 adjacent categories
  total_category <- (imp_n+1) * (k0-1) + 1
  imp_ind1 <- rep(imp_n, (k0-1))
  # 割り切れない数にあたった場合，低いカテゴリから順にぶち込んでいく
  if(total_category != k1){
    diff_n <- abs(total_category - k1)
    for(k in 1:diff_n){
      imp_ind1[k] <- imp_ind1[k] + 1
    }
  }
  imp_ind2 <- numeric(k1)
  t <- 0 # counter
  imp_ind2[1] <- 1
  for(i in imp_ind1){
    t <- t + 1
    imp_ind2[t] <- 1
    t <- t + i
  }
  imp_ind2[k1] <- 1
  list(ind1 = imp_ind1, ind2 = imp_ind2 %>% as.logical())
}

imp_ind <- cal_imp_n(3, 5)

# calculate imputed category parameter, based on uniform distribution
imp_category_sub <- function(para, imp_ind){
  K <- length(imp_ind$ind2) # total length
  res <- numeric(K+1) # to involve discrimination parameter
  res[1] <- para[1] # a para
  i <- 1
  for(k in 1:K){
    if(imp_ind$ind2[k]){
      i <- i+1
      res[k+1] <- para[i]
    } else if((i+1) <= length(para)) {
      b1 <- para[i]
      b2 <- para[i+1]
      distance <- (b2-b1)/(imp_ind$ind1[i-1] + 1)
      res[k+1] <- res[k+1] + b1
      res[(k+1):(k+imp_ind$ind1[i-1])] <- res[(k+1):(k+imp_ind$ind1[i-1])] + distance
    }
  }
  res
}

# calculate imputed category parameter, based on empirical distribution
imp_category_sub2 <- function(para, imp_ind, distribution){
  K <- length(imp_ind$ind2) # total length
  res <- numeric(K+1) # to involve discrimination parameter
  res[1] <- para[1] # a para
  # increasement of b parameter
  i <- 1
  increasement <- numeric(length(imp_ind$ind1))
  for(j in 1:length(imp_ind$ind1)){
    denominator <- sum(distribution[seq(i, length.out = imp_ind$ind1[j] + 1)])
    numerator <- abs(para[j+1] - para[j+2])
    increasement[j] <- numerator/denominator
  }
  # imputation
  i <- 1
  for(k in 1:K){
    if(imp_ind$ind2[k]){
      i <- i+1 # index for original b parameter
      res[k+1] <- para[i] # insert original b parameter to 'res' vector(res[1] is a parameter)
    } else if((i+1) <= length(para)) { # not to reach K + 1(not k + 1!) in original parameter vector
      b1 <- para[i]
      b2 <- para[i+1]
      inc <- increasement[i-1]
      # このあとは適切なヒストグラムの度数を求めて，それをincrasementに乗じて，b1に足せばOK
      distance <- (b2-b1)/(imp_ind$ind1[i-1] + 1)
      res[k+1] <- res[k+1] + b1
      res[(k+1):(k+imp_ind$ind1[i-1])] <- res[(k+1):(k+imp_ind$ind1[i-1])] + distance
    }
  }
  res
}

test_cate <- c(rep(1,100), rep(2,500), rep(3,1500), rep(4,500), rep(5,200), rep(6, 50), rep(7, 150))
test_table <- table(test_cate)
test_para <- c(1, -1, 0, 1)
imp_ind <- cal_imp_n(3, 6)

# test <- est_para_7[[1]][-6]
# test <-  test$Item001[-2,]
#
# imp_category_sub(test, imp_ind)

# imputation function
imp_category <- function(coef_object, K = 16){
  J <- length(coef_object)
  para <- coef_object[-J] %>% # remove estimated group parameter in list
    map(function(x)x["pars",]) # remove SE row
  # the number of category parameter of pre imputation(merged)
  k0 <- para %>% map(function(x){length(x)-1})
  # the target number of category parameter
  k1 <- para %>% map(function(x){K}) # K in not the number of categories but category difficulty parameters
  imp_ind <- map2(k0, k1, cal_imp_n)
  map2(para, imp_ind, imp_category_sub)
}

# imputation function that is based on subject frequency observed in merged category
imp_category <- function(coef_object, reponse_data, K = 16){
  #----NOW CODING
}

#----MAP theta estimation algorithm----

# Log likelihood
LLgrm_sub <- function(u, theta, para, D){
  a <- para[1]
  b <- para[-1]
  p1 <- log(pgrm(theta = theta, b = b, a = a, k = u, D = D))
}
LLgrm <- function(u, theta, para, D = 1.0, mu = 0, sigma = 1){
  ul <- as.list(u)
  lp <- purrr::map2_df(ul, para, LLgrm_sub, theta = theta, D = D)
  sum(lp) + log(dnorm(theta, mu, sigma))
}
LLgrm_apply <- function(u, para, D = 1.0, mu = 0, sigma = 1,  optimise_range = c(-8,8)){
  opt <- optimise(LLgrm, interval = optimise_range, u = u, para = para, maximum = T, mu = mu, sigma = sigma)
  opt$maximum
}

ext_par <- function(x){x["par",]} # extract only item parameter
# test <- est_para_17[[1]][-6] %>% map(pluck, ext_par) # map & pluck
modif_mirt_para <- function(para, n_of_items){
  tmp <- para[-(1+n_of_items)]
  tmp %>% purrr::map(purrr::pluck, ext_par)
}

# # optimize
# theta_sample[[1]] %>% head
# map1 <- apply(resp[[1]][,-1], 1, LLgrm_apply, para = test, D = 1.0, sigma = 3) # map: prior is N(0, 1)
#
# plot(theta_sample_list[[1]], map1, pch = 20, cex = 0.1)
# cor(theta_sample[[1]], map1)

# # validation
# # set.seed(old_seed[1])
# theta_sample <- rnorm(5000, 1, 1) # scale metric
# item_cate <- rep(5, 100)
# item_sample <- grm_para_gen(item_cate, a_dist = "beta", args_b_dist = list("mean" = 0.5, "sd" = 2),
#                             args_a_dist = list("shape1" = 2,"shape2" = 2))
# dat1 <- sim_grm_gen(theta_sample, item_sample, min_resp = 1) # a little bit slow
# dat1[,-1] %>% map(table) %>% map(length) %>% unlist %>% min
# mod <- mirt.model('B = 1-100
#                    START = (GROUP, MEAN_1, 1)')
# fit <- mirt(dat1[,-1], mod, itemtype = "graded", SE = T)
# para1 <- coef(fit, IRTpars = T, printSE = T)
# para2 <- modif_mirt_para(para1, 100)
# # map estimate
# map2 <- apply(dat1[,-1], 1, LLgrm_apply, para = para2, D = 1.0, sigma = 3)
# # visualization
# plot(theta_sample, map2, pch = 20, cex = 0.1)
# cor(theta_sample, map2)
#
# # fpd
# fpdgrm <- function(){
#
# }
# fpdgrm_sub <- function(u, theta, para, D){
#   a <- para[1]
#   b <- para[-1]
#   p1 <- pgrm(theta, b, a, u, D = D)
#   p0 <- pgrm(theta, b, a, u-1, D = D)
#
# }
# # spd?
#
# # fi?


# wrapper function
estgrmtheta <- function(dat, para, fc = 2, gc = 0, IDc = 1, D = 1.0, mu = 0, sigma = 1, optimise_range = c(-8,8) ){
  if(IDc == 0){
    id <- c(1:nrow(dat1))
  } else {
    id <- dat[IDc,]
  }
  xall <- dat[,fc:ncol(dat)]
  if(gc == 0){
    group <- rep(1, nrow(xall))
  } else {
    group <- dat[,gc]
  }
  cat("START MAP ESTIMATION!\n")
  map1 <- apply(xall, 1, LLgrm_apply, para = para, D = D, mu = mu, sigma = sigma, optimise_range = optimise_range)
  tibble::tibble(ID = id, GROUP = group, SCORE = rowSums(xall, na.rm = T), MAP = map1)
}
# estgrmtheta(resp[[1]][,-1], test, fc = 1, IDc = 0)
