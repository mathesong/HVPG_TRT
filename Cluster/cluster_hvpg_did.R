args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)

simset <- args[1]

set.seed(simset)
nsims = 50

print(simset)

trt_all <- readRDS("trt_all.rds")


HVPG_difindif_sim <- function(n, delta1, delta2, cv_delta, 
                              wscv=trt_all$wscv, 
                              mean=trt_all$mean, 
                              cv=trt_all$cv, 
                              icc=trt_all$icc,
                              decomp ) {
  
  wscv <- wscv[decomp]
  mean <- mean[decomp]
  cv   <- cv[decomp]
  icc   <- icc[decomp]
  
  var <- (cv*mean)^2
  sd_true <- sqrt(var * icc)
  
  pre_true1 <- rnorm(n, mean, sd_true)
  pre_true2 <- rnorm(n, mean, sd_true)
  
  pre_meas1 <- pre_true1 + rnorm(n, 0, mean*wscv)
  pre_meas2 <- pre_true2 + rnorm(n, 0, mean*wscv)
  
  post_true1 <- pre_true1 - rnorm(n, delta1, cv_delta*delta1)
  post_true2 <- pre_true2 - rnorm(n, delta2, cv_delta*delta2)
  
  post_meas1 <- post_true1 + rnorm(n, 0, mean*wscv)
  post_meas2 <- post_true2 + rnorm(n, 0, mean*wscv)
  
  measured <- tibble::tibble(
    ID = rep(1:n, times=2),
    Pre = c(pre_meas1, pre_meas2),
    Post = c(post_meas1, post_meas2),
    Diff = Post - Pre,
    Group = rep(c("A", "B"), each=n)
  )
  
  d <- effsize::cohen.d(measured$Diff, measured$Group, 
                        paired=FALSE)$estimate
  
  mod <- lm(Post ~ Pre + Group, data=measured)
  
  testout <- broom::tidy(mod) %>% 
    filter(term=="GroupB") %>% 
    select(-term) %>% 
    mutate(p.value = pt(statistic, mod$df, lower.tail = FALSE))
  # Note: using a one-sided p value
  
  out <- mutate(testout, d=d)
  
  return(out)
  
}


difindifsimpars <- tidyr::crossing(
  n = c( seq(6, 200, by = 2)),
  delta1 = seq(1,3, by = 0.5),
  delta2 = seq(0,2, by = 0.5),
  cv_delta = c(0, 0.5),
) %>% 
  mutate(deltadif = delta1 - delta2) %>% 
  filter(delta1 > delta2)


head(difindifsimpars)

difindifsims <- difindifsimpars %>% 
  mutate(sim = 1:nrow(difindifsimpars)) %>% 
  group_by(sim) %>% 
  nest(params = c(n, delta1, delta2, cv_delta)) %>% 
  mutate(output = map(params, .progress = TRUE,
                      ~bind_rows(purrr::rerun(nsims, 
                                               HVPG_difindif_sim(.x$n, .x$delta1, 
                                                                 .x$delta2, .x$cv_delta, 
                                                                 decomp = 1)))))

head(difindifsims)
  
saveRDS(difindifsims, 
        paste0("difindifsims_decomp_", simset, ".rds"))
