MA 589 Project 3
================
Your Name
2024-03-22

# 1

## (a)

## (b)

``` r
library(stats)
n_samples <- 100
U <- runif(n_samples)
#inverse CDF quantile function
qasin <- function(U) {
  (sin(pi*U) + 1) /2
}
samples_MC <-qasin(U)
f_quantiles <- function (x, Q, ...) {
  n <- length(x)
  Q(seq(n) / (n + 1), ...)[order(x)]
}

# QQ plot against quantile function `Q`, with confidence bands and line
plot_qq <- function (x, Q, conf.level = .95) {
  qqplot(x, f_quantiles(x, Q), conf.level = .95)
  abline(0,1)
  #qqline(x, dist = Q, datax = TRUE)
}
plot_qq(samples_MC, qasin)
```

![](project3_files/figure-gfm/unnamed-chunk-1-1.png)<!-- --> \## (c)

``` r
s <- function(x){
  (1-2*x)^2
}
#z applied to the numbers sampled from the arcsine distribution
z <- s(samples_MC)
mean_MC <- mean(z)
var_MC <- var(z)

#applying the cos(pix)^2 function instead of using the change of variable
cos_px <- function(x) {
  cos(pi*x)^2
}
n <- cos_px(U)
mean_n <- mean(n)
var_n <- var(n)
```

When using the same 100 samples from the uniform distribution, we var_n
is the variance from the cos(pi\*x)^2 distribution var_MC is the
variance from the change of variable distribution var_n = var_MC,
however mean_n != mean_MC, although they are very close. \## (d)

``` r
#inverse CDF for wedge distribution
qwedge <- function (x) {
  #returns (1 + 1*sqrt(2*x - 1))/2 if x >= 0.5, else returns (1 - 1*sqrt(1-2*x))/2 )
  ifelse(x >= 0.5, (1 + 1*sqrt(2*x - 1))/2,
    (1 - 1*sqrt(1-2*x))/2 )
}
```

If S is bernoulli(0.5) we have 0.5 chance of dealing with (1+1sqrt(U))/2
0.5 chance of dealing with (1-1sqrt(U))/2 U ~ U(0,1) this is the same
distribution as the wedge distribution, which takes the sqrt of l = 1-2U
or l = 2U-1 but only over the ranges of 0-0.5 and 0.5-1 respectively, so
these l’s are ~ U(0,1), therefore they have the same distribution.

## (e)

``` r
samples_wedge <- qwedge(U)
```

    ## Warning in sqrt(2 * x - 1): NaNs produced

    ## Warning in sqrt(1 - 2 * x): NaNs produced

``` r
w <- function(x) {
  ifelse (x <= 0.5, cos(pi*x)^2 / (2-4*x), cos(pi*x)^2 /(4*x-2) )
}
w_intervals <- w(samples_wedge)
w_mean <- mean(w_intervals)
```

# 2

## (a)

If U is sampled from a uniform distribution 0 to 1, this is equivalent
to the given probability distribution because given a categorical
distribution, the cumulative sum of the probabilities will be a discrete
list that goes from 0 to 1, so for a given U, the minimum k such that U
\<= Sum (pi) from i to k, has a probability pk. If the pi’s are
0.1,0.3,0.2,0.4 then the cumulative probabilities are 0.1,0.4,0.6,1.
This maps to the uniform distribution so if we sample from a uniform it
has a p\[i\] - p\[i-1\] probability of being in any given category.

``` r
num_samples <- 100
u2 <- runif(num_samples)
log_u <- log(u2)
#log pi list
#log_pi <- [#list of log(pi)]
pi <- c(0.1,0.2,0.1,0.3,0.15,0.1,0.05)

lse2 <- local({
LOGEPS <- log(.Machine$double.eps / 2) # cache 
function (x, y) {
m <- pmax(x, y); d <- -abs(x - y)
ifelse(d < LOGEPS, m, m + log(1 + exp(d))) # m + log1pe(d) 
}})

#calculate the log of the cumulative sum of the e^log_pi's
cum_sum <- function(log_pi) {
  len <- length(log_pi)
  cum_sums <- numeric(len) #initialize an empty array
  cum_sums[1] <- log_pi[1]
  for (i in 2: length(log_pi)) {
    cum_sums[i] <- lse2 (cum_sums[i-1], (log_pi[i]))
  }
  cum_sums
}
#avoid underflows by working in log space !!!
hyp_geo_sampler <-function(U, cum_sums){
  log_U <- log(U)
  log_U <= cum_sums
  min(which(log_U <= cum_sums))
}
```

## (b)

``` r
m1 <- 200
m2 <- 775
s1 <-205
theta <- 3
x_j <- function(theta,m1,m2,s1,j) {
(lchoose(m1, j) + lchoose(m2, s1 - j) + theta * j) 
 }
a <- max(0,s1-m2)
b <- min(m1,s1)

result_list <- list()
for (i in a:b) {
  result_list[i] <- x_j(3,200,775,205, i)
}

lse2 <- local({
LOGEPS <- log(.Machine$double.eps / 2) # cache
function (x, y) {
m <- pmax(x, y); d <- -abs(x - y)
ifelse(d < LOGEPS, m, m + log(1 + exp(d)))# m + log1pe(d) 
}
})
#log sum of exponents
lse <- function (x) Reduce(lse2, x)

lprobfunc_3 <- lse(c(result_list)) #762.9187

#function to calculate log(f-theta(y)) for theta = 3
log_pi_3 <- list()
  for (i in a:b) {
    log_pi_3[i] <- lchoose(m1,i) + lchoose(m2,s1-i) + theta*i - lprobfunc_3
  }
log_pi_3 <- as.numeric(unlist(log_pi_3))
#sampler
u3<- runif(n_samples)

cum_sums_hyp_3 <- cum_sum(log_pi_3)
hyp_3_sampler <- lapply(u3,hyp_geo_sampler, cum_sums=cum_sums_hyp_3)
hyp_3_sampler <- as.numeric(hyp_3_sampler)
mean(hyp_3_sampler) #134.55
```

    ## [1] 134.56

## (c)

``` r
m1 <- 20
m2 <- 75
s1 <- 21
y <-10
theta <- 0

a <- max(0,s1-m2)
b <- min(m1,s1)

result_list_2 <- list()
for (i in a:b) {
  result_list_2[i] <- x_j(0,20,75,21,i)
}
#P0(theta=0)
lprobfunc_0 <- lse(c(result_list_2))
#log of the pi's for the hyp geo distribution
log_pi_0 <- list()
  for (i in a:b) {
    log_pi_0[i] <- lchoose(m1,i) + lchoose(m2,s1-i) + theta*i - lprobfunc_0
  }
#convert to a numeric vector from a list
log_pi_0 <- as.numeric(unlist(log_pi_0))

#sampler
u4<- runif(500)
#cumulative sums to sample from uniform to hyp geometric distribution
cum_sums_hyp_0 <- cum_sum(log_pi_0)

hyp_0_sampler <- lapply(u4,hyp_geo_sampler, cum_sums=cum_sums_hyp_0)
hyp_0_sampler <- as.numeric(hyp_0_sampler)

#P(Y>=y) aka p-value
p_value_y10 <-sum(hyp_0_sampler >= 10) /sum(hyp_0_sampler > 0) #0.002, got 0 with 100 trials
```

## (d)

``` r
variance_est <- var(hyp_0_sampler) #1000 samples
#2.711679
```

## (e)

``` r
#HG Sampler in One Function
hg_sampler <- function(m1,m2,s1,theta,num_samples) {
a <- max(0,s1-m2)
b <- min(m1,s1)
#calculate the P_0(theta)
#function to calculate the exponent of each element
x_j <- function(theta,m1,m2,s1,j) {
(lchoose(m1, j) + lchoose(m2, s1 - j) + theta * j) 
} 
#vector of exponents
xj_vec <- list()
for (i in a:b) {
  xj_vec[i] <- x_j(theta,m1,m2,s1,i)
}
#log P_0 = summation of exponents in log space using lse
lprobfunc <- lse(c(xj_vec))

#log of pi's for hyp-geo distribution
log_pi <- list()
  for (i in a:b) {
    log_pi[i] <- lchoose(m1,i) + lchoose(m2,s1-i) + theta*i - lprobfunc
  }
#change var type to a vector of numbers
log_pi <- as.numeric(unlist(log_pi))
#create a vector u that samples from random uniform
u <- runif(num_samples)
#create cumulative sums of probabilities to sample from
cum_sums <- cum_sum(log_pi)
#sample by getting the minimum cum_sum s.t. U <= cum_sum IN LOG SPACE
#don't need cum_sums as an input bc it is defined in the HG function directly above so this function inherits it. 
hyp_geo_sampler <-function(U){
  log_U <- log(U)
  log_U <= cum_sums
  min(which(log_U <= cum_sums))
}
#apply 
hyp_samples <- as.numeric(lapply(u,hyp_geo_sampler))

hyp_samples
}
#samples for m1=20,m2=75,s1=21,theta=1.5, n_samples = 100
hg_15_100 <- hg_sampler(20,75,21,1.5,100)
#n_samples = 1000
hg_15_1000 <- hg_sampler(20,75,21,1.5,1000)
var(hg_15_100) #3.646465
```

    ## [1] 2.989495

``` r
var(hg_15_1000) #3.402839
```

    ## [1] 3.383848

``` r
#the variance of the samples when theta =1.5 is significantly larger, results vary every time this is run from var = 3.15 to 3.5

#p_value for Y>=y, y=10
p_value_y10 <-sum(hg_15_1000 >= 10) /1000 #p-val = 0.422
```
