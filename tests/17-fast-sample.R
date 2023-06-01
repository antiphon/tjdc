# check the speed of our sample

devtools::load_all()


library(microbenchmark)


# rcat
if(0) {

K <- 11
prob <- exp((1:K))
prob <- prob / sum(prob)

mb <- microbenchmark(
  regu = sample(K, 1, prob = prob),
  reguv = sample(1:K, 1, prob = prob),
  regu.int = sample.int(K, 1, prob = prob, useHash = FALSE),
  fast = tj_fast_sample(prob),
  fast_c = sample_fast_c(prob))

# check results
y <- sapply(1:100000, \(i) sample_fast_c(prob))


par(mfrow=c(2,1))
boxplot(mb)
plot(table(y)/length(y))
points(1:11, prob)
mb
}

#rtnorm
if(1) {
  # current mine
  mu <- 0
  sigma <- 10
  a <- 0
  b <- 100
  rtnorm <- function() qnorm( pnorm( (a-mu)/sigma ) + runif(1) * (pnorm((b-mu)/sigma)-pnorm((a-mu)/sigma)) ) *sigma + mu
  mb <- microbenchmark(truncnorm =   truncnorm::rtruncnorm(n = 1, a = a, b = b, mean = mu, sd = sigma),
                       myrtnorm = rtnorm(),
                       rnorm = rnorm(1, mu, sigma))


  print(mb)
  boxplot(mb)

}



