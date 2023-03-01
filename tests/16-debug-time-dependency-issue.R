# Fix the likely bug of early-time loss of power

library(devtools)
load_all()

if(!exists("s01")) {
  o  <- load("data/testcubes_easystripes.rda") |> get()
  r1 <- o$strip
  j1 <- o$jump
  r2 <- o$strip_rev
  j2 <- o$jump_rev

  rd1 <- tj_stars_to_data(r1)
  jd1  <- as_tibble(j1) |> transmute(c_x = x, c_y = y, jump_true = values, has_jump = jump_true < 11)
  rd2 <- tj_stars_to_data(r2)
  jd2  <- as_tibble(j2) |> transmute(c_x = x, c_y = y, jump_true = values, has_jump = jump_true < 11)

  ## The starting point:
  prior_theta <- list(m = c(a=100, b=0),
                      S = diag(1e5*c(.001, .000001)))
  fitter1 <- function(d) tj_fit_m0.6_v1(dat = d,
                                    gamma = 0.1,
                                    prior_k = 0.5,
                                    niter = 500,
                                    prior_theta = prior_theta,
                                    prior_sigma2 = c(10, 500),
                                    verbose = TRUE)
  # parse relevant summaries:
  sumit <- function(m, j, nam) {
    ss <- tj_summarise_m0.3(m)
    ss |> left_join(j) |>
      mutate(name = nam, jump_pred = pred_jump_mode)
  }

  set.seed(1)
  t01 <- system.time(m01 <- fitter1(rd1))
  s01 <- sumit(m01, jd1, "0f")

  #t02 <- system.time(m02 <- fitter1(rd2))
  #s02 <- sumit(m02, jd2, "0r")
}

## Then the new code.
fitter2 <- function(d, ...) tj_fit_m0.6(dat = d,
                                   gamma = 0.1,
                                   prior_k = .5,
                                   niter = m01$niter,
                                   prior_theta = prior_theta,
                                   prior_sigma2 = c(10, 500),
                                   ...)

if(!exists("s11")) {
  set.seed(1)
  t11 <- system.time(m11 <- fitter2(rd1, verbose = TRUE))
  s11 <- sumit(m11, jd1, "1f")
#  t12 <- system.time(m12 <- fitter2(rd2))
#  s12 <- sumit(m12, jd2, "1r")
}

print(rbind(t01, t11))



################################################################################################
#
## Results plots and tables
##
##
# Check the time dependence
res <- bind_rows(s01,
  #              s02,
                 s11,
   #              s12
                 ) |>
  filter(c_x > min(c_x), c_x < max(c_x))
qtab <- res |>
  group_by(name, jump_true) |>
  summarise( TP = sum( jump_pred == jump_true) / n(),
             FP = sum(jump_pred) /n() )

ftab <- res |> group_by(name) |> summarise(TP = sum( jump_pred[jump_true <11] < 11)/sum(jump_true< 11),
                                           FP = sum( jump_pred[jump_true==11] < 11)/sum(jump_true==11))

print(ftab)

library(ggplot2)
library(patchwork)

qt <- qtab |>
  ggplot() +
  geom_point(aes(jump_true, TP)) + ylim(0:1) +
  facet_wrap(~name)

# maps
qm <- res |> ggplot() +
  geom_raster(aes(c_x, c_y, fill = factor(jump_pred)), show.legend = FALSE) +
  scale_fill_manual(values = c(pals::kelly(12)[-(1:2)], "black")) +
  facet_wrap(~name, ncol = 1) +
  coord_fixed(exp=F)

qmc <- with(m01,
            tibble(sigma2= hist_sigma2, iter = 1:niter)) |>
              ggplot() + geom_line(aes(iter, sigma2))



# example cells and their estimates:
cells <- c(2, 50, 80)  +  dim(r1)[["x"]]
p01 <- tj_predict_m0.3(m01, cells, s = s01)
p11 <- tj_predict_m0.3(m11, cells, s = s11)
qmb <- rd1 |> filter(cell %in% cells) |>
  ggplot() + geom_point(aes(time, value)) +
  geom_line(aes(time, value), data = p01, col = 2) +
  geom_line(aes(time, value), data = p11, col = 3) +
  facet_wrap(~cell)



print( qt/qm/qmc/qmb )



### mis
if(0) {
  plot(j2)
}
