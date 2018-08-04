#' # Preamble
#' 
#' This exercise is based on the [Appendix of Miller et al 2013](http://distancesampling.org/R/vignettes/mexico-analysis.html). In this example we're ignoring all kinds of important things like detectability and availability. This should not be treated as a serious analysis of these data! For a more complete treatment of detection-corrected abundance estimation via distance sampling and generalized additive models, see Miller et al. (2013).
#' 
#' From that appendix:
#' 
#' *The analysis is based on a dataset of observations of pantropical dolphins in the Gulf of Mexico (shipped with Distance 6.0 and later). For convenience the data are bundled in an `R`-friendly format, although all of the code necessary for creating the data from the Distance project files is available [on github](http://github.com/dill/mexico-data). The OBIS-SEAMAP page for the data may be found at the [SEFSC GoMex Oceanic 1996](http://seamap.env.duke.edu/dataset/25) survey page.*
#' 
#' 
#' 
#' # Data format
#' 
#' The data are provided in the `data/mexdolphins` folder as the file `mexdolphins.RData`. Loading this we can see what is provided:
#' 
#' 
#' If you haven't already installed these packages, do so:
install.packages(c(
  "mgcv",
  "plyr",
  "ggplot2",
  "viridis",
  "dplyr",
  "tidyr",
  "mvtnorm",
  "statmod"
))

## ----loaddata------------------------------------------------------------
load("data/mexdolphins/mexdolphins.RData")
ls()

#' 
#' - `mexdolphins` the `data.frame` containing the observations and covariates, used to fit the model.
#' - `pred_latlong` an `sp` object that has the shapefile for the prediction grid, used for fancy graphs
#' - `preddata` prediction grid without any fancy spatial stuff
#' 
#' Looking further into the `mexdolphins` frame we see:
#' 
## ----frameinspect--------------------------------------------------------
str(mexdolphins)

#' 
#' A brief explanation of each entry:
#' 
#' - `Sample.Label` identifier for the effort "segment" (approximately square sampling area)
#' - `Transect.Label` identifier for the transect that this segment belongs to
#' - `longitude`, `latitude` location in lat/long of this segment
#' - `x`, `y` location in projected coordinates (projected using the [North American Lambert Conformal Conic projection](https://en.wikipedia.org/wiki/Lambert_conformal_conic_projection))
#' - `Effort` the length of the current segment
#' - `depth` the bathymetry at the segment's position
#' - `count` number of dolphins observed in this segment
#' - `segment.area` the area of the segment (`Effort` multiplied by the width of the segment
#' - `off.set` the logarithm of the `segment.area` multiplied by a correction for detectability (see link to appendix above for more information on this)
#' 
#'
#' Let's plot some of this data to get an idea what we're dealing with: 
## ------------------------------------------------------------------------

library(ggplot2)
library(viridis)
ggplot(mexdolphins, aes(x=depth, y=count)) + geom_point()
ggplot(mexdolphins, aes(x=x, y=y, size=count)) + geom_point() + coord_equal()
ggplot() +
  stat_summary_2d(data= preddata, mapping=aes(x=x, y=y, z=depth), binwidth = c(20000, 20000)) +
  scale_fill_viridis() +
  geom_point(data=mexdolphins, mapping=aes(x=x,y=y, size=count), col="orange") +
  coord_equal()


#' 
#' 
#' # Modelling
#' 
#' Our objective here is to build of where and how many dolphins there are In some sense this is a kind of species distribution model. Our possible covariates to model abundance are location and depth.
#' 
#' Here is an example of a simple GAM we could build with *mgcv* 
## ----simplemodel---------------------------------------------------------
library(mgcv)
d_depth <- gam(count ~ s(depth) + offset(off.set),
                      data = mexdolphins,
                      family = poisson(),
                      method = "REML")


summary(d_depth)

plot(d_depth)

#'  We could build a more complex gam where we model population as a function
#'  of space.  Here we use a two dimensional term of `x, y`, and use 
#'  `vis.gam`
## ------------------------------------------------------------------------
d_xy <- gam(count ~ s(x, y) + offset(off.set),
                 data = mexdolphins,
                 family=poisson(),
                 method="REML")

summary(d_xy)

par(mfrow=c(2,2))
vis.gam(d_xy, view=c("x","y"), phi=45, theta=20, asp=1)
vis.gam(d_xy, view=c("x","y"), phi=45, theta=60, asp=1)
vis.gam(d_xy, view=c("x","y"), phi=45, theta=160, asp=1)
dev.off()


#' Since this is an _additive_ model, we can of course include multiple
#' terms.  
## ------------------------------------------------------------------------
d_dxy <- gam(count ~ s(depth) + s(x, y) + offset(off.set),
                 data = mexdolphins,
                 family=poisson,
             method="REML")

summary(d_dxy)
plot(d_dxy, scale=0, pages=1, scheme=2)


#' We can generate predictions from our model at new data points.  Here
#' we generate 'response' predictors - one the scale of the data - and
#' plot the outputs
## ------------------------------------------------------------------------
d_preds <- predict(d_dxy, newdata=preddata, type="response")
str(d_preds)
preddata$prediction <- d_preds

ggplot() +
  stat_summary_2d(data= preddata, mapping=aes(x=x, y=y, z=prediction), binwidth = c(20000, 20000)) +
  coord_equal() +
  scale_fill_viridis()


#' So what just happened? Let's go back to the slides!
#'---

#' Inspecting the GAM object
## ----coefs------------------------------------------------------------
coef(d_dxy)
dim(model.matrix(d_dxy))


#' ---
#' 
#' # Exercise!
#' 
#' -  We set $k$ in `s()` terms with `s(variable, k=n)`
#' -  Re-fit models with small, medium and large $k$ values.
#' -  Look at coefficients, model.matrix, summaries and plots
#' -  How do deviance explained, EDF values, smooth shapes change?
#' -  What is default $k$?
#' 
#' ---


#' Prediction
#' 
#' Let's look at the variance of our parameters, and how it changes
#' when we account for uncertainty in the smooth
vcov(d_dxy)
vcov(d_dxy, unconditional = TRUE)
plot(diag((vcov(d_dxy, unconditional = TRUE) - vcov(d_dxy))/vcov(d_dxy)))



#' We make predictions on the response scale with type="response"
preddata2 <- within(preddata, {x = mean(x); y = mean(y)})
pred <- predict(d_dxy, newdata=preddata2, se.fit=TRUE, type="response")
str(pred)
preddata2$pmean <- pred$fit
ggplot(preddata2, aes(x=depth, y = pmean)) +
  geom_line()

#' However, to look at errors on the response scale, we need to use the 'link' type
#' and convert ourselves
pred2 <- predict(d_dxy, preddata2, se.fit=TRUE, type="link", unconditional=TRUE)
preddata2$pmean <- exp(pred2$fit)
preddata2$lo <- exp(pred2$fit - 2* pred2$se.fit)
preddata2$hi <- exp(pred2$fit + 2* pred2$se.fit)
ggplot(preddata2, aes(x=depth, y = pmean)) +
  geom_line() +
  geom_line(mapping=aes(y=lo), lty=2) +
  geom_line(mapping=aes(y=hi), lty=2)

#' An exercise
#' Make predictions of all the points in preddata, not just at the mean
#' location
#' 

#' Varying smooth shapes
library(mvtnorm)
lp <- predict(d_dxy, preddata2, se.fit=TRUE, type="lpmatrix", unconditional=TRUE)
str(lp)
coef_samps <- rmvnorm(20, mean = coef(d_dxy), sigma = vcov(d_dxy, unconditional = TRUE))
str(coef_samps)
pred_samps <- lp %*% t(coef_samps)
str(pred_samps)
rownames(pred_samps) <- preddata2$depth
curves <- reshape2::melt(pred_samps)
ggplot(curves, aes(x = Var1, y = value, group=Var2)) + geom_line(alpha=0.1)
ggplot(curves, aes(x = Var1, y = exp(value), group=Var2)) + geom_line(alpha=0.1)
