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
#' 
#' 
#' 
#' 
#' That is we fit the counts as a function of depth, using the offset to take into account effort. We use a quasi-Poisson response (i.e., modelling just the mean-variance relationship such that the variance is proportional to the mean) and use REML for smoothness selection.
#' 
#' We can check the assumptions of this model by using `gam.check`:
#' 
## ----simplemodel-check---------------------------------------------------
par(mfrow=c(2,2))
gam.check(dolphins_depth_xy)
par(mfrow=c(1,1))

#' 
#' As is usual for count data, these plots are a bit tricky to interpret. For example the residuals vs linear predictor plot in the top right has that nasty line through it that makes looking for pattern tricky. We can see easily that the line equates to the zero count observations:
#' 
## ----zeroresids, fig.width=7, fig.height=7-------------------------------
# code from the insides of mgcv::gam.check
resid <- residuals(dolphins_depth, type="deviance")
linpred <- napredict(dolphins_depth$na.action, dolphins_depth_xy$linear.predictors)
plot(linpred, resid, main = "Resids vs. linear pred.",
     xlab = "linear predictor", ylab = "residuals")

# now add red dots corresponding to the zero counts
points(linpred[mexdolphins$count==0],resid[mexdolphins$count==0],
       pch=19, col="red", cex=0.5)

#' 
#' We can use randomised quantile residuals instead of deviance residuals to get around this in some cases (though not quasi-Poisson, as we don't have a proper likelihood!).
#' 
#' Ignoring the plots for now (as we'll address them in the next section), let's look at the text output. It seems that the `k` value we set (or rather the default of 10) seems to have been adequate.
#' 
#' We could increase the value of `k` by replacing the `s(...)` with, for example, `s(depth, k=25)` (for a possibly very wiggly function) or `s(depth, k=3)` (for a much less wiggly function). Making `k` big will create a bigger design matrix and penalty matrix.
#' 
#' 
#' **Exercise**
#' 
#' Look at the differences in the size of the design and penalty matrices by using `dim(odel.matrix(...))` and `dim(model$smooth[[1]]$S[[1]])`, replacing `...` and `model` appropriately for models with `k=3` and `k=30`.
#' 
## ----simplemodel-bigsmall------------------------------------------------


#' 
#' (Don't worry about the many square brackets etc to get the penalty matrix!)
#' 
#' ### Plotting
#' 
#' We can plot the smooth we fitted using `plot`.
#' 
#' **Exercise**
#' 
#' Compare the first model we fitted with the two using different `k` values above. Use `par(mfrow=c(1,3))` to put them all in one graphic. Look at `?plot.gam` and plot the confidence intervals as a filled "confidence band". Title the plots appropriately so you can check which is which.
#' 
## ----plotk---------------------------------------------------------------
par(mfrow=c(1,3))


#' 
#' ## Count distributions
#' 
#' In general quasi-Poisson doesn't seem to do too great a job at modelling data with many zeros. Luckily we have a few tricks up our sleeves...
#' 
#' 
#' ### Tweedie
#' 
#' Adding a smooth of `x` and `y` to our model with `s(x,y)`, we can then switch the `family=` argument to use `tw()` for a Tweedie distribution.
#' 
## ----tw------------------------------------------------------------------
dolphins_xy_tw <- gam(count ~ s(x,y) + s(depth) + offset(off.set),
                         data = mexdolphins,
                         family = tw(),
                         method = "REML")

#' 
#' More information on Tweedie distributions can be found in Foster & Bravington (2012) and Shono (2008).
#' 
#' 
#' 
#' ### Negative binomial
#' 
#' **Exercise**
#' 
#' Now do the same using the negative binomial distribution (`nb()`).
#' 
## ----nb------------------------------------------------------------------


#' 
#' Looking at the quantile-quantile plots only in the `gam.check` output for these two models, which do you prefer? Why?
#' 
#' *Looks like Tweedie is better here as the points are closer to the x=y line in the Q-Q plot. Also the histogram of residuals looks more (though not very) normal.*
#' 
#' Look at the results of calling `summary` on both models and note that there are differences in the resulting models, due to the differing mean-variance relationships.
#' 
#' ## Smoothers
#' 
#' Now let's move onto using different bases for the smoothers. We have a couple of different options here.
#' 
#' 
#' ### Thin plate splines with shrinkage
#' 
#' By default we use the `"tp"` basis. This is just plain thin plate regression splines (as defined in Wood, 2003). We can also use the `"ts"` basis, which is the same but with extra shrinkage on the usually unpenalised parts of model. In the univariate case this is the linear slope term of the smooth.
#' 
#' **Exercise**
#' 
#' Compare the results from one of the models above with a version using the thin plate with shrinkage using the `bs="ts"` argument to `s()` for both terms.
#' 
## ----tw-ts---------------------------------------------------------------

#' 
#' What are the differences (use `summary`)?
#' 
#' What are the visual differences (use `plot`)?
#' 
#' 
#' ### Soap film smoother
#' 
#' We can use a soap film smoother (Wood, 2008) to take into account a complex boundary, such as a coastline or islands.
#' 
#' Here I've built a simple coastline of the US states bordering the Gulf of Mexico (see the `soap_pred.R` file for how this was constructed). We can load up this boundary and the prediction grid from the following `RData` file:
#' 
## ----soapy---------------------------------------------------------------
load("data/mexdolphins/soapy.RData")

#' 
#' Now we need to build knots for the soap film, for this we simply create a grid, then find the grid points inside the boundary. We don't need too many of them.
#' 
## ----soapknots-----------------------------------------------------------
soap_knots <- expand.grid(x=seq(min(xy_bnd$x), max(xy_bnd$x), length.out=10),
                          y=seq(min(xy_bnd$y), max(xy_bnd$y), length.out=10))
x <- soap_knots$x; y <- soap_knots$y
ind <- inSide(xy_bnd, x, y)
rm(x,y)
soap_knots <- soap_knots[ind, ]
## inSide doesn't work perfectly, so if you get an error like:
## Error in crunch.knots(ret$G, knots, x0, y0, dx, dy) :
##  knot 54 is on or outside boundary
## just remove that knot as follows:
soap_knots <- soap_knots[-8, ]
soap_knots <- soap_knots[-54, ]

#' 
#' We can now fit our model. Note that we specify a basis via `bs=` and the boundary via `xt` (for e`xt`ra information) in the `s()` term. We also include the knots as a `knots=` argument to `gam`.
#' 
## ----soapmodel-----------------------------------------------------------
dolphins_soap <- gam(count ~ s(x,y, bs="so", xt=list(bnd=list(xy_bnd))) +
                     offset(off.set),
                     data = mexdolphins,
                     family = tw(),
                     knots = soap_knots,
                     method = "REML")


#' 
#' **Exercise**
#' 
#' Look at the `summary` output for this model and compare it to the other models.
#' 
## ----soap-summary--------------------------------------------------------

#' 
#' The plotting function for soap film smooths looks much nicer by default than for other 2D smooths -- try it out.
#' 
#' 
## ----soap-plot-----------------------------------------------------------

#' 
#' # Predictions
#' 
#' As we saw in the intro to GAMs slides, `predict` is your friend when it comes to making predictions for the GAM.
#' 
#' We can do this very simply, calling predict as one would with a `glm`. For example:
#' 
## ----pred-qp-------------------------------------------------------------
pred_qp <- predict(dolphins_depth, preddata, type="response")

#' 
#' Now, this just gives a long vector of numbers for the predicted number of animals per cell. We can find the total abundance using `sum`.
#' 
#' **Exercise**
#' 
#' How many dolphins are there in the total area? What is the maximum in a given cell? What is the minimum?
#' 
#' ## Plotting predictions
#' 
#' *Note that this section requires quite a few additional packages to run the examples, so may not run the first time. You can use* `install.packages` *to grab the packages you need.*
#' 
#' Plotting predictions in projected coordinate systems is tricky. I'll show to methods here but not go into too much detail, as that's not the aim of this workshop.
#' 
#' For the non-soap models, we'll use the below helper function to put the predictions into a bunch of squares and then return an appropriate `ggplot2` object for us to plot:
#' 
## ----gridplotfn----------------------------------------------------------
library(plyr)
# fill must be in the same order as the polygon data
grid_plot_obj <- function(fill, name, sp){

  # what was the data supplied?
  names(fill) <- NULL
  row.names(fill) <- NULL
  data <- data.frame(fill)
  names(data) <- name

  spdf <- SpatialPolygonsDataFrame(sp, data)
  spdf@data$id <- rownames(spdf@data)
  spdf.points <- fortify(spdf, region="id")
  spdf.df <- join(spdf.points, spdf@data, by="id")

  # seems to store the x/y even when projected as labelled as
  # "long" and "lat"
  spdf.df$x <- spdf.df$long
  spdf.df$y <- spdf.df$lat

  geom_polygon(aes_string(x="x",y="y",fill=name, group="group"), data=spdf.df)
}

#' 
#' Then we can plot the predicted abundance:
#' 
## ----plotpred------------------------------------------------------------
library(ggplot2)
library(viridis)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)

# projection string
lcc_proj4 <- CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
# need sp to transform the data
pred.polys <- spTransform(pred_latlong, CRSobj=lcc_proj4)

p <- ggplot() +
      # abundance
      grid_plot_obj(pred_qp, "N", pred.polys) +
      # survey lines
      geom_line(aes(x, y, group=Transect.Label), data=mexdolphins) +
      # observations
      geom_point(aes(x, y, size=count),
                 data=mexdolphins[mexdolphins$count>0,],
                 colour="red", alpha=I(0.7)) +
      # make the coordinate system fixed
      coord_fixed(ratio=1, ylim = range(mexdolphins$y),
                  xlim = range(mexdolphins$x)) +
      # use the viridis colourscheme, which is more colourblind-friendly
      scale_fill_viridis() +
      # labels
      labs(fill="Predicted\ndensity",x="x",y="y",size="Count") +
      # keep things simple
      theme_minimal()
print(p)

#' 
#' If you have the `maps` package installed you can also try the following plot the coastline too:
#' 
## ----maptoo--------------------------------------------------------------
library(maps)
map_dat <- map_data("usa")
map_sp <- SpatialPoints(map_dat[,c("long","lat")])

# give the sp object a projection
proj4string(map_sp) <-CRS("+proj=longlat +datum=WGS84")
# re-project
map_sp.t <- spTransform(map_sp, CRSobj=lcc_proj4)
map_dat$x <- map_sp.t$long
map_dat$y <- map_sp.t$lat

p <- p + geom_polygon(aes(x=x, y=y, group = group), fill = "#1A9850", data=map_dat)

print(p)

#' 
#' ## Comparing predictions from soap and thin plate splines
#' 
#' The `soap_preddata` frame has a much larger prediction grid that covers a much wider area.
#' 
#' Predicting using the soap film smoother is exactly the same as for the other smoothers:
#' 
## ----soap-pred-----------------------------------------------------------
soap_pred_N <- predict(dolphins_soap, soap_preddata, type="response")
soap_pred_N <- cbind.data.frame(soap_preddata, N=soap_pred_N, model="soap")

#' 
#' Use the above as a template (along with `ggplot2::facet_wrap()`) to build a comparison plot with the predictions of the soap film and another model.
#' 
## ----xy-soap-pred--------------------------------------------------------
dolphins_xy_q <- gam(count ~ s(x,y, bs="ts") + offset(off.set),
                     data = mexdolphins,
                     family = tw(),
                     method = "REML")
xy_pred_N <- predict(dolphins_xy_q, soap_preddata, type="response")
xy_pred_N <- cbind.data.frame(soap_preddata, N=xy_pred_N, model="xy")
all_pred_N <- rbind(soap_pred_N, xy_pred_N)

p <- ggplot() +
      # abundance
      geom_tile(aes(x=x, y=y, width=sqrt(area), height=sqrt(area), fill=N), data=all_pred_N) +
      # facet it!
      facet_wrap(~model) +
      # make the coordinate system fixed
      coord_fixed(ratio=1, ylim = range(mexdolphins$y),
                  xlim = range(mexdolphins$x)) +
      # use the viridis colourscheme, which is more colourblind-friendly
      scale_fill_viridis() +
      # labels
      labs(fill="Predicted\ndensity", x="x", y="y") +
      # keep things simple
      theme_minimal() +
      geom_polygon(aes(x=x, y=y, group = group), fill = "#1A9850", data=map_dat)
print(p)

#' 
#' ## Plotting uncertainty
#' 
#' **Exercise**
#' 
#' We can use the `se.fit=TRUE` argument to get the per-cell standard errors, then divide these through by the abundances to get a coefficient of variation (CV) per cell. We can then plot that using the same technique as above. Try this out below.
#' 
#' Note that the resulting object from setting `se.fit=TRUE` is a list with two elements.
#' 
## ----CVmap---------------------------------------------------------------


#' 
#' Compare the uncertainties of the different models you've fitted so far.
#' 
#' ## `lpmatrix` magic
#' 
#' Now for some `lpmatrix` magic. We said before that the `lpmatrix` maps the parameters onto the predictions. Let's show that's true:
#' 
## ----lppred--------------------------------------------------------------
# make the Lp matrix
lp <- predict(dolphins_depth, preddata, type="lpmatrix")
# get the linear predictor
lin_pred <- lp %*% coef(dolphins_depth)
# apply the link and multiply by the offset as lpmatrix ignores this
pred <- preddata$area * exp(lin_pred)
# all the same?
all.equal(pred[,1], as.numeric(pred_qp), check.attributes=FALSE)

#' 
#' What else can we do? We can also grab the uncertainty for the sum of the predictions:
#' 
## ----lpNvar--------------------------------------------------------------
# extract the variance-covariance matrix
vc <- vcov(dolphins_depth)

# reproject the var-covar matrix to be on the linear predictor scale
lin_pred_var <- tcrossprod(lp %*% vc, lp)

# pre and post multiply by the derivatives of the link, evaluated at the
# predictions -- since the link is exp, we just use the predictions
pred_var <- matrix(pred,nrow=1) %*% lin_pred_var %*% matrix(pred,ncol=1)

# we can then calculate a total CV
sqrt(pred_var)/sum(pred)

#' 
#' As you can see, the `lpmatrix` can be very useful!
#' 
#' # Extra credit
#' 
#' - Experiment with `vis.gam` (watch out for the `view=` option) and plot the 2D smooths you've fitted. Check out the `too.far=` agument.
#' - Use the randomized quantile residuals plot to inspect the models you fitted above. What are the differences you see between them and the deviance residuals you see above? Which model would you choose?
#' - Redo the side-by-side soap film and another smoother plot, but showing the coefficient of variation. What difference does the soap film make?
#' 
#' # References
#' 
#' - Foster, S. D., & Bravington, M. V. (2012). A Poisson–Gamma model for analysis of ecological non-negative continuous data. Environmental and Ecological Statistics, 20(4), 533–552. http://doi.org/10.1007/s10651-012-0233-0
#' - Miller, D. L., Burt, M. L., Rexstad, E. A., & Thomas, L. (2013). Spatial models for distance sampling data: recent developments and future directions. Methods in Ecology and Evolution, 4(11), 1001–1010. http://doi.org/10.1111/2041-210X.12105
#' - Shono, H. (2008). Application of the Tweedie distribution to zero-catch data in CPUE analysis. Fisheries Research, 93(1-2), 154–162. http://doi.org/10.1016/j.fishres.2008.03.006
#' - Wood, S. N. (2003). Thin plate regression splines. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 65(1), 95–114.
#' - Wood, S. N., Bravington, M. V., & Hedley, S. L. (2008). Soap film smoothing. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(5), 931–955. http://doi.org/10.1111/j.1467-9868.2008.00665.x
#' 
