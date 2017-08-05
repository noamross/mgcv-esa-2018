mgcv miscellanea 
=================
author:Noam Ross
date: 2017-08-04
css: custom.css
transition: none


We've done a lot today. But there's more!
=========================================

Some other topics we can touch on:

- More nuanced prediction: variation in smooth shape and derived values
- A beastiary of smooth types
- Random effects and heirarchical GAMs (GAMMs!): Three different ways
- Different slopes for different folks: slope-random effect interactions
- New approaches to full Bayesian GAMs with JAGS or Stan

Exercises
=========
In the project folder you will find a series of `.Rmd` files with additional
excercises:

- `example-spatial-mexdolphins.Rmd`: An extended version of this morning's
dolphin excercise. Includes using _quantile residuals_, another checking tool,
which you'll find in the `code_snippets/` folder.

- `example-spatio-temporal-data.Rmd`: Using _tensor_ smooths to model and 
separate interactions between smooths that may operate at different scales.
Also random effects.

- `example-linear-functional.Rmd`: For when your $y$ outcome variable is
dependent on a nonlinear or weighted average of multiple $x$ variables.


Exercises (2)
=============

- `example-nonlinear-timeseries.Rmd`:  Time series analysis with decomposition
into components. Cyclic (seasonal) smooths, temporal autocorrelation using
`gamm()`

- `example-bivariate-timeseries-and-ti.Rmd`: Extending time series analysis
to interaction between the components at different time scales

-  `exampled-forest-health.Rmd`: Modeling ordered categorical outcomes, and
discrete spatial units with Markov random fields


