# CJSboost: machine-learning for capture-recapture

This repository contains a tutorial and code (R,C++) to accompany the pre-print publication ["EM and component-wise boosting for Hidden Markov Models: a machine-learning approach to capture-recapture", Robert W. Rankin (2016), bioRxiv pre-print, doi:10.1101/052266](http://www.biorxiv.org/content/early/2016/05/09/052266).  Interested users should study the preprint article [here](http://www.biorxiv.org/content/early/2016/05/09/052266). If you wish to run the tutorial, see below for how to download the files and a summary of the kinds of analyses you can run. To dive straight into the tutorial, see the file `R_TUTORIAL_dipper.R`.

![foo_you](https://github.com/faraway1nspace/HMMboost/blob/master/img/ANIMATION_stabselect.gif)

_Demonstration of selection probabilities for 30 high-dimensional simulations. As the boosting iteration (<i>m</i>) get large, regularization gets weaker, and all covariates have a higher selection probability S (estimated from a bootstrap). Lines in <b>red</b> are truly influential covariates. Lines in <b>gray</b> are non-influential covariates. Lines in <b>pink</b> are not-influential for THETA, but are influential in the other parameter (not THETA). Lines in <b>blue</b> represent the time-as-a-categorical-variable base-learner representing p(t) or phi(t), which in these simulations was non-influential._
## Getting Started
1. Download the CJSboost Github source code and dummy data
2. Install dependencies (see below)
3. Thoroughly read the [CJSboost manuscript](http://www.biorxiv.org/content/early/2016/05/09/052266), as well as documentation about mboost's "base-learners" and the R formula-interface. CJSboost draws heavily on `mboost`. See ["Model-based Boosting in R: A Hands-on Tutorial Using the R Package mboost"](https://epub.ub.uni-muenchen.de/12754/).
4. Make sure you can compile the source code `source("/path/to/HMMboost/R_CJSboost_SOURCE.R")`
5. Step through the tutorial file `R_TUTORIAL_dipper.R`

### Dependencies
* Working installation of [R](www.r-project.org) (available GNU/Linux, Mac, Windows)
* C++ libraries [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)
* install the R packages [mboost](https://cran.r-project.org/web/packages/mboost/index.html) (and [party](https://cran.r-project.org/web/packages/party/index.html), if you want to use boosted regression trees)
* (_optional_) working installation of [Program MARK](http://www.phidot.org/software/mark/) and [RMark](https://cran.r-project.org/web/packages/RMark/index.html) (only if you want to compare CJSboost to Mark's estimates). Linux and Mac users can find binaries [here](http://www.phidot.org/software/mark/rmark/linux/).

### Files
* `R_TUTORIAL_dipper.R` main tutorial file for running CJSboost analyses (_start here!_)
* `R_CJSboost_SOURCE.R` code for CJSboost functions;
* (optional) `R_simcjsdataset_SOURCE.R` source file for simulating fake datasets (used in tutorial);
* (optional) `data_dipper_lebreton.txt` dipper data, used in the tutorial (from [Lebreton et al 1992](http://dx.doi.org/10.2307/2937171/));

# Background
## What is CJSboost?

CJSboost uses a statistical learning technique called "boosting" for inference under a Cormack-Jolly-Seber (CJS) capture-recapture model. It iteratively builds an _ensemble_ of weak _base-learners_ to yield a strong prediction function. In univariate settings, it a powerful prediction technique (see [Schmid et al 2009](http://link.springer.com/article/10.1007/s11222-009-9162-7) for a good overview; [free pdf here](https://epub.ub.uni-muenchen.de/7788/1/TR.pdf)). The iterative weighting of base-learners can be thought of as a type of "multi-model inference", not-unlike model-averaging by AICc weights. Boosting and AIC model-averaging have some theorectical connections, both being motivated by predictive preformance and minimizing a loss criteria. However, unlike AIC model-averaging, boosting can incoporate highly flexible learners like "boosted regression trees" and random effects.

Until now, boosting could only be used for univariate regression and classification tasks. The key contribution of the [Rankin manuscript](https://drive.google.com/open?id=0BxeoeRy1g2juWVl4Yk8yZ2J2V1E) is make boosting possible for time-series and capture-histories data, by embedding a Expectation-Maximization sub-algorithm within the boosting gradient descent. 

## Benefits and Drawbacks

* automatic variable selection, detection of interactions, and non-linear effects;
* avoids overfitting through regularization and minimization of the _prediction error_ (aka the expected loss);
* shrinkage of unimportant covariates to zero;
* a type of multi-model inference, built-up in a stepwise manner. This is easier than having to fit and average over dozens or thousands of plausible fixed-effects models (like in model-averaging);
* avoids dubious boundary-value estimates (like 100% survival or 100% detection) because of shrinkage of extreme values;
* negotiates the "bias-variance trade-off" by incurring a little bias, but estimates are more stable over multiple realizations of the data and outliers;
* inference based on prediction, rather than dubious p-values or 95%CI;
* unified framework for addressing many hot-topics in capture-recapture, such as spatial capture-recapture, individual heterogeneity, non-linear effects, by leveraging the variety of base-learners in the `mboost` R packages.

The main drawback of boosting is that it depends on regularization hyperparameters (`mstop`, `nu`) which constrain the model-complexity. These cannot be estimated by the data. Instead, one must use bootstrap or cross-validation to optimize them. Typically, one must perform about 7 or 10 seperate k-fold cross-validations, which may take hours on a modern computer. 

# TUTORIAL
Interested users should set-through the R code in the file `R_TUTORIAL_dipper.R`. The tutorial does the following:
* Parts 1-4 analyses the Dipper dataset via Least-Squares, Splines and boosted regression trees.
* Part 5 compares boosted estimates to AICc model-averaging (requires RMark)
* Part 6 demonstrates CJSboost's ability to weed out unimportant variables, but adding a bunch of fake data to the dipper dataset
* Part 7 is a high-dimensional example from the manuscript: it demonstrates CJSboost's ability to find a sparse true model out of 21 candidate covariates. It compares PLS base-learners and boosted-regression trees.
* Part 8 is a strictly academic example, comparing two CJSboosting algorithms (Expectation-Maximization vs. Monte-Carlo approximation), to show they produce (approximately) the same estimates.

See below for some example plots from this tutorial.

## R formula inference and base-learners

CJSboost uses the model-specification formula of the R package `mboost`. This requires an understanding of "base-learner" formula syntax. See ["Model-based Boosting in R: A Hands-on Tutorial Using the R Package mboost"](https://epub.ub.uni-muenchen.de/12754/), and the R `mboost` help.files.  For example, open R, are run the commands `library(mboost); ?bols` to read all about specifying models with the various base-learners. 

Some examples in the above table: lets say we have two covariates called `sex` and `time`. 

|Model         | RMark | CJSboost                                                    |
|--------------|:---------:|-----------------------------------------------------------|
| intercept only| `~1` | `~ bols(interc,intercept=FALSE)`|
| main effects | `~sex+time` | `~bols(interc,intercept=FALSE)+bols(timefactor,df=1)+bols(sex,df=1)`|
| interaction  | `~sex*time` | `~bols(interc,intercept=FALSE)+bols(timefactor,df=1)+bols(sex,df=1)+bols(sex,timefactor,df=1)+bols(timefactor,by=sex,df=1)`|


* The first row is the phi(dot) model.
* The second row includes both sex and time, but the sex effect is the same over all capture periods. Notice how in CJSboost, we must explicitly include a base-learner for the intercept, whereas in RMark, this is always implied. Also, notice that the Penalized Least Squares base-learners have their flexibility/complexity controlled with the `df=1` argument. 
* In the third row, there is a two way interaction between sex and time, such that the effect of time is different for each sex. Notice that in boosting, we must include the lower-level main effects, as well as the two way interaction (specified with the `by=` argument).

We can also do more exotic semi-parametric base-learners, like regression trees (`btree`) and splines (`bbs`).


|Model         | RMark | CJSboost                                                    |
|--------------|:---------:|-----------------------------------------------------------|
| time as a continous variable| NA | `~ bols(interc,intercept=FALSE)+bols(time,df=1)+bbs(time,center=TRUE,df=1)`|
| regression rees | `~sex*time` | `~bols(interc,intercept=FALSE)+btree(sex,timefactor)`|


Notice that for the spline, we decompose the spline into its intercept, a linear trend (2nd term) and finally the wiggly-component. This decomposition of higher-order effects into their lower-order constituents is necessary only if we are particularly interested in the degree to which covariates are linear vs. non-linear or have interactions vs. no interactions. 

## Typical CJSboost analysis

### Part 1: 
