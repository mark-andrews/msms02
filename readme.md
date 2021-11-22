# Model selection and model simplification

This two day course covers the important and general topics of statistical model building, model evaluation, model selection, model comparison, model simplification, and model averaging.
These topics are vitally important to almost every type of statistical analysis, yet these topics are often poorly or incompletely understood.
We begin by considering the fundamental issue of how to measure model fit and a model's predictive performance, and discuss a wide range of other major model fit measurement concepts like likelihood, log likelihood, deviance, residual sums of squares etc.
We then turn to nested model comparison, particularly in general and generalized linear models, and their mixed effects counterparts.
We then consider the key concept of out-of-sample predictive performance, and discuss over-fitting or how excellent fits to the observed data can lead to very poor generalization performance.
As part of this discussion of out-of-sample generalization, we introduce leave-one-out cross-validation and Akaike Information Criterion (AIC).
We then cover general concepts and methods related to variable selection, including stepwise regression, ridge regression, Lasso, and elastic nets.
Following this, we turn to model averaging, which is an arguably always preferable alternative to model selection.
Finally, we cover Bayesian methods of model comparison.
Here, we describe how Bayesian methods allow us to easily compare completely distinct statistical models using a common metric.
We also describe how Bayesian methods allow us to fit all the candidate models of potential interest, including cases were traditional methods fail.

## Intended Audience

This course is aimed at anyone who is interested in advanced statistical modelling as it is practiced widely throughout academic scientific research, as well as widely throughout the public and private sectors.

## Teaching Format

This course will be largely practical, hands-on, and workshop based. For each topic, there will first be some lecture style presentation, i.e., using slides or blackboard, to introduce and explain key concepts and theories. Then, we will cover how to perform the various statistical analyses using R. Any code that the instructor produces during these sessions will be uploaded to a publicly available GitHub site after each session. 

The course will take place online using Zoom. On each day, the live video broadcasts will occur between (UK local time) at:

* 10am-12pm
* 1pm-3pm
* 4pm-6pm
 
All the sessions will be video recorded, and made available immediately on a private video hosting website as soon as possible after each 2hr session. 

Although not strictly required, using a large monitor or preferably even a second monitor will make the learning experience better, as you will be able to see my RStudio and your own RStudio simultaneously. 

## Assumed quantitative knowledge

We will assume familiarity with general statistical concepts, linear models, statistical inference (p-values, confidence intervals, etc). Anyone who has taken undergraduate (Bachelor's) level introductory courses on (applied) statistics can be assumed to have sufficient familiarity with these concepts.

## Assumed computer background

Prior experience with R and RStudio is required. Attendees should be familiar with some basic R syntax and commands, how to write code in the RStudio console and script editor, how to load up data from files, etc. 

## Equipment and software requirements

Attendees of the course will need to use RStudio. Most people will want to use their own computer on which they install the RStudio desktop software. This can be done Macs, Windows, and Linux, though not on tablets or other mobile devices. Instructions on how to install and configure all the required software, which is all free and open source, will be provided before the start of the course. We will also provide time at the beginning of the workshops to ensure that all software is installed and configured properly. An alternative to using a local installation of RStudio is to use RStudio cloud (https://rstudio.cloud/). This is a free to use and full featured web based RStudio. It is not suitable for computationally intensive work but everything done in this class can be done using RStudio cloud. 

In addition, we will use a set of additional R packages. 
Instructions on how to install all the software is [here](software.md).

An RStudio server session with Stan and `brms` installed and ready to use is available by clicking this button.
[![Binder](https://notebooks.gesis.org/binder/badge_logo.svg)](https://notebooks.gesis.org/binder/v2/gh/mark-andrews/intro_bda_qub/HEAD?urlpath=rstudio)

# Course programme 

## Day 1 

* Topic 1: *Measuring model fit*. In order to introduce the general topic of model evaluation, selection, comparison, etc., it is necessary to understand the fundamental issue of how we measure model fit. Here, the concept of conditional probability of the observed data, or of future data, is of vital importance. This is intimately related, though distinct, to concept of likelihood and the likelihood function, which is in turn related to the concept of the log likelihood or deviance of a model. Here, we also show how these concepts are related to concepts of residual sums of squares, root mean square error (rmse), and deviance residuals.
* Topic 2: *Nested model comparison*. In this section, we cover how to do nested model comparison in general linear models, generalized linear models, and their mixed effects (multilevel) counterparts. First, we precisely define what is meant by a nested model. Then we show how nested model comparison can be accomplished in general linear models with F tests, which we will also discuss in relation to $R^2$ and adjusted $R^2$. In generalized linear models, and mixed effects models, we can accomplish nested model comparison using deviance based chi-square tests via Wilks's theorem.
* Topic 3: *Out of sample predictive performance: cross validation and information criteria*. In the previous sections, the focus was largely on how well a model fits or predicts the *observed* data. For reasons that will be discussed in this section, related to the concept of overfitting, this can be a misleading and possibly even meaningless means of model evaluation. Here, we describe how to measure *out of sample* predictive performance, which measures how well a model can generalize to new data. This is arguably the gold-standard for evaluating any statistical models. A practical means to measure out of sample predictive performance is cross-validation, especially leave-one-out cross-validation. Leave-one-out cross-validation can, in relatively simple models, be approximated by *Akaike Information Criterion* (AIC), which can be exceptionally simple to calculate. We will discuss how to interpret AIC values, and describe other related information criteria, some of which will be used in more detail in later sections.

## Day 2

* Topic 4: *Variable selection*. Variable selection is a type of nested model comparison. It is also one of the most widely used model selection methods, and variable selection of some kind is almost always done routinely in all data analysis. Although we will also have discussed variable selection as part of Topic 2 above, we discuss the topic in more detail here. In particular, we cover stepwise regression (and its limitations), all subsets methods, ridge regression, Lasso, and elastic nets.
* Topic 5: *Model averaging*. Rather than selecting one model from a set of candidates, it is arguably always better perform model averaging, using all the candidates models, weighted by the predictive performance. We show how to perform model average using information criteria.
* Topic 6: *Bayesian model comparison methods*. Bayesian methods afford much greater flexibility and extensibility for model building than traditional methods. They also allow us to easily directly compare completely unrelated statistical models of the same data using information criteria such as WAIC and LOOIC. Here, we will also discuss how Bayesian methods allow us to fit all models of potential interest to us, including cases where model fitting is computationally intractable using traditional methods (e.g., where optimization convergence fails). This allows us therefore to consider all models of potential interest, rather than just focusing on a limited subset where the traditional fitting algorithms succeed.
