\name{MF.beta4-package}
\alias{MF.beta4-package}
\title{Measuring ecosystem multifunctionality and assessing BEF relationships}

\description{

\code{MF.beta4} is a R package for measuring ecosystem multifunctionality and assessing BEF relationships. It was originally developed for the Beta4 project (Müller et al. 2022) on the effect of enhancing the beta diversity between forest patches on ecosystem multifunctionality and forest resilience across spatial scales.    \cr


Based on a framework of Hill-Chao numbers of orders q = 0, 1 and 2, MF.beta4 features the following multifunctionality measures for a single and multiple ecosystems; see Chao et al. (2023) for pertinent methodology and decomposition theory. \cr

(1)	Multifunctionality measures in a single ecosystem: MF.beta4 computes a class of weighted multifunctionality measures for given function weights. Multifunctionality measures that correct for strong correlations between ecosystem functions to avoid redundancy are also provided. \cr

(2)	Multifunctionality measures in multiple ecosystems: MF.beta4 computes gamma multifunctionality of pooled ecosystems, within-ecosystem component (alpha multifunctionality) and among-ecosystem component (beta multifunctionality).    \cr

Based on biodiversity and function data from ecosystems, this package also provides graphics for assessing the biodiversity-ecosystem functioning (BEF) relationships across scales by `ggMF`.   \cr


This package contains four main functions:   \cr

1. \code{func_normalized} normalizing ecosystem functions data and turn all values between 0 and 1. \cr 

2. \code{MF_single} computing quantification of multi-fuctionality measures in a single ecosystem when all functions are assumed to be independent, and when the correlation between any two functions is adjusted for. In the latter case, our MF value is obtained based on an integrated measure by considering all possible threshold levels. \cr

3. \code{MF_multiple} computing decomposition of multi-fuctionality measures in multiple ecosystems when all functions are assumed to be independent, and when the correlation between any two functions is adjusted for. In the latter case, our MF value is obtained based on an integrated measure by considering all possible threshold levels. \cr

4. \code{ggMF} visualizing the output from the function \code{MF_single} or \code{MF_multiple}.

}

\author{
Anne Chao, Chun-Yu Liu  \cr
Maintainer: Anne Chao <chao@stat.nthu.edu.tw> 
}

\references{
* Chao, A., Chiu, C. H., Hu, K. H., van der Plas, F., Cadotte, M. W., Mitesser, O., et al. (2023). Hill-Chao numbers in multifunctionality allows decomposing gamma multifunctionality into alpha and beta components. Under revision, Ecology Letters.  

* Chao, A. and Hu, K.-H. (2023) MF.beta4: Software for measuring ecosystem multifunctionality and assessing BEF relationships. Code and User’s guide available from \cr 
http://chao.stat.nthu.edu.tw/wordpress/software_download/MF-beta4/.

}
