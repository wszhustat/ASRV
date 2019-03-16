 \name{ASRV}
 \alias{ASRV}
 \alias{ASRV.SSD.OneSet}
 \alias{ASRV.SSD.OneSet_SetIndex}
 \title{Association Test for Rare Variants Based on Algebraic Statistics}
 \description{
     Test for associations between rare variants and a dichotomous phenotype based on an exact test (a generalization of Fisher's exact test) in a two-way contingency table. The p-value of the exact test is computed based on Monte Carlo approach by sampling observations on a fiber from a Markov chain constructed based on Markov bases in the framework of algebraic statistics.
 }
 \usage{

ASRV(genoTypeData, phenoTypeStatus, nDropOut = 100000, nReplicates = 1000000, isCheckParam = TRUE)


 }
\arguments{
      \item{genoTypeData}{a numeric matrix or data frame with each row as a different individual and each column as a separate SNP. The genotype on each SNP should be coded as 0, 1, 2. No missing data allowed.}
      \item{phenoTypeStatus}{a numeric vector with phenotype status: 0=controls, 1=cases. No missing data allowed.}
      \item{nDropOut}{an integer specifying the number of drop-out samples , which is out of interest in the Markov chain Monte Carlo sampling.}
      \item{nReplicates}{an integer specifying the number of replicates used in the Markov chain Monte Carlo sampling. }
      \item{isCheckParam}{a logical value indicating whether to check the validity of the phenotype status and genotype data (default= TRUE). If TRUE, it checks if genoTypeData and phenoTypeStatus are within its given range.}

}


\value{
	\item{p.value}{p-value of ASRV. }

}


\author{Jingbo Meng, Wensheng Zhu, Canhui Li, and KyongSon Jon}

\references{
Agresti, A. (1992). A survey of exact inference for contingency tables. \emph{Statistical science}, 7(1):131-153.


Darroch, J. N. (1980). Markov fields and log-linear interaction models for contingency tables. \emph{The Annals of Statistics}, 8(3):522-539.


Diaconis, P. and Sturmfels, B. (1998). Algebraic algorithms for sampling from conditional
distributions. \emph{The Annals of Statistics}, 26:363-397.


Dinwoodie, I. H. (1998). The diaconis-sturmfels algorithm and rules of succession. \emph{Bernoulli},
4(3):401-410.

Pachter, L. and Sturmfels, B. (2005). Algebraic Statistics for Computational Biology. \emph{Cambridge University Press}.


Rapallo, F. (2003). Algebraic markov bases and mcmc for two-way contingency tables.
\emph{Scandinavian journal of statistics}, 30(2):385-397.
}

\examples{


 ## Not run:

  # number of cases
  cases = 500

  # number of controls
  controls = 500

  # total (cases + controls)
  total = cases + controls

  # phenotype vector
  phenotype = c(rep(1, cases), rep(0, controls))

  # genotype matrix with 10 variants (random data)
  set.seed(1234)
  genotype = matrix(rbinom(total * 10, 2, 0.03), nrow = total, ncol = 10)

  # apply ASRV test
  myasrv = ASRV(genotype, phenotype, 1000, 100000, TRUE)
  myasrv
## End(Not run)




}


