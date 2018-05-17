# PERMIT
**PE**nalized multivariate **R**esponse regression with **M**issingness **I**n response ma**T**rix.

## About
PERMIT is a penalized regression for multivariate-response data with missingness in response matrix. By imposing group-lasso penalty on each row of coefficient matrix, PERMIT reaches more accurate and robust prediction over penalized regression trained for each response separately. The approximation procedure proposed for coordinate descent iteration makes training faster and yields small approximation error. Details of PERMIT could be found in our [UTMOST manuscript](https://www.biorxiv.org/content/early/2018/03/21/286013).

The repo also contains codes for replicating analysis in [UTMOST manuscript](https://www.biorxiv.org/content/early/2018/03/21/286013) and training cross-tissue gene expression prediction model with GTEx data as a demonstration.

## Tutorial
1. Clone repo
```bash
git clone https://github.com/yiminghu/PERMIT.git
```
Compile the optim.c
```bash
cd PERMIT
R CMD SHLIB optim.c
```