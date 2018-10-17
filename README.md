# CTIMP
**C**ross **T**issue gene expression **IMP**utation

## About
CTIMP uses a multivariate-response penalized regression to predict cross-tissue gene expression. By imposing group-lasso penalty on each row of coefficient matrix, CTIMP reaches higher accuracy and robust prediction over penalized regression trained for each response separately. The approximation procedure proposed for coordinate descent iteration makes training faster and yields small approximation error. Details of CTIMP could be found in our [UTMOST manuscript](https://www.biorxiv.org/content/early/2018/03/21/286013).

The repo also contains codes for replicating analysis in [UTMOST manuscript](https://www.biorxiv.org/content/early/2018/03/21/286013) and training cross-tissue gene expression prediction model with GTEx data as a demonstration.

## Tutorial
Clone repo
```bash
git clone https://github.com/yiminghu/CTIMP.git
```
Compile the optim.c
```bash
cd PERMIT
R CMD SHLIB optim.c
```

### Input format
X: covariate matrix with each row being an observation and each column being a variable. The first column is the ID of each observation. No header line.

Y_1, Y_2, ..., Y_P: response vectors, each file contains two columns: ID and response values. Due to the missingness in reponse matrix, each file may have different number of rows (number of observations).

info: files contain column infomation of X, see example/ for details

### Usage
```bash
X=example/X.txt
Y_folder=example/Y_folder/
info=example/info.txt
ntune=50 #number of grids for each tuning parameter
output_path=example/output/
mkdir ${output_path}
output_predix=test # prefix of output files

Rscript --vanilla main.R ${X} ${info} ${Y_folder} ${ntune} ${output_prefix} ${output_path}
```
main.R will perform 5-fold cross-validation to select tuning parameters and generate estimation of the coefficient matrix.
