args = commandArgs(trailingOnly=TRUE)
### optimization part ###
## pre-calculate some metrics for gradient
## args
## X: a list of covariate matrices corresponding to each response
## Y: a list of response vectors
## value
## XY: a list of matrices X^TY for each response
grad_prep <- function(X, Y){
	ll = length(Y)
	P = ncol(X[[1]])
	XY = matrix(0,P,ll)
	for(i in 1:ll){
		XY[,i] = t(X[[i]])%*%Y[[i]]/nrow(X[[i]])
	}
	XY
}

## helper function for generating cross-validation sets
## args
## N: number of sample size
## fold: number of folds
## values
## perm: a permutation of 1 to N
## idx: matrix of fold by 2 with first col being starting index and second col being ending index
cv_helper <- function(N, fold){
	valid_num = floor(N/fold)
	set.seed(123)
	perm = sample(1:N, size = N)
	idx1 = seq(1,N,valid_num)
	idx2 = c(idx1[-1]-1,N)
	list(perm=perm, idx=cbind(idx1,idx2))
}

## get the minimum and maximum of lambda searched in cross-validation of an elastic net model
## args
## lst: an object returned by glmnet
## value
## min_lam: smallest lambda searched in glmnet cross-validation
## max_lam: largest lambda searched in glmnet cross-validation
minmax_lambda <- function(lst){
	max_lam = max(unlist(lapply(lst, function(x){max(x$lambda)})))
	min_lam = min(unlist(lapply(lst, function(x){min(x$lambda)})))
	c(min_lam, max_lam)
}

## evaluate the performance of elastic net on each response
## args
## lst: a list of glmnet object (fitted elastic net model for each response)
## X_tune: a list of covariate matrices corresponding for each response (for tuning lambda)
## Y_tune: a list of response vectors (for tuning lambda)
## X_test: a list of covariate matrices corresponding for each response (for testing performance)
## Y_test: a list of response vectors (for testing performance)
## value
## lam: best performing lambda (on (X_tune,Y_tune)) for each response
## mse: list of matrices with each element being a matrix of predicted vs observed response
## est: estimated effect sizes for each response (B matrix)
elastic_net_mse <- function(lst, X_tune, Y_tune, X_test, Y_test){
	P = length(lst)
	M = ncol(X_tune[[1]])
	lam_V = rep(0, P)
	test_res = list()
	test_beta = matrix(0, M, P)
	for(t in 1:P){
		ncv = length(lst[[t]]$lambda)
		tmp_mse = rep(0, ncv)
		for(k in 1:ncv){
			tmp_mse[k] = mean((Y_tune[[t]] - X_tune[[t]]%*%lst[[t]]$glmnet.fit$beta[,k])^2)
		}
		ss = which.min(tmp_mse)
		test_beta[,t] = lst[[t]]$glmnet.fit$beta[,ss]
		lam_V[t] = lst[[t]]$lambda[ss]
		predicted = X_test[[t]]%*%lst[[t]]$glmnet.fit$beta[,ss]
		test_res[[t]] = cbind(Y_test[[t]], predicted)
	}
	list(lam = lam_V, mse = test_res, est = test_beta)
}

## create a list of Y_true and Y_predicted for analysis
multi_mse <- function(theta_est, X_test, Y_test){
  answer = list()
	P = ncol(theta_est)
	for(t in 1:P){
		predicted = X_test[[t]]%*%theta_est[,t]
		answer[[t]] = cbind(Y_test[[t]], predicted)
	}
	answer
}

## average prediction metrics (rsq, mse, adjusted mse on validation sets)
avg_perm <- function(mse_lst){
	fd = length(mse_lst)
	P = length(mse_lst[[1]])
	rsq = mse = adj_mse = matrix(0, fd, P)
	for(f in 1:fd){
		for(t in 1:P){
			rsq[f,t] = (cor(mse_lst[[f]][[t]])[1,2])^2
			mse[f,t] = mean((mse_lst[[f]][[t]][,1]-mse_lst[[f]][[t]][,2])^2)
			adj_mse[f,t] = mse[f,t]/var(mse_lst[[f]][[t]][,1])
		}
	}
	cbind(apply(rsq, 2, mean), apply(mse, 2, mean), apply(adj_mse, 2, mean))
}

## test Y_true ~ Y_predicted
pred_test <- function(Y){
	if(sum(Y[,2]==0)==nrow(Y)|var(Y[,2])==0){
		return(2)
	}else{
		summary(lm(Y[,1]~Y[,2]))$coefficients[2,4]
	}
}

## group lasso on data with missing covariates, use validation data for stopping criterion
glasso <- function(X, Y, X1, Y1, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3, verbose = FALSE){
	Bgt = Sys.time()
	M = nrow(XY)
	P = length(X)
	NN = unlist(lapply(X, nrow))
	old_objV1 = rep(0,P)
	for(t in 1:P){
		old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
	}
  if (verbose) {
    cat("Training error: ", old_objV1, '\n')	
  }
	old_objV2 = rep(0,P)
	for(t in 1:P){
		old_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
	}
	if (verbose) {
    cat("Testing error: ", old_objV2, '\n')
  }
	beta_j_lasso = rep(0, P)
	tmp_XYj = 0
	if(!is.loaded("wrapper")){
		dyn.load("optim.so") # change this to the abs path to optim.so
	}
	for(i in 1:maxiter){
		bgt = Sys.time()
		res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
		edt = Sys.time()
		
		new_objV1 = new_objV2 = rep(0,P)
		for(t in 1:P){
			new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
		}
		if (verbose) {cat("Training error: ", new_objV1, '\n')}
		for(t in 1:P){
			new_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
		}
		if (verbose) {cat("Testing error: ", new_objV2, '\n')}
		if(mean(new_objV2) > mean(old_objV2)|mean(new_objV1) > mean(old_objV1)){
			break
		}else{
			old_objV2 = new_objV2
		}
		if(max(abs(new_objV1-old_objV1)) < eps){
			break
		}else{
			old_objV1 = new_objV1
		}
	}
	Edt = Sys.time()
	cat("total training time: ", Edt-Bgt, "\n")
	list(est = theta, avg_tune_err = mean(new_objV2), tune_err=new_objV2)
}

## simpler version of glasso, train model until converges 
glasso_no_early_stopping <- function(X, Y, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3, verbose = FALSE){
  cat("running glasso_no_early_stopping\n")
	M = nrow(XY)
	P = length(X)
	NN = unlist(lapply(X, nrow))
	old_objV1 = rep(0,P)
	for(t in 1:P){
		old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
	}
	if (verbose) {cat("Training error: ", mean(old_objV1), '\n')}
	beta_j_lasso = rep(0, P)
	tmp_XYj = 0
	if(!is.loaded("wrapper")){
		dyn.load("optim.so")
	}
	for(i in 1:maxiter){
		res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
		new_objV1 = rep(0,P)
		for(t in 1:P){
			new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
		}
		if (verbose) {cat("Training error: ", mean(new_objV1), '\n')}
		if(max(abs(new_objV1-old_objV1)) < eps|mean(new_objV1) > mean(old_objV1)){
			break
		}else{
			old_objV1 = new_objV1
		}
	}
	list(est = theta, avg_train_err = mean(new_objV1), train_err = new_objV1)
}

#----------------------- command line input ----------------------###
### data import ###
options(stringsAsFactors=F)
library(glmnet)
library(foreach)

dose_path = normalizePath(args[1])
info_path = normalizePath(args[2])
Yt_path = normalizePath(args[3])
Yt = dir(args[3])
ntune = as.numeric(args[4])
gene_id = args[5]
outdir = normalizePath(args[6])
if (is.na(args[7])) {
  if_verbose = FALSE
} else {
  if_verbose = as.logical(args[7])
}
sprintf("path to feature matrix %s; (id, ftr1, ftr2, ..., ftrM)", dose_path)
sprintf("path to info file: %s; ()", info_path)
sprintf("path to response files: %s (one response one file, two cols: id, value)", Yt_path)
sprintf("ntune = %s; (number of tuning parameters grid)", ntune)
sprintf("prefix for all output files: %s", gene_id)
sprintf("output dir: %s", outdir)
sprintf("print detailed log info: %s", if_verbose)

P = length(Yt)
fold = 5
sprintf("total number of responses = %s", P)
sprintf("number of cross-validation fold: %s", fold)
sprintf("using %s fold-data for training, %s fold-data for tuning (validation data) 
        and %s fold-data for evaluating performance (test data)", fold - 2, 1, 1)

if(P){
	dir.create(outdir, showWarnings = FALSE)
	#setwd(outdir)
  #-------------------------- load response files ---------------------------#
  cat("loading response files ...\n")
  # each response corresponds to a single file with two columns: id and value #
	Y = list()
	for(t in 1:P){
		Y[[t]] = read.table(paste0(Yt_path, '/', Yt[t]), header=F)
	}
	ssize = unlist(lapply(Y, nrow))
	T_num = length(Yt)
  cat("print head of one Yt:\n")
  print(head(Y[[1]]))

	#-------------------------- load feature files  ---------------------------#
  # i.e. dosage files for GTEx data
  cat("loading X matrix ...\n")
	dose = read.table(dose_path, header = F)
  cat("print dim of one Xt:\n")
  print(dim(dose))
  cat("centering each column of X\n")
	for(j in 2:ncol(dose)){ 
    ## GTEx dosage file has an extra 'dose' column, change 2 to 3 if you load GTEx
    ## dosage file for X
		dose[,j] = dose[,j] - mean(dose[,j])
	}
	N = nrow(dose)
	#-------------------------- computing covariance matrix ---------------------------#
	cat("computing covariance matrix using all X ...\n")
  tmp = as.matrix(dose[ , -1])
	XX = t(tmp)%*%as.matrix(tmp)/N
  print(dim(XX))
	Xnorm = diag(XX)
	remove(tmp); remove(XX)
  # gtex dosage file has subject_id->sample_id format
	# sub_id = matrix(unlist(strsplit(dose[,1], "->")), ncol=2, byrow=T)[,1] 
  sub_id = dose[, 1]
	M = ncol(dose) - 1
  cat('M = ', M, '\n')
	sub_id_map = list()
	for(t in 1:T_num){
		tmp = rep(0, nrow(Y[[t]]))
		for(j in 1:length(tmp)){
			tmp[j] = which(sub_id == Y[[t]][j,1])
		}
		sub_id_map[[t]] = tmp
	}
	cv_config = cv_helper(N, fold)
	cv_perm = cv_config$perm
	cv_idx = cv_config$idx
	
	single_res_test = list()
	single_lam = matrix(0,fold,P)
	single_theta_est = list()
	
	multi_res_test = list()
	multi_lam = matrix(0,fold,2)
	multi_theta_est = list()

	multi_res_test2 = list()
	multi_lam2 = array(0, dim=c(fold, P, 2))
	multi_theta_est2 = list()
  
  #-------------------------- dividing cross-validation sets ---------------------------#
	cat('dividing data into train-validate-test sets\n')
  res_tune = list()
	rec_lamv = matrix(0, fold, ntune)
	for(f in 1:fold){
		bgt = Sys.time()
		test_index = cv_perm[cv_idx[f,1]:cv_idx[f,2]]
		test_id = sub_id[test_index]
		tuning_index = cv_perm[cv_idx[f%%fold+1,1]:cv_idx[f%%fold+1,2]]
		tuning_id = sub_id[tuning_index]
	
		X_test = list()
		Y_test = list()
		X_tune = list()
		Y_tune = list()
		X_train = list()
		Y_train = list()
		for(t in 1:T_num){
			X_train_tmp = sub_id_map[[t]][!(sub_id_map[[t]]%in%c(test_index,tuning_index))]
			Y_train_tmp = !(sub_id_map[[t]]%in%c(test_index,tuning_index))
			X_tuning_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%tuning_index)]
			Y_tuning_tmp = (sub_id_map[[t]]%in%tuning_index)
			X_test_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%test_index)]
			Y_test_tmp = (sub_id_map[[t]]%in%test_index)
			X_train[[t]] = apply(as.matrix(dose[X_train_tmp,-c(1)]),2,as.numeric)
			Y_train[[t]] = Y[[t]][Y_train_tmp, 2]
			X_tune[[t]] = apply(as.matrix(dose[X_tuning_tmp,-c(1)]),2,as.numeric)
			Y_tune[[t]] = Y[[t]][Y_tuning_tmp, 2]
			X_test[[t]] = apply(as.matrix(dose[X_test_tmp,-c(1)]),2,as.numeric)
			Y_test[[t]] = Y[[t]][Y_test_tmp, 2]
		}
		
		#------------- get initial est by elasticNet on single tissue --------------#
    ## model training ##	
		## train elastic net and used average lambda as tuning parameters ##
		cat("training elastic net on single response to get initial est\n")
    single_initial_est = matrix(0, ncol(X_train[[1]]), T_num)
		single_summary = list()
		for(t in 1:T_num){
			tt = cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5)
			single_summary[[t]] = tt
			single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
		}
		## performance of Elastic net on tuning and testing data with various tuning parameters
		els_output = elastic_net_mse(single_summary, X_tune, Y_tune, X_test, Y_test)
		single_res_test[[f]] = els_output$mse
		single_lam[f,] = els_output$lam
		single_theta_est[[f]] = els_output$est
		remove(els_output)
		
		## use elastic net ests row norm as weights ##
		lam_range = minmax_lambda(single_summary)
		sig_norm = apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
		sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
		sig_norm = sig_norm/sum(sig_norm)
		weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);

		tis_norm = apply(single_initial_est, 2, function(x){sum(abs(x))})
		tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
		tis_norm = tis_norm/sum(tis_norm)
		weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);
		lam_V = seq(lam_range[1], lam_range[2], length.out = ntune)
		#lam_V = seq(lam_range[1], lam_range[2], length.out = ntune)
		

		initial_numeric = as.numeric(single_initial_est)
		remove(single_summary); remove(single_initial_est);
	
		#-------------------------- train - validate - test ---------------------------#
    cat('starting train-validate-test\n')
    ## preparation
		XY = grad_prep(X_train, Y_train)
		XX_train = lapply(X_train, function(x){t(x)%*%x/nrow(x)})
		spsz = unlist(lapply(X_train,nrow))
		#res_tune = rep(0, ntune)
		res_tune[[f]] = array(-1, dim=c(ntune, ntune, P))
		#best.lam = 0
		rec_lamv[f,] = lam_V
		for(lam1 in 1:ntune){
			for(lam2 in 1:ntune){
				single_est = matrix(initial_numeric, M, P)
				ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[lam1]/spsz, lambda2=lam_V[lam2], theta=single_est, verbose = if_verbose)
				if(sum(ans$est!=0)>0){
					res_tune[[f]][lam1,lam2, ] = ans$tune_err
					if (if_verbose) { cat("lambda1=",lam_V[lam1], "; lambda2=", lam_V[lam2], "; avg tune err=", ans$avg_tune_err, '\n') }
					remove(single_est); remove(ans);
				}else{
					remove(single_est); remove(ans);
					break
				}			
			}
		}

    #-------------------------- save results on test set for evaluation ---------------------------#
    cat("saving cross-validation results for evaluation and analysis\n")
		avg_tune_res = apply(res_tune[[f]], c(1,2), mean)
		best.lam = which(avg_tune_res == min(avg_tune_res[avg_tune_res>=0]), arr.ind = TRUE)[1,]
		single_est = matrix(initial_numeric, M, P)
		ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam[1]]/spsz, lambda2=lam_V[best.lam[2]], theta=single_est, verbose = if_verbose)
		multi_res_test[[f]] = multi_mse(ans$est, X_test, Y_test)
		multi_lam[f,] = lam_V[best.lam]
		multi_theta_est[[f]] = ans$est
		remove(single_est); remove(ans);
		
    ## tune by each tissue ##
		#tmp_est = matrix(0, M, P)
		#for(t in 1:P){
		#	sig_tune_res = res_tune[[f]][,,t]
		#	best.lam.single = which(sig_tune_res == min(sig_tune_res[sig_tune_res>=0]), arr.ind = TRUE)[1,]
		#	multi_lam2[f,t,] = lam_V[best.lam.single]
		#	ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam.single[1]]/sig_norm, lambda2=lam_V[best.lam.single[2]]/sig_norm, theta=matrix(initial_numeric, M, P), verbose = if_verbose)
		#	tmp_est[,t] = ans$est[,t]
		#	remove(ans);
		#}
		#multi_res_test2[[f]] = multi_mse(tmp_est, X_test, Y_test)
		#multi_theta_est2[[f]] = tmp_est
		edt = Sys.time()
		cat("cross-validation time for fold ", f, " is ", edt-bgt, '\n')
	}
	save(single_res_test, single_lam, single_theta_est, multi_res_test, multi_lam, multi_theta_est, res_tune, rec_lamv, file = paste0(outdir, '/', gene_id, ".cv.evaluation.RData"))
	#res_single = avg_perm(single_res_test)
	#res_multi = avg_perm(multi_res_test)
	#cat("Elastic net average testing error (all): ", apply(res_single, 2, mean), '\n')
	#cat("glasso averge testing error (all): ", apply(res_multi, 2, mean), '\n')
	#cat("Number of all zero tissues in elastic net is ", sum(is.na(res_single[,1])), '\n')
	#cat("Number of all zero tissues in glasso is ", sum(is.na(res_multi[,1])), '\n')
	#cat("Elastic net average testing error (non-zero): ", apply(res_single[!is.na(res_multi[,1]),], 2, mean), '\n')
	#cat("glasso averge testing error (non-zero): ", apply(res_multi[!is.na(res_multi[,1]),], 2, mean), '\n')

  #------------ use tuning parameter chosen above to train model on entire dataset -------------#
	## generate an estimate with whole data ##
  cat('training a model on entire data with parameters chosen from cv\n')
	X_all = list()
	Y_all = list()
	for(t in 1:T_num){
		X_all_tmp = sub_id_map[[t]]
		X_all[[t]] = apply(as.matrix(dose[X_all_tmp,-c(1)]),2,as.numeric)
		Y_all[[t]] = Y[[t]][,2]
	}
	# initial values 
	single_initial_est = matrix(0, ncol(X_train[[1]]), T_num)
	for(t in 1:T_num){
		tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5)
		single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
	}

	sig_norm = apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
	sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
	sig_norm = sig_norm/sum(sig_norm)
	weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);

	tis_norm = apply(single_initial_est, 2, function(x){sum(abs(x))})
	tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
	tis_norm = tis_norm/sum(tis_norm)
	weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);

	spsz = unlist(lapply(X_all,nrow))
	initial_numeric = as.numeric(single_initial_est)
	#remove(single_initial_est)
	XY = grad_prep(X_all, Y_all)
	XX_all = lapply(X_all, function(x){t(x)%*%x/nrow(x)})
	tmp_res = rep(0, fold)
	for(f in 1:fold){
		ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=multi_lam[f,1]/spsz, lambda2=multi_lam[f,2], theta=matrix(initial_numeric,M,P), verbose = if_verbose)
		tmp_res[f] = ans$avg_train_err
	}
	final.lam = multi_lam[which.min(tmp_res),]
	ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=final.lam[1]/spsz, lambda2=final.lam[2], theta=matrix(initial_numeric,M,P), verbose = if_verbose)
	info = read.table(info_path, header=T, sep='\t')
	downstream_est = data.frame(info[,1:3], ans$est)
	multi_all_res = multi_mse(ans$est, X_all, Y_all)
	single_all_res = multi_mse(single_initial_est, X_all, Y_all)
  cat('writing final estimates\n')
	write.table(downstream_est, paste0(outdir, '/', gene_id, ".est"), quote = F, row.names = F, col.names = c("SNP", "REF.0.", "ALT.1.", Yt))

	#tmp_res = matrix(0, fold, P)
	#for(f in 1:fold){
	#	for(t in 1:P){
	#		ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=multi_lam2[f,t,1]/sig_norm, lambda2=multi_lam2[f,t,2]/sig_norm, theta=matrix(initial_numeric,M,P))
	#		tmp_res[f,t] = ans$train_err[t]
	#	}
	#}
	#final.lam2 = matrix(0, P, 2)
	#tmp_est2 = matrix(0, M, P)
	#for(t in 1:P){
	#	final.lam2[t,] = multi_lam2[which.min(tmp_res[,t]),t,]
	#	ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=final.lam2[t,1]/sig_norm, lambda2=final.lam2[t,2]/sig_norm, theta=matrix(initial_numeric,M,P))
	#	tmp_est2[,t] = ans$est[,t]
	#}
	#downstream_est2 = data.frame(info[,1:3], tmp_est2)
	#write.table(downstream_est2, paste0(gene_id, "2.est"), quote=F, row.names=F, col.names=c("SNP", "REF.0.", "ALT.1.", Yt))

	#multi_all_res2 = multi_mse(tmp_est2, X_all, Y_all)
  cat('saving the prediction on all data for future analysis\n')
	save(multi_all_res, single_all_res, final.lam, ans, file = paste0(outdir, '/', gene_id, ".prediction_on_all_data.RData"))
  cat('done!\n')
}

