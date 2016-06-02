#' This function computes the sfpc when the response variable is continuous. 
#'
#' @param xmat a matrix. Each column represents a sample. 
#' @param y a vector. The response variable for each sample. 
#' @param npc_select an integer. The number of FPCs required. 
#' @param xmat_new a matrix. New dataset wish to predict on, each col corresponds to one sample. 
#' @param lambda a positive number. The smoothing parameter.
#' @param theta between 0 and 1. The weight parameter.
#' @export
#' @import fda
#' @import dplyr
#' @examples
#' \dontrun{
#' data(medfly)
#' xmat = medfly$eggcount  # 26 days and 50 flies 
#' y = medfly$lifetime
#' res = sfpca(xmat,y,0.1,1e3,npc_select=5,xmat_new = xmat[,1:10])
#' plot(res$beta_fd,ylab="beta(t)")
#' res$fitted
#' res$predicted 
#' res$sfpcs
#' }

sfpca_con = function(xmat,y,theta,lambda,npc_select = 3,timepts=NULL,xmat_new=NULL)
{
stopifnot(length(unique(y))<=2) 
stopifnot(all(c("fda","dplyr")%in%rownames(installed.packages()))) 
train_matrix =xmat;rowM = rowMeans(train_matrix);train_matrix =sweep(train_matrix,1,rowM)
train_y = y; meany = mean(train_y);train_y = (train_y-meany)
if (is.null(timepts)) timepts=1:nrow(xmat)
norder=4 ## cubic B-spline
nbasis=norder+length(timepts)-2; (nbasis)
spline.basis=create.bspline.basis(rangeval=range(timepts),nbasis,norder,timepts)
D2Lfd <- int2Lfd(m=2)
D2fdPar     <- fdPar(spline.basis, D2Lfd, 1e2)
train.fd <- Data2fd(y=train_matrix, argvals=timepts,D2fdPar,nderiv=2,lambda=1e2)


if(!is.null(xmat_new)) 
	{
		test_matrix = xmat_new;test_matrix =sweep(test_matrix,1,rowM)
		test.fd <- Data2fd(y=test_matrix, argvals=timepts,D2fdPar,nderiv=2,lambda=1e2)

	}



W = inprod(spline.basis,spline.basis)
S = train.fd$coef 
maty = matrix(rep(train_y,each=nrow(S)),nrow=nrow(S))
M = rowSums(maty*W%*%S)
MM = as.matrix(M)%*%t(as.matrix(M))
sqrM = function (X) 
{
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
}

fpca = function(theta,lambda,npc)
{

Q =theta/length(train_y)*W%*%S%*%t(S)%*%W+(1-theta)*MM/(length(train_y)^2)

D = inprod(spline.basis,spline.basis,Lfdobj1=2,Lfdobj2=2)
G = W + lambda*D

halfG_inv  = sqrM(G)

tM = t(halfG_inv)%*%Q%*%halfG_inv
eigen_res = eigen(tM)
varprp  = round(head(eigen_res$values/sum(eigen_res$values),npc)*100,4)

fd_list = lapply(1:npc,function(ipc)
{
cat(ipc)
coef_pc= halfG_inv%*%as.matrix(eigen_res$vectors[,ipc])
fd= fd(coef=coef_pc,spline.basis)
fd
	})
list(har = fd_list,varprp = varprp)
}

fpca_res = fpca(theta,lambda,npc=npc_select)
## train scores
scores_res_train = do.call(cbind.data.frame,lapply(fpca_res$har,function(x) as.numeric(inprod(x,train.fd))))
names(scores_res_train) = paste0("score",1:npc_select);scores_res_train$theta=theta;scores_res_train$lambda=lambda;scores_res_train$y=train_y
## test scores

fpcs = fpca_res$har
fpcs_varprp = fpca_res$varprp

betas = scores_res_train%>%data.frame%>%select(matches("score"),y)%>%lm(y~.,data=.)%>%coef
beta0 = as.vector(betas[1])
beta_func = as.vector(betas[-1])

beta_fd = beta_func[1]*fpcs[[1]]
for (i in 2:length(beta_func))
{
	beta_fd  = beta_fd + beta_func[i]*fpcs[[i]]
}

intercept = beta0 + meany
fitted_y =  inprod(beta_fd,train.fd) + intercept
predicted_y = NULL
if(!is.null(xmat_new)) {
	predicted_y =  inprod(beta_fd,test.fd) + intercept
}


return(list(intercept =  intercept, beta_fd = beta_fd, fitted = fitted_y, predicted = predicted_y, sfpcs=fpcs,sfpcs_varprp=fpcs_varprp,theta=theta,lambda=lambda))
}



