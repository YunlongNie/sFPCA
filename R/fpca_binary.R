#' This function computes the sfpc when the response variable is binary. 
#'
#' @param xmat a matrix. Each column represents a sample. 
#' @param y a vector either 0 or 1. The binary response variable for each sample. 
#' @param npc_select an integer. The number of FPCs required. 
#' @param xmat_new a matrix. New dataset wish to predict on, each col corresponds to one sample.
#' @param lambda a positive number. The smoothing parameter.
#' @param theta between 0 and 1. The weight parameter.
#' @export
#' @import fda
#' @import dplyr
#' @examples
#' \dontrun{
#' data(binary_dat);xmat = binary_dat$x;y=binary_dat$y
#' res = sfpcs_binary(xmat,y,npc_select=2,theta=1,xmat_new = xmat)
#' res$fitted
#' mean(y==res$fitted) # fitted accuracy
#' res$fitted # prediction on new data 
#' res$sfpcs # sFPCs 
#' res$beta_fd # coefficient function
#' }

sfpcs_binary = function(xmat,y,npc_select=2,theta= 0.5, lambda = 10,timepts=NULL,xmat_new = NULL)
{

stopifnot(length(unique(y))==2) 
stopifnot(all(c("fda","dplyr")%in%rownames(installed.packages()))) 
train_matrix = xmat;
rowM = rowMeans(train_matrix);train_matrix = sweep(train_matrix,1,rowM)
train_y = y

if (is.null(timepts)) timepts=1:nrow(train_matrix);
norder=4 ## cubic B-spline
nbasis=norder+length(timepts)-2; (nbasis)
spline.basis=create.bspline.basis(rangeval=c(1,nrow(train_matrix)),nbasis,norder,timepts)
D2Lfd <- int2Lfd(m=2)
D2fdPar     <- fdPar(spline.basis, D2Lfd,1)
train.fd <- Data2fd(y=train_matrix, argvals=timepts,D2fdPar,nderiv=2)
if(!is.null(xmat_new)) 
{
 		test_matrix = xmat_new;test_matrix =sweep(test_matrix,1,rowM)
		test.fd <- Data2fd(y=test_matrix, argvals=timepts,D2fdPar,nderiv=2,lambda=1e2)
}


W = inprod(spline.basis,spline.basis)
### matrix S the coef matrix for the observed curves 
S = train.fd$coef 
M1 = W%*%S%*%as.matrix(train_y)
M2 = W%*%S%*%(matrix(1,nrow=length(train_y),ncol=1)-as.matrix(train_y))
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
MM = M1 %*% t(M1)/sum(train_y) + M2 %*% t(M2)/(length(train_y) - 
            sum(train_y))
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

formula_pc = as.formula(paste0("y~",paste0(paste0("score",1:npc_select),collapse="+")))

glm_fit = glm(formula_pc,scores_res_train,family="binomial")
fitted = ifelse(predict(glm_fit,newdata=scores_res_train,type="link")>0,1,0)


predicted_y  = NULL
if(!is.null(xmat_new)) {
	scores_res_test = do.call(cbind.data.frame,lapply(fpca_res$har,function(x) as.numeric(inprod(x,test.fd))))
	names(scores_res_test) = paste0("score",1:npc_select);
	predicted_y = ifelse(predict(glm_fit,newdata=scores_res_test,type="link")>0,1,0)
}


fpcs = fpca_res$har
betas = glm_fit%>%coef
beta0 = as.vector(betas[1])
beta_func = as.vector(betas[-1])
beta_fd = beta_func[1]*fpcs[[1]]
for (i in 2:length(beta_func))
{
    beta_fd  = beta_fd + beta_func[i]*fpcs[[i]]
}


return(list( intercept = beta0,  beta_fd=beta_fd,fitted=fitted,predicted = predicted_y, sfpcs=fpcs,theta=theta,lambda=lambda))
        

}



