#' This function computes the sfpc when the response variable is binary. 
#'
#' @param xmat a matrix. Each column represents a sample. 
#' @param y a vector either 0 or 1. The binary response variable for each sample. 
#' @param npc_select an integer. The number of FPCs required. 
#' @param lambda a positive number. The smoothing parameter.
#' @param theta between 0 and 1. The weight parameter.
#' @export
#' @import fda
#' @import dplyr
#' @examples
#' \dontrun{
#' data(binary_dat);xmat = binary_dat$x;y=binary_dat$y
#' res = sfpcs_binary(xmat,y,npc=2,theta=1)
#' res$fitted
#' mean(y==res$fitted)
#' res$sfpcs
#' res$beta0
#' }
sfpcs_binary = function(xmat,y,npc=1,theta= 0.5, lambda = 10,timepts=NULL)
{
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

scores = function(theta,lambda,npc)
{

MM = M1%*%t(M1)/sum(train_y) + M2%*%t(M2)/(length(train_y)-sum(train_y))
U =theta/length(train_y)*W%*%S%*%t(S)%*%W+(1-theta)*MM

D = inprod(spline.basis,spline.basis,Lfdobj1=2,Lfdobj2=2)

G = W + lambda*D

halfG_inv  = sqrM(G)

tM = t(halfG_inv)%*%U%*%halfG_inv
eigen_res = eigen(tM)

temp_score = lapply(1:npc,function(ipc)
{
coef_pc= halfG_inv%*%as.matrix(eigen_res$vectors[,ipc])
fd= fd(coef=coef_pc,spline.basis)
train_score = inprod(fd,train.fd)
as.numeric(train_score)
    });

fpcs = lapply(1:npc,function(ipc)
{
coef_pc= halfG_inv%*%as.matrix(eigen_res$vectors[,ipc])
fd= fd(coef=coef_pc,spline.basis)
fd
    });

traind = do.call(cbind.data.frame,temp_score)
names(traind) = paste0("score",1:npc)
traind$theta=theta;traind$y=train_y;traind$lambda=lambda
list(traind=traind,fpcs=fpcs)
}

res= scores(theta,lambda,npc)
formula_pc = as.formula(paste0("y~",paste0(paste0("score",1:npc),collapse="+")))
glm_fit = glm(formula_pc,res$traind,family="binomial")
predy = ifelse(predict(glm_fit,newdata=res$traind,type="link")>0,1,0)
fpcs = res$fpcs

betas = glm_fit%>%coef
beta0 = as.vector(betas[1])
print(betas)
beta_func = as.vector(betas[-1])

beta_fd = beta_func[1]*fpcs[[1]]
for (i in 2:length(beta_func))
{
    beta_fd  = beta_fd + beta_func[i]*fpcs[[i]]
}


return(list(fitted=predy,beta_fd=beta_fd,beta0 = beta0,sfpcs=fpcs,theta=theta,lambda=lambda))
        

}



