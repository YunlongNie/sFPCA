# install the package 
# library("devtools")
# devtools::install_github("YunlongNie/sFPCA")

library(sFPCA)
# a continuous response case 
data(medfly)
xmat = medfly$eggcount  # 26 days and 50 flies 
y = medfly$lifetime
res = sFPCA::sfpca_con(xmat,y,0.1,1e3,npc_select=5,xmat_new = xmat[,1:10])
plot(res$beta_fd,ylab="beta(t)")
res$fitted 
plot(res$sfpcs[[1]], main="the first sFPC")



# a binary response case 
data(binary_dat);xmat = binary_dat$x;y=binary_dat$y
res = sfpcs_binary(xmat,y,npc_select=2,theta=1,xmat_new = xmat)
res$fitted
mean(y==res$fitted)
res$sfpcs
res$beta0