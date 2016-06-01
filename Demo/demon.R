# install the package 
# library("devtools")
# install_github("YunlongNie/sFPCA")

library(sFPCA)
# a continuous response case 
data(medfly)
xmat = medfly$eggcount  # 26 days and 50 flies 
y = medfly$lifetime
res = sfpca(xmat,y,0.1,1e3,npc_select=5)
plot(res$beta_fd,ylab="beta(t)")
res$fitted 
res$sfpcs



# a binary response case 
data(binary_dat);xmat = binary_dat$x;y=binary_dat$y
res = sfpcs_binary(xmat,y,npc=2,theta=1)
res$fitted
mean(y==res$fitted)
res$sfpcs
res$beta0