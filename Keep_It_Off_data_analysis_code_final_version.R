#######################
require(readxl)
mydata = read_excel("/secure/project/Volpp_W2H_WLossR01/Yong_Chen/weigh-in data.xlsx")
dim(mydata)
# [1] 68796    27
num_p <- seq(1:length(unique(mydata$Participant)))
length(num_p)
# 189

mydata = transform(mydata,id=as.numeric(factor(mydata$Participant)))
# mydata = mydata[which(mydata$id%%5==0),]
mydata_6 = mydata[which(mydata$month < 7),]
##########################



##########################
baseline <- read_excel("/secure/project/Volpp_W2H_WLossR01/Yong_Chen/baseline data.xlsx")
up_baseline <- baseline[baseline$Participant %in% mydata$Participant,]
dim(up_baseline)
#[1] 189 112
##########################



###############################
##  R CMD SHLIB 
dyn.load("~/project/keep_it_off/mylikC.so")
dyn.load("~/project/keep_it_off/covC.so")
library(mvtnorm)
library(geepack)
set.seed(5)
###############################



###################### functions ###################################
mylik <- function(mypar, mydat)
{
  res <- .C("mylikC", as.integer(mydat$p), as.integer(mydat$n.t), as.double(mydat$yvec), as.integer(mydat$index), as.integer(mydat$count),
            as.double(mydat$xall), as.integer(mydat$neach), as.double(mypar), as.double(mydat$sigma2), result=double(1))
  return(-res[["result"]])
}

covcal <- function(mydat)
{
  matdim <- mydat$p*mydat$p
  res <- .C("covcalC", as.integer(mydat$p), as.integer(mydat$n.t), as.double(mydat$yvec), as.integer(mydat$index), as.integer(mydat$count),
            as.double(mydat$xall), as.integer(mydat$neach), as.double(mydat$sigma2), as.double(mydat$beta.t), mat1=double(matdim), mat2=double(matdim), phi=double(mydat$p))
  mat1 <- matrix(res[["mat1"]], p,p, byrow=T)
  mat2 <- matrix(res[["mat2"]], p,p, byrow=T)
  phi <- res[["phi"]]
  print(phi)
  print(mat1)
  print(mat2)
  print(res[[1]])
  matcov <- solve(mat1)%*%mat2%*%solve(mat1)
  print(matcov)
  return(matcov)
}
####################################################################



###############################
mydata <- mydata_6
# mydata <- mydata[which(mydata$Arm != "A"),]  # B&C = 112
# up_baseline <- up_baseline[which(up_baseline$Arm != "L"),]
# mydata <- mydata[which(mydata$Arm != "B"),] #A&C = 115
# up_baseline <- up_baseline[which(up_baseline$Arm != "J"),]
# arm assignment with baseline
# C: K (38) control group
# B: J (74) direct payment group
# A: L (77) lottery group
###############################




###############################
##tocount is used to record how many points each subject has
tocount = numeric(0)
for (i in unique(mydata$id)){
  subset = mydata[which(mydata$id==i),]
  tocount_temp = rep(sum((subset$weight_daily != "."),na.rm=TRUE),183)
  tocount = c(tocount,tocount_temp)
}
mydata$tocount = tocount
###############################


###############################
# exclude the rows with NA in weight
mydata_wo_NA <- mydata[which(mydata$weight_daily != "."),]
wo_NA_id <- which(mydata$weight_daily != ".")
###############################



###############################
## count is used to record the order of the point
l = 1
for (i in unique(mydata_wo_NA$id)){
  for (j in 1:sum(mydata_wo_NA$id==i)){
    mydata_wo_NA$count[l] = j
    l = l+1
  }
}
###############################

######################################################################
################# ======  analysis ================
######################################################################
###############################
## index is used to record the starting point of each subject if recorded all outcomes in vector
mydata_wo_NA$row_id <- c(1:dim(mydata_wo_NA)[1])
index <- mydata_wo_NA[which(mydata_wo_NA$count == 1),]$row_id-1
print(paste("The length of index is", length(index)))
###############################

###############################
##yvec is a vector which records all outcomes 
yvec <- as.numeric(mydata_wo_NA$Wchange_daily)
# yvec <- as.numeric(mydata_wo_NA$Wchange_daily)
print(paste("The length of yvec is", length(yvec)))
###############################

###############################
######## predictors ###########
# gender
x1vec_t <- rep(up_baseline$gender,each=183)
x1vec <- x1vec_t[wo_NA_id]
print(paste("The length of x1vec is", length(x1vec)))

# BMI at study start
x2vec_t <- mydata$BMI_study_start[wo_NA_id]
x2vec <- x2vec_t-mean(x2vec_t)
print(paste("The length of x2vec is", length(x2vec)))

# age
x3vec_t <- rep(up_baseline$age,each=183)[wo_NA_id]
x3vec <- x3vec_t-mean(x3vec_t)
print(paste("The length of x3vec is", length(x3vec)))

# time
# mydata$day_W <-(mydata$day - min(mydata$day))/7
x5vec_t <- mydata$day[wo_NA_id]
x5vec <- x5vec_t/max(x5vec_t)
print(paste("The length of x5vec is", length(x5vec)))


# indicator * time
# # I(group)
x4vec_1 <- rep(0,length(x1vec))
for (i in 1:length(x1vec)){
  if (mydata_wo_NA$Arm[i] == "A"){
    x4vec_1[i] = 1
  }
}
x4vec_1 <- x4vec_1 * x5vec
print(paste("The length of x4vec_1 (77 participants) is", length(x4vec_1)))


x4vec_2 <- rep(0,length(x1vec))
for (i in 1:length(x1vec)){
  if (mydata_wo_NA$Arm[i] == "B"){
    x4vec_2[i] = 1
  }
}
x4vec_2 <- x4vec_2 * x5vec
print(paste("The length of x4vec_2 (74 participants) is", length(x4vec_2)))
###############################


###############################
# simple linear regression
linear_fit <- lm(yvec ~ 0+x1vec+x2vec+x3vec+x5vec+x4vec_1+x4vec_2)
print(summary(linear_fit))
###############################

########## Analysis with control vs Arm 1 and 2 #############
idvec <- mydata_wo_NA$Participant
neach <- length(x1vec)
p <- 6
ini.val <- rep(0, p)

xall <- c(x1vec, x2vec, x3vec, x5vec,x4vec_1, x4vec_2)
par.geeglm <- geeglm(yvec~0+x1vec+x2vec+x3vec+x5vec+x4vec_1+x4vec_2, id=idvec, family=gaussian, corstr = "independence")
totalsigma2 <- var(par.geeglm$res[mydata_wo_NA$count==1])
print(par.geeglm)
print(summary(par.geeglm))

count = mydata_wo_NA$tocount[mydata_wo_NA$count==mydata_wo_NA$tocount]
res1<- optim(ini.val, mylik, method = "BFGS",control = list(trace=TRUE,REPORT=1),
             hessian=T, mydat=list(p=p, n.t=length(unique(mydata_wo_NA$id)), yvec=yvec, index=index,
                                   count=count, xall=xall, neach=neach,
                                   sigma2=totalsigma2))
# point estimate
par.est <- res1$par
print(par.est)

# se
n.t=length(unique(mydata_wo_NA$id))
mydat=list(p=p, n.t=n.t, yvec=yvec, index=index, count=count, 
           xall=xall, neach=neach, sigma2=totalsigma2, beta.t=res1$par)		
res2 <- covcal(mydat)/n.t
par.ase <- sqrt(res2[col(res2)==row(res2)]) 




