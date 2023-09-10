#---------------------------------------------------------------------------------------------
#Importing librariesk
#---------------------------------------------------------------------------------------------
library(gmm)

#---------------------------------------------------------------------------------------------
#Setting working directory
#---------------------------------------------------------------------------------------------
setwd("C:\\Users\\Laptop\\Downloads\\thesis\\Bivariate")




#---------------------------------------------------------------------------------------------
#Function to compute expectations of integrated variance
#---------------------------------------------------------------------------------------------

IV1_acov = function(xi_true,lambda_true,nu_true,H_true,l){
  
  
  # cat("IV acf for lag: ",l,"\n")
  # cat("Params :",xi_true,lambda_true,nu_true,H,"\n")
  acov = (xi_true^2) * integrate(function(y) {
    sapply(y, function(y) {
      (1-y)*exp((nu_true^2)/(2*(lambda_true^(2*H_true)))*(0.5*integrate(function(x) 
        exp(-abs(x))*abs(lambda_true*(l+y) + x)^(2*H_true), -Inf, Inf)$value 
        - abs(lambda_true*(l+y))^(2*H_true)))
    })
  }, 0, 1)$value 
  + (xi_true^2) * integrate(function(y) {
    sapply(y, function(y) {
      (1-y)*exp((nu_true^2)/(2*(lambda_true^(2*H_true)))*(0.5*integrate(function(x) 
        exp(-abs(x))*abs(lambda_true*abs(l-y) + x)^(2*H_true), -Inf, Inf)$value 
        - abs(lambda_true*abs(l-y))^(2*H_true)))
    })
  }, 0, 1)$value  
  
  
  #cat("acov:",acov,"\n")
  return(acov)
  
}

IV2 = function(xi_true2,lambda_true2,nu_true2,H_true2,rho_true2,H_true1,l)
{
  
  
  # cat("IV acf for lag: ",l,"\n")
  # cat("Params :",xi_true,lambda_true,nu_true,H,"\n")
  acov = (xi_true2^2) * integrate(function(y) {
    sapply(y, function(y) {
      (1-y)*exp((rho_true2*nu_true2^2)/(2*lambda_true2^(2*H_true1))*(0.5*integrate(function(x)
        exp(-abs(x))*abs(lambda_true2*(l+y) + x)^(2*H_true1), -Inf, Inf)$value
        - abs(lambda_true2*(l+y))^(2*H_true1))
        
        + (sqrt(1-rho_true2^2)*nu_true2^2)/(2*lambda_true2^(2*H_true1))*(0.5*integrate(function(x)
          exp(-abs(x))*abs(lambda_true2*(l+y) + x)^(2*H_true1), -Inf, Inf)$value
          - abs(lambda_true2*(l+y))^(2*H_true1))
        
      )
    })
  }, 0, 1)$value
  + (xi_true2^2) * integrate(function(y) {
    sapply(y, function(y) {
      (1-y)*exp((rho_true2*nu_true2^2)/(2*lambda_true2^(2*H_true2))*(0.5*integrate(function(x)
        exp(-abs(x))*abs(lambda_true2*abs(l-y) + x)^(2*H_true2), -Inf, Inf)$value
        - abs(lambda_true2*abs(l-y))^(2*H_true2))
        +
          (sqrt(1-rho_true2^2)*nu_true2^2)/(2*lambda_true2^(2*H_true2))*(0.5*integrate(function(x)
            exp(-abs(x))*abs(lambda_true2*abs(l-y) + x)^(2*H_true2), -Inf, Inf)$value
            - abs(lambda_true2*abs(l-y))^(2*H_true2)))
    })
  }, 0, 1)$value
  
  
  
}

#---------------------------------------------------------------------------------------------
#Function to compute bias correction - Approx Method(CLT)
#---------------------------------------------------------------------------------------------

bias_CLT = function(xi_true,lambda_true,nu_true,H_true,n=78){
  
  
  k0= ((nu_true^2)*gamma(1+2*H_true))/(2*lambda_true^(2*H_true))
  bias_cl = ((2*xi_true^2)*exp(k0))/n
  
  return(bias_cl)
  
}

ICoV = function(xi_true1,lambda_true1,nu_true1,H_true1,
                xi_true2,lambda_true2,nu_true2,H_true2,rho_true2,rho_true1){
  K10 = K0(nu_est = nu_true1,lambda_est = lambda_true1,H_est = H_true1)
  K20 = K0(nu_est = nu_true2,lambda_est = lambda_true2,H_est = H_true2)
  
  icov = rho_true1*sqrt(xi_true1*xi_true2)*exp(-(K10+K20+2*rho_true2*sqrt(K10*K20))/8)
  
  return(icov)
  #mean(log(as.numeric(RCoV[1,])))/(sqrt(df_bias_clt$xi1[1]*df_bias_clt$xi2[1])*exp(-(K10+K20)/8))
  
  
}




#---------------------------------------------------------------------------------------------
#Function to define the moment condition
#---------------------------------------------------------------------------------------------

g1 <- function(tet,x)
{
  
  #RV_acov = acf(RV1_mean, lag = 50, type = "covariance", plot = FALSE)
  
  m1 <- (tet[1]-x)
  
  m2 <- (IV1_acov(tet[1],tet[2],tet[3],tet[4],0)+
           bias_CLT(tet[1],tet[2],tet[3],tet[4]) - x^2)#RV_acov$acf[1])
  m3 <- (IV1_acov(tet[1],tet[2],tet[3],tet[4],l=1) - 
           x[-length(x)]*x[-1])#RV_acov$acf[2])
  m4 <- (IV1_acov(tet[1],tet[2],tet[3],tet[4],l=2) - 
           x[-1:-2]*x[-(length(x)-1):-length(x)])#RV_acov$acf[3])
  m5 <- (IV1_acov(tet[1],tet[2],tet[3],tet[4],l=3) - 
           x[-1:-3]*x[-(length(x)-2):-length(x)])#RV_acov$acf[4])
  m6 <- (IV1_acov(tet[1],tet[2],tet[3],tet[4],l=5) - 
           x[-1:-5]*x[-(length(x)-4):-length(x)])#RV_acov$acf[6])
  m7 <- (IV1_acov(tet[1],tet[2],tet[3],tet[4],l=20) - 
           x[-1:-20]*x[-(length(x)-19):-length(x)])#RV_acov$acf[21])
  m8 <- (IV1_acov(tet[1],tet[2],tet[3],tet[4],l=50) - 
           x[-1:-50]*x[-(length(x)-49):-length(x)])#RV_acov$acf[51])
  
  #print(m2)
  f <- cbind(m1,m2,m3,m4,m5,m6,m7,m8)
  return(f)
}

g2 <- function(tet,x)
{
  
  #RV_acov = acf(RV1_mean, lag = 50, type = "covariance", plot = FALSE)
  
  m1 <- (tet[1]-x[,1])
  
  m2 <- (IV2(tet[1],tet[2],tet[3],tet[4],tet[5],0.05,0)+
           bias_CLT(tet[1],tet[2],tet[3],tet[4]) - x[,1]^2)#RV_acov$acf[1])
  m3 <- (IV2(tet[1],tet[2],tet[3],tet[4],tet[5],0.05,l=1) - 
           x[,1][-length(x[,1])]*x[,1][-1])#RV_acov$acf[2])
  m4 <- (IV2(tet[1],tet[2],tet[3],tet[4],tet[5],0.05,l=2) - 
           x[,1][-1:-2]*x[,1][-(length(x[,1])-1):-length(x[,1])])#RV_acov$acf[3])
  m5 <- (IV2(tet[1],tet[2],tet[3],tet[4],tet[5],0.05,l=3) - 
           x[,1][-1:-3]*x[,1][-(length(x[,1])-2):-length(x[,1])])#RV_acov$acf[4])
  m6 <- (IV2(tet[1],tet[2],tet[3],tet[4],tet[5],0.05,l=5) - 
           x[,1][-1:-5]*x[,1][-(length(x[,1])-4):-length(x[,1])])#RV_acov$acf[6])
  m7 <- (IV2(tet[1],tet[2],tet[3],tet[4],tet[5],0.05,l=20) - 
           x[,1][-1:-20]*x[,1][-(length(x[,1])-19):-length(x[,1])])#RV_acov$acf[21])
  m8 <- (IV2(tet[1],tet[2],tet[3],tet[4],tet[5],0.05,l=50) - 
           x[,1][-1:-50]*x[,1][-(length(x[,1])-49):-length(x[,1])])#RV_acov$acf[51])
  m9 <- (ICoV(0.0225,0.005,1.25,0.05,tet[5],tet[6],tet[7],tet[8],tet[9],tet[10]) -
             x[,2])#[-1:-50]*x[,2][-(length(x[,2])-49):-length(x[,2])])#RV_acov$acf[51])
  # 
  #print(m2)
  f <- cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9 )
  return(f)
}
#---------------------------------------------------------------------------------------------
#Initial parameter values
#---------------------------------------------------------------------------------------------

#Setting xi

set_init_vals <- function(RV,t_final=1000){
  
  
  
  
  xi0 = mean(RV) 
  
  #print(xi0)
  
  #Application of the scaling law - Setting H and nu
  
  q = 2
  m = 6
  h = c(1,2,3,4,5,6)  
  Kq = 2^(q/2)*gamma((q+1)/2)/sqrt(pi)
  
  
  gamma1 =  sum(abs(diff(log(RV),lag = 1))[1:(t_final-m)]^q)/(t_final-m)
  gamma2 =  sum(abs(diff(log(RV),lag = 2))[1:(t_final-m)]^q)/(t_final-m)
  gamma3 =  sum(abs(diff(log(RV),lag = 3))[1:(t_final-m)]^q)/(t_final-m)
  gamma4 =  sum(abs(diff(log(RV),lag = 4))[1:(t_final-m)]^q)/(t_final-m)
  gamma5 =  sum(abs(diff(log(RV),lag = 5))[1:(t_final-m)]^q)/(t_final-m)
  gamma6 =  sum(abs(diff(log(RV),lag = 6))[1:(t_final-m)]^q)/(t_final-m)
  
  gamma_vec = c(gamma1,gamma2,gamma3,gamma4,gamma5,gamma6)
  
  model = lm(log(gamma_vec) ~ log(abs(h)))
  
  #write.csv(summary(model),file =  "lm_summary.csv")
  #print(summary(model))
  coeffs = unname(model$coefficients)
  
  nu0 = exp((coeffs[1]-log(Kq))/q)
  H0 = coeffs[2]/q
  
  #Setting lambda0
  
  sample_var_RV = var(log(RV))
  lambda0 = (((nu0^2)*gamma(1+2*H0))/(2*sample_var_RV))^(1/(2*H0))
  
  return(c(xi0,lambda0,nu0,H0))
}

#cat("Intial xi, lambda, H and nu:", xi0,lambda0,H0,nu0)



K0 <- function(nu_est,lambda_est,H_est){
  
  
  k = ((nu_est^2)/(2*lambda_est)^(2*H_est))*gamma(1+2*H_est)
  return(k)
  
}

# ICoV = function(xi_true1,lambda_true1,nu_true1,H_true1,
#                 xi_true2,lambda_true2,nu_true2,H_true2,rho_true2,rho_true1){
#   K10 = K0(nu_est = nu_true1,lambda_est = lambda_true1,H_est = H_true1)
#   K20 = K0(nu_est = nu_true2,lambda_est = lambda_true2,H_est = H_true2)
#   
#   icov = rho_true1*sqrt(xi_true1*xi_true2)*exp(-(K10+K20+2*rho_true2*sqrt(K10*K20))/8)
#   
#   return(icov)
#   #mean(log(as.numeric(RCoV[1,])))/(sqrt(df_bias_clt$xi1[1]*df_bias_clt$xi2[1])*exp(-(K10+K20)/8))
#   
#   
# }
#---------------------------------------------------------------------------------------------
#GMM estimation step
#---------------------------------------------------------------------------------------------


GMM05_estimation <-function(RV){
  
  initial_vals = set_init_vals(RV = RV)
  
  
  res_bias_clt <- gmm(g1,RV,c(xi_1 = initial_vals[1],lambda_1 = initial_vals[2],
                              nu_1=initial_vals[3],H_1 =initial_vals[4]), optfct = 'nlminb', 
                      type = "iterative", kernel = "Parzen",approx = "ARMA(1,1)",
                      lower =c(0.01,0.001,1.22,0.01), upper=c(0.03,0.009,1.4,0.07),itermax = 1000)
  
  summ_res_bias_clt = summary(res_bias_clt )
  
  #print(summ_res_bias_clt)
  return(c(initial_vals,summ_res_bias_clt$coefficients[,1],
           as.list(summ_res_bias_clt$stest[2]$test)[[1]]))
}

GMM3_estimation <-function(RV){
  
  initial_vals = set_init_vals(RV = RV)
  
  
  res_bias_clt <- gmm(g1,RV,c(xi_1 = initial_vals[1],lambda_1 = initial_vals[2],
                              nu_1=initial_vals[3],H_1 =initial_vals[4],rho2=0,rho1=0), optfct = 'nlminb', 
                      type = "iterative", kernel = "Parzen",approx = "ARMA(1,1)",
                      lower =c(0.02,0.009,0.5,0.25), upper=c(0.03,0.02,0.8,0.35),itermax = 1000)
  
  summ_res_bias_clt = summary(res_bias_clt )
  
  #print(summ_res_bias_clt)
  return(c(initial_vals,summ_res_bias_clt$coefficients[,1],
           as.list(summ_res_bias_clt$stest[2]$test)[[1]]))
}




#---------------------------------------------------------------------------------------------
#Loading data, setting simulation parameters and running GMM
#---------------------------------------------------------------------------------------------


RV1= read.csv('RV1_rhox0.5_rhoy_0.0.csv',nrows=1) 
RV2= read.csv('RV2_rhox0.5_rhoy_0.0.csv',nrows=1) 
RCoV=read.csv('RCoV_rhox0.5_rhoy_0.0.csv',nrows=1) 

t_final = 1000        #250#4000#4000#500
N = 390              #23400

delta_t = 1/N      #Timestep of integration #0.001#1/23400#0.001
sim = 10000        #Number of simulations     

# lower1 =c(0.01,0.001,1.22,0.01)
# upper1 =c(0.03,0.009,1.4,0.07)
# 
# lower2 =c(0.01,0.009,0.7,0.07)
# upper2 = c(0.03,0.02,0.8,0.15)

start_time <- Sys.time()

df1_bias_clt = as.data.frame(t(apply(RV1, 1, function(x) GMM05_estimation(RV = as.numeric(x)))))
colnames(df1_bias_clt)=c("xi01","lamda01","nu01","H01","xi1", "lambda1", "nu1","H1","J1")



df2_bias_clt = as.data.frame(t(apply(RV2, 1, function(x) GMM3_estimation(RV = as.numeric(x)))))
colnames(df2_bias_clt)=c("xi02","lamda02","nu02","H02","xi2", "lambda2", "nu2","H2","J2")

# initial_vals= set_init_vals(RV =as.numeric(RV2))
# res1 = gmm(g1,RV2,c(xi_1 = initial_vals[1],lambda_1 = initial_vals[2],
#                    nu_1=initial_vals[3],H_1 =initial_vals[4]), optfct = 'nlminb', 
#            type = "iterative", kernel = "Parzen",approx = "ARMA(1,1)",
#            lower =c(0.02,0.009,0.7,0.07,0.07), upper=c(0.03,0.02,0.8,0.15),itermax = 1000)


end_time <- Sys.time()
cat("Execution time:", end_time - start_time)

df_bias_clt = merge(df1_bias_clt,df2_bias_clt)

df_bias_clt[nrow(df_bias_clt) + 1,] = colMeans(df_bias_clt)
df_bias_clt[nrow(df_bias_clt) + 1,] = sapply(df_bias_clt,sd)
write.csv(df_bias_clt,file='GMM_RV1_rx7_ry0_clt_1.csv', row.names=FALSE)


K10 = K0(nu_est = df_bias_clt$nu1[1],lambda_est = df_bias_clt$lambda1[1],H_est = df_bias_clt$H1[1])
K20 = K0(nu_est = df_bias_clt$nu2[1],lambda_est = df_bias_clt$lambda2[1],H_est = df_bias_clt$H2[1])
# 
#mean(as.numeric(log(RCoV[1,]))/(sqrt(df_bias_clt$xi1[1]*df_bias_clt$xi2[1])*exp(-(K10+K20)/8)))
      
     