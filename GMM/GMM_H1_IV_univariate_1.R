#---------------------------------------------------------------------------------------------
#Importing librariesk
#---------------------------------------------------------------------------------------------
library(gmm)

#---------------------------------------------------------------------------------------------
#Setting working directory
#---------------------------------------------------------------------------------------------
#setwd("C:\\Users\\Laptop\\Downloads\\thesis\\Univariate\\Simulation")




#---------------------------------------------------------------------------------------------
#Function to compute expectations of integrated variance
#---------------------------------------------------------------------------------------------

IV_acov = function(xi_true,lambda_true,nu_true,H_true,l){
  
  
  # cat("IV acf for lag: ",l,"\n")
  # cat("Params :",xi_true,lambda_true,nu_true,H,"\n")
  acov = (xi_true^2) * integrate(function(y) {
    sapply(y, function(y) {
      (1-y)*exp((nu_true^2)/(2*lambda_true^(2*H_true))*(0.5*integrate(function(x) 
        exp(-abs(x))*abs(lambda_true*(l+y) + x)^(2*H_true), -Inf, Inf)$value 
        - abs(lambda_true*(l+y))^(2*H_true)))
    })
  }, 0, 1)$value 
  + (xi_true^2) * integrate(function(y) {
    sapply(y, function(y) {
      (1-y)*exp((nu_true^2)/(2*lambda_true^(2*H_true))*(0.5*integrate(function(x) 
        exp(-abs(x))*abs(lambda_true*abs(l-y) + x)^(2*H_true), -Inf, Inf)$value 
        - abs(lambda_true*abs(l-y))^(2*H_true)))
    })
  }, 0, 1)$value  
  
  
  #cat("acov:",acov,"\n")
  return(acov)
  
}



#---------------------------------------------------------------------------------------------
#Function to define the moment condition
#---------------------------------------------------------------------------------------------

g1 <- function(tet,x)
{
  
  #RV_acov = acf(RV1_mean, lag = 50, type = "covariance", plot = FALSE)
  
  m1 <- (tet[1]-x)
  
  m2 <- (IV_acov(tet[1],tet[2],tet[3],tet[4],0)- x^2)#RV_acov$acf[1])
  m3 <- (IV_acov(tet[1],tet[2],tet[3],tet[4],l=1) - 
           x[-length(x)]*x[-1])#RV_acov$acf[2])
  m4 <- (IV_acov(tet[1],tet[2],tet[3],tet[4],l=2) - 
           x[-1:-2]*x[-(length(x)-1):-length(x)])#RV_acov$acf[3])
  m5 <- (IV_acov(tet[1],tet[2],tet[3],tet[4],l=3) - 
           x[-1:-3]*x[-(length(x)-2):-length(x)])#RV_acov$acf[4])
  m6 <- (IV_acov(tet[1],tet[2],tet[3],tet[4],l=5) - 
           x[-1:-5]*x[-(length(x)-4):-length(x)])#RV_acov$acf[6])
  m7 <- (IV_acov(tet[1],tet[2],tet[3],tet[4],l=20) - 
           x[-1:-20]*x[-(length(x)-19):-length(x)])#RV_acov$acf[21])
  m8 <- (IV_acov(tet[1],tet[2],tet[3],tet[4],l=50) - 
           x[-1:-50]*x[-(length(x)-49):-length(x)])#RV_acov$acf[51])
  
  #print(m2)
  f <- cbind(m1,m2,m3,m4,m5,m6,m7,m8)
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


#---------------------------------------------------------------------------------------------
#GMM estimation step
#---------------------------------------------------------------------------------------------


GMM_estimation <-function(RV){
  
  initial_vals = set_init_vals(RV = RV)
  
  res <- gmm(g1,RV,c(xi_1 = initial_vals[1],lambda_1 = initial_vals[2],
                                nu_1=initial_vals[3],H_1 =initial_vals[4]), optfct = 'nlminb', 
                        type = "iterative", kernel = "Parzen",approx = "ARMA(1,1)",
                        lower =c(0.0222,0.011,0.76,0.09), upper=c(0.03,0.02,0.8,0.15),itermax = 1000)
  
  summ_res = summary(res )
  
  #print(summ_res)
  return(c(initial_vals,summ_res$coefficients[,1],
           as.list(summ_res$stest[2]$test)[[1]]))
}


#---------------------------------------------------------------------------------------------
#Loading data, setting simulation parameters and running GMM
#---------------------------------------------------------------------------------------------


RV_daily = read.csv('IV1_uni_1.csv',nrows=100) 

t_final = 1000        #250#4000#4000#500
N = 390              #23400

delta_t = 1/N      #Timestep of integration #0.001#1/23400#0.001
#sim = 10000        #Number of simulations     


start_time <- Sys.time()

df = as.data.frame(t(apply(RV_daily, 1, function(x) GMM_estimation(RV = as.numeric(x)) )))
colnames(df)=c("xi0","lamda0","nu0","H0","xi", "lambda", "nu","H","J")


end_time <- Sys.time()
cat("Execution time:", end_time - start_time)


df[nrow(df) + 1,] = colMeans(df)
df[nrow(df) + 1,] = sapply(df,sd)
write.csv(df,file='GMM_H1_univariate_IV_1.csv', row.names=FALSE)


