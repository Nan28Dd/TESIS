library(dplyr)
library(tidyverse)
library(reshape)
library(ggplot2)
library(purrr)
library(car)
library(tictoc)


#Datos
#Los datos y las Hipótesis se usaran tanto para la Simulación 
#del Fitness como para la Simulación de las Probabilidades de
#Fijación

experimentos1=200
N1=100 #Tamaño de la muestra de la población a ocupar cada día
gamma1=2  #Parámetro para obtener la proporción de la población a rebasar al final de cada día
K_01=1    #Tamaño de la población del tipo mutante
No_K01=length(K_01)+1:N1   #Tamaño de la población del tipo no mutante
r0_1=1.05   #Tasa de reproducción para los no mutantes
t_max_1=4
t_aum_1=t_max_1+1
dia_adrian=10000

# Parámaetros a ajustar para obtener las mejores aproximaciones
set.seed(500)
b1= runif(1,min=0,max=1/2)
a1= 4*b1#(a>3*b)
a11=runif(1,min=b1,max=a1)#(a>b)
rho_1=N1^{-b1}
mu_1=N1^{-a1}
mu_11=N1^{-a11}

b2= runif(1,min=0,max=2/3)
a2= 3*b2
rho_2=N1^{-b2}
mu_2=N1^{-a2}

b3= runif(1,min=0,max=1)
a3= 3*b3

rho_3=N1^{-b3}
mu_3=N1^{-a3}

ep1= runif(1,min=0,max=1/2)
d1= runif(1,min=0.75)


#Reescalamiento del Fitness 1
s=rho_1^{-2}*mu_1^{-1}
dias=dia_adrian/s
y=dias*s
eje_x=seq(from=0, to=dias, by=1/(dias*s))# Para Fitnes f
eje_x1=seq(from=0, to=s*dias-1, by=1)# Para Fitness Adrian


t1=1
eje_y=((r0_1+rho_1)*t1)/(r0_1*t1)
eje_y3=((r0_1+rho_2)*t1)/(r0_1*t1)
eje_y4=((r0_1+rho_3)*t1)/(r0_1*t1)

#OTROS ARGUEMENTOS IMPORTANTES
C_gamma=(gamma1*log(gamma1))/(gamma1-1)
#Valor al que se debe de aproximar las probabilidades de fijación de acuerdo con el Teorema 4
P_N1=rho_1*(C_gamma/(r0_1))
P_N2=rho_2*(C_gamma/(r0_1))
P_N3=rho_3*(C_gamma/(r0_1))

t_f1=(7/8)^{rho_1^{-d1}}
t_f2=(7/8)^{rho_2^{-d1}}
t_f3=(7/8)^{rho_3^{-d1}}

const1=rho_1^{-1-3*d1}
const2=rho_2^{-1-3*d1}
const3=rho_3^{-1-3*d1}

##################################
#Funciones Generales

#Función que mueve de conjunto al elemento seleccionado para adquirir una mutación
#(es decir, que se vea afectado por la ventaja selectiva y que ahora se reproduzca más rápido)
#eliminandolo del conjunto de No Mutantes y agregandolo al conjunto de Mutantes.
Mut=function(n_mut_i,K,No_K){
  no_k_i=which(No_K==n_mut_i)
  K=append(K,n_mut_i)
  No_K=No_K[-no_k_i]
  Y=list("K"=K,"NoK"=No_K)
  return(Y)
}

#Función que asigna la tasa  de reproducción a cada individuo dependiendo del 
#conjunto en donde se encuentre (Mutantes o No Mutantes)
R_ind=function(r0,rho,K,No_K){
  R_m<-K
  R_nm<-No_K
  R_j<-c()
  for(j in R_m){
    R_m<-rep(r0 + rho, times=length(K))
  }  
  for (j in R_nm){
    R_nm <-rep(r0, times=length(No_K))
  }
  R_j=c(R_m,R_nm)
  return(R_j)
}

#Función del desarrollo de un día en el experimento de Lenski 
#para el tiempo de paro 1.
Lenski_dia_i1 <- function(r0,rho,mu,gamma,N,K,No_K,Ind){
  #inicio del día
  #Ind=1 pernetece a la parte de la función para el Estudio del Fitness
  #Ind =0 pertenece a la parte de la función para el estudio de las probabilidades de Fijación
  if(length(K)>0){
    K=1:length(K)
  }
  if(length(No_K)>0){
    No_K=(length(K)+1):N
  }#aquí se determina si habría o no una mutación beneficiosa
  if(length(No_K)>0){
    aux=rbinom(1,1,mu)
    if(aux==1){#si la hay, se muestrea sobre la población no mutante y el individuo se agrega a la población Mutante
      n_mut_i=sample(No_K,1,replace=F)
      Y_i=Mut(n_mut_i,K,No_K)
    }else{
      K=K
      No_K=No_K
      Y_i=list("K"=K,"NoK"=No_K)
    }
  }else{
    K=K
    No_K=No_K
    Y_i=list("K"=K,"NoK"=No_K)
  }
  Ki=Y_i$K
  No_Ki=Y_i$NoK
  Y=unlist(Y_i)
  r<-R_ind(r0,rho,Ki,No_Ki)
  T1<-0
  if(Ind==1){
    while((length(Y)<gamma*N)&(length(Ki)<N)){
      tiempo<-rep(0,times=length(Y))
      for(l in 1:length(Y)){#se guardan las exponenciales, que representan los tiempos en el que van "naciendo" los individuos
        tiempo[l] <-rexp(n=1,rate=r[l])
      }#nos fijamos en el tiempo mínimo y en que tipo de individuo es y en base a eso se va agregando al conjunto al que pertenece
      mi=min(tiempo)
      T1=T1+mi
      t=which(tiempo==mi)[1]
      for (v in No_Ki) {
        if(v==t){
          No_Ki=c(No_Ki,length(Y)+1)
        }else{
          No_Ki=No_Ki
        }
      }
      for (c in Ki) {
        if(c==t){
          Ki=c(Ki,length(Y)+1)
        }else{
          Ki=Ki  
        }
      }
      r<-R_ind(r0,rho,Ki,No_Ki)
      Y=c(Ki,No_Ki)
    }
  }else{
    while((length(Y)<gamma*N)&((length(Ki)<N) |(length(No_Ki)<N))){
      tiempo<-rep(0,times=length(Y))
      for(l in 1:length(Y)){#se guardan las exponenciales, que representan los tiempos en el que van "naciendo" los individuos
        tiempo[l] <-rexp(n=1,rate=r[l])
      }
      mi=min(tiempo)
      T1=T1+mi#nos fijamos en el tiempo mínimo y en que tipo de individuo es y en base a eso se va agregando al conjunto al que pertenece
      t=which(tiempo==mi)[1]
      for (v in No_Ki) {
        if(v==t){
          No_Ki=c(No_Ki,length(Y)+1)
        }else{
          No_Ki=No_Ki
        }
      }
      for (c in Ki) {
        if(c==t){
          Ki=c(Ki,length(Y)+1)
        }else{
          Ki=Ki  
        }
      }
      r<-R_ind(r0,rho,Ki,No_Ki)
      Y=c(Ki,No_Ki)
    }
  }
  #Aqui es cuando ya supero la propporción dada por gamma*N y se hace el 
  #muestreo de los N individuos con los que se empezará el día siguiente.
    if(length(No_Ki)>0){
      K<-NULL
      No_K<-NULL
      nombres1<-rep(1,times=length(Ki))
      nombres2<-rep(0,times=length(No_Ki))
      Ytipo<-c(nombres1,nombres2)
      Pob<-data.frame("Y"=Y,"Tipo"=Ytipo)
      Y_n<-sample(1:nrow(Pob), N,replace = FALSE)
      Pob[Y_n,]
      for(z in Y_n){
        if(Pob$Tipo[z]==1){
          K<-c(K,Pob$Y[z])
        }else{
          No_K<-c(No_K,Pob$Y[z])
        }
      }
    }else{
      K=Ki
      No_K=No_Ki
    }
  POB_dia<-list("Mut"=K, "No_Mut"=No_K, "Num_Mutantes"=length(K), "Tiempos"=T1)
}

#Función del desarrollo de un día en el experimento de Lenski 
#para el tiempo de paro 2.
Lenski_dia_i2 <- function(r0,rho,mu,gamma,N,K,No_K,Ind){
  #inicio del día
  if(length(K)>0){
    K=1:length(K)
  }
  if(length(No_K)>0){
    No_K=(length(K)+1):N
  }#aquí se determina si habrá o no una mutación beneficiosa
  if(length(No_K)>0){
    aux=rbinom(1,1,mu)
    if(aux==1){#si la hay, se muestrea sobre la población no mutante y el individuo se agrega a la población Mutante
      n_mut_i=sample(No_K,1,replace=F)
      Y_i=Mut(n_mut_i,K,No_K)
    }else{
      K=K
      No_K=No_K
      Y_i=list("K"=K,"NoK"=No_K)
    }
  }else{
    K=K
    No_K=No_K
    Y_i=list("K"=K,"NoK"=No_K)
  }
  Y_i
  Ki=Y_i$K
  No_Ki=Y_i$NoK
  r<-R_ind(r0,rho,Ki,No_Ki)
  Y=unlist(Y_i)
  T2<-0
  if(Ind==1){
    while((T2<(log(gamma)/r0))&(length(Ki)<N)){
      tiempo2<-rep(0,times=length(Y))
      for(m in 1:length(Y)){#se guardan las exponenciales, que representan los tiempos en el que van "naciendo" los individuos
        tiempo2[m]=rexp(1,r[m])
      }#nos fiajamos en el tiempo minimo y en que tipo de individuo es y en base a eso se va agragando al conjunto al que pertenece
      tiempo2
      mi=min(tiempo2)
      t=which(tiempo2==mi)[1]
      for (v1 in No_Ki) {
        if(v1==t){
          No_Ki=c(No_Ki,length(Y)+1)
        }else{
          No_Ki=No_Ki
        }
      }
      for (v1 in Ki) {
        if(v1==t){
          Ki=c(Ki,length(Y)+1)
        }else{
          Ki=Ki
        }
      }
      r<-R_ind(r0,rho,Ki,No_Ki)
      Y=c(Ki,No_Ki)
      T2=T2+mi
    }
  }else{
    while((T2<(log(gamma)/r0))&((length(Ki)<N) |(length(No_Ki)<N))){
      tiempo2<-rep(0,times=length(Y))
      for(m in 1:length(Y)){#se guardan las exponenciales, que representan los tiempos en el que van "naciendo" los individuos
        tiempo2[m]=rexp(1,r[m])
      }
      tiempo2#nos fiajamos en el tiempo minimo y en que tipo de individuo es y en base a eso se va agragando al conjunto al que pertenece
      mi=min(tiempo2)
      t=which(tiempo2==mi)[1]
      for (v1 in No_Ki) {
        if(v1==t){
          No_Ki=c(No_Ki,length(Y)+1)
        }else{
          No_Ki=No_Ki
        }
      }
      for (v1 in Ki) {
        if(v1==t){
          Ki=c(Ki,length(Y)+1)
        }else{
          Ki=Ki
        }
      }
      r<-R_ind(r0,rho,Ki,No_Ki)
      Y=c(Ki,No_Ki)
      T2=T2+mi
    }
  }
  #Aqui es cuando ya supero la proporción dada por gamma*N y se hace el 
  #muestreo de los N individuos con los que se empezará el día siguiente.
    if(length(No_Ki)>0){
      K<-NULL
      No_K<-NULL
      nombres1<-rep(1,times=length(Ki))
      nombres2<-rep(0,times=length(No_Ki))
      Ytipo<-c(nombres1,nombres2)
      Pob<-data.frame("Y"=Y,"Tipo"=Ytipo)
      Y_n<-sample(1:nrow(Pob), N,replace = FALSE)
      Pob[Y_n,]
      for(z in Y_n){
        if(Pob$Tipo[z]==1){
          K<-c(K,Pob$Y[z])
        }else{
          No_K<-c(No_K,Pob$Y[z])
        }
      }
    }else{
      K=Ki
      No_K=No_Ki
    }
  POB_dia<-list("Mut"=K, "No_Mut"=No_K,"Num_Mutantes"=length(K),"Tiempos"=T2)
}

#Función de la repetición de varios días del experimento. En donde se obtienen
#el número de elemntos mutantes y no mutantes cada d\'ia observado (para el tiempo de paro 1).
Lenski_Exp1=function(r0,rho,mu,gamma,N,K,No_K,t_aum,Ind){
  f<-list()
  f[[1]]=list("Mut"=K,"No_Mut"=No_K,"Num_Mutantes"=length(K),"Tiempos"=0) # Corresponde  al dia cero
  for(i in 2:t_aum){
    tic()
    f[[i]]=Lenski_dia_i1(r0,rho,mu,gamma,N,f[[i-1]]$Mut,f[[i-1]]$No_Mut,Ind)
    exectime <- toc()
    exectime <- exectime$toc - exectime$tic
    print(exectime)
  }
  
  Mut<-list()
  for(q in 1:t_aum){
    if(is.null(f[[q]]$Mut)){
      Mut[[q]]=0
    }else{
      Mut[[q]]=f[[q]]$Mut
    }
  }
  
  No_Mut<-list()
  for(p in 1:t_aum){
    if(is.null(f[[p]]$No_Mut)){
      No_Mut[[p]]=0
    }else{
      No_Mut[[p]]=f[[p]]$No_Mut
    }
  }
  
  N_Mut<-c()
  for(q in 1:t_aum){
    N_Mut[q]=f[[q]]$Num_Mutantes
  }
  
  Tiempos<-c()
  for(q in 1:t_aum){
    Tiempos[q]=f[[q]]$Tiempos
  }
  
  Tasas<-vector("list",length(t_aum))
  for (x in 1:t_aum) {
    Tasas[[x]]=R_ind(r0,rho,f[[x]]$Mut,f[[x]]$No_Mut)
  }
  
  POB<-list("Dias"=1:t_aum,"Mutantes"=Mut, "NoMutantes"=No_Mut, "Num_Mut"=N_Mut,"Tasas"=Tasas, "Tiempos"=Tiempos)
  POB
}

#Función de la repetición de varios días del experimento. En donde se obtienen
#el numero de elemntos mutantes y no mutantes cada día observado (para el tiempo de paro 2).
Lenski_Exp2=function(r0,rho,mu,gamma,N,K,No_K,t_aum,Ind){
  f<-list()
  f[[1]]=list("Mut"=K,"No_Mut"=No_K,"Num_Mutantes"=length(K),"Tiempo"=0)
  for(i in 2:t_aum){
    tic()
    f[[i]]=Lenski_dia_i2(r0,rho,mu,gamma,N,f[[i-1]]$Mut,f[[i-1]]$No_Mut,Ind)
    exectime <- toc()
    exectime <- exectime$toc - exectime$tic
    print(exectime)
  }
  Mut<-list()
  for(q in 1:t_aum){
    if(is.null(f[[q]]$Mut)){
      Mut[[q]]=0
    }else{
      Mut[[q]]=f[[q]]$Mut
    }
  }
  
  No_Mut<-list()
  for(p in 1:t_aum){
    if(is.null(f[[p]]$No_Mut)){
      No_Mut[[p]]=0
    }else{
      No_Mut[[p]]=f[[p]]$No_Mut
    }
  }
  
  N_Mut<-c()
  for(q in 1:t_aum){
    N_Mut[q]=f[[q]]$Num_Mutantes
  }
  Tiempos<-c()
  for(q in 1:t_aum){
    Tiempos[q]=f[[q]]$Tiempos
  }
  
  Tasas<-vector("list",length(t_aum))
  for (x in 1:t_aum) {
    Tasas[[x]]=R_ind(r0,rho,f[[x]]$Mut,f[[x]]$No_Mut)
  }
  
  POB<-list("Dias"=1:t_aum,"Mutantes"=Mut, "NoMutantes"=No_Mut,  "Num_Mut"=N_Mut, "Tasas"=Tasas,"Tiempos"=Tiempos)
  POB
}



#################################
#Funciones para la Simulación del Fitness

#Función de l aproximación de la curva logarítmica f(t)
Fitness1 = function(gamma,r0,rho,mu,t){
  s=rho^{-2}*mu^{-1}
  dias=t/s
  #i=seq(from=0, to=t, by=dias)
  i=1:t
  aux_t=(2*gamma*log(gamma))/((gamma-1)*r0^{2})
  ft=sqrt(1+(aux_t*i))
}

#Función que represnta a la aprimación del Fitness de Lenski Fi
#con el rescalando
Fitness_Adrian = function(N,lenski,t_aum,r0,rho,mu,gamma){
  t=lenski$Tiempos
  r<-array(unlist(lenski$Tasas),dim=c(N,length(lenski$Tasas)))
  denominador<-r0*t
  numerador<-c()
  fit<-c()
  for(i in 1:t_aum){
    numerador[i]<-log((1/N)*sum(exp(t[i]*r[,i])))
    fit[i]=numerador[i]/denominador[i]
  }
  fit
}

##################################
#SIMULACIÓN DE LA APROXIMACIÓN DEL FITNESS
#Tiempo de paro 1
Lenski_1=Lenski_Exp1(r0_1,rho_1,mu_1,gamma1,N1,K_01,No_K01,dia_adrian,1)
Lenski_2=Lenski_Exp1(r0_1,rho_1,mu_11,gamma1,N1,K_01,No_K01,dia_adrian,1)
Lenski_3=Lenski_Exp1(r0_1,rho_2,mu_2,gamma1,N1,K_01,No_K01,dia_adrian,1)
Lenski_4=Lenski_Exp1(r0_1,rho_3,mu_3,gamma1,N1,K_01,No_K01,dia_adrian,1)

Fit1=Fitness_Adrian(N1,Lenski_1,dia_adrian,r0_1,rho_1,mu_1,gamma1)
Fit2=Fitness_Adrian(N1,Lenski_2,dia_adrian,r0_1,rho_1,mu_11,gamma1)
Fit3=Fitness_Adrian(N1,Lenski_3,dia_adrian,r0_1,rho_2,mu_2,gamma1)
Fit4=Fitness_Adrian(N1,Lenski_4,dia_adrian,r0_1,rho_3,mu_3,gamma1)

#Vaciado de los ajustes realizados con los distintos parámetros
#de rho_n y mu_N, en una tabla para un mejor uso y las gáficas correctas
#de las simulaciones. 
Ft1=Fitness1(gamma1,r0_1,rho_1,mu_1,dia_adrian)
Fit_1<-data.frame(Ft1,Fit1,Fit2,Fit3,Fit4)

Fit1lm=lm(Ft1~sqrt(Fit1),data=Fit_1[,1:2])
Fit2lm=lm(Ft1~sqrt(Fit2),data=Fit_1[,c(1,3)])
Fit3lm=lm(Ft1~sqrt(Fit3),data=Fit_1[,c(1,4)])
Fit4lm=lm(Ft1~sqrt(Fit4),data=Fit_1[,c(1,5)])

#Generación del valor para el criterio AIC y con el objetivo de determinar que ajuste es mejor.
AIC_Crit_1<-c(AIC(Fit1lm),AIC(Fit2lm),AIC(Fit3lm),AIC(Fit4lm))
Ajuste<-c("Ajuste 1","Ajuste 2", "Ajuste 3", "Ajuste 4")
Selec_Model<-data.frame(Ajuste,AIC_Crit_1)

#Tiempo de paro 2
Lenski_11=Lenski_Exp2(r0_1,rho_1,mu_1,gamma1,N1,K_01,No_K01,dia_adrian,1)
Lenski_22=Lenski_Exp2(r0_1,rho_1,mu_11,gamma1,N1,K_01,No_K01,dia_adrian,1)
Lenski_33=Lenski_Exp2(r0_1,rho_2,mu_2,gamma1,N1,K_01,No_K01,dia_adrian,1)
Lenski_44=Lenski_Exp2(r0_1,rho_3,mu_3,gamma1,N1,K_01,No_K01,dia_adrian,1)

Fit11=Fitness_Adrian(N1,Lenski_11,dia_adrian,r0_1,rho_1,mu_1,gamma1)
Fit22=Fitness_Adrian(N1,Lenski_22,dia_adrian,r0_1,rho_1,mu_11,gamma1)
Fit33=Fitness_Adrian(N1,Lenski_33,dia_adrian,r0_1,rho_2,mu_2,gamma1)
Fit44=Fitness_Adrian(N1,Lenski_44,dia_adrian,r0_1,rho_3,mu_3,gamma1)

Ft2=Fitness1(gamma1,r0_1,rho_1,mu_1,dia_adrian)
Fit_2<-data.frame(Ft2,Fit11,Fit22,Fit33,Fit44)

Fit11lm=lm(Ft2~sqrt(Fit11), data = Fit_2[,1:2])
Fit22lm=lm(Ft2~sqrt(Fit22), data = Fit_2[,c(1,3)])
Fit33lm=lm(Ft2~sqrt(Fit33), data = Fit_2[,c(1,4)])
Fit44lm=lm(Ft2~sqrt(Fit44), data = Fit_2[,c(1,5)])


AIC_Crit_2<-c(AIC(Fit11lm),AIC(Fit22lm),AIC(Fit33lm),AIC(Fit44lm))
Ajuste<-c("Ajuste 1","Ajuste 2", "Ajuste 3", "Ajuste 4")
Selec_Model_1<-data.frame(Ajuste,AIC_Crit_2)

#Gráficas de las Simulaciones del Tiempo de Paro 1
par(mfcol = c(2, 2))

# Gráficos Profundizados en el Crecimiento de la población mutante
plot(eje_x1[1:195],Fit1[1:195],col="skyblue3",lwd=2,ylab="Fitness",xlab="Días",las=1,col.axis="black")
par(new=TRUE)
plot(eje_x1,Ft1,col="maroon4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 1")

plot(eje_x1[612:736],Fit2[612:736],col="goldenrod1",lwd=2,ylab="Fitness",xlab="Días" )
par(new=TRUE)
plot(eje_x1,Ft1,col="maroon4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 2")

plot(eje_x1[5593:5734],Fit3[5593:5734],col="indianred2",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft1,col="maroon4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 3")

plot(eje_x1[2889:3049],Fit4[2889:3049],col="olivedrab3",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft1,col="maroon4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 4")

par(mfcol = c(2, 2))
# Gráficos de todos los días
plot(eje_x1[0:195],Fit1[0:195],col="skyblue3",lwd=2,ylab="Fitness",xlab="Días",las=1,col.axis="black")
par(new=TRUE)
plot(eje_x1,Ft1,col="maroon4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 1")

plot(eje_x1[0:736],Fit2[0:736],col="goldenrod1",lwd=2,ylab="Fitness",xlab="Días" )
par(new=TRUE)
plot(eje_x1,Ft1,col="maroon4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 2")

plot(eje_x1[0:5735],Fit3[0:5735],col="indianred2",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft1,col="maroon4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 3")

plot(eje_x1[0:4050],Fit4[0:4050],col="olivedrab3",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft1,col="maroon4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 4")

##################################################################
#Gráficas de las Simulaciones del Tiempo de Paro 2
par(mfcol = c(2, 2))

# Gráficos Profundizados en el Crecimiento de la población mutante
plot(eje_x1[9402:9763],Fit11[9402:9763],col="chocolate",lwd=2,ylab="Fitness",xlab="Días",las=1,col.axis="black")
par(new=TRUE)
plot(eje_x1,Ft2,col="thistle4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 1")

plot(eje_x1[1:313],Fit22[1:313],col="palegreen3",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft2,col="thistle4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 2")

plot(eje_x1[7350:7822],Fit33[7350:7822],col="wheat3",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft2,col="thistle4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 3")

plot(eje_x1[9350:9822],Fit44[9350:9822],col="steelblue2",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft2,col="thistle4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 4")

# Gráficos de todos los días
par(mfcol = c(2, 2))

plot(eje_x1[0:9763],Fit11[0:9763],col="chocolate",lwd=2,ylab="Fitness",xlab="Días",las=1,col.axis="black")
par(new=TRUE)
plot(eje_x1,Ft2,col="thistle4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 1")

plot(eje_x1[0:313],Fit22[0:313],col="palegreen3",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft2,col="thistle4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 2")

plot(eje_x1[0:7822],Fit33[0:7822],col="wheat3",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft2,col="thistle4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 3")

plot(eje_x1[0:9822],Fit44[0:9822],col="steelblue2",lwd=2,ylab="Fitness",xlab="Días")
par(new=TRUE)
plot(eje_x1,Ft2,col="thistle4",lwd=2,ylab="Fitness",xlab="Días",type="l",axes = FALSE,main="Fitness del Ajuste 4")

##################################
#Funciones Importantes para las simulaciones 
#de las Probabilidades de Fijación

#Función que sirve para obtener las aproximaciones de las probabilidades de
#Fijación repitiendo el experimento (con varios días) m cantidad de veces.
Prob_Fit1=function(r0,rho,mu,gamma,N,K,No_K,t_aum,experimentos,Ind){
  E<-list()
  for(i in 1:experimentos){
    tic()
    E[[i]]=Lenski_Exp1(r0,rho,mu,gamma,N,K,No_K,t_aum,Ind)
    exectime <- toc()
    exectime <- exectime$toc - exectime$tic
    print(exectime)
  }
  
  PROB_FIJ<-list()
  for (j in 1:experimentos) {
    PROB_FIJ[[j]]=E[[j]]$Num_Mut
  }
  
  P<-array(unlist(PROB_FIJ),dim=c(t_aum,length(PROB_FIJ)))
  dimnames(P)<-list(c(1:t_aum),c(1:experimentos))
  Prob.Ki.dado.K0<-array(dim=c(t_aum,length(PROB_FIJ)))
  for(g in 1:t_aum){
    for(s in 1:length(PROB_FIJ)){
      ifelse(P[g,s]==N,Prob.Ki.dado.K0[g,s]<-1,Prob.Ki.dado.K0[g,s]<-0)
    }
  }
  S<-colSums(Prob.Ki.dado.K0)
  P_fijacion<-c()
  for(g in 1:length(PROB_FIJ)){
    if(S[g]>0){
      P_fijacion[g]<-1
    }else{
      P_fijacion[g]<-0
    }
  }
  
  Proba_fij<-0
  Proba_fij=sum(P_fijacion)
  Proba_fij=Proba_fij/experimentos
}

#Función que sirve para obtener las aproximaciones de las probabilidades de
#Fijación repitiendo el experimento (con varios días) m cantidad de veces.
Prob_Fit2=function(r0,rho,mu,gamma,N,K,No_K,t_aum,experimentos,Ind){
  E<-list()
  for(i in 1:experimentos){
    tic()
    E[[i]]=Lenski_Exp2(r0,rho,mu,gamma,N,K,No_K,t_aum,Ind)
    exectime <- toc()
    exectime <- exectime$toc - exectime$tic
    print(exectime)
  }
  
  PROB_FIJ<-list()
  for (j in 1:experimentos) {
    PROB_FIJ[[j]]=E[[j]]$Num_Mut
  }
  
  P<-array(unlist(PROB_FIJ),dim=c(t_aum,length(PROB_FIJ)))
  dimnames(P)<-list(c(1:t_aum),c(1:experimentos))
  Prob.Ki.dado.K0<-array(dim=c(t_aum,length(PROB_FIJ)))
  for(g in 1:t_aum){
    for(s in 1:length(PROB_FIJ)){
      ifelse(P[g,s]==N,Prob.Ki.dado.K0[g,s]<-1,Prob.Ki.dado.K0[g,s]<-0)
    }
  }
  S<-colSums(Prob.Ki.dado.K0)
  P_fijacion<-c()
  for(g in 1:length(PROB_FIJ)){
    if(S[g]>0){
      P_fijacion[g]<-1
    }else{
      P_fijacion[g]<-0
    }
  }
  
  Proba_fij<-0
  Proba_fij=sum(P_fijacion)
  Proba_fij=Proba_fij/experimentos
}

####################################################
#SIMULACIóN DE LAS PROBABILIDADES DE FIJACIÓN

Prob1=Prob_Fit1(r0_1,rho_1,0,gamma1,N1,K_01,No_K01,dia_adrian,experimentos1,0)
Prob2=Prob_Fit1(r0_1,rho_1,0,gamma1,N1,K_01,No_K01,dia_adrian,experimentos1,0)
Prob3=Prob_Fit1(r0_1,rho_2,0,gamma1,N1,K_01,No_K01,dia_adrian,experimentos1,0)
Prob4=Prob_Fit1(r0_1,rho_3,0,gamma1,N1,K_01,No_K01,dia_adrian,experimentos1,0)

Prob11=Prob_Fit2(r0_1,rho_1,mu_1,gamma1,N1,K_01,No_K01,dia_adrian,experimentos1,0)
Prob22=Prob_Fit2(r0_1,rho_1,mu_11,gamma1,N1,K_01,No_K01,dia_adrian,experimentos1,0)
Prob33=Prob_Fit2(r0_1,rho_2,mu_2,gamma1,N1,K_01,No_K01,dia_adrian,experimentos1,0)
Prob44=Prob_Fit2(r0_1,rho_3,mu_3,gamma1,N1,K_01,No_K01,dia_adrian,experimentos1,0)


