# ///////////////////////////////////////////////////////////////////////////////////////////
# -------------------------------------------------------------------------------------------
# Calculo del VaR para un portafolio de dos activos usando Teoria de Cópulas
# -------------------------------------------------------------------------------------------
# ///////////////////////////////////////////////////////////////////////////////////////////

#####
# Cargamos las bibliotecas necesarias
library(quantmod)
library(rriskDistributions)  
library(fitdistrplus)
library(MASS) 
library(actuar) 
library(univariateML)
library(VGAM)
library(lcopula)
library(VineCopula)
library(VC2copula)
library(scatterplot3d)
library(cvar)  



#####
# Cargamos la serie de datos
cartera=c("WALMEX.MX","KOFUBL.MX")
getSymbols(cartera,src = "yahoo",from="2020-01-01", to="2022-12-31")


#####
# Guardamos en unas variables los precios de cierre de cada emisora
p_WALMEX<-WALMEX.MX$WALMEX.MX.Close
p_FEMSA<-KOFUBL.MX$KOFUBL.MX.Close


#####
# Creamos una tabla donde juntamos ambos precios de cierre
tabla_precios<-cbind(p_WALMEX,p_FEMSA)
colnames(tabla_precios)<-c("WALMEX","FEMSA")


#####
# Definimos la tabla como un dataframe para su mejor manejo
tabla_precios<-as.data.frame(tabla_precios)


#####
# Calculamos los rendimientos  
tabla_rendimientos<-data.frame()

  for(i in 1:2){
      for (j in 1:length(tabla_precios$WALMEX)){
      tabla_rendimientos[j,i]<-tabla_precios[j+1,i]/tabla_precios[j,i]-1
      }  
  } 


#####
# Cambiamos el nombre de las columnas al dataframe   
colnames(tabla_rendimientos)<-c("WALMEX","FEMSA")


#####
# Eliminamos la ultima observacion pues es un dato NA  
tabla_rendimientos<-tabla_rendimientos[-c(nrow(tabla_rendimientos)),]


#####
# Guardamos el ultimo valor observado de los precios de cada una de las emisoras 
ultimo_precio<-tabla_precios[c(nrow(tabla_precios)),]
ultimo_precio<-as.data.frame(ultimo_precio)


#####
# Con los rendimientos calculados previamente y el ultimo precio observado de cada emisora, hacemos
# la revaluacion.

tabla_revaluacion<-data.frame()
  for(i in 1:2){
    for( j in 1:nrow(tabla_rendimientos)){
      tabla_revaluacion[j,i]<-ultimo_precio[,i]*(1+tabla_rendimientos[j,i])
    } 
  }

colnames(tabla_revaluacion)<-c("WALMEX","FEMSA")


#####
# P&L individual
# Creamos la funcion P&L de cada emisora y las juntamos en un dataframe llamado PL_Portafolio

PL_Portafolio<-data.frame()  

  for (i in 1:2){
      for (j in 1:nrow(tabla_revaluacion)){
        PL_Portafolio[j,i]<-ultimo_precio[i]-tabla_revaluacion[j,i]
      }
  }  


hist(PL_Portafolio$WALMEX) # Revisamos la distribucion de la P&L


# -------------------------------------------------------------------------------------------
# Buscamos que distribucion de probabilidad se ajusta a los rendimientos de WALMEX
  hist(tabla_rendimientos$WALMEX, freq = F)
  
  fit.cont(tabla_rendimientos$WALMEX)  
      # La mejor distribucion es la Logistica

#####
# Hallamos los parametros de la distribucion logistica  
  m1<-model_select(tabla_rendimientos$WALMEX,models="logis",criterion="aic")
  coef(m1) 


#####
# Calculamos los valores de la cdf de la logistica, estos valores corresponden a los valores de u1
# de la copula  
    u1<-plogis(tabla_rendimientos$WALMEX,location=coef(m1)[1],scale=coef(m1)[2])

  
  
# -------------------------------------------------------------------------------------------
# Buscamos que distribucion de probabilidad se ajusta a los rendimientos de FEMSA
  hist(tabla_rendimientos$FEMSA, freq = F)
  fit.cont(tabla_rendimientos$FEMSA)  
  # La mejor distribucion es la Logistica
  

# Hallamos los parametros de la distribucion logistica  
  m2<-model_select(tabla_rendimientos$FEMSA,models="logis",criterion="aic")
  coef(m2) 
  

#####
# Calculamos los valores de la cdf de la logistica, estos valores corresponden a los valores de u2
# de la copula  

  u2<-plogis(tabla_rendimientos$FEMSA,location=coef(m2)[1],scale=coef(m2)[2])
  

#####
# Grafico de los puntos u1,u2
  plot(u1,u2)  
  

##### Para checar dependencia
  
  cor.test(x=tabla_rendimientos$WALMEX,
           y=tabla_rendimientos$FEMSA,method = "spearman")
    # Comentario: comprobamos que efectivamente existe dependencia entre los valores de u1 y u2

  valores<-cbind(tabla_rendimientos$WALMEX,tabla_rendimientos$FEMSA)
  K.plot(valores)
    # Comentario: existe una dependencia positiva, pero no es muy fuerte
  
  
# -------------------------------------------------------------------------------------------
# Buscamos determinar cual sera la mejor copula que describa el comportamiento de nuestros datos (u1,u2)
  
  BiCopSelect(u1,u2,familyset=NA,selectioncrit = "AIC")
  # La mejor copula es la Survival Gumbel
  # par=1.15, tau=0.13
  
  
##### Construccion de la copula 
  copula<-surGumbelCopula(param = 1.15)
  

##### Construccion de la funcion de distribucion conjunta F(x,y)
  
  DistConj<-mvdc(copula,margins=c("logis","logis"),
                 paramMargins=list(
                   list(location=coef(m1)[1],scale=coef(m1)[2]),
                   list(location=coef(m2)[1],scale=coef(m2)[2])
                 )
  )
  
##### Graficos de mvdc

# Simulamos 1000 valores de X e Y
  z <- rMvdc(n=1000,DistConj)
  
# Creamos el grafico de la pdf
  pdf<-dMvdc(z,DistConj) 
  scatterplot3d(z[,1],z[,2],pdf,highlight.3d=T)
  
# Creamos el grafico de la cdf
  cdf<-pMvdc(z,DistConj) 
  scatterplot3d(z[,1],z[,2],cdf,highlight.3d=T)

  
  valores<-cbind(tabla_rendimientos$WALMEX,tabla_rendimientos$FEMSA)
  gofCopula(copula,valores)    
    # "valores" seran transformados a pseudo-obs
  
# Prueba con ua copula random
  copulaprueba<-claytonCopula(param=1.5,dim=2)
  gofCopula(copulaprueba,valores) # p-value de 0.03047
  
  
# -------------------------------------------------------------------------------------------
# Simulacion de rendimientosy obtención del VaR
  
VaR95<-data.frame()
tVaR95<-data.frame()  


  for(i in 1:5){ # Este "for" corresponde al numero de simulaciones, tarda mucho tiempo, se puede correr una sola simulacion para ver los resultados rapidos
    
    rendimientos_sim<-data.frame()
    
   
    for(k in 1:2){
      for(j in 1:nrow(tabla_rendimientos)){
        rendimientos_sim[j,k]<-rMvdc(n=nrow(tabla_rendimientos),DistConj)[j,k]
      }
    }# Simula rendimientos de X e Y, ie, de ambas emisoras 
    
    
    
    #######################
    # Revaluacion
    # Una vez obtenido los rendimientos simulados, procedemos a calcular la revaluacion para cada emisora
    
    tabla_revaluacionSIM<-data.frame()
    for(k in 1:2){
      
      for( j in 1:nrow(rendimientos_sim)){
        tabla_revaluacionSIM[j,k]<-ultimo_precio[,k]*(1+rendimientos_sim[j,k])
      }
      
    }
    
    #######################
    # P&L indivual
    # Construimos la P&L de cada emisora
    PL_EmisorasSIM<-data.frame()
    
    for (k in 1:2){
      
      for(j in 1:nrow(tabla_revaluacionSIM)){
        
        PL_EmisorasSIM[j,k]<-ultimo_precio[k]-tabla_revaluacionSIM[j,k]
      }
    }
    
    
    #######################
    # Para cada emisora, calculamos el VaR y los guardamos en un dataframe, de modo que
    # tendremos 10 valores de VaR al 95% de confianza y se guardan todos en una fila
    # del dataframe VARSM95 (para los otros dos es analogo)
    
    
#####################################################################################  
# Todo lo que esta abajo se hizo para calcular el VaR y ES de forma manual
# Calculamos el tVaR

    
    for(j in 1:2){
      VaR95[i,j]<-quantile(PL_EmisorasSIM[,j],probs=0.95)
      tVaR95[i,j]<-ES(dist=PL_EmisorasSIM[,j],p_loss = 0.05)
    }
    

        
    
}
colnames(VaR95)<-c("WALMEX","FEMSA")
colnames(tVaR95)<-c("WALMEX","FEMSA")  
  

# VaR y tVaR
colMeans(VaR95)
colMeans(tVaR95)
  
#table<-rbind("VaR al 95%", "tVaR al 95%")  
#table<-cbind(colMeans(VaR95),colMeans(tVaR_95))  




  
  
