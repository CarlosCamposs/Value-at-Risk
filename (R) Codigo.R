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
library(knitr)
library(DescTools)




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
# Creamos la funcion P&L de cada emisora y cada una aparece en un dataframe llamado PL_Portafolio

PL_Portafolio<-data.frame()  

  for (i in 1:2){
      for (j in 1:nrow(tabla_revaluacion)){
        PL_Portafolio[j,i]<-ultimo_precio[i]-tabla_revaluacion[j,i]
      }
  }  

#####
# P&L portafolio diversificado
# Creamos la tabla que representa la suma de cada una de las P&L de las emisoras, a esto se le
# conoce como un "neteo"

  PL_Portafolio2 <-data.frame()
      for(i in 1:nrow(PL_Portafolio)){
        PL_Portafolio2[i,1]<-rowSums(PL_Portafolio)[i]
      }

View(PL_Portafolio2)

#####
# Pegamos esta tabla con la tabla PL_Portafolio, esta ultima tiene las P&L de cada emisora de forma individual
  PL_Portafolio<-cbind(PL_Portafolio,PL_Portafolio2)

  colnames(PL_Portafolio)<-c("WALMEX","FEMSA","Diversificado")
  View(PL_Portafolio)



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

  
#  valores<-cbind(tabla_rendimientos$WALMEX,tabla_rendimientos$FEMSA)
    # gofCopula(copula,valores)    
    # "valores" seran transformados a pseudo-obs
  
# Prueba con ua copula random
#  copulaprueba<-claytonCopula(param=1.5,dim=2)
#  gofCopula(copulaprueba,valores) # p-value de 0.03047
  
  
# -------------------------------------------------------------------------------------------
# Simulacion de rendimientos y obtención del VaR usando copulas
  
VaRCopulas<-data.frame()
tVaRCopulas<-data.frame()  
Diversificado_Copulas <- data.frame()  

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
    # P&L del portafolio diversificado
    # Construimos la P&L que consta del neteo de ambas emisoras
    PL_DiversificadoSIM<-data.frame()
    
    for(j in 1:nrow(PL_Portafolio)){
      PL_DiversificadoSIM[j,1]<-rowSums(PL_EmisorasSIM)[j]
    }

            
#####################################################################################  
  
    #######################
    # Calculamos el VaR tVaR para cada una de las emisoras

  
    for(j in 1:2){
      VaRCopulas[i,j]<-quantile(PL_EmisorasSIM[,j],probs=0.95)
      tVaRCopulas[i,j]<-ES(dist=PL_EmisorasSIM[,j],p_loss = 0.05)
    } # Estos dos dataframes que estamos llenando se crearon antes del "for" que cuenta el numero de simulaciones
    
    #######################
    # Calculamos el VaR tVaR del portafolio diversificado
        
    Diversificado_Copulas[i,1]<-quantile(PL_DiversificadoSIM$V1,probs=0.95) # Calculamos el VaR
    Diversificado_Copulas[i,2]<-ES(dist=PL_DiversificadoSIM$V1,p_loss = 0.05) # Calculamos el tVaR
    
    
}


#####
# VaR y tVaR de cada emisora
  colnames(VaRCopulas)<-c("WALMEX","FEMSA")
  colnames(tVaRCopulas)<-c("WALMEX","FEMSA")  # Cada uno de las filas de estos dataframes representan una simulacion del valor del VaR y tVaR, respectivamente
  
  table1_Copulas<-rbind(colMeans(VaRCopulas),colMeans(tVaRCopulas))  
  colnames(table1_Copulas)<-c("WALMEX", "FEMSA") 
  row.names(table1_Copulas)<-c("VaRCopulas", "tVaRCopulas")
  View(table1_Copulas)

#####
# VaR y tVaR del portafolio diversificado
  View(Diversificado_Copulas)  
  
  table2_Copulas<-cbind(colMeans(Diversificado_Copulas))
  colnames(table2_Copulas) <-"Diversificado" 
  row.names(table2_Copulas) <- c("VaRCopulas","tVaRCopulas")
  View(table2_Copulas)
  
  
# -------------------------------------------------------------------------------------------
# Simulacion Historica

View(PL_Portafolio)

######
# VaR - No diversificado (Simulacion Historica - 1 día)
VaR_SH<-data.frame()
tVaR_SH<-data.frame()
  
  for(i in 1:length(cartera)){
      VaR_SH[1,i]<-quantile(PL_Portafolio[,i],0.95)
      tVaR_SH[1,i]<-ES(dist=PL_Portafolio[,i],p_loss = 0.05)
  }

colnames(VaR_SH)<-c("WALMEX","FEMSA")
colnames(tVaR_SH)<-c("WALMEX","FEMSA")

View(VaR_SH)
View(tVaR_SH)

table1_SH<-rbind(VaR_SH,tVaR_SH)
rownames(table1_SH)<-c("VaR_SH","tVaR_SH")
View(table1_SH)

######
# VaR - Diversificado (Simulacion Historica - 1 día)
View(PL_Portafolio)

table2_SH<-data.frame()
    for(i in 1:length(cartera)){
      table2_SH[1,1]<-quantile(PL_Portafolio[,3],0.95)
      table2_SH[2,1]<-ES(dist=PL_Portafolio[,3],p_loss = 0.05)
    }

  colnames(table2_SH)<-c("Diversificado")
  rownames(table2_SH)<-c("VaR_SH","tVaR_SH")
  View(table2_SH)


# -------------------------------------------------------------------------------------------
# Metodo de Simulacion Montecarlo

#####
# Calculamos las medias y sd de cada columna de los rendimientos de cada emisora

means<-colMeans(tabla_rendimientos) 
sd<-c(sd(tabla_rendimientos$WALMEX),sd(tabla_rendimientos$FEMSA))


#####
# Simulaciones

VaRSM<-data.frame()  
tVaRSM<-data.frame()  
Diversificado_SM<-data.frame() #Servira para almacenar el VaR y tVaR del portafolio diversificado

for(i in 1:5){ # Este "for" corresponde al numero de simulaciones, tarda mucho tiempo, se puede correr una sola simulacion para ver los resultados rapidos
 
  rendimientos_simMC<-data.frame()
  
  for(k in 1:length(cartera)){
    for(j in 1:nrow(tabla_rendimientos)){
      rendimientos_simMC[j,k]<-rnorm(nrow(tabla_rendimientos),mean = means[k],sd = sd[k])[j]
    }
    
  } 
  
  
  
  #######################
  # Revaluacion
  # Una vez obtenido los rendimientos simulados, procedemos a calcular la revaluacion para cada emisora
  
  tabla_revaluacionSM<-data.frame()
  for(k in 1:length(cartera)){
    
    for( j in 1:nrow(rendimientos_simMC)){
      tabla_revaluacionSM[j,k]<-ultimo_precio[,k]*(1+rendimientos_simMC[j,k])
    }
    
  }
  
  #######################
  # P&L indivual
  # Construimos la P&L de cada emisora
  PL_EmisorasSM<-data.frame()
  
  for (k in 1:ncol(tabla_revaluacionSM)){
    
    for(j in 1:nrow(tabla_revaluacionSM)){
      
      PL_EmisorasSM[j,k]<-ultimo_precio[k]-tabla_revaluacionSM[j,k]
    }
  }
  
  #######################
  # P&L Portafolio diversificado
  # Construimos la P&L del portafolio diversificado, la cual resulta de realizar un neteo
  
  PL_SM<-data.frame()
  
  for(k in 1:nrow(PL_EmisorasSM)){
    PL_SM[k,1]<-rowSums(PL_EmisorasSM)[k]
  }
  
  
  #######################
  # VaR y tVaR de cada emisora (Simulacion Montecarlo - 1 día)

  for(j in 1:length(cartera)){
    VaRSM[i,j]<-quantile(PL_EmisorasSM[,j],probs=0.95)
    tVaRSM[i,j]<-ES(dist=PL_EmisorasSM[,j],p_loss = 0.05)
  }
  
  #######################
  # VaR y tVaR del portafolio diversificado (Simulacion Montecarlo - 1 día)
  
    Diversificado_SM[i,1]<-quantile(PL_SM$V1,probs=0.95) # VaR
    Diversificado_SM[i,2]<-ES(dist=PL_SM$V1,p_loss = 0.05) #tVaR
    
}
View(VaRSM)
View(tVaRSM)
View(Diversificado_SM)

# El VaR por metodo de simulacion de Monte Carlo de cada emisora es el promedio de cada columna

######
# VaR - No Diversificado (Simulacion Monte Carlo - 1 día)


table1_SM<-data.frame()

  for (j in 1:length(cartera)){
    table1_SM[1,j]<- colMeans(VaRSM)[j]
    table1_SM[2,j]<- colMeans(tVaRSM)[j]
  }

  colnames(table1_SM)<-c("WALMEX","FEMSA")
  rownames(table1_SM)<-c("VaR_SM", "tVaR_SM")
  View(table1_SM)

######
# VaR - Diversificado (Simulacion Monte Carlo - 1 día)

table2_SM<-data.frame()
    for (j in 1:length(cartera)){
      table2_SM[j,1]<- colMeans(Diversificado_SM)[j]
    }

  colnames(table2_SM)<-c("Diversificado")
  rownames(table2_SM)<-c("VaR_SM", "tVaR_SM")
  View(table2_SM)


  
  
# -------------------------------------------------------------------------------------------
# Metodo de Simulacion Bootstrap

VaRBoots<-data.frame()  
tVaRBoots<-data.frame()  
Diversificado_Boots <- data.frame()  


for (i in 1:5){ # Este "for" es para las simulaciones, se tarda mucho tiempo. Se puede establecer que haga una simulacion para ver los resultados rapidos
  
  # Hacemos un remuestreo del P&L de cada emisora y lo metemos en un dataframe llamado PL_EmisorasBoots  
  # el remuestro de cada emisora se mete en una columna del de dataframe "PL_EmisorasBoots"

  PL_EmisorasBoots<-data.frame()
  PL_DiversificadoBoots <- data.frame()
  
  #####
  # Realizamos el remuestro para cada una de las emisoras
      for(k in 1:length(cartera)){
        for (j in 1:nrow(PL_Portafolio)){
      
          PL_EmisorasBoots[j,k]<-sample(PL_Portafolio[,k], size=nrow(PL_Portafolio), replace = TRUE)[j]
          
        }
      }
  
  #####
  # Realizamos el remuestro del portafolio diversificado 
        for(j in 1:nrow(PL_Portafolio)){
        PL_DiversificadoBoots[j,1] <- sample(PL_Portafolio[,3], size=nrow(PL_Portafolio), replace = TRUE)[j]
      } 
  
  #####
  # Calculamos el VaR y tVaR para cada emisora  
  
    for (j in 1:length(cartera)){
      VaRBoots[i,j]<-quantile(PL_EmisorasBoots[,j],0.95) # VaR
      tVaRBoots[i,j] <- ES(dist=PL_EmisorasBoots[,j],p_loss = 0.05) #tVaR 
    }
  
  #####
  # Calculamos el VaR y tVaR del portafolio diversificado  
  
    Diversificado_Boots[i,1] <- quantile(PL_DiversificadoBoots$V1,0.95) # VaR del portafolio diversificado
    Diversificado_Boots[i,2] <- ES(dist=PL_DiversificadoBoots$V1,p_loss = 0.05) # tVaR del portafolio diversificado
  
  
}
View(VaRBoots)
View(tVaRBoots)
View(Diversificado_Boots)


# El VaR por simulacion Bootstrap de cada emisora es el promedio de cada columna


######
# VaR - No Diversificado (Simulacion Bootstrapping - 1 día)

table1_Boots<-data.frame()

    for (j in 1:length(cartera)){
      table1_Boots[1,j]<- colMeans(VaRBoots)[j]
      table1_Boots[2,j]<- colMeans(tVaRBoots)[j]
    }

  colnames(table1_Boots)<-c("WALMEX","FEMSA")
  rownames(table1_Boots)<-c("VaR_Boots","tVaR_Boots")
  View(table1_Boots)

######
# VaR - Diversificado (Simulacion Bootstrapping - 1 día)
  
table2_Boots<-data.frame()  
  
  for (j in 1:length(cartera)){
    table2_Boots[j,1]<- colMeans(Diversificado_Boots)[j]
  }
  
  colnames(table2_Boots)<-c("Diversificado")
  rownames(table2_Boots)<-c("VaR_SM", "tVaR_SM")
  View(table2_Boots)
  
  
  
# -------------------------------------------------------------------------------------------
# Recapitulacion de resultados
View(table1_SH)
View(table2_SH)
  
View(table1_SM)
View(table2_SM)

View(table1_Boots)
View(table2_Boots)

View(table1_Copulas) # VaR y tVaR no diversificado
View(table2_Copulas)

# table1_X[1,] Contiene los valores del VaR de cada emisora
# table1_X[2,] Contiene los valores del tVaR de cada emisora

# table2_X[1,1] Contiene los valores del VaR del portafolio diversificado
# table1_X[2,1] Contiene los valores del tVaR del portafolio diversificado

Copulas <- cbind(table1_Copulas,table2_Copulas)
SH <- cbind(table1_SH,table2_SH)
SM <- cbind(table1_SM,table2_SM)
Boots <- cbind(table1_Boots,table2_Boots)



#####
# VaR 

VaR_vector<-as.data.frame(
  rbind(Copulas[1,],SH[1,],SM[1,],Boots[1,]))
rownames(VaR_vector)<-c("Copulas","Simulacion Historica","Simulacion Monte-Carlo","Bootstrap")

kable(VaR_vector,digits = 4,caption = "VaR al 95%",col.names = c("WALMEX","FEMSA","Diversificado"))


#####
# tVaR 

tVaR_vector<-as.data.frame(
  rbind(Copulas[2,],SH[2,],SM[2,],Boots[2,]))
rownames(tVaR_vector)<-c("Copulas","Simulacion Historica","Simulacion Monte-Carlo","Bootstrap")

kable(tVaR_vector,digits = 4,caption = "tVaR al 95%",col.names = c("WALMEX","FEMSA","Diversificado"))


# -------------------------------------------------------------------------------------------
