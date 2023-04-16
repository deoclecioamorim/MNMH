
######################################################################
##--Modelagem da hormesis por regressão não linear multivariada
######################################################################
##--Chamando os pacotes e o conjunto de dados
library(readxl)
library(MVN)
library(minpack.lm)
library(MASS)
library(car)
Dados<-read_excel("Dados/Dados.xlsx",sheet =1)
##--Ajustando o modelo de Modelo de Brain-Cousens (1989)
##--Etapa I: Ajuste individual dos modelos de regressão não linear

##--Altura de planta (AP)
c=0 #Parâmetros fixado
mode1.AP <-nlsLM(plant_height~
                (c+(((d-c)+(f*Dose))/(1+exp(b*log(Dose/E))))),
                start =list(d=54,E=100,b=2,f=11),
                control = nls.lm.control(maxiter = 600),
                algorithm="LM",trace= T,data =Dados)

##--Área folia (AF)
c=0
mode1.AF <- nlsLM(leaf_area~
                      (c+(((d-c)+(f*Dose))/(1+exp(b*log(Dose/E))))),
                       start =list(d=416,E=100,b=2,f=100),
                       control = nls.lm.control(maxiter = 600),
                       algorithm="LM",trace= T,data =Dados)

##--Massa de matéria seca (MMS)
mode1.MMS <- nlsLM(dry_matter_mass~ 
                     (c+(((d-c)+(f*Dose))/(1+exp(b*log(Dose/E))))),
                      start =list(d=6,E=40,c=1.5,b=2,f=2),
                      control = nls.lm.control(maxiter = 600),
                      algorithm="LM",trace= T,data =Dados)

##--Etapa II: Obtenção da estimativa da matriz de variâncias
##--e covariâncias dos resíduos

N<-length(Dados$plant_height) #comprimento dos dados
M<-3 #Três modelos univariados

#Construindo a matriz Sigma para as variáveis AP, AF e MMS
resid_m1<-as.matrix(residuals(mode1.AP))
(resid_m1<-t(resid_m1))
resid_m2<-as.matrix(residuals(mode1.AF))
(resid_m2<-t(resid_m2))
resid_m3<-as.matrix(residuals(mode1.MMS))
(resid_m3<-t(resid_m3))

Matrix.erros_m<-matrix(c(resid_m1,resid_m2,resid_m3),
                       nrow = M, ncol = N,byrow = T)

##--Matriz Sigma
SIGMA<-((Matrix.erros_m)%*%t(Matrix.erros_m))/N

##--Etapa III: Composição do modelo não linear multivariado
##--e obtenção da otimização única
##--Inversa de Moore Penrose
inv_Sigma<-ginv(SIGMA) 
cf <- chol(inv_Sigma) #Fator de Cholesky

j<-matrix(1,N,1) #Vetor de 1's
In<-diag(1,N,N)  #Matriz identidade
p1<-kronecker(diag(1, M), j) %*% cf[,1] 
p2<-kronecker(diag(1, M), j) %*% cf[,2]
p3<-kronecker(diag(1, M), j) %*% cf[,3]

#Y empilhado
y_emp<-c(Dados$plant_height,Dados$leaf_area,Dados$dry_matter_mass)
Yr<-kronecker(cf, In)%*%y_emp

#####################################################################
# Modelo não linear multivariado: Resultado apresentado na Tabela 4.1
######################################################################
##--Parâmetros fixados
c1=0;
c2=0;
mult1<-nlsLM(Yr~p1*(c1+(((d1-c1)+(f1*Dose))/(1+exp(b1*log(Dose/E1)))))+
               p2*(c2+(((d2-c2)+(f2*Dose))/(1+exp(b2*log(Dose/E2)))))+
               p3*(c3+(((d3-c3)+(f3*Dose))/(1+exp(b3*log(Dose/E3))))),
             start =list(d1=44, f1=0.7, b1=2, E1=55,
                         d2=300,f2=5, b2=2, E2=60,
                         c3=3, d3=4,f3=0.08, b3=2, E3=40),
             algorithm="LM",trace= T,data =Dados)
summary(mult1)

##--Teste de hormesis no contexto multivariado empregando 
##--o teste da razão de verossimilhanças

#--Modelo reduzido
mult2<-nlsLM(Yr~p1*(c1+(((d1-c1))/(1+exp(b1*log(Dose/E1)))))+
               p2*(c2+(((d2-c2))/(1+exp(b2*log(Dose/E2)))))+
               p3*(c3+(((d3-c3))/(1+exp(b3*log(Dose/E3))))),
             start =list(d1=44, b1=2, E1=55,
                         d2=300, b2=2, E2=60,
                         c3=3,   d3=4, b3=2, E3=40),
             algorithm="LM",trace= T,data =Dados)

summary(mult2)

##--Teste LRT
anova(mult2,mult1)
qf(0.95,3,227)

##--Ajuste para características quantitativas sub-NOAEL
#######################################################################
# Ajustes univariados: Dose-M 
#######################################################################
##--Dose M
##--Altura de planta (AP)
c1=0;
M_AP<-nlsLM(plant_height~
              (c1+((d1-c1)+(f1*Dose))/(1+(f1*M1/(((d1-c1)*b1)
              -f1*M1*(1-b1)))*exp(b1*log(Dose/M1)))),
              start =list(d1=44, M1=20,b1=2,f1=0.7),
              algorithm="LM",trace= T,
              control = nls.lm.control(maxiter = 600),data =Dados)

##--Área foliar (AF)
c2=0
M_AF<-nlsLM(leaf_area~
            (c2+((d2-c2)+(f2*Dose))/(1+(f2*M2/(((d2-c2)*b2)
            -f2*M2*(1-b2)))*exp(b2*log(Dose/M2)))),
            start =list(d2=300,M2=30,b2=2,f2=5),
            algorithm="LM",trace= T,
            control = nls.lm.control(maxiter = 600),data =Dados)

##--Massa de matéria seca (MSS)
M_MSS <-nlsLM(dry_matter_mass~(c3+((d3-c3)+(f3*Dose))
              /(1+(f3*M3/(((d3-c3)*b3)-f3*M3*(1-b3)))*
              exp(b3*log(Dose/M3)))),
              start =list(d3=4,  M3=20,b3=2,f3=0.08,c3=2),
              algorithm="LM",trace= T,
              control = nls.lm.control(maxiter = 600),data =Dados)
##--Extraindo informações para construção da Tabela S4
summary(M_AP)
summary(M_AF)
summary(M_MSS)
######################################################################
#         Dose-M: Modelo não linear multivariado 
######################################################################
doseM<-nlsLM(Yr~p1*(c1+((d1-c1)+(f1*Dose))/(1+(f1*M1/
            (((d1-c1)*b1)-f1*M1*(1-b1)))*exp(b1*log(Dose/M1))))+
            p2*(c2+((d2-c2)+(f2*Dose))/(1+(f2*M2/
            (((d2-c2)*b2)-f2*M2*(1-b2)))*exp(b2*log(Dose/M2))))+
            p3*(c3+((d3-c3)+(f3*Dose))/(1+(f3*M3/
            (((d3-c3)*b3)-f3*M3*(1-b3)))*exp(b3*log(Dose/M3)))),
             start =list(d1=44, M1=20,b1=2,f1=0.7,
                         d2=300,M2=30,b2=2,f2=5,
                         d3=4,  M3=20,b3=2,f3=0.08,c3=2),
                         algorithm="LM",trace= T,
          control = nls.lm.control(maxiter = 600),data =Dados)

summary(doseM)
##--Dose M média
car::deltaMethod(doseM,"(M1+M2+M3)/3",vcov. = vcov(doseM))




