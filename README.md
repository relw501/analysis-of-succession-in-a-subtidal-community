# analysis-of-succession-in-a-subtidal-community
# define transition matrix
PP <- matrix(c(.771, .102,.017,.004,.015,.001,.018,.012,.002,.014,.003,.002,.005,.003,.029, 
              .145, .609,.031,.011,.028,.005,.022,.025,.011,.015,.012,.008,.005,.004,.069,
              .052, .061,.710,.004,.020,.004,.008,.008,.025,.003,.005,.007,.002,.008,.084,
              .017, .054,.006,.839,.005,.000,.004,.006,.008,.004,.006,.011,.000,.003,.036,
              .117, .218,.035,.004,.404,.008,.033,.032,.013,.007,.006,.005,.006,.005,.108,
              .009, .024,.012,.000,.016,.863,.001,.007,.016,.003,.004,.007,.000,.000,.036,
              .241, .223,.051,.016,.080,.024,.105,.041,.014,.033,.025,.005,.014,.012,.115,
              .199, .235,.038,.018,.089,.007,.044,.154,.015,.027,.016,.020,.009,.009,.122,
              .056, .147,.026,.011,.020,.006,.011,.026,.586,.021,.006,.005,.001,.005,.074,
              .309, .228,.031,.010,.027,.006,.042,.031,.010,.165,.013,.008,.012,.006,.104,
              .056, .222,.028,.008,.036,.000,.025,.020,.007,.007,.507,.002,.005,.003,.076,
              .025, .068,.018,.030,.016,.000,.010,.016,.004,.003,.001,.537,.003,.003,.266,
              .321, .179,.023,.000,.063,.000,.030,.020,.003,.020,.017,.000,.248,.000,.076,      
              .158, .448,.018,.018,.085,.006,.030,.018,.018,.030,.006,.006,.000,.030,.127, 
              .101, .320,.025,.009,.062,.005,.048,.034,.013,.031,.017,.017,.011,.013,.294),
            nrow=15, ncol=15, byrow=TRUE)


x<- c(1/15,1/15,1/15,1/15,1/15,1/15,1/15,1/15,1/15,1/15,
      1/15,1/15,1/15,1/15,1/15)

#find stationary dist
n<-0 #current step number
n <- n+1
x <- x %*% PP
barplot(x,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Probabilities after step n=", n))

eigen(PP)

P <- matrix(c(.773, .102,.017,.004,.015,.001,.018,.012,.002,.014,.003,.002,.005,.003,.029, 
              .145, .609,.031,.011,.028,.005,.022,.025,.011,.015,.012,.008,.005,.004,.069,
              .052, .061,.709,.004,.020,.004,.008,.008,.025,.003,.005,.007,.002,.008,.084,
              .017, .054,.006,.840,.005,.000,.004,.006,.008,.004,.006,.011,.000,.003,.036,
              .117, .218,.035,.004,.403,.008,.033,.032,.013,.007,.006,.005,.006,.005,.108,
              .009, .024,.012,.000,.016,.865,.001,.007,.016,.003,.004,.007,.000,.000,.036,
              .242, .223,.051,.016,.080,.024,.105,.041,.014,.033,.025,.005,.014,.012,.115,
              .199, .233,.038,.018,.089,.007,.044,.154,.015,.027,.016,.020,.009,.009,.122,
              .056, .147,.026,.011,.020,.006,.011,.026,.585,.021,.006,.005,.001,.005,.074,
              .307, .228,.031,.010,.027,.006,.042,.031,.010,.165,.013,.008,.012,.006,.104,
              .056, .222,.028,.008,.036,.000,.025,.020,.007,.007,.505,.002,.005,.003,.076,
              .025, .068,.018,.030,.016,.000,.010,.016,.004,.003,.001,.537,.003,.003,.266,
              .321, .179,.023,.000,.063,.000,.030,.020,.003,.020,.017,.000,.248,.000,.076,      
              .158, .450,.018,.018,.085,.006,.030,.018,.018,.030,.006,.006,.000,.030,.127, 
              .101, .320,.025,.009,.062,.005,.048,.034,.013,.031,.017,.017,.011,.013,.294),
            nrow=15, ncol=15, byrow=TRUE)

for (i in 1:15)
  {cumsum(P[i,][15])}
   
n<-0 #current step number
x<- c(1/15,1/15,1/15,1/15,1/15,1/15,1/15,1/15,1/15,1/15,
      1/15,1/15,1/15,1/15,1/15)
n <- n+1
x <- x %*% P
barplot(x,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Probabilities after step n=", n))
w<-x
w
#species level successions
#disturbances
D<-c(0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0)
for (i in 1:15)
  {D[i] = P[i,15]}
print(D)
  
barplot(D,
          names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
          main = paste("Disturbance")) 

#colonisation
C<-c(0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0)
for (i in 1:15)
{C[i] = P[15,i]}
print(C)

barplot(C,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Colonisation")) 

#persistence
PE<-c(0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0)
for (i in 1:15)
{PE[i] = P[i,i]}
print(PE)

barplot(PE,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Persistence")) 

#replacement of
RO<-c(0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0)
for (i in 1:14)
{  RO[i] = 1- P[i,i] - P[i,15]}
RO[15] = 1- P[15,15]
print(RO)

barplot(RO,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Replacement of")) 

#replacement by
RB<-c(0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0)

for (i in 1:14)
{   RB[i] = (1/13)*(P[1,i]+P[2,i]+P[3,i]+P[4,i]+
                    P[5,i]+P[6,i]+P[7,i]+P[8,i]+
                    P[9,i]+P[10,i]+P[11,i]+P[12,i]+
                    P[13,i]+P[14,i]-P[i,i])}
RB[15] = (1/14)*(P[1,15]+P[2,15]+P[3,15]+P[4,15]+
  P[5,15]+P[6,15]+P[7,15]+P[8,15]+
  P[9,15]+P[10,15]+P[11,15]+P[12,15]+
  P[13,15]+P[14,15])
print(RB)

barplot(RB,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Replacement by")) 

# Community level transitions
#colonisation
CC=1-P[15,15]
CC

#disturbance
CDL<-c(0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0)
ww = w
ww[15]=0
for (i in 1:14)
{CDL[i]=w[i]*P[i,15]}
CD=(cumsum(CDL)[14]/cumsum(ww)[15])
CD

#replacement of
CROI<-c(0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0)
for (i in 1:14)
{CROI[i]=w[i]*RO[i]}
CRO = cumsum(CROI)[14]/cumsum(ww)[15]
CRO

#replacement by
CRBI<-c(0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0)
for (i in 1:14)
{CRBI[i]=w[i]*RB[i]}
CRB=cumsum(CRBI)[14]/cumsum(ww)[15]
CRB


#Successional dynamics

#Turnover rate
T<-c(0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0)
for (i in 1:15)
{T[i]=1-P[i,i]}
T

#Turnover time
t<-c(0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0)
for (i in 1:15)
{t[i]=1/T[i]}
t
barplot(t,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Turnover times")) 

#mean biotic turnover time
t2<-c(0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0)
for (i in 1:14)
{t2[i]=t[i]*w[i]}
tbar=cumsum(t2)[14]/cumsum(ww)[15]
tbar

#recurrence time
RT<-c(0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0)
for (i in 1:15)
{RT[i]=(1-w[i])/(w[i]*(1-P[i,i]))}
barplot(RT,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Recurrence times")) 

RT[15] #average time between disturbances
RT[15]/tbar #point changes its biotic state this many time between disturbances on average.
#mean recurrence time
MRTI<-c(0,0,0,0,0,
         0,0,0,0,0,
         0,0,0,0,0)
for (i in 1:15)
{MRTI[i]=w[i]*RT[i]}
MRT=cumsum(MRTI)
MRT[15]

#mean biotic recurrence time
BRT=cumsum(MRTI)[14]/cumsum(w)[14]
BRT

#entropy, complexity and predictability fo succession
# entropy of P[j,.] 
H<-c(0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0)
HI<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
           nrow=15, ncol=15, byrow=TRUE)

for (i in 1:15)
{for (j in 1:15) 
{if (P[i,j]==0)
{HI[i,j]=0}
  else {HI[i,j] = -P[i,j]*log(P[i,j])}}}

for (i in 1:15)
{H[i]=cumsum(HI[i,])[15]}


#max value is -log(1/15) so rescale
for (i in 1:15)
{H[i]=H[i]/(-log(1/15))}
H

barplot(H,
        names.arg =c("HY1","CRI","MXY","MYC","FIL","URT","ASC","APL","HY2","IDM","COR","MET","PAR","SPI","BR"),
        main = paste("Entropy")) 

#entropy of whole markov chain
HH <- c(0,0,0,0,0,
        0,0,0,0,0,
        0,0,0,0,0)
for (i in 1:15)
{HH[i]=w[i]*H[i]}
 HHH=cumsum(HH)[15]
 HHH
 
 #damping ratio
 eig<-eigen(P)
 A<-eig$values[1]
 B<-eig$values[2]
 A/B

 C<-log(A/B) #convergence rate per year
 
 #half life
 log(2)/C
 
 #dobrushkins coeff
 D<-matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
            nrow=15, ncol=15, byrow=TRUE)
 I15<-diag(15)
 for (i in 1:15)
 {for (j in 1:15)
 {D[i,j]=norm(P[j,]*I15-P[i,]*I15)}}
 D
 
 
