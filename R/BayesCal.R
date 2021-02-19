#####################################Function BayesCal################################################
`BayesCal` <-
function(Lum,alpha=FALSE,supra=FALSE,Temp=600,debug=FALSE) #appel Bug 
	{
	filename<-file.path("BayesCal.bug")

	#rev 8-dec-2015 ajout possibilité supra - en cours

	require(R2WinBUGS)
	require(rv)


	### version alternative avec supra - en cours ##################
	ModelBayesCalSup <- function(){
	#(c) Antoine Zink 2-fev-2016
	#calibration TL (1st growth curves+2nd growth)

	#d'après calibration TL (1st growth) rev 28-jan-2015 

	#prior	
	#N+beta
	X0 ~ dnorm(mu.X0,t.X0)%_%I(,0)
	mbeta~ dnorm(mu.mB0,t.mB0)
	nbeta~dnorm(mu.nB0,t.nB0)
	mu.Y0~ dnorm(Y0,t.Y0) 
	for (i in 1:9){
		t.Delta0[i]~dgamma(a,b)
		}
	
	#supra
	#Xs ~ dnorm(mu.Xs,t.Xs)%_%I(,0)
	msup~ dnorm(mu.mS0,t.mS0)
	nsup~dnorm(mu.nS0,t.nS0)
	mu.Ys~ dnorm(Ys,t.Ys) 
	for (i in 1:9){
		t.DeltaS[i]~dgamma(a,b)
		}

	for (n in 1:600) {
  		prior[n] <- 1/300*step(n-250.5)*step(550.5-n);   # Uniform prior on 250 to 550°C
		}
	Temp2 ~ dcat(prior[]);	
	for (n in 1:600) {
  		prior1[n] <- step(Temp2-n-0.5)*step(n-249.5);   # Uniform prior on 250 to Temp2-10 ; step(>0)=1 ; step(<0)=0		
		}
	l<-sum(prior1[]);
	for (n in 1:600) {
  		prior2[n] <- prior1[n]/l;   # Uniform prior boucle à 1 		
		}

	Temp1~ dcat(prior2[]);
	
	#logical link

		#naturel et N+bêta
		for (i in 1:3){ # i irradiation		
			for (g in 1:3){#g groupe de disques				
				for (j in 1:600){ #j température				
					v[j,i,g]<-step(j-Temp1+0.5)*step(Temp2-j+0.5)			
					YY[j,(g-1)*3+i]<-L[j,i,g]*v[j,i,g]
					YYn[j,(g-1)*3+i]<-LN[j,(g-1)*3+i]*v[j,i,g]
					}
				u[(g-1)*3+i]<-sum(YY[,(g-1)*3+i])
				w[(g-1)*3+i]<-sum(YYn[,(g-1)*3+i])
				y1[(g-1)*3+i]<-u[(g-1)*3+i]/w[(g-1)*3+i]
				D[(g-1)*3+i]<-(y2[i]-y1[(g-1)*3+i])

				}
			y2[i]<-nbeta+mbeta*mu.X[i]
			}

		Y0<- nbeta+mbeta*X0

	#supra en cours
		for (i in 1:3){ # i irradiation		
			for (g in 1:3){#g groupe de disques				
				for (j in 1:600){ #j température				
					vs[j,i,g]<-step(j-Temp1+0.5)*step(Temp2-j+0.5)			
					YYs[j,(g-1)*3+i]<-Ls[j,i,g]*vs[j,i,g]
					YYns[j,(g-1)*3+i]<-LNs[j,(g-1)*3+i]*vs[j,i,g]
					}
				us[(g-1)*3+i]<-sum(YYs[,(g-1)*3+i])
				ws[(g-1)*3+i]<-sum(YYns[,(g-1)*3+i])
				y1s[(g-1)*3+i]<-us[(g-1)*3+i]/ws[(g-1)*3+i]
				Ds[(g-1)*3+i]<-(y2s[i]-y1s[(g-1)*3+i])

				}
			y2[i]<-nsup+msup*mu.X[i]
			}

		Ys<- nsup+msup*Xs


	#data

	p ~ dunif(0,1)
	for (i in 1:9){ 
		group[i] ~ dbern(p)
		index[i] <- group[i] + 1
		mu.D[i,1]<-0
		mu.D[i,2]~dnorm(1,4)
		Delta0[i]<-D[i]-mu.D[i,index[i]]
		mu.Delta0[i]~dnorm(Delta0[i],t.Delta0[i])

		}
	}
	##### fin version alternative ######################################################


	ModelBayesCal <- function(){
	#(c) Antoine Zink 20-Mars-2014
	#calibration TL (1st growth curves)

	#rev. 26-juin-2014
	#rev 26-sept-2014 extrait de test 4 en se limitant au beta
	#rev 1-oct-2014 version correcte pour 3 courbes normalisée de la première chauffe
	#rev 3-oct-2014 version pour l'ensemble des courbes normalisée bêta de la première chauffe
	#rev 10-oct-2014 distribution de t.D selon une loi gamma
	#rev 27-jan-2015 prior basé sur l'analyse de la pente entre 370-400°C
	#rev 28-jan-2015 introduction d'une erreur supplémentaire Delta1 avec une probabilité p (mixture model)


	#prior	
	X0 ~ dnorm(mu.X0,t.X0)%_%I(,0)
	mbeta~ dnorm(mu.mB0,t.mB0)
	nbeta~dnorm(mu.nB0,t.nB0)
	mu.Y0~ dnorm(Y0,t.Y0) 
	for (i in 1:9){
		t.Delta0[i]~dgamma(a,b)
		}
	
	for (n in 1:600) {
  		prior[n] <- 1/300*step(n-250.5)*step(550.5-n);   # Uniform prior on 250 to 550°C
		}
	Temp2 ~ dcat(prior[]);	
	for (n in 1:600) {
  		prior1[n] <- step(Temp2-n-0.5)*step(n-249.5);   # Uniform prior on 250 to Temp2-10 ; step(>0)=1 ; step(<0)=0		
		}
	l<-sum(prior1[]);
	for (n in 1:600) {
  		prior2[n] <- prior1[n]/l;   # Uniform prior boucle à 1 		
		}

	Temp1~ dcat(prior2[]);
	
	#logical link

		#naturel et N+bêta
		for (i in 1:3){ # i irradiation		
			for (g in 1:3){#g groupe de disques				
				for (j in 1:600){ #j température				
					v[j,i,g]<-step(j-Temp1+0.5)*step(Temp2-j+0.5)			
					YY[j,(g-1)*3+i]<-L[j,i,g]*v[j,i,g]
					YYn[j,(g-1)*3+i]<-LN[j,(g-1)*3+i]*v[j,i,g]
					}
				u[(g-1)*3+i]<-sum(YY[,(g-1)*3+i])
				w[(g-1)*3+i]<-sum(YYn[,(g-1)*3+i])
				y1[(g-1)*3+i]<-u[(g-1)*3+i]/w[(g-1)*3+i]
				D[(g-1)*3+i]<-(y2[i]-y1[(g-1)*3+i])

				}
			y2[i]<-nbeta+mbeta*mu.X[i]
			}

		Y0<- nbeta+mbeta*X0



	#data

	p ~ dunif(0,1)
	for (i in 1:9){ 
		group[i] ~ dbern(p)
		index[i] <- group[i] + 1
		mu.D[i,1]<-0
		mu.D[i,2]~dnorm(1,4)
		Delta0[i]<-D[i]-mu.D[i,index[i]]
		mu.Delta0[i]~dnorm(Delta0[i],t.Delta0[i])

		}

	}

	B<-Lum$b
	A<-Lum$a
	N<-Lum$n
	S<-Lum$bsup
	Ns<-Lum$nsup

	if (alpha){
		N<-N[,,seq(1,13)]
		}
	else {
		N<-N[,,seq(1,9)]
		}

	supra<-Lum$supra
	Doseb<-Lum$beta
	Dosea<-Lum$alpha

	write.model(ModelBayesCal,filename)#model.bug)
	model.file<-ModelBayesCal

	if (supra) {	
		write.model(ModelBayesCalSup,filename)#model.bug)
		model.file<-ModelBayesCalSup
		}


	if (alpha) {	
		write.model(ModelBayesCala,filename)#model.bug)
		model.file<-ModelBayesCala
		}

	#data

	mu.X <<-Doseb;
	t.X<<-rep(0.000000000000001,3);

	mu.Delta0<-rep(0,9)
	mu.Delta1<-rep(1,9)

	mu.X0<<- Doseb[2]
	sd.X0<<-mu.X0/2
	t.X0<<-1/sd.X0^2

	mu.Y0<<-0
	t.Y0<<-10000

	mu.Ys<<-0
	t.Ys<<-10000

	L<<-B[,,seq(1,3)] 
	Ls<<-S
	LA<<-A
	LN<<-N
	LNs<<-Ns

	my<-rep(0,9)
	mx<-rep(0,9)
	sigy<-rep(0,3)
	for (j in 1:3){ 
	for (i in 1:3)
		{
		my[(j-1)*3+i]<-sum(B[seq(370,400),i,j])/sum(N[seq(370,400),(j-1)*3+i])  
		mx[(j-1)*3+i]<-mu.X[i]
		} 
	}

	mu.mB0=coef(lm(my~mx))[2] #slope 
	t.mB0=6000000
	mu.nB0=coef(lm(my~mx))[1]  #intercept
	t.nB0=10


	sigy<-(my-mu.nB0-mu.mB0*mx)

	n0<<-2*mean(sigy^2)^2/var(sigy^2)+4
	S0<<-(n0-2)*mean(sigy^2)/n0
	a<<-n0/2
	b<<-n0*S0/2


	data <- c("L","LN","mu.Delta0","mu.X","mu.Y0","t.Y0","mu.mB0","t.mB0","mu.nB0","t.nB0","mu.X0","t.X0","a","b")


	if (supra){
			data <- c(data,"Ls","LNs","mu.Ys","t.Ys","mu.mS0","t.mS0","mu.nS0","t.nS0","mu.Xs","t.Xs")
			}


	#parameters
	inits <- function(){
    			list(Temp1=251,Temp2=549)

		}

	parameters <- c("X0","mbeta","nbeta","Temp2","Temp1","y1","y2","D","Delta0","p","mu.D")

	if (alpha) {
		data <- c(data,"LA")
		parameters <- c(parameters,"Xba")
		}


	#simulation
	Cal.sim <<- bugs(data, parameters, inits=inits, model.file,
    	n.chains=3, n.iter=10,codaPkg=FALSE, n.sims=10,n.burnin=0,
    	bugs.directory="d:/Program Files/WinBUGS14/",debug=debug,DIC=TRUE)

	Cal.sim
	#detach(Cal.sim)

	}

