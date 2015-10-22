# R code used to generate Figs 2–6 in Harris et al., Nat Chem Biol 2015
# primary code written by Darren Tyson
#
#
source('NatMethFxns.r')

#########################################################################
# 						  
# 							FIGURES 2 & 3
# 						  
#########################################################################

#########################################################################
# 						  
# 							Set parameter values 
# 						  
#########################################################################
drug.conc		<-	signif(10^((-110:-40)/10),3)		# concentrations of drug (1e-11 to 1e-4 M)
times			<-	0:120								# hours
states			<-	c(Cell=100,CellPrime=0)

k_on	<-	c(1e8,1e8,1e5)
k_off	<-	c(1,1,1e-3)
DIP0	<-	c(0.06,0.01,0.06)*log(2)					# parameters in base e to solve ODEs
DIPmax	<-	c(-0.03,-0.005,-0.03)*log(2)				
param	<-	rbind(k_on=k_on,k_off=k_off,DIP0=DIP0,DIPmax=DIPmax)
colnames(param)	<-	c('fast','slow','fast.ne')		# ne = non-equilibrium
param	<-	as.matrix(param)


drugCols	<-	colorpanel(n=length(drug.conc),low='blue',mid=grey(0.5),high='red')
timeCols	<-	topo.colors(120)
timeCols	<-	sub(topo.colors(120)[72], 'red', timeCols)

#########################################################################
# 						  
# 	Generate cell counts for each cell type and drug concentration
# 						  
#########################################################################

cell.counts	<-	list()
for(i in colnames(param))
	cell.counts[[i]]	<-	getCellCount(param[,i])		# call ODE solver with relevant parameters


makeFig2	<-	function(toFile=FALSE)
{
	if(toFile)	pdf(file='Theor time bias.pdf', width=10, height=5)
	if(!toFile)	dev.new(width=10, height=5)
	par(mfrow=c(3,6), mar=c(0.5,3,1,1), oma=c(6,1,2,1))
	for(i in colnames(param))
	{
		dtp	<-	cell.counts[[i]]
		
		# linear growth curve
		plot(1, xlim=c(0,120),ylim=c(25,10000), type='n', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
		axis(side=1, at=c(0:5)*24, padj=-.5, labels=FALSE)
		axis(side=2, at=seq(from=0, to=10, by=2)*1000, labels=c(0,NA,4000,NA,8000,NA), padj=.5)
		for(c in colnames(dtp)[-1])
			lines(dtp$time,dtp[,c],col=drugCols[match(c,drug.conc)], type='l')
		if(i==colnames(param)[1])
		{
			mtext(side=3,'Growth curve', line=1.25, font=2, cex=.75)
			mtext(side=3,'linear scale', line=0.25, font=3, cex=.75)
		}
		if(i==colnames(param)[2])
			mtext(side=2, 'Number of cells', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3])
		{
			axis(side=1, at=c(0:5)*24, padj=-.5)
			lines(0:120,dtp[,'6.31e-07'], col='orange', lwd=3)
			mtext('Time (h)', side=1, font=2, line=2, cex=0.75)
		}
		
		# log growth curve
		plot(1, xlim=c(0,120),ylim=c(-4,6), type='n', xlab=NA, ylab=NA, xaxt='n', yaxt='n')
		axis(side=1, at=c(0:5)*24, padj=-.5, labels=FALSE)
		axis(side=2, at=seq(from=-4, to=6, by=2), padj=.5)
		for(c in colnames(dtp)[-1])
			lines(dtp$time,log(dtp[,c])-log(dtp[,c])[1],col=drugCols[match(c,drug.conc)], type='l')
		if(i==colnames(param)[1])
		{
			mtext(side=3,'Growth curve', line=1.25, font=2, cex=.75)
			mtext(side=3,'log scale', line=0.25, font=3, cex=.75)
		}
		if(i==colnames(param)[2])
			mtext(side=2, 'Population doublings', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3]) 
		{
			lines(0:120,log(dtp[,'6.31e-07'])-log(dtp[,'6.31e-07'])[1], col='orange', lwd=3)
			mtext('Time (h)', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(0:5)*24, padj=-.5)
		}
		
		
		# static response ratio DRC
		dfsdrc	<-	apply(dtp[-1], 1, FUN=function(x) x/x[1])
		plot(1, xlim=c(1e-11,1e-4),ylim=c(-0.6,1.2), type='n', xlab=NA, ylab=NA, log='x', xaxt='n', yaxt='n')
		axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), padj=-.5, labels=FALSE)
		axis(side=2, at=seq(from=-0.5, to=1, by=0.5), padj=.5)
		ec50.static	<-	numeric()
		aa.static	<-	numeric()
		emax.static	<-	numeric()
		
		# dynamic response ratio DRC (add to static plot)

		if(i==colnames(param)[3])
		{
			DIP		<-	1/apply(dtp[48:120,-1],2,
				FUN=function(x) coef(lm(48:120 ~ log2(x), na.action='na.omit')))['log2(x)',]
			stable.range	<-	48:120	
		} else {
			myargs	<-	append(as.list(param[,i]),list(drug=drug.conc))
			DIP	<-	do.call(dipPEA, args=myargs)
			stable.range	<-	1:120	
		}
			
		m.dyn.norm	<-	drm(DIP/DIP[1] ~ drug.conc, fct=LL.4(), na.action='na.omit')
		for(ti in 12:120)
		{
			m.static	<-	drm(dfsdrc[,ti] ~ drug.conc, fct=LL.4(), na.action='na.omit')
			ec50.static	<-	append(ec50.static,signif(ED(m.static, 50, display=FALSE),3)[1])
			emax.static	<-	append(emax.static,coef(m.static)[2])
			aa.static		<-	append(aa.static,getAA(m.static))
			plot(m.static, type='none', add=TRUE, col=ifelse(ti==72,'red',timeCols[ti]), lwd=ifelse(ti==72,2.5,1))
		}
		plot(m.dyn.norm, type='none', add=TRUE)
		abline(h=0, lwd=1.5, lty=2, col=grey(.5))
		if(i==colnames(param)[1])
		{
			mtext(side=3,'Dose-response curve', line=1.25, font=2, cex=.75)
			mtext(side=3,'Static vs. DIP rate', line=0.25, font=3, cex=.75)
		}
		if(i==colnames(param)[2])
			mtext(side=2, 'Response ratio', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3]) 	
		{
			mtext('[drug], log10 M', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), labels=c(-11,-9,-7,-5), padj=-.5)
		}
		
		# Time bias on EC50
		plot(12:120, ec50.static, xlim=c(0,120), ylim=c(1e-9,1e-5), type='n', xlab=NA, ylab=NA, log='y', xaxt='n', yaxt='n')
		points(12:120, ec50.static, col=timeCols[12:120], pch=19, cex=0.5)
		axis(side=1, at=c(0:5)*24, padj=-.5, labels=FALSE)
		axis(side=2, at=c(10^(-9:-5)), labels=c(-9,NA,-7,NA,-5), padj=.5)
		points(stable.range[1],signif(ED(m.dyn.norm, 50, display=FALSE),3)[1], pch=2, cex=1.5)
		points(stable.range[-1],rep(signif(ED(m.dyn.norm, 50, display=FALSE),3)[1], length(stable.range[-1])), pch=3, cex=.25)
		
		if(i==colnames(param)[1])
		{
			mtext(side=3,'Drug-exposure-time bias', line=1.25, font=2, cex=.75)
			mtext(side=3,'EC50', line=0.25, font=3, cex=.75)
		}
		if(i==colnames(param)[2])
			mtext(side=2, '[drug], log10 M', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3])
		{
			mtext('Time (h)', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(0:5)*24, padj=-.5)
		}
		
		# Time bias on AA
		plot(12:120,aa.static, type='l', xlim=c(0,120), ylim=c(0,1), xlab=NA, ylab=NA, xaxt='n')
		points(12:120, aa.static, col=timeCols[12:120], pch=19, cex=0.5)
		axis(side=1, at=c(0:5)*24, padj=-.5, labels=FALSE)
		points(stable.range[1],getAA(m.dyn.norm), pch=2, cex=1.5)
		points(stable.range[-1],rep(getAA(m.dyn.norm), length(stable.range[-1])), pch=3, cex=.25)
		if(i==colnames(param)[1])
		{
			mtext(side=3,'Drug-exposure-time bias', line=1.25, font=2, cex=.75)
			mtext(side=3,'AA', line=0.25, font=3, cex=.75)
		}
		if(i==colnames(param)[2])
			mtext(side=2, 'Activity area, AU', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3])
		{
			mtext('Time (h)', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(0:5)*24, padj=-.5)	
		}
		
		# dynamic effect DRC
		plot(1, xlim=c(1e-11,1e-4),ylim=c(-0.036,0.072), type='n', xlab=NA, ylab=NA, log='x', xaxt='n', yaxt='n')
		m.dyn	<-	drm(DIP ~ drug.conc, fct=LL.4(), na.action='na.omit')
		plot(m.dyn, type='none', add=TRUE)
		axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), padj=-.5, labels=FALSE)
		axis(side=2, at=c(seq(-4,6,2)*.01), labels=c(-0.04,NA,0,NA,0.04,NA), padj=.5)
		abline(h=0, lwd=1.5, lty=2, col=grey(.5))
		if(i==colnames(param)[1])
		{
			mtext(side=3,'Dose-response curve', line=1.25, font=2, cex=.75)
			mtext(side=3,'DIP rate', line=0.25, font=3, cex=.75)
		}
		if(i==colnames(param)[2])
			mtext(side=2, 'DIP rate', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3])
		{
			mtext('[drug], log10 M', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), labels=c(-11,-9,-7,-5), padj=-.5)
		}
	}

	# color legend for drug conc
	li	<- as.raster(matrix(drugCols, nrow=1))
	grid.raster(x=0.018, y=0.05, li, height=.01, width=.1, just=c('left','top'), interpolate=FALSE)
	mtext(paste(seq(-11,-5,2),collapse="     "), outer=TRUE, side=1, line=3, cex=0.7, adj=0)

	# color legend for times
	li2 <- as.raster(matrix(timeCols[12:120], nrow=1))
	grid.raster(x=0.495, y=0.05, li2, height=.01, width=.1, just=c('center','top'), interpolate=FALSE)
	mtext(paste(seq(24,120,24),collapse="  "), outer=TRUE, side=1, line=3, cex=0.7)

	if(toFile)	dev.off()
}

makeFig3	<-	function(toFile=FALSE, actual.DTP=TRUE)
{
	fn	<-	ifelse(actual.DTP,
		'Theor drug effect delay (actual NCI DTP).pdf',
		'Theor drug effect delay (ll4 fit to data).pdf')
	if(!toFile)	dev.new(width=3, height=6)
	if(toFile)	pdf(file=fn, width=3, height=6)
	par(mar=c(3,3,1,1), oma=c(0,0,1,0), mfrow=c(2,1))
	dtp	<-	cell.counts[['fast.ne']][,c('1e-11','6.31e-07')]

	#
	# growth curves over time and predicted change from baseline from measurements at distinct times (24, 48, 72, 96 and 120 h)
	#
	plot(0:120,dtp[,2], type='l', xlab=NA, ylab=NA, col='orange', lwd=3,
		ylim=c(0,300), xaxt='n', yaxt='n')
	axis(side=1, at=c(0:5)*24, padj=-.5)
	axis(side=2, at=c(0:3)*100, padj=.5)
	lines(0:120,dtp[,1], lwd=3)
	abline(h=100, lty=2, col=grey(.5))
	for(i in c(24,48,72,96,120))
	{
		curve(x*(dtp[,2][i]-dtp[,2][1])/i+100, from=0, to=120, col=timeCols[i], lwd=2, add=TRUE, lty=4)
		points(i,dtp[,2][i])
	}
	mtext('Time (h)', side=1, line=2, font=2)
	mtext('Number of cells', side=2, line=2, font=2)
	mtext('Effects of single drug conc over time', outer=TRUE, side=3, line=1)
	legend('bottomright', bty='n', legend=c('no drug','630 nM'), col=c('black','orange'), cex=0.75, lty=1, lwd=3)
	legend('bottomleft', bty='n', legend=seq(from=24,to=120, by=24), col=timeCols[seq(from=24,to=120, by=24)], 
		cex=0.75, lty=4, lwd=3, y.intersp=.75)

	#
	# dose–response curves illustrating time-dependent bias
	#
	dtp	<-	cell.counts[['fast.ne']]
	DIP		<-	1/apply(dtp[61:120,-1],2,FUN=function(x) coef(lm(61:120 ~ log2(x))))['log2(x)',]
	conc	<-	as.numeric(names(DIP))		
	m.dyn.norm	<-	drm(DIP/DIP[1] ~ conc, fct=LL.4(), na.action='na.omit')
	dfgi	<-	apply(dtp[dtp$time %in% 1:120,-1],1,
		FUN=function(x)	as.numeric(getGInci(t0val=dtp[1,2],tnval=x, cnval=x[1])))
	plot(1, xlim=c(1e-11,1e-4),ylim=c(-1.2,1.2), type='n', xlab=NA, ylab=NA, log='x', xaxt='n', yaxt='n')
	axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), padj=-.5, labels=c(-11,-9,-7,-5))
	axis(side=2, at=seq(from=-1, to=1, by=0.5), padj=.5)
	# NOTE: plotting log-logistic model fit to the simulated outputs here, not the actual sim outputs.
	if(actual.DTP)	
	{
	# plot the actual NCI DTP values
		for(ti in 12:120)
			lines(drug.conc, dfgi[,ti], col=timeCols[ti], lwd=ifelse(ti==72,3,1))
	} else {
	# plot log-logistic curve fits to the data
		for(ti in 12:120)
			plot(drm(dfgi[,ti] ~ drug.conc, fct=LL.4(), na.action='na.omit'), 
				type='none', add=TRUE, col=timeCols[ti], lwd=ifelse(ti==72,3,1))
	}
	plot(m.dyn.norm, type='none', add=TRUE, lwd=2)
	abline(h=0, lwd=1.5, lty=2, col=grey(.5))
	legend('bottomleft', bty='n', legend='DIP rate', col='black', cex=0.75, lty=1, lwd=3)
	abline(v=6.31e-7, col='orange', lwd=2)
	# add points at times matching points shown in growth curve above
	points(rep(6.31e-7,length(seq(24,120,24))), dfgi[match(6.31e-7,drug.conc),seq(24,120,24)])
	mtext('[drug], log10 M', side=1, line=2, font=2)
	mtext('Response ratio', side=2, line=2, font=2)
	li2 <- as.raster(matrix(timeCols[12:120], nrow=1))
	grid.raster(x=0.41, y=0.18, li2, height=.0075, width=.3, just=c('center','top'), interpolate=FALSE)
	text(2e-11,-0.85, paste(seq(24,120,24),collapse="  "), cex=0.65, pos=4)
	if(toFile)	dev.off()
}	

#########################################################################
# 						  
# 							FIGURE 4
# 						  
#########################################################################


PEAdata		<-	read.csv('../data for figs/BrCa_PEA_data.csv', as.is=TRUE)
timeCols	<-	topo.colors(120)
timeCols	<-	sub(topo.colors(120)[72], 'red', timeCols)

makeFig4	<-	function(d=PEAdata, toFile=FALSE,scale.bars=TRUE)
{
	id2plot	<-	unique(d$UID)
	h		<-	5.5
	nr		<-	2
	if(!toFile)	dev.new(width=6,height=h)
	if(toFile)	pdf(file='../R-generated figs/PEA data fig.pdf', width=6,height=h)

	par(mfrow=c(nr,3), mar=c(1,3,2,1), oma=c(5,1,2,1))

	for(uid in rev(id2plot))
	{
		dtp	<-	d[d$UID==uid,]
		cl	<-	unique(dtp$cellLine)
		drug	<-	c('phenformin','rotenone','rotenone','2DG')[match(unique(dtp$drug),c("phen","rot","rotenone","2DG"))]
	# aggregating data by summing nuclei counts in each condition and time point
		a	<-	aggregate(NumNuclei ~ time + factor(conc), data=dtp, FUN=sum)
		colnames(a)[colnames(a)=='factor(conc)']	<-	'conc'
		a$nl2	<-	log2norm(a$NumNuclei, norm_id=a$time, ids=a$conc)
		a$conc	<-	as.numeric(as.character(a$conc))
		conc	<-	sort(unique(a$conc))
		conc.range	<-	range(conc[conc!=0])
		conc.col	<-	c('#000000',colorpanel(n=length(conc)-1,low='blue',mid=grey(0.5),high='red'))

	# PLOT 1: Growth curves of experimental data (log scale)
	# make basic plot for pop growth curves
		plot(a$time,a$nl2, xlim=c(0,100), ylim=c(-1,5), xlab=NA, ylab=NA, main=NA, type='n', xaxt='n', yaxt='n')
		axis(side=1, at=c(0:5)*24, padj=-.5)
		axis(side=2, at=c(-1:5), padj=.5)
		if(uid==tail(id2plot,1))		mtext('Time (h)', side=1, line=1.5, font=2, cex=0.75)
		mtext('Population doublings', side=2, line=1.5, font=2, cex=0.75)
		if(grepl('rot',uid) & scale.bars)
		{
			text(0,4.1,'[rotenone], log10 M', font=2, pos=4, cex=0.75)
			text(0,4.8,paste(seq(-11,-5,2),collapse='    '), pos=4, cex=0.75)
		}
		if(grepl('phen',uid) & scale.bars)
		{
			text(0,4.1,'[phenformin], log10 M', font=2, pos=4, cex=0.75)
			text(0,4.8,paste(seq(-8,-2,2),collapse='    '), pos=4, cex=0.75)
		}
		for(c in unique(a$conc))
		{
			dfl	<-	a[a$conc==c,]
			lines(dfl$time,dfl$nl2, col=conc.col[match(c,conc)], lwd=2)
		}
		abline(h=0, lwd=1.5, lty=2, col=grey(.5))
		
		
	# aggregating data by averaging the nuclei counts in each condition and time point
	# must average or the variation of cells per condition is too high when normalizing to control
		a2	<-	aggregate(NumNuclei ~ time + factor(conc), data=dtp, FUN=mean)
		dfsdrc	<-	data.frame(t(sapply(1:9,FUN=function(x) a2$NumNuclei[seq(x, length(a2$NumNuclei), 9)])))
		colnames(dfsdrc)	<-	as.character(conc)
		rownames(dfsdrc)	<-	unique(a2$time)
		
		xl	<-	c(10^floor(log10(conc.range)[1]),10^ceiling(log10(conc.range)[2]))

	# PLOT 2: static vs dynamic DRC
	# make basic plot for DRC (static metric)
		plot(1, xlim=xl,ylim=c(-0.2,1.2), type='n', xlab=NA, ylab=NA, log='x', xaxt='n', yaxt='n')
		axis(side=1, at=10^seq(from=log10(xl[1]), to=log10(xl[2])), padj=-.5, labels=seq(from=log10(xl[1]), to=log10(xl[2])))
		axis(side=2, at=seq(from=-0.2, to=1.2, by=0.2), padj=.5)

		static.effect	<-	t(apply(dfsdrc[-1,],1,FUN=function(x) x/mean(as.numeric(dfsdrc[1,]))))	# response ratio at each time point after 0
		
		# dynamic response ratio DRC (add to static plot)
		# converting DIP rates to log2 scale (doublilngs per hour)
		linMod		<-	coef(lm(nl2 ~ time * factor(conc), data=a))
		linMod		<-	linMod[grep('time',names(linMod))]
		linMod		<-	c(linMod[1],linMod[-1]+linMod[1])
		dyn.effect	<-	linMod/linMod[1]
		m.dyn.norm	<-	drm(linMod/linMod[1] ~ conc, fct=LL.4(), na.action='na.omit')

		# dynamic response ratio DRC (add to static plot)
		# converting DIP rates to log2 scale (doublilngs per hour)
		DIP		<-	coef(lm(nl2 ~ time * factor(conc), data=a[a$time>24,]))
		DIP		<-	DIP[grep('time',names(DIP))]
		DIP		<-	c(DIP[1],DIP[-1]+DIP[1])
		DIP.ratio	<-	DIP/DIP[1]
		
		m.DIP.ratio	<-	drm(DIP/DIP[1] ~ conc, fct=LL.4(), na.action='na.omit')
		for(ti in unique(a2$time)[-1])
		{
			m.static	<-	drm(as.numeric(dfsdrc[match(ti,rownames(dfsdrc)),])/as.numeric(dfsdrc[match(ti,rownames(dfsdrc)),])[1] ~ conc, fct=LL.4(), na.action='na.omit')
			plot(m.static, type='none', add=TRUE, col=timeCols[floor(ti)], lwd=2)
			if(floor(ti)==72) cat(paste(uid,'EC50=',signif(ED(m.static, 50, interval='delta',display=FALSE),3)[1],'\n'))
		}
		plot(m.dyn.norm, type='none', add=TRUE, lwd=2)
		cat(paste(uid,'DIP EC50=',signif(ED(m.dyn.norm, 50, interval='delta',display=FALSE)[1],3),'\n'))

		abline(h=0, lwd=1.5, lty=2, col=grey(.5))
		if(uid==tail(id2plot,1))		mtext('[drug], log10 M', side=1, font=2, line=2, cex=0.75)
		mtext('Response ratio', side=2, line=1.5, font=2, cex=0.75)
		text(xl[2],1.2,'Static effect',pos=2,font=2)		
		
	# PLOT 3: NCI DTP dynamic DRC versis DIP-based DRC
	# make basic plot for DRC (NCI dynamic metric)
		plot(1, xlim=xl,ylim=c(-0.2,1.2), type='n', xlab=NA, ylab=NA, log='x', xaxt='n', yaxt='n')
		axis(side=1, at=10^seq(from=log10(xl[1]), to=log10(xl[2])), padj=-.5, labels=seq(from=log10(xl[1]), to=log10(xl[2])))
		axis(side=2, at=seq(from=-0.2, to=1.2, by=0.2), padj=.5)
		
	# obtain DRC using NCI DTP approach (dynamic metric but uses cell number vs log cell number)
		t0	<-	as.numeric(dfsdrc[1,])
		c0	<-	t0[1]
		cn	<-	dfsdrc[,1]
		refTable	<-	as.matrix(dfsdrc)
		dimnames(refTable)	<-	NULL
		whichCol	<-	function(x) which(apply(refTable,2,identical,x))
		whichRow	<-	function(x) which(apply(refTable,1,identical,x))
		NCI.DTP		<-	t(apply(dfsdrc[-1,],1,FUN=function(x) getGInci(t0val=t0, tnval=x, c0val=c0, cnval=cn[whichRow(as.numeric(x))])))

		for(ti in unique(a2$time)[-1])
		{
			m.NCI.DTP	<-	drm(NCI.DTP[match(ti,rownames(NCI.DTP)),] ~ conc, fct=LL.4(), na.action='na.omit')
			plot(m.NCI.DTP, type='none', add=TRUE, col=timeCols[floor(ti)], lwd=2)
		}
		plot(m.dyn.norm, type='none', add=TRUE, lwd=2)
		abline(h=0, lwd=1.5, lty=2, col=grey(.5))
		if(uid==tail(id2plot,1))		mtext('[drug], log10 M', side=1, font=2, line=2, cex=0.75)
		mtext('Response ratio', side=2, line=1.5, font=2, cex=0.75)
		text(xl[2],1.2,'NCI DTP',pos=2,font=2)		
	}

	# color legend for drug conc
	li	<- as.raster(matrix(conc.col, nrow=1))
	# for rotenone growth curve
	if(scale.bars)	grid.raster(x=0.11, y=0.87, li, height=.01, width=.1, just=c('left','top'), interpolate=FALSE)
	# for phenformin growth curve
	if(scale.bars)	grid.raster(x=0.11, y=0.455, li, height=.01, width=.1, just=c('left','top'), interpolate=FALSE)

	# color legend for times
	li2 <- as.raster(matrix(timeCols[12:120], nrow=1))
	grid.raster(x=0.495, y=0.055, li2, height=.01, width=.1, just=c('center','top'), interpolate=FALSE)
	mtext(paste(seq(24,120,24),collapse=" "), outer=TRUE, side=1, line=3, cex=0.5)

	if(toFile)	dev.off()
}


#########################################################################
# 						  
# 							FIGURE 5
# 						  
#########################################################################


# must ensure sufficient digits are available
options(digits=12)
PC9			<-	read.csv('../data for figs/PC9DoseResp.csv', as.is=TRUE)
DS.corr		<-	read.csv('../data for figs/DS_corr_vals.csv', as.is=TRUE)
allDRM		<-	getDRM()
post72		<-	read.csv('../data for figs/post72hCounts.csv', as.is=TRUE)
DS.dip		<-	assembleDIP(post72)
# add DIP rate values into DS.corr data.frame
DS.corr[DS.corr$var=='DIP',][sapply(DS.dip$all.rates$UID, FUN=function(x) 
	grep(x, DS.corr[DS.corr$var=='DIP','ID'])),'val']	<- DS.dip$all.rates$DIP
DS345_3erl		<-	read.csv('../data for figs/DS345 growth curve data.csv', as.is=TRUE)
DS345_erlDR	<-	read.csv('../data for figs/DS345+erl in HTS core.csv',as.is=TRUE)
colnames(DS345_erlDR)	<-	c('plate','ID','Subline','conc','Time_h','Cell.count')

brCa		<-	read.csv('../data for figs/BrCa+EGFRi_data.csv', as.is=TRUE)

brCa_static	<-	read.csv('../data for figs/BrCa+EGFRi_static_data.csv', as.is=TRUE)
# estimated cell numbers after 72 hours of drug treatment
aggDynData	<-	function(d)	{
	out	<-	data.frame()
	for(tx in unique(d$Treatment))
	{
		dta	<-	d[d$Conc==0 | d$Treatment==tx,]
		a	<-	aggregate(nl2 ~ Conc + Time, data=dta, FUN=mean)
		a	<-	cbind(Treatment=tx,a)
		a	<-	a[order(a$Conc,a$Time),]
		ifelse(nrow(out)==0, out <- a, out <- rbind(out,a))
	}
	rownames(out)	<-	NULL
	out$Treatment	<-	as.character(out$Treatment)
	out
}

brCa_dyn	<-	aggDynData(d=brCa)
maxConcData	<-	subset(brCa_dyn, (Treatment=='BIBW2992' & Conc %in% c(0, 2e-6)) | 
	(Treatment %in% c('Erlotinib','Lapatinib') & Conc == 8e-6))
maxConcData[maxConcData$Conc==0,'Treatment']	<-	'control'

stats	<-	generateStats(DS.corr)

makeFig5a	<-	function(toFile=FALSE)
{
	if(toFile)	pdf('DS growth curves.pdf', width=4.5, height=3)
	if(!toFile)	dev.new(width=4.5, height=3)
	par(mfrow=c(1,2), mar=c(5,4,1,0), cex=0.5, cex.axis=2, oma=c(1,1,1,5))
	plot(nl2 ~ Time_h, data=DS345_3erl[DS345_3erl$Subline=='PC9',], xlab=NA, ylab=NA, xlim=c(0,150), 
		ylim=c(-0.5,5), axes=FALSE, frame.plot=TRUE)
	axis(side=1, at=c((0:2)*72), padj=0.5)
	axis(side=2, at=0:5)
	mtext(side=2, 'Population doublings', font=2, line=3)
	abline(h=0, col=grey(0.2),lty=2)
	for(ds in c('DS3','DS4','DS5'))
	{
		points(DS345_3erl[DS345_3erl$Subline==ds,'Time_h'],DS345_3erl[DS345_3erl$Subline==ds,'nl2'], 
			col=c('green','red','orange')[match(ds,c('DS3','DS4','DS5'))])
		dip.curve	<-	coef(lm(nl2 ~ Time_h, data=DS345_3erl[DS345_3erl$Subline==ds & DS345_3erl$Time_h >= 72,]))
		curve(x*dip.curve[2]+dip.curve[1], from=72, to=150, lwd=4, add=TRUE, 
			col=c('green','red','orange')[match(ds,c('DS3','DS4','DS5'))])
	}
	ctrl.lm.coef	<-	coef(lm(nl2 ~ Time_h, data=DS345_3erl[DS345_3erl$Subline=='PC9',]))
	curve(x*ctrl.lm.coef[2]+ctrl.lm.coef[1], from=0, to=140, lwd=4, add=TRUE)
	
	p72	<-	DS345_3erl[DS345_3erl$Time >= 72,]
	plot(nl2 ~ Time_h, data=p72, xlab=NA, ylab=NA, xlim=c(72,150), ylim=c(-1,0.5), type='n', 
		axes=FALSE, frame.plot=TRUE)
	axis(side=1, at=c((3:6)*24), padj=0.5)
	axis(side=4, at=c(-1,-0.5,0,0.5), padj=0.5)
	abline(h=0, col=grey(0.2),lty=2)
	for(ds in c('DS3','DS4','DS5'))
	{
		points(p72[p72$Subline==ds,'Time_h'],p72[p72$Subline==ds,'nl2']-p72[p72$Subline==ds,'nl2'][1], 
			col=c('green','red','orange')[match(ds,c('DS3','DS4','DS5'))])
		dip.curve2	<-	coef(lm(nl2-nl2[1] ~ Time_h, data=p72[p72$Subline==ds,]))
		curve(x*dip.curve2[2]+dip.curve2[1], from=72, to=150, lwd=4, add=TRUE, 
			col=c('green','red','orange')[match(ds,c('DS3','DS4','DS5'))])
	}
	mtext(side=4, 'Population doublings', font=2, line=3, padj=0.5)
	mtext(side=1, outer=TRUE,text='Time (h)', font=2, line=-1)
	if(toFile)	dev.off()
}

makeFig5b	<-	function(toFile=FALSE,fn=paste0('Fig 5b DS345 DIP DRC ',' (',Sys.Date(),')','.pdf'))
{
	if(toFile) 
	{
		pdf(file=fn,width=3,height=3)
	}	else {dev.new(width=3,height=3)}
	par(mar=c(4,3,1,1), oma=c(0,2,0,0))
	plot(1,1, ylim=c(-0.25, 1), xlim=c(1e-9,1e-5),log='x', type='n', xlab=NA, ylab=NA, axes=FALSE, frame.plot=TRUE)
	axis(side=1, at=10^(-9:-5), labels=(-9:-5), padj=0)
	axis(side=2, at=seq(from=-0.2, to=1.2, by=0.2), padj=.5)
	mtext('[erlotinib], log10 M', side=1, font=2, line=2)
	mtext('DIP rate', side=2, font=2, line=3)
	mtext('response ratio', side=2, font=2, line=2)
	a	<-	aggregate(DS345_erlDR[,'Cell.count'], by=list(time=DS345_erlDR$Time_h, cellLine=DS345_erlDR$Subline, conc=DS345_erlDR$conc), FUN=sum)
	a$ID	<-	paste(a$cellLine,a$conc,sep='_')
	a$nl2	<-	log2norm(a$x, a$ID)
	a$drug	<-	'erl'
	a$conc	<-	a$conc*1e-6
	for(cl in c('DS4','DS5','DS3'))
	{
		m	<-	drcDIP(dtf=a, cl=cl, drg='erl', PIP.range=c(64.4,140))	
		plot(m, add=TRUE, pch=NA, lwd=2, col=c('red','orange','green')[match(cl,c('DS4','DS5','DS3'))])
	}
	abline(h=0, lty=2, col=gray(.2))
	if(toFile)	dev.off()
}

makeFig5c	<-	function(toFile=FALSE)
{
	if(toFile)	pdf('10d corr.pdf', width=6, height=3)
	if(!toFile)	dev.new(width=6, height=3)
	par(mfrow=c(1,2), mar=c(3,3,1,1))
	for(ty in c('DIP','norm72'))
	{
		if(ty=='DIP')
		{
			xl	<-	c(-0.02,0.01)
			add.line	<-	TRUE
			xlab.txt	<-	bquote(bold(DIP~rate~','~doublings~h^-1))
			ylab.txt	<-	bquote(bold('10day fold change'))
		}	else	{
	
			xl	<-	c(.25,.75)
			add.line	<-	FALSE
			xlab.txt	<-	'Relative cell number @ 72 h'
			ylab.txt	<-	NA
		}
		scatterPlotErr(stats,ty,'fold10d', col='grey',sfrac=0.01,
			gap=0.3,xlim=xl,add.line=add.line,xlab=xlab.txt,ylab=ylab.txt)
	}
	if(toFile)	dev.off()
}

makeFig5d	<-	function(myData=maxConcData[maxConcData$Time >= 48,], toFile=FALSE, xl=c(48,144), yl=c(0,2))	{
	if(!toFile)	dev.new(width=3, height=3)
	if(toFile)	pdf(file='HCC1954 max drug growth curves.pdf', width=3, height=3)
	par(font.lab=2, cex.lab=1.2, mar=c(3,3,1,1))
	plot(	myData$Time-48, myData$nl2, type='n',
			xlim=xl, ylim=yl,
				xlab=NA, 
				ylab=NA,
				main=NA,
				xaxt='n',
				yaxt='n'
	)
	axis(side=1, at=c(0:5)*24, padj=-.5)
	axis(side=2, at=c(0:4), padj=.5)
	for(tx in unique(myData$Treatment))	{
		dtp	<-	myData[myData$Treatment==tx,]
		lines(dtp$Time,dtp$nl2-dtp$nl2[1], lwd=3, col=c('black','blue','green','red')[match(tx, unique(myData$Treatment))])
		lab	<-	c('control','afatinib','lapatinib','erlotinib')[match(tx,unique(myData$Treatment))]
		text(	96, dtp[dtp$Time==96,]$nl2-dtp[dtp$Time==48,]$nl2, 
				lab, pos=4, cex=0.85)
		DIP	<-	round(coef(lm(nl2 ~ Time, data=dtp))['Time'],4)
	}
	mtext('Time (h)', side=1, line=1.5, font=2)
	mtext('Population doublings', side=2, line=1.5, font=2)
	if(toFile) dev.off()
}

makeFig5e	<-	function(a=brCa_dyn, toFile=FALSE)	{
	if(!toFile)	dev.new(width=7.5, height=2.5)
	if(toFile)	pdf(file='../R-generated figs/HCC1954 conc effect curves w DIP.pdf', width=7.5, height=2.5)
	par(mfrow=c(1,3),font.lab=2, mar=c(3,1,1,1), oma=c(1,3,0,0))
	for(tx in sort(unique(a$Treatment)))	{
		dtp		<-	a[a$Time>=48 & (a$Treatment==tx | a$Conc==0),]
		dip		<-	coef(lm(nl2 ~ Time * factor(Conc), data=dtp))
		dip		<-	dip[grep('Time',names(dip))]
		dip		<-	c(dip[1],(dip[-1] + dip[1]))
		myConc	<-	sort(unique(dtp$Conc))
		m		<-	drm(dip ~ myConc, fct=LL.4(), na.action=na.omit)
		ec50	<-	signif(ED(m, 50, interval = "delta", display=FALSE)[1],3)
		emax	<-	signif(min(dip),3)
		plot(m, ylim=c(-.02,0.04), xlim=c(1e-11,1e-5), xlab=NA, ylab=NA, type='confidence', axes=FALSE)
		plot(m, type='obs', add=TRUE)
		axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), labels=c(-11,-9,-7,-5), padj=-0.1,cex.axis=1.25)
		lab	<-	sub('bibw2992','afatinib',tolower(tx))
		mtext(paste('[',lab,'], log10 M', sep=""), side=1, line=2, font=2)
		if(tx==sort(unique(a$Treatment))[1])
		{
			mtext(bquote(bold(DIP~'rate,'~doublings~h^-1)), side=2, line=2)
			axis(side=2, at=c(seq(-4,6,2)*.01), labels=c(-0.04,-0.02,0,0.02,0.04,NA),cex.axis=1.25)
		}
		abline(h=0)
		abline(h=emax, col='grey', lty=4, lwd=2)
		abline(v=ec50, col='red', lty=2, lwd=2)
		text(ec50,-0.01,'EC50 =',pos=2)
		text(ec50,-0.015,paste(round(ec50*1e9,0),'nM'),pos=2)
		text(3e-10,emax-0.0005,paste('Emax =',emax),pos=1)
	}
	if(toFile) dev.off()
}

makeFig5f	<-	function(a=brCa_static, toFile=FALSE)	{
	if(!toFile)	dev.new(width=7.5, height=2.5)
	if(toFile)	pdf(file='../R-generated figs/HCC1954 standard DRC.pdf', width=7.5, height=2.5)
	par(mfrow=c(1,3),font.lab=2, mar=c(3,1,1,1), oma=c(1,3,0,0))
	for(tx in sort(unique(a$Treatment)))	{
		dtp	<-	a[a$Treatment==tx,]
		plot(	0, type='n', 
				xlim=c(1e-11,1e-4), ylim=c(0,1.2), log='x',
				xlab=NA, ylab=NA, xaxt='n', yaxt='n')
		points(dtp$Conc,dtp$respRatio)
		m	<-	drm(respRatio ~ Conc, data=dtp, fct=LL.4())
		ec50	<-	signif(ED(m, 50, interval = "delta", display=FALSE)[1],3)
		emax	<-	signif(coef(m)['c:(Intercept)'],3)
		plot(m, type='confidence', add=TRUE)
		axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), labels=c(-11,-9,-7,-5), padj=-0.1,cex.axis=1.25)
		abline(h=emax, col='grey', lty=4, lwd=2)
		abline(v=ec50, col='red', lty=2, lwd=2)
		text(ec50,1.15,'EC50 =',pos=4)
		text(ec50,1.05,paste(round(ec50*1e9,0),'nM'),pos=4)
		text(1e-10,emax,paste('Emax =',emax),pos=1)
		lab	<-	sub('bibw2992','afatinib',tolower(tx))
		mtext(paste('[',lab,'], log10 M', sep=""), side=1, line=2, font=2)
		if(tx==sort(unique(a$Treatment))[1])
		{
			mtext('Response ratio', font=2, side=2, line=2)
			axis(side=2, at=c(seq(0,1.2,0.2)),cex.axis=1.25)
		}
	}
	if(toFile) dev.off()
}


#########################################################################
# 						  
# 							OUTPUT FIGURES
# 						  
#########################################################################

makeFig2()
makeFig3()
makeFig4()
makeFig5a()
makeFig5b()
makeFig5c()
makeFig5d()
makeFig5e()
makeFig5f()
