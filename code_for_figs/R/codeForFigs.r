# R code used to generate figures in Harris et al., Nature Methods 2016
# primary code written by Darren Tyson
#
#
source('NatMethFxns.r')
source('../../dipDRC/dipDRC.r')
output.to.file = FALSE

 
#########################################################################
# ___  ___      _        ______ _                           
# |  \/  |     (_)       |  ___(_)                          
# | .  . | __ _ _ _ __   | |_   _  __ _ _   _ _ __ ___  ___ 
# | |\/| |/ _` | | '_ \  |  _| | |/ _` | | | | '__/ _ \/ __|
# | |  | | (_| | | | | | | |   | | (_| | |_| | | |  __/\__ \
# \_|  |_/\__,_|_|_| |_| \_|   |_|\__, |\__,_|_|  \___||___/
#                                 __/ /                    
#                                |___/                     
#
#        Thanks to: http://www.kammerl.de/ascii/AsciiSignature.php
#########################################################################


#########################################################################
# 						  
# 					FUNCTIONS USED TO PRODUCE FIGURES
# 						  
#########################################################################

#########################################################################
# 						FUNCTION FOR FIGURE 1
#########################################################################

makeFig1	<-	function(toFile=FALSE)
{
	if(toFile)	pdf(file='Theor time bias.pdf', width=10, height=5)
	if(!toFile)	dev.new(width=10, height=5)
	par(mfrow=c(3,6), mar=c(0.5,3,1,1), oma=c(6,1,2,1))
	for(i in colnames(param))
	{
		dtp	<-	cell.counts[[i]]
		
		# linear growth curve
		plot(1, xlim=c(0,120),ylim=c(25,10000), type='n', xlab=NA, ylab=NA, xaxt='n', yaxt='n', bty='n')
		axis(side=1, at=c(0,24,72,120), padj=-.5, labels=FALSE)
		axis(side=2, at=seq(from=0, to=8, by=4)*1000, labels=c(0,'4,000','8,000'), padj=.5)
		for(c in colnames(dtp)[-1])
			lines(dtp$time,dtp[,c],col=drugCols[match(c,drug.conc)], type='l')
		if(i==colnames(param)[2])
			mtext(side=2, 'Number of cells', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3])
		{
			axis(side=1, at=c(0:5)*24, padj=-.5, tick=FALSE)
			lines(0:120,dtp[,'6.31e-07'], col='orange', lwd=3)
			mtext('Time (h)', side=1, font=2, line=2, cex=0.75)
		}
		
		# log growth curve
		plot(1, xlim=c(0,120),ylim=c(-4,6), type='n', xlab=NA, ylab=NA, xaxt='n', yaxt='n', bty='n')
		axis(side=1, at=c(0,24,72,120), padj=-.5, labels=FALSE)
		axis(side=2, at=seq(from=-4, to=6, by=2), padj=.5)
		for(c in colnames(dtp)[-1])
			lines(dtp$time,log(dtp[,c])-log(dtp[,c])[1],col=drugCols[match(c,drug.conc)], type='l')
		if(i==colnames(param)[2])
			mtext(side=2, 'Population doublings', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3]) 
		{
			lines(0:120,log(dtp[,'6.31e-07'])-log(dtp[,'6.31e-07'])[1], col='orange', lwd=3)
			mtext('Time (h)', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(0,24,72,120), padj=-.5, tick=FALSE)
		}
		
		
		# static response ratio DRC
		dfsdrc	<-	apply(dtp[-1], 1, FUN=function(x) x/x[1])
		plot(1, xlim=c(1e-11,1e-4),ylim=c(-0.6,1.2), type='n', xlab=NA, ylab=NA, log='x', xaxt='n', yaxt='n', bty='n')
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
		plot(m.dyn.norm, type='none', add=TRUE, lwd=2)
		abline(h=0, col=grey(.5))
		if(i==colnames(param)[2])
			mtext(side=2, 'Response ratio', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3]) 	
		{
			mtext('[drug], log10 M', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), labels=c(-11,-9,-7,-5), padj=-.5, tick=FALSE)
		}
		
		# Time bias on EC50
		plot(12:120, ec50.static, xlim=c(0,120), ylim=c(1e-9,1e-5), type='n', xlab=NA, ylab=NA, log='y', xaxt='n', yaxt='n', bty='n')
		points(12:120, ec50.static, col=timeCols[12:120], pch=19, cex=0.5)
		axis(side=1, at=c(0,24,72,120), padj=-.5, labels=FALSE)
		axis(side=2, at=c(10^c(-9,-7,-5)), labels=c(-9,-7,-5), padj=.5)
		points(stable.range[1],signif(ED(m.dyn.norm, 50, display=FALSE),3)[1], pch=2, cex=1.5)
		points(stable.range[-1],rep(signif(ED(m.dyn.norm, 50, display=FALSE),3)[1], length(stable.range[-1])), pch=3, cex=.25)
		
		if(i==colnames(param)[2])
			mtext(side=2, expression(bold('EC'[50]*', log10 M')), line=2, font=2, cex=0.75)
		if(i==colnames(param)[3])
		{
			mtext('Time (h)', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(0,24,72,120), padj=-.5, tick=FALSE)
		}
		
		# Time bias on AA
		plot(12:120,aa.static, type='l', xlim=c(0,120), ylim=c(0,1), xlab=NA, ylab=NA, xaxt='n', yaxt='n', bty='n')
		points(12:120, aa.static, col=timeCols[12:120], pch=19, cex=0.5)
		axis(side=1, at=c(0,24,72,120), padj=-.5, labels=FALSE)
		axis(side=2, at=seq(0,1,0.5), padj=.5)	
		points(stable.range[1],getAA(m.dyn.norm), pch=2, cex=1.5)
		points(stable.range[-1],rep(getAA(m.dyn.norm), length(stable.range[-1])), pch=3, cex=.25)
		if(i==colnames(param)[2])
			mtext(side=2, 'Activity area, AU', line=2, font=2, cex=0.75)
		if(i==colnames(param)[3])
		{
			mtext('Time (h)', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(0,24,72,120), padj=-.5, tick=FALSE)	
		}
		
		# DIP DRC
		plot(1, xlim=c(1e-11,1e-4),ylim=c(-0.036,0.072), type='n', xlab=NA, ylab=NA, log='x', xaxt='n', yaxt='n', bty='n')
		m.dyn	<-	drm(DIP ~ drug.conc, fct=LL.4(), na.action='na.omit')
		plot(m.dyn, type='none', add=TRUE, lwd=2)
		axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), padj=-.5, labels=FALSE)
		axis(side=2, at=c(seq(-4,6,4)*.01), labels=c(-0.04,0,0.04), padj=.5)
		abline(h=0, col=grey(.5))
		if(i==colnames(param)[2])
			mtext(side=2, expression(bold('DIP rate, doublings h'^-1)), line=2, font=2, cex=0.75)
		if(i==colnames(param)[3])
		{
			mtext('[drug], log10 M', side=1, font=2, line=2, cex=0.75)
			axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), labels=c(-11,-9,-7,-5), padj=-.5, tick=FALSE)
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

#########################################################################
# 						FUNCTIONS FOR FIGURE 2
#########################################################################

makeFig2ab	<-	function(d=PEAdata, toFile=FALSE,scale.bars=TRUE)
{
	id2plot	<-	unique(d$UID)
	h		<-	6
	nr		<-	2
	if(!toFile)	dev.new(width=4.5,height=h)
	# PEA = partial equilibrium assumption
	if(toFile)	pdf(file='PEA data fig.pdf', width=4.5,height=h)

	par(mfcol=c(nr,2), mar=c(1,3,2,1), oma=c(5,1,2,1))

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
		plot(a$time,a$nl2, xlim=c(0,100), ylim=c(-1,5), xlab=NA, ylab=NA, main=NA, type='n', xaxt='n', yaxt='n', bty='n')
		axis(side=1, at=c(0:5)*24, padj=-.5)
		axis(side=2, at=c(-1:5), padj=.5)
		mtext('Time (h)', side=1, line=1.5, font=2)
		mtext('Population doublings', side=2, line=2, font=2)
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
		plot(1, xlim=xl,ylim=c(-0.2,1.2), type='n', xlab=NA, ylab=NA, log='x', xaxt='n', yaxt='n', bty='n')
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
		abline(h=0, lwd=1.5, lty=2, col=grey(.5))
		mtext('[Drug], log10 M', side=1, font=2, line=2)
		mtext('Response ratio', side=2, line=2, font=2)
		
	}

	# color legend for drug conc
	li	<- as.raster(matrix(conc.col, nrow=1))
	# for rotenone growth curve
	if(scale.bars)	grid.raster(x=0.155, y=0.86, li, height=.01, width=.2, just=c('left','top'), interpolate=FALSE)
	# for phenformin growth curve
	if(scale.bars)	grid.raster(x=0.625, y=0.86, li, height=.01, width=.2, just=c('left','top'), interpolate=FALSE)

	# color legend for times
	li2 <- as.raster(matrix(timeCols[12:120], nrow=1))
	grid.raster(x=0.48, y=0.065, li2, height=.01, width=.19, just=c('center','top'), interpolate=FALSE)
	mtext(paste(seq(24,120,24),collapse=" "), outer=TRUE, side=1, line=3, cex=0.75)

	if(toFile)	dev.off()
}


makeFig2c	<-	function(toFile=FALSE)
{
	if(toFile)	pdf('DS growth curves.pdf', width=2.5, height=6)
	if(!toFile)	dev.new(width=2.5, height=6)
	par(mfrow=c(2,1), mar=c(1,3,2,1), oma=c(5,1,2,0),cex.axis=.9)
	plot(nl2 ~ Time_h, data=DS345_3erl[DS345_3erl$Subline=='PC9',], xlab=NA, ylab=NA, xlim=c(0,150), 
		ylim=c(0,5), axes=FALSE, frame.plot=FALSE, type='l')
	axis(side=1, at=c((0:2)*72), padj=-.75)
	axis(side=2, at=0:5, padj=.75)
	mtext(side=2, 'Population doublings', font=2, line=2)
	mtext(side=1, 'Time (h)', font=2, line=1.5)
	for(ds in c('DS3','DS4','DS5'))
	{
		lines(DS345_3erl[DS345_3erl$Subline==ds,'Time_h'],DS345_3erl[DS345_3erl$Subline==ds,'nl2'], 
			col=c('green','red','orange')[match(ds,c('DS3','DS4','DS5'))], lwd=2)
		dip.curve	<-	coef(lm(nl2 ~ Time_h, data=DS345_3erl[DS345_3erl$Subline==ds & DS345_3erl$Time_h >= 72,]))
		curve(x*dip.curve[2]+dip.curve[1], from=72, to=150, lwd=4, add=TRUE, 
			col=c('green','red','orange')[match(ds,c('DS3','DS4','DS5'))])
	}
	ctrl.lm.coef	<-	coef(lm(nl2 ~ Time_h, data=DS345_3erl[DS345_3erl$Subline=='PC9',]))
	curve(x*ctrl.lm.coef[2]+ctrl.lm.coef[1], from=0, to=140, lwd=4, add=TRUE)
	legend('topleft', legend=c('PC9 parental','DS4 + 3µM erl','DS5 + 3µM erl','DS3 + 3µM erl'), 
		col=c('black','red','orange','green'), lwd=4, bty='n', cex=0.5)
	
	plot(1,1, ylim=c(-0.25, 1), xlim=c(1e-9,1e-5),log='x', type='n', xlab=NA, ylab=NA, axes=FALSE, frame.plot=FALSE)
	axis(side=1, at=10^(-9:-5), labels=(-9:-5), padj=-0.5)
	axis(side=2, at=seq(from=-0.2, to=1.2, by=0.2), padj=.5)
	mtext('[erlotinib], log10 M', side=1, font=2, line=2)
	mtext('DIP rate', side=2, font=2, line=3)
	mtext('response ratio', side=2, font=2, line=2)
	a	<-	aggregate(DS345_erlDR[,'Cell.count'], by=list(rel.time=DS345_erlDR$Time_h, 
			cellLine=DS345_erlDR$Subline, conc=DS345_erlDR$conc), FUN=sum)
	a$ID	<-	paste(a$cellLine,a$conc,sep='_')
	a$nl2	<-	log2norm(a$x, a$ID)
	a$drug	<-	'erl'
	a$conc	<-	a$conc*1e-6
	a$expt.date <- 1
	
	m	<-	dipDRC(a,xName='rel.time',yName='x',var=c('cellLine','drug','conc','expt.date'), add=TRUE, norm=TRUE, 
			showEC50=FALSE, type='none')	
	for(i in names(m)) plot(m[[i]],type='none',col=c('green','red','orange')[match(i,names(m))], add=TRUE, lwd=2)
	if(toFile)	dev.off()
	invisible(m)
}


makeFig2d	<-	function(myData=maxConcData, toFile=FALSE, xl=c(0,144), yl=c(0,4))	{
	if(!toFile)	dev.new(width=3, height=3)
	if(toFile)	pdf(file='HCC1954 max drug growth curves.pdf', width=3, height=3)
	par(font.lab=2, cex.lab=1.2, mar=c(3,3,1,1))
	plot(	myData$Time, myData$nl2, type='n', bty='n',
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
		lines(dtp$Time,dtp$nl2-dtp$nl2[1], lwd=3, col=c('black','blue','orange')[match(tx, unique(myData$Treatment))])
		DIP	<-	round(coef(lm(nl2 ~ Time, data=dtp))['Time'],4)
	}
	mtext('Time (h)', side=1, line=1.5, font=2)
	mtext('Population doublings', side=2, line=1.5, font=2)
	legend('topleft', legend=c('control','8µM erl','8µM lap'), 
		col=c('black','blue','orange'), lwd=2, bty='n', cex=0.9)
	if(toFile) dev.off()
}

makeFig2e1	<-	function(a=brCa_dyn[brCa_dyn$Treatment%in%c('Erlotinib','Lapatinib'),], toFile=FALSE)	{
	if(!toFile)	dev.new(width=5, height=2.5)
	if(toFile)	pdf(file='HCC1954 DIP DRC.pdf', width=5, height=2.5)
	par(mfrow=c(1,2),font.lab=2, mar=c(3,1,1,1), oma=c(1,3,0,0))
	for(tx in sort(unique(a$Treatment)))	{
		dtp		<-	a[a$Time>=48 & (a$Treatment==tx | a$Conc==0),]
		dip		<-	coef(lm(nl2 ~ Time * factor(Conc), data=dtp))
		dip		<-	dip[grep('Time',names(dip))]
		dip		<-	c(dip[1],(dip[-1] + dip[1]))
		myConc	<-	sort(unique(dtp$Conc))
		m		<-	drm(dip ~ myConc, fct=LL.4(), na.action=na.omit)
		ec50	<-	signif(ED(m, 50, interval = "delta", display=FALSE)[1],3)
		emax	<-	signif(min(dip),3)
		plot(m, ylim=c(-.02,0.04), xlim=c(1e-11,1e-5), xlab=NA, ylab=NA, type='confidence', axes=FALSE, bty='n')
		plot(m, type='obs', add=TRUE, cex=0.75)
		axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), labels=c(-11,-9,-7,-5), padj=-0.1)
		lab	<-	sub('bibw2992','afatinib',tolower(tx))
		mtext(paste('[',lab,'], log10 M', sep=""), side=1, line=2, font=2)
		if(tx==sort(unique(a$Treatment))[1])
		{
			mtext(bquote(bold(DIP~'rate,'~doublings~h^-1)), side=2, line=2)
			axis(side=2, at=c(seq(-4,6,2)*.01), labels=c(-0.04,-0.02,0,0.02,0.04,NA))
		}
		abline(h=0)
		abline(v=ec50, col='red', lty=2, lwd=2)
		text(ec50,-0.008,'EC50 =',pos=2)
		text(ec50,-0.015,paste(round(ec50*1e9,0),'nM'),pos=2)
	}
	if(toFile) dev.off()
}

makeFig2e2	<-	function(a=brCa_static[brCa_static$Treatment!='BIBW2992',], toFile=FALSE)	{
	if(!toFile)	dev.new(width=5, height=2.5)
	if(toFile)	pdf(file='HCC1954 standard DRC.pdf', width=5, height=2.5)
	par(mfrow=c(1,2),font.lab=2, mar=c(3,1,1,1), oma=c(1,3,0,0))
	for(tx in sort(unique(a$Treatment)))	{
		dtp	<-	a[a$Treatment==tx,]
		plot(	0, type='n',  bty='n',
				xlim=c(1e-11,1e-5), ylim=c(0,1.2), log='x',
				xlab=NA, ylab=NA, xaxt='n', yaxt='n')
		points(dtp$Conc,dtp$respRatio, cex=0.75)
		m	<-	drm(respRatio ~ Conc, data=dtp, fct=LL.4())
		ec50	<-	signif(ED(m, 50, interval = "delta", display=FALSE)[1],3)
		emax	<-	signif(coef(m)['c:(Intercept)'],3)
		plot(m, type='confidence', add=TRUE)
		axis(side=1, at=c(1e-11,1e-9,1e-7,1e-5), labels=c(-11,-9,-7,-5), padj=-0.1)
		abline(v=ec50, col='red', lty=2, lwd=2)
		abline(h=0)
		text(ec50,0.3,'EC50 =',pos=2)
		text(ec50,0.15,paste(round(ec50*1e9,0),'nM'),pos=2)
		lab	<-	sub('bibw2992','afatinib',tolower(tx))
		mtext(paste('[',lab,'], log10 M', sep=""), side=1, line=2, font=2)
		if(tx==sort(unique(a$Treatment))[1])
		{
			mtext('Response ratio', font=2, side=2, line=2)
			axis(side=2, at=c(seq(0,1.2,0.2)))
		} else {
			axis(side=2, at=c(seq(0,1.2,0.2)), labels=NA)

		}
	}
	if(toFile) dev.off()
}

#########################################################################
# 					FUNCTION FOR FIGURE 3
#########################################################################

makeFig3	<-	function(toFile=FALSE, bias='IC50', after=48)
{
	if(!toFile)	dev.new(width=6.5,height=6)
	if(toFile)	pdf(file='Melanoma + BRAFi time bias.pdf',width=6.5,height=6)
	par(mfcol=c(3,4), mar=c(3,2,.5,.5), oma=c(0,2,.5,0))

	cln			<-	unique(mel$cellLine)
	timeCols	<-	topo.colors(120)
	drug.conc	<-	unique(mel$Conc)
	drugCols	<-	colorpanel(n=length(drug.conc),low='blue',mid=grey(0.5),high='red')

	# growth curves
	for(cl in cln)
	{
		dtp			<-	mel[mel$cellLine==cl,]
		a	<-	aggregate(cell.count ~ rel.time + Conc, data=dtp, sum)
		a$nl2	<-	log2norm(a$cell.count,a$Conc)
		plot(a$rel.time, a$nl2, ylim=c(-0.5,5), type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, bty='n')
		abline(h=0, lwd=1.5, lty=2, col=grey(.5))
		for(co in unique(a$Conc))	lines(a[a$Conc==co,'rel.time'],a[a$Conc==co,'nl2'], col=drugCols[match(co,drug.conc)], lwd=2)
		axis(side=1, at=0:5*24, padj=-.5)
		ifelse(cl == cln[1], axis(side=2, at=0:5, padj=.5), axis(side=2, at=0:5, padj=.5, labels=NA))
		text(0,4.5, cl, font=2, pos=4, cex=1.5)
		if(cl == cln[1])	mtext('Population doublings', font=2, line=2, side=2)

	# dose-response curves
	
		m.dyn		<-	list()
		m.dyn.norm	<-	list()
		params		<-	data.frame()
		addCCLE		<-	FALSE
		conc		<-	unique(dtp$Conc)
		dip			<-	getDIP(dtp)
		resp.ratio	<-	dip/dip['0']
		ec50		<-	numeric()
		ic50		<-	numeric()
		aa			<-	numeric()
		m.dyn[[cl]]	<-	drm(dip ~ conc, fct=LL.4(), na.action='na.omit')
		m.dyn.norm[[cl]]	<-	drm(resp.ratio ~ conc, fct=LL.4(), na.action='na.omit')
		plot(m.dyn.norm[[cl]], type='none', ylim=c(0,1.25), 
			xtlab=c(0,-7,-6,-5), xlab=NA, ylab=NA, lwd=2, axes=FALSE, bty='n')			
		axis(side=1, at=c(0,1e-7,1e-6,1e-5), labels=c(0,-7,-6,-5), padj=-.5)
		if(cl == cln[1])	mtext('Response ratio', font=2, line=2, side=2)
		ifelse(cl == cln[1], axis(side=2, at=seq(from=0, to=1, by=0.5),  padj=.5), 
			axis(side=2, at=seq(from=0, to=1, by=0.5),  padj=.5, labels=NA))
		abline(h=0, lwd=1.5, lty=2, col=grey(.5))
		for(ti in unique(dtp$rel.time)[-1])
		{
			dfsdrc	<-	dtp[dtp$rel.time==ti,]
			if(length(unique(dfsdrc$Conc)) <= 2) next
			ctrl	<-	mean(dfsdrc[dfsdrc$Conc==0,'cell.count'])
			resp.ratio	<-	dfsdrc$cell.count/ctrl
			m	<-	drm(resp.ratio ~ dfsdrc$Conc, fct=LL.4())
			plot(m, add=TRUE, type='none', col=timeCols[ti], lwd=2)
			ec50	<-	append(ec50,coef(m)['e:(Intercept)'])
			# suppress errors in obtaining IC50 values when response is insufficient to achieve
			ic50	<-	append(ic50,tryCatch({ED(m,0.5,'delta',type="absolute", 
				display=FALSE)[1]},warning=function(cond){return(NA)}))
			aa		<-	append(aa,getAA(m))
		}
		plot(m.dyn.norm[[cl]], type='none', lwd=2, add=TRUE)			# 
		if(inCCLE(cl) & addCCLE)
		{
			m2	<-	getDRMfromCCLE(ccle.raw[grepl(cl,ccle.raw$CCLE.Cell.Line.Name),])
			plot(m2, add=TRUE, col='red',lwd=2)
		}
		times	<-	unique(dtp$rel.time)[-1]
		params	<-	rbind(params,data.frame(rel.time=times,EC50=ec50,IC50=ic50,AA=aa,cellLine=cl))


	# time-bias in IC50 
		plot(1, 1e-8, xlim=c(0,120), ylim=c(1e-8,1e-4), type='n', xlab=NA, ylab=NA, log='y', xaxt='n', yaxt='n', bty='n')
		stable.range	<-	seq(after,tail(params[params$cellLine==cl,'rel.time'],1))
		points(params[params$cellLine==cl,'rel.time'], params[params$cellLine==cl,bias], 
			col=timeCols[params[params$cellLine==cl,'rel.time']], pch=19, cex=2)
		if(cl %in% ccle$cellLine)	points(72,ccle[ccle$cellLine==cl,bias]*1e-6,pch=8, cex=1.5)
		if(cl %in% sData$cell_line_name)	points(72,sData[sData$cell_line_name==cl,bias],pch=12, cex=1.5)
		axis(side=1, at=c(0:5)*24, padj=-.5)
		ifelse(cl == 'A2058',
			axis(side=2, at=10^(-8:-4), labels=c(-8:-4), padj=.5),
			axis(side=2, at=10^(-8:-4), labels=FALSE, padj=.5))
		if(cl == 'A2058')
		{
			legend('bottomleft',c('DIP','static','CCLE','GDSC'),pch=c(95,19,8,12), col=c('black',timeCols[48],'black','black'), bty='n',border='n')
			mtext(side=2,expression(bold('PLX4720 IC'[50]~', log10 M')), line=2, font=2)
		}
		ic50	<-	ED(m.dyn.norm[[cl]],0.5,'delta',type="absolute", display=FALSE)[1]

		points(after,signif(ic50,3), pch=2, cex=1.5)
		points(stable.range[-1],rep(signif(ic50,3), length(stable.range[-1])), pch=3, cex=.25)
	}
	if(toFile)	dev.off()
}

#########################################################################
# 						  
# 					 LOAD/INITIALIZE REQUIRED DATA
# 						  
#########################################################################

#########################################################################
# 						LOAD DATA FOR FIGURE 1
#########################################################################
#########################################################################					  
# 						  Set parameter values 					  
#########################################################################
drug.conc		<-	signif(10^((-110:-40)/10),3)		# concentrations of drug (1e-11 to 1e-4 M)
times			<-	0:120								# hours
states			<-	c(Cell=100,CellPrime=0)				# initial conditions (all 100 cells are in Cell state)

k_on	<-	c(1e8,1e8,1e5)
k_off	<-	c(1,1,1e-3)
DIP0	<-	c(0.06,0.01,0.06)*log(2)					# parameters in base e to solve ODEs
DIPmax	<-	c(-0.03,-0.005,-0.03)*log(2)				
param	<-	rbind(k_on=k_on,k_off=k_off,DIP0=DIP0,DIPmax=DIPmax)
colnames(param)	<-	c('fast','slow','fast.ne')			# ne = non-equilibrium
param	<-	as.matrix(param)


drugCols	<-	colorpanel(n=length(drug.conc),low='blue',mid=grey(0.5),high='red')
timeCols	<-	topo.colors(120)
timeCols	<-	sub(topo.colors(120)[72], 'red', timeCols)

#########################################################################					  
# 	Generate cell counts for each cell type and drug concentration					  
#########################################################################

cell.counts	<-	list()
for(i in colnames(param))
	cell.counts[[i]]	<-	getCellCount(param[,i])		# call ODE solver with relevant parameters


#########################################################################
# 					LOAD DATA FOR FIGURE 2
#########################################################################
PEAdata		<-	read.csv('../../data_for_figs/BrCa_PEA_data.csv', as.is=TRUE)
timeCols	<-	topo.colors(120)
timeCols	<-	sub(topo.colors(120)[72], 'red', timeCols)

DS345_3erl		<-	read.csv('../../data_for_figs/DS345 growth curve data.csv', as.is=TRUE)
DS345_erlDR	<-	read.csv('../../data_for_figs/DS345+erl in HTS core.csv',as.is=TRUE)
colnames(DS345_erlDR)	<-	c('plate','ID','Subline','conc','Time_h','Cell.count')

brCa		<-	read.csv('../../data_for_figs/BrCa+EGFRi_data.csv', as.is=TRUE)

brCa_static	<-	read.csv('../../data_for_figs/BrCa+EGFRi_static_data.csv', as.is=TRUE)
# estimated cell numbers after 72 hours of drug treatment

brCa_dyn	<-	aggDynData(d=brCa)

maxConcData	<-	subset(brCa_dyn, (Treatment %in% c('Erlotinib','Lapatinib') & Conc==8e-6) | Conc==0)
maxConcData[maxConcData$Conc==0,'Treatment']	<-	'control'
maxConcData	<-	maxConcData[!duplicated(makeUCond(maxConcData,c('Treatment','Time'))),]



#########################################################################
# 					LOAD DATA FOR FIGURE 3
#########################################################################


mel	<-	read.csv('../../data_for_figs/mel+plx time course.csv')

# Data downloaded from CCLE: Raw data dose-response data generated by the GDSCP used to obtain model fits
# http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_27/CCLE_NP24.2009_Drug_data_2015.02.24.csv?downloadff=true&fileId=20777
# Data for PLX4720-treated melanoma cell lines was saved in separate file, maintaining exact data structure
ccle.raw	<-	read.csv('../../data_for_figs/CCLE raw data (skin+PLX4720).csv', as.is=TRUE)

# Data and dose-response curve parameters downloaded from CCLE
# http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_27/CCLE_GNF_data_090613.xls?downloadff=true&fileId=13199
# Data for PLX4720-treated melanoma cell lines was saved in separate file, maintaining exact data structure
ccle		<-	read.csv('../../data_for_figs/BRAF mut DRC param from CCLE.csv', as.is=TRUE)
ccle$cellLine	<-	gsub('-','',ccle$cellLine)

# Data and dose-response curve parameters downloaded from GDSCP: 
# ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_drug_sensitivity_fitted_data_w5.zip
# (VERY LARGE FILE! >135 MB)
# Data for PLX4720-treated melanoma cell lines was saved in separate file, maintaining exact data structure
sData	<-	read.csv('../../data_for_figs/Sanger PLX in melanoma cell lines.csv', as.is=TRUE)
sData	<-	sData[,c('cell_line_name','drug_id','max_conc','ic_50_est','alpha_est','beta_est','i_0_est','i_max_est','e_est')]
sData$drug_name	<-	'PLX4720'
sData$EC50	<-	exp(sData$e_est)*1e-6				# convert from µM to M
sData$IC50	<-	exp(sData$ic_50_est)*1e-6			# convert from µM to M


#########################################################################
# 						  
# 							OUTPUT MAIN FIGURES
# 						  
#########################################################################

makeFig1()
makeFig2ab()
makeFig2c()
makeFig2d()
makeFig2e1()
makeFig2e2()
makeFig3()



#########################################################################
#
#  _____                   _                           _                   
# /  ___|                 | |                         | |                  
# \ `--. _   _ _ __  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _ _ __ _   _ 
#  `--. \ | | | '_ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | '__| | | |
# /\__/ / |_| | |_) | |_) | |  __/ | | | | |  __/ | | | || (_| | |  | |_| |
# \____/ \__,_| .__/| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|_|   \__, |
#             | |   | |                                               __/ |
#             |_|   |_|                                              |___/ 
#
#        Thanks to: http://www.kammerl.de/ascii/AsciiSignature.php
#########################################################################


#########################################################################
# 						  
# 					FUNCTIONS FOR SUPPLEMENTARY FIGURES
# 						  
#########################################################################

#########################################################################
# FUNCTION FOR SUPPLEMENTARY FIGURE 4 - ILLUSTRATION OF THEORETICAL BIAS
#########################################################################

makeSuppFigTheorDelay	<-	function(toFile=FALSE, actual.DTP=TRUE)
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
# 					FUNCTION FOR SUPPLEMENTARY FIGURE 5
#########################################################################

# DRC of PC9 and sublines treated with erlotinib using 72 h cell count (static) metric

makeDRC	<-	function(mydrm=allDRM[c('DS3','DS4','DS5','PC9')], toFile=FALSE, w=8, h=2.5, ...)
{
	if(toFile) {pdf(file='PC9 subline DRC.pdf',width=w, height=h)} else {
		dev.new(width=w, height=h)}
	par(mfrow=c(1,4), mar=c(3.5,2,.5,.5), oma=c(0,2,0,0))
	for(n in names(mydrm)) 
	{
		plot(	0, type='n',  bty='n',
				xlim=c(1e-11,1e-5), ylim=c(0,1.2), log='x',
				xlab=NA, ylab=NA, xaxt='n', yaxt='n')
		plot(mydrm[[n]], add=TRUE, ...)
		plot(mydrm[[n]], type='confidence', add=TRUE, ...)
		axis(side=1, at=10^(seq(-12,-6,2)), labels=seq(-12,-6,2), cex.axis=1.25)
		legend('topright',n,pch='', bty='n', text.font=2, cex=1.5)
		if(n==names(mydrm)[1]) {
			mtext('Response ratio', side=2, line=2, font=2)
			axis(side=2, at=seq(from=0, to=1, by=0.2), cex.axis=1.25)
		} else {
			axis(side=2, at=seq(from=0, to=1, by=0.2), labels=NA)
		}	
		ec50	<-	signif(ED(mydrm[[n]], 50, interval = "delta", display=FALSE)[1],3)
		abline(v=ec50, col='red', lty=2, lwd=2)
		text(ec50,0.35,'EC50 =',pos=2, cex=1.5)
		text(ec50,0.25,paste(round(ec50*1e9,0),'nM'),pos=2, cex=1.5)
	}
	mtext('[erlotinib], log10 M', side=1, line=-1, font=2, outer=TRUE)
	if(toFile) dev.off()
}


#########################################################################
# FUNCTION FOR SUPPLEMENTARY FIGURE 6 - AUTOMATED ESTIMATION OF DIP RATE
#########################################################################
findDIPgraphs	<-	function(dat=ds.data, toFile=FALSE, met='ar2',...)
{
	m <- list()
	if(!toFile)	dev.new(width=6, height=6)
	if(toFile)	pdf(file='findDIPrate.pdf', width=6, height=6)
	par(mfrow=c(3,3), oma=c(0,1,0,0), mar=c(4,4,1,0.5))
	for(sl in unique(dat$Subline))
	{
		dtp <- dat[dat$Subline==sl,c('Time_h','l2')]
		m[[sl]] <-	plotGC_DIPfit(dtp=dtp,tit=sl, metric=met, newDev=FALSE, toFile=FALSE,...)
	}
	if(toFile)	dev.off()
	invisible(m)
}

#########################################################################
# FUNCTION FOR SUPPLEMENTARY FIGURE 7 - EFFECTS OF SAMPLING FREQUENCY
#########################################################################

makeSamplingFig	<-	function(toFile=FALSE)
{
	if(!toFile)	dev.new(width=7,height=12)
	if(toFile)	pdf(file='Sampling.pdf', width=7,height=12)
	par(mfrow=c(5,3), mar=c(3.5,3.5,1,1))
	out <- subsamp(m$DS3$data,metric='ar2',tit='DS3',toFile=FALSE, newDev=FALSE)
	if(toFile) dev.off()
	invisible(out)
}

#########################################################################
# FUNCTION FOR SUPPLEMENTARY FIGURE XX - CORRELATIONS WITH 10-DAY OUTCOME
#########################################################################

makeDIPcorrFig	<-	function(toFile=FALSE)
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

#########################################################################
# FUNCTION FOR SUPPLEMENTARY FIGURE 9 - EFFECTS OF CELL SEEDING DENSITY 
#########################################################################


makeDensDepGCfig	<-	function(toFile=FALSE)
{
	if(!toFile)	dev.new(width=10,height=2.5)
	if(toFile)	pdf(file='SeedDensGC.pdf', width=10,height=2.5)
	par(mfrow=c(1,6), oma=c(1,2,0,0), mar=c(3,2,1,0.5), cex=0.75)
	for(s in unique(p8$Seeding.density))
	{
		dtp <- p8[p8$Seeding.density==s,]
		plot(log2(Cell.count) ~ Time, main=unique(dtp$Seeding.density), data=dtp, xlab=NA,ylab=NA,
			ylim=c(5,11.5), type='n')
		for(w in unique(dtp$Well))	lines(dtp[dtp$Well==w,]$Time, log2(dtp[dtp$Well==w,]$Cell.count))
	}
	mtext('Time (h)', side=1, font=2, line=-.5, outer=TRUE)
	mtext('log2(cell number)', side=2, font=2, line=0.5, outer=TRUE)
	if(toFile)	dev.off()
}

makeDensDep_nl2Fig	<-	function(toFile=FALSE)
{
	mycol	<-	topo.colors(length(unique(seed.data$Seeding.density)))
	if(!toFile)		dev.new(width=4,height=4)
	if(toFile)	pdf(file='SeedDens_nl2.pdf', width=4,height=4)
	par(font.lab=2)
	plot(nl2 ~ Time, data=p8, type='n', ylim=c(-1,2), ylab='Population doublings')
	for(w in unique(p8$Well)) lines(p8[p8$Well==w,'Time'],p8[p8$Well==w,'nl2'], 
		col=mycol[match(p8[p8$Well==w,'Seeding.density'],unique(p8$Seeding.density))])
	legend('bottomright',legend=unique(p8$Seeding.density),lwd=1, col=mycol, cex=0.6, bty='n')
	if(toFile) dev.off()
}

makeDensDIPfig	<-	function(toFile=FALSE)
{
	if(!toFile)	dev.new(width=5, height=2.5)
	if(toFile)	pdf(file='SeedDensOnDIP.pdf', width=6, height=3)
	par(mar=c(4,4,1,1),font.lab=2)
	plot(DIP ~ Seed.dens.fac, data=p8.rates, ylab=NA, xlab=NA, ylim=c(0,0.03))
	mtext(expression(bold(DIP~rate~','~doublings~h^-1)), side=2, font=2, line=2)
	mtext("Cells per well", side=1, font=2, line=2)
	if(toFile) dev.off()
}


#########################################################################
# 						  
# 					LOAD DATA FOR SUPPLEMENTARY FIGURES
# 						  
#########################################################################

# PC9 subline static dose-response curves and 10-day treatment correlates
# must ensure sufficient digits are available
options(digits=12)
PC9			<-	read.csv('../../data_for_figs/PC9_72hCellCountDoseResp.csv', as.is=TRUE)
DS.corr		<-	read.csv('../../data_for_figs/DS_corr_vals.csv', as.is=TRUE)

# This generates dose-response curve models from 72 hour cell count measurements
# of PC9 parental and DS sublines and stores them as a list of drm objects 
# (see DRC library in R for more information)
allDRM		<-	getDRM()


post72		<-	read.csv('../../data_for_figs/post72hCounts.csv', as.is=TRUE)
DS.dip		<-	assembleDIP(post72)
# add DIP rate values into DS.corr data.frame
DS.corr[DS.corr$var=='DIP',][sapply(DS.dip$all.rates$UID, FUN=function(x) 
	grep(x, DS.corr[DS.corr$var=='DIP','ID'])),'val']	<- DS.dip$all.rates$DIP

										# statistical tests for correlations between various
stats	<-	generateStats(DS.corr)		# metrics and number of cells after 10 days
										# of drug treatment

# Used for findDIP functions
ds.data		<-	read.csv('../../data_for_figs/DS345 growth curve data.csv', as.is=TRUE)
ds.data		<- ds.data[ds.data$Subline!='PC9'& ds.data$Time_h<=120,]

seed.data	<-	read.csv('../../data_for_figs/CellSeedingData.csv')
rates <- getWellRates(seed.data, c(70,200))

p8 <- seed.data[seed.data$conc==8,]
p8.rates <- getWellRates(p8, c(70,200))
p8.rates$Seeding.density	<-	p8[p8$Time==0,'Seeding.density']
p8.rates$Seed.dens.fac <- factor(p8.rates$Seeding.density, levels=unique(p8.rates$Seeding.density))

# Statistical tests of seeding density effects on DIP rate
#m.p8 <- lm(DIP ~ factor(Seeding.density), data=p8.rates)
#summary(m.p8)
#anova(m.p8)
leveneTest(lm(DIP ~ Seed.dens.fac, data=p8.rates))

#########################################################################
# 						  
# 				    	OUTPUT SUPPLEMENTARY FIGURES
# 						  
#########################################################################

# Supp Fig 5 - Static DRC of PC9 sublines
makeDRC()

m	<-	findDIPgraphs(met='rmse',o=0.001,add.line.met='ar2',dat.type='l2')
m2	<-	makeSamplingFig()


makeDensDepGCfig()
makeDensDep_nl2Fig()
makeDensDIPfig()

