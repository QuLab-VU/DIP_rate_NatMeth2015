# DIP rate paper functions and all required libraries
# 

require(deSolve)		# ODE solver usedfor generating outputs from 2-state model
						# without partial equilibrium assumption
require(gplots)			# needed for colors
require(grid)
require(drc)			# for dose-response curves; version 2.5-12


twoStateModel	<-	function(times, states, parameters)
{
	with(as.list(c(states, parameters)),{
		# rates
		dCell		<-	(DIP0 - k_on * Drug) * Cell + (k_off * CellPrime)
		dCellPrime	<-	(DIPmax - k_off) * CellPrime + (k_on * Drug * Cell)
		
		# return list of rates
		list(c(dCell,dCellPrime))
	})
}

dipPEA		<-	function(k_on, k_off,  DIP0, DIPmax, drug)
# analytical solution for value of DIP rate under PEA
{
	( 1 / log(2) ) * ( k_off * DIP0 + (k_on * drug * DIPmax) ) / (k_on * drug + k_off )
}

biasF	<- function(kx,ky,t)	exp((kx-ky)*t)


#
# 4-parameter log-logistic function
#
#                           d - c           
# f(x)  =  c  +  ---------------------------
#                1 + exp(b(log(x) - log(e)))
#
#	b: Hill coefficient
#	c: lower limit (Emax)
#	d: upper limit (Emin)
#	e: EC50
#

getAA	<-	function(drmod, conc=drug.conc, max=1)
{
	# drmod is dose–response model (drm)
	# normalize to max of 0 by subtracting 1 from all values
	# and divide by number of measurements
	vals	<-	PR(drmod, conc)
	-sum(vals-1) / length(vals)
}

getGI	<-	function(t0val,tnval,cnval)
# t0val = starting cell #
# tnval = cell # at time n (in specific condition)
# cnval = control cell # at time n
{
	(tnval-t0val) / (cnval-t0val)
}

getGInci	<-	function(t0val,tnval,c0val=t0val,cnval)
# t0val = starting cell #
# tnval = cell # at time n (in specific condition)
# c0val = starting cell # of control
# cnval = control cell # at time n
{
	ifelse(	tnval > t0val, 
			(tnval-t0val) / (cnval-c0val),
			(tnval-t0val) / c0val
	)
}


getCellCount	<-	function(parameters)
{
	for(dc in drug.conc)
	{	
		p	<-	append(parameters,dc)
		names(p)[5]	<-	'Drug'
		temp		<- as.data.frame(ode(y = states, times = times, func = twoStateModel, parms = p))
		Cell.tot	<-	as.numeric(temp$Cell+temp$CellPrime)
		ifelse(dc==drug.conc[1], out <- data.frame(times,Cell.tot), out <- cbind(out,Cell.tot))
	}
	colnames(out)	<-	c('time',as.character(drug.conc))
	out
}


getDRM	<-	function(d=PC9)
{
	allDRM	<-	list()
	for(cl in unique(d$cellLine))	
	{
		m				<-	drm(resp.ratio ~ conc, data=d[d$cellLine==cl,], fct=LL.4())		# to see curves from each expt add: , ECID
		name			<-	ifelse(cl=='par','PC9',cl)
		allDRM[[name]]	<-	m
	}
	allDRM
}

assembleDIP	<-	function(dat)
{
	mean.rates	<-	data.frame(subline=sort(unique(dat$Subline)),DIP=0,DIP.dev=0)
	all.rates	<-	data.frame()
	for(sl in sort(unique(dat$Subline)))
	{
		dtf	<-	dat[dat$Subline==sl,]
# uncomment if you want to see the raw data to which the model is being fit
#		plot(log2(Cell.count) ~ Time_h, data=dtf[dtf$Time_h>72,], main=sl)
		m	<-	lm(log2(Cell.count) ~ Time_h * factor(UID), data=dtf[dtf$Time_h>72,])
		dip	<-	coef(m)[grep('Time_h',names(coef(m)))]
		dip	<-	c(dip[1],dip[-1]+dip[1])
		names(dip)	<-	gsub('Time_h:factor\\(UID\\)','',names(dip))
		names(dip)[1]	<-	setdiff(unique(dtf$UID),names(dip))
		all.rates	<-	rbind(all.rates,data.frame(Subline=sl,DIP=dip, UID=names(dip)))
		mean.val	<-	t.test(dip)[['estimate']]
		dev.val		<-	( t.test(dip)['conf.int'][[1]][2]-t.test(dip)['conf.int'][[1]][1] ) / 2
		mean.rates[mean.rates$subline==sl,2:3]	<-	c(mean.val,dev.val)
	}
	list(all.rates=all.rates, mean.rates=mean.rates)
}

generateStats	<-	function(dat=DS.corr)
{
	stats	<-	aggregate(val ~ subl + var, data=dat, mean)
	colnames(stats)[3]	<-	'mean'
	stats	<-	cbind(stats,aggregate(val ~ subl + var, data=dat, length)['val'])
	colnames(stats)[4]	<-	'n'

	temp	<-	as.data.frame(aggregate(val ~ subl + var, data=dat, FUN=function(x) as.numeric(t.test(x)$conf.int))[,3])
	dev		<-	(temp[,2]-temp[,1])/2
	stats	<-	cbind(stats,dev)

	# add EC50 values
	ec50.mean	<-	sapply(allDRM,FUN=function(x) ED(x,50,interval = "delta", display=FALSE)[1])
	ec50.dev	<-	(sapply(allDRM,FUN=function(x) ED(x,50,interval = "delta", display=FALSE)[3])-
		sapply(allDRM,FUN=function(x) ED(x,50,interval = "delta", display=FALSE)[2]))/2
	ec50.n		<-	sapply(allDRM,FUN=function(x) length(unique(x$origData$exptDate)))
	ec50		<-	cbind(mean=ec50.mean,dev=ec50.dev,n=ec50.n)
	ec50[,1:2]	<-	apply(ec50[,1:2], c(1,2), FUN=function(x) signif(x*1e9,3))		# convert to nM

	stats		<-	rbind(stats,data.frame(subl=rownames(ec50),var='EC50',ec50))
	stats		<-	stats[order(stats$var,stats$subl),]
	stats[,3:4]	<-	apply(stats[,3:4], c(1,2), FUN=function(x) signif(x, digits=3))
	stats		<-	stats[!stats$subl %in% c('DS8','DS9','PC9'),]
	stats$dev	<-	sapply(stats$dev,FUN=function(x) signif(x,3))
	rm(list=c('ec50.mean','ec50.dev','ec50.n'))

	type	<-	as.character(unique(stats$var))
	stats
}

scatterPlotErr	<-	function(dtp,x.var,y.var,add.line=FALSE,xlab.txt=x.var,ylab.txt=y.var,...)
{
	x		<-	dtp[dtp$var==x.var,]
	y		<-	dtp[dtp$var==y.var,]
			
	# plots with y-axis errors
	plotCI(
		x$mean,
		y$mean,
		uiw=y$dev,
		err='y',
		pch=NA,
		...
	)

	# plots with x-axis errors
	plotCI(
		x$mean,
		y$mean,
		uiw=x$dev,
		err='x',
		pch=NA,
		add=TRUE,
		...
	)
	points(x$mean,y$mean,pch=19)
	mtext(side=1, xlab.txt, font=2, line=2)
	mtext(side=2, ylab.txt, font=2, line=2)	
	my.fit<- lm(y$mean ~ x$mean)
	ar2	<-	round(summary(my.fit)$adj.r.squared,2)
	mtext(as.expression(bquote(Adj~R^2~'='~.(ar2))), side=1, line=-1.25, cex=0.75,adj=.95)
	if(add.line)	abline(coef(my.fit),lty=2,lwd=2)	
}

getDIP		<-	function(dtf)	{
	myMod	<-	lm(nl2 ~ time, data=dtf)
	coef(myMod)['time']
}

getDRM	<-	function(dat=PC9)
{
	allDRM	<-	list()
	for(cl in unique(dat$cellLine))	
	{
		m				<-	drm(resp.ratio ~ conc, data=dat[dat$cellLine==cl,], fct=LL.4())		# to see curves from each expt add: , ECID
		name			<-	ifelse(cl=='par','PC9',cl)
		allDRM[[name]]	<-	m
	}
	allDRM
}

# Function to extract DIP rate across multiple condititons 
# and calculate a 4-param logistic fit

drcDIP	<-	function(dtf,cl=cl,drg=drg,PIP.range=PIP.range,norm=T)
{	
	rates	<-	data.frame(cellLine=character(), conc=integer(),
		DIPrate=numeric(),drug=character())
	dtp		<-	dtf[dtf$cellLine==cl & dtf$drug==drg,]
	Uconc	<-	unique(dtp$conc)
	for(co in 2:length(Uconc))
	{
		rates	<-	rbind(rates, 
			data.frame(cellLine=cl,conc=Uconc[co],
			DIPrate=getDIP(dtp[dtp$conc==Uconc[co] & dtp$time>=PIP.range[1] & dtp$time<=PIP.range[2],]),drug=drg))
	}
	rates$normDIP	<-	rates$DIPrate/rates[rates$conc==min(rates$conc),'DIPrate']
	if(norm)
	{	
		drm(normDIP~conc,data=rates,fct=LL.4())
	} else
	{
		drm(DIPrate~conc,data=rates,fct=LL.4())	
	}
}

log2norm <- function(count, ids, norm_type=c('idx','ref')[1], 
	norm_id="0", norm_vals, norm_idx=1, zero=log2(0.999))
{
# this function will normalize in one of two distinct ways:
# 1) using the position of each vector specified by 'idx'
# 2) using the position identified from matching 'norm_id'
# to its position of the vector from 'norm_val'
# e.g. norm_id = 3 could refer to a vector 'day' containing a time variable
# made this modification on 2014-09-23, but should still work with older code
    
    # l2
    l2 <- log2(count)
		# finds time points with no cells (0), and replace it
		# with zero (log2(0.999)), so that the data can be displayed 
		# in log scale, yet easily found.
    l2[is.infinite(l2)] <- zero
	norm <- numeric()
	group <- as.character(unique(ids))

    if(norm_type=='idx')
    {
		for(i in group)
		{
			d <- l2[ids == i]
			norm <- append(norm, d - d[norm_idx])
		}
	
		norm
		
	}	else	{
	
		for(i in group)
		{
			norm_pos		<-	match(norm_id, as.character(norm_vals)[ids == i])
			d <- l2[ids == i]
			norm <- append(norm, d - d[norm_pos])
		}
		norm
	}

}
