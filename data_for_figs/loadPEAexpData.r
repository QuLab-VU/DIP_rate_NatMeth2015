loadPEAexpData	<-	function(data.fol='../data for figs/')
{
	start.dir	<-	getwd()
	setwd(data.fol)
	source('loadMetabolBCdata.R')

	getMapInfo	<-	function(mapName)
	{
		library(gdata)
		ifelse(grepl('.xl',mapName),
			mapColNames	<-	as.character(read.xls(mapName, nrow=1, header=FALSE, as.is=TRUE)),
			mapColNames	<-	as.character(read.csv(mapName, nrow=1, header=FALSE, as.is=TRUE)))
	
		mapColNames	<-	gsub('\xb5','micro',mapColNames)
		mapColNames	<-	gsub('TR:', "", mapColNames)
	
		ifelse(grepl('.xl',mapName),
			map	<-	read.xls(mapName, head=FALSE, skip=1, as.is=TRUE),
			map	<-	read.csv(mapName, head=FALSE, skip=1, as.is=TRUE))
		colnames(map)	<-	mapColNames
		map
	}

	addMapInfo		<-	function(dfa,mapFile,drugName='PLX')
	{
		# dfa = data for annotation
		# may need to modifiy to find well column regardless of case (Well or well)
		map	<-	getMapInfo(mapFile)
		dfa$exptDate		<-	map[match(dfa$Well,map$Well),'Date']
		dfa$cellLine	<-	map[match(dfa$Well,map$Well),'Description']
		dfa$drug		<-	drugName
		dfa$conc		<-	0
		dfa$conc		<-	map[match(dfa$Well,map$Well),grep(drugName,colnames(map))]
		dfa
	}

	convertDrugConc	<-	function(dtc, drugCol='drug', concCol='conc', conc.str='mM')
	{
		# obtain correction factor based on concentration string (conc.str)
		corrFactor	<-	switch(conc.str,
			mM = 1e-3,
			uM = 1e-6,
			µM = 1e-6,
			nM = 1e-9,
			pM = 1e-12
		)
	
		# find rows with drug name containing '.' and conc.str (default='mM')
		dfc	<-	dtc[sapply(strsplit(dtc[,drugCol],'[\\.]'), FUN=function(x) any(x==conc.str)),]
	
		dfc[sapply(strsplit(dfc[,drugCol],'[\\.]'), FUN=function(x) any(x==conc.str)),concCol]	<-	dfc[sapply(strsplit(dfc[,drugCol],'[\\.]'), FUN=function(x) any(x==conc.str)),concCol]*corrFactor
		dfc[sapply(strsplit(dfc[,drugCol],'[\\.]'), FUN=function(x) any(x==conc.str)),drugCol]	<-	sapply(strsplit(dfc[sapply(strsplit(dfc[,drugCol],'[\\.]'), FUN=function(x) any(x==conc.str)),drugCol],'[\\.]'), FUN=function(x) x[1])
		dtc[which(sapply(strsplit(dtc[,drugCol],'[\\.]'), FUN=function(x) any(x==conc.str))),]	<-	dfc
		dtc
	}	
	

	bca	<-	loadMetabolBCdata()			# 20130901 MDAMB231 + Phen/Rot/Oligo
	bca	<-	bca[bca$drug %in% c('rot','phen'),]
	
	mel	<-	getCVCounts(read.csv('../data for figs/20130729 SKMel28 Rot_NACPhen_NACRot_CSVdata.csv'))
	mel	<-	addMapInfo(mel,'../data for figs/20130729 SKMEL28 Rot Phen plate map.xlsx','rotenone.nM')
	mel	<-	mel[mel$Row %in% c('B','C'),]
	rownames(mel)	<-	NULL
	colnames(mel)[match('Cell.Count', colnames(mel))]	<-	'NumNuclei'
	colnames(mel)[match('Time', colnames(mel))]	<-	'time'
	
	pc9	<-	getCVCounts(read.csv('../data for figs/20131011 PC9clone10 2DG_oligo_17AAG_CSVData.csv'))
	pc9	<-	pc9[pc9$Row %in% c('B','C') | pc9$Column == 2,]
	pc9	<-	addMapInfo(pc9,'../data for figs/20131011 PC9clone10 2DG_oligo_17AAG plate map.xlsx','2DG.mM')
	colnames(pc9)[colnames(pc9)=='Time']	<-	'time'
	colnames(pc9)[colnames(pc9)=='Cell.Count']	<-	'NumNuclei'
	
	
	setwd(start.dir)
	
# put all data in list or data.frame named <out>
	out		<-	data.frame()
	for(i in 1:3)	out		<-	rbind(out,list(bca,mel,pc9)[[i]])
	
	out	<-	convertDrugConc(out)
	out	<-	convertDrugConc(out,conc.str='nM')
	out$time	<-	floor(out$time)

	out$UID	<-	paste(out$exptDate,out$cellLine,out$drug, sep='_')
	out
	
}