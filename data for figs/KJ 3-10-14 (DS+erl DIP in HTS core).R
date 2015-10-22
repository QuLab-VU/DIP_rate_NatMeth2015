# reorganize raw data from Katie's HTS core PC9 DS + erl experiment

getCellLine		<-	function(myData)	{
	map		<-	read.csv('KJ 3-10-14 map.csv')
	id		<-	paste(map$plate,map$row,sep='_')
	map[match(substr(myData$ID, 1, 3), id),'cellLine']
}

d <- read.csv('KJ 3-10-14 analysis.csv')
d <- d[,1:9]
colnames(d)[colnames(d)=='Total.Cells..MultiWaveScoring.']	<-	'NumNuclei'
colnames(d)[colnames(d)=='Positive.W2..MultiWaveScoring.']	<-	'NumFUCCIpos'
colnames(d)[colnames(d)=='Well.Name']						<-	'Well'

conc		<-	c(0,round(3*4^(1:-6),4),0)
keepCol		<-	c('Plate.ID','Well','Site.ID','NumNuclei','NumFUCCIpos')

d1	<-	d[,keepCol]
plateIDs		<-	unique(d1$Plate.ID)
d1$myPlateID	<-	(match(d1$Plate.ID, plateIDs)-1) %% 3 + 1
d1$UID	<-	paste(d1$myPlateID,d1$Well,sep='_')

a	<-	aggregate(d1[,c('NumNuclei','NumFUCCIpos')], by=list(plate=d1$Plate.ID, ID=d1$UID), FUN=sum)
a$cellLine	<-	getCellLine(a)
a$conc		<-	conc[as.integer(as.character(substr(a$ID,4,5)))-1]
a$time		<-	rep(c(0,24,(0:16*4+48)))				# Guessing about when drug was added

a			<-	a[,c('plate','ID','cellLine','conc','time','NumNuclei')] 

# write.csv(a[a$conc != 12,], file='DS1-9+erl in HTS core.csv')
# write.csv(a[a$cellLine %in% c('DS3','DS4','DS5'),], file='DS345+erl in HTS core.csv')