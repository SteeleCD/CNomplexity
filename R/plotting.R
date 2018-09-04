# ========================================================================
#		PLOT A SEGMENT PROFILE
# ========================================================================

# plot a single chromosome
segChromPlot = function(seg.mean,num.mark=NULL,seg.start,seg.end,seg.extra=NULL,
		YLIM,colours,
		highlightStarts1,highlightEnds1,highlightCol1=rgb(1,0,0,0.5),
		highlightStarts2,highlightEnds2,highlightCol2=rgb(0,0,1,0.1),
		fusionPos,chrom,offset=0.0125)
	{
	seg.mean = as.numeric(seg.mean); SEG.MEAN<<-seg.mean
	if(!is.null(num.mark)) num.mark = as.numeric(num.mark); NUM.MARK<<-num.mark
	seg.start = as.numeric(seg.start); SEG.START<<-seg.start
	seg.end = as.numeric(seg.end); SEG.END<<-seg.end
	if(!is.null(seg.extra)) seg.extra = as.numeric(seg.extra); SEG.EXTRA<<-seg.extra
	COLOURS<<-colours
	# empty plot
	plot(NA,ylim=c(0,max(seg.mean,na.rm=TRUE)),xlim=range(c(seg.start,seg.end),na.rm=TRUE),col=colours[num.mark],xaxt="n",main=chrom)
	abline(h=0)
	abline(h=1,lty=2)
	# offset if two alleles	
	offset = offset*max(seg.mean)
	if(!is.null(seg.extra)) seg.mean=seg.mean+offset
	# plot CN segments
	if(!is.null(num.mark))
		{
		# colour segments by number of markers
		sapply(1:length(seg.mean),FUN=function(x) lines(x=c(seg.start[x],seg.end[x]),y=rep(seg.mean[x],2),col=colours[num.mark][x],lwd=3))
		} else {
		# single colour for all segments
		sapply(1:length(seg.mean),FUN=function(x) lines(x=c(seg.start[x],seg.end[x]),y=rep(seg.mean[x],2),col=colours,lwd=3))
		}
	if(!is.null(seg.extra)) sapply(1:length(seg.extra),FUN=function(x) lines(x=c(seg.start[x],seg.end[x]),y=rep(seg.extra[x]-offset,2),col="green",lwd=3))
	# higlight regions of interest
	if(length(highlightStarts1)>0)
		{
		for(i in 1:length(highlightStarts1))
			{
			polygon(x=c(highlightStarts1[i],
					highlightEnds1[i],
					highlightEnds1[i],
					highlightStarts1[i]),
				y=c(YLIM[2]+1,YLIM[2]+1,YLIM[1]-1,YLIM[1]-1),
				col=highlightCol1,border=NA)
			}
		}
# higlight regions of interest
	if(length(highlightStarts2)>0)
		{
		for(i in 1:length(highlightStarts2))
			{
			polygon(x=c(highlightStarts2[i],
					highlightEnds2[i],
					highlightEnds2[i],
					highlightStarts2[i]),
				y=c(YLIM[2]+1,YLIM[2]+1,YLIM[1]-1,YLIM[1]-1),
				col=highlightCol2,border=NA)
			}
		}
	# plot vertical lines at fusions
	if(length(fusionPos)>0)
		{
		abline(v=fusionPos,lty=2)
		}
	}

# plot a segment profile
segPlot = function(segFile=NULL,segObj=NULL,dataDir,
		outDir=NULL,fileName=NULL,
		sampleCol=2,chromCol=1,startCol=3,endCol=4,nMarkCol=NULL,segMeanCol=6,segBCol=NULL,
		highlightChroms1=NULL,highlightStarts1=NULL,highlightEnds1=NULL,
		highlightChroms2=NULL,highlightStarts2=NULL,highlightEnds2=NULL,
		fusionChroms=NULL,fusionPos=NULL,segColour=NULL)
	{
	# checking/getting seg profile
	if(is.null(segFile)&is.null(segObj)) stop("Must specify a segment profile")
	if(!is.null(segFile)) data = read.table(paste0(dataDir,"/",segFile),head=TRUE,sep="\t")
	if(!is.null(segObj)) data = segObj
	# plot by chms
	chms = unique(data[,chromCol])
	individual = unique(data[,sampleCol])
	# plot colours
	if(!is.null(nMarkCol)) 
		{
		colours = colorRampPalette(c("red","green"))(max(data[,nMarkCol]))
		} else {
		colours=ifelse(is.null(segColour),"black",segColour)
		}
	# plot
	if(!is.null(fileName)) pdf(paste0(outDir,'/',fileName))
	# loop over individuals
	sapply(individual,FUN=function(x) {
		# set plotting parameters
		par(mfrow=c(4,6),mar=c(1,2,2,0))
		# loop over chromosomes
		sapply(chms,FUN=function(y) {
			print(paste0(x,":",y))
			# any highlighted regions in this chms?
			if(!is.null(highlightChroms1))
				{
				highlightIndex1 = which(highlightChroms1==y)
				} else {
				highlightIndex1 = c()
				} 
			if(!is.null(highlightChroms2))
				{
				highlightIndex2 = which(highlightChroms2==y)
				} else {
				highlightIndex2 = c()
				} 
			# any fusions in this chms?
			if(!is.null(fusionChroms))
				{
				fusionIndex = which(fusionChroms==y)
				} else {
				fusionIndex = c()
				}
			# which rows are this sample and this chms
			index = which(data[,sampleCol]==x&data[,chromCol]==y)
			# optional inputs
			if(is.null(nMarkCol))
				{
				nMark=NULL
				} else {
				nMark=data[index,nMarkCol]
				}
			if(is.null(segBCol))
				{
				segB=NULL
				} else {
				segB=data[index,segBCol]
				}
			# plotting
			# coloured by number of markers
			segChromPlot(seg.mean=data[index,segMeanCol],
				num.mark=nMark,
				seg.start=data[index,startCol],
				seg.end=data[index,endCol],
				seg.extra=segB,
				YLIM=c(0,max(data[,segMeanCol])),
				colours=colours,
				highlightStarts1=highlightStarts1[highlightIndex1],
				highlightEnds1=highlightEnds1[highlightIndex1],
				highlightStarts2=highlightStarts2[highlightIndex2],
				highlightEnds2=highlightEnds2[highlightIndex2],
				fusionPos=fusionPos[fusionIndex],
				chrom=y)

			})
		})
	if(!is.null(fileName)) dev.off()
	}
