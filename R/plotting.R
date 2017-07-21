# ========================================================================
#		PLOT A SEGMENT PROFILE
# ========================================================================

# plot a single chromosome
segChromPlot = function(seg.mean,num.mark=NULL,seg.start,seg.end,
		YLIM,colours,highlightStarts,highlightEnds,highlightCol=rgb(0.1,0.1,0.1,0.1),fusionPos,chrom)
	{
	# wmpty plot
	plot(NA,ylim=range(seg.mean),xlim=range(c(seg.start,seg.end)),col=colours[num.mark],xaxt="n",main=chrom)
	abline(h=0)
	# plot CN segments
	if(!is.null(num.mark))
		{
		# colour segments by number of markers
		sapply(1:length(seg.mean),FUN=function(x) lines(x=c(seg.start[x],seg.end[x]),y=rep(seg.mean[x],2),col=colours[num.mark][x],lwd=3))
		} else {
		# single colour for all segments
		sapply(1:length(seg.mean),FUN=function(x) lines(x=c(seg.start[x],seg.end[x]),y=rep(seg.mean[x],2),col=colours,lwd=3))
		}
	# higlight regions of interest
	if(length(highlightStarts)>0)
		{
		for(i in 1:length(highlightStarts))
			{
			polygon(x=c(highlightStarts[i],
					highlightEnds[i],
					highlightEnds[i],
					highlightStarts[i]),
				y=c(YLIM[2]+1,YLIM[2]+1,YLIM[1]-1,YLIM[1]-1),
				col=highlightCol)
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
		sampleCol=2,chromCol=1,startCol=3,endCol=4,nMarkCol=NULL,segMeanCol=6,
		highlightChroms=NULL,highlightStarts=NULL,highlightEnds=NULL,
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
			if(!is.null(highlightChroms))
				{
				highlightIndex = which(highlightChroms==y)
				} else {
				highlightIndex = c()
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
			# plotting
			if(!is.null(nMarkCol))
				{
				# coloured by number of markers
				segChromPlot(seg.mean=data[index,segMeanCol],
					num.mark=data[index,nMarkCol],
					seg.start=data[index,startCol],
					seg.end=data[index,endCol],
					YLIM=range(data[,segMeanCol]),
					colours=colours,
					highlightStarts=highlightStarts[highlightIndex],
					highlightEnds=highlightEnds[highlightIndex],
					fusionPos=fusionPos[fusionIndex],
					chrom=y)
				} else {
				# coloured black
				segChromPlot(seg.mean=data[index,segMeanCol],
					seg.start=data[index,startCol],
					seg.end=data[index,endCol],
					YLIM=range(data[,segMeanCol]),
					colours=colours,
					highlightStarts=highlightStarts[highlightIndex],
					highlightEnds=highlightEnds[highlightIndex],
					fusionPos=fusionPos[fusionIndex],
					chrom=y)
				}
			})
		})
	if(!is.null(fileName)) dev.off()
	}
