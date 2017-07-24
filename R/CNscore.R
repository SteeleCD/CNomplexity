# function to score severity of copy number changes
CNscore = function(segFiles,threshLow=0.275,
                      threshHigh=0.8,doRel=TRUE,
                      doCat=FALSE,doScale=FALSE,
		      searchPattern=NULL,
			tumTotCol=7,normTotCol=5,
			startCol=3,endCol=4,chromCol=2,
			sampleCol=1,head=TRUE)
	{
	# seg info
	if(length(segFiles)==1)
		{
		seg = readFile(segFiles,head=head)
		samples = unique(seg[,sampleCol])
		} else {
		samples = segFiles
		}
	outScores = NULL
	# loop over samples
	for(i in 1:length(samples))
		{
		print(samples[i])
		# read seg file
		if(length(segFiles)>0)
			{
			data = seg[which(seg[,sampleCol]==samples[i]),]
			} else {
			data = readFile(samples[i],head=head)
			}
		# copy number ratios
		ratios = data[,tumTotCol]/data[,normTotCol] # copy number ratios
		# scores for each segment
		scores = vector(length=length(ratios))
		scores[which(ratios>=threshHigh)] = 2
		scores[which(ratios>=threshLow&ratios<threshHigh)] = 1
		scores[which(ratios<=c(-threshHigh))] = -2
		scores[which(ratios<=c(-threshLow)&ratios>c(-threshHigh))] = -1
		scores[which(is.na(scores))] = 0
		# length of segment as percentage of chms
		relLengths = apply(data,MARGIN=1,FUN=function(x)
			{
			(as.numeric(x[endCol])-as.numeric(x[startCol]))/
				(max(data[which(data[,chromCol]==x[chromCol]),endCol])-
					min(data[which(data[,chromCol]==x[chromCol]),startCol]))
			})
		# get categories
		if(doCat)
			{
			# get chromosome arm limits
			allArmInfo = getArmInfo()
			# get categories
			categories = unlist(sapply(unique(data[,chromCol]),
			   FUN=function(x) categoriseCN(allArmInfo[,x],
			               data[which(data[,chromCol]==x),,drop=FALSE])))
			# get into right format
			chromIndex = which(categories=="Chrom")
			focalIndex = which(categories=="Focal")
			armIndex = which(categories=="Arm")
			splitIndex = which(categories=="Large Split")
			outScores = rbind(outScores,
				c(sum(abs(scores[chromIndex])),
				sum(abs(scores[armIndex])),
				sum(abs(scores[focalIndex])),
				sum(abs(scores[splitIndex]))))
			} else {
			# otherwise score*relLength
			if(doRel) scores = scores*relLengths
			outScores = c(outScores,sum(abs(scores)))
			}  
		}
	# scale
	if(doScale) outScores = scale(outScores)
	# set names
	if(doCat)
		{
		rownames(outScores) = samples
		colnames(outScores) = c("Chrom","Arm","Focal","SplitLarge")
		} else {
		names(outScores) = samples
		}
	# return
	return(outScores)
	}

# get arm level chromosome lengths
getArmInfo = function(cytoFile=NULL)
	{
	data = loadCytoBand(cytoFile)
	data[,1] = gsub("chr","",data[,1])
	sapply(c(1:22,"X","Y"),FUN=function(x) 
		sapply(c("p","q"),FUN=function(y)
			c(min(data[which(data[,1]==x&grepl(y,data[,4])),2]),
			max(data[which(data[,1]==x&grepl(y,data[,4])),3]))))
	}

# categorise CN changes into focal/chromosome/arm/large split
categoriseCN = function(armInfo,CNinfo,startCol=3,endCol=4,chromCol=2)
	{
	armLengths = abs(diff(armInfo))[-2]
	chromLength = armInfo[4]-armInfo[1]
	lengths = CNinfo[,endCol]-CNinfo[,startCol]
	pFlag = !CNinfo[,endCol]<armInfo[2]
	qFlag = !CNinfo[,startCol]>armInfo[3]
	relArmLengths = sapply(1:length(lengths),
	         FUN=function(x) lengths[x]/armLengths[pFlag[x]+1])
	relChromLengths = lengths/chromLength
	mapply(FUN=decision,armLength=relArmLengths,
		chromLength=relChromLengths,
		bothArms=pFlag==qFlag)
	}

# categorise single CN segment
decision = function(armLength,chromLength,bothArms)
	{
	if(armLength<0.5)
		{
		return("Focal")
		} else {
		if(chromLength>=0.8)
			{
			if(bothArms)
				{
				return("Chrom")
				} else {
				return("Arm")
				}
			} else {
			if(bothArms)
				{
				return("Large split")
				} else {
				return("Arm")
				}
			}
		}
	}

# purity/ploidy adjusted CN threshold 
getCNthresh = function(purity=0.7,ploidy=2.89,ratio=1.27)
	{
	constant = 2*(1-purity)
	num = (purity*(ploidy*ratio))+constant
	den = (purity*ploidy)+constant
	return(log2(num/den))
	}

