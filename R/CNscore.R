# function to score severity of copy number changes
CNscore = function(dir,threshLow=0.275,
                      threshHigh=0.8,doRel=TRUE,
                      doCat=FALSE,doScale=FALSE,
		      searchPattern=NULL,
			tumTotCol=7,normTotCol=5,
			startCol=3,endCol=4,chromCol=2)
	{
	outScores = NULL
	files = list.files(dir)
	if(!is.null(searchPattern)) files = files[grep(searchPattern,files)]
	for(i in 1:length(files))
		{
		print(files[i])
		# read seg file
		data = read.csv(paste0(dir,"/",files[i],
		                       "/",files[i],".ascat_ngs.summary.csv")
		                ,head=FALSE,as.is=TRUE)
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
		rownames(outScores) = files
		colnames(outScores) = c("Chrom","Arm","Focal","SplitLarge")
		} else {
		names(outScores) = files
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

