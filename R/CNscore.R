CNscore = function(dir,threshLow=0.275,
                      threshHigh=0.8,doRel=TRUE,
                      doCat=FALSE,doScale=FALSE)
	{
	outScores = NULL
	files = list.files(dir)
	files = files[grep("a",files)]
	for(i in 1:length(files))
		{
		print(files[i])
		data = read.csv(paste0(dir,"/",files[i],
		                       "/",files[i],".ascat_ngs.summary.csv")
		                ,head=FALSE,as.is=TRUE)
		# copy number ratios
		ratios = data[,7]/data[,5] # copy number ratios
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
			#x = as.numeric(x)
			(as.numeric(x[4])-as.numeric(x[3]))/
				(max(data[which(data[,2]==x[2]),4])-
					min(data[which(data[,2]==x[2]),3]))
			})
		if(doCat)
			{
			# get categories
			allArmInfo = getArmInfo()
			categories = unlist(sapply(unique(data[,2]),
			   FUN=function(x) categoriseCN(allArmInfo[,x],
			               data[which(data[,2]==x),,drop=FALSE])))
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
			if(doRel) scores = scores*relLengths
			outScores = c(outScores,sum(abs(scores)))
			}  
		}
	# scale
	if(doScale) outScores = scale(outScores)
	if(doCat)
		{
		rownames(outScores) = files
		colnames(outScores) = c("Chrom","Arm","Focal","SplitLarge")
		} else {
		names(outScores) = files
		}
	return(outScores)
	}

getArmInfo = function(cytoFile="/nfs/users/nfs_c/cs32/cytoBand.txt")
	{
	data = read.table(cytoFile,head=FALSE)
	data[,1] = gsub("chr","",data[,1])
	sapply(c(1:22,"X","Y"),FUN=function(x) 
		sapply(c("p","q"),FUN=function(y)
			c(min(data[which(data[,1]==x&grepl(y,data[,4])),2]),
			max(data[which(data[,1]==x&grepl(y,data[,4])),3]))))
	}

categoriseCN = function(armInfo,CNinfo)
	{
	armLengths = abs(diff(armInfo))[-2]
	chromLength = armInfo[4]-armInfo[1]
	lengths = CNinfo[,4]-CNinfo[,3]
	pFlag = !CNinfo[,4]<armInfo[2]
	qFlag = !CNinfo[,3]>armInfo[3]
	relArmLengths = sapply(1:length(lengths),
	         FUN=function(x) lengths[x]/armLengths[pFlag[x]+1])
	relChromLengths = lengths/chromLength
	mapply(FUN=decision,armLength=relArmLengths,
		chromLength=relChromLengths,
		bothArms=pFlag==qFlag)
	}

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

