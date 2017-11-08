# ==================================================================================
#			DEFINE FUNCTIONS
# ==================================================================================
# setup windows
setupBFB = function(windowSize=10000000,windowGap=100000,runChrom="7")
	{
	print("setup")
	# get chromosome info
	chromInfo = CNomplexity:::getChromInfo()
	lims = chromInfo[,runChrom]
	# windows to run	
	library(GenomicRanges)
	end = lims[2]-windowSize
	if(end<0) end = lims[1]+windowGap
	windows = seq(from=lims[1],to=end,by=windowGap)
	windowGR = as(paste0(runChrom,":",windows,"-",windows+windowSize),"GRanges")
	return(windowGR)
	}

# get fold back inversion score
foldBackScore = function(chms,start,end,chromIndex,invIndex,windowGR,runChrom="7",minRearr=5)
	{
	print("foldback")
	# empty count objects
	counts = matrix(0,nrow=2,ncol=length(windowGR))
	# windowed proportion  of ++ or --
	# counts inverted rearrangements in bins
	index = which(invIndex&chromIndex)
	if(length(index)>0)
		{
		invGR = as(paste0(chms[index],":",start[index],"-",end[index]),"GRanges")
		overlapsInv = findOverlaps(invGR,windowGR)
		info = table(overlapsInv@to)
		counts[1,as.numeric(names(info))] = info
		}
	# counts non-inverted rearrangements in bins
	index = which(!invIndex&chromIndex)
	if(length(index)>0)
		{
		noninvGR = as(paste0(chms[index],":",start[index],"-",end[index]),"GRanges")
		overlapsNoninv = findOverlaps(noninvGR,windowGR)
		info = table(overlapsNoninv@to)
		counts[2,as.numeric(names(info))] = info
		}
	# rearrangement fold-back inversion score
	f = counts[1,]/colSums(counts)
	f[which(colSums(counts)<minRearr)] = NA
	return(f)
	}

# get CN score for a single window
singleCNscore = function(CN,minLength=8,minMatch=8,maxError=0,softDir="/home/chris/software/BFB-Zakov/mckinsel-bfb-7127b9339b2b/bin",tmpFile="~/tmp.txt")
	{
	print("singleCN")
	if(length(CN)<minLength) return(list(d=NA,matchString=NA))
	# string for Zakov
	rearrString = paste(c("[",paste(CN,collapse=","),"]"),collapse="")
	print(rearrString)
	# run Zakov
	#www.bitbucket.org/mckinsel/bfb
	system(paste0("java -cp ",softDir," bfb.BFB_Algo counts:",rearrString," mode:substring model:Poisson maxError:",maxError," > ",tmpFile)) 
	# read Zakov results 
	bfbOut = read.table("~/tmp.txt",sep="\t",head=FALSE,as.is=TRUE)
	#system("rm ~/tmp.txt")
	# data munge
	info = strsplit(bfbOut[2,],split="[:]")[[1]]
	matchString = strsplit(info[1],split="[(]")[[1]][1]
	n = strsplit(info[2],split="/")[[1]][1]
	n = gsub(" ","",n)
	if(n<minMatch) return(list(d=NA,matchString=matchString))
	# distance
	d = rev(info)[1]
	d = gsub(" ","",d)
	d = gsub(")","",d)
	d = as.numeric(d)
	return(list(d=d,matchString=matchString))
	}

# get CN score for all windows
CNscore = function(chrom,start,end,CN,windowGR,minLength=8,minMatch=8,maxError=0,softDir="/home/chris/software/BFB-Zakov/mckinsel-bfb-7127b9339b2b/bin",tmpFile="~/tmp.txt")
	{
	print("CN")
	if(length(chrom)<1) return(NA)
	# genomic ranges for CN
	CNGR = as(paste0(chrom,":",start,"-",end),"GRanges")
	# overlaps
	overlaps = findOverlaps(CNGR,windowGR)
	indices = sapply(1:length(windowGR),
		FUN=function(x) overlaps@from[which(overlaps@to==x)],
		simplify=FALSE)
	# get score
	print("d")
	d = sapply(indices,FUN=function(x) singleCNscore(CN[x],minLength=minLength,minMatch=minMatch,maxError=maxError,softDir=softDir,tmpFile=tmpFile))
	return(d["d",])
	}

# combine CN score and foldback score
combineFun = function(d,f,w=0.5)
	{
	(w*d)+((1-w)*(1-f))
	}




# function to identify BFB on a single chromosome
identifyBFBwindow = function(chrom1bedpe,chrom2bedpe,
			strand1bedpe,strand2bedpe,
			start1bedpe,start2bedpe,
			end1bedpe,end2bedpe,
			chromSeg,startSeg,
			endSeg,CNseg,runChrom="7",
			windowSize=50000000,windowGap=5000000,
			threshold=0.18,
			minLength=12,minMatch=12,maxError=0,
			minRearr=5,
			softDir="/home/chris/software/BFB-Zakov/mckinsel-bfb-7127b9339b2b/bin")
	{
	print("identify BFB")
	# get data into right format 
	chms = c(chrom1bedpe,chrom2bedpe)
	orientation = rep(paste0(strand1bedpe,strand2bedpe),2)
	start = c(start1bedpe,start2bedpe)
	end = c(end1bedpe,end2bedpe)
	invIndex = orientation%in%c("++","--")
	chromIndex = chms==runChrom
	# get windows
	windowGR = setupBFB(windowSize=windowSize,
			windowGap=windowGap,
			runChrom=runChrom)
	print(paste0("Windows: ",length(windowGR)))
	# foldback inversion proportion
	f=foldBackScore(chms=chms,
			start=start,
			end=end,
			chromIndex=chromIndex,
			invIndex=invIndex,
			windowGR=windowGR,
			runChrom=runChrom,
			minRearr=minRearr)
	# copy number
	d = CNscore(chrom=chromSeg,start=startSeg,
		end=endSeg,CN=CNseg,
		windowGR=windowGR,
		minLength=minLength,
		minMatch=minMatch,
		maxError=maxError,
		softDir=softDir)
	# combine score
	combined = combineFun(d,f)
	# return
	return(list(decision=any(combined<threshold,na.rm=TRUE),scores=combined,threshold=threshold,d=d,f=f))
	}



# foldback score not counting rearrs twice
foldBackScoreBP = function(chrom1,chrom2,start1,start2,end1,end2,strand1,strand2,windowGR,runChrom,minRearr=5)
	{
	chromIndex = chrom1==runChrom|chrom2==runChrom
	if(sum(chromIndex)==0) return(NA)
	invIndex = (chrom1==chrom2)&paste0(strand1,strand2)%in%c("++","--")
	overlaps1 = findOverlaps(as(paste0(chrom1,":",start1,"-",end1),"GRanges"),windowGR)
	overlaps2 = findOverlaps(as(paste0(chrom2,":",start2,"-",end2),"GRanges"),windowGR)
	overlaps = unique(c(overlaps1@from,overlaps2@from))
	if(length(overlaps)<minRearr) return(NA)
	sum(invIndex[overlaps])/length(invIndex[overlaps])
	}

# function to identify BFB on a single chromosome - using set CN seg lengths
identifyBFBsetlength = function(chrom1bedpe,chrom2bedpe,
			strand1bedpe,strand2bedpe,
			start1bedpe,start2bedpe,
			end1bedpe,end2bedpe,
			segChrom,segStart,segEnd,segCN,
			runChrom="7",
			maxError=0,
			Length=8,
			minRearr=5,
			threshold=0.18,
			softDir="/home/chris/software/BFB-Zakov/mckinsel-bfb-7127b9339b2b/bin",
			rearrMethod="sep",
			tmpFile="~/tmp.txt")
	{
	if(length(segChrom)<Length) return(list(decision=NA,
						score="Too few segments: length(segChrom)<Length",
						d=NA,
						f=NA,
						threshold=threshold))
	# get data into right format 
	chms = c(chrom1bedpe,chrom2bedpe)
	orientation = rep(paste0(strand1bedpe,strand2bedpe),2)
	start = c(start1bedpe,start2bedpe)
	end = c(end1bedpe,end2bedpe)
	invIndex = orientation%in%c("++","--")
	chromIndex = chms==runChrom
	# set window starts
	windowStarts = 1:(length(segChrom)-(Length-1))
	scores = sapply(windowStarts,FUN=function(x) 
		{
		# which CN segs in this window
		index = x:(x+(Length-1))
		# distances
		d = singleCNscore(CN=segCN[index],
				minLength=Length,
				minMatch=Length,
				maxError=maxError,
				softDir=softDir,
				tmpFile=tmpFile)
		# genomic range of window
		gr = as(paste0(unique(segChrom[index]),":",min(segStart[index]),"-",max(segEnd[index])),"GRanges")
		# fold back score
		if(rearrMethod=="sep")
			{
			f = foldBackScore(chms=chms,
				start=start,
				end=end,
				chromIndex=chromIndex,
				invIndex=invIndex,
				windowGR=gr,
				runChrom=runChrom,
				minRearr=minRearr)
			} else {
			f = foldBackScoreBP(chrom1=chrom1bedpe,chrom2=chrom2bedpe,
				start1=start1bedpe,start2=start2bedpe,
				end1=end1bedpe,end2=end2bedpe,
				strand1=strand1bedpe,strand2=strand2bedpe,
				windowGR=gr,runChrom=runChrom,minRearr=minRearr)
			}
		score = combineFun(d$d,f)
		window=paste0(unique(segChrom[index]),":",min(segStart[index]),"-",max(segEnd[index]))
		segString=paste(segCN[index],collapse=",")
		return(list(score=score,d=d$d,f=f,window=window,segString=segString,matchString=d$matchString))
		})
	return(list(decision=any(unlist(scores["score",])<threshold,na.rm=TRUE),
			score=unlist(scores["score",]),
			window=unlist(scores["window",]),
			d=unlist(scores["d",]),
			f=unlist(scores["f",]),
			segString=unlist(scores["segString",]),
			matchString=unlist(scores["matchString",]),
			threshold=threshold))
	}
