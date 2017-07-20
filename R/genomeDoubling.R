# ======================================================================================
#				SIMULATION FUNCTIONS
# ======================================================================================

# simulate a single sample
simSample = function(N, # number of copy number alterations
                     Ps, # probabiltiies for each chromosome arm to be altered
	                   sepP=FALSE # whether Ps are provided per arm
                     )
  {
	if(!sepP)
	  {
	  # generate neutral allele counts
	  allelecounts = matrix(1,ncol=length(Ps),nrow=2)
	  # loop over alterations
	  for(i in 1:N)
  		{
		  # which chromosome to alter (ignoring entirely deleted chromosomes)
		  zeroIndex = which(!colSums(allelecounts)==0)
		  chrom = sample((1:length(Ps))[zeroIndex],size=1,prob=Ps[zeroIndex])
		  print(chrom)
		  # which allele to alter
		  allele = rbinom(1,1,allelecounts[2,chrom]/sum(allelecounts[,chrom]))
		  print(allele)
		  # whether to delete or amplify
		  change = rbinom(1,1,0.5)
		  if(change==0) change=c(-1)
		  print(change)
		  # make change
		  allelecounts[allele+1,chrom] = allelecounts[allele+1,chrom]+change
		  print(allelecounts)
	    }
	 } else {
	   Ps = t(Ps)
	   # generate neutral allele counts
	   allelecounts = matrix(1,ncol=ncol(Ps),nrow=2)
	   # loop over alterations
	   for(i in 1:N)
	   {
	     # which chromosome & allele to alter (ignoring entirely deleted chromosomes)
	     zeroIndex = which(!allelecounts==0)
	     chromAllele = sample((1:length(Ps))[zeroIndex],size=1,prob=Ps[zeroIndex])
	     print(chromAllele)
	     # whether to delete or amplify
	     change = rbinom(1,1,0.5)
	     if(change==0) change=c(-1)
	     print(change)
	     # make change
	     allelecounts[chromAllele] = allelecounts[chromAllele]+change
	     print(allelecounts)
	   } 
	 }
	return(allelecounts)
	}


# get score
sampleSimScore = function(allelecounts,  # allelecounts from simSample()
                          testVal)	# prop chrom arms with major allele >2 CN in real sample 
	{
	score = sum(pmax(allelecounts[1,],allelecounts[2,])>2)>testVal
	return(score)
	}

# simulate single replicate
getSingleScore = function(N,Ps,testVal,testFUN=sampleSimScore,diag=FALSE,sepP=FALSE)
	{
	allelecounts = simSample(N,Ps,sepP=sepP)
	score = testFUN(allelecounts,testVal)
	if(diag) return(list(score=score,allelelcounts=allelecounts))
	return(mean(score))
	}

getChromArmAbberations = function(seg,armLims)
	{
	abberations = abs(1-seg[,c("tumCNa","tumCNb")])
	apply(armLims,MARGIN=1,FUN=function(x)
		{
		index = which(grepl(abberations))
		})
	}

# single sample simulation
singleSamp = function(N,Ps,testProps,nReps=10000,testFUN=sampleSimScore,diag=FALSE,sepP=FALSE,doParallel=FALSE,nCores=NULL)
	{
  if(doParallel)
    {
	  res = unlist(mclapply(1:nReps,FUN=function(x) getSingleScore(N,Ps,testProps,testFUN=testFUN,diag=diag,sepP=sepP),
	                        mc.cores=nCores))
    } else {
    res = replicate(nReps,getSingleScore(N,Ps,testProps,testFUN=testFUN,diag=diag,sepP=sepP))
    }
  res
  }

# ======================================================================================
#				DATA MUNGING FUNCTIONS
# ======================================================================================


# split segments that cross arm boundaries
splitSeg = function(x,split,startCol=4,endCol=5)
	{
	splitseg =  rbind(x,x)
	splitseg[1,endCol] = split
	splitseg[2,startCol] = split+1
	return(splitseg)
	}

# function to get CN of chromsome arms from seg file
getArmCN = function(seg,armLims,
                    aCol = 9,
                    bCol = 10,
                    startCol = 4,
                    endCol = 5,
                    chromCol = 3)
	{
	# get which segments cross arm boundaries
	index = ((1:(nrow(armLims)/2))*2)-1
	armMiddles = armLims[index,2]
	names(armMiddles) = gsub("q","",gsub("p","",gsub("chr","",names(armMiddles))))
	crossMiddle = apply(seg,MARGIN=1,FUN=function(x) 
		{
		index = which(names(armMiddles)==paste0(x[chromCol]))
		as.numeric(x[startCol])<armMiddles[index]&
		as.numeric(x[endCol])>armMiddles[index]
		})
	tmp = seg[which(crossMiddle),]
	if(nrow(tmp)>0)
		{
		# split segments that cross middle of chromosome 
		segSub = seg[-which(crossMiddle),]
		tmp = sapply(1:nrow(tmp),
			FUN=function(x) splitSeg(tmp[x,,drop=FALSE],
				armMiddles[tmp[x,chromCol]],
				startCol=startCol,
				endCol=endCol),simplify=FALSE)
		tmp = do.call(rbind,tmp)
		segNew = rbind(segSub,tmp)
		crossMiddleCheck = apply(segNew,MARGIN=1,FUN=function(x) 
			{
			index = which(names(armMiddles)==paste0(x[chromCol]))
			as.numeric(x[startCol])<armMiddles[index]&
			as.numeric(x[endCol])>armMiddles[index]
			})
		print(sum(crossMiddleCheck))
		} else {
		segNew = seg
		}
	# add chromosome arm to seg files
	arm = apply(segNew,MARGIN=1,FUN=function(x) as.numeric(x[endCol])<=
		as.numeric(armMiddles[paste0(unlist(x[chromCol]))]))
	segNew$arm = c("q","p")[arm+1]
	segNew$chromArm = paste0(segNew[,chromCol],c("q","p")[arm+1])
	armCol = ncol(segNew)-1
	# relative lengths of each segment
	relLength = apply(segNew,MARGIN=1,FUN=function(x)
		{
		index = which(rownames(armLims)==
			paste0("chr",paste0(x[chromCol]),x[armCol]))
		segLength = as.numeric(x[endCol])-as.numeric(x[startCol])
		armLength = armLims[index,2]-armLims[index,1]
		segLength/armLength
		})
  
	index = which(relLength>1)
	print(length(index>0))
 	# arm-level copy number
	CN = sapply(unique(segNew$chromArm),FUN=function(x) 
		{
		index = which(segNew$chromArm==x)
		c(sum(relLength[index]*segNew[index,aCol]),
		sum(relLength[index]*segNew[index,bCol]))
		})
	return(CN)
	}

# function to get arm limits in right format
getArmLims = function(file)
	{
	armLims = read.table(file,sep="\t")
	chrom = unique(armLims[,1])
	armLims = sapply(chrom,FUN=function(x) 
		sapply(c("p","q"),FUN=function(y)
			range(armLims[which(armLims[,1]==x&
				grepl(y,armLims[,4])),2:3]),
			simplify=FALSE),
		simplify=FALSE)
	chrom = rep(chrom,each=2)
	armLims = do.call(rbind,sapply(armLims,
		FUN=function(x) do.call(rbind,x),simplify=FALSE))
	rownames(armLims) = paste0(chrom,rownames(armLims))
	return(armLims)
	}

# ======================================================================================
#			SIMULATION VARIABLE FUNCTIONS
# ======================================================================================



# function to get variables for simulations
getVars = function(CN, # CN from getArmCN()
                   armLims, # limits of chromosome arms
                   testFUN=ABSOLUTEscore, # score function
                   samples,       # samples
                   sepP=FALSE      # sep P and N for each copy? 
                  )
	{
	rownames(armLims) = gsub("chr","",rownames(armLims))
	armLengths = armLims[,2]-armLims[,1]
	# CN changes
	CNchanges = round(abs(1-CN)) # how different to 1 is CN
	# total events
	if(!sepP)
	  {
	  Ns = sapply(samples,FUN=function(x) 
		  {
		  index = grep(x,rownames(CNchanges))
		  sum(CNchanges[index,],na.rm=TRUE)
		  })
	  names(Ns) = samples
	  } else {
	  Ns = apply(CNchanges,MARGIN=1,FUN=function(x) sum(x))
	  names(Ns) = rownames(CN)
	  }
	# rates for each arm (per base)
	if(!sepP)
	  {
	  #rates = sapply(samples,FUN=function(x) 
  	#	{
	#	  index = grep(x,rownames(CN))
	#	  sapply(colnames(CN),FUN=function(y)
  	#		{
	#		  sum(CNchanges[index,y],na.rm=TRUE)/armLengths[y]
	#		  })
	#	  })
	  rates = sapply(samples,FUN=function(x) 
  		{
		  index = grep(x,rownames(CN))
		  sapply(colnames(CN),FUN=function(y)
  			{
			  sum(CNchanges[index,],na.rm=TRUE)/armLengths[y]
			  })
		  })
	  rates[which(is.na(rates))] = 0
	  rates = t(rates)
	  } else {
	 #rates = sapply(rownames(CN),FUN=function(x) 
	 #   {
	 #     sapply(colnames(CN),FUN=function(y)
	 #     {
	 #       sum(CNchanges[x,y],na.rm=TRUE)/armLengths[y]
	 #     })
	 #   })
	 rates = sapply(rownames(CN),FUN=function(x) 
	    {
	      sapply(colnames(CN),FUN=function(y)
	      {
	        sum(CNchanges[x,],na.rm=TRUE)/armLengths[y]
	      })
	    })  
	    rates[which(is.na(rates))] = 0
	    rates = t(rates)
	  }
	# normalise to probabilities
	probs = apply(rates,MARGIN=1,FUN=function(x) x/sum(x))
	if(!sepP)
	  {
	  colnames(probs) = samples
	  } else {
	  colnames(probs) = rownames(CN)
	  }
	rownames(probs) = colnames(CN)
	# testVals
	CN = round(CN)
	testVals = sapply(samples,FUN=function(x)
		{
		index = grep(x,rownames(CN))
		testFUN(CN[index,])
		})
	names(testVals) = samples
	# return for simulation
	return(list(Ps=probs,Ns=Ns,testVals=testVals))
	}



# ======================================================================================
#				SCORE FUNCTIONS
# ======================================================================================

# function for ABSOLUTE score (even high)
ABSOLUTEscore = function(x)
	{
	sampMax = pmax(x[1,],x[2,])
	mean(sampMax%%2==0,na.rm=TRUE)
	}

# function to test ABSOLUTE score
ABSOLUTEtest = function(allelecounts,testVal)
	{
	ABSOLUTEscore(allelecounts)>=testVal
	}

# Swanton method score
swantonScore = function(x)
	{
	sampMax = pmax(x[1,],x[2,])
	mean(sampMax>=2,na.rm=TRUE)
	}

# Swanton method test
swantonTest = function(allelecounts,testVal)
	{
	swantonScore(allelecounts)>=testVal
	}

# ======================================================================================
#				WRAPPER FUNCTION
# ======================================================================================

# function to run full genome doubling analysis
genomeDoubling = function(segFile,	# segment file
                       armFile,	# cytoband file
                       nReps=1000,	# number of simulations
                       scoreFUN=ABSOLUTEscore,	# score to get p value
                       doSim=TRUE,	# whether to simulate
                       withCNb=TRUE,	# does seg file have b CN?
                       sampleCol = 1,	# seg sample column
                       aCol = 9,	# seg a CN column
                       bCol = 10,	# seg b CN column
                       startCol = 4,	# seg start column
                       endCol = 5,	# seg end column
                       chromCol = 3,	# seg chrom column
                       totCol = 11,	# seg total CN column
                       head=TRUE,	# does seg file have headers?
                       diag=FALSE,
                       sepP=FALSE, # vars separate for each allele?
                       doParallel=FALSE, # run in parallel 
                       nCores = NULL # number of cores
                       ) 
	{
  if(doParallel&is.null(nCores)) nCores=detectCores()
	# get arm limits
	print("load arm lims")
	armLims = getArmLims(armFile)
	# get segs
	if(is.character(segFile))
		{
		print("load seg")
		if(grepl("[.]csv",segFile))
			{
			seg = read.csv(segFile,head=head,as.is=TRUE)
			} else {
			seg = read.table(segFile,head=head,sep="\t",as.is=TRUE)  
			}
		} else {
		seg = segFile
		}
	# reformat sex chromosome names
	print("formatting")
	if(!any(grepl("X",seg[,chromCol])))
		{
		seg[grep("23",seg[,chromCol]),chromCol] = 
			gsub("23","X",seg[grep("23",seg[,chromCol]),chromCol])
		seg[grep("24",seg[,chromCol]),chromCol] = 
			gsub("24","Y",seg[grep("24",seg[,chromCol]),chromCol])
		}
	# reformat hg/grch chromsome names
	if(any(grepl("chr",seg[,chromCol])))
		{
		seg[,chromCol] = gsub("chr","",seg[,chromCol])
		}
	# if no CN b
	if(!withCNb)
		{
		seg = cbind(seg,seg[,totCol]-seg[,aCol])
		bCol = ncol(seg)
		}
	# if no sample column
	if(is.na(sampleCol))
		{
		seg = cbind(seg,rep("sample",nrow(seg)))
		sampleCol = ncol(seg)
		}
	# get arm CNs
	print("sample CN")
	samples = unique(seg[,sampleCol])
	CNs = sapply(samples,FUN=function(x) getArmCN(seg[which(seg[,sampleCol]==x),],
		armLims,
		aCol = aCol,
		bCol = bCol,
		startCol = startCol,
		endCol = endCol,
		chromCol = chromCol),simplify=FALSE)
	names(CNs) = samples
	CNcomb = do.call(smartbind,CNs)
	# remove sex chromsome copy numbers
	xIndex = grep("X",colnames(CNcomb))
	if(length(xIndex)>0) CNcomb = CNcomb[,-xIndex]
	yIndex = grep("Y",colnames(CNcomb))
	if(length(yIndex)>0) CNcomb = CNcomb[,-yIndex]
	# get simulation variables
	print("get simulation variables")
	simVars = getVars(CNcomb,armLims,testFUN=scoreFUN,samples=samples,sepP=sepP)
	# make test function
	funtest = function(allelecounts,testVal)
		{
		scoreFUN(allelecounts)>=testVal
		}
	if(doSim)
		{
		print("run simulation")
		# run simulations
	  if(!sepP)
	    {
	    # single P per arm
		  res = sapply(names(simVars$Ns),
  			FUN=function(x) singleSamp(N=simVars$Ns[x],
			    testProps=simVars$testVals[x],
			    Ps=simVars$Ps[,x],
			    nReps=nReps,
			    testFUN=funtest,
			    diag=diag,sepP=sepP,
			    doParallel=doParallel,
			    nCores=nCores))
	    } else {
	   # separate Ps per allele
	    res = sapply(samples,
	      FUN=function(x)
	        {
	        index = grep(x,names(simVars$Ns))
	        singleSamp(N=sum(simVars$Ns[index]),
	          testProps=simVars$testVals[x],
	          Ps=simVars$Ps[,index],
	          nReps=nReps,
	          testFUN=funtest,
	          diag=diag,sepP=sepP,
	          doParallel=doParallel,
	          nCores=nCores)})
	    }
		return(res)
		} else {
		return(list(CN=CNcomb,vars=simVars))
		}
	}
