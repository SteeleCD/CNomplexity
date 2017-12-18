# ============================================================================
#		REARRANGEMENT SIGNATURES SETUP
# ============================================================================
# function to get whether rearrangments are clustered or not
getClustered = function(chrom1,pos1,chrom2,pos2,threshold=NULL)
	{
	pos1 = as.numeric(pos1)
	pos2 = as.numeric(pos2)
	# get segments per chromosome
	dists = sapply(unique(c(chrom1,chrom2)),FUN=function(x)
		{
		bp = matrix(NA,ncol=2,nrow=0)
		# positions separate for each member of breakpoint
		pos=c()
		index1 = which(chrom1==x)
		if(length(index1)>0)
			{
			pos = c(pos,pos1[index1])
			bp = rbind(bp,cbind(index1,1))
			}
		index2 = which(chrom2==x)
		if(length(index2)>0)
			{
			pos = c(pos,pos2[index2])
			bp = rbind(bp,cbind(index2,2))
			} 
		# distances between breakpoint
		dists = diff(sort(pos))
		# segment
		forCN = data.frame(chrom=x,pos=sort(pos),dist=dists[order(pos)])
		return(list(info=forCN,seg=pcf(forCN,gamma=25,kmin=10)))
		},simplify=FALSE)
	# threshold below which to call clustered
	if(is.null(threshold)) threshold = 0.1*mean(unlist(sapply(dists,FUN=function(x) x$info[,3])),na.rm=TRUE)
	# which segments are below threshold
	regions = do.call(rbind,sapply(dists,FUN=function(x)
		{
		x$seg[which(x$seg$mean<threshold),]
		},simplify=FALSE))
	if(nrow(regions)==0)
		{
		return(rep("unclustered",length(chrom1)))
		}
	# which rearrangements are in clustered regions
	print(paste0(regions$chrom,":",regions$start.pos,"-",regions$end.pos))
	regionGR = as(paste0(regions$chrom,":",regions$start.pos,"-",regions$end.pos),"GRanges")
	clustered = rep("unclustered",length=length(pos1))
	print(paste0(chrom1,":",pos1,"-",pos1))
	index1 = findOverlaps(as(paste0(chrom1,":",pos1,"-",pos1),"GRanges"),regionGR)@from
	if(length(index1)>0) clustered[index1] = "clustered"
	print(paste0(chrom2,":",pos2,"-",pos2))
	index2 = findOverlaps(as(paste0(chrom2,":",pos2,"-",pos2),"GRanges"),regionGR)@from
	if(length(index2)>0) clustered[index2] = "clustered" 
	return(clustered)
	}

# get rearrangement classes
getRearrClass = function(sample,
	chrom1,
	start1,
	end1,
	chrom2,
	start2,
	end2,
	type)
	{
	clustered = unlist(sapply(unique(sample),FUN=function(x)
		{
		index = which(sample==x)
		getClustered(chrom1=chrom1[index],
			pos1=start1[index],
			chrom2=chrom2[index],
			pos2=start2[index])
		}))
	size =  sapply(1:nrow(data), FUN=function(x) 
		{
		if(type[x]=="translocation") return(NA)
		length = end2[x]-start1[x]
		if(length<1000) return("<1kb")
		if(length>=1000&length<10000) return("1-10kb")
		if(length>=10000&length<100000) return("10-100kb")
		if(length>=100000&length<1000000) return("100kb-1Mb")
		if(length>=1000000&length<=10000000) return("1Mb-10Mb")
		if(length>10000000) return(">10Mb")
		})
	return(paste0(type,":",size,":",clustered))
	}

# ============================================================================
#		RUN NMF
# ============================================================================

# bootstrap with serenas method
bootNMF = function(data,rank=2:10)
	{
	require(NMF)
	d = apply(data,MARGIN=1,FUN=function(x)
		{
		sampled = sample(length(x),sum(x),replace=TRUE,prob=x/sum(x))
		sapply(1:length(x),FUN=function(x) sum(sampled==x))
		})
	colnames(d) = rownames(data)
	rownames(d) = colnames(data)
	nmfres = nmf(x=d,
		rank=rank,
		method="brunet",
		nrun=30)
	sil = sapply(nmfres$fit,FUN=function(x) silhouette(x)[,"sil_width"])
	return(list(nmf=nmfres,sil=sil))
	}


# ============================================================================
#		FIT METRICS TO DETERMINE N
# ============================================================================
# NMF fit metrics
nmfFitMetrics = function(N,res,data,plotsil=FALSE,maxiter=1000,nstarts=20,algo="Hartigan-Wong",doMetric=TRUE,clusterMeth="kmeans",hclustInput="sigs",kmeansDist="cosine")
	{
	# get matrices
	Ps = sapply(res,FUN=function(x) x$nmf$fit[[N]]@fit@W,simplify=FALSE)
	Es = sapply(res,FUN=function(x) x$nmf$fit[[N]]@fit@H,simplify=FALSE)
	# combine P signatures
	sigs = t(do.call(cbind,Ps))
	exps = do.call(rbind,Es)
	# get distances
	dists = matrix(NA,ncol=nrow(sigs),nrow=nrow(sigs))
	for(x in 1:nrow(sigs)) for(y in x:nrow(sigs))
		{
		dist = 1-cosine(sigs[x,],sigs[y,])
		dists[x,y] = dist
		dists[y,x] = dist
		}
	distsObj = as.dist(dists)
	# cluster signatures into right number
	if(clusterMeth=="kmeans")
		{
		if(kmeansDist=="cosine")
			{
			clusters = pam(distsObj,k=as.numeric(N))
			clusters = clusters$clustering
			} else {
			clusters = kmeans(sigs,
				centers=as.numeric(N),
				iter.max=maxiter,
				nstart=nstarts,
				algorithm=algo)
			clusters = clusters$cluster
			}
		} else {
		if(hclustInput=="sigs")
			{
			clusters = hclust(dist(sigs),method="ward.D2")
			clusters = cutree(clusters,as.numeric(N))
			} else {
			clusters = hclust(distsObj,method="complete")
			clusters = cutree(clusters,as.numeric(N))
			}
		}	
	# pca
	pca = prcomp(sigs)
	plot(pca$x,col=clusters)
	# get group centroids
	centroids = sapply(unique(clusters),FUN=function(x) {
		index = which(clusters==x)
		colSums(sigs[index,])/length(index)
		})
	# get group exposures
	exposures = t(sapply(unique(clusters),FUN=function(x) {
		index = which(clusters==x)
		colSums(exps[index,])/length(index)
		}))
	if(doMetric)
		{
		# get frobenius distance
		frobenius = sqrt(sum((data-(centroids%*%exposures))^2))
		# get silhouette metric
		groups = unique(clusters)
		sil = sapply(1:nrow(sigs),FUN=function(x)
			{
			index = which(clusters==clusters[x])
			index = index[-which(index==x)]
			a = mean(dists[x,index])
			b = min(sapply(groups[-which(groups==clusters[x])],FUN=function(y)
				{
				mean(dists[x,which(groups==y)])
				}))
			(b-a)/max(c(a,b))
			})
		if(plotsil) barplot(sil[order(clusters)],col=sort(clusters))
		avsil = sapply(groups,FUN=function(x) mean(sil[which(clusters==x)]))
		return(list(sil=mean(avsil),frob=frobenius,P=centroids,E=exposures))
		} else {
		return(list(P=centroids,E=exposures))
		}
	}
