# ======================================================================================
#			UTILITY FUNCTIONS
# ======================================================================================

# read in a file
readFile = function(file,head)
	{
	ending = rev(strsplit(file,split="[.]")[[1]])[1]
	if(ending=="csv") return(read.csv(file,head=head))
	if(ending%in%c("txt","tsv")) return(read.table(file,sep="\t",head=head))
	return(read.table(file,sep="\t",head=head))	
	}


# load cytoBand file
loadCytoBand = function(file=NULL)
	{
	if(is.null(file))
		{
		tmpEnv = new.env()
		data(list="cytoBand", package='CNomplexity',envir=tmpEnv)
		return(tmpEnv[["cytoBand"]])
		} else {
		return(readFile(file))
		}
	}
