
get_nullprop <-function(input_pvalues)
{
  input_pvalues <- as.matrix(input_pvalues)
  nullprop <-null_estimation(input_pvalues)
  return(nullprop)
}

Corrected_qqplot <- function(input_pvalues,nullprop,exact=1)
{
  input_pvalues <- as.matrix(input_pvalues)
  
  pnull <- adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                                     nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=exact)
  pmax <- apply(input_pvalues,1,max)
  
  nmed <- length(pmax)
  pindex=1:nmed
 
  xmax <- max(c(-log(pnull[order(pnull,decreasing=T)],base=10)[pindex],-log10(1/nmed)))
  ymax <- max(c(-log(pmax[order(pmax,decreasing=T)],base=10)[pindex],-log(pmax[order(pmax,decreasing=T)],base=10)[pindex]))
  if (xmax>0.8*ymax & xmax<1.25*ymax)
  {
    xmax <- max(xmax,ymax)
    ymax <- max(xmax,ymax)
  }
  
  plot((-log(pnull[order(pnull,decreasing=T)],base=10))[pindex],(-log(pmax[order(pmax,decreasing=T)],base=10))[pindex],xlab="log base 10 (expected null p-values)", ylab="log base 10 (observed p-values)",col=2,xlim=c(0,xmax),ylim=c(0,ymax))
  #points((-log((nmed:1)/nmed,base=10))[pindex],(-log(pmax[order(pmax,decreasing=T)],base=10))[pindex],pch=2,col=2)
  #legend(0.1,max(-log(pmax,base=10)),c("Uniform null","Mixture null"),pch=2:1,col=2:3,bty="n")
  abline(0,1,lty=2)
  title("q-q plot for testing mediation")
}

get_fdr_est <- function(input_pvalues,nullprop,exact=1,cutoff=0.05)
{
  input_pvalues <- as.matrix(input_pvalues)
  #nullprop <-null_estimation(input_pvalues)
  fdr <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                 nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=exact)
  nsig <- sum(fdr<=cutoff)
  if (sum(fdr<=cutoff)>0)
  {
    #idx <- order(fdr)
    #idx1 <- idx[1:min(c(maxn,sum(fdr<cutoff)))]
    if(is.null(rownames(input_pvalues)))
    {
      rownames(input_pvalues)=1:nrow(input_pvalues)
    }
    #print(rownames(input_pvalues)[idx1])
    output_p = as.data.frame(cbind(input_pvalues,fdr))
    colnames(output_p) <- c("p1","p2","FDR")
    output_p = output_p[fdr<=cutoff,]
    output_p = output_p[order(output_p[,3]),]
    return(list(nsig=nsig,output_p=output_p))
  } else return(list(nsig=nsig,output_p=NULL))
}

get_fwer_sig <- function(input_pvalues,nullprop,exact=1,cutoff=0.05)
{
  input_pvalues <- as.matrix(input_pvalues)
  #nullprop <-null_estimation(input_pvalues)
  fwercut <- fwer_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                 nullprop$alpha1,nullprop$alpha2,input_pvalues,alpha=cutoff,exact=exact)
  pmax <- apply(input_pvalues,1,max)
  nsig <- sum(pmax<=fwercut)
  if (nsig>0){
    output_p = as.data.frame(cbind(input_pvalues[pmax<=fwercut,]))
    colnames(output_p) <- c("p1","p2")
    return(list(nsig=nsig,output_p=output_p))
  } else return(list(nsig=nsig,output_p=NULL))
    
}
