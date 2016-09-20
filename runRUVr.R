

##  RUVr function; runRUVr(set, res,k1, k2)
#set = input SeqExpression set for normalization, same as used in edgeR glmFit to obtain residuals
#res = residual output from edgeR glmFit
#k1 = to loop through increasing vlaues of k, this is starting point
#k2 = highest k to test
#n.s = number of samples
###################################
## - testing runRUVr(set.uq, res.UQ, k1,k2)
#         set <- set.raw
#         res <- res.uq
#         i <-1
#         k2 =2


require(sva)
require(Biobase)    
require(RUVSeq)
require(EDASeq)
require(ShortRead)
##source(paste(.libPaths()[1], "mzw","runEdgeR_1f.R", sep = "/"))


runRUVr<- function(set,res,k1,k2,n.s) {
  rr <- k2+1; pv <- matrix(dimnames = list(c(k1:rr),c("K","PC1","PC2", "PC3", "PC4", "Cum1.3")), nrow=rr, ncol = 6)
  set <- betweenLaneNormalization(set, which="upper", offset = FALSE, round = TRUE)
  nc <- normCounts(set); data.prc <- prcomp(nc, scale = TRUE); p <- (data.prc$sdev ^2)/sum(data.prc$sdev ^2); pv[rr,] <- c(0,p[1:4], sum(p[1:3]));
  controls <- rownames(set)
  
  for(i in k1:k2) {
    set.r <- RUVr(set,controls, k=i, res)
    assign(paste("set_r", i, sep="."), set.r, envir = .GlobalEnv)
    
    #--run edgeR for norm cpm
    a <- DGEList(counts=normCounts(set.r))
    cpm <- cpm(a, normalized.lib.sizes=FALSE, log=FALSE, prior.count=1)
    assign(paste("dge",i,sep="."), cpm, envir = .GlobalEnv)
    
    #--diagnostics
    nc <- normCounts(set.r); data.prc <- prcomp(nc, scale = TRUE); p <- (data.prc$sdev ^2)/sum(data.prc$sdev ^2); pv[i,] <- c(i,p[1:4], sum(p[1:3]));
    colnames(set.r) <- paste(d.factor, seq(1:n.s), sep="")
    fn <- paste("rle_",i,"_",n.s,".png",sep=""); png(fn); plotRLE(set.r,outline=FALSE,ylim=c(-2, 2),col=cc); dev.off();
    fn <- paste("pca_",i,"_",n.s,".png",sep=""); png(fn); plotPCA(log2(normCounts(set.r)+1),col=cc,cex=1.2); dev.off();
    fn <- paste("pca_Log2CPM_",i,".png",sep=""); png(fn); plotPCA(log2(cpm+1), k=2, col=cc, main = paste("RUVr Log2CPM k", i, sep = "")) ; dev.off();

  }
  
  assign("pca.e", pv, envir = .GlobalEnv)
}

