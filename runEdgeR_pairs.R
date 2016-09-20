

require(sva)
require(Biobase)    
require(RUVSeq)
require(EDASeq)
require(ShortRead)
require(edgeR)

##edgeRUsersGuide(view=TRUE)


runEdgeR.p <-
function(x1,x2,y,z,m){

   
    
    ##  EdgeR function; runEdgeR(x,y,z,m)
    #y = data set
    #x = factors for the model, design is pairs so x is c("donor", "treatment")
    #z = name to append
    #m = method for calculating norm factors , "TMM" or "none"
    # check the glmLRT for coeficcient or contrasts before running
    
    
          set <- y
          
          design <-  model.matrix(~x1+x2, data = pData(set))
        if(m == "normCounts") {a <- DGEList(counts=normCounts(set));print("normCounts")}
        if(m != "normCounts") {a <- DGEList(counts=counts(set)); a <- calcNormFactors(a, method=m);print("m")}
        
          a <- estimateGLMCommonDisp(a, design)
          a <- estimateGLMTagwiseDisp(a, design)
          cpm <- cpm(a, normalized.lib.sizes=TRUE, log=FALSE, prior.count=1)
          
          # fit the models
          fit <- glmFit(a, design)
          lrt <- glmLRT(fit)
          print(topTags(lrt, n=10, adjust.method="BH", sort.by="PValue"))
          
          # calculate residuals here to have available for RUVr
          res <- residuals(fit, type="deviance")
          
          assign(paste("cpm",z,sep="."), cpm, envir = .GlobalEnv)
          assign(paste("y",z,sep="."), a, envir = .GlobalEnv)
          assign(paste("fit",z,sep="."), fit, envir = .GlobalEnv)
          assign(paste("lrt",z,sep="."), lrt, envir = .GlobalEnv)
          assign(paste("res",z,sep="."), res, envir = .GlobalEnv)
          
        }
