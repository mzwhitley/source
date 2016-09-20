require(Biobase)    
require(RUVSeq)
require(EDASeq)
require(ShortRead)
require(edgeR)


runEdgeR.1f.norm <-
  function(x,y,z){
    
    
    
    ##  EdgeR function; runEdgeR(x,y,z)
    #y = data set
    #x = factors for the model, design is simple one factor, ~x
    #z = name to append
 
    # check the glmLRT for coeficcient or contrasts before running
    
    set <- y
    
    w1 <- pData(set)
    design <-  model.matrix(~x, data = pData(set))
    a <- DGEList(counts=normCounts(set), group = x)
    ##a <- calcNormFactors(a, method=m)
    a <- estimateGLMCommonDisp(a, design)
    a <- estimateGLMTagwiseDisp(a, design)
    cpm <- cpm(a, normalized.lib.sizes=TRUE, log=FALSE, prior.count=1)
    
    # fit the models
    fit <- glmFit(a, design)
    lrt <- glmLRT(fit, coef=2:length(levels(x)))
    
    # calculate residuals here to have available for RUVr
    res <- residuals(fit, type="deviance")
    
    ## assign(paste("cpm",z,sep="."), cpm, envir = .GlobalEnv) # can derive from dge
    assign(paste("dge",z,sep="."), a, envir = .GlobalEnv)
    assign(paste("fit",z,sep="."), fit, envir = .GlobalEnv)
    assign(paste("lrt",z,sep="."), lrt, envir = .GlobalEnv)
    assign(paste("res",z,sep="."), res, envir = .GlobalEnv)
    assign(paste("set.norm",z,z,sep="."), set, envir = .GlobalEnv)
    
  }