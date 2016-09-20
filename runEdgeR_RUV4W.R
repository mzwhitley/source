require(sva)
require(Biobase)    
require(RUVSeq)
require(EDASeq)
require(ShortRead)
require(edgeR)


runEdgeR.ruv <-
  function(x,y,z,m,k){
    
    ##runEdgeR.ruv(d.factor, set_r.4, "4W_uq", "upperquartile", 4) ## with Ws 1-4, using origninal values and UQ
    
    
    ##  EdgeR function; runEdgeR(x,y,z,m)
#     y = set_r.4
#     x = d.factor
#     z = "4W_uq"
#     m = "upperquartile"
#     k = 4
    # check the glmLRT for coeficcient or contrasts before running
    
    ##  EdgeR function; runEdgeR(x,y,z,m)
    #y = data set
    #x = factors for the model, design is simple one factor, ~x
    #z = name to append
    #m = method for calculating norm factors , "TMM" or "none"
    #k = number of RUVr weights from pData to use
    # check the glmLRT for coeficcient or contrasts before running
    
    
    set <- y
    if (k==1) {design <-  model.matrix(~ x + W_1, data = pData(set))}
    if (k==2) {design <-  model.matrix(~ x + W_1 + W_2, data = pData(set))}
    if (k==3) {design <-  model.matrix(~ x + W_1 + W_2 + W_3, data = pData(set))}
    if (k==4) {design <-  model.matrix(~ x + W_1 + W_2 + W_3 + W_4, data = pData(set))}
    
  
    a <- DGEList(counts=counts(set), group = x)
    a <- calcNormFactors(a, method=m)
    a <- estimateGLMCommonDisp(a, design)
    a <- estimateGLMTagwiseDisp(a, design)
    cpm <- cpm(a, normalized.lib.sizes=TRUE, log=FALSE, prior.count=1)
    
    # fit the models
    fit <- glmFit(a, design)
    lrt <- glmLRT(fit, coef=2:length(levels(x)))
    
    # calculate residuals here to have available for RUVr
    res <- residuals(fit, type="deviance")
    
    assign(paste("cpm",z,sep="."), cpm, envir = .GlobalEnv)
    assign(paste("y",z,sep="."), a, envir = .GlobalEnv)
    assign(paste("fit",z,sep="."), fit, envir = .GlobalEnv)
    assign(paste("lrt",z,sep="."), lrt, envir = .GlobalEnv)
    assign(paste("res",z,sep="."), res, envir = .GlobalEnv)
    
  }