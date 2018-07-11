#### Functions

plot.2.file <- function(name, type="pdf", res.f=1){
  full.name <- paste(name,type,sep=".")
  if (type == "pdf"){
    pdf(full.name)
  }else if(type == "png"){
    png(full.name,width = 480*res.f, height = 480*res.f, pointsize = 12*res.f)
  }else if (type == "tiff"){
    tiff(full.name,width = 480*res.f, height = 480*res.f, pointsize = 12*res.f)
  }else if (type == "jpeg"){
    jpeg(full.name,width = 480*res.f, height = 480*res.f, pointsize = 12*res.f)
  }
}
chop_chr<-function  (name, splitcharacter = "__") {   strsplit(name, splitcharacter)[[1]][1]}

cellgenes.hist2<-function(data,name,col){
  Genes<-apply(data,2,function(x) sum(x>=1))
  # pdf(paste(getwd(),"/plots/","genes_cell",name,".pdf",sep=""))
  hist(Genes,breaks=50, col=col, main=paste("Genes per cell", paste("(",name,")",sep=""), sep=" "), sub=paste("median",median(Genes)))
  # dev.off()
}
name2id <- function(x,id) id[sub("\\_\\_chr\\w+","",id) %in% x]

SCseq <- setClass("SCseq", slots = c(expdata = "data.frame", ndata = "data.frame", fdata = "data.frame", distances = "matrix", tsne = "data.frame", cluster = "list", background = "list", out = "list", cpart = "vector", fcol = "vector", filterpar = "list", clusterpar = "list", outlierpar ="list" ))
setValidity("SCseq",function(object) {
              msg <- NULL
              if ( ! is.data.frame(object@expdata) ){
                msg <- c(msg, "input data must be data.frame")
              }else if ( nrow(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }else if (sum( apply( is.na(object@expdata),1,sum ) ) > 0 ){
                msg <- c(msg, "NAs are not allowed in input data")
              }else if (sum( apply( object@expdata,1,min ) ) < 0 ){
                msg <- c(msg, "negative values are not allowed in input data")
              }
              if (is.null(msg)) TRUE
              else msg
            })
setGeneric("filterdata", function(object, mintotal=3000, minexpr=5, minnumber=1, maxexpr=Inf, downsample=TRUE, dsn=1, rseed=17000) standardGeneric("filterdata"))
setMethod("filterdata",
          signature = "SCseq",
          definition = function(object,mintotal,minexpr,minnumber,maxexpr,downsample,dsn,rseed) {
            if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
            if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
            if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
            if ( ! ( is.numeric(downsample) | is.logical(downsample) ) ) stop( "downsample has to be logical (TRUE/FALSE)" )
            if ( ! is.numeric(dsn) ) stop( "dsn has to be a positive integer number" ) else if ( round(dsn) != dsn | dsn <= 0 ) stop( "dsn has to be a positive integer number" )
            object@filterpar <- list(mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, maxexpr=maxexpr, downsample=downsample, dsn=dsn)
            object@ndata <- object@expdata[,apply(object@expdata,2,sum,na.rm=TRUE) >= mintotal]
            if ( downsample ){
              set.seed(rseed)
              object@ndata <- downsample(object@expdata,n=mintotal,dsn=dsn)
            }else{
              x <- object@ndata
              object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE)) + .1 )
            }
            x <- object@ndata
            object@fdata <- x[apply(x>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
            x <- object@fdata
            object@fdata <- x[apply(x,1,max,na.rm=TRUE) < maxexpr,]
            return(object)
          })
setMethod("initialize",
          signature = "SCseq",
          definition = function(.Object, expdata ){
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@fdata <- expdata
            validObject(.Object)
            return(.Object)
          }
)

downsample <- function(x,n,dsn){
  x <- round( x[,apply(x,2,sum,na.rm=TRUE) >= n], 0)
  nn <- min( apply(x,2,sum) )
  for ( j in 1:dsn ){
    z  <- data.frame(GENEID=rownames(x))
    rownames(z) <- rownames(x)
    initv <- rep(0,nrow(z))
    for ( i in 1:dim(x)[2] ){
      y <- aggregate(rep(1,nn),list(sample(rep(rownames(x),x[,i]),nn)),sum)
      na <- names(x)[i]
      names(y) <- c("GENEID",na)
      rownames(y) <- y$GENEID
      z[,na] <- initv
      k <- intersect(rownames(z),y$GENEID)
      z[k,na] <- y[k,na]
      z[is.na(z[,na]),na] <- 0
    }
    rownames(z) <- as.vector(z$GENEID)
    ds <- if ( j == 1 ) z[,-1] else ds + z[,-1]
  }
  ds <- ds/dsn + .1
  return(ds)
}
