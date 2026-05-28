library(dartR.base)
test<-bandicoot.gl
validpops<-popNames(test)[1:3]
test<-gl.keep.pop(test, validpops)
validinds<-indNames(test_subset)[1:15]
validloci<-locNames(test_subset)[1:200]
test_subset<-gl.keep.ind(test, validinds)
rm(test);gc()
test_subset<-gl.keep.loc(test_subset,  loc.list=validloci)
test_subset@other$loc.metrics$TrimmedSequence <-
    test_subset@other$loc.metrics$AlleleSequence
outdf<-data.frame()
for (i in indNames(test_subset)){
    print(i)
    latlon<-(gl.keep.ind(test_subset,i))[1] 
    latlon<-latlon@other$ind.metrics
    latlon$indName<-i
    outdf<-rbind.data.frame(outdf,latlon)
}
npoints<-nrow(outdf)
minlon <- -72
maxlon <- -70
minlat <- -33.95
maxlat <- -32.03
outdf$lat <- runif(npoints, min = minlat, max = maxlat)
outdf$lon <- runif(npoints, min = minlon, max = maxlon)
envardf<-outdf[,c("id", "lat","lon")]
for (i in 1:5){
    meanr <- runif(1, -10, 100)
    sdrex <- runif(1, -3, 3)
    sdr<-10^sdrex
    print(sdr)
    varname <- paste0("envar_", i)
    values <- rnorm(
        nrow(envardf),
        mean = meanr,
        sd = sdr
    )
    envardf[,varname] <- values
}

save_test_path<-function(x = NULL){
    if (is.null(x)) {
        outvect <- paste(getwd(), "test-data",sep="/")
        return(outvect)
    }
    return(paste(getwd(), "test-data", x,sep="/"))
    
}


write.csv(envardf, file=save_test_path("envars.csv"))
write.csv(outdf,file=save_test_path("indmetrics.csv"))

path_plink<-paste0( getwd(),"/plink")
keyword<-"test-data"
savepath<-save_test_path()
test_subset<-gl.impute(test_subset,method="random")

saveRDS(test_subset, file = save_test_path(paste0(keyword, ".rds")))
gl2structure(test_subset, outfile = paste0(keyword, ".str"),        outpath = savepath, ploidy = 2)
#gl2structure(test_subset, outfile = paste0(keyword, "_ploid1.str"), outpath = savepath, ploidy = 1)
gl2vcf(test_subset,     outfile = basename(keyword),              outpath = savepath, plink.bin.path = path_plink)
gl2fasta(test_subset,  outfile = paste0(keyword, ".fasta"),      outpath = savepath, method = 3)
gl2plink(test_subset,     outfile = keyword,                        outpath = savepath, plink.bin.path = path_plink, bed.files = TRUE)
gl2gds(test_subset,       outfile = paste0(keyword, ".gds"),        outpath = savepath)