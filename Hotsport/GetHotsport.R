USAGE <- function(){
         message("USAGE:\n\tRscript GetHotsport.R <snpfile>  <cutoff>  <cf.length>  <tolerance.length>  <to_up_len> <output>\n")
         message("\tsnpfile must have two cols(chrom posi) and header\n")
         message("\tcutoff is the  max Pvalue usually use 0.05\n")
         message("\tcf.length is at least  how many pvalue <= cutoff\n")
         message("\ttolerance.lenght tolerance for False within the hotsport region\n")
         message("\tto_up_len is the minimum value for continuous True within the region\n")
         message("\toutpath should have filepath and filename \n")
}

RegionAsBed <- function(marker, cf.length, tolerance.length,to_up_len ,chrom) {
        r <- rle(marker)
        if(is.null(tolerance.length)){
                end <- with(r,cumsum(lengths)[values & lengths>=cf.length])
                start <- end - with(r,lengths[values & lengths>=cf.length]) + 1
                if(length(start) == 0 || length(end) == 0){
                        df <- data.frame(chrom=character(), start=integer(), end=integer())
                } else {
                        df <- data.frame(chrom, start, end)
                }
        } else {
                end<-cumsum(r$lengths)
                start<- end - r$lengths + 1
                df.tmp <- data.frame(start = start  , end = end , value = r$values)
                df.tmp$len <- df.tmp$end - df.tmp$start + 1
                df.tolerance <- df.tmp[!df.tmp$value & df.tmp$len <= tolerance.length,]
                tolerance_rowname <- as.numeric(rownames(df.tolerance))
                tolerance_data <- data.frame(tolerance_rowname = tolerance_rowname)
                tolerance_data$tolerance_up <- tolerance_rowname - 1
                tolerance_data$tolerance_down <- tolerance_rowname + 1
                tolerance_data_sel <- tolerance_data[with(tolerance_data, tolerance_up > 0 & tolerance_down <= nrow(df.tmp)),]
                tolerance_mark <- tolerance_data_sel$tolerance_rowname[df.tmp$len[tolerance_data_sel$tolerance_up] >= to_up_len & df.tmp$len[tolerance_data_sel$tolerance_down] >= to_up_len]
                dt.to <- data.table(df.tmp[tolerance_mark,])
		dt.to <- dt.to[complete.cases(dt.to),]
                if (nrow(dt.to)!=0) {
                        dt.to$id <- 1:nrow(dt.to)
                        dt.to.index <- dt.to[,.(index = start:end), by = id]
                        marker[dt.to.index$index] <- TRUE
                }

                r <- rle(marker)
                end <- with(r,cumsum(lengths)[values & lengths>=cf.length])
                start <- end - with(r,lengths[values & lengths>=cf.length]) + 1
                if(length(start) == 0 || length(end) == 0){
                        df <- data.frame(chrom=character(), start=integer(), end=integer())
                } else {
                        df <- data.frame(chrom, start, end)
                }
        }
    return(df)
}



args <- commandArgs(T)
USAGE()
library(data.table)
library(zoo)
snpfile <- fread(args[1])
colnames(snpfile) <- c("chrom","posi")
snpfile$chrom <- as.character(snpfile$chrom)
cutoff <- as.numeric(args[2])
cflength <- as.integer(args[3])
tolerancelength  <-as.integer(args[4])
to_up_len <- as.integer(args[5])
outpath <- args[6]
chrom <- unique(snpfile$chrom)
result <- data.frame()
for(chr in chrom){
    message(chr,"  is running  ",date())
    split <- snpfile[snpfile$chrom == chr,]
    distance <- diff(split$posi)
    split <- split[2:nrow(split),]
    split$distance <- distance
    split$id <- 1:nrow(split)
    lamta  <- 1/mean(split$distance)
    split[,Pvalue:= (1 - exp(-lamta * distance)),by=id]
    split$id <- NULL
    split <- as.data.frame(split)
    m <- split$Pvalue < cutoff
    hotsport <- RegionAsBed(m,cflength,tolerancelength,to_up_len,chr)
    hotsport$start <- split$posi[hotsport$start]
    hotsport$end <- split$posi[hotsport$end]
    result <- rbind(result,hotsport)
}
result$length <- result$end - result$start + 1
write.table(result,file=outpath,quote=F,col.names=T,row.names=F)
