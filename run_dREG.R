require(dREG)
options("scipen"=100, "digits"=4)

## Process command arguments
args <- commandArgs(trailingOnly=TRUE)

## Read arguments from thwe web page.
ps_plus_path  <- args[1]
ps_minus_path <- args[2]
## Read arguments from default parameters in run_dREG.sh
outfile <- args[3]

ncores <- as.integer(args[5])
if (is.na(ncores)) ncores <- 1;
gpu_cores <- toupper(as.integer(args[6]))
if (is.na(ncores)) ncores <- 1;

cat("Bigwig(plus):", ps_plus_path, "\n");
cat("Bigwig(minus):", ps_minus_path, "\n");
cat("Output:", outfile, "\n");
cat("dREG model:", args[4], "\n");
cat("ncores:", ncores, "\n");
cat("GPU:", gpu_cores, "\n");


cat("1) -------- Loading model\n");
## Load the dRGE model including two ojects 'asvm' and 'gdm'.
## Do this before loading ps_plus_path, just in case those are saved in the model file.
## Should have (by default) gdm and asvm.
load(args[4])

## Now scan all positions in the genome ...
cat("2) -------- Peak calling\n");
r <- peak_calling( svm, gdm, ps_plus_path, ps_minus_path, cpu_cores=ncores, gpu_cores=gpu_cores, use_rgtsvm=(gpu_cores!=0) )

make_index_gz<-function( df_bed, out_file)
{
	file.tmp <- tempfile(fileext=".bed");
	write.table( df_bed, file=file.tmp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
	system(paste( "sort-bed ", file.tmp, " | bgzip -f > ", out_file, ".gz", sep="") );
	system(paste( "tabix -f -p bed ", out_file, ".gz", sep="") );
	unlink(file.tmp)
}

make_bw<-function( df_bed, out_file)
{
	## decrease the size
	df_bed[,4] <- round(df_bed[,4], 2);

	file.tmp <- tempfile(fileext=".bed");
	write.table( df_bed, file=file.tmp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");

	system(paste( "sort-bed ", file.tmp, " | bedGraphToBigWig -  ", out.chrom.info, out_file,  sep=" ") );
	unlink(file.tmp)
}

make.chrom.info <- function( file.bw.plus, file.bw.minus, file.chrom.info)
{
	bw.plus <- load.bigWig(file.bw.plus);
    bw.minus <- load.bigWig(file.bw.minus);

    chrom <- rbind( cbind( bw.plus$chroms, bw.plus$chromSizes), cbind( bw.minus$chroms, bw.minus$chromSizes) );
    chr.size <- unlist( lapply( unique(chrom[,1]), function(chr){max( as.numeric( chrom[which(chrom[,1]==chr),2])) } ) );

	df.bed <- data.frame( V1=unique(chrom[,1]), V2=chr.size );
	write.table( df_bed, file=file.chrom.info, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
}


out.infp.bed <- paste(outfile, "dREG.infp.bed", sep=".")
out.infp.bw  <- paste(outfile, "dREG.infp.bw", sep=".")
out.score.bed<- paste(outfile, "dREG.peak.score.bed", sep=".")
out.score.bw <- paste(outfile, "dREG.peak.score.bw", sep=".")
out.prob.bed <- paste(outfile, "dREG.peak.prob.bed", sep=".")
out.prob.bw  <- paste(outfile, "dREG.peak.prob.bw", sep=".")
out.full.bed <- paste(outfile, "dREG.peak.full.bed", sep=".")
out.chrom.info <- paste(outfile, "chrom.info", sep=".")

cat("3) -------- outputting result\n");
make.chrom.info(ps_plus_path, ps_minus_path, out.chrom.info);

make_index_gz( r$peak_bed, out.full.bed );

make_index_gz( r$infp_bed, out.info.bed );
make_bw( r$infp_bed, out.info.bw );

make_index_gz( r$peak_bed[,c(1:4)], out.score.bed );
make_bw( r$peak_bed[,c(1:4)], out.score.bw );

r$peak_bed[,5] <- 1 - r$peak_bed[,5];
make_index_gz( r$peak_bed[,c(1,2,3,5)], out.prob.bed );
make_bw( r$peak_bed[,c(1,2,3,5)], out.prob.bw );

system( paste("tar -cvzf ", outfile, ".dREG.tar.gz", " ", outfile, ".dREG", sep="") );
cat("Result:", paste(outfile, ".dREG.tar.gz", sep=""), "\n");
cat("Done!\n");

