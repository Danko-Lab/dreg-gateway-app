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
gpu_cores <- as.integer(args[6])
if (is.na(gpu_cores)) gpu_cores <- 1;

out.log  <- paste(outfile, "dREG.log", sep=".")

cat("Bigwig(plus):", ps_plus_path, "\n");
cat("Bigwig(minus):", ps_minus_path, "\n");
cat("Output:", outfile, "\n");
cat("dREG model:", args[4], "\n");
cat("ncores:", ncores, "\n");
cat("GPU:", gpu_cores, "\n");

cat("1) -------- Checking bigWig files\n");
b1 <- check_bigwig(ps_plus_path, strand="+", out.file = out.log );
b2 <- check_bigwig(ps_minus_path, strand="-", out.file = out.log );
if( !b1 || !b2 )
{
    cat("Warning: bigWig files maybe not meet the requirements. See dREG requirement in https://github.com/Danko-Lab/dREG");
}

cat("2) -------- Loading model\n");
## Load the dRGE model including two ojects 'asvm' and 'gdm'.
## Do this before loading ps_plus_path, just in case those are saved in the model file.
## Should have (by default) gdm and asvm.
load(args[4])

## Now scan all positions in the genome ...
cat("3) -------- Peak calling\n");
t.run <- system.time( r <- peak_calling( asvm, gdm, ps_plus_path, ps_minus_path, cpu_cores=ncores, gpu_cores=gpu_cores, use_rgtsvm=(gpu_cores!=0) ) )

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

	file.tmp2 <- tempfile(fileext=".bed");
	system(paste( "sort-bed ", file.tmp, ">", file.tmp2) );
	system(paste( "bedGraphToBigWig ", file.tmp2, out.chrom.info, out_file,  sep=" ") );

	unlink(file.tmp)
	unlink(file.tmp2)
}

make.chrom.info <- function( file.bw.plus, file.bw.minus, file.chrom.info)
{
	bw.plus <- load.bigWig(file.bw.plus);
 	bw.minus <- load.bigWig(file.bw.minus);

 	chrom <- rbind( cbind( bw.plus$chroms, bw.plus$chromSizes), cbind( bw.minus$chroms, bw.minus$chromSizes) );
   	chr.size <- unlist( lapply( unique(chrom[,1]), function(chr){max( as.numeric( chrom[which(chrom[,1]==chr),2])) } ) );

	df.bed <- data.frame( V1=unique(chrom[,1]), V2=chr.size );
	write.table( df.bed, file=file.chrom.info, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t");
}

make_summary<-function( r, out.log )
{
	cat("[summary]\n", file=out.log);
	cat("Running time:", t.run[3], "\n", file=out.log, append=TRUE);
	cat("Informative loci:", NROW(r$infp_bed), "\n", file=out.log, append=TRUE);
	cat("Broad peaks:", NROW(r$peak_broad), "\n", file=out.log, append=TRUE);
	cat("Threshold:", r$min_score, "\n", file=out.log, append=TRUE);
	cat("Narrow peaks:", NROW(r$raw_peak), "\n", file=out.log, append=TRUE);
	cat("Peaks(p<=0.05 ):", r$peak_sum$adjust.none.0.05, "\n", file=out.log, append=TRUE);
	cat("Peaks(p<=0.05 with FDR correction):", r$peak_sum$adjust.fdr.0.05, "\n", file=out.log, append=TRUE);
	cat("Peaks(p<=0.05 with Bonferroni correction):", r$peak_sum$adjust.bonferroni.0.05, "\n", file=out.log, append=TRUE);
	cat("Peaks(width<50):", r$peak_sum$peak.narrow100, "\n", file=out.log, append=TRUE);
	cat("Peaks(width<50 and p<=0.05):", r$peak_sum$peak.narrow100.sig, "\n", file=out.log, append=TRUE);

}

cat("4) -------- outputting result\n");

make_summary( r, out.log );

out.infp.bed <- paste(outfile, "dREG.infp.bed", sep=".")
out.infp.bw  <- paste(outfile, "dREG.infp.bw", sep=".")
out.score.bed<- paste(outfile, "dREG.peak.score.bed", sep=".")
out.score.bw <- paste(outfile, "dREG.peak.score.bw", sep=".")
out.prob.bed <- paste(outfile, "dREG.peak.prob.bed", sep=".")
out.prob.bw  <- paste(outfile, "dREG.peak.prob.bw", sep=".")
out.full.bed <- paste(outfile, "dREG.peak.full.bed", sep=".")
out.chrom.info <- paste(outfile, "chrom.info", sep=".")

make.chrom.info(ps_plus_path, ps_minus_path, out.chrom.info);

make_index_gz( r$infp_bed, out.infp.bed );
make_bw( r$infp_bed[,1:4], out.infp.bw );

make_index_gz( r$peak_bed, out.full.bed );

make_index_gz( r$peak_bed[,c(1:4)], out.score.bed );
make_bw( r$peak_bed[,c(1:4)], out.score.bw );

r$peak_bed[,5] <- 1 - r$peak_bed[,5];
make_index_gz( r$peak_bed[,c(1,2,3,5)], out.prob.bed );
make_bw( r$peak_bed[,c(1,2,3,5)], out.prob.bw );

system( paste("tar -cvzf ", outfile, ".tar.gz", " ", outfile, ".dREG.*", sep="") );
cat("Result:", paste(outfile, ".tar.gz", sep=""), "\n");

cat("Done!\n");
