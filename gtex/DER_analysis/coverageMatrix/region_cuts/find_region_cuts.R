## Load required libraries
library('derfinder')
library('GenomicRanges')
library('RColorBrewer')
library('BiocParallel')
library('devtools')

chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
cuts <- seq(from = 0.2, to = 7, by = 0.1)

## Parallel environment to use
bp <- SnowParam(workers = 25, outfile = Sys.getenv('SGE_STDERR_PATH'))

## Find the regions for all the chromosomes given a specific cutoff
getRegChr <- function(cutoff, chr, meanCov) {
    suppressPackageStartupMessages(library('derfinder'))
    suppressPackageStartupMessages(library('IRanges'))
    message(paste(Sys.time(), 'processing', chr, 'with cutoff', cutoff))
    regs <- findRegions(position = Rle(TRUE, length(meanCov$coverage[[1]])), fstats = meanCov$coverage[[1]], chr = chr, maxClusterGap = 300L, cutoff = cutoff, verbose = FALSE)
    if(is.null(regs)) {
        message(paste(Sys.time(), 'found no regions for', chr, 'with cutoff', cutoff))
        regs <- GRanges()
        return(regs)
    }
    else if(!is(regs, 'GRanges')) {
        message(paste(Sys.time(), 'processing failed for', chr, 'with cutoff', cutoff))
        regs <- GRanges()
        return(regs)
    }
    names(regs) <- NULL
    return(regs)
}
getRegs <- function(chr, chrlen, cutoffs, param) {
    message(paste(Sys.time(), 'loading meanCov for', chr))
    meanCov <- loadCoverage('/dcl01/leek/data/gtex_work/gtex_mean_coverage.bw', chr, chrlen = chrlen)
    res <- bplapply(cutoffs, getRegChr, chr = chr, meanCov = meanCov, BPPARAM = param)
    names(res) <- as.character(cutoffs)
    return(res)
}

region_cuts_raw <- mapply(getRegs, chrInfo$chr, chrInfo$length, MoreArgs = list(cutoffs = cuts, param = bp), SIMPLIFY = FALSE)
names(region_cuts_raw) <- chrInfo$chr

message(paste(Sys.time(), 'saving region_cuts_raw'))
save(region_cuts_raw, file = 'region_cuts_raw.Rdata')


region_cuts <- lapply(as.character(cuts), function(cutoff, chromosomes = chrInfo$chr) {
    regs <- lapply(chromosomes, function(chr) { region_cuts_raw[[chr]][[cutoff]] })
    names(regs) <- chromosomes
    regs <- unlist(GRangesList(regs))
    
    seqlen <- chrInfo$length[match(names(seqlengths(regs)), chrInfo$chr)]
    names(seqlen) <- names(seqlengths(regs))
    seqlengths(regs) <- seqlen
    
    return(regs)
})
names(region_cuts) <- as.character(cuts)
message(paste(Sys.time(), 'saving region_cuts'))
save(region_cuts, file = 'region_cuts.Rdata')

## Extract information of interest
info <- lapply(region_cuts, function(regs) {
    res <- list(mean = mean(width(regs)), n = length(regs), quantile = quantile(width(regs), seq(0, 1, by = 0.1)), sd = sd(width(regs)))
    return(res)
})

n <- unlist(lapply(info, '[[', 'n'))
mean <- unlist(lapply(info, '[[', 'mean'))
sd <- unlist(lapply(info, '[[', 'sd'))
quantile <- do.call(rbind, lapply(info, '[[', 'quantile'))

## Arrange information into a single data frame
regInfo <- data.frame(cutoff = cuts, n = n, mean = mean, sd = sd, stringsAsFactors = FALSE)
regInfo <- cbind(regInfo, quantile)
save(regInfo, file = 'regInfo.Rdata')


## Make plots exploring the relationship between the cutoff and the number of regions
pdf(file = 'regInfo-plots.pdf', width = 10, height = 10)
plot(x = regInfo$cutoff, y = regInfo$n, xlab = 'Cutoff', ylab = 'Number of regions', type = 'o')
plot(x = regInfo$cutoff, y = regInfo$mean, xlab = 'Cutoff', ylab = 'Mean region width', type = 'o')
plot(x = regInfo$cutoff, y = regInfo$sd, xlab = 'Cutoff', ylab = 'SD of region width', type = 'o')
plot(x = regInfo$sd, y = regInfo$mean, xlab = 'SD of region width', ylab = 'Mean region width', type = 'o')
plot(x = regInfo$cutoff, y = regInfo$mean * regInfo$n, xlab = 'Cutoff', ylab = 'Total bp in regions', type = 'o')

percent <- regInfo$mean * regInfo$n / sum(as.numeric(chrInfo$length)) * 100
names(percent) <- cuts
plot(x = regInfo$cutoff, y = percent, xlab = 'Cutoff', ylab = 'Percent of genome in regions', type = 'o')
percent

matplot(x = regInfo$cutoff, data.frame(m_minusSD = regInfo$mean - regInfo$sd, m = regInfo$mean, m_plusSD = regInfo$mean + regInfo$sd), type = 'o', lty = 1, lwd = 2, xlab = 'Cutoff', ylab = 'Mean region width (+- SD)', col = brewer.pal(n = 3, 'Set1'), pch = 21)
abline(h = 0, col = 'black')


colors <- brewer.pal(11, 'PuOr')
colors[6] <- 'magenta'

matplot(x = regInfo$cutoff, log2(as.matrix(regInfo[, 5:15]) + 1), type = 'o', xlab = 'Cutoff', col = colors, lty = 1, lwd = 2, ylab = 'log2(length + 1) quantiles by 10% increments', pch = 21)

matplot(x = regInfo$cutoff, log10(as.matrix(regInfo[, 5:15]) + 1), type = 'o', xlab = 'Cutoff', col = colors, lty = 1, lwd = 2, ylab = 'log10(length + 1) quantiles by 10% increments', pch = 21)
dev.off()


## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
session_info()
