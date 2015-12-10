## Load required libraries
library('derfinder')
library('GenomicRanges')
library('RColorBrewer')
library('BiocParallel')
library('devtools')

chrs <- paste0('chr', c(1:22, 'X', 'Y'))
cuts <- seq(from = 0.1, to = 5, by = 0.1)

## Parallel environment to use
bp <- SnowParam(workers = 10, outfile = Sys.getenv('SGE_STDERR_PATH'))

## Load the mean coverage for all the chromosomes
meanCovs <- bplapply(chrs, function(chr) {
    library('derfinder')
    message(paste(Sys.time(), 'processing', chr))
    meanCov <- loadCoverage('/dcl01/leek/data/gtex_work/gtex_mean_coverage.bw', chr)
    file.path <- file.path(Sys.getenv('TMPDIR'), paste0(chr, '-meanCov.Rdata'))
    message(paste(Sys.time(), 'saving meanCov for', chr, 'at', file.path))
    save(meanCov, file = file.path)
    return(file.path)
}, BPPARAM = bp)
names(meanCovsFiles) <- chrs
#save(meanCovs, file = 'meanCovs.Rdata')
#> print(object.size(meanCovs), units = 'Gb')
#27.9 Gb

## Find the regions for all the chromosomes given a specific cutoff
getRegs <- function(chr, meanCovFile, cutoff) {
    library('derfinder')
    library('IRanges')
    message(paste(Sys.time(), 'loading meanCov file', meanCovFile))
    load(meanCovFile)
    message(paste(Sys.time(), 'processing', chr, 'with cutoff', cutoff))
    regs <- findRegions(position = Rle(TRUE, length(meanCov$coverage[[1]])), fstats = meanCov$coverage[[1]], chr = chr, maxClusterGap = 300L, cutoff = cutoff, verbose = FALSE)
    return(regs)
}

region_cuts <- lapply(cuts, function(cutoff) {
    message(paste(Sys.time(), 'working with cutoff', cutoff))
    regs <- bpmapply(getRegs, chrs, meanCovsFiles, MoreArgs = list(cutoff = cutoff), BPPARAM = bp)
    regs <- unlist(GRangesList(regs))
    return(regs)
})
save(region_cuts, file = 'region_cuts.Rdata')

## Extract information of interest
info <- lapply(region_cuts, function(regs) {
    res <- list(mean = mean(width(regs)), n = length(regs), quantile = quantile(width(regs), seq(0, 1, by = 0.1)), sd = sd(width(regs)))
    return(res)
})

n <- unlist(lapply(info, '[[', 'n'))
mean <- unlist(lapply(info, '[[', 'mean'))
sd <- unlist(lapply(info, '[[', 'mean'))
quantile <- do.call(rbind, lapply(info, '[[', 'quantile'))

## Arrange information into a single data frame
regInfo <- data.frame(cutoff = cuts, n = n, mean = mean, sd = sd, stringsAsFactors = FALSE)
regInfo <- cbind(regInfo, quantile)
save(regInfo, file = 'regInfo.Rdata')


## Make plots exploring the relationship between the cutoff and the number of regions
pdf(file = 'regInfo-plots.pdf', width = 10, height = 10)
plot(x = regInfo$cutoff, y = regInfo$n, xlab = 'Cutoff', ylab = 'Number of regions')
plot(x = regInfo$cutoff, y = regInfo$mean, xlab = 'Cutoff', ylab = 'Mean region width')
plot(x = regInfo$cutoff, y = regInfo$sd, xlab = 'Cutoff', ylab = 'SD of region width')

matplot(x = regionInfo$cutoff, data.frame(m_minusSD = regInfo$mean - regInfo$sd, m = regInfo$mean, m_plusSD = regInfo$mean + regInfo$sd), lty = 1, lwd = 2, xlab = 'Cutoff', ylab = 'Mean region width (+- SD)', col = brewer.pal(n = 3, 'Set1'))


colors <- brewer.pal(11, 'PuOr')
colors[6] <- 'magenta'

matplot(x = regInfo$cutoff, log2(as.matrix(regInfo[, 5:15]) + 1), type = 'l', xlab = 'Cutoff', col = colors, lty = 1, lwd = 2, ylab = 'log2(length + 1)')

matplot(x = regInfo$cutoff, log10(as.matrix(regInfo[, 5:15]) + 1), type = 'l', xlab = 'Cutoff', col = colors, lty = 1, lwd = 2, ylab = 'log10(length + 1)')
dev.off()


## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
session_info()
