library('shiny')
library('DT')
library('markdown')

# cp ../metadata_web/meta_web_sra.Rdata .
# mkdir www
# cp ../genes/ucsc-knowngene-hg38-exons.Rdata www/
# cp ../genes/ucsc-knowngene-hg38-genes-bp-length.Rdata www/
load('meta_web_sra.Rdata')

meta_web$species <- factor(meta_web$species)
colnames(meta_web)[colnames(meta_web) == 'number_samples'] <- 'number of samples'
colnames(meta_web)[colnames(meta_web) == 'rse_gene'] <- 'RSE gene'
colnames(meta_web)[colnames(meta_web) == 'rse_exon'] <- 'RSE exon'
colnames(meta_web)[colnames(meta_web) == 'counts_gene'] <- 'counts gene'
colnames(meta_web)[colnames(meta_web) == 'counts_exon'] <- 'counts exon'
colnames(meta_web)[colnames(meta_web) == 'files_info'] <- 'files info'

shinyServer(function(input, output, session) {
    output$metadata = DT::renderDataTable(
        meta_web[, - which(colnames(meta_web) %in% c('genes', 'exons'))],
        escape = which(colnames(meta_web) %in% c('number of samples', 'species',
            'abstract')),
        style = 'bootstrap', rownames = FALSE, filter = 'top',
        options = list(
            columnDefs = list(
                list(className = 'dt-center', targets = 1)
            ),
            pageLength = 10,
            lengthMenu = c(5, 10, 25, 50, 100, nrow(meta_web))
        )
    )
    output$popular = DT::renderDataTable(
        meta_web[meta_web[, 2] > 500, - which(colnames(meta_web) %in%
            c('genes', 'exons'))],
        escape = which(colnames(meta_web) %in% c('number of samples', 'species',
            'abstract')),
        style = 'bootstrap', rownames = FALSE, filter = 'top',
        options = list(
            columnDefs = list(
                list(className = 'dt-center', targets = 1)
            ),
            pageLength = 10,
            lengthMenu = c(5, 10, 25, 50, 100, nrow(meta_web))
        )
    )
})
