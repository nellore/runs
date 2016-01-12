# Load libraries
library('RSQLite')
library('magrittr')
library('SRAdb')
library('stringr')

# Define functions
"%p%"  <- function(x, y) paste0(x, y)


sqlfile <- file.path('.', 'SRAmetadb.sqlite')

# Download sql if necessary:
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()

# Create connection
sra_con <- dbConnect(SQLite(),sqlfile)

# Query database
query <- "SELECT 
run_accession,
sample_accession,
experiment_accession,
study_accession,
submission_accession,
sra_ID,                                        
run_ID,
run_alias,                    
run_date,                     
updated_date,
spots,                        
bases,
run_center,                   
experiment_name,                         
run_attribute,                
experiment_ID,
experiment_alias,             
experiment_title,
study_name,
sample_name,
design_description,
library_name,
library_strategy,
library_source,
library_selection,
library_layout,
library_construction_protocol,
read_spec,
platform,
instrument_model,
platform_parameters,
experiment_url_link,
experiment_attribute,
sample_ID,
sample_alias,
taxon_id,
common_name,
description,
sample_url_link,
sample_attribute,
study_ID,
study_alias,
study_title,
study_type,
study_abstract,
center_project_name,
study_description,
study_url_link,
study_attribute,
related_studies,
primary_study,
submission_ID,
submission_comment,
submission_center,
submission_lab,
submission_date,
sradb_updated
FROM sra
WHERE platform = 'ILLUMINA' AND
library_strategy = 'RNA-Seq' AND
taxon_id = 9606;"

selected <- dbGetQuery(sra_con, query)
query <- paste0("SELECT * FROM fastq WHERE run_accession IN ('", 
                paste(selected$run_accession, collapse="', '"), "')")
fastq_namefiles <- dbGetQuery(sra_con, query)

# Disconnect from db
dbDisconnect(sra_con)

metadata <- merge(selected, fastq_namefiles, by = "run_accession")

metadata$date_download <- rep(Sys.time(),nrow(metadata)) %>% as.character()

# Rename column name for file data
names(metadata)[names(metadata) == 'sradb_updated.y'] <- 'sradb_updated_file'

get.fastq.urls <- function(df) {
    url <- rep(NA, length(df$run_accession))
    for(i in 1:length(df$run_accession)) {
        run <- df$run_accession[i]
        filename <- df$file_name[i]
        if(nchar(run) < 10) {
            url[i] <- file.path('ftp://ftp.sra.ebi.ac.uk/vol1/fastq',
                                substring(run, 1, 6), run, filename)
        } else {
            dir2 <- paste( c(rep(x='0', 12-nchar(run)), substring(run, 10, 
                                                                  nchar(run))), collapse = '' )
            url[i] <- file.path('ftp://ftp.sra.ebi.ac.uk/vol1/fastq', 
                                substring(run, 1, 6), dir2, run, filename)
        }  
    }
    return(url)
}


metadata$URL <- get.fastq.urls(metadata)
metadata$layout <- unlist(lapply(strsplit(metadata$library_layout, " - "), `[`, 1))

# Get all PAIRED-end and SINGLE runs with 2 associated files (forward and reverse)

## PAIRED _1
paired_and_forward <- grep("_1", metadata$URL)
## PAIRED _2
paired_and_reverse <- grep("_2", metadata$URL)

equal <- setequal(metadata$run_accession[paired_and_forward], 
                  metadata$run_accession[paired_and_reverse])
cat("N°studies:(paired-end AND '_1')\t", length(paired_and_forward), "\n",
    "N°studies:(paired-end AND '_2')\t", length(paired_and_reverse), "\n",
    "Are run_accessions equal for paired AND '_1' and paired AND '_2'?:\n", equal)
diff1 <- setdiff(metadata$run_accession[paired_and_forward], metadata$run_accession[paired_and_reverse])
diff2 <- setdiff(metadata$run_accession[paired_and_reverse], metadata$run_accession[paired_and_forward])
if (!equal){
    cat("Differences forward-reverse:\n")
    print(diff1)
    cat("Differences reverse-forward:\n")
    print(diff2)
}

file_names_diff2 <- subset(metadata, run_accession == diff2)[,"file_name"]

if(!(grepl("_", file_names_diff2[1]) & grepl("_", file_names_diff2[2]))){
    failed_url <- unlist(strsplit(subset(metadata, run_accession == diff2)[1,"URL"], "/"))
    positions <- rownames(subset(metadata, run_accession == diff2))
    pos_file <- grep('.fastq.gz', failed_url)
    failed_url
    f_file <- failed_url -> r_file
    f_file[pos_file[1]] <- as.character(diff2) %p% "_2.fastq.gz"
    r_file[pos_file[1]] <- as.character(diff2) %p% "_1.fastq.gz"
    f_file <- paste(f_file, collapse = "/")
    r_file <- paste(r_file, collapse = "/")
    metadata[positions[1],c("URL")] <- f_file
    metadata[positions[2],c("URL")] <- r_file
    metadata[positions[2],c("md5")] <- '0'
}

## PAIRED _1
paired_and_forward <- grep("_1", metadata$URL)
## PAIRED _2
paired_and_reverse <- grep("_2", metadata$URL)

equal <- setequal(metadata$run_accession[paired_and_forward], 
                  metadata$run_accession[paired_and_reverse])
cat("N°studies:(paired-end AND '_1')\t", length(paired_and_forward), "\n",
    "N°studies:(paired-end AND '_2')\t", length(paired_and_reverse), "\n",
    "Are run_accessions equal for paired AND '_1' and paired AND '_2'?:\n", equal)
diff1 <- setdiff(metadata$run_accession[paired_and_forward], metadata$run_accession[paired_and_reverse])
diff2 <- setdiff(metadata$run_accession[paired_and_reverse], metadata$run_accession[paired_and_forward])
if (!equal){
    cat("Differences forward-reverse:\n")
    print(diff1)
    cat("Differences reverse-forward:\n")
    print(diff2)
}



metadata_twofiles_single_F <- table(metadata[paired_and_forward,]$layout)
metadata_twofiles_single_R <- table(metadata[paired_and_reverse,]$layout)
metadata_twofiles_single_F[3] <- "    <-- Runs with 2 files | F | _1"
metadata_twofiles_single_R[3] <- "    <-- Runs with 2 files | R | _2"

print(metadata_twofiles_single_F)
print(metadata_twofiles_single_R)

# Make sure single runs have F and R associated files 
metadata_twofiles_single_F_layout <- subset(metadata[paired_and_forward,], layout == "SINGLE")
nrow(metadata_twofiles_single_F_layout)
metadata_twofiles_single_R_layout <- subset(metadata[paired_and_reverse,], layout == "SINGLE")
nrow(metadata_twofiles_single_R_layout)
setdiff(metadata_twofiles_single_F_layout$run_accession, 
        metadata_twofiles_single_R_layout$run_accession)

# Normalize metadata table
ids_filename <- c("run_accession", "study_accession", "sample_accession", 
                  "experiment_accession", "submission_accession", "file_name", "URL", "md5")

PF <- metadata[paired_and_forward, ids_filename]
PR <- metadata[paired_and_reverse, ids_filename]

PF <- cbind(PF, paired_and_forward)
PR <- cbind(PR, paired_and_reverse)
paired_experiments <- merge(PF, PR, by = head(ids_filename, -3))

colnames(paired_experiments) <- c("run_accession", "study_accession", "sample_accession", "experiment_accession",
                                  "submission_accession", "file_name.f", "URL.f", "md5.f","paired_and_forward", "file_name.r",
                                  "URL.r", "md5.r","paired_and_reverse")

# Assign urls to paired runs (reverse)
reverse_url <- rep("NA", nrow(metadata))
reverse_url[paired_experiments$paired_and_forward] <- paired_experiments$URL.r
metadata$URL_R <- reverse_url

# Assign md5 to paired runs (reverse)
reverse_md5 <- rep("NA", nrow(metadata))
reverse_md5[paired_experiments$paired_and_forward] <- paired_experiments$md5.r
metadata$md5_R <- reverse_md5

# Remove reverse paired-end runs from metadata table
metadata <- metadata[-paired_experiments$paired_and_reverse, ]


# Identify duplicated run_accessions
print(length(metadata$run_accession))
print(length(unique(metadata$run_accession)))
#metadata[duplicated(metadata$run_accession),]$run_accession
duplicated_accessions <- metadata[duplicated(metadata$run_accession),]
i <- metadata[,1] %in% duplicated_accessions$run_accession
duplicated_accessions <-metadata[i,]
print(duplicated_accessions[1:4,])

# Remove duplicated paired accessions associated with 3 files
# Keep those with two files
# Example:
# (ERR033743.fastq.gz, ERR033743_1.fastq.gz, and ERR033743_2.fastq.gz)
# Keep ERR033743_1.fastq.gz, and ERR033743_2.fastq.gz pair


paired_with_F <- intersect(grep("PAIRED", metadata$layout), grep("_1", metadata$URL))
paired_with_F_and_R <- intersect(paired_with_F, grep("_2", metadata$URL_R))


duplicated_paired_wo_R <- setdiff(grep("PAIRED", metadata$layout), paired_with_F_and_R)
print("Number of duplicated runs associated with 3 files")
print(nrow(metadata[duplicated_paired_wo_R,]))

metadata <- metadata[-duplicated_paired_wo_R,]

# Corroborations 

print("Make sure metadata does not have duplicated runs")
print(nrow(metadata))
print(length(unique(metadata$run_accession)))

print("Make sure all paired experiments have 2 files")
paired_1 <- intersect(grep("PAIRED", metadata$layout), grep("_1", metadata$URL))
paired_2 <- intersect(grep("PAIRED", metadata$layout), grep("_2", metadata$URL_R))
print(setdiff(paired_1, paired_2))
print(setdiff(paired_2, paired_1))

print("Are there single runs with 2 associated files?")
single_1 <- intersect(grep("SINGLE", metadata$layout), grep("_1", metadata$URL))
single_2 <- intersect(grep("SINGLE", metadata$layout), grep("_2", metadata$URL_R))
length(intersect(single_1, single_2))

print("Are there single runs with files ending with _2 in URL?")
length(intersect(grep("SINGLE", metadata$layout), grep("_2", metadata$URL)))

print("Are there single runs with files ending with _1 in URL_R?")
length(intersect(grep("SINGLE", metadata$layout), grep("_1", metadata$URL_R)))


search.field <- function(column, field) {
    str_split(column, "\\|\\|") %>% 
        lapply(str_trim) %>%
        lapply(str_extract, regex(field %p% ":.*", ignore_case = TRUE)) %>%
        lapply(na.omit) %>%
        lapply(as.vector) %>% ifelse(. == 'character(0)', 'NA', .) %>%
        lapply(`[[`, 1) %>%
        unlist() %>% 
        str_replace_all(regex(field %p% ": ", ignore_case = TRUE), "") %>%
        str_to_upper()
}

# Get cell type
metadata$cell_type <- search.field(metadata$sample_attribute, "cell type")

# Get tissue
metadata$tissue <- search.field(metadata$sample_attribute, "tissue")

# Get cell line
metadata$cell_line <- search.field(metadata$sample_attribute, "cell line")

# Get strain
metadata$strain <- search.field(metadata$sample_attribute, "Strain")

# Get age
metadata$age <- search.field(metadata$sample_attribute, "age")

# Get disease
metadata$disease <- search.field(metadata$sample_attribute, "disease")



# Get population
metadata$population <- str_replace_all(search.field(metadata$sample_attribute, "population"), 
                                       regex('hapmap ', ignore_case = TRUE), "")
# Get race
metadata$race <- search.field(metadata$sample_attribute, "race")

metadata$population <- ifelse((metadata$population == 'NA') & (metadata$race == 'NA'),'NA',
                              ifelse(metadata$population == 'NA', metadata$race, metadata$population))
metadata$race <- NULL

# Get sex of individuals
metadata$sex <- search.field(metadata$sample_attribute, "sex") %>%
    str_replace(regex("^male$", ignore_case = TRUE), "M") %>%
    str_replace(regex("^female$", ignore_case = TRUE), "F") %>%
    str_replace("^U$", "NA") %>%
    str_replace(regex("^missing$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^asexual$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^mixed$", ignore_case = TRUE), "B") %>%
    str_replace(regex("^mixed sex$", ignore_case = TRUE), "B") %>%
    str_replace(regex("^mixture$", ignore_case = TRUE), "B") %>%
    str_replace('1 Male, 2 Female', "M, FF") %>%
    str_replace('2 Male, 1 Female', "MM, F") %>%
    str_replace('^3 Female$', "FFF") 

# Get gender
metadata$gender <- search.field(metadata$sample_attribute, "gender") %>%
    str_replace(regex("^male$", ignore_case = TRUE), "M") %>%
    str_replace(regex("^female$", ignore_case = TRUE), "F") %>%
    str_replace(regex("^male$", ignore_case = TRUE), "M") %>%
    str_replace("^MAL$", "M") %>%
    str_replace("^FEM$", "F") %>%
    str_replace('XY', "M") %>%
    str_replace('XX', "F") %>%
    str_replace(regex("^mixture$", ignore_case = TRUE), "B") %>%
    str_replace(regex("^mixed$", ignore_case = TRUE), "B") %>%
    str_replace("^3 Female$", "FFF")


metadata$sex <- ifelse((metadata$sex == 'NA') & (metadata$gender == 'NA'),'NA',
                       ifelse(metadata$sex == 'NA', metadata$gender, metadata$sex))

metadata$gender <- NULL

# Cheking merged columns
# table(metadata$gender)[c(1,3,4)] + table(metadata$sex)[c(2,4,7)]

# Get source type of samples
metadata$source_name <- search.field(metadata$sample_attribute, "source_name")

# Removing new line characters
metadata$read_spec <- gsub("\n", "", metadata$read_spec)

# Create manifest file
labels <- c("study_accession", "sample_accession", 
            "experiment_accession", "run_accession")

sample_labels <- as.vector(apply(metadata[,labels], 1 ,
                                 function(x){paste( x ,collapse = "_")}))


s_i <- grep("_", metadata$URL, invert = TRUE)
p_i <- grep("_", metadata$URL)

single <- cbind(metadata$URL[s_i], metadata$md5[s_i], sample_labels[s_i] %p% "-1-1")
rownames(single) <- s_i
paired <- cbind(metadata$URL[p_i], metadata$md5[p_i], metadata$URL_R[p_i], metadata$md5_R[p_i], 
                sample_labels[p_i] %p% "-1-1")
rownames(paired) <- p_i
order_list <- c(p_i,s_i)
rownames(metadata) <- 1:nrow(metadata)

paired_as_single <- subset(paired, paired[,3] == 'NA' | paired[,4] == 'NA')
print("Number of studies reported as paired but just one fastq file is given:")
print(nrow(paired_as_single))
print(paired_as_single[,c(1,5)])


order_list <- c(as.numeric(rownames(paired)),as.numeric(rownames(single)))


convert_mis_to_na <- . %>% str_replace_na() %>%
    str_replace(regex("^not_documented$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^not determined$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^not applicable$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^Unknown$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^not collected$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^none provided$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^unspecified$", ignore_case = TRUE), "NA") %>%
    str_replace(regex("^none$", ignore_case = TRUE), "NA") %>%
    str_replace("^N/A$", "NA") %>%
    str_replace("^<NA>$", "NA") %>%
    str_replace("^--$", "NA") %>% 
    str_replace("$^", "NA") %>%
    str_replace("$ ^", "NA") %>%
    ifelse(. == 'NA', NA, .)
    

metadata <- lapply(metadata, convert_mis_to_na) %>% 
    lapply(str_trim) %>% as.data.frame()

# Write table with all Illumina data
write.table(metadata[match(order_list, rownames(metadata)),], "all_illumina_sra_for_human.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(paired, "manifest_file_illumina_sra_human", sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
write.table(single, "manifest_file_illumina_sra_human", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE, append = TRUE)

# Ensure reproducibility
options(width = 120)
devtools::session_info()

# MANIFEST FILE FORMAT
# <FASTQ URL>(tab)<optional MD5>(tab)<sample label>
# <FASTQ URL 1>(tab)<optional MD5 1>(tab)<FASTQ URL 2>(tab)<optional MD5 2>(tab)<sample label>
