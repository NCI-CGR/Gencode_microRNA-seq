## library is sense stranded
library(stringr)
files <- list.files(
        path = "star_align",
        recursive = TRUE,
        pattern = "ReadsPerGene.out.tab",
        full.names = TRUE)
sample_names <- list()
for (file in files) {
        basename <- basename(file)
        sample <- str_replace(basename, "ReadsPerGene.out.tab", "")
        sample_names <- append(sample_names, sample)
        if (!exists("rc")) {
                rc <- read.table(
                        file, header = FALSE, row.names = 1,
                        colClasses = c(NA,"NULL",NA,"NULL"),
                        col.names=c("gene_id", "NULL", file, "NULL"))
        } else if (exists("rc")) {
                temp_rc=read.table(
                        file,header = FALSE, row.names = 1,
                        colClasses = c(NA, "NULL", NA, "NULL"),
                        col.names=c("gene_id", "NULL", file, "NULL"))
                rc <- cbind(rc, temp_rc)
                rm(temp_rc)
        }
}

colnames(rc) <- sample_names
dir.create(file.path("reads_count"), showWarnings = FALSE)
write.csv(rc, "reads_count/reads_count.csv")
