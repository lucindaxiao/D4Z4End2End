# D4Z4End2End methylation analysis (using NanoMethViz)
# usage: Rscript D4Z4_methylation.R -s [sample] -g [group] -a [annotation file] -b [modBAM file] -o [output directory]

library(optparse)
library(NanoMethViz)
library(dplyr)

options("NanoMethViz.site_filter" = 1)

option_list <- list(
    make_option(c("-s", "--sample"), type = "character", default = "sample", help = "Sample name [default: %default]", metavar = "character"),
    make_option(c("-g", "--group"), type = "character", default = "group", help = "Group name [default: %default]", metavar = "character"),
    make_option(c("-a", "--annofile"), type = "character", default = "D4Z4_haplotyping_ref_seqs/reads_features.txt", help = "", metavar = "character"),
    make_option(c("-b", "--bamfile"), type = "character", default = "D4Z4_alleles/D4Z4_alleles.bam", help = "BAM file path [default: %default]", metavar = "character"),
    make_option(c("-o", "--outdir"), type = "character", default = "D4Z4_methylation", help = "Output directory [default: %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

sample <- opt$sample
group <- opt$group
annofile <- paste0("output/", opt$annofile)
bamfile <- paste0("output/", opt$bamfile)
outdir_plots <- paste0("output/", opt$outdir, "/5mC_plots")
outdir_stats <- paste0("output/", opt$outdir, "/5mC_stats")

if (!dir.exists(outdir_plots)) {
    dir.create(outdir_plots, recursive = TRUE)
}
if (!dir.exists(outdir_stats)) {
    dir.create(outdir_stats, recursive = TRUE)
}

## Get D4Z4 unit and DUX4 exon annotations from the reads_features.txt file

anno_data <- read.table(annofile, header = TRUE)

anno_df <- data.frame(
    chr = anno_data$read_name,
    gene_id = anno_data$element,
    transcript_id = anno_data$element,
    symbol = anno_data$element,
    strand = anno_data$strand,
    start = anno_data$start,
    end = anno_data$end
)

anno_D4Z4 <- anno_df %>%
    dplyr::filter(grepl("D4Z4_unit", symbol)) %>%
    dplyr::mutate(symbol = "D4Z4")

anno_DUX4 <- anno_df %>%
    dplyr::filter(grepl("DUX4", symbol))

## Plot D4Z4 methylation and calculate methylation statistics for each allele

combined_methy_full_array <- data.frame()

unique_haps <- unique(anno_D4Z4$chr)

for (haplo in unique_haps) {
  
    filtered_anno_D4Z4 <- anno_D4Z4 %>%
        dplyr::filter(chr == haplo)
    
    filtered_anno_DUX4 <- anno_DUX4 %>%
        dplyr::filter(chr == haplo)
    
    start <- min(filtered_anno_D4Z4$start) - 1000
    end <- max(filtered_anno_D4Z4$end) + 1000
    
    mbr <- ModBamResult(
        methy = ModBamFiles(
            samples = sample,
            paths = bamfile
        ),
        samples = data.frame(
            sample = sample,
            group = group
        ),
        exons = filtered_anno_D4Z4
    )
    
    # generate methylation plots

    plot <- NanoMethViz:::plot_region(mbr, haplo, start, end, anno_regions = filtered_anno_DUX4, heatmap=TRUE, line_size=0.5, smoothing_window=2000)
    ggsave(filename = paste0(sample, "_", haplo, ".png"), plot = plot, path = outdir_plots, width = 8.06, height = 5.89, units = "in")
    
    # calculate RU-level methylation statistics (mean methylation probability and mean methylation rate per D4Z4 unit)
    
    methy_per_RU <- region_methy_stats(mbr, filtered_anno_D4Z4, threshold = 0.5)

    if (any(methy_per_RU$strand == '-')) {
        methy_per_RU_minus <- methy_per_RU %>%
            dplyr::filter(strand == '-') %>%
            dplyr::arrange(start)
        methy_per_RU_minus <- methy_per_RU_minus %>%
            dplyr::mutate(RU = dplyr::row_number() - nrow(methy_per_RU_minus) - 1)
        methy_per_RU_plus <- methy_per_RU %>%
            dplyr::filter(strand == '+') %>%
            dplyr::arrange(start) %>%
            dplyr::slice(-1, -n()) %>%
            dplyr::mutate(RU = dplyr::row_number())
        methy_per_RU <- rbind(methy_per_RU_minus, methy_per_RU_plus)
    } else {
        methy_per_RU <- methy_per_RU %>%
            dplyr::filter(strand == '+') %>%
            dplyr::arrange(start) %>%
            dplyr::slice(-1, -n()) %>%
            dplyr::mutate(RU = dplyr::row_number())
    }
    
    methy_per_RU <- methy_per_RU %>%
        dplyr::select(chr, gene_id, strand, start, end, mean_methy_prob, prevalence, RU) %>%
        dplyr::rename(haplotype = chr, feature = gene_id, mean_methy_rate = prevalence)

    write.table(methy_per_RU, file = paste0(outdir_stats, "/", sample, "_", haplo, "_5mCperRU.tsv"), sep = "\t", row.names = FALSE)

    # calculate array-level methylation statistics (mean methylation probability and mean methylation rate over the full D4Z4 array)

    filtered_anno_D4Z4_plus <- filtered_anno_D4Z4 %>%
        dplyr::filter(strand == '+')
    
    D4Z4_array <- data.frame(
        chr = haplo,
        gene_id = "D4Z4_array",
        transcript_id = "D4Z4_array",
        symbol = "D4Z4_array",
        strand = "+",
        start = min(filtered_anno_D4Z4_plus$start),
        end = max(filtered_anno_D4Z4_plus$end)
    )
    methy_full_array <- region_methy_stats(mbr, D4Z4_array, threshold = 0.5)
    methy_full_array <- methy_full_array %>%
        dplyr::select(chr, gene_id, strand, start, end, mean_methy_prob, prevalence) %>%
        dplyr::rename(haplotype = chr, feature = gene_id, mean_methy_rate = prevalence)
    combined_methy_full_array <- rbind(combined_methy_full_array, methy_full_array)
}

write.table(combined_methy_full_array, file = paste0(outdir_stats, "/", sample, "_5mC_fullD4Z4.tsv"), sep = "\t", row.names = FALSE)