library("karyoploteR")


plot_coverage <- function(
    organism_name, genome_size, coverage_data, genes_data
    ) {

    virus <- data.frame(
        genome = organism_name,
        start = 1L,
        end = genome_size
    )
    
    coverage_data <- read.table(coverage_data, header = TRUE)
    coverage_data$name <- organism_name
    
    genes <- read.table(genes_data, header = TRUE)
    genes$chr <- organism_name

    custom.genome = regioneR::toGRanges(virus)
    custom.genes = regioneR::toGRanges(genes)

    length <- max(col(coverage_data))
    col_namews <- colnames(coverage_data)

    png(
        "Coverage.png",
        width = 5,
        height = 6,
        units = "in",
        res = 1200,
        pointsize = 6
    )

    kpp <- getDefaultPlotParams(plot.type = 1)
    kpp$data1height <- 300
    kpp$ideogramheight <- 30
    
    kp <-
        plotKaryotype(
            genome = custom.genome,
            cytobands = custom.genes,
            plot.type = 2,
            plot.params = kpp
        )
    
    kpAddCytobandLabels(
        kp,
        force.all = FALSE,
        srt = 0,
        col = "black",
        cex = .5
    )
    
    end <- max(genes$end)
    
    # Genome axis
    
    kpAddBaseNumbers(
        kp,
        tick.dist = ifelse(end <= 30000, 1000, 10000),
        tick.len = 5,
        add.units = TRUE,
        minor.tick.dist = ifelse(end <= 30000, 200, 2000),
        minor.tick.len = 3,
        cex = 0.51,
        units = "auto"
    )
    
    for (i in 3:length) {
        if ((length - 3) <= 4) {
            g <- i - 2
            j0 = -0.18 + g * 0.23
            j1 = -0.16 + g * 0.23
            
            region_df <- data.frame(
                name = coverage_data$name,
                start = coverage_data[, i],
                end = coverage_data[, i]
            )
            
            region_df_UP <- toGRanges(
                region_df,
                genome = custom.genome
            )
            
            kpDataBackground(
                kp,
                r0 = j0,
                r1 = j1 + 0.05,
                color = circlize::rand_color(
                    n = 1,
                    hue = "monochrome",
                    luminosity = "light"
                )
            )
            
            kpLines(
                kp,
                data = region_df_UP,
                chr = organism_name,
                x = 1:length(region_df$start),
                y = log10(region_df$start),
                r0 = j0,
                r1 = j1,
                col = "black"
            )
            
            kpAxis(
                kp,
                ymax = log10(signif(mean(coverage_data[, i]), digits = 1)),
                side = 2,
                r0 = j0,
                r1 = j1 + 0.05,
                data.panel = 1,
                chromosomes = organism_name,
                cex = 0.5,
                numticks = 5
            )
        
            kpAddLabels(
                kp,
                labels = col_namews[i],
                srt = 90,
                r0 = j0,
                r1 = j1 + 0.15,
                data.panel = 1
            )
        
            kpAddLabels(
                kp,
                labels = "Coverage (log10)",
                r0 = length / 12,
                r1 = length / 12,
                cex = 1,
                srt = 90,
                label.margin = 0.035
            )
        
        } else {
            g <- i - 2
            j0 = -0.1 + g * 0.1
            j1 = -0.08 + g * 0.1
            
            region_df <- data.frame(
                name = coverage_data$name,
                start = coverage_data[, i],
                end = coverage_data[, i]
            )
        
            region_df_UP <- toGRanges(
                region_df,
                genome = custom.genome
            )
        
            kpDataBackground(
                kp,
                r0 = j0,
                r1 = j1 + 0.05,
                color = circlize::rand_color(
                    n = 1,
                    hue = "monochrome",
                    luminosity = "light"
                )
            )
            
            kpLines(
                kp,
                data = region_df_UP,
                chr = organism_name,
                x = 1:length(region_df$start),
                y = log10(region_df$start),
                r0 = j0,
                r1 = j1,
                col = "black"
            )
            
            kpAxis(
                kp,
                ymax = signif(mean(coverage_data[, i]), digits = 1),
                side = 2,
                r0 = j0,
                r1 = j1 + 0.05,
                data.panel = 1,
                chromosomes = organism_name,
                cex = 0.5
            )
            
            kpAddLabels(
                kp,
                labels = col_namews[i],
                srt = 90,
                r0 = j0,
                r1 = j1 + 0.12,
                cex = 0.75,
                data.panel = 1
            )
            
            kpAddLabels(
                kp,
                labels = "Coverage (log10)",
                r0 = length / 20,
                r1 = length / 20,
                cex = 1,
                srt = 90,
                label.margin = 0.035
            )
        }
    }
    
    dev.off()
}


plot_coverage(
    organism_name = "hPIV-3",
    genome_size = 7600,
    coverage_data = "F:/Telegram Desktop/coverage.txt",
    genes_data = "F:/Telegram Desktop/genes.txt"
)
