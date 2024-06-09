library("karyoploteR")

plot_primers <- function(
    organism_name, genome_size, primers_data, genes_data
    ) {
    
    virus <- data.frame(
        genome = organism_name,
        start = 1,
        end = genome_size
    )
    
    custom.genome <- toGRanges(virus)
    
    genes <- read.table(genes_data, header = TRUE)
    genes$chr <- organism_name
    custom.genes <- toGRanges(genes)
    
    
    primers <- read.table(primers_data, header = TRUE)
    primers$genome <- organism_name
    
    png(
        "Primers.png",
        width = 5,
        height = 6,
        units = "in",
        res = 1200,
        pointsize = 6
    )

    kp <-
        plotKaryotype(
            genome = custom.genome,
            cytobands = custom.genes,
            plot.type = 2
    )

    kpAddCytobandLabels(
        kp,
        force.all = TRUE,
        srt = 0,
        col = "black",
        cex = .5
    )

    # Genome axis
    kpAddBaseNumbers(
        kp,
        tick.dist = 1000,
        tick.len = 5,
        add.units = TRUE,
        minor.tick.dist = 200,
        minor.tick.len = 3,
        cex = 0.51,
        units = "b"
    )

    for (i in 1:max(row(primers)))
        if (i  %% 2) {
            region <- primers[i, ]
            region_UP <- joinRegions(region)
            kpRect(
                kp,
                data = region_UP,
                y0 = 0,
                y1 = 0.1,
                col = c(circlize::rand_color(
                    n = i)),
                border = 1,
                r0 = 0,
                r1 = 0.5
            )
            
            marks_1 <- toGRanges(
                data.frame(
                    genome = organism_name,
                    start = c(region$start, region$end)
                )
            )

            mcols(marks_1) <- data.frame(labels = c(
                paste0("P", i, " Start: ", region$start),
                paste0("P", i, " End: ", region$end))
            )
            
            kpPlotMarkers(
                kp,
                data = marks_1,
                label.color = "#333333",
                r1 = 0.35,
                cex = 0.5,
                label.margin = 5,
                data.panel = 1,
                text.orientation = "vertical",
                adjust.label.position = FALSE
            )
        } else {
            
            region <- primers[i, ]
            region_OP <- joinRegions(region)
            
            kpRect(
                kp,
                data = region_OP,
                y0 = 0,
                y1 = 0.1,
                col = c(circlize::rand_color(
                    n = i)),
                border = 1,
                r0 = 0,
                r1 = 0.5,
                data.panel = 2
            )
            
            marks_2 <- toGRanges(
                data.frame(
                    genome = organism_name,
                    start = c(region$start, region$end)
                )
            )
            
            mcols(marks_2) <- data.frame(labels = c(
                paste0("P", i, " Start: ", region$start),
                paste0("P", i, " End: ", region$end))
            )
        
            kpPlotMarkers(
                kp,
                data = marks_2,
                label.color = "#333333",
                r1 = 0.35,
                cex = 0.5,
                label.margin = 5,
                data.panel = 2,
                text.orientation = "vertical",
                adjust.label.position = FALSE
            )
        }
    dev.off()
}

plot_primers(
    organism_name = "HSV-1",
    genome_size = 152222,
    primers_data = "C:/Users/AlexSir/Desktop/Infz/HSV-1_Genes_Ampl.tsv",
    genes_data = "C:/Users/AlexSir/Desktop/Infz/HSV-1_Genes_Ann.tsv"
)
