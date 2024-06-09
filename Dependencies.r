if (!require("BiocManager", quietly = TRUE) == FALSE) {
    message("Installing the 'BiocManager' package")
    install.packages("BiocManager")
}

if (!require("karyoploteR")) {
    message("Installing the 'karyoploteR' package")
    BiocManager::install("karyoploteR")
}

if (!require("circlize")) {
    message("Installing the 'circlize' package")
    install.packages("circlize")
}