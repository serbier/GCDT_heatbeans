rm(list = ls())
pks <- c('argparse','multtest','gplots','genetics','ape','EMMREML','compiler','GAPIT3','scatterplot3d')

pks <- suppressPackageStartupMessages(sapply(pks, require, character.only = TRUE))

if (any(!pks)) stop('The package(s) ', names(pks)[!pks], ' is/are not available for load.')

pr <- ArgumentParser()

pr$add_argument("-geno", type = "character", metavar = 'file',
                    help = "SNP information in hapmap format")

pr$add_argument("-pheno", type = "character", metavar = "file",
                    help = "Phenotypic information indexed with genotype")

pr$add_argument("-struc", type = "character", metavar = "title",
                    help = "Structure file for population structure adjust", default = NULL)
pr$add_argument("-pca", type = "integer", metavar = "title",
                    help = "Number of PCA components for population structure adjust", default = NULL)
pr$add_argument("-modelSelection", type = "logical", metavar = "title",
                    help = "Selection of best number of components for PCA", default = FALSE)
pr$add_argument("-out", type = "character", metavar = "file",
                    help = "out path")

ar = pr$parse_args()


if(!is.null(ar$struc)){
    myCV <- read.table(ar$struc, head = F)
} else {
    myCV <- NULL
}


myY = read.csv(ar$pheno, header = T)

initial = read.table(ar$geno,nrows = 100, comment.char='')
classes = sapply(initial, class)
myG     = read.table(ar$geno, colClasses=classes, header=F, comment.char='')

setwd(ar$out)

mygapit = GAPIT(
        Y=myY,
        G=myG,
        model=c("MLM","Blink"),
        CV=myCV,
        PCA.total=ar$pca,
        Model.selection = ar$modelSelection)