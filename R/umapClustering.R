#' Prepare data for UMAP
#' @export
#' @usage prepareUmapData(longDf = NULL, taxonRank = NULL, type = "taxa", 
#'     taxDB = NULL, filterVar = "both", cutoff = 0, groupLabelsBy = "taxa")
#' @param longDf input phyloprofile file in long format
#' @param taxonRank taxonomy rank for labels (e.g. "phylum")
#' @param type type of clustering, either "taxa" (default) or "genes"
#' @param taxDB path to taxonomy database
#' @param filterVar choose variable (either "var1", "var2" or "both") to filter 
#' the data. Default: "both"
#' @param cutoff cutoff to filter data values. Default: 0
#' @param groupLabelsBy group labels by the number of "taxa" (default) or 
#' "genes"
#' @return A dataframe in wide format
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' prepareUmapData(longDf, "phylum")

prepareUmapData <- function(
        longDf = NULL, taxonRank = NULL, type = "taxa", taxDB = NULL, 
        filterVar = "both", cutoff = 0, groupLabelsBy = "taxa"
) {
    if (is.null(longDf)) stop("Input data cannot be NULL")
    if (is.null(taxonRank)) stop("Taxon rank must be specified!")
    FAS_F <- FAS_B <- geneID <- ncbiID <- abbrName <- fullName <- NULL
    var1 <- var2 <- Freq <- Label <- NULL
    # filter and subset input df
    filterFlag <- 1
    if (ncol(longDf) == 6) {
        colnames(longDf) <- c(
            "geneID", "ncbiID", "orthoID", "var1", "var2", "geneName"
        )
    } else if (ncol(longDf) == 5) {
        colnames(longDf) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
    } else if (ncol(longDf) == 4) {
        colnames(longDf) <- c("geneID", "ncbiID", "orthoID", "var1")
    } else if (ncol(longDf) == 3) {
        colnames(longDf) <- c("geneID", "ncbiID", "orthoID")
        longDf$var1 <- 1
        filterFlag <- 0
    }
    if (filterFlag == 1) {
        if (filterVar == "both") {
            longDfSub <- longDf %>% filter(var1 >= cutoff & var2 >= cutoff)
        } else if (filterVar == "var1") {
            longDfSub <- longDf %>% filter(var1 >= cutoff)
        } else if (filterVar == "var2") {
            longDfSub <- longDf %>% filter(var2 >= cutoff)
        }
    } else  longDfSub <- longDf
    # calculate mean of var1 and var2 (if var2 available)
    if (all(c("var1", "var2") %in% colnames(longDf))) 
        longDfSub <- longDfSub %>% mutate(var1 = (var1 + var2) / 2) 
    longDfSub <- longDfSub %>% select(geneID, ncbiID, var1)
    
    # get taxon names for input taxa based on a selected supertaxon rank
    taxMatrix <- getTaxonomyMatrix(taxDB)
    nameList <- getNameList(taxDB)
    superTaxonDf <- taxMatrix %>% 
        filter(abbrName %in% longDfSub$ncbiID) %>% 
        select(abbrName, one_of(taxonRank))
    colnames(superTaxonDf) <- c("ncbiID", "superID")
    superTaxonDf <- dplyr::left_join(
        superTaxonDf, nameList, by = c("superID" = "ncbiID")
    )
    # transform to wide format, add taxon names and ortholog count
    # if co-orthologs present, get max value (e.g. FAS) as the representative
    if (type == "taxa") {
        wideDf <- data.table::dcast(
            setDT(longDfSub), ncbiID ~ geneID, value.var = c("var1"), 
            fun.aggregate = max, fill = -1
        )
        wideDf <- dplyr::left_join(
            wideDf, superTaxonDf %>% select(ncbiID, Label = fullName), 
            by = "ncbiID"
        )
        # add count (how many seed genes each taxon has, excluding co-orthologs)
        seedWithTax <- longDfSub %>% select(geneID, ncbiID)
        seedWithTax <- seedWithTax[!duplicated(seedWithTax),]
        countDf <- data.frame(table(seedWithTax$ncbiID))
        colnames(countDf) <- c("ncbiID", "Freq")
        wideDf <- dplyr::left_join(
            data.frame(wideDf), countDf, by = c("ncbiID")
        )
    } else {
        wideDf <- data.table::dcast(
            setDT(longDfSub), geneID ~ ncbiID, value.var = c("var1"), 
            fun.aggregate = max, fill = -1
        )
        wideDf$Label <- wideDf$geneID
        # add count (how many taxa has orthologs for each seed gene)
        seedWithTax <- longDfSub %>% select(geneID, ncbiID)
        seedWithTax <- seedWithTax[!duplicated(seedWithTax),]
        countDf <- data.frame(table(longDfSub$geneID))
        colnames(countDf) <- c("geneID", "Freq")
        wideDf <- dplyr::left_join(
            data.frame(wideDf), countDf, by = c("geneID")
        )
    }
    # calculate genes / taxa frequency for each label
    if (groupLabelsBy == "genes") {
        countFreqDf <- data.frame(
            wideDf %>% group_by(Label) %>% summarise(n = max(Freq))
        )
    } else {
        countFreqDf <- data.frame(wideDf %>% count(Label))
    }
    wideDf <- left_join(wideDf, countFreqDf, by = "Label")
    wideDf$n <- as.numeric(wideDf$n)
    return(wideDf)
}

#' Perform UMAP clustering
#' @export
#' @usage umapClustering(data4umap = NULL, by = "taxa", type = "binary", 
#'     randomSeed = 123)
#' @param data4umap data for UMAP clustering (output from prepareUmapData)
#' @param by cluster data by "taxa" (default) or "genes"
#' @param type type of data, either "binary" (default) or "non-binary"
#' @param randomSeed random seed. Default: 123
#' @import umap
#' @return A list contain UMAP cluster objects
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareUmapData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4umap <- prepareUmapData(longDf, "phylum")
#' umapClustering(data4umap)

umapClustering <- function(
        data4umap = NULL, by = "taxa", type = "binary", randomSeed = 123
) {
    if (is.null(data4umap)) stop("Input data cannot be NULL!")
    ncbiID <- Label <- Freq <- geneID <- NULL
    if (by == "taxa") {
        subsetDt <- subset(data4umap, select = -c(ncbiID, Label, Freq))
        
    } else {
        subsetDt <- subset(data4umap, select = -c(geneID, Label, Freq))
    }
    if (type == "binary") {
        subsetDt[subsetDt >= 0] <- 1
        subsetDt[subsetDt < 0] <- 0
    }
    df.umap <- umap::umap(
        subsetDt, random_state = randomSeed, preserve.seed = TRUE
    )
    return(df.umap)
}

#' Reduce the number of labels for UMAP plot based on the gene/taxon frequency
#' @export
#' @usage groupLabelUmapData(data4umap = NULL, freqCutoff = 0)
#' @param data4umap data for UMAP clustering (output from prepareUmapData)
#' @param freqCutoff gene/taxon frequency cutoff
#' @return A dataframe similar to input data4umap, but with modified Label 
#' column, where less frequent labels are grouped together as "Other"
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareUmapData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4umap <- prepareUmapData(longDf, "phylum")
#' groupLabelUmapData(data4umap, freqCutoff = 3)

groupLabelUmapData <- function( data4umap = NULL, freqCutoff = 0) {
    if (is.null(data4umap)) stop("Input data cannot be NULL!")
    if (length(data4umap) == 0) stop("Input data cannot be EMPTY!")
    
    # minFreq <- tail(sort(unique(data4umap$n)), labelNr)[1]
    # keepLabel <- unique(data4umap$Label[data4umap$n >= minFreq])
    # data4umap$Label[!(data4umap$Label %in% keepLabel)] <- "[Other]"
    data4umap$Label[data4umap$n < freqCutoff] <- "[Other]"
    return(data4umap)
}

#' Create UMAP cluster plot
#' @export
#' @usage createUmapPlotData(umapData = NULL, data4umap = NULL, freqCutoff = 0, 
#'     excludeTaxa = "None")
#' @param umapData data contains UMAP cluster (output from umapClustering())
#' @param data4umap data for UMAP clustering (output from prepareUmapData())
#' @param freqCutoff gene/taxon frequency cutoff
#' @param excludeTaxa hide taxa from plot. Default: "None"
#' @importFrom utils tail
#' @return A plot as ggplot object
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareUmapData}}, \code{\link{umapClustering}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4umap <- prepareUmapData(longDf, "phylum")
#' umapData <- umapClustering(data4umap)
#' createUmapPlotData(umapData, data4umap)

createUmapPlotData <- function(
    umapData = NULL, data4umap = NULL, freqCutoff = 0, excludeTaxa = "None"
) {
    if (is.null(umapData) | is.null(data4umap)) 
        stop("Input data cannot be NULL!")
    if (length(umapData) == 0 | length(data4umap) == 0) 
        stop("Input data cannot be EMPTY!")
    
    Label <- Freq <- NULL
    data4umap$X <- umapData$layout[,1]
    data4umap$Y <- umapData$layout[,2]
    # join less freq items into "other"
    data4umap <- groupLabelUmapData(data4umap, freqCutoff)
    # exclude taxa
    if (length(excludeTaxa) > 0) {
        data4umap$X[data4umap$Label %in% excludeTaxa] <- NA
        data4umap$Y[data4umap$Label %in% excludeTaxa] <- NA
    }
    return(data4umap)
}


#' Create UMAP cluster plot
#' @export
#' @usage plotUmap(plotDf = NULL, legendPos = "right", colorPalette = "Set2", 
#'     transparent = 0, textSize = 12)
#' @param plotDf data for UMAP plot 
#' @param legendPos position of legend. Default: "right"
#' @param colorPalette color palette. Default: "Set2"
#' @param transparent transparent level (from 0 to 1). Default: 0
#' @param textSize size of axis and legend text. Default: 12
#' @return A plot as ggplot object
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareUmapData}}, \code{\link{umapClustering}},
#' \code{\link{createUmapPlotData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' umapData <- prepareUmapData(longDf, "phylum")
#' data.umap <- umapClustering(umapData)
#' plotDf <- createUmapPlotData(data.umap, umapData)
#' plotUmap(plotDf)

plotUmap <- function(
        plotDf = NULL, legendPos = "right", colorPalette = "Set2", 
        transparent = 0, textSize = 12
) {
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    X <- Y <- Label <- Freq <- NULL
    # generate plot
    plot <- ggplot(plotDf, aes(x = X, y = Y, color = Label)) +
        geom_point(aes(size = Freq), alpha = 1 - transparent) +
        geom_rug(alpha = 1) +
        theme_minimal() +
        labs(x = "", y = "") +
        guides(
            color = guide_legend(override.aes = list(alpha = 1)),
            size = guide_legend(title = "Number of genes")
        ) +
        theme(
            legend.position = legendPos,
            legend.text = element_text(size = textSize),
            legend.title = element_text(size = textSize),
            axis.text = element_text(size = textSize), 
            axis.title = element_text(size = textSize),
        )
    if (checkColorPallete(levels(as.factor(plotDf$Label)), colorPalette))
        plot <- plot + scale_color_brewer(palette = colorPalette)
    return(plot)
}