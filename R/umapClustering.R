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
    if (type == "genes") groupLabelsBy <- "genes"
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
        countDf <- data.frame(table(seedWithTax$geneID))
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
    wideDf$Label <- as.character(wideDf$Label)
    return(wideDf)
}

#' Perform UMAP clustering 2D
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
    if ("geneID" %in% colnames(data4umap)) by <- "genes"
    if (by == "taxa") {
        subsetDt <- subset(data4umap, select = -c(ncbiID, Label, Freq, n))
    } else {
        subsetDt <- subset(data4umap, select = -c(geneID, Label, Freq, n))
    }
    if (type == "binary") {
        subsetDt[subsetDt >= 0] <- 1
        subsetDt[subsetDt < 0] <- 0
    }
    checkDt4Umap <- tryCatch(
        {
            suppressWarnings(umap::umap(
                subsetDt, random_state = randomSeed, preserve.seed = TRUE
            ))
        },
        error = function(cond) {
            message(conditionMessage(cond))
            NA
        },
        warning = function(cond) {
            message(conditionMessage(cond))
            NA
        }
    )
    if (!(is.na(checkDt4Umap[1]))) {
        df.umap <- umap::umap(
            subsetDt, random_state = randomSeed, preserve.seed = TRUE
        )
    } else {
        warning("PROBLEM: Too few samples for UMAP!!!")
        df.umap <- umap::umap(
            subsetDt, random_state = randomSeed, preserve.seed = TRUE,
            n_neighbors = max(1, nrow(subsetDt) - 1)
        )
    }
    return(df.umap)
}

#' Perform UMAP clustering 3D
#' @export
#' @usage umapClustering3D(data4umap = NULL, by = "taxa", type = "binary", 
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
#' umapClustering3D(data4umap)

umapClustering3D <- function(
        data4umap = NULL, by = "taxa", type = "binary", randomSeed = 123
) {
    if (is.null(data4umap)) stop("Input data cannot be NULL!")
    ncbiID <- Label <- Freq <- geneID <- NULL
    if ("geneID" %in% colnames(data4umap)) by <- "genes"
    if (by == "taxa") {
        subsetDt <- subset(data4umap, select = -c(ncbiID, Label, Freq, n))
    } else {
        subsetDt <- subset(data4umap, select = -c(geneID, Label, Freq, n))
    }
    if (type == "binary") {
        subsetDt[subsetDt >= 0] <- 1
        subsetDt[subsetDt < 0] <- 0
    }
    checkDt4Umap <- tryCatch(
        {
            suppressWarnings(umap::umap(
                subsetDt, random_state = randomSeed, preserve.seed = TRUE,
                n_components = 3
            ))
        },
        error = function(cond) {
            message(conditionMessage(cond))
            NA
        },
        warning = function(cond) {
            message(conditionMessage(cond))
            NA
        }
    )
    if (!(is.na(checkDt4Umap[1]))) {
        df.umap <- umap::umap(
            subsetDt, random_state = randomSeed, preserve.seed = TRUE,
            n_components = 3
        )
    } else {
        warning("PROBLEM: Too few samples for UMAP!!!")
        df.umap <- umap::umap(
            subsetDt, random_state = randomSeed, preserve.seed = TRUE,
            n_neighbors = max(1, nrow(subsetDt) - 1), n_components = 3
        )
    }
    return(df.umap)
}

#' Reduce the number of labels for UMAP plot based on the gene/taxon frequency
#' @export
#' @usage groupLabelUmapData(data4umap = NULL, freqCutoff = c(0,200))
#' @param data4umap data for UMAP clustering (output from prepareUmapData)
#' @param freqCutoff gene/taxon frequency cutoff range. Any labels that are 
#' outside of this range will be assigned as [Other]
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
#' groupLabelUmapData(data4umap, freqCutoff = c(3,5))

groupLabelUmapData <- function(data4umap = NULL, freqCutoff = c(0,200)) {
    if (is.null(data4umap)) stop("Input data cannot be NULL!")
    if (length(data4umap) == 0) stop("Input data cannot be EMPTY!")
    
    data4umap$Label[
        data4umap$n < freqCutoff[1] | data4umap$n > freqCutoff[2]
    ] <- "[Other]"
    return(data4umap)
}

#' Create UMAP cluster plot
#' @export
#' @usage createUmapPlotData(umapData = NULL, data4umap = NULL, 
#'     freqCutoff = c(0,200), excludeTaxa = "None", currentNCBIinfo = NULL)
#' @param umapData data contains UMAP cluster (output from umapClustering())
#' @param data4umap data for UMAP clustering (output from prepareUmapData())
#' @param freqCutoff gene/taxon frequency cutoff range. Any labels that are 
#' outside of this range will be assigned as [Other]
#' @param excludeTaxa hide taxa from plot. Default: "None"
#' @param currentNCBIinfo table/dataframe of the pre-processed NCBI taxonomy
#' data (/PhyloProfile/data/preProcessedTaxonomy.txt)
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
    umapData = NULL, data4umap = NULL, freqCutoff = c(0, 200), 
    excludeTaxa = "None", currentNCBIinfo = NULL
) {
    if (is.null(umapData) | is.null(data4umap)) 
        stop("Input data cannot be NULL!")
    if (length(umapData) == 0 | length(data4umap) == 0) 
        stop("Input data cannot be EMPTY!")
    Label <- Freq <- NULL
    data4umap$X <- umapData$layout[,1]
    data4umap$Y <- umapData$layout[,2]
    if (ncol(umapData$layout) == 3) data4umap$Z <- umapData$layout[,3]
    # join less freq items into "other"
    data4umap <- groupLabelUmapData(data4umap, freqCutoff)
    # exclude taxa
    if (length(excludeTaxa) > 0) {
        data4umap$X[data4umap$Label %in% excludeTaxa] <- NA
        data4umap$Y[data4umap$Label %in% excludeTaxa] <- NA
    }
    # convert tax IDs into names
    if ("ncbiID" %in% colnames(data4umap)) {
        if (is.null(currentNCBIinfo)) {
            dataPath <- system.file(
                "PhyloProfile", "data/",
                package = "PhyloProfile", mustWork = TRUE
            )
            nameFullFile <- paste0(dataPath, "/preProcessedTaxonomy.txt")
            if (file.exists(nameFullFile)) {
                currentNCBIinfo <- as.data.frame(
                    data.table::fread(nameFullFile)
                )
            }
        }
        if (!is.null(currentNCBIinfo)) {
            id2nameDf <- PhyloProfile::id2name(
                gsub("ncbi","",unique(data4umap$ncbiID)), currentNCBIinfo
            )
            id2nameDf$ncbiID <- paste0("ncbi", id2nameDf$ncbiID)
            data4umap <- merge(
                data4umap, id2nameDf, by = "ncbiID", all.x = TRUE
            )
        } else {
            data4umap$fullName <- "NA"
        }
    }
    return(data4umap)
}


#' Create UMAP cluster plot
#' @export
#' @usage plotUmap(plotDf = NULL, legendPos = "bottom", colorPalette = "Set2", 
#'     transparent = 0, textSize = 12, font = "Arial", highlightTaxa = NULL,
#'     dotZoom = 0)
#' @param plotDf data for UMAP plot 
#' @param legendPos position of legend. Default: "right"
#' @param colorPalette color palette. Default: "Set2"
#' @param transparent transparent level (from 0 to 1). Default: 0
#' @param textSize size of axis and legend text. Default: 12
#' @param font font of text. Default = Arial"
#' @param highlightTaxa list of taxa to be highlighted
#' @param dotZoom dot size zooming factor. Default: 0
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
#' plotUmap(plotDf, font = "sans")

plotUmap <- function(
    plotDf = NULL, legendPos = "bottom", colorPalette = "Set2", 
    transparent = 0, textSize = 12, font = "Arial", highlightTaxa = NULL,
    dotZoom = 0
) {
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    X <- Y <- Label <- Freq <- NULL
    # add colors
    plotDf <- addUmapTaxaColors(plotDf, colorPalette, highlightTaxa)
    if ("color" %in% colnames(plotDf)) {
        colorScheme <- structure(
            plotDf$color, .Names = plotDf$Label
        )
    } else {
        colorScheme <- structure(
            head(
                suppressWarnings(
                    RColorBrewer::brewer.pal(
                        nlevels(as.factor(plotDf$Label)), colorPalette
                    )
                ), levels(as.factor(plotDf$Label))
            ), .Names = levels(as.factor(plotDf$Label))
        )
    }
    # adapt plot height based on number of labels

    # generate plot
    plot <- ggplot(plotDf, aes(x = X, y = Y, color = Label)) +
        geom_point(aes(size = Freq), alpha = 1 - transparent) +
        geom_rug(alpha = 1) +
        theme_minimal() +
        labs(x = "", y = "")+
        scale_radius(range = c(1 + dotZoom, 6 + dotZoom))
    # change legend title
    if ("ncbiID" %in% colnames(plotDf)) {
        plot <- plot + guides(
            color = guide_legend(override.aes = list(alpha = 1), ncols = 5),
            size = guide_legend(title = "Number of genes", ncols = 2)
        )
    } else
        plot <- plot + guides(
            color = guide_legend(override.aes = list(alpha = 1), ncols = 5),
            size = guide_legend(title = "Number of taxa", ncols = 1)
        )
    plot <- plot + theme(
            legend.position = legendPos,
            legend.text = element_text(size = textSize),
            legend.title = element_text(size = textSize),
            axis.text = element_text(size = textSize), 
            axis.title = element_text(size = textSize),
        )
    plot <- plot + scale_color_manual(values = colorScheme)
    plot <- plot + theme(text = element_text(family = font))
    return(plot)
}


#' Create UMAP cluster 3D plot
#' @export
#' @usage plotUmap3D(plotDf = NULL, legendPos = "bottom", 
#'     colorPalette = "Set2", transparent = 0,highlightTaxa = NULL, 
#'     dotZoom = 0)
#' @param plotDf data for UMAP plot 
#' @param legendPos position of legend. Default: "right"
#' @param colorPalette color palette. Default: "Set2"
#' @param transparent transparent level (from 0 to 1). Default: 0
#' @param highlightTaxa list of taxa to be highlighted
#' @param dotZoom dot size zooming factor. Default: 0
#' @return A plot as ggplot object
#' @rawNamespace import(plotly, except = last_plot)
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareUmapData}}, \code{\link{umapClustering}},
#' \code{\link{createUmapPlotData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' umapData <- prepareUmapData(longDf, "phylum")
#' data.umap <- umapClustering3D(umapData)
#' plotDf <- createUmapPlotData(data.umap, umapData)
#' plotUmap3D(plotDf)

plotUmap3D <- function(
    plotDf = NULL, legendPos = "bottom", colorPalette = "Set2", 
    transparent = 0, highlightTaxa = NULL, dotZoom = 0
) {
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    X <- Y <- Z <- Label <- Freq <- color <- NULL
    if (!("Z" %in% colnames(plotDf))) stop("PlotDf seems to be 2D plot!")
    # add colors
    plotDf <- addUmapTaxaColors(plotDf, colorPalette, highlightTaxa)
    colorDf <- unique(plotDf %>% select(Label, color))
    # calculate min and max dot size
    minSize <- dotZoom
    maxSize <- (max(plotDf$Freq)*dotZoom/min(plotDf$Freq))
    # generate 3D scatter plot
    if ("ncbiID" %in% colnames(plotDf)) {
        plot <- plot_ly(
            plotDf, x = ~X, y = ~Y, z = ~Z,
            hovertemplate = ~paste0(
                "Taxon ID: ", gsub("ncbi","",ncbiID), "<br>Name: ", fullName,
                "<br>Label: ", Label, "<br>Freq: ", Freq
            )
        )
    } else {
        plot <- plot_ly(
            plotDf, x = ~X, y = ~Y, z = ~Z,
            hovertemplate = ~paste0(
                "Gene ID: ", geneID, "<br>Freq: ", Freq
            )
        )
    }
    
    if (dotZoom == 0) {
        plot <- plot %>% add_markers(
            color = ~Label, colors = colorDf$color, size = ~Freq, 
            alpha = 1 - transparent, span = I(0)
        )
    } else {
        plot <- plot %>% add_markers(
            color = ~Label, colors = colorDf$color, size = ~Freq, 
            sizes = c(minSize, maxSize), alpha = 1 - transparent, span = I(0)
        )
    }
    if (legendPos == "none") plot <- plot %>% hide_legend()
    return(plot)
} 


#' Add colors for taxa in UMAP plot
#' @usage addUmapTaxaColors(plotDf = NULL, colorPalette = "Set2", 
#'     highlightTaxa = NULL)
#' @param plotDf data for UMAP plot 
#' @param colorPalette color palette. Default: "Set2"
#' @param highlightTaxa list of taxa to be highlighted
#' @return A dataframe for UMAP plot with an additional column for the assigned
#' color to each taxon
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
#' PhyloProfile:::addUmapTaxaColors(plotDf, colorPalette = "Set2")

addUmapTaxaColors <- function(
    plotDf = NULL, colorPalette = "Set2", highlightTaxa = NULL
) {
    if (is.null(plotDf)) stop("plotDf cannot be null!")
    allTaxa <- levels(as.factor(plotDf$Label))
    allColors <- getQualColForVector(allTaxa)
    if (checkColorPalette(allTaxa, colorPalette)) {
        allColors <-
            suppressWarnings(head(
                RColorBrewer::brewer.pal(length(allTaxa), colorPalette),
                length(allTaxa)
            ))
    }
    colorScheme <- data.frame(color = allColors, Label = allTaxa)
    if (!(is.null(highlightTaxa)))
        colorScheme$color[!(colorScheme$Label %in% highlightTaxa)] <- "#d4d6d9"
    plotDf <- merge(plotDf, colorScheme, by = "Label", all.x = TRUE)
    return(plotDf)
}
