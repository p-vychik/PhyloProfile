#' Get all taxa that share a common ancestor
#' @description Identify the common ancestor for a selected taxa and return a
#' list of all taxa that have that common ancestor from an large input taxa set.
#' @export
#' @param inputTaxa ID list of all input taxa (e.g. "ncbi12345")
#' @param inGroup ID list of selected taxa used for identify the common ancestor
#' (e.g.: "ncbi55555")
#' @param taxDB Path to the taxonomy DB files
#' @return A list containing the taxonomy rank and name of the common ancestor,
#' together with a dataframe storing the full taxonomy info of all taxa that
#' share that corresponding common ancestor.
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @examples
#' inputTaxa <- c("ncbi34740", "ncbi9606", "ncbi374847", "ncbi123851",
#'     "ncbi5664", "ncbi189518", "ncbi418459", "ncbi10116", "ncbi284812",
#'     "ncbi35128", "ncbi7070")
#' inGroup <-  c("ncbi9606", "ncbi10116")
#' getCommonAncestor(inputTaxa, inGroup)

getCommonAncestor <- function(inputTaxa = NULL, inGroup = NULL, taxDB = NULL) {
    if (is.null(inputTaxa) | is.null(inGroup))
        stop("Input taxa and in-group ID list cannot be NULL!")
    # get list of pre-calculated taxonomy info
    taxMatrix <- getTaxonomyMatrix(taxDB, TRUE, inputTaxa)
    # get subset taxonomy info for selected in-group taxa
    selectedTaxMatrix <- taxMatrix[
        taxMatrix$abbrName %in% inGroup,
        which(!duplicated(t(taxMatrix)))
        ]
    # identify common ancestor
    V1 <- vapply(
        selectedTaxMatrix[,c(seq(4, ncol(selectedTaxMatrix)))], max,
        FUN.VALUE = numeric(1)
    )
    V2 <- vapply(
        selectedTaxMatrix[,c(seq(4, ncol(selectedTaxMatrix)))],
        function (x) as.integer(sum(x)/length(x)),
        FUN.VALUE = numeric(1)
    )
    checkDf <- t(rbind(V1, V2))
    commonRank <- rownames(checkDf[V1 == V2,])[1]
    commonID <- checkDf[V1 == V2,][1]
    commonTaxa <- taxMatrix[taxMatrix[, commonRank] == commonID, ]
    return(list(commonRank, commonID, commonTaxa))
}

#' Compare the score distributions between 2 taxon groups
#' @description Given the phylogenetic profiles that contains up to 2 additional
#' variables besides the presence/absence information of the orthologous
#' proteins. This function will compare the distribution of those variables
#' between 2 different taxon groups (e.g. parasitic species vs non-parasitic
#' species), which are defined as in-group and out-group. In-group is identified
#' by the user. Out-group contains all taxa in the input phylogenetic profiles
#' that are not part of the in-group.
#' @usage compareTaxonGroups(data, inGroup, useCommonAncestor, variable,
#'     significanceLevel, taxDB)
#' @export
#' @param data input phylogenetic profile in long format (see ?mainLongRaw and
#' ?createLongMatrix)
#' @param inGroup ID list of in-group taxa (e.g. "ncbi1234")
#' @param useCommonAncestor TRUE/FALSE if using all taxa that share the same
#' common ancestor with the pre-selected in-group as the in-group taxa.
#' Default = TRUE.
#' @param variable name of the variable that need to be compared
#' @param significanceLevel significant cutoff for the statistic test (between
#' 0 and 1). Default = 0.05.
#' @param taxDB Path to the taxonomy DB files
#' @return list of genes that have a significant difference in the variable
#' distributions between the in-group and out-group taxa and their corresponding
#' p-values.
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi9606", "ncbi10116")
#' variable <- colnames(data)[4]
#' compareTaxonGroups(data, inGroup, TRUE, variable, 0.05)

compareTaxonGroups <- function(
    data = NULL,
    inGroup = NULL,
    useCommonAncestor = TRUE,
    variable = NULL,
    significanceLevel = 0.05,
    taxDB = NULL
) {
    if (is.null(data) | is.null(inGroup) | is.null(variable))
        stop("Input profiles, in-group IDs and variable name cannot be NULL!")
    if (!(variable %in% colnames(data))) stop("Invalid variable")
    # add other taxa that share a common ancestor with the given in-group
    commonTaxa <- getCommonAncestor(
        levels(as.factor(data$ncbiID)), inGroup, taxDB
    )
    if (useCommonAncestor == TRUE)
        inGroup <- as.character(commonTaxa[[3]]$abbrName)

    # perform distribution comparison test
    data$geneID <- as.character(data$geneID)
    pvalues <- vapply(
        unique(data$geneID),
        function (x) {
            varIn <- data[data$geneID == x & data$ncbiID %in% inGroup, variable]
            varOut <- data[
                data$geneID == x & !(data$ncbiID %in% inGroup), variable
                ]
            pvalue <- try(distributionTest(varIn, varOut, significanceLevel))
            if (is.numeric(pvalue))
                return(distributionTest(varIn, varOut, significanceLevel))
            else return(999)
        }, FUN.VALUE = numeric(1)
    )
    return(sort(unlist(pvalues)))
}

#' Compare the distribution of 2 numeric vectors
#' @description This function tests the difference between the distributions of
#' two input numeric samples using the statistical tess. First the
#' Kolmogorov-Smirnov is used to check if 2 samples have the same distribution.
#' If yes, Wilcoxon-Mann-Whitney will be used to compare the distribution
#' difference.
#' @usage distributionTest(varIn, varOut, significanceLevel)
#' @param varIn first numeric vector
#' @param varOut second numeric vector
#' @param significanceLevel significant cutoff of the Kolmogorov-Smirnov test.
#' Default = 0.05.
#' @return p-value of the comparison test.
#' @author Carla Mölbert (carla.moelbert@gmx.de)

distributionTest <- function(
    varIn = NULL, varOut = NULL, significanceLevel = 0.05
){
    if (is.null(varIn) | is.null(varOut))
        stop("Vector of in-group and out-group cannot be NULL!")
    # remove NA values
    varIn <- varIn[!is.na(varIn)]
    varOut <- varOut[!is.na(varOut)]
    # if there is no data in one of the groups the p-value is NULL
    if (length(varIn) == 0 | length(varOut) == 0)
        stop("No data for In/Out-group taxa!")
    else {
        # * Kolmogorov-Smirnov Test
        # H0 : The two samples have the same distribution
        ks <- suppressWarnings(
            stats::ks.test(unique(varIn), unique(varOut), exact = FALSE)
        )

        if (ks$p.value <= significanceLevel) return(ks$p.value)
        else {
            # * Wilcoxon-Mann-Whitney Test
            # H0: the samples have the same location parameters
            wilcox <- suppressWarnings(
                stats::wilcox.test(
                    varIn, varOut, alternative = "two.sided", paired = FALSE
                )
            )
            return(wilcox$p.value)
        }
    }
}

#' Compare the median values of a variable between 2 taxon groups
#' @description Given the phylogenetic profiles that contains up to 2 additional
#' variables besides the presence/absence information of the orthologous
#' proteins. This function will compare the median scores of those variables
#' between 2 different taxon groups (e.g. parasitic species vs non-parasitic
#' species), which are defined as in-group and out-group. In-group is identified
#' by the user. Out-group contains all taxa in the input phylogenetic profiles
#' that are not part of the in-group.
#' @usage compareMedianTaxonGroups(data, inGroup, useCommonAncestor, variable,
#'     taxDB)
#' @export
#' @param data input phylogenetic profile in long format (see ?mainLongRaw and
#' ?createLongMatrix)
#' @param inGroup ID list of in-group taxa (e.g. "ncbi1234")
#' @param useCommonAncestor TRUE/FALSE if using all taxa that share the same
#' common ancestor with the pre-selected in-group as the in-group taxa.
#' Default = TRUE.
#' @param variable name of the variable that need to be compared
#' @param taxDB Path to the taxonomy DB files
#' @return List of genes that have a difference in the variable's median scores
#' between the in-group and out-group taxa and their corresponding delta-median.
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi9606", "ncbi10116")
#' variable <- colnames(data)[4]
#' compareMedianTaxonGroups(data, inGroup, TRUE, variable)

compareMedianTaxonGroups <- function(
    data = NULL, inGroup = NULL, useCommonAncestor = TRUE, variable = NULL,
    taxDB = NULL
) {
    if (is.null(data) | is.null(inGroup) | is.null(variable))
        stop("Input profiles, in-group IDs and variable name cannot be NULL!")
    if (!(variable %in% colnames(data))) stop("Invalid variable")
    # add other taxa that share a common ancestor with the given in-group
    commonTaxa <- getCommonAncestor(
        levels(as.factor(data$ncbiID)), inGroup, taxDB
    )
    if (useCommonAncestor == TRUE)
        inGroup <- as.character(commonTaxa[[3]]$abbrName)

    # return delta-median scores for two taxa groups
    data$geneID <- as.character(data$geneID)
    deltaMedian <- vapply(
        unique(data$geneID),
        function (x) {
            varIn <- data[data$geneID == x & data$ncbiID %in% inGroup, variable]
            varOut <- data[
                data$geneID == x & !(data$ncbiID %in% inGroup), variable]
            return(
                abs(
                    stats::median(varIn[!is.na(varIn)])
                    - stats::median(varOut[!is.na(varOut)])
                )
            )
        },
        FUN.VALUE = numeric(1)
    )
    return(sort(unlist(deltaMedian)))
}


#' Create data for variable distribution comparison plot
#' @description Create data for plotting the distribution comparison between 2
#' groups of taxa for a selected gene.
#' @usage dataVarDistTaxGroup(data, inGroup, gene, variable)
#' @export
#' @param data input phylogenetic profile in long format (see ?mainLongRaw and
#' ?createLongMatrix)
#' @param inGroup ID list of in-group taxa (e.g. "ncbi1234")
#' @param gene ID of gene that need to be plotted the distribution comparison
#' between in- and out-group taxa.
#' @param variable var1 or c(var1, var2)
#' @return Dataframe containing list of values for all available variables for
#' the selected genes in in-group and out-group taxa (max. 3 columns).
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{createLongMatrix}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi9606", "ncbi10116")
#' variable <- colnames(data)[c(4, 5)]
#' dataVarDistTaxGroup(data, inGroup, "101621at6656", variable)

dataVarDistTaxGroup <- function(
    data = NULL, inGroup = NULL, gene = NULL, variable = NULL
) {
    if (is.null(data) | is.null(inGroup) | is.null(gene) | is.null(variable))
        stop("Input profiles, in-group IDs, gene & variable ID cannot be NULL!")
    # remove "empty" variable (char "")
    variable <- variable[unlist(lapply(variable, function (x) x != ""))]
    # get 2 lists of values for in-group and out-group
    varIn <- data.frame(
        data[data$geneID == gene & data$ncbiID %in% inGroup, ][, variable]
    )
    colnames(varIn) <- variable
    varOut <- data.frame(
        data[data$geneID == gene & !(data$ncbiID %in% inGroup), ][, variable]
    )
    colnames(varOut) <- variable
    if (nrow(varIn) == 0 & nrow(varOut) == 0) stop("No data for in-/out-group!")

    varIn$type <- "In-group"
    varOut$type <- "Out-group"
    out <- rbind(
        varIn[stats::complete.cases(varIn),],
        varOut[stats::complete.cases(varOut),]
    )
    return(out[, c(variable, "type")])
}

#' Create a single violin distribution plot
#' @export
#' @param plotDf dataframe for plotting containing values for each variable in
#' in-group and out-group.
#' @param parameters plot parameters, including size of x-axis, y-axis,
#' legend and title; position of legend ("right", "bottom" or "none");
#' mean/median point; names of in-group and out-group; and plot title.
#' NOTE: Leave blank or NULL to use default values.
#' @param variable name of variable that need to be plotted (one of the column
#' names of input dataframe plotDf).
#' @return A violin plot as a ggplot object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @import ggplot2
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi9606", "ncbi10116")
#' varNames <- colnames(data)[c(4, 5)]
#' plotDf <- dataVarDistTaxGroup(data, inGroup, "101621at6656", varNames)
#' plotParameters <- list(
#'     "xSize" = 12,
#'     "ySize" = 12,
#'     "titleSize" = 15,
#'     "legendSize" = 12,
#'     "legendPosition" = "right",
#'     "mValue" = "mean",
#'     "inGroupName" = "In-group",
#'     "outGroupName" = "Out-group",
#'     "title" = "101621at6656"
#' )
#' generateSinglePlot(plotDf, plotParameters, colnames(plotDf)[1])

generateSinglePlot <- function(plotDf, parameters, variable) {
    type <- .data <- NULL
    xNames <- c(
        paste(
            parameters$inGroupName, " \n n = ",
            nrow(plotDf[plotDf$type == parameters$inGroupName,]), sep = ""
        ),
        paste(
            parameters$outGroupName, " \n n = ",
            nrow(plotDf[plotDf$type == parameters$outGroupName,]), sep = ""
        )
    )

    plot <- ggplot(plotDf, aes(x = factor(type), y = .data[[variable]])) +
        geom_violin(
            aes(fill = factor(type)), position = position_dodge(),
            scale = "width", alpha = .5) +
        geom_boxplot(width = 0.1) +
        scale_x_discrete(labels = xNames) +
        # add mean/median point and color of that point
        stat_summary(
            aes(colour = parameters$mValue), fun = parameters$mValue,
            geom = "point", size = 3, show.legend = TRUE, shape = 8
        ) +
        scale_color_manual("", values = c("red")) +
        # add title and theme
        labs(x = element_blank(), y = variable) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(size = parameters$xSize, hjust = 1),
            axis.text.y = element_text(size = parameters$ySize),
            axis.title.y = element_text(size = parameters$ySize),
            legend.position = parameters$legendPosition,
            legend.text = element_text(size = parameters$legendSize),
            legend.title = element_blank()
        )
    return(plot)
}

#' Plot Multiple Graphs with Shared Legend in a Grid
#' @export
#' @usage gridArrangeSharedLegend(...,  ncol = length(list(...)), nrow = 1,
#'     position = c("bottom", "right"), title = NA, titleSize = 12)
#' @param ...  Plots to be arranged in grid
#' @param ncol Number of columns in grid
#' @param nrow Number of rows in grid
#' @param position Gird position (bottom or right)
#' @param title Title of grid
#' @param titleSize Size of grid title
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @return Grid of plots with common legend
#' @author Phil Boileau, \email{philippe.boileau@rimuhc.ca}
#' @note adapted from https://rdrr.io/github/PhilBoileau/CLSAR/src/R/
#' gridArrangeSharedLegend.R
#' @examples
#' \dontrun{
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi9606", "ncbi10116")
#' varNames <- colnames(data)[c(4, 5)]
#' plotDf <- dataVarDistTaxGroup(data, inGroup, "101621at6656", varNames)
#' plotParameters <- list(
#'     "xSize" = 12,
#'     "ySize" = 12,
#'     "titleSize" = 15,
#'     "legendSize" = 12,
#'     "legendPosition" = "right",
#'     "mValue" = "mean",
#'     "inGroupName" = "In-group",
#'     "outGroupName" = "Out-group",
#'     "title" = "101621at6656"
#' )
#' plotVar1 <- generateSinglePlot(plotDf, plotParameters, colnames(plotDf)[1])
#' plotVar2 <- generateSinglePlot(plotDf, plotParameters, colnames(plotDf)[2])
#' g <- gridArrangeSharedLegend(
#'     plotVar1, plotVar2,
#'     position = plotParameters$legendPosition,
#'     title = plotParameters$title,
#'     size = plotParameters$titleSize
#' )
#' }

gridArrangeSharedLegend <- function(
    ...,  ncol = length(list(...)), nrow = 1, position = c("bottom", "right"),
    title = NA, titleSize = 12
) {
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[
        which(vapply(g, function(x) x$name, FUN.VALUE = "") == "guide-box")
        ]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(
        position,
        "bottom" = gridExtra::arrangeGrob(
            do.call(arrangeGrob, gl),
            legend,
            ncol = 1,
            heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight),
            top = grid::textGrob(
                title,
                vjust = 1,
                gp = grid::gpar(
                    fontface = "bold", fontsize = titleSize
                )
            )
        ),
        "right" = gridExtra::arrangeGrob(
            do.call(arrangeGrob, gl),
            legend,
            ncol = 2,
            widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth),
            top = grid::textGrob(
                title,
                vjust = 1,
                gp = grid::gpar(
                    fontface = "bold", fontsize = titleSize
                )
            )
        )
    )
    return(combined)
}

#' Create variable distribution comparison plot
#' @description Create variable distribution plots between 2 groups of taxa for
#' a selected gene.
#' @export
#' @param data dataframe for plotting. Last column
#' indicates what type of taxon group (in- or out-group). The first (or first 2)
#' column contains values of the variables. See ?dataVarDistTaxGroup
#' @param plotParameters plot parameters, including size of x-axis, y-axis,
#' legend and title; position of legend ("right", "bottom" or "none");
#' mean/median point; names of in-group and out-group; and plot title.
#' NOTE: Leave blank or NULL to use default values.
#' @return Distribution plots as a grob (gtable) object. Use grid.draw to plot.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{dataVarDistTaxGroup}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi9606", "ncbi10116")
#' variable <- colnames(data)[c(4, 5)]
#' plotDf <- dataVarDistTaxGroup(data, inGroup, "101621at6656", variable)
#' plotParameters <- list(
#'     "xSize" = 12,
#'     "ySize" = 12,
#'     "titleSize" = 15,
#'     "legendSize" = 12,
#'     "legendPosition" = "right",
#'     "mValue" = "mean",
#'     "inGroupName" = "In-group",
#'     "outGroupName" = "Out-group",
#'     "title" = "101621at6656"
#' )
#' g <- varDistTaxPlot(plotDf, plotParameters)
#' grid::grid.draw(g)

varDistTaxPlot <- function(data, plotParameters) {
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (missing(plotParameters)) stop("Plot parameters are missing!")
    # rename in-group and out-group
    data$type[data$type == "In-group"] <- plotParameters$inGroupName
    data$type[data$type == "Out-group"] <- plotParameters$outGroupName
    # return plot(s)
    if (ncol(data) == 2) {
        plotVar1 <- generateSinglePlot(data, plotParameters, colnames(data)[1])
        return(
            arrangeGrob(
                plotVar1,
                top = grid::textGrob(
                    plotParameters$title, vjust = 1,
                    gp = grid::gpar(
                        fontface = "bold", fontsize = plotParameters$titleSize
                    )
                )
            )
        )
    } else {
        plotVar1 <- generateSinglePlot(data, plotParameters, colnames(data)[1])
        plotVar2 <- generateSinglePlot(data, plotParameters, colnames(data)[2])
        if (plotParameters$legendPosition == "none") {
            return(
                arrangeGrob(
                    plotVar1, plotVar2,
                    nrow = 1,
                    top = grid::textGrob(
                        plotParameters$title, vjust = 1,
                        gp = grid::gpar(
                            fontface = "bold",
                            fontsize = plotParameters$titleSize
                        )
                    )
                )
            )
        } else {
            return(
                gridArrangeSharedLegend(
                    plotVar1, plotVar2,
                    position = plotParameters$legendPosition,
                    title = plotParameters$title
                )
            )
        }
    }
}

#' Create data for feature distribution comparison plot
#' @description Create data for plotting the distribution of the protein domain
#' features between 2 group of taxa for a selected gene (average number of
#' feature occurrency per protein/ortholog).
#' @usage dataFeatureTaxGroup(mainDf, domainDf, inGroup, gene)
#' @export
#' @param mainDf input phylogenetic profile in long format (see ?mainLongRaw
#' and ?createLongMatrix)
#' @param domainDf dataframe contains domain info for the seed and ortholog.
#' This including the seed ID, orthologs IDs, sequence lengths, feature names,
#' start and end positions, feature weights (optional) and the status to
#' determine if that feature is important for comparison the architecture
#' between 2 proteins* (e.g. seed protein vs ortholog) (optional). (see
#' ?parseDomainInput)
#' @param inGroup ID list of in-group taxa (e.g. "ncbi1234")
#' @param gene ID of gene that need to be plotted the feature distribution
#' comparison between in- and out-group taxa.
#' @return Dataframe containing all feature names, their frequencies (absolute
#' count and the average instances per protein - IPP) in each taxon group and
#' the corresponding taxa group type (in- or out-group).
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{createLongMatrix}}, \code{\link{parseDomainInput}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' mainDf <- mainLongRaw
#' gene <- "101621at6656"
#' inputFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' type <- "file"
#' domainDf <- parseDomainInput(gene, inputFile, type)
#' inGroup <- c("ncbi9606", "ncbi10116")
#' dataFeatureTaxGroup(mainDf, domainDf, inGroup, gene)

dataFeatureTaxGroup <- function(
    mainDf = NULL, domainDf = NULL, inGroup = NULL, gene = NULL
) {
    if (is.null(mainDf) | is.null(inGroup) | is.null(gene) | is.null(domainDf))
        stop("Input profiles, in-group IDs, gene & domain data cannot be NULL!")
    # get ncbiIDs for the domain data
    mainDf$orthoID <- gsub("\\|", ":", mainDf$orthoID)
    domainDfSub <- merge(
        domainDf[grep(paste0(gene, "#"), domainDf$seedID),][
            ,c("orthoID","seedID","feature")
        ],
        mainDf[, c("orthoID", "ncbiID")], by = "orthoID", all.x = TRUE
    )
    domainDfSub <- domainDfSub[stats::complete.cases(domainDfSub),]
    # identify in-group and out-group
    domainDfSub$type[domainDfSub$ncbiID %in% inGroup] <- "In-group"
    domainDfSub$type[!(domainDfSub$ncbiID %in% inGroup)] <- "Out-group"
    # count number of orthologs for each taxon group
    nInGroup <- length(
        unique(domainDfSub$seedID[domainDfSub$type == "In-group"]))
    nOutGroup <- length(
        unique(domainDfSub$seedID[domainDfSub$type == "Out-group"]))
    # get intances and count the number of intances for each taxon group
    domainIn <- domainDfSub$feature[domainDfSub$ncbiID %in% inGroup]
    domainOut <- domainDfSub$feature[!(domainDfSub$ncbiID %in% inGroup)]
    countDomainIn <- data.frame(table(unlist(as.character(domainIn))))
    countDomainIn$type <- "In-group"
    countDomainIn$ipp <- countDomainIn$Freq/nInGroup
    countDomainOut <- data.frame(table(unlist(as.character(domainOut))))
    countDomainOut$type <- "Out-group"
    countDomainOut$ipp <- countDomainOut$Freq/nOutGroup
    # calculate delta IPP
    mergedDf <- merge(countDomainIn, countDomainOut, by = "Var1", all = TRUE)
    mergedDf[is.na(mergedDf)] <- 0
    mergedDf$dIPPtmp <- mergedDf$ipp.x - mergedDf$ipp.y
    mergedDf$dIPP <- mergedDf$dIPPtmp / (mergedDf$ipp.x + mergedDf$ipp.y)
    # return
    outDf <- merge(
        rbind(countDomainIn, countDomainOut), mergedDf[, c("Var1", "dIPP")],
        by = "Var1", all.x = TRUE)
    colnames(outDf) <- c("Feature", "Count", "Taxon_group", "IPP", "dIPP")
    return(outDf)
}

#' Create feature distribution comparison plot
#' @description Create protein feature distribution plots between 2 groups of
#' taxa for a selected gene.
#' @export
#' @param data dataframe for plotting (see ?dataFeatureTaxGroup)
#' @param plotParameters plot parameters, including size of x-axis, y-axis,
#' legend and title; position of legend ("right", "bottom" or "none"); names of
#' in-group and out-group; flip the plot coordinate ("Yes" or "No").
#' NOTE: Leave blank or NULL to use default values.
#' @return Distribution plots as a ggplot2 object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{dataFeatureTaxGroup}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' gene <- "101621at6656"
#' inputFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' type <- "file"
#' domainDf <- parseDomainInput(gene, inputFile, type)
#' inGroup <- c("ncbi9606", "ncbi10116")
#' plotDf <- dataFeatureTaxGroup(data, domainDf, inGroup, gene)
#' plotParameters <- list(
#'     "xSize" = 12,
#'     "ySize" = 12,
#'     "angle" = 15,
#'     "legendSize" = 12,
#'     "inGroupName" = "In-group",
#'     "outGroupName" = "Out-group",
#'     "flipPlot" = "No"
#' )
#' featureDistTaxPlot(plotDf, plotParameters)

featureDistTaxPlot <- function(data, plotParameters) {
    Feature <- NULL
    IPP <- NULL
    Taxon_group <- NULL
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (missing(plotParameters)) stop("Plot parameters are missing!")

    data$Taxon_group[data$Taxon_group == "In-group"] <-
        plotParameters$inGroupName
    data$Taxon_group[data$Taxon_group == "Out-group"] <-
        plotParameters$outGroupName

    plot <- ggplot(data, aes(x = Feature, y = IPP, fill = Taxon_group)) +
        geom_bar(stat="identity", width=.5, position = "dodge") +
        theme_minimal() +
        theme(
            axis.title.y = element_text(size = plotParameters$ySize),
            axis.text.y = element_text(size = plotParameters$ySize),
            axis.title.x = element_text(size = plotParameters$xSize),
            axis.text.x = element_text(
                size = plotParameters$xSize,
                angle = plotParameters$angle,
                hjust = 1
            ),
            legend.text = element_text(size = plotParameters$legendSize),
            legend.title = element_blank()
        )
    if (plotParameters$flipPlot == "Yes") plot <- plot + coord_flip()
    return(plot)
}
