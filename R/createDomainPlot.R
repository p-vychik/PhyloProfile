#' Create protein's domain architecure plot
#' @description Create architecture plot for both seed and orthologous protein.
#' If domains of ortholog are missing, only architecture of seed protein will
#' be plotted. NOTE: seed protein ID is the one being shown in the profile plot,
#' which normally is also the orthologous group ID.
#' @export
#' @usage createArchiPlot(info, domainDf, labelArchiSize, titleArchiSize,
#'     showScore, domainNamePos, firstDist, nameType, nameSize, nameColor, labelPos,
#'     colorType, ignoreInstanceNo, currentNCBIinfo, featureTypeSort,
#'     featureTypeOrder, colorPallete, resolveOverlap)
#' @param info A list contains seed and ortholog's IDs
#' @param domainDf Dataframe contains domain info for the seed and ortholog.
#' This including the seed ID, orthologs IDs, sequence lengths, feature names,
#' start and end positions, feature weights (optional) and the status to
#' determine if that feature is important for comparison the architecture
#' between 2 proteins* (e.g. seed protein vs ortholog) (optional).
#' @param labelArchiSize Lable size (in px). Default = 12.
#' @param titleArchiSize Title size (in px). Default = 12.
#' @param showScore Show/hide E-values and Bit-scores. Default = NULL (hide)
#' @param domainNamePos list of positions for domain names, choose from "plot",
#' "legend" or "axis". Default: "plot"
#' @param firstDist Distance of the first domain to plot title. Default = 0.5
#' @param nameType Type of domain names, either "Texts" or "Labels" (default)
#' @param nameSize Size of domain names (for Texts only). Default = 5
#' @param nameColor Color of domain names (for Texts only). Default = "black"
#' @param labelPos Position of doman names (for Labels only). Choose from
#' @param colorType Choose to color "all", "shared", "unique" features or color
#' by "Feature type". Default = "all"
#' @param ignoreInstanceNo Ignore number of feature instances while identifying
#' shared or unique features. Default = FALSE
#' @param currentNCBIinfo Dataframe of the pre-processed NCBI taxonomy
#' data. Default = NULL (will be automatically retrieved from PhyloProfile app)
#' @param featureTypeSort Choose to arrange plot by shared features. Default =
#' "Yes"
#' @param featureTypeOrder Vector of ordered feature types
#' @param colorPallete Choose between "Paired", "Set1", "Set2", "Set3",
#' "Accent", "Dark2" for the color pallete
#' @param resolveOverlap Choose to merge non-overlapped features of a feature
#' type into one line. Default = "Yes"
#' @return A domain plot as arrangeGrob object. Use grid::grid.draw(plot) to
#' render.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{singleDomainPlotting}}, \code{\link{sortDomains}},
#' \code{\link{parseDomainInput}}, \code{\link{getQualColForVector}}
#' @examples
#' seedID <- "101621at6656"
#' orthoID <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' info <- c(seedID, orthoID)
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' plot <- createArchiPlot(info, domainDf, 9, 9, seqIdFormat = "bionf")
#' grid::grid.draw(plot)

createArchiPlot <- function(
        info = NULL, domainDf = NULL, labelArchiSize = 12, titleArchiSize = 12,
        showScore = NULL, domainNamePos = "plot", firstDist = 0.5,
        nameType = "Labels", nameSize = 5, nameColor = "#000000", labelPos = "Above",
        colorType = "Unique", ignoreInstanceNo = FALSE, currentNCBIinfo = NULL,
        featureTypeSort = "Yes", featureTypeOrder = NULL, colorPallete = "Paired",
        resolveOverlap = "Yes"
){
    if (is.null(info) | is.null(domainDf)) return(ggplot() + theme_void())
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    subdomainDf <- domainDf[domainDf$seedID == grepID, ]
    subdomainDf$feature <- as.character(subdomainDf$feature)
    orthoID <- NULL

    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
    else {
        # get minStart and maxEnd
        minStart <- min(subdomainDf$start)
        maxEnd <- max(subdomainDf$end)
        if ("length" %in% colnames(subdomainDf))
            maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        if (nrow(seedDf) == 0) seedDf <- orthoDf
        seed <- as.character(seedDf$orthoID[1])
        if (nrow(seedDf) == 0) return(paste0("No domain info available!"))

        # simplify seed/ortho seq IDs if they are in bionf format
        if (!is.null(currentNCBIinfo)) {
            if (str_count(seed, ":") >= 2 & str_count(seed, "@") >= 2) {
                seedTmp <- strsplit(as.character(seed),':', fixed = TRUE)[[1]]
                seedSpec <-
                    strsplit(as.character(seedTmp[2]),'@', fixed = TRUE)[[1]][2]
                seed <- paste0(
                    id2name(seedSpec, currentNCBIinfo)[,2], " - ", seedTmp[3]
                )
                if (ortho != seed) {
                    orthoTmp <- strsplit(as.character(ortho),':', fixed = TRUE)[[1]]
                    orthoSpec <-
                        strsplit(as.character(orthoTmp[2]),'@',fixed = TRUE)[[1]][2]
                    ortho <- paste0(
                        id2name(orthoSpec, currentNCBIinfo)[,2], " - ", orthoTmp[3]
                    )
                }
            }
        }

        # add feature colors
        featureColorDf <- addFeatureColors(seedDf, orthoDf, colorType, colorPallete, ignoreInstanceNo)
        seedDf <- featureColorDf[[1]]
        orthoDf <- featureColorDf[[2]]
        # resolve (non)overlapped features
        if (resolveOverlap == "Yes") {
            seedDf <- resolveOverlapFeatures(seedDf)
            orthoDf <- resolveOverlapFeatures(orthoDf)
        } else {
            seedDf$featureOri <- seedDf$feature
            seedDf$featureOri <- as.character(seedDf$featureOri)
            orthoDf$featureOri <- orthoDf$feature
            orthoDf$featureOri <- as.character(orthoDf$featureOri)
        }

        if (nrow(orthoDf) > 0) {
            if (all.equal(seedDf, orthoDf)[1] == TRUE) featureTypeSort <- "No"
            # sort features
            if (featureTypeSort == "Yes") {
                # change order of one df's features based on order of other df's
                if (length(orthoDf$feature) < length(seedDf$feature)) {
                    orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
                    orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
                    orderedOrthoDf <- sortDomains(orderedSeedDf, orderedOrthoDf)
                } else {
                    orderedSeedDf <- seedDf[order(seedDf$feature), ]
                    orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
                    orderedSeedDf <- sortDomains(orderedOrthoDf, orderedSeedDf)
                }
            } else {
                # change order based on list of feature types
                orderedSeedDf <- sortDomainsByList(seedDf, featureTypeOrder)
                orderedOrthoDf <- sortDomainsByList(orthoDf, featureTypeOrder)
            }

            # join weight values and feature names
            if ("weight" %in% colnames(orderedOrthoDf)) {
                NULL
                # orderedOrthoDf$yLabel <- paste0(
                #     orderedOrthoDf$feature," (",round(orderedOrthoDf$weight, 2),")")
            } else orderedOrthoDf$yLabel <- orderedOrthoDf$feature
            if ("weight" %in% colnames(orderedSeedDf)) {
                NULL
                # orderedSeedDf$yLabel <- paste0(
                #     orderedSeedDf$feature," (",round(orderedSeedDf$weight, 2),")")
            } else orderedSeedDf$yLabel <- orderedSeedDf$feature
            # plotting
            g <- pairDomainPlotting(
                seed, ortho, orderedSeedDf, orderedOrthoDf, minStart, maxEnd,
                labelArchiSize, titleArchiSize, showScore, domainNamePos, firstDist,
                nameType, nameSize, nameColor, labelPos, colorPallete)
        } else {
            # orderedSeedDf <- seedDf[order(seedDf$feature), ]
            orderedSeedDf <- sortDomainsByList(seedDf, featureTypeOrder)
            if ("weight" %in% colnames(orderedSeedDf)) {
                NULL
                # orderedSeedDf$yLabel <- paste0(
                #     orderedSeedDf$feature," (",round(orderedSeedDf$weight, 2),")")
            } else orderedSeedDf$yLabel <- orderedSeedDf$feature
            # plotting
            g <- pairDomainPlotting(
                seed, seed, orderedSeedDf, orderedSeedDf, minStart, maxEnd,
                labelArchiSize, titleArchiSize, showScore, domainNamePos, firstDist,
                nameType, nameSize, nameColor, labelPos, colorPallete)
        }
        return(g)
    }
}


#' Create protein's domain architecure plot
#' @description Create architecture plot for both seed and orthologous protein.
#' If domains of ortholog are missing, only architecture of seed protein will
#' be plotted. NOTE: seed protein ID is the one being shown in the profile plot,
#' which normally is also the orthologous group ID.
#' @export
#' @usage createArchiPlot2(info = NULL, domainDf = NULL, labelArchiSize = 12,
#'     titleArchiSize = 12, showFeature = "all", seqIdFormat = "unknown",
#'     currentNCBIinfo = NULL)
#' @param info a list contains seed and ortholog's IDs
#' @param domainDf dataframe contains domain info for the seed and ortholog.
#' This including the seed ID, orthologs IDs, sequence lengths, feature names,
#' start and end positions, feature weights (optional) and the status to
#' determine if that feature is important for comparison the architecture
#' between 2 proteins* (e.g. seed protein vs ortholog) (optional).
#' @param labelArchiSize lable size (in px). Default = 12.
#' @param titleArchiSize title size (in px). Default = 12.
#' @param showFeature choose to show all, common or unique features.
#' Default = "all"
#' @param seqIdFormat sequence ID format (either bionf or unknown).
#' Default = "unknown"
#' @param currentNCBIinfo dataframe of the pre-processed NCBI taxonomy
#' data. Default = NULL (will be automatically retrieved from PhyloProfile app)
#' @importFrom gridExtra arrangeGrob
#' @import ggplot2
#' @return A domain plot as arrangeGrob object. Use grid::grid.draw(plot) to
#' render.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{singleDomainPlotting}}, \code{\link{sortDomains}},
#' \code{\link{parseDomainInput}}, \code{\link{getQualColForVector}}
#' @examples
#' seedID <- "101621at6656"
#' orthoID <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' info <- c(seedID, orthoID)
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' plot <- createArchiPlot(info, domainDf, 9, 9, seqIdFormat = "bionf")
#' grid::grid.draw(plot)

createArchiPlot2 <- function(
    info = NULL, domainDf = NULL, labelArchiSize = 12, titleArchiSize = 12,
    showFeature = "all", seqIdFormat = "unknown", currentNCBIinfo = NULL
){
    if (is.null(info) | is.null(domainDf)) return(ggplot() + theme_void())
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
    subdomainDf$feature <- as.character(subdomainDf$feature)
    orthoID <- NULL
    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
    else {
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        # filter common features
        if (!(showFeature == "all")) {
            allFeats <- c(
                levels(as.factor(orthoDf$feature)),
                levels(as.factor(seedDf$feature))
            )
            countFeats <- as.data.frame(table(allFeats))
            commonFeats <- countFeats$allFeats[countFeats$Freq > 1]
            if (showFeature == "common") {
                orthoDf <- orthoDf[orthoDf$feature %in% commonFeats,]
                seedDf <- seedDf[seedDf$feature %in% commonFeats,]
            } else {
                orthoDf <- orthoDf[!(orthoDf$feature %in% commonFeats),]
                seedDf <- seedDf[!(seedDf$feature %in% commonFeats),]
            }
        }
        # final check
        if (nrow(seedDf) == 0) seedDf <- orthoDf
        seed <- as.character(seedDf$orthoID[1])
        if (nrow(seedDf) == 0) return(paste0("No domain info available!"))
        if (nrow(orthoDf) == 0) {
            ortho <- seed
            orthoDf <- seedDf
        }
        # change order of one df's features based on order of other df's
        if (length(orthoDf$feature) < length(seedDf$feature)) {
            if (nrow(orthoDf) > 0) {
                orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
                orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
            }
        } else {
            orderedSeedDf <- seedDf[order(seedDf$feature), ]
            orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
        }
        # modify feature (domain) names (used as x-axis label in domain plot)
        orderedOrthoDf <- modifyFeatureName(orderedOrthoDf)
        orderedSeedDf <- modifyFeatureName(orderedSeedDf)
        # simplify seq IDs if they are in bionf format
        if (seqIdFormat == "bionf") {
            seedTmp <- strsplit(as.character(seed),':', fixed = TRUE)[[1]]
            seedSpec <-
                strsplit(as.character(seedTmp[2]),'@', fixed = TRUE)[[1]][2]
            seed <- paste0(
                id2name(seedSpec, currentNCBIinfo)[,2], " - ", seedTmp[3]
            )
            if (ortho != seed) {
                orthoTmp <- strsplit(as.character(ortho),':', fixed = TRUE)[[1]]
                orthoSpec <-
                    strsplit(as.character(orthoTmp[2]),'@',fixed = TRUE)[[1]][2]
                ortho <- paste0(
                    id2name(orthoSpec, currentNCBIinfo)[,2], " - ", orthoTmp[3]
                )
            }
        }
        # plotting
        minStart <- min(subdomainDf$start)
        maxEnd <- max(subdomainDf$end)
        if ("length" %in% colnames(subdomainDf))
            maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
        g <- pairDomainPlotting(
            seed, ortho, orderedSeedDf, orderedOrthoDf, minStart, maxEnd,
            labelArchiSize, titleArchiSize)
        return(g)
    }
}

#' Create architecure plot for a single protein
#' @usage singleDomainPlotting(df, geneID, sep, labelSize, titleSize, minStart,
#'     maxEnd, colorPallete, showScore, domainNamePos, firstDist, nameType,
#'     nameSize, nameColor, labelPos)
#' @param df Domain dataframe for ploting containing the seed ID, ortholog ID,
#' ortholog sequence length, feature names, start and end positions,
#' feature weights (optional) and the status to determine if that feature is
#' important for comparison the architecture between 2 proteins* (e.g. seed
#' protein vs ortholog) (optional)
#' @param geneID ID of seed or orthologous protein
#' @param sep Separate indicator for title. Default = "|"
#' @param labelSize Lable size. Default = 12
#' @param titleSize Title size. Default = 12
#' @param minStart The smallest start position of all domains
#' @param maxEnd The highest stop position of all domains
#' @param colorPallete Color pallete. Default = Paired"
#' @param showScore Show/hide E-values and Bit-scores. Default = NULL (hide)
#' @param domainNamePos List of positions for domain names, choose from "plot",
#' "legend" or "axis". Default: "plot"
#' @param firstDist Distance of the first domain to plot title. Default = 0.5
#' @param nameType Type of domain names, either "Texts" or "Labels" (default)
#' @param nameSize Size of domain names (for Texts only). Default = 5
#' @param nameColor Color of domain names (for Texts only). Default = "black"
#' @param labelPos Position of doman names (for Labels only). Choose from
#' "Above" (default), "Below" or "Inside" the domain bar
#' @return Domain plot of a single protein as a ggplot object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{getQualColForVector}},
#' \code{\link{parseDomainInput}}
#' @import ggplot2
#' @examples
#' \dontrun{
#' # get domain data
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' seedID <- "101621at6656"
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' df <- domainDf[
#'     domainDf$orthoID == "101621at6656:AGRPL@224129@0:224129_0:001955:1",]
#' # create color scheme for all domain types
#' allFeatures <- levels(as.factor(df$feature))
#' allColors <- getQualColForVector(allFeatures)
#' colorScheme <- structure(
#'     allColors,
#'     .Names = allFeatures
#' )
#' # other parameters
#' geneID <- "AGRPL@224129@0|224129_0:001955|1"
#' sep <- "|"
#' labelSize <- 9
#' titleSize <- 9
#' minStart <- min(df$start)
#' maxEnd <- max(df$end)
#' # do plotting
#' PhyloProfile:::singleDomainPlotting(
#'     df,
#'     geneID,
#'     sep,
#'     labelSize, titleSize,
#'     minStart, maxEnd,
#'     colorScheme
#' )
#' }

singleDomainPlotting <- function(
        df = NULL, geneID = "GeneID", sep = "|", labelSize = 12, titleSize = 12,
        minStart = NULL, maxEnd = NULL, colorPallete = "Set2",
        showScore = NULL, domainNamePos = "plot", firstDist = 0.5,
        nameType = "Labels", nameSize = 5, nameColor = "#000000", labelPos = "Above"
){
    feature <- feature_id_mod <- end <- start <- NULL

    # parse parameters
    if (is.null(df)) return(ggplot() + theme_void())
    if (is.null(minStart)) minStart <- min(df$start)
    if (is.null(maxEnd)) maxEnd <- max(df$end)
    if ("color" %in% colnames(df)) {
        colorScheme <- structure(
            df$color, .Names = df$featureOri
        )
    } else {
        colorScheme <- structure(
            head(
                suppressWarnings(
                    RColorBrewer::brewer.pal(
                        nlevels(as.factor(df$feature)), colorPallete
                    )
                ),
                levels(as.factor(df$feature))
            ),
            .Names = levels(as.factor(df$feature)))
    }

    # initiate ggplot object
    gg <- ggplot(df, aes(y = feature, x = end))

    # draw lines for representing sequence length
    if ("length" %in% colnames(df))
        gg <- gg + geom_segment(
            data = df, size = 1, color = "#b2b2b2", alpha = 0.0,
            aes(x = 0, xend = length, y = feature, yend = feature))
    # draw features
    gg <- gg + geom_segment(
        data = df, aes(x = start, xend = end, y = feature, yend = feature, color = as.factor(featureOri)),
        linewidth = nameSize, lineend = "round", alpha = 0.7) +
        scale_color_manual(values = colorScheme)

    # add feature names
    if ("plot" %in% domainNamePos) {
        if (nameType == "Labels") {
            if (labelPos == "Above") {
                gg <- gg + geom_label(
                    aes(label = str_wrap(feature_id_mod), x = (start+end)/2),
                    color = "black", vjust = -0.5
                )
            } else if (labelPos == "Below") {
                gg <- gg + geom_label(
                    aes(label = str_wrap(feature_id_mod), x = (start+end)/2),
                    color = "black", vjust = 1.5
                )
            } else {
                gg <- gg + geom_label(
                    aes(label = str_wrap(feature_id_mod), x = (start+end)/2),
                    color = "black", size = nameSize - 2
                )
            }
        } else {
            gg <- gg + geom_text(
                aes(label = str_wrap(feature_id_mod),
                    x = (start+end)/2),
                color = nameColor, check_overlap = TRUE, size = nameSize - 2
            )
        }
    }
    # add scores if selected
    if ("Bit-score" %in% showScore & "E-value" %in% showScore) {
        gg <- gg + geom_label(
            aes(
                label = ifelse(evalue == "NA" & bitscore == "NA", "", str_wrap(
                    paste0(
                        "E-value: ", evalue, "; Bitscore: ", bitscore
                    )
                )),
                x = (start+end)/2
            ),
            color = "black", vjust = 1.25
        )
    } else if ("Bit-score" %in% showScore) {
        gg <- gg + geom_label(
            aes(
                label = ifelse(bitscore == "NA", "", str_wrap(
                    paste0(
                        "Bitscore: ", bitscore
                    )
                )),
                x = (start+end)/2
            ),
            color = "black", vjust = 1.25
        )
    } else if ("E-value" %in% showScore) {
        gg <- gg + geom_label(
            aes(
                label = ifelse(evalue == "NA", "", str_wrap(
                    paste0(
                        "E-value: ", evalue
                    )
                )),
                x = (start+end)/2
            ),
            color = "black", vjust = 1.25
        )
    }
    # theme format
    gg <- gg + labs(
        title = paste0(gsub(":", sep, geneID)), color = "Feature"
    )
    gg <- gg + theme_minimal() + theme(panel.border = element_blank())
    # gg <- gg + theme(axis.ticks = element_blank())
    gg <- gg + theme(plot.title = element_text(face = "bold", size = titleSize))
    gg <- gg + theme(plot.title = element_text(hjust = 0.5))
    gg <- gg + theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
    # add feature names on the axis (if required)
    if ("axis" %in% domainNamePos) {
        if ("plot" %in% domainNamePos | "legend" %in% domainNamePos) {
            gg <- gg + scale_y_discrete(
                expand = c(0.075, 0), breaks = df$feature, labels = df$feature_type)
        } else {
            gg <- gg + scale_y_discrete(
                expand = c(0.075, 0), breaks = df$feature, labels = df$feature)
        }
        gg <- gg + theme(axis.text.y = element_text(size = labelSize))
    } else {
        gg <- gg + theme(axis.text.y = element_blank())
    }
    # add legend (if required)
    if ("legend" %in% domainNamePos) {
        gg <- gg + theme(legend.position = "bottom")
    } else {
        gg <- gg + theme(legend.position = "none")
    }
    # add space on the top of the plot (for feature name)
    gg <- gg + coord_cartesian(
        clip = 'off', ylim = c(1, nlevels(as.factor(df$feature)) + firstDist)
    )
    return(gg)
}


#' Create architecure plot for a single protein
#' @usage singleDomainPlotting2(df, geneID = "GeneID", sep = "|", labelSize = 12,
#'     titleSize = 12, minStart = NULL, maxEnd = NULL, colorScheme)
#' @param df domain dataframe for ploting containing the seed ID, ortholog ID,
#' ortholog sequence length, feature names, start and end positions,
#' feature weights (optional) and the status to determine if that feature is
#' important for comparison the architecture between 2 proteins* (e.g. seed
#' protein vs ortholog) (optional).
#' @param geneID ID of seed or orthologous protein
#' @param sep separate indicator for title. Default = "|".
#' @param labelSize lable size. Default = 12.
#' @param titleSize title size. Default = 12.
#' @param minStart the smallest start position of all domains
#' @param maxEnd the highest stop position of all domains
#' @param colorScheme color scheme for all domain types
#' @return Domain plot of a single protein as a ggplot object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{getQualColForVector}},
#' \code{\link{parseDomainInput}}
#' @import ggplot2
#' @examples
#' \dontrun{
#' # get domain data
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' seedID <- "101621at6656"
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' df <- domainDf[
#'     domainDf$orthoID == "101621at6656:AGRPL@224129@0:224129_0:001955:1",]
#' # create color scheme for all domain types
#' allFeatures <- levels(as.factor(df$feature))
#' allColors <- getQualColForVector(allFeatures)
#' colorScheme <- structure(
#'     allColors,
#'     .Names = allFeatures
#' )
#' # other parameters
#' geneID <- "AGRPL@224129@0|224129_0:001955|1"
#' sep <- "|"
#' labelSize <- 9
#' titleSize <- 9
#' minStart <- min(df$start)
#' maxEnd <- max(df$end)
#' # do plotting
#' PhyloProfile:::singleDomainPlotting(
#'     df,
#'     geneID,
#'     sep,
#'     labelSize, titleSize,
#'     minStart, maxEnd,
#'     colorScheme
#' )
#' }

singleDomainPlotting2 <- function(
    df = NULL, geneID = "GeneID", sep = "|", labelSize = 12, titleSize = 12,
    minStart = NULL, maxEnd = NULL, colorScheme = NULL
){
    feature <- end <- start <- NULL
    # parse parameters
    if (is.null(df)) return(ggplot() + theme_void())
    if (is.null(minStart)) minStart <- min(df$start)
    if (is.null(maxEnd)) maxEnd <- max(df$end)
    if (is.null(colorScheme)) {
        colorScheme <- structure(
            getQualColForVector(levels(as.factor(df$feature))),
            .Names = levels(as.factor(df$feature)))}
    gg <- ggplot(df, aes(y = feature, x = end, color = as.factor(feature))) +
        geom_segment(
            data = df, color = "white", linewidth = 0,
            aes(y = feature, yend = feature, x = minStart, xend = maxEnd)) +
        scale_color_manual(values = colorScheme)
    # draw lines for representing sequence length
    if ("length" %in% colnames(df))
        gg <- gg + geom_segment(
            data = df, linewidth = 1, color = "#b2b2b2",
            aes(x = 0, xend = length, y = feature, yend = feature))
    # draw line and points
    gg <- gg + geom_segment(
        data = df, aes(x = start, xend = end, y = feature, yend = feature),
        linewidth = 1.5) +
        geom_point(data = df, aes(y = feature, x = start),
                    color = "#b2b2b2", size = 3, shape = 3) +
        geom_point(data = df, aes(y = feature, x = end),
                    color = "#edae52", size = 3, shape = 5)
    # draw dashed line for domain path
    gg <- gg + geom_segment(
        data = df[df$path == "Y", ], linewidth = 3, linetype = "dashed",
        aes(x = start, xend = end, y = feature, yend = feature))
    # theme format
    gg <- gg + scale_y_discrete(
        expand = c(0.075, 0), breaks = df$feature, labels = df$yLabel)
    gg <- gg + labs(title = paste0(gsub(":", sep, geneID)), y = "Feature")
    gg <- gg + theme_minimal() + theme(panel.border = element_blank())
    gg <- gg + theme(axis.ticks = element_blank())
    gg <- gg + theme(plot.title = element_text(face = "bold", size = titleSize))
    gg <- gg + theme(plot.title = element_text(hjust = 0.5))
    gg <- gg + theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.y = element_text(size = labelSize),
        axis.title.y = element_blank(),
        panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
    return(gg)
}


#' Create architecure plot for a pair of seed and ortholog protein
#' @usage pairDomainPlotting(seed, ortho, seedDf, orthoDf, minStart, maxEnd,
#'     labelSize, titleSize, showScore, domainNamePos, firstDist,  nameType,
#'     nameSize, nameColor, labelPos, colorPallete)
#' @param seed Seed ID
#' @param ortho Ortho ID
#' @param seedDf domain dataframe for seed domains containing the seed ID,
#' ortholog ID, sequence length, feature names, start and end positions,
#' feature weights (optional) and the status to determine if that feature is
#' important for comparison the architecture between 2 proteins* (e.g. seed
#' protein vs ortholog) (optional)
#' @param orthoDf domain dataframe for ortholog domains (same format as seedDf)
#' @param minStart the smallest start position of all domains
#' @param maxEnd the highest stop position of all domains
#' @param labelSize lable size. Default = 12
#' @param titleSize title size. Default = 12
#' @param showScore show/hide E-values and Bit-scores. Default = NULL (hide)
#' @param domainNamePos list of positions for domain names, choose from "plot",
#' "legend" or "axis". Default: "plot"
#' @param firstDist distance of the first domain to plot title. Default = 0.5
#' @param nameType type of domain names, either "Texts" or "Labels" (default)
#' @param nameSize size of domain names (for Texts only). Default = 5
#' @param nameColor color of domain names (for Texts only). Default = "black"
#' @param labelPos position of doman names (for Labels only). Choose from
#' "Above" (default), "Below" or "Inside" the domain bar
#' @param colorPallete color pallete. Default = Paired"
#' @return Domain plot of a pair proteins as a arrangeGrob object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' \dontrun{
#' seed <- "101621at6656"
#' ortho <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' ortho <- gsub("\\|", ":", ortho)
#' grepID <- paste(seed, "#", ortho, sep = "")
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seed, domainFile, "file")
#' subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
#' subdomainDf$feature <- as.character(subdomainDf$feature)
#' orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
#' seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
#' minStart <- min(subdomainDf$start)
#' maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
#' g <- pairDomainPlotting(seed,ortho,seedDf,orthoDf,minStart,maxEnd,9,9)
#' grid::grid.draw(g)
#' }

pairDomainPlotting <- function(
        seed = NULL, ortho = NULL, seedDf = NULL, orthoDf = NULL,
        minStart = 0, maxEnd = 999, labelSize = 12, titleSize = 12,
        showScore = NULL, domainNamePos = "plot", firstDist = 0.5,
        nameType = "Labels", nameSize = 5, nameColor = "#000000",
        labelPos = "Above", colorPallete = "Paired"
) {
    if(is.null(seed) | is.null(ortho) | is.null(seedDf) | is.null(orthoDf))
        stop("Seed/Ortho ID or domain dataframe is NULL!")

    sep <- "|"
    plotOrtho <- singleDomainPlotting(
        orthoDf, ortho, sep, labelSize, titleSize, minStart, maxEnd, colorPallete,
        showScore, domainNamePos, firstDist, nameType, nameSize, nameColor, labelPos)
    plotSeed <- singleDomainPlotting(
        seedDf, seed, sep, labelSize, titleSize, minStart, maxEnd, colorPallete,
        showScore, domainNamePos, firstDist, nameType, nameSize, nameColor, labelPos)
    if (ortho == seed) {
        g <- plotSeed
    } else {
        seedHeight <- length(levels(as.factor(seedDf$feature)))
        orthoHeight <- length(levels(as.factor(orthoDf$feature)))
        if ("legend" %in% domainNamePos) {
            g <- gridArrangeSharedLegend(
                plotSeed, plotOrtho,
                ncol = 1, nrow = 2, position = "bottom"
            )
        } else {
            g <- gridExtra::arrangeGrob(plotSeed, plotOrtho, ncol = 1, nrow = 2)
        }
    }
    return(g)
}

#' Create architecure plot for a pair of seed and ortholog protein
#' @usage pairDomainPlotting2(seed, ortho, seedDf, orthoDf, minStart, maxEnd,
#'     labelSize, titleSize)
#' @param seed Seed ID
#' @param ortho Ortho ID
#' @param seedDf domain dataframe for seed domains containing the seed ID,
#' ortholog ID, sequence length, feature names, start and end positions,
#' feature weights (optional) and the status to determine if that feature is
#' important for comparison the architecture between 2 proteins* (e.g. seed
#' protein vs ortholog) (optional).
#' @param orthoDf domain dataframe for ortholog domains (same format as seedDf).
#' @param minStart the smallest start position of all domains
#' @param maxEnd the highest stop position of all domains
#' @param labelSize lable size. Default = 12.
#' @param titleSize title size. Default = 12.
#' @return Domain plot of a pair proteins as a arrangeGrob object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' \dontrun{
#' seed <- "101621at6656"
#' ortho <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' ortho <- gsub("\\|", ":", ortho)
#' grepID <- paste(seed, "#", ortho, sep = "")
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seed, domainFile, "file")
#' subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
#' subdomainDf$feature <- as.character(subdomainDf$feature)
#' orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
#' seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
#' minStart <- min(subdomainDf$start)
#' maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
#' g <- pairDomainPlotting(seed,ortho,seedDf,orthoDf,minStart,maxEnd,9,9)
#' grid::grid.draw(g)
#' }

pairDomainPlotting2 <- function(
    seed = NULL, ortho = NULL, seedDf = NULL, orthoDf = NULL,
    minStart = 0, maxEnd = 999, labelSize = 12, titleSize = 12
) {
    if(is.null(seed) | is.null(ortho) | is.null(seedDf) | is.null(orthoDf))
        stop("Seed/Ortho ID or domain dataframe is NULL!")
    # create color scheme, so that the same features in seed & ortholog will
    # have the same colors
    featureSeed <- levels(as.factor(seedDf$feature))
    featureOrtho <- levels(as.factor(orthoDf$feature))
    allFeatures <- c(featureSeed, featureOrtho)
    allColors <- getQualColForVector(allFeatures)
    colorScheme <- structure(allColors, .Names = allFeatures)
    # plot
    sep <- "|"
    plotOrtho <- singleDomainPlotting(
        orthoDf, ortho, sep, labelSize, titleSize, minStart, maxEnd,colorScheme)
    plotSeed <- singleDomainPlotting(
        seedDf, seed, sep, labelSize, titleSize, minStart, maxEnd, colorScheme)
    if (ortho == seed) {
        g <- gridExtra::arrangeGrob(plotSeed, ncol = 1)
    } else {
        seedHeight <- length(levels(as.factor(seedDf$feature)))
        orthoHeight <- length(levels(as.factor(orthoDf$feature)))
        g <- gridExtra::arrangeGrob(
            plotSeed, plotOrtho, ncol = 1, heights = c(seedHeight, orthoHeight)
        )
    }
    return(g)
}

#' Sort one domain dataframe based on the other domain dataframe
#' @description Sort domain dataframe of one protein (either seed or ortholog)
#' based on the dataframe of the its paired protein, in order to bring the
#' common domain feature in the same order which make it easy for comparing.
#' @param seedDf data of seed protein
#' @param orthoDf data of ortholog protein
#' @return Dataframe contains sorted domain list.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' \dontrun{
#' # get domain data
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' # get seedDf and orthoDf
#' subDf <- domainDf[
#'     domainDf$seedID ==
#'     "101621at6656#101621at6656:AGRPL@224129@0:224129_0:001955:1",]
#' orthoDf <- subDf[subDf$orthoID == "101621at6656:DROME@7227@1:Q9VG04",]
#' seedDf <- subDf[subDf$orthoID != "101621at6656:DROME@7227@1:Q9VG04",]
#' # sort
#' PhyloProfile:::sortDomains(seedDf, orthoDf)
#' }

sortDomains <- function(seedDf, orthoDf){
    if (is.null(seedDf) | is.null(orthoDf))
        stop("Domain data for seed & ortholog cannot be NULL!")
    orderNo <- NULL
    # get list of features in seedDf
    featureList <- as.data.frame(levels(as.factor(seedDf$feature)))
    colnames(featureList) <- c("feature")
    # and add order number to each feature
    featureList$orderNo <- seq(length(featureList$feature))

    # merge those info to orthoDf
    orderedOrthoDf <- merge(orthoDf, featureList, all.x = TRUE)

    # sort orthoDf
    index <- with(orderedOrthoDf, order(orderNo))
    orderedOrthoDf <- orderedOrthoDf[index, ]

    #turn feature column into a character vector
    orderedOrthoDf$feature <- as.character(orderedOrthoDf$feature)
    #then turn it back into an ordered factor (to keep this order when plotting)
    orderedOrthoDf$feature <- factor(
        orderedOrthoDf$feature, levels = unique(orderedOrthoDf$feature)
    )
    #return sorted df
    return(orderedOrthoDf)
}

#' Sort one domain dataframe based on list of ordered feature types
#' @description Sort domain dataframe of one protein based on a given list of
#' ordered feature types
#' @param domainDf domain dataframe
#' @param featureTypeOrder vector of ordered feature types
#' @return Dataframe contains sorted domain list.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom dplyr left_join
#' @examples
#' \dontrun{
#' NEED TO ADD AN EXAMPLE HERE
#' }

sortDomainsByList <- function(domainDf = NULL, featureTypeOrder = NULL) {
    if (is.null(domainDf) | is.null(featureTypeOrder))
        stop("Domain data or feature type order is NULL!")
    featureTypeOrder <- rev(featureTypeOrder[
        featureTypeOrder %in% levels(as.factor(domainDf$feature_type))
    ])
    orderedDomainDf <- dplyr::left_join(
        data.frame(feature_type = featureTypeOrder),
        domainDf, by = "feature_type"
    )
    orderedDomainDf$feature <- factor(
        orderedDomainDf$feature, levels = unique(orderedDomainDf$feature)
    )
    return(orderedDomainDf)
}

#' Modify feature names
#' @description Simplify feature names (e.g. TM for transmembrane domain,
#' LCR for low complexity regions, remove tool names from domain name) and add
#' weight to feature names (if available)
#' @param domainDf domain data as a dataframe object
#' @return Dataframe contains simlified domain names in yLabel column
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' \dontrun{
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' seedID <- "101621at6656"
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' PhyloProfile:::modifyFeatureName(domainDf)
#' }

modifyFeatureName <- function(domainDf = NULL) {
    if (is.null(domainDf)) stop("Domain data cannot be NULL!")
    domainDf$yLabel <- domainDf$feature_id
    domainDf$yLabel[domainDf$yLabel == "transmembrane"] <- "TM"
    domainDf$yLabel[domainDf$yLabel == "low complexity regions"] <- "LCR"
    domainDf$yLabel[domainDf$yLabel == "low_complexity_regions"] <- "LCR"
    return(domainDf)
}


#' Join multiple plots and merge legends
#' @description See: https://github.com/tidyverse/ggplot2/wiki
#' /Share-a-legend-between-two-ggplot2-graphs
#' @param ... List of plots, seperated by comma
#' @param ncol number of columns for final plot
#' @param nrow number of row for final plot
#' @param position position of legend (bottom or right)

gridArrangeSharedLegend <- function (
        ..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")
) {
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(
        position,
        "bottom" = gridExtra::arrangeGrob(
            do.call(arrangeGrob, gl),
            legend,
            ncol = 1,
            heights = grid::unit.c(unit(1, "npc") - lheight, lheight)
        ),
        "right" = gridExtra::arrangeGrob(
            do.call(arrangeGrob, gl),
            legend,
            ncol = 2,
            widths = grid::unit.c(unit(1, "npc") - lwidth, lwidth)
        )
    )

    return(combined)
}

#' Identify feature type(s) containing overlapped domains/features
#' @param domainDf input domain dataframe
#' @return List of feature types that have overlapped domains
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom dplyr group_by add_count
#' @examples
#' \dontrun{
#' NEED TO ADD AN EXAMPLE HERE
#' }

checkOverlapDomains <- function(domainDf) {
    if (is.null(domainDf)) stop("Domain data cannot be NULL!")
    feature_type <- NULL
    df <- domainDf[, c("feature", "feature_type", "start", "end")]
    df <- data.frame(df %>% dplyr::group_by(feature_type) %>% dplyr::add_count(feature_type))
    df$tmp <- paste(df$feature_type, df$end, sep = "_")
    overlappedType <- lapply(
        df$tmp[df$n > 1],
        function (x) {
            ed <- df$end[df$tmp == x][1]
            st <- df$start[df$tmp == x][1]
            type <- df$feature_type[df$tmp == x][1]

            subDf <- df[
                !(df$tmp == x) & df$feature_type == type & df$start < ed & df$end > st,
            ]
            if (nrow(subDf) > 0) {
                if(!(subDf$feature[1] == levels(as.factor(df$feature[df$tmp == x]))))
                    return(subDf$feature_type)
            }
        }
    )
    return(unique(unlist(overlappedType)))
}

#' Modify domain dataframe to resolve overlapped domains/features
#' @param domainDf input domain dataframe
#' @return Domain dataframe with modified feature names that join multiple
#' domains of the same type that are not overlapped
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' \dontrun{
#' NEED TO ADD AN EXAMPLE HERE
#' }

resolveOverlapFeatures <- function(domainDf) {
    if (is.null(domainDf)) stop("Domain data cannot be NULL!")
    overlappedType <- checkOverlapDomains(domainDf)
    domainDf$featureOri <- domainDf$feature
    domainDf$featureOri <- as.character(domainDf$featureOri)
    if (length(overlappedType) > 0) {
        domainDf$feature <- ifelse(
            !(domainDf$feature_type %in% overlappedType),
            domainDf$feature_type,
            domainDf$feature
        )
    } else {
        domainDf$feature <- domainDf$feature_type
    }
    return(domainDf)
}


#' Add colors for each feature/domain
#' @description Add colors to features/domains of 2 domain dataframes. Users can
#' choose to color only the shared features, unique features, all features
#' (default) or based on feature types. Default color pallete is "Paired", but
#' it can be changed.
#' @param seedDf Domain dataframe of seed protein (protein 1)
#' @param orthoDf Domain dataframe of orthologs protein (protein 2)
#' @param colorType Choose to color "all", "shared", "unique" features or color
#' by "Feature type". Default: "all"
#' @param colorPallete Choose between "Paired", "Set1", "Set2", "Set3",
#' "Accent", "Dark2" for the color pallete
#' @param ignoreInstanceNo Ignore number of feature instances while identifying
#' shared or unique features. Default: FALSE
#' @return 2 dataframes (seedDf and orthoDf) with an additional column for the
#' assigned colors to each feature instance
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom utils head
#' @examples
#' \dontrun{
#' NEED TO ADD AN EXAMPLE HERE
#' }

addFeatureColors <- function(
        seedDf = NULL, orthoDf = NULL, colorType = "all",
        colorPallete = "Paired", ignoreInstanceNo = FALSE
) {
    if (is.null(seedDf) | is.null(orthoDf)) stop("Domain Df cannot be null!")
    feature <- NULL
    featureSeedCount <- seedDf %>% dplyr::count(feature)
    featureOrthoCount <- orthoDf %>% dplyr::count(feature)
    featureCount <- merge(featureSeedCount, featureOrthoCount, by = "feature", all = TRUE)
    featureCount$type <- "unique"
    featureCount$type[featureCount$n.x == featureCount$n.y] <- "shared"
    if (ignoreInstanceNo == TRUE)
        featureCount$type[!is.na(featureCount$n.x) & !is.na(featureCount$n.y)] <- "shared"

    featureCount$feature <- as.character(featureCount$feature)
    sharedFeatures <- unique(featureCount$feature[featureCount$type == "shared"])
    uniqueFeatures <- unique(featureCount$feature[featureCount$type == "unique"])
    allFeatures <- c(sharedFeatures, uniqueFeatures)

    if (colorType == "Unique" & length(uniqueFeatures) == 0)
        colorType = "All"
    if (colorType == "Shared" & length(sharedFeatures) == 0)
        colorType = "All"

    if (colorType == "Unique") {
        sharedFeaturesColors <- rep("#C9C9C9", length(sharedFeatures))
        uniqueFeaturesColors <- getQualColForVector(uniqueFeatures)
        if (checkColorPallete(uniqueFeatures, colorPallete) == TRUE) {
            uniqueFeaturesColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(length(uniqueFeatures), colorPallete),
                    length(uniqueFeatures)
                ))
        }
        allColors <- c(sharedFeaturesColors, uniqueFeaturesColors)
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else if (colorType == "Shared") {
        sharedFeaturesColors <- getQualColForVector(sharedFeatures)
        uniqueFeaturesColors <- rep("#C9C9C9", length(uniqueFeatures))
        if (checkColorPallete(sharedFeatures, colorPallete) == TRUE) {
            sharedFeaturesColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(length(sharedFeatures), colorPallete),
                    length(sharedFeatures)
                ))
        }
        allColors <- c(sharedFeaturesColors, uniqueFeaturesColors)
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else if (colorType == "All") {
        allColors <- getQualColForVector(allFeatures)
        if (checkColorPallete(allFeatures, colorPallete) == TRUE) {
            allColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(length(allFeatures), colorPallete),
                    length(allFeatures)
                ))
        }
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else {
        tmpDf <- data.frame(str_split_fixed(allFeatures, "_", 2))
        tmpDf$name <- paste(tmpDf$X1, tmpDf$X2, sep = "_")
        tmpDf$name[tmpDf$X2 == ""] <- tmpDf$X1[tmpDf$X2 == ""]
        tmpDf$X2[tmpDf$X2 == ""] <- tmpDf$X1[tmpDf$X2 == ""]
        typeColorDf <- data.frame(
            colors = head(
                suppressWarnings(RColorBrewer::brewer.pal(nlevels(as.factor(tmpDf$X1)), colorPallete)),
                nlevels(as.factor(tmpDf$X1))
            ),
            X1 = levels(as.factor(tmpDf$X1))
        )
        tmpDf <- merge(tmpDf, typeColorDf, all.x = TRUE)
        colorScheme <- data.frame(color = tmpDf$colors, feature = tmpDf$name)
    }

    # add color to seedDf and orthoDf
    seedDf <- merge(seedDf, colorScheme, by = "feature", all.x = TRUE)
    orthoDf <- merge(orthoDf, colorScheme, by = "feature", all.x = TRUE)

    return(list(seedDf, orthoDf))
}
