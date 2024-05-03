#' Protein domain architecture plot
#' @param pointInfo() info of clicked point
#' (from reactive fn "pointInfoDetail")
#' @param domainInfo() domain information
#' (from reactive fn "getDomainInformation")
#' @param labelArchiSize lable size (from input$labelArchiSize)
#' @param titleArchiSize title size (from input$titleArchiSize)
#' @param archiHeight plot height (from input$archiHeight)
#' @param archiWidth plot width (from input$archiWidth)
#' @param seqIdFormat sequence ID format (either bionf or unknown)
#' @param currentNCBIinfo dataframe of the pre-processed NCBI taxonomy data
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

createArchitecturePlotUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                3,
                radioButtons(
                    ns("resolveOverlap"),
                    "Merge non-overlapped features", 
                    choices = c("Yes","No"), selected = "Yes", 
                    inline = TRUE
                ),
                checkboxGroupInput(
                    ns("showName"),
                    "Display feature names",
                    choices = c(
                        "On the plot" = "plot",
                        "As a legend" = "legend",
                        "On the y-axis" = "axis"
                    ),
                    selected = c("plot","axis")
                )
            ),
            column(
                3,
                selectInput(
                    ns("feature"),
                    "Exclude features",
                    choices = c(
                        "flps","seg","coils","signalp","tmhmm",
                        "smart","pfam",
                        "without E-value" = "noEvalue",
                        "without Bit-score" = "noBitscore"
                    ),
                    multiple = TRUE
                ),
                bsButton(ns("featureOpt"), "Other feature options"),
                checkboxInput(ns("plotConfig"), "Plot configuration", value = FALSE)
            ),
            column(
                3,
                checkboxGroupInput(
                    ns("showScore"),
                    "Show information",
                    choices = c(
                        "E-value", "Bit-score"
                    )
                ),
                uiOutput(ns("filterEvalue.ui")),
                uiOutput(ns("filterBitscore.ui"))
            ),
            column(
                3,
                checkboxGroupInput(
                    ns("showInstance"),
                    "Show only instances with",
                    choices = c(
                        "Best E-value" = "evalue", 
                        "Best Bit-score" = "bitscore",
                        "Paths" = "path"
                    )
                )
            )
        ),
        br(),
        fluidRow(
            conditionalPanel(
                condition = {sprintf("input['%s'] == 1", ns("plotConfig"))},
                column(
                    3,
                    createPlotSize(ns("archiHeight"), "Plot height(px)",400, 200),
                    createPlotSize(ns("archiWidth"), "Plot width (px)", 800, 200)
                ),
                column(
                    3,
                    createTextSize(
                        ns("titleArchiSize"), "Title/Seq ID size (px)", 14, 200
                    ),
                    createTextSize(
                        ns("labelArchiSize"), "Axis label size(px)", 12, 200
                    )
                ),
                column(
                    6,
                    createTextSize(
                        ns("nameSize"), "Feature segment size (mm)", 5, 200
                    ),
                    sliderInput(
                        ns("firstDist"), "Distance between plot title and the 1st feature", 
                        min = 0, max = 5, value = 0.5, step = 0.1, width = 400
                    )
                ),
                column(
                    6,
                    radioButtons(
                        ns("colorType"),"Color feature instances", inline = TRUE,
                        choices = c("Shared","Unique","All","Feature type"), selected = "All"
                    ),
                    checkboxInput(
                        ns("ignoreInstanceNo"), "Ignore number of instances", value = FALSE
                    ),
                ),
                column(
                    6,
                    selectInput(
                        ns("colorPallete"),
                        "Color pallete",
                        choices = c("Paired", "Set1", "Set2", "Set3", "Accent", "Dark2"),
                        selected = "Paired"
                    )
                )
            )
        ),
        hr(),
        uiOutput(ns("archiPlot.ui")),
        verbatimTextOutput(ns("hover_info")),
        br(),
        downloadButton(ns("archiDownload"), "Download plot", class = "butDL"),
        hr(),
        tableOutput(ns("linkTable")),
        checkboxInput(
            ns("showDomainTable"), "Show detailed feature table", value = FALSE
        ),
        conditionalPanel(
            condition = {sprintf("input['%s'] == 1", ns("showDomainTable"))},
            DT::dataTableOutput(ns("domainTable"))
        ),
        
        bsModal(
            ns("featureConfigBs"),
            "Feature name configuration",
            ns("featureOpt"),
            size = "small",
            br(),
            radioButtons(
                ns("nameType"),"Type of feature names", inline = TRUE,
                choices = c("Labels","Texts"), selected = "Labels"
            ),
            conditionalPanel(
                condition = {sprintf("input['%s'] == 'Labels'", ns("nameType"))},
                radioButtons(
                    ns("labelPos"),"Label position", inline = TRUE,
                    choices = c("Above","Inside","Below"), selected = "Above"
                )
            ),
            conditionalPanel(
                condition = {sprintf("input['%s'] == 'Texts'", ns("nameType"))},
                colourpicker::colourInput(
                    ns("nameColor"),
                    "Feature name color",
                    value = "#000000"
                )
            ),
            selectInput(
                ns("excludeNames"),
                "Exclude feature names of",
                choices = c(
                    "flps","seg","coils","signalp","tmhmm",
                    "smart","pfam"
                ),
                selected = c("tmhmm","signalp","seg","flps","coils"),
                multiple = TRUE
            ),
            radioButtons(
                ns("featureTypeSort"),"Sort feature types by shared features", inline = TRUE,
                choices = c("Yes","No"), selected = "Yes"
            ),
            conditionalPanel(
                condition={sprintf("input['%s'] == 'No'",ns("featureTypeSort"))},
                selectInput(
                    ns("featureTypeOrder"),
                    "Feature type order",
                    choices = c(
                        "pfam", "smart", "tmhmm", "coils", "signalp", "seg", "flps"
                    ),
                    selected = c(
                        "pfam", "smart", "tmhmm", "coils", "signalp", "seg", "flps"
                    ),
                    multiple = TRUE
                )
            )
        )
        
        # column(
        #     4,
        #     style = "padding:0px;",
        #     selectInput(
        #         ns("showFeature"),
        #         "Show",
        #         choices = c(
        #             "All features" = "all",
        #             "Shared features" = "common",
        #             "Unique features" = "unique"
        #         ),
        #         selected = "all"
        #     )
        # ),
        # column(
        #     4,
        #     selectizeInput(
        #         ns("excludeFeature"),
        #         "Exclude feature type(s)",
        #         choices = c(
        #             "flps","seg","coils","signalp","tmhmm","smart","pfam"
        #         ),
        #         multiple = TRUE, options=list(placeholder = 'None')
        #     )
        # ),
        # column(4, uiOutput(ns("featureList.ui"))),
        # column(12, uiOutput(ns("archiPlot.ui"))),
        # downloadButton(ns("archiDownload"), "Download plot", class = "butDL"),
        # tags$head(
        #     tags$style(HTML(
        #         ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        # ),
        # br(),
        # br(),
        # h4(strong("LINKS TO ONLINE DATABASE")),
        # textOutput(ns("selectedDomain")),
        # tableOutput(ns("domainTable")),
        # HTML(
        #     paste0(
        #         "<p><em><strong>Disclaimer:</strong> ",
        #         "External links are automatically generated and may point to ",
        #         "a wrong target (see <a ",
        #         "href=\"https://github.com/BIONF/PhyloProfile/wiki/FAQ",
        #         "#wrong-info-from-public-databases\" ",
        #         "target=\"_blank\">FAQ</a>)</em></p>"
        #     )
        # )
    )
}

createArchitecturePlot <- function(
    input, output, session, pointInfo, domainInfo, 
    # labelArchiSize, 
    # titleArchiSize, archiHeight, archiWidth, 
    seqIdFormat, currentNCBIinfo
){
    # update excludeNames if no feature type on the y-axis =====================
    observe({
        req(input$showName)
        if (!("axis" %in% input$showName) & !("legend" %in% input$showName)) {
            updateSelectInput(
                session, "excludeNames",
                "Exclude feature names of",
                choices = c(
                    "flps","seg","coils","signalp","tmhmm",
                    "smart","pfam"
                )
            )
        } else if ("axis" %in% input$showName | "legend" %in% input$showName) {
            updateSelectInput(
                session, "excludeNames",
                "Exclude feature names of",
                choices = c(
                    "flps","seg","coils","signalp","tmhmm",
                    "smart","pfam"
                ),
                selected = c("tmhmm","signalp","flps","seg","coils")
            )
        }
    })
    
    # * render e-value / bitscore filter =======================================
    output$filterEvalue.ui <- renderUI({
        req(domainInfo())
        df <- domainInfo()
        maxEvalue = 1
        if ("evalue" %in% colnames(df))
            maxEvalue <- format(max(df$evalue[!is.na(df$evalue)]), scientific = TRUE, digits = 2)
        if ("E-value" %in% input$showScore) {
            numericInput(
                "minEvalue", "Filter E-value:",
                min = 0,
                max = maxEvalue,
                value = format(0.00001, scientific = TRUE, digits = 2)
            )
        }
    })
    
    output$filterBitscore.ui <- renderUI({
        req(domainInfo())
        df <- domainInfo()
        if ("Bit-score" %in% input$showScore) {
            numericInput(
                "minBitscore", "Filter Bit-score:",
                min = min(df$bitscore[!is.na(df$bitscore)]),
                max = 9999,
                value = min(df$bitscore[!is.na(df$bitscore)])
            )
        }
    })
    
    # * filter domain features -------------------------------------------------
    filterDomainDf <- reactive({
        outDf <- domainInfo()
        if (is.null(nrow(outDf))) stop("Domain info is NULL!")
        
        # # filter domain df by feature type
        # if (!("feature_type" %in% colnames(df))) {
        #     df[c("feature_type","feature_id")] <- 
        #         stringr::str_split_fixed(df$feature, '_', 2)
        #     df$feature_id[df$feature_type == "smart"] <-
        #         paste0(df$feature_id[df$feature_type == "smart"], "_smart")
        # }
        # df <- df[!(df$feature_type %in% input$excludeFeature),]
        # return(outDf)
        
        # filter domain df by features
        if (nrow(outDf) == 0) stop("Domain info is NULL!")
        
        outDf[c("feature_type","feature_id")] <- stringr::str_split_fixed(outDf$feature, '_', 2)
        outDf <- outDf[!(outDf$feature_type %in% input$feature),]
        # filter filters without e-value and/or bitscore
        if ("evalue" %in% colnames(outDf)) {
            if ("noEvalue" %in% input$feature)
                outDf <- outDf[!is.na(outDf$evalue),]
            if ("noBitscore" %in% input$feature)
                outDf <- outDf[!is.na(outDf$bitscore),]
        }
        # modify feature IDs
        outDf$feature_id_mod <- outDf$feature_id
        outDf$feature_id_mod <- gsub("SINGLE", "LCR", outDf$feature_id_mod)
        outDf$feature_id_mod[outDf$feature_type == "coils"] <- "Coils"
        outDf$feature_id_mod[outDf$feature_type == "seg"] <- "LCR"
        outDf$feature_id_mod[outDf$feature_type == "tmhmm"] <- "TM"
        # exclude features IDs
        if (!is.null(input$excludeNames)) {
            outDf$feature_id_mod[outDf$feature_type %in% input$excludeNames] <- NA
        }
        
        # enable/disable option for showing evalue/bitscore
        if ("evalue" %in% colnames(outDf)) {
            shinyjs::enable("showScore")
        } else {
            shinyjs::disable("showScore")
        }
        
        # Filter data by e-value, bit-score and feature path
        if ("evalue" %in% colnames(outDf)) {
            # filter by e-value and/or bit-score
            if ("E-value" %in% input$showScore) {
                # req(input$minEvalue)
                minEvalue <- format(input$minEvalue, scientific = FALSE)
                naOutDf <- outDf[is.na(outDf$evalue),]
                outDf <- outDf[!is.na(outDf$evalue) & outDf$evalue <= input$minEvalue,]
                outDf <- rbind(outDf,naOutDf)
            }   
            if ("Bit-score" %in% input$showScore) {
                # req(input$minBitscore)
                naOutDf <- outDf[is.na(outDf$bitscore),]
                outDf <- outDf[!is.na(outDf$bitscore) & outDf$bitscore >= input$minBitscore,]
                outDf <- rbind(outDf,naOutDf)
            }
            # get only best instances
            if ("evalue" %in% input$showInstance) {
                naOutDf <- outDf[is.na(outDf$evalue),]
                outDf <- outDf %>% group_by(feature, orthoID) %>% filter(evalue == min(evalue))
                outDf <- rbind(outDf,naOutDf)
            }
            if ("bitscore" %in% input$showInstance) {
                naOutDf <- outDf[is.na(outDf$bitscore),]
                outDf <- outDf %>% group_by(feature, orthoID) %>% filter(bitscore == max(bitscore))
                outDf <- rbind(outDf,naOutDf)
            }
            if ("path" %in% input$showInstance) {
                outDf <- outDf %>% group_by(feature) %>% filter(path == "Y")
            }
            # Format e-values
            outDf$evalue[!is.na(outDf$evalue)] <- 
                format(outDf$evalue[!is.na(outDf$evalue)], scientific = TRUE, digits = 2)
        }
        return(outDf[!is.na(outDf$seedID),])
    })
    
    getSeqIdFormat <- reactive({
        if (seqIdFormat() == 1) return("bionf")
        return("unknown")
    })
    
    # * render plot ------------------------------------------------------------
    output$archiPlot <- renderPlot({
        if (is.null(nrow(filterDomainDf()))) stop("Domain info is NULL!")
        # remove user specified features (from input$featureList)
        df <- filterDomainDf()
        df <- df[!(df$feature_id %in% input$featureList),]
        # generate plot
        g <- createArchiPlot(
            pointInfo(), df, input$labelArchiSize, input$titleArchiSize,
            "all", getSeqIdFormat(), currentNCBIinfo()
        )
        if (any(g == "No domain info available!")) {
            msgPlot()
        } else {
            grid.draw(g)
        }
    })

    output$archiPlot.ui <- renderUI({
        ns <- session$ns
        if (is.null(nrow(filterDomainDf()))) {
            msg <- paste0(
                "<p><em>Wrong domain file has been uploaded!
        Please check the correct format in
        <a href=\"https://github.com/BIONF/PhyloProfile/wiki/",
                "Input-Data#ortholog-annotations-eg-domains\"
        target=\"_blank\" rel=\"noopener\">our PhyloProfile wiki</a>.</em></p>"
            )
            HTML(msg)
        } else {
            shinycssloaders::withSpinner(
                plotOutput(
                    ns("archiPlot"),
                    height = input$archiHeight,
                    width = input$archiWidth,
                    click = ns("archiClick")
                )
            )
        }
    })

    output$archiDownload <- downloadHandler(
        filename = function() {
            c("domains.pdf")
        },
        content = function(file) {
            # remove user specified features (from input$featureList)
            df <- filterDomainDf()
            df <- df[!(df$feature_id %in% input$featureList),]
            # generate plot
            g <- createArchiPlot(
                pointInfo(), filterDomainDf(), input$labelArchiSize, input$titleArchiSize,
                input$showFeature, getSeqIdFormat(), currentNCBIinfo()
            )
            grid.draw(g)
            # save plot to file
            ggsave(
                file, plot = g,
                width = input$archiWidth * 0.056458333,
                height = input$archiHeight * 0.056458333,
                units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
            )
        }
    )

    # output$selectedDomain <- renderText({
    #     if (is.null(input$archiClick$y)) return("No domain selected!")
    #     y <- input$archiClick$y
    #     # paste(y, round(y), convertY(unit(y, "npc"), "px"))
    #     
    # })
    
    output$featureList.ui <- renderUI({
        ns <- session$ns
        allFeats <- getAllFeatures(pointInfo(), filterDomainDf())
        selectizeInput(
            ns("featureList"), "Exclude individual feature(s)", multiple = TRUE,
            choices = allFeats,
            options=list(placeholder = 'None')
        )
    })
    
    output$linkTable <- renderTable({
        if (is.null(nrow(filterDomainDf()))) return("No domain info available!")
        features <- getDomainLink(pointInfo(), filterDomainDf(), getSeqIdFormat())
        features
    }, sanitize.text.function = function(x) x)
}

#' plot error message
#' @return error message in a ggplot object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
msgPlot <- function() {
    msg <- paste(
        "No information about domain architecture!",
        "Please check:","if you uploaded the correct domain file/folder; or ",
        "if the selected genes (seed & ortholog) do exist in the uploaded file",
        "(please search for the corresponding seedID and hitID)",
        sep = "\n"
    )
    x <- c(1,2,3,4,5)
    y <- c(1,2,3,4,5)
    g <- ggplot(data.frame(x, y), aes(x,y)) +
        geom_point(color = "white") +
        annotate(
            "text", label = msg, x = 3.5, y = 0.5, size = 5, colour = "red"
        ) +
        theme(axis.line = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_blank()) +
        ylim(0,1)
    return(g)
}

#' Get list of all features
#' @return A list of all features
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getAllFeatures <- function(info, domainDf) {
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
    orthoID <- feature <- NULL
    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
    else {
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        feature <- c(
            levels(as.factor(orthoDf$feature_id)), 
            levels(as.factor(seedDf$feature_id))
        )
    }
    return(levels(as.factor(feature)))
}

#' get pfam and smart domain info (domain name, acc, profile HMM,...)
#' @return dataframe for each type of database
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getDomainInfo <- function(info, domainDf, type) {
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    domainDf <- domainDf[grep(grepID, domainDf$seedID),]
    domainDf <- domainDf[domainDf$feature_type == type,]
    domainInfoDf <- data.frame()
    if (nrow(domainDf) > 0) {
        # get all features of for this pair proteins
        feature <- getAllFeatures(info, domainDf)
        # filter domain info
        if ("acc" %in% colnames(domainDf)) {
            domainInfoDf <- domainDf[
                , c("orthoID", "feature_id", "acc", "evalue", "bitscore")
            ]
            domainInfoDf$evalue <- format(domainInfoDf$evalue, scientific=TRUE)
        } else {
            domainInfoDf <- domainDf[, c("orthoID", "feature_id")]
            domainInfoDf <- domainInfoDf[!(duplicated(domainInfoDf)),]
            domainInfoDf$acc <- rep(NA, nrow(domainInfoDf))
            domainInfoDf$evalue <- rep(NA, nrow(domainInfoDf))
            domainInfoDf$bitscore <- rep(NA, nrow(domainInfoDf))
        }
    }
    return(domainInfoDf)
}

#' get pfam and smart domain links
#' @return dataframe with domain IDs and their database links
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getDomainLink <- function(info, domainDf, seqIdFormat) {
    featurePfam <- getDomainInfo(info, domainDf, "pfam")
    pfamDf <- createLinkTable(featurePfam, "pfam")
    featureSmart <- getDomainInfo(info, domainDf, "smart")
    smartDf <- createLinkTable(featureSmart, "smart")
    
    featDf <- rbind(pfamDf, smartDf)
    featDf <- subset(featDf, select=c(orthoID,feature_id,link,evalue,bitscore))
    featDf <-featDf[order(featDf$orthoID),]
    if (seqIdFormat == "bionf") {
        featDf[c("groupID", "spec", "geneID", "misc")] <- 
            stringr::str_split_fixed(featDf$orthoID, ":", 4)
        featDf <- subset(featDf,select=c(geneID,feature_id,link,evalue,bitscore))
    }
    colnames(featDf) <- c("Gene ID", "Domain ID", "URL", "E-value", "Bit-score")
    return(featDf)
}

#' plot error message
#' @param featureDf dataframe contains 2 columns feature_id and acc
#' @param featureType pfam or smart
#' @return dataframe contains 2 columns feature_id and link
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
createLinkTable <- function(featureDf, featureType) {
    featureDf <- featureDf[!(duplicated(featureDf)),]
    if (nrow(featureDf) > 0) {
        if (featureType == "pfam") {
            featureDf$link[is.na(featureDf$acc)] <- paste0(
                "<a href='http://pfam-legacy.xfam.org/family/", featureDf$feature_id,
                "' target='_blank'>", "PFAM", "</a>"
            )
            featureDf$link[!(is.na(featureDf$acc))] <- paste0(
                "<a href='https://www.ebi.ac.uk/interpro/entry/pfam/", featureDf$acc, 
                "' target='_blank'>", "INTERPRO", "</a>"
            )
        } else {
            featureDf$link <- paste0(
                "<a href='http://smart.embl-heidelberg.de/smart/",
                "do_annotation.pl?BLAST=DUMMY&DOMAIN=",
                # featureDf$feature_id, "' target='_blank'>",
                gsub("_smart", "",featureDf$feature_id), "' target='_blank'>",
                "SMART", "</a>"
            )
        }
    }
    return(featureDf)
}
 

