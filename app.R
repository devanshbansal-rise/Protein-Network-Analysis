library(shiny)
library(shinycssloaders)
library(DT)
library(plotly)
library(visNetwork)
library(igraph)
library(data.table)
library(R.utils)
library(stringr)
library(corrplot)
library(ComplexHeatmap)
library(enrichR)
library(dplyr)
library(tidyr)
library(plotly)
library(tibble)
library(circlize)
library(nVennR)
degree <- igraph::degree
sd <- stats::sd
options(shiny.maxRequestSize = 30*1024^2)

ui <- fluidPage(
  tabsetPanel(
    tabPanel(title = "Input Data",
      titlePanel('Tool for Computational Protein Network Analysis'),
      'Welcome! Begin by choosing the data the tool will analyze. You may upload your own data - in which case you will have to fill some info about your table so it can be read - or choose to analyze one of our sample sets already curated for these purposes.',
      'The minimum requirements for uploading your own data are that the files are datatables/dataframes where each row has the same number of columns, every datatable is formatted the same way and uses the same ID mapping, and that format contains two columns of interacting proteins (so each row among those two columns represents an interaction). This is generally how binary interaction data is curated and how our tool generates the networks.',
      'Note that if your protein IDs are not already in gene name format, then ID mapping will be necessary to conduct enrichment (otherwise you will just get empty results). The other analyses can be conducted without mapping.',
      'Once you have made your choice, hit the Update button. The data will process and a table will show up with some basic statistics about each network.',
      'Then you can navigate to any of the other tabs to find more analyses!',
      sidebarLayout(
      sidebarPanel(
        selectInput('data','Data to Use', choices = list('Uploading My Own Data' = 1)),
                                                         # 'TissueNet Sample' = 2,
                                                         # 'PPI Database Sample' = 3)),
        conditionalPanel('input.data == 1', fileInput(inputId = "upload", label = "Upload files", multiple = TRUE)),
        conditionalPanel('input.data == 1', numericInput('column1','Which two columns of the datatable contain the data?',value=1, min=1)),
        conditionalPanel('input.data == 1', numericInput('column2',label='',value=2, min=2)),
        conditionalPanel('input.data == 1', checkboxInput('header', label = 'My datatable has a header')),
        conditionalPanel('input.data == 1', checkboxInput('convert', label = 'My data needs to be ID mapped to the respective gene names')),
        conditionalPanel("input.convert", selectInput('mapto', 'What ID are your proteins?', choices = list('UniProtKB-AC' = 1, 'UniProtKB-ID' = 2, 'GeneID (EntrezGene)' = 3, 'RefSeq' = 4, 'GI' = 5, 'PDB' = 6, 'GO' = 7, 'UniRef100' = 8, 'UniRef90' = 9, 'UniRef50' = 10, 'UniParc' = 11, 'PIR' = 12, 'NCBI-taxon' = 13, 'MIM' = 14, 'UniGene' = 15, 'PubMed' = 16, 'EMBL' = 17, 'EMBL-CDS' = 18, 'Ensembl' = 19, 'Ensembl_TRS' = 20, 'Ensembl_PRO' = 21, 'Additional PubMed' = 22))),
      actionButton('update','Update')
      ),
      mainPanel(withSpinner(tableOutput("info")))
    )
    ),
    tabPanel(title = 'Network Visualization and Clustering',
      'Here you can visualize each network in its entirety or by individual cluster. Below the cluster view, you can find analysis on each cluster individually - the cluster itself, its function determined by extracting all unique enrichments, and the interactions present. This is especially useful for very large networks to help break it down into smaller components.',
      'Clusters are determined by the Louvain method.',
      'If the enrichment result is empty, that means there were no unique enrichment results for that cluster, so it does not occupy a significant functional module.',
      'Since large networks can take a bit of time to view, that is an option you may toggle on your own. The remaining analyses will happen automatically, those may take a bit of processing too so please be patient.',
      uiOutput('choosegraph'),
      actionButton('view','View the full network in one graph'),
      withSpinner(plotOutput('fullnetwork')),
      'Below is the basic info on each cluster, followed by on option to zoom in one specific cluster and view its plot, unique function, and interactions.',
      dataTableOutput('clusterinfo'),
      uiOutput('choosecluster'),
      withSpinner(plotOutput('plotcluster')),
      plotOutput('clusterenrichment'),
      dataTableOutput('clusterints')
    ),
    tabPanel(title = 'Network Comparison',
      'There are various comparisons that can be made between the different databases, and various ways to visualize each of those.',
      titlePanel('Pair-wise Comparisons: Correlation Matrices'),
      'The first comparison below visualizes the intersection of each possible pair of databases in correlation matrices. The second two show those intersections as the number of nodes or edges shared, but that number itself may not provide as much value as a percentage of the database found in the other.',
      'Thus, the first two matrices convert the raw numeric data to percentages and plot that. The percentage in each cell should be read as how much of the database in the ROW is found in the database in the COLUMN.',
      withSpinner(plotOutput("correlation", height = 'auto')),
      titlePanel('Dendrograms'),
      'The dendrograms below help visualize the global distance between all the networks based on the percentages calculated above. One dendrogram is based on the node calculations, the other on the edges (usually they mirror each other).',
      plotOutput("dend"),
      titlePanel('Global Comparisons'),
      'Below you will find some more global comparisons of the networks.',
      'The bar graphs shows how much of each network is unique and universally shared among all networks (the tables just repeat that data to be read in its raw form).',
      'The Euler plots are just another way to visualize these overlaps. Note that these may take more time to process with a large number of datasets.',
      'The first set of charts are for intersections of nodes, the second for edges.',
      plotOutput('bar1'),
      dataTableOutput('uniqshared1'),
      htmlOutput('euler1'),
      plotOutput('bar2'),
      dataTableOutput('uniqshared2'),
      htmlOutput('euler2'),
      titlePanel('Network Specifity Measures'),
      'Below is a table listing each gene and how many network-specific interactions it is a part of.',
      dataTableOutput('tsintcount'),
      'The histograms below show the distribution of specificity of nodes and edges across the networks. They essentially answer the question of how network-specific/unique nodes and edges are throughout.',
      'More left-leaning histograms mean generally more unique, while right-leaning means generally more global.',
      plotOutput('histograms')
    ),
    tabPanel(title = 'Differential Connectivity',
      'This tab focuses on differential connectivity, the idea of how a gene\'s importance, determined by its degree (the number of genes it interacts with), varies across each network. If a gene has a low degree in one network and high in another, that could signal it as a gene/protein of interest.',
      'The differential connectivity of each gene found in all networks can be found below, ranked by the coefficient of variation (CV) of its degrees.',
      dataTableOutput('dctable'),
      'This table is visualized in a heatmap below. Each row is a gene, and the columns are its degree values in each network scaled by Z-score. Columns that stick out indicate networks with generally much lower (blue) or higher (red) degrees than the other networks.',
      withSpinner(plotOutput('dcheatmap')),
      'The next question is if a gene has high differential connectivity across networks, which network(s) is or are causing that variance in degree? If a network is causing high variance of a gene (i.e. the gene\'s degree changes heavily in that network), that gene can be associated with that network like a hub.',
      'This is visualized below in a Sankey plot mapping each network to its top \'hubs\' in this manner. This was determined by measuring the change in differential connectivity of a gene when the degree of that gene in each network was not counted toward the variance calculation.',
      'You may use the slider below to change how many genes show up for each network. Note that the number of actual genes that show up may be greater if there is a tie among genes (based on the score described before). You may also view that score for each association by hovering over the respective line.',
      sliderInput("cutoff", label = 'Top # of Associations Per Network', min = 1, max = 10, value = 3),
      withSpinner(plotlyOutput('sankey'))
    ),
    tabPanel(title = 'Isolate and Enrich Genes of Interest',
      'Here, via connection to EnrichR, you can isolate certain genes of interest by selecting from the drop down menu, and then view their annotations. At the moment, the annotations are from Gene Ontology\'s Biological Process database.',
      'The app is designed to conduct all enrichments first and then allow you to sort through the options, so it may take a little time at first to process all that at first.',
      'Another note is that hubs are calculated as the top 1% of genes based on either degree or closeness score. These are two very common measures of centrality and importance of a node in a network, and taking the top 1% helps fine-tune the most central proteins.',
      'The format of this tab is first the list of genes that are being isolated and enriched, then the resulting enrichment table where the results are ranked by p-value, then those results seen in a bar graph where the results are ranked by combined score and the top 20 or less are displayed.',
      selectInput('genes', 'Genes to Enrich', 
                  choices = list('Shared Interactions' = 1,
                                 'High DC' = 2,
                                 'Unique Proteins' = 3,
                                 'Hubs' = 4)),
      uiOutput('tochoosenet'),
      dataTableOutput('enrichedList'),
      dataTableOutput('enrichedDT'),
      plotOutput('enrichedPlot')
    ),
    tabPanel(title = 'Close-up on Select Proteins',
      'Here you can type in any protein found in any of the networks throughout your data and view a close-up on all its immediate interactions in your data. Hover over an edge line to view the network in which that interaction was found.',
      'The first toggle box allows you to view the neighbors of the neighbors as well. The second toggle allows you to view the associations the input gene has with each network (essentially switching the titles on each edge and the nodes).',
      'This plot is done through an interactive tool called visNetwork, which can be very slow for large networks. Thus, a warning will be given if the network that would be plotted is very large, and it is advised you only update the datatable in those cases and view the interactions in their raw form.',
      'This can be especially useful to further study genes of interest calculated from the previous two panel!',
      textInput('gene','Type in a gene'),
      checkboxInput('expand','Expand interactions of Neighbors'),
      checkboxInput('switch', 'View interactions by network rather than by gene'),
      textOutput('warningsize'),
      actionButton('update3','Update the table'),
      actionButton('update2','Update the graph'),
      fluidRow(column(3, dataTableOutput('dtints')), column(9, visNetworkOutput('string')))
    ),
    tabPanel(title = 'Predicting Functions',
      'This panel deals less with network comparison and more with one of the best utilities behind network biology: predicting the function of understudied genes. Protein interactomes allow us to associate genes with other genes by their interactions, which is how we can make predictions about genes without annotations.',
      'Below, all genes without any annotations for biological process in GO are listed in the drop-down menu (this will show up after a little processing). Select any and view its predicted functions based on enrichment of its neighbors throughout the networks.',
      uiOutput('nonann'),
      plotOutput('nonannenriched')
    )
  )
)

server <- function (input,output) {
  processing <- function (data,col1=1,col2=2,convert=TRUE,mapto=NULL, mapping=fullmap) {
    #condense data
    data <- subset(data, select = c(col1,col2))
    #map IDs if necessary
    #if (convert) {
      map <- subset(mapping, select = c(2,mapto))
      splvec <- str_split(map[[2]],"; ")
      veccounts <- lengths(splvec)
      a1 = rep.int(map[[1]],veccounts)
      a2 = unlist(splvec)
      map <- data.frame(V1 = a1,
                        V2 = a2)
      dict <- setNames(map[,1], map[,2])
      data[[1]] <- c(dict[unlist(data[[1]])])
      data[[2]] <- c(dict[unlist(data[[2]])])
      data <- data[which(complete.cases(data)),]
    #} else {}
    #remove _HUMAN parts of IDs
    data[[1]] <- str_remove(data[[1]],"_HUMAN")
    data[[2]] <- str_remove(data[[2]],"_HUMAN")
    #generate the igraph object
    vec <- as.vector(t(data))
    g <- graph(edges = vec, directed = FALSE)
    g <- simplify(g)
    isolated = which(degree(g)==0);
    g <- delete.vertices(g, isolated)
    output = g
  }
  # ensgmap <- fread(file= "uniprot_human_mapping.csv", header = TRUE)
  # v <- unzip("HPA-Protein.zip")
  # tissuefiles <- lapply(v[c(5,12,22,42,43)],fread,header=FALSE)
  # tissuegraphs <- lapply(tissuefiles,processing,mapto=3, mapping=ensgmap)
  # tissuenames <- str_sub(str_remove(str_extract(v[c(5,12,22,42,43)], "(?<=/).+(?=.tsv)")," "),1,15)

  fullmap <- fread(file = "mapping.tab", header = FALSE)
  # huri <- fread(file = "huri.tsv", header = FALSE)
  # hippie <- fread(file = "hippie.txt", header = FALSE)
  # string <- fread(file = "string.gz", header = TRUE)
  # bioplex1 <- fread(file = "bioplex1.tsv", header = TRUE)
  # bioplex2 <- fread(file = "bioplex2.tsv", header = TRUE)
  # biogrid <- fread(file = "biogrid.txt", header = FALSE)
  # corum <- fread(file = "corum.txt", header = TRUE)
  # huri_g <- processing(huri,1,2,TRUE,19)
  # hippie_g <- processing(hippie,1,3,TRUE,2)
  # biogrid[[3]] <- str_extract(biogrid[[3]], "(?<=-prot:)[[:alnum:]]+(?=|)" )
  # biogrid[[4]] <- str_extract(biogrid[[4]], "(?<=-prot:)[[:alnum:]]+(?=|)" )
  # biogrid_g <- processing(biogrid,3,4,TRUE,1)
  # bioplex1[[3]] <- substr(bioplex1[[3]],1,6)
  # bioplex1[[4]] <- substr(bioplex1[[4]],1,6)
  # bioplex2[[3]] <- substr(bioplex2[[3]],1,6)
  # bioplex2[[4]] <- substr(bioplex2[[4]],1,6)
  # bioplex_g <- processing(rbind(bioplex1,bioplex2),3,4,TRUE,1)
  # extrcorum <- function (data) {
  #   l <- str_split((data[c(which(nchar(data[[6]])>6)),6])[[1]],";")
  #   l <- lapply(l,sort)
  #   test <- unlist(lapply(l,combn,m=2))
  #   data <- as.data.frame(matrix(test, ncol = 2, byrow = TRUE))
  # }
  # corum_g <- processing(extrcorum(corum),1,2,TRUE,1)
  # string <- string[which(string[[3]]>=700),]
  # string[[1]] <- str_remove(string[[1]],"9606.")
  # string[[2]] <- str_remove(string[[2]],"9606.")
  # string_g <- processing(string,1,2,TRUE,21)
  # databasegraphs <- list(huri_g,hippie_g,biogrid_g,bioplex_g,corum_g,string_g)
  # databasenames <- c('HuRI', 'HIPPIE', 'BioGrid', 'BioPlex', 'CORUM', 'STRING')
  netdata <- eventReactive(input$update, {
    choice = as.numeric(input$data)
    if (input$data == 1) {
      #req(input$upload)
      files = input$upload$datapath
      files <- lapply(files, fread, header=input$header)
      graphs = lapply(files, processing, col1 = input$column1, col2 = input$column2, convert = input$convert, mapto = as.numeric(input$mapto))
      names = input$upload$name
      output=list(graphs,names)
    } else {
      if (input$data == 2) {
        output=list(tissuegraphs,tissuenames)
      } else {
        output=list(databasegraphs,databasenames)
      }
    }
  })
  lengthdata = reactive(length(netdata()[[1]]))
  createtable <- function (l, n) {
    extractdata <- function (g) {
      info <- vector ()
      info[1] <- NA
      info[2] <- length(V(g))
      info[3] <- length(E(g))
      info[4] <- round(mean(degree(g)), digits = 2)
      # number of components
      info[5] <- components(g)[[3]]
      # clusters by louvain, includes all components
      clustered = cluster_louvain(g)
      info[6] <- length(clustered)
      info[7] <- round(modularity(clustered), digits = 4)
      output=info
    }
    t <- t(sapply(l,extractdata))
    colnames(t) <- c('Databases', 'Nodes', 'Edges', 'Avg. Degree',
                     '# of Components', '# of Clusters (by Louvain)','Modularity')
    t[,1] = n
    output = t
  }
  output$info = renderTable(createtable(netdata()[[1]],netdata()[[2]]))
  corrmatrix <- function (a, n) {
    f <- function (v) {
      i <- intersection(a[[v[1]]],a[[v[2]]],keep.all.vertices = FALSE)
      v[3] = round((length(V(i))/length(V(a[[v[1]]])))*100, digits=2)
      v[4] = round((length(V(i))/length(V(a[[v[2]]])))*100, digits=2)
      v[5] = round((length(E(i))/length(E(a[[v[1]]])))*100, digits=2)
      v[6] = round((length(E(i))/length(E(a[[v[2]]])))*100, digits=2)
      v[7] = length(V(i))
      v[8] = length(E(i))
      output = v
    }
    l <- (combn(1:length(a),2,f))
    mvp = mep = matrix(100,nrow=length(a),ncol=length(a))
    mv = me = matrix(NA,nrow=length(a),ncol=length(a))
    for (i in 1:ncol(l)) {
      mvp[l[1,i],l[2,i]] = l[3,i]; mvp[l[2,i],l[1,i]] = l[4,i] 
      mep[l[1,i],l[2,i]] = l[5,i]; mep[l[2,i],l[1,i]] = l[6,i]
      mv[l[1,i],l[2,i]] = l[7,i]; me[l[1,i],l[2,i]] = l[8,i]}
    colnames(mvp) = rownames(mvp) = colnames(mep) = rownames(mep) = 
      colnames(mv) = rownames(mv) = colnames(me) = rownames(me) = n
    output = list(mvp,mep,mv,me)
  }
  colorforcm = colorRampPalette(c('blue', 'red'))
  cm <- reactive(corrmatrix(netdata()[[1]],netdata()[[2]]))
  output$correlation = renderPlot({
    par(mfrow=c(2,2))
    corrplot(cm()[[1]], is.corr=FALSE, title = "Nodes", diag = FALSE, mar=c(0,0,2,0),
             col = colorforcm(20), addCoef.col = "white", method = 'color', cl.pos = "n", 
             number.cex = 0.75, number.digits = 1)
    corrplot(cm()[[2]], is.corr=FALSE, title = "Edges", diag = FALSE, mar=c(0,0,2,0),
             col = colorforcm(20), addCoef.col = "white", method = 'color', cl.pos = "n", 
             number.cex = 0.75, number.digits = 1)
    corrplot(cm()[[3]], is.corr=FALSE, diag = FALSE, type = 'upper',
             col = colorforcm(20), addCoef.col = "white", method = 'color', cl.pos = "n", 
             number.cex = 0.75, number.digits = 1)
    corrplot(cm()[[4]], is.corr=FALSE, diag = FALSE, type = 'upper',
             col = colorforcm(20), addCoef.col = "white", method = 'color', cl.pos = "n", 
             number.cex = 0.75, number.digits = 1)
  }, height = function() {150*ncol(cm()[[1]])})
  output$dend = renderPlot({
    par(mfrow=c(1,2))
    plot(hclust(as.dist(100-cm()[[1]])), main = "Nodes")
    plot(hclust(as.dist(100-cm()[[2]])), main = "Edges")
  })
  comparisons <- function (l, n) {
    l = setNames(l, n)
    m = make_comb_mat(l)
    x = length(l)
    v <- matrix(c(
      comb_size(m,1), 
      rep(comb_size(m, x), times = x), 
      (set_size(m)-comb_size(m,1)-comb_size(m,x)),
      set_size(m)
    ), byrow = TRUE, nrow = 4, ncol = x)
    rownames(v) <- c('Unique','Shared by all','Remaining','Total')
    colnames(v) <- n
    output = v
  }
  f <- function (g) {
    data = as_edgelist(g)
    want2 = which(data[,1] > data[,2])
    if (length(want2) > 0) {
      data <- rbind(t(apply(data[want2,], 1, sort)),data[-want2,])
    }
    i <- paste(data[,1],data[,2],sep="-")
    output = i
  }
  getPage<-function(l, n) {
    l = setNames(l, n)
    myv = plotVenn(l, labelRegions = FALSE, systemShow=FALSE, outFile = 'testsvg')
    if (length(l) == 6) {plotVenn(l, myv, labelRegions = FALSE, systemShow=FALSE, outFile = 'testsvg')} else {}
    return(includeHTML("testsvg"))
  }
  output$euler1<-renderUI({getPage(lapply(netdata()[[1]],function(g) V(g)$name), netdata()[[2]])})
  comps <- reactive(comparisons(lapply(netdata()[[1]],function(g) V(g)$name), netdata()[[2]]))
  output$bar1 = renderPlot(barplot(comps()[1:3,], beside = FALSE, legend.text=TRUE,
          border = NA, main = "Node Comparison"))
  output$uniqshared1 = renderDataTable(t(comps()))
  alpha <- reactive(lapply(netdata()[[1]], f))
  comps2 <- reactive(comparisons(alpha(), netdata()[[2]]))
  output$euler2<-renderUI({getPage(alpha(), netdata()[[2]])})
  output$bar2 = renderPlot(barplot(comps2()[1:3,], beside = FALSE, legend.text=TRUE,
                                   border = NA, main = "Edge Comparison"))
  output$uniqshared2 = renderDataTable(t(comps2()))
  
  retrieveunique <- function (l) {
    e = l
    v = l2 = lapply(l,function(g) V(g)$name)
    n = length(l)
    for (i in 1:n) {
      for (j in 1:n) {
        if (j != i) {
          e[[i]] = difference(e[[i]],l[[j]])
          v[[i]] = setdiff(v[[i]],l2[[j]])} else {}
      }
      e[[i]] = delete.vertices(e[[i]], which(degree(e[[i]])==0))
    }
    output = list(e,v)
  }
  tsinteractioncount <- function (l) {
    t = as.matrix(table(unlist(lapply(l,function(g) c(as_edgelist(g)[,1],as_edgelist(g)[,2])))))
    output = t[order(-t[,1]),]
  }
  uniques = reactive(retrieveunique(netdata()[[1]]))
  tsintcount = reactive(tsinteractioncount(uniques()[[1]]))
  output$tsintcount = renderDataTable(
    data.frame(tsintcount())
  )
  
  tissuecount <- function (l, n, edges = FALSE) {
    nodes <- function(g) {data.frame(V(g)$name, rep(1,times=length(V(g))))}
    if (edges) {
      nodes <- function(v) {data.frame(v, rep(1, times=length(v)))}
    }
    l <- lapply(l, nodes)
    for (i in 1:length(l)) {colnames(l[[i]])=c('Names',n[i])}
    merged <- Reduce(function(x,y) merge(x = x, y = y, by = "Names", all=TRUE), l)
    count <- apply(merged[,-1],1,sum,na.rm=TRUE)
    names(count) = merged[,1]
    output = count
  }
  output$histograms = renderPlot({
    par(mfrow = c(1,2))
    hist(tissuecount(netdata()[[1]],netdata()[[2]]), breaks=lengthdata(), xlab = '# of networks found in', main = 'By Node')
    hist(tissuecount(alpha(),netdata()[[2]], edges=TRUE), breaks=lengthdata(), xlab = '# of networks found in', main = 'By Edge')
  })
  
  cv <- function (v) {100*sd(v)/mean(v)}
  dc <- function (l, n) {
    l <- lapply(l, function(g) data.frame(names(igraph::degree(g)), igraph::degree(g)))
    for (i in 1:length(l)) {colnames(l[[i]])=c('Names',n[i])}
    merged <- Reduce(function(x,y) merge(x = x, y = y, by = "Names", all=FALSE), l)
    cv <- apply(merged[,-1],1,function (a) cv(a))
    r <- apply(merged[,-1],1,function (a) max(a)-min(a))
    merged[,ncol(merged)+1] = cv
    merged[,ncol(merged)+1] = r
    colnames(merged) = c('Names',n,'CV','Range')
    output = merged[order(-merged$CV,-merged$Range),]
  }
  dctable <- reactive(dc(netdata()[[1]], netdata()[[2]]))
  dcmatrix <- reactive(as.matrix(dctable()[,-c(1,ncol(dctable()),ncol(dctable())-1)][which(dctable()[,ncol(dctable())]>0),]))
  #converts matrix to respective Z-scores to plot in heatmap
  dczs = reactive(t(apply(dcmatrix(), 1, function(v) unlist(lapply(1:length(v), function(i) if(sd(v) == 0) {0} else {(v[i] - mean(v))/sd(v)})))))
  output$dctable = renderDataTable(dctable())
  output$dcheatmap = renderPlot({
    Heatmap(dczs(), cluster_columns = FALSE, show_row_names=FALSE,
            column_names_gp = gpar(fontsize = 8), heatmap_legend_param = list(title = 'Z-score'),
            column_labels = colnames(dcmatrix()),
            col=colorRamp2(c(-2,0,2), c("blue", "white", "red"))
    )
  })
  dcps = reactive({
    z = t(apply(dcmatrix(), 1, function(v) unlist(lapply(1:length(v), function(i) abs(cv(v)-cv(v[-i]))))))
    rownames(z) = dctable()[,1]
    colnames(z) = colnames(dcmatrix())
    output=z
  })
  makesankey <- function (data, cutoff, l) {
    links <- data %>% 
      as.data.frame() %>% 
      rownames_to_column(var="source") %>% 
      gather(key="target", value="value", -1)
    links = split(links,rep(1:l,each=nrow(links)/l))
    links = lapply(links, function (d) top_n(d, n=cutoff, wt=d[,3]))
    links = do.call('rbind',links)
    nodes <- data.frame(
      name=c(as.character(links$source), as.character(links$target)) %>% 
        unique()
    )
    links$IDsource <- match(links$source, nodes$name)-1 
    links$IDtarget <- match(links$target, nodes$name)-1
    #output = list(links,nodes)
    fig <- plot_ly(type = "sankey", orientation = "h",
                   node = list(label = nodes[[1]]),
                   link = list(source = links[[4]], target = links[[5]], value =  links[[3]])
    )
    fig <- layout(fig, font = list(size = 7))
    output=fig
  }
  sankeyp = reactive(makesankey(dcps(), input$cutoff, lengthdata()))
  output$sankey = renderPlotly(sankeyp())
  dbs <- "GO_Biological_Process_2021"
  plotE <- function (d) {plotEnrich(d, orderBy = 'Combined.Score', y = 'Ratio')}
  shared <- reactive({
    shared = do.call('intersection', c(netdata()[[1]], keep.all.vertices = FALSE)) 
    shared = delete.vertices(shared, which(degree(shared)==0))
    output = shared
  })
  highdc <- reactive(dctable()[1:round(0.05*nrow(dctable())),1])

  uniqueints <- reactive({
    output = retrieveunique(netdata()[[1]])[[2]]
  })
  hubsby <- function (g) {
    d = sort(degree(g), decreasing=TRUE)
    c = sort(closeness(g), decreasing=TRUE)
    #b = estimate_betweenness(g, cutoff=5)
    top = round(0.01*length(V(tissuegraphs[[1]])))
    output = union(names(d[1:top]),names(c[1:top]))
  }
  hubs <- reactive({
    lapply(netdata()[[1]],hubsby)
  })
  graphsaslist = reactive({
    l = 1:lengthdata()
    names(l) = netdata()[[2]]
    output = l
  })
  output$tochoosenet = renderUI({
    if ((input$genes == 3) | (input$genes == 4)) {output = selectInput('nettoshow', 'Choose Graph', choices = graphsaslist())}
  })
  net <- reactive({
    if (is.null(input$nettoshow)) {output = 1} else {output = as.numeric(input$nettoshow)}
  })
  egenes <- reactive(list(names(V(shared())),highdc(),uniqueints()[[net()]],hubs()[[net()]]))
  #note: the [[1]] at the end of each of these means it's only printing the biological process one
  egenes1 <- reactive(enrichr(names(V(shared())),dbs)[[1]])
  egenes2 <- reactive(enrichr(highdc(),dbs)[[1]])
  egenes3 <- reactive(
    lapply(1:length(uniqueints()), function(i) enrichr(uniqueints()[[i]],dbs)[[1]])
  )
  egenes4 <- reactive(
    lapply(1:length(hubs()), function(i) enrichr(hubs()[[i]],dbs)[[1]])
  )
  enriched <- reactive(list(egenes1(),egenes2(),egenes3()[[net()]],egenes4()[[net()]])[[as.numeric(input$genes)]])
  output$enrichedList = renderDataTable(as.matrix(egenes()[[as.numeric(input$genes)]]))
  output$enrichedDT = renderDataTable(enriched())
  output$enrichedPlot = renderPlot(plotE(enriched()))

  
  output$choosegraph = renderUI(selectInput('graph','Choose Graph',choices=graphsaslist()))
  clusters <- reactive({
    g = netdata()[[1]][[choseng()]]
    clustered = cluster_louvain(g)
    output <- lapply(1:length(clustered), function (i) 
      induced_subgraph(g,clustered[[i]]))
  })
  chosencl <- reactive({
    if (is.null(input$cluster)) {output = 1} else {output = as.numeric(input$cluster)}
  })
  choseng <- reactive({
    # if (is.null(input$graph)) {output = 1} else {output = as.numeric(input$graph)}
    output = as.numeric(input$graph)
  })
  output$choosecluster = renderUI(numericInput('cluster','Choose Cluster', value=1, min=1, max=length(clusters())))
  output$plotcluster = renderPlot({
    clusterg = clusters()[[chosencl()]]
    output = plot(clusterg, vertex.label = NA, vertex.size = 2)
  })
  netplot <- eventReactive(input$view, {
    netplot() <- plot(netdata()[[1]][[choseng()]], vertex.label=NA, vertex.shape = 'none')
  })
  output$fullnetwork = renderPlot(netplot())
  
  enrichedclusters <- reactive({
    g = netdata()[[1]][[choseng()]]
    clustered = cluster_louvain(g)  
    enriched = lapply(1:length(clustered), function (i)
      e = enrichr(clustered[[i]],'GO_Biological_Process_2021')[[1]])
    retrieveuniqueenrichments <- function (v) {
      l2 <- v
      n = length(v)
      for (i in 1:n) {
        for (j in 1:n) {
          if (j != i) {
            v[[i]] = setdiff(v[[i]],l2[[j]])} else {}
        }
      }
      output = v
    }
    x = retrieveuniqueenrichments(lapply(1:length(enriched),function(i) enriched[[i]][[1]]))
    for (i in 1:length(enriched)) {enriched[[i]]=enriched[[i]][which(enriched[[i]][,1] %in% x[[i]]),]}
    output = enriched
  })
  output$clusterenrichment = renderPlot({
    plotE(enrichedclusters()[[chosencl()]])
  })
  output$clusterints = renderDataTable({
    as_edgelist(clusters()[[chosencl()]])
  })
  output$clusterinfo = renderDataTable({
    createtable(clusters(),1:length(clusters()))
  })
  
  #vector with names of proteins in all networks
  allproteins = reactive(V(do.call('intersection',netdata()[[1]]))$name)
  getints <- function (gene) {
    x = lapply(netdata()[[1]], function(g) if (gene %in% V(g)$name) {neighbors(g,gene)$name} else {})
    y = rep.int(netdata()[[2]], lengths(x))
    output = cbind(rep(gene,times=length(y)),unlist(x),y)
  }
  getallints <- function (gene) {
    neighs = lapply(netdata()[[1]],function(g) if (gene %in% V(g)$name) {neighbors(g,gene)$name} else {})
    neighs = unique(unlist(neighs))
    neighints = lapply(neighs, getints)
    neighints = do.call('rbind', neighints)
    output = rbind(getints(gene), neighints)
  }
  links = reactive({
    if (input$expand) {output = getallints(input$gene)} else {output = getints(input$gene)}
  })
  output$warningsize = renderText({
    if (nrow(links()) > 500) {
      output = 'Warning! This network has over 500 interactions and may be extremely slow or crash the program when you update the graph. It is better to just update the table of interactions below.'
    }
  })
  dtints <- eventReactive(input$update3, {output = links()})
  output$dtints = renderDataTable(dtints())
  stringstyle <- eventReactive(input$update2,{
    if (input$switch) {a = c(3,2)} else {a = c(2,3)}
    nodes = data.frame(id = unique(c(links()[,1],links()[,a[1]])), label = unique(c(links()[,1],links()[,a[1]])))
    edges = data.frame(from = links()[,1], to = links()[,a[1]], title = links()[,a[2]])
    output = visNetwork(nodes, edges)
  })
  output$string = renderVisNetwork(stringstyle())
  nonann <- reactive({
    enriched = enrichr(allproteins(),'GO_Biological_Process_2021')[[1]][[9]]
    splvec <- unlist(str_split(enriched,";"))
    output = allproteins()[!(allproteins() %in% splvec)]
  })
  
  
  output$nonann = renderUI(selectInput('choosenonann', 'Choose a gene', choices = nonann()))
  output$nonannenriched = renderPlot({
    ns = unlist(lapply(netdata()[[1]], function (g) neighbors(g,input$choosenonann)$name))
    enrichedns = enrichr(ns,'GO_Biological_Process_2021')[[1]]
    plotE(enrichedns)
  })
}

shinyApp(ui=ui, server=server)