library(shiny)
library(shinycssloaders)
library(DT)
library(plotly)
options(shiny.maxRequestSize = 30*1024^2)

ui <- fluidPage(
  tabsetPanel(
    tabPanel(title = "Input Data",
      titlePanel('Tool for Computational Protein Network Analysis'),
      'Welcome! Begin by uploading the files with the data that you want to analyze.',
      'For now, the app is set to process data from TissueNet, so please upload that.',
      'Once the data is processed, a table will show up with some basic statistics about each network.',
      'Then you can navigate to any of the other tabs to find more analyses!',
      fileInput(inputId = "upload", label = "Upload files", multiple = TRUE),
      withSpinner(tableOutput("info"))
    ),
    tabPanel(title = 'Network Visualization (by Cluster)',
      'Here you can visualize each graph by clusters, which are pre-determined by the Louvain method.',
      'This is done because the whole network altogether is quite large and will thus look messy and take a long time to plot.',
      'By the way, an error may pop up momentarily before the initial processing, that can be ignored.',
      uiOutput('choosegraph'),
      uiOutput('choosecluster'),
      withSpinner(plotOutput('plotcluster'))
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
      plotOutput('euler1'),
      plotOutput('bar2'),
      dataTableOutput('uniqshared2'),
      plotOutput('euler2'),
      titlePanel('Network Specifity Measures'),
      'Below is a table listing each gene and how many network-specific interactions it is a part of.',
      dataTableOutput('tsintcount'),
      'The histograms below show the distribution of specificity of nodes and edges across the networks. They essentially answer the question of how network-specific/unique nodes and edges are throughout.',
      'More left-leaning histograms mean generally more unique, while right-leaning means generally more global.',
      plotOutput('histograms')
      #plotOutput('neighanal')
    ),
    tabPanel(title = 'Differential Connectivity',
      'This tab focuses on differential connectivity, the idea of how a gene\'s importance, determined by its degree (the number of genes it interacts with), varies across each network. If a gene has a low degree in one network and high in another, that could signal it as a gene/protein of interest.',
      'The differential connectivity of each gene found in all networks can be found below, ranked by the coefficient of variation (CV) of its degrees.',
      dataTableOutput('dctable'),
      'This table is visualized in a heatmap below. Each row is a gene, and the columns are its degree values in each network scaled by Z-score. Columns that stick out indicate networks with generally much lower (blue) or higher (red) degrees than the other networks.',
      withSpinner(plotOutput('dcheatmap')),
      #radioButtons('disttype', label = 'Choose type of distribution to filter by', 
      #             choiceNames = c('Poisson','Normal'), choiceValues = c(1,2)),
      'The next question is if a gene has high differential connectivity across networks, which network(s) is or are causing that variance in degree? If a network is causing high variance of a gene (i.e. the gene\'s degree changes heavily in that network), that gene can be associated with that network like a hub.',
      'This is visualized below in a Sankey plot mapping each network to its top \'hubs\' in this manner. This was determined by measuring the change in differential connectivity of a gene when the degree of that gene in each network was not counted toward the variance calculation.',
      'You may use the slider below to change how many genes show up for each network. Note that the number of actual genes that show up may be greater if there is a tie among genes (based on the score described before). You may also view that score for each association by hovering over the respective line.',
      sliderInput("cutoff", label = 'Top # of Associations Per Network', min = 1, max = 10, value = 3),
      withSpinner(plotlyOutput('sankey'))
    ),
    tabPanel(title = 'Isolate and Enrich Genes of Interest',
      selectInput('genes', 'Genes to Enrich', 
                  choices = list('Shared Interactions' = 1,
                                 'High DC' = 2,
                                 'Unique Interactions' = 3)),
      uiOutput('tochoosenet'),
      # conditionalPanel("input.genes == 3",
      #   numericInput('nettoshow', 'Which Network?', value = 1)),
      checkboxGroupInput('test', label = 'Test',
                         choices = list('GO_Biological_Process_2021', 
                                        'GO_Molecular_Function_2021',
                                        'GO_Cellular_Component_2021'),  selected = 1),
      dataTableOutput('enrichedList'),
      dataTableOutput('enrichedDT'),
      plotOutput('enrichedPlot')
    )
  )
)

server <- function (input,output) {
  setwd("~/R")
  library(igraph)
  library(data.table)
  library(R.utils)
  #library(VennDiagram)
  library(stringr)
  library(corrplot)
  #library(UpSetR)
  library(ComplexHeatmap)
  library(enrichR)
  library(dplyr)
  #library(readxl)
  #library(rmarkdown)
  #library(devtools)
  library(networkD3)
  library(tidyr)
  library(plotly)
  library(tibble)
  library(circlize)
  library(distr)
  library(visNetwork)
  library(eulerr)
  library(viridis)
  degree <- igraph::degree
  sd <- stats::sd
  processing <- function (data,col1=1,col2=2,convert=TRUE,mapto=NULL) {
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
  mapping <- fread(file= "uniprot_human_mapping.csv", header = TRUE)
  netdata <- reactive({
    req(input$upload)
    files = input$upload$datapath
    files <- lapply(files, fread, header=FALSE)
    graphs = lapply(files, processing, mapto=3)
    names = str_sub(str_extract(input$upload$name, ".+(?=.tsv)"),1,15)
    #graphs = 1:3
    output=list(graphs,names)
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
    z = str_remove_all(comb_name(m, readable=TRUE), "(?<=\\()[[^&]]+(?=\\))")
    z = str_remove_all(z, '&\\(\\)|\\(\\)&')
    x = comb_size(m); names(x) = z
    e = euler(x, shape = 'ellipse')
    x = length(l)
    v <- matrix(c(
      comb_size(m,1), 
      rep(comb_size(m, x), times = x), 
      (set_size(m)-comb_size(m,1)-comb_size(m,x)),
      set_size(m)
    ), byrow = TRUE, nrow = 4, ncol = x)
    rownames(v) <- c('Unique','Shared by all','Remaining','Total')
    colnames(v) <- n
    output = list(e,v)
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
  up <- reactive(comparisons(lapply(netdata()[[1]],function(g) V(g)$name), netdata()[[2]]))
  output$euler1 = renderPlot({
    par(mfrow=c(1,1))
    plot(up()[[1]], shape='ellipse', labels = FALSE, legend=TRUE, 
         fills = viridis(lengthdata()), edges = viridis(lengthdata()), main = "Euler of Nodes")
  })
  output$bar1 = renderPlot(barplot(up()[[2]][1:3,], beside = FALSE, legend.text=TRUE,
          border = NA, main = "Node Comparison"))
  output$uniqshared1 = renderDataTable(t(up()[[2]]))
  alpha <- reactive(lapply(netdata()[[1]], f))
  up2 <- reactive(comparisons(alpha(), netdata()[[2]]))
  output$euler2 = renderPlot({
    par(mfrow=c(1,1))
    plot(up2()[[1]], shape='ellipse', labels = FALSE, legend=TRUE, 
         fills = viridis(lengthdata()), edges = viridis(lengthdata()), main = "Euler of Edges")
  })
  output$bar2 = renderPlot(barplot(up2()[[2]][1:3,], beside = FALSE, legend.text=TRUE,
                                   border = NA, main = "Edge Comparison"))
  output$uniqshared2 = renderDataTable(t(up2()[[2]]))
  
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
    hist(tissuecount(netdata()[[1]],netdata()[[2]]), breaks=lengthdata(), xlab = '# of tissues found in', main = 'By Node')
    hist(tissuecount(alpha(),netdata()[[2]], edges=TRUE), breaks=lengthdata(), xlab = '# of tissues found in', main = 'By Edge')
  })
  
  # retrieveallneighbors <- function (n,g) {
  #   neighbor = function (g) {if (n %in% V(g)$name) {output = neighbors(g,n)$name} else {}}
  #   output = unique(unlist(lapply(g,neighbor)))
  # }
  # neighanal = reactive(unlist(lapply(names(tsintcount()), function(n) 
  #   if (n %in% unlist(uniques()[[2]])) {NA} else {tsintcount()[n]/length(retrieveallneighbors(n, g=netdata()[[1]]))}
  #   ))
  # )
  # output$neighanal = renderPlot(hist(neighanal(), breaks=100, xlab = 'Fraction of Interactions that are TS',
  #      main='Distribution of Interactions that are TS'))
  
  cv <- function (v) {100*sd(v)/mean(v)}
  dc <- function (l, n) {
    l <- lapply(l, function(g) data.frame(names(igraph::degree(g)), igraph::degree(g)))
    for (i in 1:length(l)) {colnames(l[[i]])=c('Names',n[i])}
    merged <- Reduce(function(x,y) merge(x = x, y = y, by = "Names", all=FALSE), l)
    #merged <- Reduce(function(x,y) merge(x = x, y = y, by = "Names", all=TRUE), l)
    #merged[is.na(merged)] <- 0
    cv <- apply(merged[,-1],1,function (a) cv(a))
    r <- apply(merged[,-1],1,function (a) max(a)-min(a))
    merged[,ncol(merged)+1] = cv
    merged[,ncol(merged)+1] = r
    colnames(merged) = c('Names',n,'CV','Range')
    #first item is sorted by coeff. of variation, second sorted by range
    output = merged[order(-merged$CV,-merged$Range),]
  }
  dctable <- reactive(dc(netdata()[[1]], netdata()[[2]]))
  dcmatrix <- reactive(as.matrix(dctable()[,-c(1,ncol(dctable()),ncol(dctable())-1)]))
  #converts matrix to respective Z-scores to plot in heatmap
  dczs = reactive(t(apply(dcmatrix(), 1, function(v) unlist(lapply(1:length(v), function(i) if(sd(v) == 0) {0} else {(v[i] - mean(v))/sd(v)})))))
  output$dctable = renderDataTable(dctable())
  output$dcheatmap = renderPlot({
    Heatmap(dczs(), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names=FALSE,
            column_names_gp = gpar(fontsize = 8), heatmap_legend_param = list(title = 'Z-score'),
            column_labels = colnames(dcmatrix()),
            col=colorRamp2(c(-2,0,2), c("blue", "white", "red"))
    )
  })
  #dcps = reactive(t(apply(dcmatrix(), 1, function(v) ppois(v, lambda = mean(v)))))
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
  # dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
  #          "GO_Biological_Process_2021")
  #dbs = reactive(input$test)
  dbs <- "GO_Biological_Process_2021"
  plotE <- function (d) {plotEnrich(d, orderBy = 'Combined.Score', y = 'Ratio')}
  shared <- reactive({
    shared = do.call('intersection', c(netdata()[[1]], keep.all.vertices = FALSE)) 
    shared = delete.vertices(shared, which(degree(shared)==0))
    output = shared
  })
  highdc <- reactive(dctable()[1:round(0.05*nrow(dctable())),1])

  uniqueints <- reactive({
    u = retrieveunique(netdata()[[1]])[[1]]
    output=lapply(1:length(u), function (i) V(u[[i]])$name)
  })
  hubsbydcandb <- function (g) {
    d = degree(g); c = closeness(g); b = estimate_betweenness(g, cutoff=10)
    m = as.data.frame(cbind(d, c = c[names(d)], b = b[names(d)]))
    topfive = round(0.05*nrow(m))
    output = Reduce(union, list(rownames(m[order(-m$d)[1:topfive],]), 
                                rownames(m[order(-m$c)[1:topfive],]),
                                rownames(m[order(-m$b)[1:topfive],])))
  }
  hubs <- reactive({
    lapply(netdata()[[1]],hubsbydandc)
  })
  output$tochoosenet = renderUI({
    if (input$genes == 3) {output = numericInput('nettoshow', 'Which Network?', value = 1, min = 1, max = lengthdata())}
  })
  net <- reactive({
    if (is.null(input$nettoshow)) {output = 1} else {output = as.numeric(input$nettoshow)}
  })
  egenes <- reactive(list(names(V(shared())),highdc(),uniqueints()[[net()]]))
  #note: the [[1]] at the end of each of these means it's only printing the biological process one
  egenes1 <- reactive(enrichr(names(V(shared())),dbs)[[1]])
  egenes2 <- reactive(enrichr(highdc(),dbs)[[1]])
  egenes3 <- reactive(
    lapply(1:length(uniqueints()), function(i) enrichr(uniqueints()[[i]],dbs)[[1]])
  )
  enriched <- reactive(list(egenes1(),egenes2(),egenes3()[[net()]])[[as.numeric(input$genes)]])
  # genechoice <- reactive(as.numeric(input$genes))
  # e = list (NULL, NULL)
  # enriched <- reactive({
  #   if (genechoice() == 1) {if (is.null(e[[1]])) {e[[1]] = enrichr(names(V(shared())),dbs)}}
  #   if (genechoice() == 2) {e[[2]] = enrichr(highdc(),dbs)}
  #   output = e[[genechoice()]][[1]]
  # })
  output$enrichedList = renderDataTable(as.matrix(egenes()[[as.numeric(input$genes)]]))
  output$enrichedDT = renderDataTable(enriched())
  output$enrichedPlot = renderPlot(plotE(enriched()))
  
  
  graphsaslist = reactive({
    l = 1:lengthdata()
    names(l) = netdata()[[2]]
    output = l
  })
  output$choosegraph = renderUI(selectInput('graph','Choose Graph',choices=graphsaslist()))
  clusters <- reactive({
    g = netdata()[[1]][[as.numeric(input$graph)]]
    clustered = cluster_louvain(g)
    output <- lapply(1:length(clustered), function (i) 
      induced_subgraph(g,clustered[[i]]))
  })
  makeVis <- function (g) {
    nodes = data.frame(id = V(g)$name, title = V(g)$name)
    t = as_edgelist(g)
    edges = data.frame(from = t[,1], to = t[,2])
    output = visNetwork(nodes, edges)
  }
  output$choosecluster = renderUI(numericInput('cluster','Choose Cluster', value=1, min=1, max=length(clusters())))
  output$plotcluster = renderPlot({
    clusterg = clusters()[[as.numeric(input$cluster)]]
    output = plot(clusterg, vertex.label = NA, vertex.size = 2)
  })
  # output$plotcluster = renderPlotly({
  #   clusterg = clusters()[[as.numeric(input$cluster)]]
  #   if (length(E(clusterg)) > 500) {output = plot(clusterg, vertex.label = NA, vertex.size = 2)}
  #     else {output = makeVis(clusterg)}
  # })
  # output$plotcluster = renderPlot({
  #   par(mfrow = c(ceiling(length(clusters())/3),3))
  #   for (i in 1:length(clusters())) {plot(clusters()[[i]], vertex.label = NA, vertex.shape = 'none')}
  # })
}

shinyApp(ui=ui, server=server)