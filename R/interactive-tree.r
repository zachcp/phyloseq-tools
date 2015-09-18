#' create an interactive (html5) phylogenetic tree 
#'
#' interactive phlogenetic tree that returns a phylogenetic tree as an htmlwidget
#'
#' @import leaflet
#' @import dplyr
#' @import phyloseq
#' @importFrom htmltools htmlEscape
#'
#' @param  physeq a phyloseg object
#' @param ladderize boolean. ladderize or not
#' @param aspectratio control the dimensions of the tree
#' @param linecolor color lines
#' @param lineweight
#' @param lineopacity
#' @param size size of circles
#' @param abundancescale a sclae for circle size
#' @param method
interactivetree <- function(physeq,
                            ladderize=TRUE,
                            aspectratio=1,
                            #lineoptions
                            linecolor= "black",
                            lineweight = 2,
                            lineopacity = 1,
                            #circle options
                            size = 3,
                            abundancescale=1,
                            color = "blue",
                            stroke = TRUE,
                            fill = TRUE,
                            circlestrokeweight = 3,
                            circlestrokeopacity = 0.5,
                            fillColor = color,
                            fillOpacity = 0.2,
                            palette = "Blues",
                            numerictype = "numeric",
                            method="sampledodege",
                            legend=TRUE){
  
  #get tree layout data
  treeSegs <- phyloseq::tree_layout(physeq, ladderize = ladderize)
  
  # use the x and y max/min values to find an x and yscle factor to length 1
  # and modify these scale sby the aspect ratio
  x_max <- max(treeSegs$edgeDT$xright)
  y_max <- max(treeSegs$edgeDT$y)
  xscale = 1/x_max
  yscale = 1/y_max
  xscale = xscale*aspectratio
  
  #' helper function to create a leaflet-compatable two-column matrix of lines
  make_edgematrix <-function(){
    #get the edges
    edges <- treeSegs$edgeDT
    edges <- edges %>% add_rownames()
    start <- edges %>% mutate(lng = xleft*xscale, lat = y*yscale)
    end   <- edges %>% mutate(lng = xright*xscale, lat = y*yscale)
    space <- edges
    space$lng = NA
    space$lat = NA
    edges2 <- rbind(start, end, space) %>%
      arrange(rowname) %>%
      select(lng,lat) %>%
      as.matrix()
  }
  #' helper function to create a leaflet-compatable two-column matrix of lines
  make_vertexmatrix <- function(){
    #and the vertices
    verts <- treeSegs$vertDT
    verts <- verts %>% add_rownames()
    start <- verts %>% mutate(lng=x*xscale, lat=vmin*yscale)
    end   <- verts %>% mutate(lng=x*xscale, lat=vmax*yscale)
    space <- verts
    space$lng = NA
    space$lat = NA
    verts2 <- rbind(start, end, space) %>%
      arrange(rowname) %>%
      select(lng,lat) %>%
      as.matrix()
  }
  
  ################################################################
  ################################################################
  # get tree edges and vertices and create the basemap
  edges <- make_edgematrix()
  vertices <- make_vertexmatrix()
  m <- leaflet() %>%
    addPolylines(data=edges,    color=linecolor, weight=lineweight, opacity=lineopacity) %>%
    addPolylines(data=vertices, color=linecolor, weight=lineweight, opacity=lineopacity)
  
  # add the tree tips. tip metadata will come from either the tax_table (method=="tree")
  # or from the taxtable + sampel_data (method=="sampledodege"). in the latter case,
  # there will be multiple points correspondin to each OTU-sample combination
  points <- treeSegs$edgeDT %>%
    data.frame() %>%
    filter(!is.na(OTU)) %>%
    mutate(lng = xright*xscale, lat = y*yscale) %>%
    select(OTU,lng,lat)
  
  if(method=="tree"){
    if (size=="Abundance") stop("Abundance plotting can only be done with the sampledodge method")
    
    # merge point location with tax_table information and draw the tree
    taxdata <- tax_table(physeq) %>%
      data.frame() %>%
      mutate(OTU = rownames(.))
    pointdata <- merge(points, taxdata, by="OTU") %>%
      mutate(treelat = lat, treelng=lng) %>%
      select(-lat,-lng)
    
  } else if(method=="sampledodge"){
    warning("This method is not yet implemented")
    points <- points %>%
      mutate(treelat = lat, treelng=lng) %>%
      select(-lat,-lng)
    melted <- phyloseq::psmelt(physeq)
    pointdata <- merge(points,melted, by="OTU")
    pointdata <- pointdata %>% group_by(OTU) %>%
      mutate(rownum = row_number(treelng)) %>%
      mutate(treelng = treelng+(rownum*0.003))
    
  } else {
    stop("'sampledoge' and 'tree' are the only acceptable method options")
  }
  
  ################################################################
  ################################################################
  #check size input
  if (is.numeric(size)){
    pointdata$pointsize <- size*100
  } else if (size == "Abundance"){
    pointdata$pointsize <- pointdata$Abundance *abundancescale
    #pointdata$pointsize <- lapply(pointdata$Abundance, function(x) sqrt(x) * abundancescale)
  } else {
    stop("Invalid size option: numeric value or 'Abundance")
  }
  
  ################################################################
  ################################################################
  # check color varaialbe and add a color column to the points df
  if(!color %in% names(pointdata)){
    pointdata$color <- color
  }else{
    #checkfactor
    if(is.factor(pointdata[[color]])){
      colorfn <- leaflet::colorFactor(palette = palette,
                                      domain = pointdata[[color]])
      pointdata$color <- colorfn(pointdata[[color]])
    }else if(is.numeric(pointdata[[color]]) & numerictype=="numeric"){
      colorfn <- leaflet::colorNumeric(palette = palette,
                                       domain = pointdata[[color]])
      pointdata$color <- colorfn(pointdata[[color]])
    }else if(is.numeric(pointdata[[color]]) & numerictype=="bin"){
      colorfn <- leaflet::colorBin(palette = palette,
                                   domain = pointdata[[color]])
      pointdata$color <- colorfn(pointdata[[color]])
    }else if(is.numeric(pointdata[[color]]) & numerictype=="quantile"){
      colorfn <- leaflet::colorQuantile(palette = palette,
                                        domain = pointdata[[color]])
      pointdata$color <- colorfn(pointdata[[color]])
    }
  }
  
  #set fillcolor
  if (fillColor == color){
    pointdata$fillcolor <- pointdata$color
  } else {
    pointdata$fillcolor <- fillColor
  }
  
  
  ################################################################
  ################################################################
  #function for useful popups
  createPopup <- function(){
    
  }
  ################################################################
  ################################################################
  #add points to the map
  m <- m %>% addCircles(data=pointdata,
                        lng=~treelng,
                        lat=~treelat,
                        radius=~pointsize,
                        color = ~color,
                        stroke = stroke,
                        fill = fill,
                        weight = circlestrokeweight,
                        opacity = circlestrokeopacity,
                        fillColor = ~fillcolor,
                        fillOpacity = fillOpacity,
                        popup=~sprintf("<b>OTU Data</b><br/>
                                       Sample: %s </br>
                                       Location:  %s, %s",
                                       htmlEscape(Sample),
                                       htmlEscape(treelng), 
                                       htmlEscape(treelat))
  )
  
  print(pointdata)
  print(names(pointdata))
  ################################################################
  ################################################################
  #add legend
  if (legend & exists("colorfn")){
    m <- m %>% addLegend(position="bottomright",
                         pal = colorfn,
                         values = pointdata[[color]],
                         opacity = 1)
  }
  
  return(m)
}

