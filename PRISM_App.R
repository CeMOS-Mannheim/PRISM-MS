library(shiny)
library(dplyr)
library(msiImporter)
library(xml2)
library(stringr)
library(Cardinal)
library(shinyjs)
library(EBImage)

customCSS <- "
body { 
  background-image: url('background.png'); 
  background-size: cover; 
  color: #333; 
  font-family: 'Arial', sans-serif; 
}
.h4, .h3, .h2, .h1, h1, h2, h3, h4 {
  color: #337ab7;
  font-size: 200%; /* Increase title size */
  font-weight: bold; /* Make title */
}
.btn {
  background-color: #337ab7; 
  color: #fff; 
  border-radius: 4px; 
  transition: background-color 0.5s ease; /* Smooth transition for color change */
}
.btn:hover {
  background-color: #285f8f; /* Darker shade on hover */
}
.btn:active {
  background-color: #1d4f73; /* Even darker shade to provide feedback on click */
  border-color: #122b40; /* Optional: change border color on click for more feedback */
}
.nav-tabs > li > a {
  color: #337ab7; 
  font-size: 16px; 
  border: 1px solid #ddd;
  border-bottom-color: transparent;
}
.nav-tabs > li.active > a {
  font-weight: bold;
}
.sidebar-panel {
  background-color: #e9ecef;
  padding: 20px;
  border-radius: 5px;
}
.titlePanel {
  color: #fff; 
  background-color: #337ab7; 
  padding: 20px;
  border-radius: 0;
  border-bottom: 5px solid darken(#337ab7, 10%);
}
.main-panel {
  padding: 20px;
}
/* Custom styles for progress bars */
.shiny-progress-container {
  height: 30px; /* Make the progress bar taller */
  border: 1px solid #337ab7; /* Border color matches button color */
  border-radius: 5px; /* Rounded corners for the progress bar */
  background-color: #f7f7f7; /* Light background color */
}
.shiny-progress-bar {
  background-color: #337ab7; /* Progress bar color matches button color */
  height: 100%; /* Ensure the progress bar fills the container */
  border-radius: 4px; /* Rounded corners for the progress bar */
}
.shiny-progress-text {
  position: absolute;
  width: 100%;
  text-align: center;
  color: #fff; /* Text color */
  font-weight: bold; /* Make text bold */
}
.responsive-img-container img {
  width: 100%;
  height: auto;
}
"

ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML(customCSS))),
  titlePanel("PRISM-MS"),
  uiOutput("dynamicUI")
)


server <- function(input, output, session) {
  
  #increase maximum file uploade size to 48 GB
  options(shiny.maxRequestSize=48000*1024^2)
  
  
  loaded_data <- reactiveVal()
  step <- reactiveVal(1)
  
  values <- reactiveValues(
    m_range = c(200, 1000),
    tolerance = 10,
    goalmz = 786.5,
    thresh = 1
  )
  
  observe({
    # Update reactive values when inputs change
    
    if (!is.null(input$m_range)) {
      values$m_range <- input$m_range
    }
    if (!is.null(input$tolerance)) {
      values$tolerance <- input$tolerance
    }

    debouncedGoalmz <- debounce(reactive({input$goalmz}), 2000)
    observe({
      if (!is.null(debouncedGoalmz())) {
        values$goalmz <- debouncedGoalmz()
      }
    })
    debouncedThresh <- debounce(reactive({ input$thresh }), 1000)
    observe({
      if (!is.null(debouncedThresh())) {
        values$thresh <- debouncedThresh()
      }
    })

    
  })
  
  output$dynamicUI <- renderUI({
    if (step() == 1) {
      sidebarLayout(
        sidebarPanel(
          fileInput("imzmlFile", "imzML & ibd File", multiple = TRUE, buttonLabel = "Browse Files", accept = c(".imzML", ".ibd")),
          fileInput("mis_file", "Choose MIS File", accept = c(".mis")),
          selectInput("normalize", "Normalization", choices = c("none", "TIC", "RMS", "Reference")),
          numericInput("feat", "Reference Peak",value=208.1140),
          sliderInput("m_range", "Mass Range", min = 1, max = 5000, value = c(100, 900)),
          numericInput("tolerance", "Tolerance (ppm)", value = 10),
          actionButton("load", "Load", style = "display: none;")
        ),
        mainPanel()
      )
    } else if (step() != 1) {
      sidebarLayout(
        sidebarPanel(
          numericInput("goalmz", "Goal m/z", values$goalmz, min = 1, step = 0.01),
          tags$h4("Threshold Adjustment"),
          numericInput("thresh", "Threshold Level", values$thresh,step=0.25),
          actionButton("otsu_thresh", "Set Otsu Threshold"),
          tags$hr(),
          textInput("new_method", "Choose New Method Path"),
          numericInput("resulting_res", "Resulting Spatial Resolution (μm)", 20),
          numericInput("spot_dilation", "Spot Dilation", 1),
          checkboxInput("fill_holes", "Fill Holes", TRUE),
          actionButton("write", "Write MIS File")
        ),
        mainPanel(
          fluidRow(
            column(12, plotOutput("threshPlot"))
          ),
          fluidRow(uiOutput("downloadButtonUI")),
          fluidRow(
            column(6, div(class = "responsive-img-container", imageOutput("realImageDisplay"))),
            column(6, div(class = "responsive-img-container", imageOutput("thresholdedImageDisplay")))
          )
        )
      )
    }
  })
  
  observeEvent(input$normalize, {
    if (input$normalize == "Reference") {
      shinyjs::show("feat")
    } else {
      shinyjs::hide("feat")
    }
  })
  
  xmlContent <- reactiveVal()
  
  
  observe({
    validFiles <- reactive({
      req(input$imzmlFile)
      req(input$mis_file)
      
      imzmlFiles <- input$imzmlFile$name
      misFile <- input$mis_file$name
      
      # Extract base names without extension
      imzmlBaseNames <- tools::file_path_sans_ext(imzmlFiles)
      
      # Check if there are exactly 2 imzML files (one .imzML and one .ibd) with the same base name
      validImzml <- length(unique(imzmlBaseNames)) == 1 && length(imzmlFiles) == 2 &&
        any(grepl(".imzML$", imzmlFiles)) && any(grepl(".ibd$", imzmlFiles))
      
      # Check if the .mis file matches the imzML base name
      validMis <- grepl(paste0("^", unique(imzmlBaseNames), "\\.mis$"), misFile)
      
      validImzml && validMis
    })
    
    observe({
      if(validFiles()) {
        shinyjs::show("load")
      } else {
        shinyjs::hide("load")
        output$mainPanel <- renderUI({
          "Files are not correctly uploaded or do not match the required naming conventions. Please upload again."
        })
      }
    })
  })
  
  # Step 1: Load Data & Set Parameters
  observeEvent(input$load, {
    
    
    withProgress(message = 'Loading data...', {
      incProgress(0.1, detail = "Preparing files")
      
      
      # Check if all inputs are available
      if (is.null(input$imzmlFile) | is.null(input$mis_file)) {
        print(input)
        print("something wrong with de paths")
        return()
      }
      
      # Copy imzML files to tempdir
      for(i in 1:length(input$imzmlFile$name)){
        file.copy(input$imzmlFile$datapath[i], paste0(tempdir(),"/", input$imzmlFile$name[i]))
      }
      
      # Find and set the path for the imzML file
      imzmlidx <- which(str_detect(input$imzmlFile$datapath, ".imzML"))[1]
      imzmlpath <- paste0(tempdir(),"/", input$imzmlFile$name[imzmlidx])
      
      # Copy the MIS file to tempdir and set its path
      mis_file_path <- paste0(tempdir(),"/", input$mis_file$name)
      file.copy(input$mis_file$datapath, mis_file_path)
      
      # Use the paths in msiImporter::getRegions
      regions <- msiImporter::getRegions(mis_file_path, imzmlPath = imzmlpath)
      region_names <- regions[["coordsSummary"]][["ROI_name"]]
      
      
      incProgress(0.5, detail = "Processing imzML and MIS files")
      
      msData_cardinal <- Cardinal::readMSIData(imzmlpath, attach.only=TRUE, mass.range = input$m_range,resolution = input$tolerance, units = "ppm")
      
      if(input$normalize=="none"){
        msData_cardinal <- msData_cardinal
      }else if(input$normalize=="TIC"){
        msData_cardinal <- process(normalize(msData_cardinal, method="tic"))
      }else if(input$normalize=="RMS"){
        msData_cardinal <- process(normalize(msData_cardinal, method="rms"))
      }else if(input$normalize=="Reference"){
        msData_cardinal <- process(normalize(msData_cardinal, method="reference",feature = input$feat))
      }
      
      
      incProgress(1, detail = "Finalizing")
      
      mis_file <- xml2::read_xml(mis_file_path)
      centr_pos <- xml_attr(xml_find_first(mis_file, "//View"), "CenterPos")
      centr_pos_x <- as.numeric(str_split(centr_pos,",")[[1]][1])
      centr_pos_y <- as.numeric(str_split(centr_pos,",")[[1]][2])
      
      mis_list <- as_list(mis_file)
      
      len_mis <- length(mis_list[[1]])
      
      mis_coords_orig <- list()
      
      reg_center <- list()
      for(reg in 1:length(region_names)){
        
        idx_reg <- len_mis-length(region_names)+reg
        
        p1 <- mis_list[[1]][[idx_reg]][[3]][[1]]
        p2 <- mis_list[[1]][[idx_reg]][[4]][[1]]
        
        p1 <- as.numeric(strsplit(p1,",")[[1]])
        p2 <- as.numeric(strsplit(p2,",")[[1]])
        
        
        reg_center[[length(reg_center)+1]] <- c(round((p1[1]+p2[1])/2),round((p1[2]+p2[2])/2))
        mis_coords_orig[[length(mis_coords_orig)+1]] <- c(p1,p2)
        
        
      }
      
      
      if(values$goalmz > min(Cardinal::mz(msData_cardinal)) && values$goalmz < max(Cardinal::mz(msData_cardinal))) {
        # goalmz is in the range, no action needed
      } else {
        values$goalmz <- median(Cardinal::mz(msData_cardinal))
      }
      
      
      loaded_data(list(msData_cardinal,mis_coords_orig,mis_list,region_names,reg_center,centr_pos_x,centr_pos_y))
      
    })
    # Move to step 2
    step(2)
  })
  
  im_mat_reactive <- eventReactive(values$goalmz, {
    
    req(step() >= 2, loaded_data())
    
    withProgress(message = 'Generating preview...', {
      
      incProgress(0.3, detail = "Calculating image matrix")
      
      goal_mz <- values$goalmz
      im_tmp <- Cardinal::image(loaded_data()[[1]], mz=goal_mz, plusminus=0.01)
      
      incProgress(0.7, detail = "Applying thresholds")
      
      x_vals <- im_tmp[["facets"]][[1]][[1]][["x"]]
      y_vals <- im_tmp[["facets"]][[1]][[1]][["y"]]
      im_mat <- t(im_tmp[["facets"]][[1]][[1]][["values"]])
      
      # im_mat <- scale(im_mat)
      im_mat <- (im_mat - mean(im_mat,na.rm=T)) / sd(im_mat,na.rm=T)
      
      
      incProgress(1, detail = "Ready")
    })
    
    return(im_mat)
  }, ignoreNULL = FALSE)
  
  output$realImageDisplay <- renderImage({
    req(im_mat_reactive())  # Ensure im_mat is ready
    im_mat <- EBImage::flip(t(im_mat_reactive()))
    
    col_palette <- colorRampPalette(c("blue", "white", "red"))(256)
    
    width <- 250
    
    x.values.x <- nrow(im_mat)
    x.values.y <- ncol(im_mat)
    
    # Create a PNG file
    png_filename <- "real_image.png"
    png(png_filename, width = width, height = (x.values.y * width / x.values.x), units = "px", res = 72)
    
    # Plot the image on the PNG device
    par(mar = c(0, 0, 0, 0))  # Set margins to 0 to remove axes
    image(im_mat, col = col_palette, axes = FALSE, asp = NA)
    
    # Close the PNG device
    dev.off()
    
    # Return the PNG file as binary data
    list(src = png_filename, contentType = "image/png")
  }, deleteFile = TRUE)
  
  output$thresholdedImageDisplay <- renderImage({
    req(im_mat_reactive())  # Ensure im_mat is ready
    im_mat <- EBImage::flip(t(im_mat_reactive()))
    thresh <- values$thresh  # Make sure you have the threshold value reactive or input
    im_mat <- im_mat > thresh
    col_palette <- colorRampPalette(c("black", "white"))(2)
    
    width <- 250
    
    x.values.x <- nrow(im_mat)
    x.values.y <- ncol(im_mat)
    
    # Create a PNG file
    png_filename <- "real_image.png"
    png(png_filename, width = width, height = (x.values.y * width / x.values.x), units = "px", res = 72)
    
    # Plot the image on the PNG device
    par(mar = c(0, 0, 0, 0))  # Set margins to 0 to remove axes
    image(im_mat, col = col_palette, axes = FALSE, asp = NA)
    
    # Close the PNG device
    dev.off()
    
    # Return the PNG file as binary data
    list(src = png_filename, contentType = "image/png")
  }, deleteFile = TRUE)
  
  output$threshPlot <- renderPlot({
    req(im_mat_reactive()) # Ensure im_mat is computed
    im_mat <- im_mat_reactive()
    im_mat_vector <- as.vector(im_mat)
    log_im_mat_vector <- im_mat_vector
    log_im_mat_vector <- na.omit(log_im_mat_vector)
    min_val <- min(log_im_mat_vector)
    max_val <- max(log_im_mat_vector)
    
    thresh_log <- (input$thresh)
    col_palette <- colorRampPalette(c("blue", ifelse(thresh_log > min_val & thresh_log < max_val, "white", "red"), "red"))(256)
    
    dens <- density(log_im_mat_vector)
    plot(dens, main = "Density Plot of Z-Score Scaled Intensities", xlab = "Z-Score", ylab = "Density", 
         col = "black", lwd = 2)
    abline(v = thresh_log, col = "green", lwd = 4)
    rect(min(dens$x), 0, thresh_log, max(dens$y), col = rgb(1, 0, 0, 0.2), border = NA)
    
    
    cols <- col_palette[findInterval(dens$x, seq(min_val, max_val, length.out = length(col_palette)))]
    for (i in 1:(length(dens$x)-1)) {
      rect(dens$x[i], 0, dens$x[i+1], dens$y[i], col = cols[i], border = NA)
    }
  })
  
  observeEvent(input$otsu_thresh, {
    req(im_mat_reactive()) # Ensure im_mat is ready
    otsu_t <- EBImage::otsu(im_mat_reactive(), range = range(im_mat_reactive(), na.rm = TRUE))
    values$thresh <- otsu_t
  })
  
  # Step 3: Write MIS File
  observeEvent(input$write, {
    
    req(im_mat_reactive()) 
    
    withProgress(message = 'Writing MIS file...', {
      incProgress(0.2, detail = "Preparing data")
      
      req(im_mat_reactive())  
      thresh <- input$thresh  
      im_mat <-  (im_mat_reactive())
      im_thresh <- im_mat>thresh
      print(sum(im_thresh,na.rm=T))
      
      mis_coords_orig <- loaded_data()[[2]]
      mis_list <- loaded_data()[[3]]
      region_names <- loaded_data()[[4]]
      reg_center <- loaded_data()[[5]]
      centr_pos_x <- loaded_data()[[6]]
      centr_pos_y <- loaded_data()[[7]]
      
      len_mis <- length(mis_list[[1]])
      
      top_left_x = Inf
      top_left_y = Inf
      bottom_right_x = -Inf
      bottom_right_y = -Inf
      
      # Iterate over the list
      for (coords in mis_coords_orig) {
        top_left_x = min(top_left_x, coords[1], coords[3])
        top_left_y = min(top_left_y, coords[2], coords[4])
        bottom_right_x = max(bottom_right_x, coords[1], coords[3])
        bottom_right_y = max(bottom_right_y, coords[2], coords[4])
      }
      
      x_range  <- c(top_left_x,bottom_right_x)
      y_range  <- c(top_left_y,bottom_right_y)
      
      im_cimg <- imager::as.cimg(im_thresh)
      
      
      if(input$spot_dilation!=1){
        canvas <- as.matrix(imager::resize(im_cimg, size_x = diff(y_range)/10, size_y = diff(x_range)/10))
        fac_x <- round(diff(x_range)/ncol(im_thresh)*input$spot_dilation/10)
        fac_y <- round(diff(y_range)/nrow(im_thresh)*input$spot_dilation/10)
        
        start.time <- Sys.time()
        k <- mmand::shapeKernel(c(fac_x,fac_y), type="box")
        canvas_dil <- mmand::dilate(canvas,k)
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        time.taken
        
        # Convert EBImage object to imager cimg object
        im_cimg <- imager::as.cimg(canvas_dil)
        
      }else{
        im_cimg <- imager::as.cimg(im_thresh)
      }
      
      # Resize the image using nearest-neighbor interpolation
      canvas_dil2_im <- imager::resize(im_cimg, size_x = diff(y_range), size_y = diff(x_range))
      
      # Convert back to EBImage object
      canvas_dil2 <- base::as.matrix(canvas_dil2_im)
      
      canvas_dil2[is.na(canvas_dil2)] <- 0
      canvas <- canvas_dil2
      
      
      incProgress(0.5, detail = "Processing contours")
      
      if(input$fill_holes){
        canvas <- (EBImage::fillHull(canvas))
      }
      
      
      options("max.contour.segments"=999999)
      
      canvas_conts <- contourLines(canvas,nlevels = 1)
      
      nth_element <- function(vector, starting_position, n) {
        vector[seq(starting_position, length(vector), n)]
      }
      
      x_max <- dim(canvas)[1]
      y_max <- dim(canvas)[2]
      
      elem_lengths <- c()
      canvas_conts_transformed <- canvas_conts
      
      for(u in 1:length(canvas_conts)){
        
        x_new <- round(canvas_conts[[u]]$x*x_max)
        y_new <- round(canvas_conts[[u]]$y*y_max)
        
        elem_lengths <- c(elem_lengths,length(x_new))
        
        index <- 0
        
        if(length(x_new)>=80){
          
          x_new <- nth_element(x_new,1,20)
          y_new <- nth_element(y_new,1,20)
          
          index <- index+1
          
        } else{
          
          x_new <- nth_element(x_new,1,10)
          y_new <- nth_element(y_new,1,10)
          
        }
        # print(c(x_new,y_new))
        
        
        canvas_conts_transformed[[u]]$x <- y_new
        canvas_conts_transformed[[u]]$y <- x_new
        
      }
      
      # Create an empty binary matrix with the same dimensions as canvas
      binary_matrix <- matrix(0, nrow = nrow(canvas), ncol = ncol(canvas))
      
      # Loop through each element in canvas_conts_transformed
      for (contour in canvas_conts_transformed) {
        x_coords <- contour$x
        y_coords <- contour$y
        
        # Set the corresponding positions in the binary matrix to 1
        for (i in 1:length(x_coords)) {
          binary_matrix[y_coords[i], x_coords[i]] <- 1
        }
      }
      
      
      incProgress(0.8, detail = "Finalizing new MIS file")
      
      
      mis_list_new <- mis_list[[1]][1:(len_mis-length(region_names))]
      
      mis_list_new$Method[[1]] <- input$new_method
      
      if(length(region_names)==1){
        
        
        body_base <- (mis_list[[1]][[len_mis]])#-1]])
        
        
      }else{
        
        
        body_base <- (mis_list[[1]][[len_mis-1]])
        
        
      }
      
      
      attributes(body_base)$Type <- "3"
      
      raster_new <- body_base[1]
      method_new <- body_base[[2]]
      name_new <- attributes(body_base)[3]
      upper_new <- body_base[3]
      
      #take a different range, not pos pixels but all pixels.
      x_centr <- centr_pos_x#round(mean(x_range))#centr_pos_x#
      y_centr <- centr_pos_y#round(mean(y_range))#centr_pos_y#r
      
      
      #this only works for one region  so far and is for the distance to the top left corner
      x_shift <- x_range[1]
      y_shift <- y_range[1]
      
      ##new resolution!!!!
      raster_new[[1]] <- list(paste(as.character(c(input$resulting_res,input$resulting_res)),collapse=","))
      method_new[[1]] <- list(input$new_method)
      
      
      dilation_factor_x <- 0#resulting_res
      dilation_factor_y <- 0#resulting_res
      for(reg in 1:length(canvas_conts_transformed)){
        
        body_new <- body_base
        
        reg_dist <- c()
        #find which center was closest to this new cluster:
        for(reg_c in reg_center){
          
          x_spot <- canvas_conts_transformed[[reg]]$x[1]+x_shift
          y_spot <- canvas_conts_transformed[[reg]]$y[1]+y_shift
          
          distance <- sqrt((x_spot-reg_c[1])^2+(y_spot-reg_c[2])^2)
          
          reg_dist <- c(reg_dist,distance)
        }
        
        
        orig_region <- region_names[which.min(reg_dist)]
        
        name_new[[1]] <- paste0(orig_region," ", reg)
        
        
        for(point_idx in 1:length(canvas_conts_transformed[[reg]]$x)){
          
          
          x_point <- canvas_conts_transformed[[reg]]$x[point_idx]
          x_dist <- (x_centr-x_point)/x_centr
          x_point <- as.character(round(x_point+x_shift))#+x_dist*+round(dilation_factor_x/2)))
          
          y_point <- canvas_conts_transformed[[reg]]$y[point_idx]
          y_dist <- (y_centr-y_point)/y_centr
          y_point <- as.character(round(y_point+y_shift))#+y_dist*+round(dilation_factor_y/2)))
          
          # print(paste(x_point,y_point))
          
          upper_new[[1]] <- list(paste(c(x_point,y_point),collapse=","))
          if(point_idx<=2){
            body_new[(2+point_idx)] <- upper_new
          }
          else{
            body_new <- append(body_new, upper_new)
          }
          
        }
        
        
        body_new[1] <- raster_new
        body_new[2] <- method_new
        
        attr_tmp <- attributes(body_new)
        names(attr_tmp) <- NULL
        attr_tmp <- unlist(attr_tmp)
        attributes(body_new) <- attributes(body_base)
        attributes(body_new)[3] <- name_new
        attributes(body_new)$names <- attr_tmp
        
        
        mis_list_new <- append(mis_list_new,list(Area=body_new))
      }
      
      mis_list_out <- mis_list
      mis_list_out[[1]] <- mis_list_new
      print("Preparing to write XML")
      
      mis_out <- as_xml_document(mis_list_out)
      
      xmlContent(mis_out)
      
      incProgress(1, detail = "Complete")
    })
    
  })
  
  output$downloadButtonUI <- renderUI({
    if (!is.null(xmlContent()) && xmlContent() != "") {
      downloadButton(
        "download", 
        "DOWNLOAD Guidance",
        style = "color: #fff; background-color: #FF0000; border-color: #CC0000; font-size: 24px; padding: 20px 32px; border-radius: 8px; display: block; margin: 20px auto; width: 60%;"  # Custom styles for large red button
      )
    }
  })
  
  
  output$download <- downloadHandler(
    filename = function() {
      # Ensure there's at least one file and access its name
      if (!is.null(input$imzmlFile) && nrow(input$imzmlFile) > 0) {
        # Extract the name of the first file, removing the extension
        base_name <- tools::file_path_sans_ext(input$imzmlFile$name[1])
        # Concatenate with the other inputs and append ".mis"
        paste0(base_name, "_", input$resulting_res, "_", input$spot_dilation, ".mis")
      } else {
        # Default filename if no file is selected
        paste0("downloaded_file", ".mis")
      }
    },
    content = function(file) {
      # Assuming xmlContent() returns an XML document object
      xml2::write_xml(xmlContent(), file)  # Directly write the XML object to the file
    },
    contentType = "text/xml"
  )
  
}

# Run the app
shinyApp(ui, server)
