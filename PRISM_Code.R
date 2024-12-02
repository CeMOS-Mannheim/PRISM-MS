#devtools::install_github("CeMOS-Mannheim/PlaquePicker")

library(dplyr)
# library(tidyverse)
library(msiImporter)
library(xml2)
library(stringr)
library(rayshader)

path <- ".imzML"

mis_path <- ".mis"


tolerance <- 10

resulting_res <- 20

spot_dilation <- 1#.25

m_range <- c(100,900)
# m_range <- c(650,800)

new_method <- "D:\\James\\241016_DAN_Sweeping_20.m"

fill_holes <- T

regions <- msiImporter::getRegions(mis_path,imzmlPath = path)
region_names <- regions[["coordsSummary"]][["ROI_name"]]

msData_cardinal <- Cardinal::readMSIData(path, attach.only=TRUE, mass.range = m_range,resolution = tolerance, units = "ppm")

mis_file <- xml2::read_xml(mis_path)

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


goal_mz <- 281.2485

im_tmp <- Cardinal::image(msData_cardinal,mz=goal_mz,plusminus=0.01)

x_vals <- im_tmp[["facets"]][[1]][[1]][["x"]]
y_vals <- im_tmp[["facets"]][[1]][[1]][["y"]]
im_mat <- t(im_tmp[["facets"]][[1]][[1]][["values"]])
mmand::display(im_mat)
im_mat <- (im_mat - mean(im_mat,na.rm=T)) / sd(im_mat,na.rm=T)


mean(t(im_mat),na.rm=T)+1*sd(im_mat,na.rm=T)
EBImage::otsu(im_mat,range=range(im_mat,na.rm=T))
thresh <- EBImage::otsu(im_mat,range=range(im_mat,na.rm=T))

thresh1 <- 1
thresh2 <- mean(im_mat,na.rm=T)+2*sd(im_mat,na.rm=T)
thresh3 <- mean(im_mat,na.rm=T)+3*sd(im_mat,na.rm=T)
im_mat_clean <- im_mat[!is.na(im_mat)]
# Create the data frame
data <- data.frame(x = 1:length(im_mat_clean), y = sort(im_mat_clean))

# Add a new column to store the color information
data$color <- ifelse(data$y > thresh1, "above", "below")
library(ggplot2)


im_thresh1 <- matrix(im_mat,nrow=nrow(im_mat),ncol=ncol(im_mat))
im_thresh1[is.na(im_mat)] <- 0
im_thresh1[im_mat<=thresh1] <- 0
im_thresh1[im_mat>thresh1] <- 1
im_thresh2 <- matrix(im_mat,nrow=nrow(im_mat),ncol=ncol(im_mat))
im_thresh2[is.na(im_mat)] <- 0
im_thresh2[im_mat<=thresh2] <- 0
im_thresh2[im_mat>thresh2] <- 1
im_thresh3 <- matrix(im_mat,nrow=nrow(im_mat),ncol=ncol(im_mat))
im_thresh3[is.na(im_mat)] <- 0
im_thresh3[im_mat<=thresh3] <- 0
im_thresh3[im_mat>thresh3] <- 1

#show possible thresholds
mmand::display(cbind(im_thresh1,matrix(rep(1,dim(im_thresh1)[1])),im_thresh2,matrix(rep(1,dim(im_thresh1)[1])),im_thresh3))


#choose which one fit best
im_thresh <- im_thresh1
im_save <- im_thresh1

im_thresh <- im_thresh1+im_save
im_thresh[im_thresh>1] <- 1
im_save <-im_thresh1

im_thresh_bin <- im_thresh > 0
im_save_bin <- im_save > 0
overlay <- im_thresh_bin + 2*im_save_bin
combined_img <- cbind(im_thresh_bin, overlay, im_save_bin)
# colors <- c("black", "cyan", "magenta", "white")
colors <- c("white", "black")
text <- colorRampPalette(colors)(4)

# plot_obj <- height_shade(t(combined_img), range = c(0, 3), texture = text)
plot_obj <- height_shade(t(im_thresh1), range = c(0, 1), texture = text)
text <- colorRampPalette(colors)(256)
im_mat_2 <- im_mat
im_mat_2[im_mat<=thresh1] <- 0
plot_obj <- height_shade(t(im_mat_2), range = range(im_mat,na.rm=T), texture = text)

# mmand::display(plot_obj)

sum(im_thresh,na.rm=T)

# Initialize variables to store the top left and bottom right coordinates
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


if(spot_dilation!=1){
  canvas <- as.matrix(imager::resize(im_cimg, size_x = diff(y_range)/10, size_y = diff(x_range)/10))
  fac_x <- round(diff(x_range)/ncol(im_thresh)*spot_dilation/10)
  fac_y <- round(diff(y_range)/nrow(im_thresh)*spot_dilation/10)
  
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
canvas_dil2 <-base::as.matrix(canvas_dil2_im)
# canvas_dil2 <- canvas_dil

dim(canvas_dil2)
# mmand::display(canvas_dil2)

canvas <- canvas_dil2
print(range(canvas_dil2))
print(sum(canvas))
if(fill_holes){
  canvas <- (EBImage::fillHull(canvas))
}

mmand::display(canvas)

options("max.contour.segments"=999999)

canvas_conts <- contourLines(canvas,nlevels = 1)

canvas_conts[[1]]

nth_element <- function(vector, starting_position, n) {
  vector[seq(starting_position, length(vector), n)]
}

x_max <- dim(canvas)[1]
y_max <- dim(canvas)[2]

elem_lengths <- c()
canvas_conts_transformed <- canvas_conts
index <- 0
for(u in 1:length(canvas_conts)){

  x_new <- round(canvas_conts[[u]]$x*x_max)
  y_new <- round(canvas_conts[[u]]$y*y_max)

  elem_lengths <- c(elem_lengths,length(x_new))

  
  if(length(x_new)>=80){
    
    x_new <- nth_element(x_new,1,20)
    y_new <- nth_element(y_new,1,20)
    
    index <- index+1
    
  } else{
    
    x_new <- nth_element(x_new,1,10)
    y_new <- nth_element(y_new,1,10)
    
  }
  # print(c(x_new,y_new))

  # break
  
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

# mmand::display(t(binary_matrix))


print(sum(elem_lengths))

mis_path_new <- paste0(strsplit(mis_path,"[.]")[[1]][1],"_",resulting_res,"_",spot_dilation,".mis")

mis_list_new <- mis_list[[1]][1:(len_mis-length(region_names))]
#

mis_list_new$Method[[1]] <- new_method

if(length(region_names)==1){
  
  
  body_base <- (mis_list[[1]][[len_mis]])#-1]])
  
  
}else{
  
  
  body_base <- (mis_list[[1]][[len_mis-1]])
  
  
}


attributes(body_base)$Type <- "3"

raster_new <- body_base[1]
#add method in the future (body_new[2])
method_new <- body_base[[2]]
name_new <- attributes(body_base)[3]
upper_new <- body_base[3]

#take a different range, not pos pixels but all pixels.
x_centr <- centr_pos_x
y_centr <- centr_pos_y

x_shift <- x_range[1]
y_shift <- y_range[1]




##new resolution!!!!
raster_new[[1]] <- list(paste(as.character(c(resulting_res,resulting_res)),collapse=","))
method_new[[1]] <- list(new_method)

print(length(canvas_conts_transformed))

dilation_factor_x <- resulting_res
dilation_factor_y <- resulting_res

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
    x_point <- as.character(round(x_point+x_shift+x_dist*+round(dilation_factor_x/2)))
    
    y_point <- canvas_conts_transformed[[reg]]$y[point_idx]
    y_dist <- (y_centr-y_point)/y_centr
    y_point <- as.character(round(y_point+y_shift+y_dist*+round(dilation_factor_y/2)))
    

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

mis_out <- as_xml_document(mis_list_out)

write_xml(mis_out,mis_path_new)


