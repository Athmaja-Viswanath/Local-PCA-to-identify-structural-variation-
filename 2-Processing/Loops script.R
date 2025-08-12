library(vcfR)

# Folder containing VCF files
vcf_folder <- "1-Input"

# Generate vector of file paths for ch1.vcf to ch19.vcf
vcf_files <- file.path(vcf_folder, paste0("chr", 17:19, ".vcf"))

# Initialize lists to store results
all_sample_names <- list()
all_vcf_pos_list <- list()
all_gt_matrix <- list()

for (vcf_file in vcf_files) {
  cat("Processing:", vcf_file, "\n")
  
  # Read VCF
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  
  # Extract sample names (excluding FORMAT column)
  samples <- colnames(vcf@gt)[-1]
  all_sample_names[[basename(vcf_file)]] <- samples
  
  # Extract chromosome and position
  chroms <- getCHROM(vcf)
  positions <- getPOS(vcf)
  # Create named list by chromosome
  vcf_pos_list <- split(positions, chroms)
  all_vcf_pos_list[[basename(vcf_file)]] <- vcf_pos_list
  
  # extract genotype matrix
  gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  head(gt_matrix)
  
  #recode genotype matrix
  gt_recode <- apply(gt_matrix, 2, function(x) {
    x <- gsub("0\\|0|0/0", "0", x) 
    x <- gsub("0\\|1|1\\|0|0/1|1/0", "1", x)
    x <- gsub("1\\|1|1/1", "2", x)
    x[x %in% c(".", "./.", ".|.")] <- NA
    return(as.numeric(x))
  }) 
  
  #making sure the extracted matrix is numeric
  gt_recode_num <- matrix(
    as.numeric(gt_recode),
    nrow = nrow(gt_matrix),
    ncol = ncol(gt_matrix))
  
  #Adding the correct column names
  colnames(gt_recode_num) <- colnames(gt_matrix)
  
  # Create named list by chromosome
  all_gt_matrix[[basename(vcf_file)]] <- gt_recode_num #list of genotype matrix
  
}

# Check results
str(all_sample_names)
str(all_vcf_pos_list)

#converting to genotype matrix
#Extracting genotype matrix
all_gt_matrix <- list()
for (vcf_file in vcf_files) {
  cat("Processing:", vcf_file, "\n")
  
  # Read VCF
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  
  gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  head(gt_matrix)
  
  #recode genotype matrix
  gt_recode <- apply(gt_matrix, 2, function(x) {
    x <- gsub("0\\|0|0/0", "0", x) 
    x <- gsub("0\\|1|1\\|0|0/1|1/0", "1", x)
    x <- gsub("1\\|1|1/1", "2", x)
    x[x %in% c(".", "./.", ".|.")] <- NA
    return(as.numeric(x))
  }) 
  
  #making sure the extracted matrix is numeric
  gt_recode_num <- matrix(
    as.numeric(gt_recode),
    nrow = nrow(gt_matrix),
    ncol = ncol(gt_matrix))
  
  #Adding the correct column names
  colnames(gt_recode_num) <- colnames(gt_matrix)
  
  # Create named list by chromosome
  all_gt_matrix[[basename(vcf_file)]] <- gt_recode_num #list of genotype matrix
  
  
}
#checking if i got matrices
str(all_gt_matrix)


##sliding windows

#definign a function
make_vcf_windows <- function(vcf_obj, window_size = 100, step_size = 50, min_snps = 10) {
  chroms <- getCHROM(vcf_obj)
  positions <- getPOS(vcf_obj)
  
  # Combine chrom and pos
  snp_df <- data.frame(chrom = chroms, pos = positions)
  windows_list <- list()
  
  for (chr in unique(snp_df$chrom)) {
    chr_pos <- snp_df[snp_df$chrom == chr, "pos"]
    n_snps <- length(chr_pos)
    
    i <- 1
    while (i + window_size - 1 <= n_snps) {
      pos_window <- chr_pos[i:(i + window_size - 1)]
      if (length(pos_window) >= min_snps) {
        windows_list[[length(windows_list) + 1]] <- data.frame(
          chrom = chr,
          start = min(pos_window),
          end = max(pos_window)
        )
      }
      i <- i + step_size
    }
  }
  
  windows_df <- do.call(rbind, windows_list)
  rownames(windows_df) <- NULL
  return(windows_df)
}

#calculating and saving sliding windows for chr17-19

all_regions <- list()

for (vcf_file in vcf_files) {
  cat("Processing:", vcf_file, "\n")
  
  # Read VCF
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  regions <- make_vcf_windows(vcf, window_size = 1000, step_size = 500)
  # Create named list by chromosome
  all_regions[[basename(vcf_file)]] <- regions #list of regions
}

str(all_regions)

#gives a dataframe with chrom, start, end as columns that can be used to "map" windows on the genome


###################
#custom window extracter function
make_custom_winfun <- function(gt_matrix, regions, variant_chroms, variant_positions) {
  # Returns a function that extracts genotype matrix for window n
  winfun <- function(n) {
    if(any(n < 1 | n > nrow(regions))) stop("Window index out of range")
    region <- regions[n, ]
    # Find variants in this chrom and position range
    idx <- which(variant_chroms == region$chrom & 
                   variant_positions >= region$start & 
                   variant_positions <= region$end)
    if(length(idx) == 0) {
      warning(sprintf("Window %d has no variants.", n))
      return(NULL)  # or matrix(nrow=0, ncol=ncol(gt_matrix))
    }
    return( gt_matrix[idx, , drop=FALSE] )
  }
  attr(winfun, "max.n") <- nrow(regions)
  attr(winfun, "region") <- function(n) { regions[n, , drop=FALSE] }
  attr(winfun, "samples") <- colnames(gt_matrix)
  class(winfun) <- c("winfun", "function")
  return(winfun)
}

all_winfuns <- list()

variant_chroms_all = list()
variant_positions_all = list()

for (vcf_file in vcf_files) {
  cat("Processing:", vcf_file, "\n")
  
  # Read VCF
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  variant_chroms <- getCHROM(vcf)
  variant_positions <- getPOS(vcf)
  variant_chroms_all[[basename(vcf_file)]] <- variant_chroms
  variant_positions_all[[basename(vcf_file)]] <- variant_positions #list of regions
}



for (chrom in names(all_gt_matrix)) {
  cat("Processing chromosome:", chrom, "\n")
  
  gt_matrix <- all_gt_matrix[[chrom]]
  regions <- all_regions[[chrom]]
  variant_chroms = variant_chroms_all[[chrom]]
  variant_positions = variant_positions_all[[chrom]]
  
  winfun <- make_custom_winfun(gt_matrix, regions, variant_chroms, variant_positions)
  
  all_winfuns[[chrom]] <- winfun
}

head(all_winfuns)



winfun <- make_custom_winfun(all_gt_matrix$chr17.vcf, all_regions$chr17.vcf, variant_chroms_all$chr17.vcf,
                             variant_positions_all$chr17.vcf)

eigen = eigen_windows(data = winfun, k = 2, win = NULL )   #uses eigen_windows() function from lostruct
windist_allsamples <- pc_dist( eigen, npc=2 )  #uses pc_dist() function from lostruct
fit2d_allsamples <- cmdscale( windist_allsamples, eig=TRUE, k=2 )
xy_coords <- fit2d_allsamples$points
stopifnot(nrow(xy_coords) == nrow(all_regions$chr17.vcf))
mds_df <- cbind(all_regions$chr17.vcf, 
                MDS1 = xy_coords[,1],
                MDS2 = xy_coords[,2])
mds_df$midpoint = (mds_df$start + mds_df$end)/2
ggplot(mds_df, aes(x = midpoint, y = MDS1)) +
  geom_point() +
  #facet_wrap(~ chrom, scales = "free_x") +
  labs(title = "MDS1 vs Genomic Position",
       x = "Genomic Midpoint",
       y = "MDS1") +
  scale_color_manual(values = c("black", "red")) +
  theme_bw()

#creating the MDS1 vs Genomic location plot
p <- ggplot(mds_df, aes(x = midpoint, y = MDS1)) +
  geom_point() +
  
  ggplotly(p) #better to run in console

winfun_18 <- make_custom_winfun(all_gt_matrix$chr18.vcf, all_regions$chr18.vcf, variant_chroms_all$chr18.vcf,
                             variant_positions_all$chr18.vcf)

eigen_18 = eigen_windows(data = winfun, k = 2, win = NULL )   #uses eigen_windows() function from lostruct
windist_allsamples <- pc_dist( eigen, npc=2 )  #uses pc_dist() function from lostruct
fit2d_allsamples <- cmdscale( windist_allsamples, eig=TRUE, k=2 )
xy_coords <- fit2d_allsamples$points
stopifnot(nrow(xy_coords) == nrow(all_regions$chr17.vcf))
mds_df <- cbind(all_regions$chr17.vcf, 
                MDS1 = xy_coords[,1],
                MDS2 = xy_coords[,2])
mds_df$midpoint = (mds_df$start + mds_df$end)/2
ggplot(mds_df, aes(x = midpoint, y = MDS1)) +
  geom_point() +
  #facet_wrap(~ chrom, scales = "free_x") +
  labs(title = "MDS1 vs Genomic Position",
       x = "Genomic Midpoint",
       y = "MDS1") +
  scale_color_manual(values = c("black", "red")) +
  theme_bw()

#creating the MDS1 vs Genomic location plot
p <- ggplot(mds_df, aes(x = midpoint, y = MDS1)) +
  geom_point() +
  
  ggplotly(p) #better to run in console





library(ggplot2)
library(plotly)

chromosomes <- c("chr17.vcf", "chr18.vcf", "chr19.vcf")
results_list <- list()
plot_list <- list() 

for (chrom in chromosomes) {
  cat("Processing:", chrom, "\n")
  
  # Create window function
  winfun <- make_custom_winfun(all_gt_matrix[[chrom]],
                               all_regions[[chrom]],
                               variant_chroms_all[[chrom]],
                               variant_positions_all[[chrom]])
  
  # Run eigen_windows
  eigen_res <- eigen_windows(data = winfun, k = 2, win = NULL)
  
  # Calculate distance matrix using pc_dist
  windist <- pc_dist(eigen_res, npc = 2)
  
  # Run MDS
  fit2d <- cmdscale(windist, eig = TRUE, k = 2)
  
  # Extract MDS coordinates
  xy_coords <- fit2d$points
  
  # Check dimension matches
  stopifnot(nrow(xy_coords) == nrow(all_regions[[chrom]]))
  
  # Create data frame with MDS + genomic info
  mds_df <- cbind(all_regions[[chrom]],
                  MDS1 = xy_coords[,1],
                  MDS2 = xy_coords[,2])
  
  mds_df$midpoint <- (mds_df$start + mds_df$end) / 2
  
  # Save in list for later use/plotting
  results_list[[chrom]] <- mds_df
  
  p <- ggplot(mds_df, aes(x = midpoint, y = MDS1)) +
    geom_point() +
    labs(title = paste("MDS1 vs Genomic Position for", chrom),
         x = "Genomic Midpoint",
         y = "MDS1") +
    theme_bw()
  
  plot_list[[chrom]] <- p  # store plot in list
  
  print(p) 
}

# Now you have all results and plots saved in results_list

# Now you can access a plot for a specific chromosome like this:
plot_list[["chr17.vcf"]]  # shows plot for chr17
plot_list[["chr18.vcf"]]  # shows plot for chr18
plot_list[["chr19.vcf"]]


