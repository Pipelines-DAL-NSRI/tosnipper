#' Convert files to SNIPPER-analysis ready files
#' @description
#' This package takes in xlsx, csv, or vcf files and their metadata to convert them to a SNIPPER-analysis-ready file. 
#' @param input can be an xlsx, csv, or vcf file. If an xlsx or csv file, the first column should be "Sample" while the succeeding columns are genotypes with marker names as headers.
#' @param reference should be an xlsx or csv file. It should contain three columns 1) Sample, 2) Population, and 3) Superpopulation/Continental region
#' @param target.pop is described by the user. This is the target population they want to classify in SNIPPER.
#' @examples
#' tosnipper("excelfile.xlsx", "reference.csv", target.pop = "Filipino")
#' @export
tosnipper <- function(input, reference, target.pop = population){
   
   # load in input file
   if (tools::file_ext(input) == "csv"){
      input.file <- readr::read_csv(input)
   } else if (tools::file_ext(input) == "xlsx"){
      input.file <- readxl::read_excel(input)
   } else if (tools::file_ext(input) == "vcf"){
      vcf_file <- vcfR::read.vcfR(input)
      vcf_gt <- vcfR::extract.gt(vcf_file, return.alleles = TRUE)
      vcf_cols <- vcfR::getFIX(vcf_file)
      vcf_to_df <- data.frame(vcf_cols, vcf_gt)   
      names(vcf_to_df) <- sub('^X', '', names(vcf_to_df))
      
      drops  <- c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER") # remove unnecessary columns
      vcf_to_df <- vcf_to_df[ ,!(names(vcf_to_df) %in% drops)] 
      
      # transposed (rsID as column headers)
      vcf_to_excel <- data.frame(t(vcf_to_df))
      vcf_to_excel <- vcf_to_excel[-1,]
      Sample <- rownames(vcf_to_excel)
      input.file <- data.frame(Sample, vcf_to_excel)
   } else {
      stop("Not an xlsx, csv, or vcf file.")
   } # TO ADD: VCF FILE AS INPUT
   
   ## correct format
   tosnipper <- lapply(
      input.file, 
      function(x){
         gsub(pattern = "/", replacement = "", x = x, fixed = TRUE)
      }
   )
   
   tosnipper <- as.data.frame(tosnipper)
   
   
   # load in the reference
   if(tools::file_ext(references) == "xlsx"){
      reference <- readxl::read_excel(references) 
   } else if(tools::file_ext(references) == "csv"){
      reference <- readr::read_csv(references)
   } else {
      stop("Not an xlsx or csv file")
   }
   
   if("Sample" %in% names(reference)){
      library(dplyr) #in case there's an error loading the library
      matched <- tosnipper %>% dplyr::left_join(reference, by = "Sample")
      last.col <- as.integer(ncol(matched)) # refer to the total number of columns
      sec.last <- last.col -1
      Superpop <- as.data.frame(matched[,last.col])
      Population <- as.data.frame(matched[,sec.last])
      
      Sample <- as.data.frame(matched$Sample)
      
      #this removes the new col added at the end of the df
      data <- as.data.frame(matched[,2:ncol(matched)-1])
      
      #to ensure that the Sample column is not included in the data df
      drops <- "Sample"
      data <- data[,!(names(data) %in% drops)]
      
      #merge in order
      to_excel <- dplyr::bind_cols(Population, Superpop, Sample, data)
      names(to_excel)[names(to_excel) == "matched[, last.col]"] <- "Superpop"
      names(to_excel)[names(to_excel) == "matched[, sec.last]"] <- "Population"
      names(to_excel)[names(to_excel) == "matched$Sample"] <- "Sample"
      
   } else {
      report::report("No Sample in file")
   }
   
   # split per population?
   tosnpr_split <- split(to_excel, to_excel$Population)
   
   library(purrr)
   tosnpr_split <- tosnpr_split %>% map(`rownames<-`, NULL)
   
   # add the no. to the label
   tosnpr_split <- lapply(
      tosnpr_split, 
      function(x){
         x$No <- rownames(x)
         as.data.frame(x)
         data <- as.data.frame(x[,2:ncol(x)-2])
         Sample <- x[,1]
         x <- dplyr::bind_cols(x$No, Sample, data)
      }
   )
   
   # merges all to a list
   merged <- plyr::ldply(tosnpr_split, data.frame)
   merged <- merged[,-c(1,3)]
   
   
   # Subset the target pop
   target <- merged[merged$Population == target.pop,]
   target$snpr <- "0"
   
   non.target <- merged[merged$Population != target.pop,]
   non.target$snpr <- "1"
   
   merged2 <- dplyr::bind_rows(target, non.target)
   
   # add value
   # count number of samples
   sample.count <- as.integer(nrow(merged2))
   
   # get the total number of pops
   pop.only <- as.data.frame(merged2$Superpop)
   pop.only <- pop.only[!duplicated(pop.only), ]
   pop.only <- as.data.frame(pop.only)
   pop.count <- as.integer(nrow(pop.only))
   
   # to add the columns
   merged3 <- merged2[,-2]
   merged3 <- as.data.frame(merged3)
   names(merged3)[names(merged3) == "...1"] <- sample.count
   names(merged3)[names(merged3) == "Superpop"] <- markers
   names(merged3)[names(merged3) == "Sample"] <- pop.count
   
   # merged3[nrow(merged3)+4,] <- NA #or
   merged3 <- rbind(NA, merged3)
   merged3 <- rbind(NA, merged3)
   merged3 <- rbind(NA, merged3)
   merged3 <- rbind(NA, merged3)
   
   openxlsx::write.xlsx(merged3, "snipper.xlsx")
}
