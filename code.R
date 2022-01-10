blast2rearangements <- function(x, minAlignmentLength = 10, minPercentID = 95, CPUs = 20){
  library(dplyr)
  library(data.table)
  
  # Return an empty tibble if an empty data frame was provided.
  if(nrow(x) == 0) return(tibble())
  
  # Return an empty tibble if no rows remain after filtering alignments.
  x <- subset(x, alignmentLength >= minAlignmentLength & gapopen <= 1 & pident >= minPercentID)
  if(nrow(x) == 0) return(tibble())
  
  # Create a virtual cluster for parallel processing.
  cluster <- parallel::makeCluster(CPUs)
  
  # Subset the blast result data and split on query name to create a list.
  z <- dplyr::select(x,  qname, qstart, qend, sstart, send, evalue) %>% group_split(qname)
  
  # Create a splitting vector to split the data across n CPUs while making 
  # sure that query ids are not split across CPUs.
  z <- rbindlist(
          mapply(function(x, n){ x$n <- n; x }, 
                         z, 
                         dplyr::ntile(1:length(z), CPUs), 
                         SIMPLIFY = FALSE))
  
  
  # Process the BLAST result in parallel.
  # Each worker function will return a data.table object which will be collated and'
  # bound into a single data.table object.
  
  r <- rbindlist(parallel::parLapply(cluster, split(z, z$n), function(b){
  #r <- rbindlist(lapply(split(z, z$n), function(b){
    library(dplyr)
    library(IRanges)
    library(data.table)
    
    # Loop over each query name and build a read model
    rbindlist(lapply(split(b, b$qname), function(b2){
      
      # Bin start and end positions to help mitigate sequencing and aligner error.
      
      # Bin alignment start positions which are shifted to the lower end of 3 NT intervals.
      breaks <- seq(1, max(b2$qstart), by = 3)
      b2$qstart_binned <- as.integer(as.character(cut(b2$qstart, breaks = c(breaks, Inf), include.lowest = TRUE, labels = breaks)))
      
      # Bin alignment end positions which are shifted to the upper end of 3 NT intervals.
      breaks <- seq(min(b2$qend), max(b2$qend), by = 3)
      b2$qend_binned <- as.integer(as.character(cut(b2$qend, breaks = c(breaks, Inf), include.lowest = TRUE, labels = breaks, right = FALSE)))
      
      # Sort BLAST results by query start position and evalue (low to high).
      b2 <- arrange(b2, qstart_binned, evalue)
      
      # Alignment to the negative strand will result in the subject end to come before the start.
      # Switch it back so that they are sequential.
      b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
      b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
      b2$strand  <- ifelse(b2$sstart > b2$send, '-', '+')
      
      # Create IRanges.
      ir <- IRanges(start = b2$qstart_binned, end = b2$qend_binned)
      
      # Name the ranges with the binned query positions followed by the actual subject positions.
      names(ir) <- paste0(b2$qstart, '..', b2$qend, '[', b2$sstart2, b2$strand, b2$send2, ']')
      
      # The ranges are ordered by start position and significance. 
      # 10 -------------           (start)
      # 10 ---                     (don't add)
      #      15 --------           (don't add)
      #                    50 ---  (add)
      #
      # Start with the first range, loop through other ranges and add new 
      # fragments if they do not overlap with the previous fragment while overlooking 2 NT overlaps.
      
      o <- ir[1]
      invisible(lapply(split(ir, 1:length(ir)), function(a){
        if(all(! countOverlaps(o, a, minoverlap = 2) > 0)){
          o <<- c(o, a)
        }
      }))
      
      r <- paste0(unique(names(o)), collapse = ';')
      
      data.table(qname = b2$qname[1], rearrangement = r)
    }))
  }))
  
  parallel::stopCluster(cluster)
  as_tibble(r)
}


r <- blast2rearangements(readRDS('blastResult.rds'))

# Show reads with rearrangements.
r[grepl(';', r$rearrangement),]

