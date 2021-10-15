# This document contains the functions that can be performed 
# on all data tables as well as any associated helper functions


# This function is used by all other functions to get the correct tables using the user's
# data pull for that table.
# takes a table and code to filter it
# returns a table that is filtered

dataPull <- function(table, dp) {
  if (table == "mutations"){
    if (dp$m == 0){
      return(mutations)
    }
    
    else{
      t <- filter(mutations, eval(dp$m))
    }
  }
  
  else if (table == "sequences"){
    if (dp$s == 0){
      return(sequences)
    }
    
    else{
      t <- filter(sequences, eval(dp$s))
    }
  }
  return(t)
  
}

# Finds the number of mutations in each row of table. if duplicates are
# not allowed, then mutations are sorted by id.
# takes a table along with its filter information and wether or not to
# count duplicates more than once
# returns a table with a column added

Mutation_Count <- function(tableName, table, dp, factor, duplicates){
  f <- parse(text=factor) # needed for proper typing in count
  filtered <- dataPull("mutations", dp)

  if(duplicates == "Yes"){
    m <-count(filtered, colname= eval(f), name="mutations") #counts number of muts per row
  }
  
  else{
    pre <- filtered %>% count(colname= eval(f), Mutation_Id) #sorts unique mutations in each country 
    m <- count(pre, colname, name="mutations") # counts num of columns in table
    
  }

  #merges list of columns by factor into existing table
  merged <- merge(table, m , by.x =factor , by.y="colname", 
                  all.x = TRUE, all.y = FALSE)
  
  merged["mutations"] [is.na(merged["mutations"])] <- 0  # needs to work if names are changed
  return(merged)
}



# calculates rate of mutations per unit time or chromosomal distance for each row
# takes table info and a unit of measurement
# returns modified table

Mutation_Rate <- function(tableName, table, dp, factor, measurement, duplicates){
  filtered <- dataPull("mutations", dp)
  filteredSeq <- dataPull("sequences", dp)
  f <- parse(text=factor)

  
  # Gathers number of mutations
  table <- Mutation_Count(tableName, table, dp, factor, duplicates) # number of mutations
 
  
  if(measurement == "Mutations per sequence"){
    m <-count(filteredSeq, colname= eval(f), name="unitNum") # number of sequences

    merged <- merge(table, m , by.x =factor , by.y="colname", 
                    all.x = TRUE, all.y = FALSE)

  }
  
  
  # In the case of chromosomal poosition, if the factor is a mutation type then
  # all sequences must be used. otherwise they musted be sorted
  else if(measurement == "Mutations per 1000 bases"){
    if (factor %in% names(filteredSeq)){
      # THIS LINE REPLACES NA WITH '' TO MATCH DATA TABLE AND ALLOW TAPPLY TO WORK 
      filteredSeq[[factor]][is.na(filteredSeq[[factor]])] <- ''
      
      temp <-tapply(filteredSeq$Alignment_Length, filteredSeq[[factor]], sum, simplify = TRUE)
      temp <- as.vector(temp)
      table$unitNum <- temp /1000
    } 
    else{
      table$unitNum <- sum(filteredSeq$Alignment_Length)/ 1000
      
    }
    merged <- table
  }
  
  else if(measurement == "Mutations per million bases"){
    
    if (factor %in% names(filteredSeq)){
      # THIS LINE REPLACES NA WITH '' TO MATCH DATA TABLE AND ALLOW TAPPLY TO WORK 
      filteredSeq[[factor]][is.na(filteredSeq[[factor]])] <- ''
      
      temp <-tapply(filteredSeq$Alignment_Length, filteredSeq[[factor]], sum, simplify = TRUE)
      temp <- as.vector(temp)
      table$unitNum <- temp /1000000
    } 
    else{
      table$unitNum <- sum(filteredSeq$Alignment_Length)/ 1000000
      
    }
    merged <- table
  }
  
  else if (measurement == "Mutations per week"){
    minWeek <- min(na.omit(filtered$Week))
    maxWeek <- max(na.omit(filtered$Week))
    table$unitNum <- (maxWeek - minWeek) + 1
    merged <- table
  }
  
  else if (measurement == "Mutations per month"){
    minMonth <- min(na.omit(filtered$Month %% 12))  # modulus is used until month attrubute is corrected
    maxMonth <- max(na.omit(filtered$Month))
    table$unitNum <- (maxMonth - minMonth) + 1
    merged <- table
  }
  
  
  merged$Mutation_Rate <- merged$mutations / merged$unitNum
  merged$unitNum <- NULL
  return(merged)
}


# helper for consensus sequence that finds most common mutation at base
ch2 <- function(df,seqlst){
  r <- df[which.max(df$n), ]
  nonMut <- seqlst[[r$colname]] - sum(df$n)
  if (r$n > nonMut){
    return(r)
    
  }
}


# HELPER FUNCTION THAT ADDS CONSESUS BASE TO EMPTY CELLS
consensusFiller <- function(position,col){
  base <- as.character(ref[as.integer(position)])
  col <- replace(col, is.na(col), base)
  return(col)
}


# finds any changes to the consensus sequence for each row, and retuns
# a table of the most common nucleotide for each position in each row
Consensus_Sequence <- function(tableName, table, dp, factor){
  filteredMut <- dataPull("mutations", dp)
  filteredSeq <- dataPull("sequences", dp)
  f <- parse(text=factor)
  t <- count(filteredSeq, colname= eval(f), name="unitNum") # number of sequences
  tlist <- as.list(c(t$unitNum)) # number of sequences as a list to allow ease of use
  names(tlist) <- t$colname

  
  # gets a count of mutation by position, new base and country, splits the table by country,
  # finds most common mutation of each position in the country and keeps it if greater than
  # number of non mutated sequences
  temp <- filteredMut %>% count(colname= eval(f),Position, New, Mutation_Id) 
  tableLst <- split(temp,list(temp$colname, temp$Position), drop = TRUE) #
  table2 <- sapply(tableLst, ch2,tlist)
  table2 <- rbindlist(table2, fill= TRUE)

  
  # constructs new data table ho hold results
  rows <- unique(table2$colname)
  cols <- unique(table2$Position)
  t <- data.frame(matrix(nrow =length(rows), ncol = length(cols)))
  row.names(t) <- rows
  names(t) <- cols
  
  #adds results to table
  for (x in 1:nrow(table2)){
    row <- table2[x,]
    col <- row[["colname"]]
    pos <- as.character(row[["Position"]])
    new <- row[["New"]]
    t[col,pos] <- new
  }

  
  # adds reference base to all unfilled cells
  final <-sapply(colnames(t),function(x) consensusFiller(x,t[,x]))
  final <- as.data.frame(final)
  rownames(final) <- rows
  names(final) <- as.character(names(final))
  
  return(final)
  
}

# works in same way as consensus sequence
Consensus_Protein <- function(tableName, table, dp, factor, protein){
  filteredMut <- dataPull("mutations", dp)
  filteredSeq <- dataPull("sequences", dp)
  f <- parse(text=factor)
  filteredMut <- filter(filteredMut, Protein == protein)
  t <- count(filteredSeq, colname= eval(f), name="unitNum") # number of sequences
  tlist <- as.list(c(t$unitNum))
  names(tlist) <- t$colname
  protTable <- filteredMut %>% count(colname=eval(f),Amino_Acid, Effect)
  pList <- split(protTable, list(protTable$colname, protTable$Amino_Acid), drop = TRUE)
  pTable2 <- sapply(pList, ch2, tlist)
  if( is.character(pTable2[[1]])){
    return(as.data.table(pTable2))
  }
  pTable2 <- rbindlist(pTable2, fill= TRUE)
  if(length(names(pTable2)) == 0){
    show_alert(
      title = "Function Complete!",
      text = tags$h3("All consensus protien sequences are identical to
                        the reference."),
      type = "warning",
      html = TRUE
    )
    return()
  }
  rows <- unique(pTable2$colname)
  cols <- unique(pTable2$Amino_Acid)
  df <- data.frame(matrix(nrow =length(rows), ncol = length(cols)))
  row.names(df) <- rows
  names(df) <- cols
  for (x in 1:nrow(pTable2)){
    row <- pTable2[x,]
    col <- row[["colname"]]
    AA <- as.character(row[["Amino_Acid"]])
    effect <- row[["Effect"]]
    df[col,AA] <- effect
  }
  
}

# finds the number of sequences for row
Sequence_Count <- function(tableName, table, dp, factor){
  f <- parse(text=factor)
  filteredSeq <- dataPull("sequences", dp)
  filteredMut <- dataPull("mutations", dp)
  # get number by sorting sequences by factor if factor is in sequences table
  if (factor %in% names(sequences)){
    temp <- filteredSeq %>% count(colname=eval(f))
  }
  
  # sorts by sequence ID if factor is in mutations
  else{
    temp <- filteredMut %>% count(Id, colname=eval(f))
    temp <- tapply(temp$colname, temp$Id, sum)
  }
  merged <- merge(table, temp, by.x =factor , by.y="colname", 
                  all.x = TRUE, all.y = FALSE)
  return(merged)

}

# takes a characteristic of a mutation (new, type, effect, etc) and for each row
# computes the number of mutations with that criteria
Mutation_Ratio <- function(tableName, t, dp, factor, cat){
  f <- parse(text=factor)
  c <- parse(text=cat)
  filteredSeq <- dataPull("sequences", dp)
  filteredMut <- dataPull("mutations", dp)
  
  # checks if new category is in mutations or sequences and pulls each possibility
  if (cat %in% names(filteredMut)){
    cols <- unique(filteredMut[[cat]])
  }
  else{
    cols <- unique(filteredSeq[[cat]])
  }
  
  
  # add columns for each possibility of the category ( ex. C, T, A, G, -)
  for(name in cols){
    t[[name]] <- 0
  }
  #View(t)
  old <- row.names(t)
  w <- t[[factor]]
  w[is.na(w)] <- ''
  row.names(t) <- w 
  
  nums <- filteredMut %>% count(colname=eval(f), cat=eval(c)) # counts mutations with each factor and category

  
  # adds results to table
  for (x in 1:nrow(nums)){
    row <- nums[x,]
    col <- row[["colname"]]
    cat <- as.character(row[["cat"]])
    count <- row[["n"]]
    t[col, cat] <- count
  }
  row.names(t) <- old
  return(t)
}

# FINDS THE UNIQUE MUTATIONS FOR EACH OBSERVATION IN THE TABLE
Unique_Mutations <- function(tableName, t, dp, factor){
  f <- parse(text=factor)
  filteredSeq <- dataPull("sequences", dp)
  filteredMut <- dataPull("mutations", dp)

  
  # THIS LINE REPLACES NA WITH '' TO MATCH DATA TABLE AND ALLOW TAPPLY TO WORK 
  filteredSeq[[factor]][is.na(filteredSeq[[factor]])] <- ''
  filteredMut[[factor]][is.na(filteredMut[[factor]])] <- ''
  
  temp <- tapply(filteredMut[[factor]],filteredMut$Mutation_Id, unique)
  temp <- as.data.frame(temp)
  temp <- temp %>% rename(colname = temp)
  temp <- rownames_to_column(temp, var= "mut")
  t$Unique_Mutations <- ''
  old <- row.names(t)
  w <- t[[factor]]
  w[is.na(w)] <- ''
  row.names(t) <- w
  for (x in 1:NROW(temp)){
    r <- temp[x,]
    row <- r[["colname"]]
    row <- as.character(row)
    print(paste(row, typeof(row)))
    print(old)
    if(row %in% t[[factor]]){
      print("Hi")
      t[row, "Unique_Mutations"] <- paste0(t[row, "Unique_Mutations"], paste0(r[["mut"]], ", "))
    }
  }
  t$Unique_Mutations[t$Unique_Mutations == ""] <- "None"
  row.names(t) <- old
  return(t)

}

# HELPER FOR BASE 
Base_Helper <- function(df,seqlst){
  r <- df[which.max(df$n), ]
  nonMut <- seqlst[[r$colname]] - sum(df$n)
  if (r$n > nonMut){
    return(r)
    
  }
}

# BROKEN :( BUT SHOULD RETURN PIE CHARTS OF BASES AT EACH POSITION THAT HAS
# AT LEAST ONE MUTATION
Base_Frequency <- function(tableName, t, dp, factor){
  show_alert(
    title = "Nice work!",
    text =tags$div(
      tags$h3("Here are the results"),
      renderPlot({
        filteredMut <- dataPull("mutations", dp)
        filteredSeq <- dataPull("sequences", dp)
        f <- parse(text=factor)
        t <- count(filteredSeq, colname= eval(f), name="unitNum") # number of sequences
        tlist <- as.list(c(t$unitNum)) # number of sequences as a list to allow ease of use
        names(tlist) <- t$colname
        
        # gets a count of mutation by position, new base and country, splits the table by country,
        # finds most common mutation of each position in the country and keeps it if greater than
        # number of non mutated sequences
        temp <- filteredMut %>% count(colname= eval(f),Position, New) 
        # THIS DATA IS ADDED TO TEMP FOR TESTING PURPOSES!!!!!!!
        fake <- data.frame("colname"= c("JAPAN", "JAPAN", "JAPAN", "JAPAN"),
                           "Position" = c(100, 100, 100, 100),
                           "New" = c("A","C","G", "T"),
                           "n" = c(60,20,13,7))
        temp <- rbind(temp, fake)
        temp %>% group_by(colname, Position) %>% mutate(Freq = ifelse(is.na(n), NA, n / sum(n))) %>%
        ungroup %>%
        ggplot(aes("", Freq, fill=factor(New))) + 
          geom_bar(width = 1, stat = "identity") +
          coord_polar("y") +       # Make it a pie chart
          facet_wrap(~colname+Position) + # Break it down into 9 charts
          # Below is just aesthetics
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.grid = element_blank(),
                axis.title = element_blank()) +
          guides(fill = FALSE)
      })
    ),
    type = "success",
    html = TRUE
  )
  
}

# HELPER FUNCTION FINDS THE GENE OF EACH POSITION
genefinder <- function(pos){
  gene = "Non-coding"
  if (pos <= 21555){
    if (pos >= 266){
      gene = "ORF1ab"
    }
  }
  else{
    if (pos <= 25384){
      if (pos >= 21563){
        gene = "S"
      }
    }
    else if (pos <= 26220){
      if (pos >= 25393){
        gene = "ORF3a"
      }
    }
    else if (pos <= 26472){
      if (pos >= 26245){
        gene = "E"
      }
    }
    else if (pos <= 27191){
      if (pos >= 26523){
        gene = "M"
      }
    }
    else if (pos <= 27387){
      if (pos >= 27202){
        gene = "ORF6"
      }
    }
    else if (pos <= 27759){
      if (pos >= 27756){
        gene = "Multiple"
      }
      else if(pos >= 27394){
        gene = "ORF7a"
      }
    }
    else if (pos <= 27887){
      if (pos >= 27759){
        gene = "ORF7b"
      }
    }
    else if (pos <= 28259){
      if (pos >= 27894){
        gene = "ORF8"
      }
    }
    else if (pos <= 29533){
      if (pos >= 28274){
        gene = "N"
      }
    }
    else if (pos <= 29674){
      if (pos >= 29558){
        gene = "ORF10"
      }
    }
  }
  return(gene)
}

# FINDS GENE FOR EACH POSITION (WORK IN PROGRESS)
Gene <- function(tableName, table, dp, factor){
  if (factor == "Position"){
    #print(lapply(table$Position, genefinder))
    lst <-lapply(table$Position, genefinder)
    table$Gene <-unlist(lst, use.names=FALSE)
    #View(table)
  }
  return(table)
}

# STUB FOR IS CODING FUNCTION
is_coding <- function(pos){
  
}

# WORK IN PROGRESS  
Coding <- function(tableName, table, dp, factor){
  vect <- vector()
  for (x in table$Position){
    filteredMut <- dataPull("mutations", dp)
    vals <- tapply(filteredMut$Gene, filteredMut$Position, ifelse(first(n) == "Non-coding", "Non-coding", "Coding"))

    
    
  }
  
}

# FINDS AN ASSOCIATION BETWEEN TWO MUTATIONS ANS IS HELPER FOR OTHER
# ASSOCIATION FUNCTIONS
association_finder <- function(mut1, mut2, table,func){

  
  seqtab1 <- table %>% filter(Mutation_Id == mut1)
  seqs1 <- seqtab1$Seq_Id
  seqtab2 <- table %>% filter(Mutation_Id == mut2)
  seqs2 <- seqtab2$Seq_Id
  if (func == "Intersection proportion"){
    return(length(intersect(seqs1,seqs2)) / length(union(seqs1,seqs2)))
  }
  else if (func == "Probability(A | B)"){
    return(length(intersect(seqs1,seqs2)) / length(seqs2))
  }
  else if (func == "Probability(B | A)"){
    return(length(intersect(seqs1,seqs2)) / length(seqs1))
  }
  else if (func == "First known date"){
    lst <- intersect(seqs1,seqs2)
    combined <- rbind(seqtab1,seqtab2)
    both <- combined %>% filter(Seq_Id %in% lst)
    return(min(both$Date, na.rm=TRUE))
  }
  else if (func == "First known sequence"){
    lst <- intersect(seqs1,seqs2)
    combined <- rbind(seqtab1,seqtab2)
    both <- combined %>% filter(Seq_Id %in% lst)
    if(nrow(both)== 0){
      return(NA)
    }
    else{
      return(both$Seq_Id[which.min(both$Date)])
    }
    
    
  }
  else if (func == "First known country"){
    lst <- intersect(seqs1,seqs2)
    combined <- rbind(seqtab1,seqtab2)
    both <- combined %>% filter(Seq_Id %in% lst)
    if(nrow(both)== 0){
      return(NA)
    }
    else{
      return(both$Country[which.min(both$Date)])
    }
    
  }
  
}

# FINDS THE ASSOCIATION BETWEEN A TABLE AND ONE MUTATION    
Association <- function(tableName, table, dp, factor, mutation, func){
  print(func)
  filteredSeq <- dataPull("sequences", dp)
  filteredMut <- dataPull("mutations", dp)
  #add option to ask if associations should be in filtered pool or from all mutations
  
  # advanced naming protocol would be nice
  lst <- lapply(table$Mutation_Id,FUN=function(x) association_finder(x,mutation,filteredMut,func))
  table[[mutation]] <- do.call("c", lst)
  return(table)
}  

# CREATES MATRIX
Association_Matrix <- function(tableName, table, dp, factor, func){
  size <- nrow(table)
  df <- data.frame(matrix(NA,nrow=size, ncol=1, dimnames= list(table[[factor]],c("Mutation"))))
  filteredSeq <- dataPull("sequences", dp)
  filteredMut <- dataPull("mutations", dp)
  for (x in table[[factor]]){
    for (y in table[[factor]]){
      df[x,y] <- association_finder(x,y,filteredMut,func)
    }
  }
  df$Mutation <- NULL


  return(df)
}  





