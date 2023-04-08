scale_variables <- function(data, mns, sds){
  data = data[,as.character(names(mns))]
  data = sweep(data, 2, mns)
  data = sweep(data, MARGIN=2, FUN="/", sds)
  data[is.na(data)] <- 0
  return(data)
}

add_adjacent_ints <- function(x, prefix){
  factors = paste("AF", seq(1, 5), sep="")
  factor_grid = expand.grid(factors, factors)
  for (i in 2:18){
    pos2 = paste(prefix, i, sep="")
    if (!(paste(pos2, "AF1", sep="_") %in% colnames(x))) { next }
    pos1 = paste(prefix, (i-1), sep="")
    for (j in 1:nrow(factor_grid)){
      ind2 = which(colnames(x)==paste(pos2, factor_grid$Var1[j], sep="_"))
      ind1 = which(colnames(x)==paste(pos1, factor_grid$Var2[j], sep="_"))
      new = x[,ind1] * x[,ind2]
      x = cbind(x, new)
      colnames(x)[ncol(x)] = paste(colnames(x)[ind1], colnames(x)[ind2], sep="_by_")
    }
  }
  return(x)
}

reformat_gene <- function(string, gene, brds){
  string = gsub("\\*", "-", string)
  pref = strsplit(string, gene)[[1]][1]
  info = strsplit(string, gene)[[1]][2]
  fam = strsplit(info, "-")[[1]][1]
  mem = strsplit(info, "-")[[1]][2]
  if (substr(fam, 1, 1)=="0"){
    fam = substr(fam, 2, nchar(fam))
  }
  if (!(is.na(mem))){
    if (substr(mem, 1, 1)=="0"){
      mem = substr(mem, 2, nchar(mem))
    }
  }
  base = paste(paste(gsub("C", "", pref), gene, sep=""), fam, sep="")
  if (!(base %in% brds) & !(is.na(mem))){
    return(paste(paste(gsub("C", "", pref), gene, sep=""), paste(fam, mem, sep="-"), sep=""))
  } else {
    return(base)
  }
}

create_hashmap <- function(key_vector, val_vector){
  h = hash::hash()
  for (i in 1:length(key_vector)){
    h[[key_vector[i]]] = val_vector[i]
  }
  return(h)
}

str_index <- function(str_length, pos, scheme, max_length){
  if (scheme=="LR" & pos<=ceiling(str_length/2)){
    return(pos)
  }
  if (scheme=="LR" & pos>ceiling(str_length/2)){
    diff = max_length - str_length
    if ((pos-diff) > ceiling(str_length/2)){
      return(pos - diff)
    } else {
      return(".")
    }
  }
  if (scheme=="mid"){
    tab = floor((max_length - str_length)/2)
    if (tab>=pos){
      return(".")
    } else {
      if ((pos-tab)<=str_length){
        return(pos - tab)
      } else {
        return(".")
      }
    }
  }
}

get_feat_score <- function(x, amap){
  if (is.na(x)) { return(NA) }
  if (x=="."){
    return(NA)
  }
  sum = 0
  for (i in 1:nchar(x)){
    sum = sum + amap[[substr(x, i, i)]]
  }
  return(sum/nchar(x))
}

featurize_tcrs <- function(data, chain, cdr3_align="mid", cdr_only = TRUE, add_ints52 = TRUE, return_seq_grid=FALSE, do_jgenes = TRUE, restrict_length=TRUE){
  brd_BV = TRBVgrid$gene[!(grepl("-", TRBVgrid$gene))]
  brd_BJ = TRBJgrid$gene[!(grepl("-", TRBJgrid$gene))]
  brd_AV = TRAVgrid$gene[!(grepl("-", TRAVgrid$gene))]
  brd_AJ = TRAJgrid$gene[!(grepl("-", TRAJgrid$gene))]
  brds <- c(brd_BV, brd_BJ, brd_AV, brd_AJ)

  if (chain %in% c("b", "ab")){
    data$TCRB_vgene = as.character(data$TCRB_vgene)
    data$TCRB_jgene = as.character(data$TCRB_jgene)
    data = data[data$TCRB_vgene!="unresolved",]
    data = data[data$TCRB_vgene!="",]
    data$TCRB_vgene = sapply(as.character(data$TCRB_vgene), function(x) reformat_gene(x, gene="V", brds))
    data = data[data$TCRB_vgene %in% TRBVgrid$gene,]
    if (do_jgenes){
      data = data[data$TCRB_jgene!="unresolved",]
      data$TCRB_jgene = sapply(as.character(data$TCRB_jgene), function(x) reformat_gene(x, gene="J", brds))
      data = data[data$TCRB_jgene %in% TRBJgrid$gene,]
    }
    data$TCRB_cdr3aa = as.character(data$TCRB_cdr3aa)
    if (restrict_length){
      nc = sapply(data$TCRB_cdr3aa, function(x) nchar(x))
      data = data[nc >=11 & nc <= 18,]
    }
  }
  if (chain %in% c("a", "ab")){
    data$TCRA_vgene = as.character(data$TCRA_vgene)
    data$TCRA_jgene = as.character(data$TCRA_jgene)
    data = data[data$TCRA_vgene!="unresolved",]
    data = data[data$TCRA_vgene!="",]
    data$TCRA_vgene = sapply(as.character(data$TCRA_vgene), function(x) reformat_gene(x, gene="V", brds))
    data = data[data$TCRA_vgene %in% TRAVgrid$gene,]
    if (do_jgenes){
      data = data[data$TCRA_jgene!="unresolved",]
      data$TCRA_jgene = sapply(as.character(data$TCRA_jgene), function(x) reformat_gene(x, gene="J", brds))
      data = data[data$TCRA_jgene %in% TRAJgrid$gene,]
    }
    data$TCRA_cdr3aa = as.character(data$TCRA_cdr3aa)
    if (restrict_length){
      nc = sapply(data$TCRA_cdr3aa, function(x) nchar(x))
      data = data[nc >=10 & nc <= 17,]
    }
  }
  pos = data
  if (chain %in% c("a", "ab")){
    pos = dplyr::left_join(pos, TRAVgrid, by=c("TCRA_vgene"="gene"))
    pos = dplyr::left_join(pos, TRAJgrid, by=c("TCRA_jgene"="gene"))
  }
  if (chain %in% c("b", "ab")){
    pos = dplyr::left_join(pos, TRBVgrid, by=c("TCRB_vgene"="gene"))
    pos = dplyr::left_join(pos, TRBJgrid, by=c("TCRB_jgene"="gene"))
  }

  if (chain %in% c("a", "ab")){
    nc = sapply(data$TCRB_cdr3aa, function(x) nchar(x))
    for (i in 1:17){
      pos$new = suppressWarnings(sapply(data$TCRA_cdr3aa, function(y) substr(y, str_index(nchar(y), i, cdr3_align, 17), str_index(nchar(y), i, cdr3_align, 17))))
      pos$new[is.na(pos$new)] = "."
      colnames(pos)[ncol(pos)] = paste("TRAcdr3_p", i, sep="")
    }
  }

  if (chain %in% c("b", "ab")){
    for (i in 1:18){
      pos$new = suppressWarnings(sapply(data$TCRB_cdr3aa, function(y) substr(y, str_index(nchar(y), i, cdr3_align, 18), str_index(nchar(y), i, cdr3_align, 18))))
      pos$new[is.na(pos$new)] = "."
      colnames(pos)[ncol(pos)] = paste("TRBcdr3_p", i, sep="")
    }
  }
  if (return_seq_grid){
    ind = which(colnames(pos)=="TRA_p1")
    pos = pos[,ind:ncol(pos)]
    rownames(pos) = as.character(data[,1])
    return(pos)
  }

  loops = data
  segs = vector(mode="character")
  positions = vector(mode="character")
  if (chain %in% c("a", "ab")){
    loops = dplyr::left_join(loops, TRAVref, by=c("TCRA_vgene"="gene"))
    loops = dplyr::left_join(loops, TRAJref, by=c("TCRA_jgene"="TRA_Gene"))
    ind = which(colnames(loops)=="TRA_FR1")
    segs = c(segs, colnames(loops)[ind:ncol(loops)], "TCRA_cdr3aa")
    positions = c(positions, colnames(pos)[grepl("TRA_", colnames(pos))], colnames(pos)[grepl("TRAcdr3", colnames(pos))])
  }
  if (chain %in% c("b", "ab")){
    loops = dplyr::left_join(loops, TRBVref, by=c("TCRB_vgene"="gene"))
    loops = dplyr::left_join(loops, TRBJref, by=c("TCRB_jgene"="TRB_Gene"))
    ind = which(colnames(loops)=="TRB_FR1")
    segs = c(segs, colnames(loops)[ind:ncol(loops)], "TCRB_cdr3aa")
    positions = c(positions, colnames(pos)[grepl("TRB_", colnames(pos))], colnames(pos)[grepl("TRBcdr3", colnames(pos))])
  }

  if (cdr_only){
    segs = segs[!(grepl("_FR", segs))]
    segs = segs[!(grepl("_Jtail", segs))]
    cdr = c("TRA_p27", "TRA_p28", "TRA_p29", "TRA_p30", "TRA_p36",
      "TRA_p37", "TRA_p38", "TRA_p56", "TRA_p57", "TRA_p58", "TRA_p63",
      "TRA_p64", "TRA_p65", "TRB_p27", "TRB_p28", "TRB_p29", "TRB_p36",
      "TRB_p37", "TRB_p38", "TRB_p56", "TRB_p57", "TRB_p58", "TRB_p59",
      "TRB_p63", "TRB_p64","TRB_p65")
    positions = positions[grepl("cdr3", positions) | positions %in% cdr]
  }
  pos = pos[,as.character(positions)]

  prc_fts = data.frame(id = as.character(data[,1]))
  ##we don't calculate G because this is our reference
  aminos <- c("A", "C", "D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

  for (i in 1:length(segs)){
    ind = which(colnames(loops)==segs[i])
    loops[,ind] = as.character(loops[,ind])
    prc_fts$new = sapply(loops[,ind], function(x) nchar(x) - stringr::str_count(x, "\\."))
    colnames(prc_fts)[ncol(prc_fts)]  = paste(segs[i], "length", sep="_")
    for (j in 1:length(aminos)){
      prc_fts$new = sapply(loops[,ind], function(x) stringr::str_count(x, aminos[j])/(nchar(x)-stringr::str_count(x, "\\.")))
      colnames(prc_fts)[ncol(prc_fts)] = paste(segs[i], paste("prc", aminos[j], sep=""), sep="_")
    }
  }
  atch_maps <- list()
  for (i in 2:ncol(atch)){
    atch_maps[[i-1]] <- create_hashmap(c(as.character(atch$symbol), "."), c(as.numeric(as.character(atch[,i])), NA))
  }

  atch_fts = data.frame(id = as.character(data[,1]))

  for (i in 1:ncol(pos)){
    for (j in 1:5){
      atch_fts$new = sapply(pos[,i], function(x) get_feat_score(x, atch_maps[[j]]))
      colnames(atch_fts)[ncol(atch_fts)] = paste(colnames(pos)[i], paste("AF", j, sep=""), sep="_")
    }
  }

  res = cbind(prc_fts, atch_fts[,2:ncol(atch_fts)])
  if (add_ints52){
    if (chain %in% c("a", "ab")){ res = add_adjacent_ints(res, "TRAcdr3_p") }
    if (chain %in% c("b", "ab")){ res = add_adjacent_ints(res, "TRBcdr3_p") }
  }

  return(res)
}


