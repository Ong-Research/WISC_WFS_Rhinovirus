#! Rscript
# functions.R

getProteinStart<-function(probe.ids) {
  probe.list = strsplit(probe.ids, ";");
  protein.start.list = lapply(
    probe.list,
    function(x) {
      
      return(x[length(x)]);
    }
  );
  return(as.integer(unlist(protein.start.list)));
}

getProteinLabel<-function(probe.ids) {
  
  probe.list = strsplit(probe.ids,";");
  protein.list = lapply(
    probe.list, function(x){
      n = length(x)-1;
      return(paste(x[1:n],sep=";",collapse=";"))
    }
  )
  protein.labels = unlist(protein.list);  
  
  return(protein.labels);
}

getProteinTiling<-function(probe_ids, return.vector=TRUE) {
  Pos = getProteinStart(probe_ids);
  Protein = getProteinLabel(probe_ids);
  tiling.df = aggregate(
    as.integer(Pos), 
    by = list(Protein = Protein), 
    function(l) {
      if (length(l) > 1) {
        l = l[order(l, decreasing = FALSE)]
        return(l[2] - l[1])
      }
      return(-1)
    }
  )
  ans = tiling.df
  colnames(ans)[2] = "Tiling";
  if (return.vector) {
    ans.vec = ans$Tiling;
    names(ans.vec) = ans$Protein;
    ans = ans.vec;
  } 
  return(ans);
}

getEpitopeID<-function(protein, start, stop) {
  return(paste(protein,start,stop, sep="_"))
}

getMetricStats<-function(blocks, probes, metric, label, proteins=getProteinLabel(probes), pos = getProteinStart(probes)) {
  
  #cat("blocks:",length(blocks)," probes:", length(probes)," metric:",length(metric), " proteins:",length(proteins), " pos", length(pos))

  keep = proteins %in% blocks$Protein;
  
  metric.df = data.frame(
    Protein = proteins[keep],
    Pos = pos[keep],
    Metric = metric[keep],
    stringsAsFactors=FALSE
  );
  
  metric_list = split(metric.df, metric.df$Protein);
  
  #cat("getMetricStats - start\n");  
  min_value = rep(NA, nrow(blocks))
  max_value = rep(NA, nrow(blocks))
  mean_value = rep(NA, nrow(blocks))
  
  for (idx in 1:nrow(blocks)) {
    temp.df = metric_list[[blocks$Protein[idx]]];
    metric_values = temp.df$Metric[temp.df$Pos >= blocks$Start[idx] & temp.df$Pos <= blocks$Stop[idx]];
    min_value[idx] = min(metric_values, na.rm=TRUE);
    max_value[idx] = max(metric_values, na.rm=TRUE);
    mean_value[idx] = mean(metric_values, na.rm=TRUE);    
  }
  
  ans = data.frame(
    min_value = min_value,
    max_value = max_value,
    mean_value = mean_value
  );
  colnames(ans) = c(
    paste0("Min.",label),
    paste0("Max.",label),
    paste0("Mean.",label)
  )
  return(ans);
}

smoothProbeMat<-function(probe_mat, 
                         probe_meta, 
                         dist_method = "aa", 
                         w=2, weighted=FALSE, 
                         weighted.alpha=0.5, 
                         debug=FALSE) {
  smoothing_mat = createProbeSmoothMat(probe_meta=probe_meta, dist_method=dist_method, w=w, weighted=weighted, weighted.alpha = weighted.alpha);
  return(smoothProbeMatInternal(probe_mat = probe_mat, smoothing_mat = smoothing_mat, debug=debug));
}

smoothProbeMatInternal<-function(probe_mat, smoothing_mat, eps=1e-6, debug=FALSE) {
  probe_mat2 = as.matrix(probe_mat);
  NAs = is.na(probe_mat2);
  
  smoothing_mat2 = smoothing_mat[rownames(probe_mat2),rownames(probe_mat2)];
  if (debug) {print(dim(smoothing_mat2))}
  if (sum(NAs) == 0) {
    probe_smooth_mat = as.matrix(
      smoothing_mat2 %*% probe_mat2
    );
  } else {
    message("smoothProbeMat, handling NAs")
    probe_smooth_mat = probe_mat2;
    for (col_idx in 1:ncol(probe_mat2)) {
      values = probe_mat2[,col_idx];
      svalues = NULL;
      vnas = is.na(values);
      smoothing_mat3 = smoothing_mat2;
      #Set the probes with the missing values to zero, then renormalize the weights in the matrix.
      #Then multiply the new smoothing matrix by the vector of with missing values set to 0.
      #Then set missing value back into result. Copy results back to result matrix.
      if (sum(vnas) > 0) {
        if(debug){message("col_idx:",col_idx, " ", sum(vnas));}
        smoothing_mat3[,vnas] = 0;
        if (debug){cat(".")}
        smoothing_mat3[vnas,] = 0;
        if (debug){cat(".")}

        rm = Matrix::rowSums(smoothing_mat3);
        if (debug){cat(".")}
        rindices = which(rm >= eps & rm <= 1-eps) #Only modify those that need to be modified. 0 - prevent division by zero, 1 - already normalized. Use eps to prevent issues with exact comparison with real numbers
        if (debug){cat(length(rindices));}
        if (debug){cat(".")}
        smoothing_mat3[rindices,] = smoothing_mat3[rindices,] / rm[rindices]

        if (debug){cat(".")}
        values[vnas] = 0;
        if (debug){cat(".")}
        svalues = smoothing_mat3 %*% values;
        if (debug){cat(".")}
        svalues[vnas && rm < eps] = NA; #These are the values those *should* be NA since the smoothing weight values are all zero.
        if (debug){cat(".\n");}
      } else {
        svalues = smoothing_mat3 %*% values;
      }
      probe_smooth_mat[,col_idx] = svalues[,1]
    }
  }
  return(probe_smooth_mat);
  
}

createProbeSmoothMat<-function(probe_meta, 
                               file_path, 
                               dist_method="aa", 
                               w=2, 
                               update=100, 
                               weighted=FALSE, 
                               weighted.alpha=0.5) {
  
  umeta = unique(probe_meta[,c("PROBE_ID","SEQ_ID","POSITION")])
  
  if (w <=0) {
    probe_probe_mat = Diagonal(nrow(umeta), x=1);
    rownames(probe_probe_mat) = umeta$PROBE_ID;
    colnames(probe_probe_mat) = umeta$PROBE_ID;
    
  } else {
    lr = w/2;
    
    split_list = split(umeta, umeta$SEQ_ID);
    
    probe_probe_mat_list = lapply(
      split_list,
      function(l) {
        sub_meta = l;
        sub_meta = sub_meta[order(sub_meta$POSITION, decreasing=FALSE),];
        sub_probe_mat = Diagonal(nrow(sub_meta), x=0);
        rownames(sub_probe_mat) = sub_meta$PROBE_ID;
        colnames(sub_probe_mat) = sub_meta$PROBE_ID;
        for (idx1 in 1:nrow(sub_meta)) {
          probe = sub_meta$PROBE_ID[idx1];
          if (dist_method == "aa") {
            pos = sub_meta$POSITION[idx1];
            pos_left = pos - lr;
            pos_right = pos + lr;
            probes = sub_meta$PROBE_ID[sub_meta$POSITION >= pos_left & sub_meta$POSITION <= pos_right];
            weights = rep(1.0,length(probes))
            
            if (weighted) {
              dists = abs(getProteinStart(probes) - pos);
              
              weights = weighted.alpha ^ dists;
            }
          } else {
            if (dist_method == "peptide") {
              idx_left = max(1,idx1-lr);
              idx_right = min(idx1+lr, nrow(sub_meta))
              probes = sub_meta$PROBE_ID[idx_left:idx_right];
            } else {
              stop("Unknown distance method:",dist_method);  
            }
          }
          sub_probe_mat[probe,probes] = weights / sum(weights);
        }
        return(sub_probe_mat);
      }
    )
    cat("Creating final matrix\n");
    probe_probe_mat = bdiag(probe_probe_mat_list);
    cat("Assigning column/row names\n");
    cols = unlist(lapply(probe_probe_mat_list, function(l){return(colnames(l))}))
    colnames(probe_probe_mat) = cols;
    rownames(probe_probe_mat) = cols;
  }
  
  return(probe_probe_mat);
}

makeEpitopeCalls<-function(probe_sample_padj, 
                           probe_cutoff = 0.05, 
                           probes = rownames(probe_sample_padj),
                           proteins = getProteinLabel(probes),
                           pos = getProteinStart(probes),
                           protein_tiling = getProteinTiling(probes),
                           one_hit_filter = TRUE) {
  
  message("makeEpitopeCalls - start");
  sample_epitopes = list();
  
  all_epitopes = NULL;
  message("Epitopes: Find sample epitopes")
  for (col_idx in 1:ncol(probe_sample_padj)) {
    sample_name = colnames(probe_sample_padj)[col_idx];
    probes_f = rownames(probe_sample_padj)[probe_sample_padj[,col_idx] < probe_cutoff];
    if (length(probes) > 0) {
      epitopes = findBlocksProbeT(probes_f,protein_tiling = protein_tiling);
      epitopes$Epitope_ID = getEpitopeID(epitopes$Protein,epitopes$Start,epitopes$Stop);
      all_epitopes = rbind(all_epitopes, epitopes);
      
      epitopes$Sample = sample_name;
      epitopes_padj_metrics = getMetricStats(
        blocks = epitopes, 
        metric = probe_sample_padj[,col_idx],
        label = "padj",
        probes = probes,
        proteins = proteins,
        pos = pos
      )
      epitopes$Epitope.padj = epitopes_padj_metrics$Max.padj;
      sample_epitopes[[sample_name]] = epitopes;
    }
  }
  
  #Find the unique set of epitopes.
  message("Epitopes: Get unique epitopes and rescore");
  all_epitopes = unique(all_epitopes);
  rownames(all_epitopes) = all_epitopes$EpitopeID;
  #For every epitope region/sample, calculate the maxFDR.
  epitope_fdrs = matrix(1, nrow=nrow(all_epitopes), ncol=ncol(probe_sample_padj));
  rownames(epitope_fdrs) = all_epitopes$Epitope_ID;
  colnames(epitope_fdrs) = colnames(probe_sample_padj)
  
  for (col_idx in 1:ncol(probe_sample_padj)) {
    epitopes_padj_metrics = getMetricStats(
      blocks = all_epitopes, 
      metric = probe_sample_padj[,col_idx],
      label = "padj",
      probes = probes,
      proteins = proteins,
      pos = pos
    )
    epitope_fdrs[,col_idx] = epitopes_padj_metrics$Max.padj;
  }
  
  message("Epitopes: K of N");  
  #Now estimate the K of N by first finding the min fdr per K.
  minFDRs = protStat_calcMinFDR(as.matrix(epitope_fdrs), additional_stats=FALSE);
  
  k_of_n = minFDRs;
  colnames(k_of_n) = paste0("K",1:ncol(minFDRs),".padj");
  k_of_n = cbind(rowSums(minFDRs < probe_cutoff),k_of_n);
  colnames(k_of_n)[1] = "K";
  
  
  
  if (one_hit_filter) {
    message("Epitopes: applying one-hit filter");
    #Find all epitopes that have 1 probe and only 1 sample call.  Remove those epitopes from the results
    k1_epitopes = rownames(k_of_n)[k_of_n$K <= 1]
    n1_epitopes = all_epitopes$Epitope_ID[all_epitopes$Number.Of.Probes == 1]
    onehit_epitopes = intersect(k1_epitopes,n1_epitopes);
    message("Removing ",length(onehit_epitopes)," epitopes from ",length(k1_epitopes), " k1 epitopes");
    for (idx in 1:length(sample_epitopes)) {
      sample_epitopes[[idx]] = sample_epitopes[[idx]][!(sample_epitopes[[idx]]$Epitope_ID %in% onehit_epitopes),];
    }
    epitope_fdrs = epitope_fdrs[!(rownames(epitope_fdrs) %in% onehit_epitopes),]
    minFDRs = minFDRs[!(rownames(minFDRs) %in% onehit_epitopes),];
    all_epitopes = all_epitopes[!(all_epitopes$Epitope_ID %in% onehit_epitopes),]
    k_of_n = k_of_n[!(rownames(k_of_n) %in% onehit_epitopes),]
    
  }
  message("Epitopes: Build list");
  ans = list();
  ans$sample_epitopes = sample_epitopes;
  ans$sample_epitopes_fdr = epitope_fdrs;
  ans$sample_epitopes_minfdr = minFDRs;
  ans$all_epitopes = all_epitopes;
  ans$k_of_n_epitopes = k_of_n;
  
  message("makeEpitopeCalls - Done ", "#epitopes:",length(all_epitopes));
  return(ans);
  
}

findBlocksProbeT<-function(
    probes,
    protein_tiling, 
    proteins = getProteinLabel(probes), 
    starts = getProteinStart(probes)) {
  
  protein.df = data.frame(
    Protein = proteins,
    Pos = starts,
    stringsAsFactors=FALSE
  );
  
  protein_list = split(protein.df, protein.df$Protein)
  
  ans_list = lapply(protein_list, findBlocksT, protein_tiling = protein_tiling);
  ans_dt = data.table::rbindlist(ans_list);
  ans_df = as.data.frame(ans_dt, stringsAsFactors=FALSE);
  
  ans_df$ProbeSet.AA.Span = ans_df$Stop - ans_df$Start + 1;
  rownames(ans_df) = getEpitopeID(ans_df$Protein, ans_df$Start, ans_df$Stop);
  return(ans_df);
  
}

findBlocksT<-function(protein.df, protein_tiling) {
  
  if (nrow(protein.df)  == 1) {
    return(
      data.frame(
        Protein = protein.df$Protein,
        Start = protein.df$Pos,
        Stop = protein.df$Pos,
        Number.Of.Probes = 1,
        stringsAsFactors = FALSE
      )
    );
  }
  tiling = protein_tiling[protein.df$Protein[1]];
  #cat("protein:",protein.df$Protein[1]," tiling:",tiling,"\n");
  protein.df = protein.df[order(protein.df$Pos, decreasing=FALSE),];
  ans.df = NULL;
  start_idx=1;
  start_pos = protein.df$Pos[start_idx];
  for (idx in 2:nrow(protein.df)) {
    current_idx = idx;
    current_pos = protein.df$Pos[current_idx];
    prev_idx = idx-1;
    prev_pos = protein.df$Pos[prev_idx];
    d = current_pos - prev_pos;
    
    if (d > tiling) {
      #Block is start_pos to prev_pos
      #cat("Add ",start_idx," ",prev_idx,"\n");
      #print(pvalues)
      ans.df = rbind(
        ans.df,
        data.frame(
          Protein=protein.df$Protein[1],
          Start = start_pos,
          Stop = prev_pos,
          Number.Of.Probes = length(seq(from=start_pos,to=prev_pos,by=tiling)),
          stringsAsFactors = FALSE
        )
      )
      #Update start_pos
      start_idx = current_idx;
      start_pos = current_pos;
    }
  }
  if (start_idx != nrow(protein.df)) {
    #Add in last block
    ans.df = rbind(
      ans.df,
      data.frame(
        Protein=protein.df$Protein[1],
        Start = start_pos,
        Stop = protein.df$Pos[nrow(protein.df)],
        Number.Of.Probes = length(seq(from=start_pos, to=protein.df$Pos[nrow(protein.df)], by=tiling)),
        stringsAsFactors = FALSE
      )
    )
  } else {
    #Last probe is by itself
    ans.df = rbind(ans.df, data.frame(
      Protein = protein.df$Protein[1],
      Start = protein.df$Pos[nrow(protein.df)],
      Stop = protein.df$Pos[nrow(protein.df)],
      Number.Of.Probes = 1,
      stringsAsFactors=FALSE
      
    ))
  }
  
  return(ans.df); 
}

protStat_calcMinFDR<-function(fdrs, additional_stats=TRUE) {
  fdrs2 = as.data.frame(fdrs,stringsAsFactors=FALSE);
  ncols = ncol(fdrs);
  
  fdrs2 = t(apply(fdrs2, 1, function(l) { return(l[order(l,decreasing=FALSE)])} ))
  fdrs2 = as.data.frame(fdrs2, stringsAsFactors=FALSE)
  colnames(fdrs2) = paste0("n",1:ncol(fdrs2))
  if (additional_stats) {
    fdrs2$meanFDR = rowMeans(fdrs);
    fdrs2$medFDR = matrixStats::rowMedians(fdrs);
  }
  n_col = paste0("n",ncol(fdrs))
  
  fdrs2 = fdrs2[order(fdrs2[,n_col],decreasing=FALSE),]
  return(fdrs2);
  
}

# END
