# PRQ (Pre-processing for Relative  Quantification in LC-MS/MS)
#
# Algorithm design and strategy by Dr David Baker (Janssen)
# Written by Dr Steven Kiddle (KCL) - steven.kiddle (at) kcl.ac.uk
#
# To check for updates see http://core.brc.iop.kcl.ac.uk/software/
#
# Copyright (c) Steven Kiddle (2013). This software may be 
# freely used and distributed for non-commercial purposes.
#
# paths should be of the form (mac example) "/Users/username/folder/subfolder_or_file"
#
# Example usage: prq("/Users/stevenkiddle/mass_spec_pipeline/AIBL/Raw Data/","/Users/stevenkiddle/mass_spec_pipeline/AIBL/Results/","/Users/stevenkiddle/mass_spec_pipeline/AIBL/six_plex_lookup.csv")


prq<-function(data_path,results_path,table_path,ref_tag = 126,cutoff = 50,verbose = T){
# Main function in PRQ (Pre-processing for Relative 
# Quantification in LC-MS/MS)

	dir.create(results_path, showWarnings = FALSE)
	
	# Perform median normalisation, relative to reference tag
	prq_norm(data_path,results_path,ref_tag,cutoff,verbose)
	
	# Calculate scores from normalised data
	prq_score(data_path,results_path,ref_tag,cutoff,verbose)
	
	# Load and order subject_table for later
	subject_table <- read.table(table_path,sep=",",header=T,as.is=T)
	subject_table <- subject_table[order(subject_table[,2],subject_table[,3]),]
	
	# Extract protein level data
	prq_protein(data_path,results_path,ref_tag,subject_table,verbose)
	
	# Extract peptide level data
	prq_peptide(data_path,results_path,ref_tag,subject_table,verbose)	
	
	if (verbose){
		print("PRQ completed")
		
	}	
	
}

prq_norm<-function(data_path,results_path,ref_tag,cutoff,verbose = T){
# Normalisation function in PRQ (Pre-processing for Relative 
# Quantification in LC-MS/MS)
#
# Performs median ratio normalisation relative to the 
# reference tag

	if (verbose){
		print("Normalising data")
		print(paste("cutoff",cutoff))
	}
	
	files<-list.files(data_path)
	
	med_ratio<-data.frame()
	prop_zeros<-data.frame()
	num_inst<-data.frame()
	max_intens<-data.frame()

	for (i in files){
	
		tmp_file<-file.path(data_path,i)
	
		tmp_table<-read.csv(tmp_file,header=TRUE,as.is=TRUE)
	
		tmp_table_mr_norm<-tmp_table
	
		ms_runs<-unique(tmp_table[,1])
	
		for (j in ms_runs){
		
			ind<-which(tmp_table[,1]==j)
		
			#if (verbose){print(j)} got annoying
			
			
		
		
		
			for (k in 126:131){
				
				if (k != ref_tag){
					
					tag_str<-paste("X",k,"_cor",sep="")
					
					ref_str<-paste("X",ref_tag,"_cor",sep="")
			
					above_cutoff<-((tmp_table[ind,tag_str]>cutoff)&(tmp_table[ind,ref_str]>cutoff))

			
					tmp_med_ratio<-median(tmp_table[ind[above_cutoff],tag_str]/tmp_table[ind[above_cutoff],ref_str])
			
					med_ratio[j,tag_str]<-tmp_med_ratio
			
					prop_zeros[j,tag_str]<-length(which(tmp_table[ind,tag_str]==0))/length(ind)
					num_inst[j,tag_str]<-length(which(tmp_table[ind,tag_str]>0))
					max_intens[j,tag_str]<-max(tmp_table[ind,tag_str])
						
					tmp_table_mr_norm[ind,tag_str]<-tmp_table[ind,tag_str]/tmp_med_ratio
					
				}
			
			}	
		
		}
	
		tmp <- nchar(results_path)
	
		if (substr(results_path,tmp,tmp) == "\\" | substr(results_path,tmp,tmp) == "/"){
			
			tmp_path<-paste(results_path, "normalised data",sep="")
			
		} else {
			
			tmp_path<-file.path(results_path, "normalised data")
	
		}
		
		dir.create(tmp_path, showWarnings = FALSE)

		tmp_str<-file.path(tmp_path, paste(substr(i,1,nchar(i)-4),"_mr_norm.csv",sep=""))
	
		write.table(tmp_table_mr_norm,file=tmp_str,row.names=FALSE,sep=",")	
	
	}
	
	if (verbose){print("normalisation done")}
	
}

prq_score<-function(data_path,results_path,ref_tag,cutoff,verbose = T){
# Scoring function in PRQ (Pre-processing for Relative 
# Quantification in LC-MS/MS)
#
# Calculates peptide ratios based on normalised data and rolls up to protein level

	#setwd("../TMT normalised")
	# go by path instead

	# Load list of files, each summarising proteins found in a single sixplex
	
	if (verbose){
		print("Scoring peptides and proteins")
	}
	
	tmp <- nchar(results_path)
	
	if (substr(results_path,tmp,tmp) == "\\" | substr(results_path,tmp,tmp) == "/"){
			
		tmp_path<-paste(results_path, "normalised data",sep="")
			
	} else {
			
		tmp_path<-file.path(results_path, "normalised data")
	
	}
	
	list_of_files<-list.files(tmp_path)


	# Loop over files in folder, gathering data from each and processing into output files

	for (l in list_of_files){
		
		if (verbose){
			
			print(l)
			
		}
				
		tmp_str<-file.path(tmp_path, l)
	
		# Current file loaded into temporary table
		tmp_table<-read.csv(tmp_str,header=TRUE,as.is=TRUE)
	
		unq_fractions<-sort(unique(tmp_table[,"fraction"]))

		# peptide level output, still seperated by sixplex.

		peptide_sums<-data.frame(fraction=NA,prot_acc=NA,pep_seq=NA,num_instances=NA,tmt_126=NA,tmt_126_instances=NA,tmt_127=NA,tmt_127_instances=NA,tmt_128=NA,tmt_128_instances=NA,tmt_129=NA,tmt_129_instances=NA,tmt_130=NA,tmt_130_instances=NA,tmt_131=NA,tmt_131_instances=NA)

		# Current row index for peptide level output
		m<-1

		# Protein ratio output, gathered together.

		protein_score<-data.frame(fraction=NA,prot_acc=NA)

		n<-1
	


		for (i in 1:length(unq_fractions)){
	
			ind<-which(tmp_table[,"fraction"]==unq_fractions[i])
			prot_acc<-as.character(tmp_table[ind,"prot_acc"])
			unq_prot_acc<-unique(prot_acc)
	
			for (j in 1:length(unq_prot_acc)){
		
				ind2<-which(tmp_table[ind,"prot_acc"]==unq_prot_acc[j])
				unq_pep_seqs<-unique(as.character(tmp_table[ind[ind2],"pep_seq"]))
		
				for (k in 1:length(unq_pep_seqs)){
					ind3<-which(as.character(tmp_table[ind[ind2],"pep_seq"])==unq_pep_seqs[k])
				
					peptide_sums[m,"fraction"]<-unq_fractions[i]
					peptide_sums[m,"prot_acc"]<-unq_prot_acc[j]
					peptide_sums[m,"pep_seq"]<-as.character(unq_pep_seqs[k])
			
					peptide_sums[m,"num_instances"]<-length(ind3)
			
					peptide_sums[m,"tmt_126"]<-sum(tmp_table[ind[ind2[ind3]],"X126_cor"])
			
					peptide_sums[m,"tmt_126_instances"]<-sum(as.numeric(tmp_table[ind[ind2[ind3]],"X126_cor"]>0))
			
					peptide_sums[m,"tmt_127"]<-sum(tmp_table[ind[ind2[ind3]],"X127_cor"])
			
					peptide_sums[m,"tmt_127_instances"]<-sum(as.numeric(tmp_table[ind[ind2[ind3]],"X127_cor"]>0))
			
					peptide_sums[m,"tmt_128"]<-sum(tmp_table[ind[ind2[ind3]],"X128_cor"])
			
					peptide_sums[m,"tmt_128_instances"]<-sum(as.numeric(tmp_table[ind[ind2[ind3]],"X128_cor"]>0))
			
					peptide_sums[m,"tmt_129"]<-sum(tmp_table[ind[ind2[ind3]],"X129_cor"])
				
					peptide_sums[m,"tmt_129_instances"]<-sum(as.numeric(tmp_table[ind[ind2[ind3]],"X129_cor"]>0))
			
					peptide_sums[m,"tmt_130"]<-sum(tmp_table[ind[ind2[ind3]],"X130_cor"])
			
					peptide_sums[m,"tmt_130_instances"]<-sum(as.numeric(tmp_table[ind[ind2[ind3]],"X130_cor"]>0))
			
					peptide_sums[m,"tmt_131"]<-sum(tmp_table[ind[ind2[ind3]],"X131_cor"])		
			
					peptide_sums[m,"tmt_131_instances"]<-sum(as.numeric(tmp_table[ind[ind2[ind3]],"X131_cor"]>0))
					
					if (ref_tag!=126){peptide_sums[m,paste("X126_",ref_tag,sep="")] <- (peptide_sums[m,"tmt_126"] / peptide_sums[m,paste("tmt_",ref_tag,sep="")])}
			
					if (ref_tag!=127){peptide_sums[m,paste("X127_",ref_tag,sep="")] <- (peptide_sums[m,"tmt_127"] / peptide_sums[m,paste("tmt_",ref_tag,sep="")])}
			
					if (ref_tag!=128){peptide_sums[m,paste("X128_",ref_tag,sep="")] <-  (peptide_sums[m,"tmt_128"] / peptide_sums[m,paste("tmt_",ref_tag,sep="")])}
			
					if (ref_tag!=129){peptide_sums[m,paste("X129_",ref_tag,sep="")] <- (peptide_sums[m,"tmt_129"] / peptide_sums[m,paste("tmt_",ref_tag,sep="")])}
			
					if (ref_tag!=130){peptide_sums[m,paste("X130_",ref_tag,sep="")] <- (peptide_sums[m,"tmt_130"] / peptide_sums[m,paste("tmt_",ref_tag,sep="")])}
			
					if (ref_tag!=131){peptide_sums[m,paste("X131_",ref_tag,sep="")] <- (peptide_sums[m,"tmt_131"] / peptide_sums[m,paste("tmt_",ref_tag,sep="")])}
			
					m<-m+1	
			
				

				}
		
				protein_score[n,"fraction"]<-unq_fractions[i]
				protein_score[n,"prot_acc"]<-unq_prot_acc[j]
		
		
				if (ref_tag!=126){
					
					protein_score[n,paste("mean_126_",ref_tag,sep="")] <- mean(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X126_",ref_tag,sep="")])
					
					protein_score[n,paste("median_126_",ref_tag,sep="")] <- median(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X126_",ref_tag,sep="")])
					
					protein_score[n,paste("sd_126_",ref_tag,sep="")] <- sd(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X126_",ref_tag,sep="")])
					
					print("check2")
					
				}
				
				if (ref_tag!=127){
					
					protein_score[n,paste("mean_127_",ref_tag,sep="")] <- mean(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X127_",ref_tag,sep="")])
					
					protein_score[n,paste("median_127_",ref_tag,sep="")] <- median(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X127_",ref_tag,sep="")])
					
					protein_score[n,paste("sd_126_",ref_tag,sep="")] <- sd(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X127_",ref_tag,sep="")])
					
				}
				
				if (ref_tag!=128){
					
					protein_score[n,paste("mean_128_",ref_tag,sep="")] <- mean(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X128_",ref_tag,sep="")])
					
					protein_score[n,paste("median_128_",ref_tag,sep="")] <- median(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X128_",ref_tag,sep="")])
					
					protein_score[n,paste("sd_128_",ref_tag,sep="")] <- sd(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X128_",ref_tag,sep="")])
					
				}
				
				if (ref_tag!=129){
					
					protein_score[n,paste("mean_129_",ref_tag,sep="")] <- mean(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X129_",ref_tag,sep="")])
					
					protein_score[n,paste("median_129_",ref_tag,sep="")] <- median(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X129_",ref_tag,sep="")])
					
					protein_score[n,paste("sd_129_",ref_tag,sep="")] <- sd(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X129_",ref_tag,sep="")])
					
				}		
				
				if (ref_tag!=130){
					
					protein_score[n,paste("mean_130_",ref_tag,sep="")] <- mean(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X130_",ref_tag,sep="")])
					
					protein_score[n,paste("median_130_",ref_tag,sep="")] <- median(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X130_",ref_tag,sep="")])
					
					protein_score[n,paste("sd_130_",ref_tag,sep="")] <- sd(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X130_",ref_tag,sep="")])
					
				}
				
				if (ref_tag!=131){
					
					protein_score[n,paste("mean_131_",ref_tag,sep="")] <- mean(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X131_",ref_tag,sep="")])
					
					protein_score[n,paste("median_131_",ref_tag,sep="")] <- median(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X131_",ref_tag,sep="")])
					
					protein_score[n,paste("sd_131_",ref_tag,sep="")] <- sd(peptide_sums[(m-1-length(unq_pep_seqs)):(m-1),paste("X131_",ref_tag,sep="")])
					
				}					
		
				n<-n+1
		
			}
		}

		# setwd("../peptide_sums")
		# use path instead
		
		tmp <- nchar(results_path)
	
		if (substr(results_path,tmp,tmp) == "\\" | substr(results_path,tmp,tmp) == "/"){
			
			tmp_path2<-paste(results_path, "protein scores",sep="")
			tmp_path3<-paste(results_path, "peptide sums",sep="")
			
		} else {
			
			tmp_path2<-file.path(results_path, "protein scores")
			tmp_path3<-file.path(results_path, "peptide sums")
	
		}
		
		dir.create(tmp_path2, showWarnings = FALSE)
		dir.create(tmp_path3, showWarnings = FALSE)

		tmp_str<-file.path(tmp_path3, paste(substr(l,1,nchar(l)-4),"_peptide_sums.csv",sep=""))	
	
		write.table(peptide_sums,file=tmp_str,row.names=FALSE,sep=",")
	
		#setwd("../protein_scores")
		# use paths instead
		
		tmp_str<-file.path(tmp_path2, paste(substr(l,1,nchar(l)-4),"_protein_scores.csv",sep=""))
	
		write.table(protein_score,file=tmp_str,row.names=FALSE,sep=",")	
	
		#setwd("../TMT normalised")
		# use path instead	
		
	}
		
}

prq_protein<-function(data_path,results_path,ref_tag,subject_table,verbose = T){
# Protein level data extraction in PRQ (Pre-processing for Relative 
# Quantification in LC-MS/MS)
#
# Arranges protein level data into one file	
	
	
	if (verbose){
		print("Extracting protein data")
	}
	
	# Load list of files, each summarising proteins found in a single sixplex
	
	tmp <- nchar(results_path)
	
	if (substr(results_path,tmp,tmp) == "\\" | substr(results_path,tmp,tmp) == "/"){
			
		tmp_path<-paste(results_path, "protein scores",sep="")
			
	} else {
			
		tmp_path<-file.path(results_path, "protein scores")
	
	}
	
	list_of_files<-list.files(tmp_path)

	# Save the current directory so we can return to it
	#sixplex_folder<-getwd()
	# now just use tmp_path

	# Go back one directory to load lookup table, which shows which samples are loaded on each sixplex
	#setwd("..")
	#six_plex_lookup<-read.table("six_plex_lookup.csv",sep=",",header=T)
	# instead use subject_table

	# Prepare empty results spreadsheet
	SubjectIds<-paste("X",as.character(subject_table[,1]),sep="")
	tmp_nas<-data.frame(SubjectIds,NA)
	tmp_nas2<-t(tmp_nas)
	colnames(tmp_nas2)<-SubjectIds
	tier_pooled<-data.frame(fraction=NA,protein_acc=NA,proportion_sixplex=NA,t(tmp_nas2[-1,]))

	# Reformat most columns of the spreadsheet to expect numeric data
	for (i in 4:(3+length(SubjectIds))){
		tier_pooled[,i]<-as.numeric(tier_pooled[,i])
	}

	# 'm' will act as a 'current row' variable, this will be helpful later
	m<-1

	# Return to the sixplex folder
	#setwd(sixplex_folder)
	# now using paths instead

	# Loop over fractions, this is only designed for gel10 experiments, but is easy to modify
	for (i in 1:10){

		# generate temporary list of proteins for later
		tmp_prot_list<-c()
	
		# generate temporary table to store the counts, i.e. how many sixplexs has each protein in this fraction been seen in
		protein_counts<-data.frame(count=NA,proportion=NA)
	
		# Loop over files, incrementing the counts in protein_counts as appropriate
		for (j in list_of_files){
		
			tmp_str<-file.path(tmp_path, j)
		
			# Load a file
			tmp<-read.table(tmp_str,sep=",",header=T)
		
			# Extract protein list for current fraction
			tmp_prot_list<-as.character(tmp[tmp[,"fraction"]==i,"prot_acc"])
		
			# Loop over the list of proteins in this fraction and file
			for (l in tmp_prot_list){
			
				# If this has been seen before, increment the count, otherwise set to one
				if (!is.na(protein_counts[l,"count"])){
				
					protein_counts[l,"count"]<-protein_counts[l,"count"]+1
				
				} else {
				
					protein_counts[l,"count"]<-1
				
				} 
			
			}
		
		}
	
		# At this stage protein counts is fully populated, although it has an artefactual first row, but this gets dropped so doesn't cause a problem
	
		# Use counts to calculate the proportion of sixplexs each protein has been seen in
		protein_counts[,"proportion"]<-protein_counts[,"count"]/length(list_of_files)
	
		# Use call rate threshold to drop some of the proteins, i.e. those not observed sufficiently often
		current_list<-unique(rownames(protein_counts)[which(protein_counts[,"proportion"]>=0)])
	
		# Use this list to populate first three columns of results spreadsheet
		for (k in current_list){
			
			tier_pooled[m,"fraction"]<-i # First column, put current fraction number
			tier_pooled[m,"protein_acc"]<-k # Second column put uniprot
			tier_pooled[m,"proportion_sixplex"]<-protein_counts[k,"proportion"] # Third column put proportion of sixplexs this is observed in
			
			m<-m+1 # Move to next row
			
			# if (verbose){print(paste(i,j,k))} #Print current progress GOT TOO ANNOYING
		}	
	}

	# To make the next loop easier, we name the rows using a combination of rows 1 and 2
	rownames(tier_pooled)<-paste(tier_pooled[,1],tier_pooled[,2],sep="")

	# Define a current six plex variable, because we loop seperately over files and current_sixplexs the order of list_of_files must be right. This can be improved in future versions	
	current_six_plex<-0
	if (verbose){print("check the numbers to the left and right match")}

	# Now we want to put in the data for each fraction, protein and subject

	# So we loop over each file
	for (j in list_of_files){
	
		current_six_plex<-current_six_plex+1
	
		# Print progress, if numbering of these don't match it can cause subjects data to be jumbled up
		if (verbose){print(paste(current_six_plex,j))}

		tmp_str<-file.path(tmp_path, j)
		
		# Load sixplex protein score data
		tmp<-read.table(tmp_str,sep=",",header=T)
	
		# Loop over fractions (gel10)
		for (i in 1:10){
		
			# Find proteins that pass the call rate AND are in this file
			tier_proteins<-as.character(tier_pooled[tier_pooled["fraction"]==i,2])[as.character(tier_pooled[tier_pooled["fraction"]==i,2]) %in% as.character(tmp[tmp[,"fraction"]==i,"prot_acc"])]
		
			# Then loop over this list
			for (k in tier_proteins){
			
				# Find the index of one of these proteins in current sixplex file
				ind<-which(as.character(tmp[tmp[,"fraction"]==i,"prot_acc"])==k)
			
				# Which subject does each tag correspond to in this sixplex
				six_plex_table<-subject_table[(subject_table[,"six_plex"]==current_six_plex),c(1,3)]
			
			
				# IF REFERENCE TAG IS NOT 126, E.G. 131, you need to delete the 131 tag section, change all "mean_num_126" to "mean_num_131" and uncomment the following:
			
				# 126 tmt
				#
				# Add an X to subject ID, as can't have column name starting with a number
				# tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==126,1],sep="")
				#
				# Put right data in right place
				if (ref_tag!=126){
					
					tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==126,1],sep="")
					
					tier_pooled[paste(as.character(i),k,sep=""),tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("mean_126_",ref_tag,sep="")][ind])
					
				}
			
			
				# 127 tmt
				#
			
				if (ref_tag!=127){
					
					tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==127,1],sep="")
					
					tier_pooled[paste(as.character(i),k,sep=""),tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("mean_127_",ref_tag,sep="")][ind])
					
				}
			

			
				# 128 tmt
			
				if (ref_tag!=128){
					
					tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==128,1],sep="")
					
					tier_pooled[paste(as.character(i),k,sep=""),tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("mean_128_",ref_tag,sep="")][ind])
					
				}			
			

			
				# 129 tmt
			
				if (ref_tag!=129){
					
					tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==129,1],sep="")
					
					tier_pooled[paste(as.character(i),k,sep=""),tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("mean_129_",ref_tag,sep="")][ind])
					
				}
			
		
			
				# 130 tmt
			
				if (ref_tag!=130){
					
					tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==130,1],sep="")
					
					tier_pooled[paste(as.character(i),k,sep=""),tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("mean_130_",ref_tag,sep="")][ind])
					
				}
			
				# 131 tmt
			
				if (ref_tag!=131){
					
					tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==131,1],sep="")
					
					tier_pooled[paste(as.character(i),k,sep=""),tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("mean_131_",ref_tag,sep="")][ind])
					
				}	
	

			}
		
		}	
	
	}	

	# setwd("..") 
	# instead use paths
	
	tmp_str<-file.path(results_path,"protein_scores.csv")

	write.table(tier_pooled,file=tmp_str,sep=",",row.names=F)
	
}

prq_peptide<-function(data_path,results_path,ref_tag,subject_table,verbose = T){
# Peptide level data extraction in PRQ (Pre-processing for Relative 
# Quantification in LC-MS/MS)
#
# Arranges peptide data into one file
	
	if (verbose){
		print("Extracting peptide data",verbose = T)
	}
	
	tmp <- nchar(results_path)
	
	if (substr(results_path,tmp,tmp) == "\\" | substr(results_path,tmp,tmp) == "/"){
			
		tmp_path<-paste(results_path, "peptide sums",sep="")
			
	} else {
			
		tmp_path<-file.path(results_path, "peptide sums")
	
	}
	
	list_of_files<-list.files(tmp_path)

	# Save the current directory so we can return to it
	#sixplex_folder<-getwd()
	# use paths instead

	# Go back one directory to load lookup table, which shows which samples are loaded on each sixplex
	#setwd("..")
	# use paths instead
	#six_plex_lookup<-read.table("six_plex_lookup.csv",sep=",",header=T)
	# use function arguement subject_table

	# Prepare empty results spreadsheet
	SubjectIds<-paste("X",as.character(subject_table[,1]),sep="")
	tmp_nas<-data.frame(SubjectIds,NA)
	tmp_nas2<-t(tmp_nas)
	colnames(tmp_nas2)<-SubjectIds
	tier_pooled<-data.frame(fraction=NA,protein_acc=NA,pep_seq=NA,proportion_sixplex=NA,t(tmp_nas2[-1,]))

	# Reformat most columns of the spreadsheet to expect numeric data
	for (i in 5:(4+length(SubjectIds))){
		tier_pooled[,i]<-as.numeric(tier_pooled[,i])
		#print(i) probably would get annoying
	}

	# 'm' will act as a 'current row' variable, this will be helpful later
	m<-1

	# Return to the sixplex folder
	#setwd(sixplex_folder)
	# use paths instead


	# Loop over fractions, this is only designed for gel10 experiments, but is easy to modify
	for (i in 1:10){
	
		# print(i) probably would get annoying

		# generate temporary list of peptides for later
		tmp_protein_list<-c()
		tmp_peptide_list<-c()
	
		# generate temporary table to store the counts, i.e. how many sixplexs has each peptide in this fraction been seen in
		peptide_counts<-data.frame(protein_id=NA,pep_seq=NA,count=NA,proportion=NA)
	
		n<-1
	
		# Loop over files, incrementing the counts in peptide_counts as appropriate
		for (j in list_of_files){
		
			# print(j) probably would get annoying
		
			tmp_str<-file.path(tmp_path, j)
		
			# Load a file
			tmp<-read.table(tmp_str,sep=",",header=T,as.is=T)
		
			# Extract protein list for current fraction
			tmp_peptide_list<-as.character(tmp[tmp[,"fraction"]==i,"pep_seq"])
		
			tmp_protein_list<-as.character(tmp[tmp[,"fraction"]==i,"prot_acc"])
		
			if (length(tmp_peptide_list)>0){
		
				# Loop over the list of proteins in this fraction and file
				for (l in 1:length(tmp_peptide_list)){
			
					count_ind<-which(peptide_counts[,"pep_seq"]==tmp_peptide_list[l] & peptide_counts[,"protein_id"]==tmp_protein_list[l])
			
					# If this has been seen before, increment the count, otherwise set to one
					if (length(count_ind)>0){
				
						peptide_counts[count_ind,"count"]<-peptide_counts[count_ind,"count"]+1
				
					} else {
				
						peptide_counts[n,"protein_id"]<-tmp_protein_list[l]
				
						peptide_counts[n,"pep_seq"]<-tmp_peptide_list[l]
				
						peptide_counts[n,"count"]<-1
				
						n<-n+1
				
					} 
				
				}
				
			}
		
		}
	
		# At this stage peptide counts is fully populated, although it has an artefactual first row, but this gets dropped so doesn't cause a problem
	
		# Use counts to calculate the proportion of sixplexs each protein has been seen in
		peptide_counts[,"proportion"]<-peptide_counts[,"count"]/length(list_of_files)
	
		# Use call rate threshold to drop some of the proteins, i.e. those not observed sufficiently often
		current_list<-peptide_counts[which(peptide_counts[,"proportion"]>=0),]
	
		# Use this list to populate first three columns of results spreadsheet
		for (k in 1:dim(current_list)[1]){
			
			tier_pooled[m,"fraction"]<-i # First column, put current fraction number
		
			tier_pooled[m,"protein_acc"]<-current_list[k,"protein_id"] # Second column put uniprot
		
			tier_pooled[m,"pep_seq"]<-current_list[k,"pep_seq"] # Second column put peptide_sequence
				
			tier_pooled[m,"proportion_sixplex"]<-current_list[k,"proportion"] # Third column put proportion of sixplexs this is observed in
			
			m<-m+1 # Move to next row
			
			#print(paste(i,j,k)) #Print current progress

		}	
	}

	# To make the next loop easier, we name the rows using a combination of rows 1 and 2
	#rownames(tier_pooled)<-paste(tier_pooled[,1],tier_pooled[,3],sep="")

	# Define a current six plex variable, because we loop seperately over files and current_sixplexs the order of list_of_files must be right. This can be improved in future versions
	current_six_plex<-0
	if (verbose){print("check the numbers to the left and right match")}

	# Now we want to put in the data for each fraction, protein and subject

	# So we loop over each file
	for (j in list_of_files){
	
		current_six_plex<-current_six_plex+1
	
		# Print progress, if numbering of these don't match it can cause subjects data to be jumbled up
		if (verbose){print(paste(current_six_plex,j))}

		tmp_str<-file.path(tmp_path, j)

		# Load sixplex protein score data
		tmp<-read.table(tmp_str,sep=",",header=T)
	
		# Loop over fractions (gel10)
		for (i in 1:10){
		
		
		
			# Find peptides that pass the call rate AND are in this file
		
		
			tmp_ind<-which(paste(tmp[tmp[,"fraction"]==i,"prot_acc"],tmp[tmp[,"fraction"]==i,"pep_seq"],sep="") %in% paste(tier_pooled[tier_pooled[,"fraction"]==i,"protein_acc"],tier_pooled[tier_pooled[,"fraction"]==i,"pep_seq"],sep=""))
		
			if (length(tmp_ind)>0){
		
				# Then loop over this list
				for (k in 1:length(tmp_ind)){
			
					#print(paste(i,k))
			
					# Find the index of one of these proteins in current sixplex file
					#ind<-which(as.character(tmp[tmp[,"fraction"]==i,"pep_seq"])==k)
			
					# Which subject does each tag correspond to in this sixplex
					six_plex_table<-subject_table[(subject_table[,"six_plex"]==current_six_plex),c(1,3)]
			
					tier_ind<-which((tier_pooled[,"protein_acc"]==tmp[tmp[,"fraction"]==i,"prot_acc"][tmp_ind[k]])&(tier_pooled[,"pep_seq"]==tmp[tmp[,"fraction"]==i,"pep_seq"][tmp_ind[k]])&(tier_pooled[,"fraction"]==i))
					
					
			
					# IF REFERENCE TAG IS NOT 126, E.G. 131, you need to delete the 131 tag section, change all "mean_num_126" to "mean_num_131" and uncomment the following:
			
					# 126 tmt
					#
					# Add an X to subject ID, as can't have column name starting with a number
					#tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==126,1],sep="")
			
					# Put right data in right place
					#tier_pooled[tier_ind,tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,"X126_131"][tmp_ind[k]])
			
					if (ref_tag!=126){
					
						tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==126,1],sep="")

						tier_pooled[tier_ind,tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("X126_",ref_tag,sep="")][tmp_ind[k]])	
			
					}
							
					# 127 tmt
					#
			
					if (ref_tag!=127){
					
						tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==127,1],sep="")

						tier_pooled[tier_ind,tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("X127_",ref_tag,sep="")][tmp_ind[k]])	
			
					}			

			
					# 128 tmt
			
					if (ref_tag!=128){
					
						tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==128,1],sep="")

						tier_pooled[tier_ind,tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("X128_",ref_tag,sep="")][tmp_ind[k]])	
			
					}		
			
					# 129 tmt
			
					if (ref_tag!=129){
					
						tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==129,1],sep="")

						tier_pooled[tier_ind,tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("X129_",ref_tag,sep="")][tmp_ind[k]])	
			
					}		
		
			
					# 130 tmt
				
					if (ref_tag!=130){
					
						tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==130,1],sep="")

						tier_pooled[tier_ind,tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("X130_",ref_tag,sep="")][tmp_ind[k]])	
			
					}		
			
					# 131 tmt
			
					if (ref_tag!=131){
					
						tmp_subject_id<-paste("X",six_plex_table[six_plex_table[,2]==131,1],sep="")

						tier_pooled[tier_ind,tmp_subject_id]<-as.numeric(tmp[tmp[,"fraction"]==i,paste("X131_",ref_tag,sep="")][tmp_ind[k]])	
			
					}		
				
				}
	
			}
		
		}	
	
	}	

	#setwd("..")
	# use paths instead

	tmp_str<-file.path(results_path,"peptide_scores.csv")

	write.table(tier_pooled,file=tmp_str,sep=",",row.names=F)
	
}