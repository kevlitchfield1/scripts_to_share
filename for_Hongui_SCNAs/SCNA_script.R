library(copynumber)
library(GenomicRanges)
library(reshape2)
require(data.table)
library("ggplot2")
library(cowplot)

rm(list=ls(all=TRUE))
case_list_for_discovery<-scan("./final_case_list_Aug2020.txt",what=character())
lymphoma.res <- pcf(data=lymphoma,gamma=12,verbose=FALSE)

segs_resp<-read.table("./responder_segs.txt",header=T,sep="\t",stringsAsFactors=F)
segs_resp<-segs_resp[segs_resp$sample %in% case_list_for_discovery,]
count_resp<-length(unique((segs_resp$sample)))
index1<-(segs_resp$cnTotal < round(segs_resp$Ploidy,0))
index2<-(segs_resp$cnTotal > round(segs_resp$Ploidy,0))
segs_resp$mean<-0
segs_resp$mean[index1]=-0.5
segs_resp$mean[index2]=0.5
segs_resp$arm<-"p"
segs_resp2<-segs_resp[,c(1,2,15,3,4,5,14)]

names(segs_resp2)<-names(lymphoma.res)
segs_resp2<-segs_resp2[segs_resp2$chrom!="23",]
segs_resp2<-segs_resp2[segs_resp2$chrom!="24",]

insertSource("./plot_line_resp.R",package="copynumber");plotFreq(segments=segs_resp2,thres.gain=0.2,thres.loss=-0.1)

segs_resp<-read.table("./non_responder_segs.txt",header=T,sep="\t",stringsAsFactors=F)
segs_resp<-segs_resp[segs_resp$sample %in% case_list_for_discovery,]
count_non_resp<-length(unique((segs_resp$sample)))
index1<-(segs_resp$cnTotal < round(segs_resp$Ploidy,0))
index2<-(segs_resp$cnTotal > round(segs_resp$Ploidy,0))
segs_resp$mean<-0
segs_resp$mean[index1]=-0.5
segs_resp$mean[index2]=0.5
segs_resp$arm<-"p"
segs_resp2<-segs_resp[,c(1,2,15,3,4,5,14)]

names(segs_resp2)<-names(lymphoma.res)
segs_resp2<-segs_resp2[segs_resp2$chrom!="23",]
segs_resp2<-segs_resp2[segs_resp2$chrom!="24",]

insertSource("./plot_line_nonresp.R",package="copynumber");plotFreq(segments=segs_resp2,thres.gain=0.2,thres.loss=-0.1)

# Do cytoband analysis
dels_resp<-read.table("./resp_del_freqs.txt",stringsAsFactors=F,header=T)
dels_non_resp<-read.table("./non_resp_del_freqs.txt",stringsAsFactors=F,header=T)
amps_resp<-read.table("./resp_amp_freqs.txt",stringsAsFactors=F,header=T)
amps_non_resp<-read.table("./non_resp_amp_freqs.txt",stringsAsFactors=F,header=T)
bands<-read.table("./cytobands_b37.txt",stringsAsFactors=F,header=T)
#bands$cyto<-NULL;names(bands)[8]<-"cyto"
bands_GR <- with(bands, GRanges(cumulative_chr, IRanges(start=(cumulative_start/1000),end=(cumulative_end/1000),names=paste0("covid:",seq(1:nrow(bands)))),cyto=cyto))

dels_resp$cumulative_start<-round(dels_resp$xleft*1000,0);dels_resp$cumulative_end<-round(dels_resp$xleft*1000,0);dels_resp$cumulative_chr<-"chrAll"
dels_resp_GR <- with(dels_resp, GRanges(cumulative_chr, IRanges(start=(cumulative_start),end=(cumulative_end)),freq.del=freq.del))
matches <- suppressWarnings(subsetByOverlaps(dels_resp_GR,bands_GR))
hits <- findOverlaps(dels_resp_GR,bands_GR,select=c("first"))
idx<-hits
values <- DataFrame(covid=names(bands_GR)[idx], cyto=bands_GR$cyto[idx])

mcols(matches) <- c(mcols(matches), values)
matches_del_resp2<-matches
out_del_resp2<-data.frame(matches)
#out_del_resp<-aggregate(out_del_resp2$freq.del, by=list(out_del_resp2$cyto), FUN=min)

dels_non_resp$cumulative_start<-round(dels_non_resp$xleft*1000,0);dels_non_resp$cumulative_end<-round(dels_non_resp$xleft*1000,0);dels_non_resp$cumulative_chr<-"chrAll"
dels_non_resp_GR <- with(dels_non_resp, GRanges(cumulative_chr, IRanges(start=(cumulative_start),end=(cumulative_end)),freq.del=freq.del))
matches <- suppressWarnings(subsetByOverlaps(dels_non_resp_GR,bands_GR))
hits <- findOverlaps(dels_non_resp_GR,bands_GR,select=c("first"))
idx<-hits
values <- DataFrame(covid=names(bands_GR)[idx], cyto=bands_GR$cyto[idx])
mcols(matches) <- c(mcols(matches), values)
matches_del_non_resp2<-matches
out_del_non_resp2<-data.frame(matches)
#out_del_non_resp<-aggregate(out_del_non_resp2$freq.del, by=list(out_del_non_resp2$cyto), FUN=min)

#find matching breakpoints for deletions and deduplicate data
idx<-nearest(matches_del_resp2,matches_del_non_resp2,ignore.strand=T)
out_del_non_resp3<-out_del_non_resp2[idx,]
out_del_comb<-cbind(out_del_resp2,out_del_non_resp3)
index<-out_del_comb[8]==out_del_comb[16]
out_del_comb<-out_del_comb[index,]
out_del_comb<-out_del_comb[,c(16,6,14,2)]
names(out_del_comb)<-c("cyto","resp_del_freq","non_resp_del_freq","start")
out_del_comb$freq_diff<-abs(out_del_comb$resp_del_freq-out_del_comb$non_resp_del_freq)
out_del_comb<-as.data.table(out_del_comb)
#remove duplicate start positions
out_del_comb<-out_del_comb[out_del_comb[, .I[freq_diff == min(freq_diff)], by=start]$V1]
out_del_comb_max_diff<-data.frame(out_del_comb[out_del_comb[, .I[freq_diff == max(freq_diff)], by=cyto]$V1])
out_del_comb_max_diff$start<-NULL
out_del_comb_max_diff<-out_del_comb_max_diff[!duplicated(out_del_comb_max_diff),]

#amps
amps_resp$cumulative_start<-round(amps_resp$xleft*1000,0);amps_resp$cumulative_end<-round(amps_resp$xleft*1000,0);amps_resp$cumulative_chr<-"chrAll"
amps_resp_GR <- with(amps_resp, GRanges(cumulative_chr, IRanges(start=(cumulative_start),end=(cumulative_end)),freq.amp=freq.amp))
matches <- suppressWarnings(subsetByOverlaps(amps_resp_GR,bands_GR))
hits <- findOverlaps(amps_resp_GR,bands_GR,select=c("first"))
idx<-hits
values <- DataFrame(covid=names(bands_GR)[idx], cyto=bands_GR$cyto[idx])
mcols(matches) <- c(mcols(matches), values)
matches_amps_resp2<-matches
out_amps_resp2<-data.frame(matches)
#out_amps_resp<-aggregate(out_amps_resp2$freq.amp, by=list(out_amps_resp2$cyto), FUN=min)

amps_non_resp$cumulative_start<-round(amps_non_resp$xleft*1000,0);amps_non_resp$cumulative_end<-round(amps_non_resp$xleft*1000,0);amps_non_resp$cumulative_chr<-"chrAll"
amps_non_resp_GR <- with(amps_non_resp, GRanges(cumulative_chr, IRanges(start=(cumulative_start),end=(cumulative_end)),freq.amp=freq.amp))
matches <- suppressWarnings(subsetByOverlaps(amps_non_resp_GR,bands_GR))
hits <- findOverlaps(amps_non_resp_GR,bands_GR,select=c("first"))
idx<-hits
values <- DataFrame(covid=names(bands_GR)[idx], cyto=bands_GR$cyto[idx])
mcols(matches) <- c(mcols(matches), values)
matches_non_amps_resp2<-matches
out_amps_non_resp2<-data.frame(matches)
#out_amps_non_resp<-aggregate(out_amps_non_resp2$freq.amp, by=list(out_amps_non_resp2$cyto), FUN=min)

#find matching breakpoints for amps and deduplicate data
idx<-nearest(matches_amps_resp2,matches_non_amps_resp2,ignore.strand=T)
out_amps_non_resp3<-out_amps_non_resp2[idx,]
out_amps_comb<-cbind(out_amps_resp2,out_amps_non_resp3)
index<-out_amps_comb[8]==out_amps_comb[16]
out_amps_comb<-out_amps_comb[index,]
out_amps_comb<-out_amps_comb[,c(16,6,14,2)]
names(out_amps_comb)<-c("cyto","resp_amp_freq","non_resp_amp_freq","start")
out_amps_comb$freq_diff<-abs(out_amps_comb$resp_amp_freq-out_amps_comb$non_resp_amp_freq)
out_amps_comb<-as.data.table(out_amps_comb)
out_amps_comb<-out_amps_comb[out_amps_comb[, .I[freq_diff == min(freq_diff)], by=start]$V1]
out_amps_comb_max_diff<-data.frame(out_amps_comb[out_amps_comb[, .I[freq_diff == max(freq_diff)], by=cyto]$V1])
out_amps_comb_max_diff$start<-NULL
out_amps_comb_max_diff<-out_amps_comb_max_diff[!duplicated(out_amps_comb_max_diff),]

#merge outputs
merged_dat<-merge(out_del_comb_max_diff,out_amps_comb_max_diff,by.x="cyto",by.y="cyto")
merged_dat$freq_diff.x<-NULL
merged_dat$freq_diff.y<-NULL
names(merged_dat)<-c("cytoband","del_freq_resp","del_freq_non_resp","amp_freq_resp","amp_freq_non_resp")

#do Fishers Exact test per cytoband
merged_dat$del_pval<-1
merged_dat$amp_pval<-1
for(i in 1:nrow(merged_dat)){
	res_del <- fisher.test(matrix(c(count_resp*merged_dat[i,2]/100,count_non_resp*merged_dat[i,3]/100,count_resp-count_resp*merged_dat[i,2]/100,count_non_resp-count_non_resp*merged_dat[i,3]/100),nrow = 2))
	merged_dat[i,6]<-res_del$p.value
	res_amp <- fisher.test(matrix(c(count_resp*merged_dat[i,4]/100,count_non_resp*merged_dat[i,5]/100,count_resp-count_resp*merged_dat[i,4]/100,count_non_resp-count_non_resp*merged_dat[i,5]/100),nrow = 2))
	merged_dat[i,7]<-res_amp$p.value
			}
merged_dat$del_qval<-p.adjust(merged_dat$del_pval, "fdr")
merged_dat$amp_qval<-p.adjust(merged_dat$amp_pval, "fdr")

del_hits<-head(merged_dat[order(merged_dat$del_pval),],10)
amp_hits<-head(merged_dat[order(merged_dat$amp_pval),],10)
del_hits$freq_diff<-(del_hits$del_freq_resp-del_hits$del_freq_non_resp)
amp_hits$freq_diff<-(amp_hits$amp_freq_resp-amp_hits$amp_freq_non_resp)

p1<-ggplot(data=del_hits, aes(x=reorder(cytoband,freq_diff), y=freq_diff)) +geom_bar(stat="identity", position=position_dodge(),fill="darkred")+ theme_minimal()+coord_flip()
p2<-ggplot(data=amp_hits, aes(x=reorder(cytoband,freq_diff), y=freq_diff)) +geom_bar(stat="identity", position=position_dodge(),fill="darkblue")+ theme_minimal()+coord_flip()
plot_grid(p1,p2)
ggsave("./cytoband_freqs.pdf",dpi=600)


del_hits$labels2<-paste0(round(del_hits$del_pval,4),",",round(del_hits$del_qval,2))
amp_hits$labels2<-paste0(round(amp_hits$amp_pval,4),",",round(amp_hits$amp_qval,2))

p3<-ggplot(data=del_hits, aes(x=reorder(labels2,freq_diff), y=freq_diff)) +geom_bar(stat="identity", position=position_dodge(),fill="darkred")+ theme_minimal()+coord_flip()
p4<-ggplot(data=amp_hits, aes(x=reorder(labels2,freq_diff), y=freq_diff)) +geom_bar(stat="identity", position=position_dodge(),fill="darkblue")+ theme_minimal()+coord_flip()
plot_grid(p3,p4)

ggsave("./cytoband_freqs_pvals.pdf",dpi=600)

## Only needed if you want to make a focal plot of a particular cytoband, to see if there is a strong peak over a particular gene
## Example one for 9q34 plot
out_del_to_plot<-data.frame(out_del_comb)
index<-out_del_to_plot[1]=="chr9_q34"
out_del_to_plot<-out_del_to_plot[index,]
#out_del_to_plot$freq_diff[index]<-0
out_del_to_plot$xleft<-out_del_to_plot$start/1000
#out_del_to_plot<-out_del_to_plot[,c(6,5)]
#names(out_del_to_plot)[2]<-"freq.del"
#write.table(out_del_to_plot,"CPI1000_analysis/CN_analysis/9q34_deletion_freqs.txt",quote=F,sep="\t")
out_del_to_plot$chr9_Mb<- out_del_to_plot$xleft-(1539159712/1000000)
out_del_to_plot_m<-melt(out_del_to_plot,measure.vars = c("resp_del_freq","non_resp_del_freq"))

ggplot(out_del_to_plot_m) + geom_line(aes(y = -freq_diff, x = chr9_Mb),color="black",size=0.6,data = out_del_to_plot, stat="identity")+scale_y_continuous()+theme_bw()
ggsave("./focal_plot.pdf",dpi=600)
