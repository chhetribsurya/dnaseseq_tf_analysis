library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)

### Stacked barplot for merged site, that is merging all the features within 100bp regions:
input_file <- "~/Dropbox/encode_3/tf_hotspots_dnase/files/merged_d_100/final_dnase_tf_hotspots_single_barplot_data.bed"
output_dir <-  "~/Dropbox/encode_3/tf_hotspots_dnase/plots"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)
final_hotspots_df <- tf_hotspot_df %>% select(state_anno, tf_counts, total_sites) 
final_hotspots_df$log2_site_count <-  log2(final_hotspots_df$total_sites)+1

final_hotspots_df %>% head(10)
prom_col <- colorRampPalette(c("red"))
strong_enh_col <- colorRampPalette(c("orange3"))
weak_enh_col <- colorRampPalette(c("yellow3"))

custom_col <- c(prom_col(1), strong_enh_col(1), weak_enh_col(1))

output_file_name <- file.path(output_dir, "Tf_hotspots_regulatory_region_dist_for_dnase.pdf")				
pdf(output_file_name)

barplot <- ggplot(final_hotspots_df, aes(x=tf_counts, y=log2_site_count, fill=state_anno)) + 
	geom_bar(stat="identity") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("(Log2+1) site counts") + 
    theme_bw() + 
    ggtitle("TF Hotspots distribution across Dnase regulatory sites") + 
    theme(
    axis.text.y = element_text(size=8, face="bold" ),
	plot.title=element_text(size=12, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + scale_x_reverse() + coord_flip() +
	scale_fill_manual(values=custom_col )+
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() 

print(barplot)
dev.off()


##############################
# Typical joint plots:

### For merged site, that is merging all the features within 100bp regions:
input_file <- "~/Dropbox/encode_3/tf_hotspots_dnase/files/merged_d_100/final_tf_hotspots_double_panels_barplot_data.bed"
output_dir <-  "~/Dropbox/encode_3/tf_hotspots_dnase/plots"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)


tf_hotspot_df <- tf_hotspot_df %>% select(tf_counts, total_sites, tf_counts_per_kb_segment)
tf_hotspot_df$log2_total_sites <-  log2(tf_hotspot_df$total_sites) + 1
tf_hotspot_df_1 <-  tf_hotspot_df %>% select(tf_counts, total_sites) %>% mutate(Annotation="Site Counts")
tf_hotspot_df_2 <-  tf_hotspot_df %>% select(tf_counts, tf_counts_per_kb_segment) %>% mutate(Annotation="TF counts per kb")

### have same colnames to make rbind compatible; else rbind gets upset with different colnames:
names(tf_hotspot_df_1) <- c("tf_counts", "value", "Annotation")
names(tf_hotspot_df_2) <- c("tf_counts", "value", "Annotation")
combined_hotspots_df <- rbind(tf_hotspot_df_1, tf_hotspot_df_2)

# change the order of factor in TF Category:
# meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))

output_file_name <- file.path(output_dir, "Tf_hotspots_regulatory_region_dist_for_dnase_1.pdf")				
pdf(output_file_name)

#test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(combined_hotspots_df, aes(x=tf_counts, y=value, fill=Annotation)) + 
	geom_bar(stat="identity") + 
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("Site counts                            TF counts per kb segment") + 
    theme_bw() + 
    ggtitle("TF Hotspots distribution across Dnase regulatory sites") + 
    theme(
    axis.text.y = element_text(size=3, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + scale_x_reverse() +  
	coord_flip() + 
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() +
	facet_wrap(~Annotation)

print(barplot)

dev.off()


### The following codes are not really being used:
################################
################################



output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_dist.pdf")       
pdf(output_file_name, width=6, height=6)

this_plot <- ggplot(tf_hotspot_df, aes(x=as.numeric(tf_hotspot_df$tf_counts))) + 
			geom_histogram(aes(fill = ..count..)) +
			scale_fill_gradient("Count", low = "green", high = "red")
this_plot
dev.off()


#hist(tf_hotspot_df$tf_counts)

input_file <- "~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_hotspots_prom_enh/final_tf_hotspots_sites_distribution_sorted.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist.pdf")       
pdf(output_file_name, width=6, height=6)

this_plot <- ggplot(tf_hotspot_df, aes(x=as.numeric(tf_hotspot_df$tf_counts))) + 
			geom_histogram(aes(fill = ..count..)) +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			xlab("HepG2 sites count") +
			ylab("TF counts")
this_plot
dev.off()


### Barplot for tf hotspots promoter and enhancer associated:
input_file <- "~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_hotspots_prom_enh/final_tf_hotspots_barplot_data.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_25.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts<=25]
this_plot_1 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="black", size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_1
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_25_75.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=26 & tf_hotspot_df$tf_counts<=75]
this_plot_2 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="black", size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_2
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_75_125.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=76 & tf_hotspot_df$tf_counts<=125]
this_plot_3 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="black", size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_3
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_125_max.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=126]
this_plot_4 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="black", size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_4
dev.off()



#require(gridExtra)
#plot_grid(this_plot_4, this_plot_3, labels=c("A", "B"), ncol = 2, nrow = 1)
#grid.arrange(this_plot_4, this_plot_3, ncol=2)


### merged data -d 1000bp, Barplot for tf hotspots promoter and enhancer associated merged within 1000bp:
input_file <- "~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_hotspots_prom_enh/merged_d_1000/final_tf_hotspots_barplot_data.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_25_merged_1000bp.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts<=25]
this_plot_1 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="black", size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_1
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_25_75_merged_1000bp.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=26 & tf_hotspot_df$tf_counts<=75]
this_plot_2 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="black", size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_2
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_75_125_merged_1000bp.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=76 & tf_hotspot_df$tf_counts<=125]
this_plot_3 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="black", size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_3
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_125_max_merged_1000bp.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=126]
this_plot_4 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="black", size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_4
dev.off()



#############################
#############################
#############################

# df$tf_counts <- as.factor(df$tf_counts)
# df$total_sites <- as.numeric(as.character(df$total_sites))
# summary(df)

""" Final data for 208 unique TFs """

input_file <- "~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_hotspots_prom_enh/merged_d_100/final_tf_hotspots_barplot_data.bed"
tf_hotspot_df_read <- fread(input_file, sep="\t", header= TRUE)
df_1 <- tf_hotspot_df_read[,c(2,3)]
df_2 <-  data.frame("tf_counts"=c(0), "total_sites"=c(225103))

new_df <- rbind(df_1, df_2)
tf_hotspot_df <- new_df[order(new_df$tf_counts)]


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_25_merged_100bp.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts<=25]
this_plot_1 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_y_continuous(expand=c(0,0), limits=c(0,250000)) +
			#scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="red", hjust=-0.5, size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_1
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_25_75_merged_100bp.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=26 & tf_hotspot_df$tf_counts<=75]
this_plot_2 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_y_continuous(expand=c(0,0), limits=c(0,250000)) +
			#scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="red", hjust=-0.5, size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_2
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_75_125_merged_100bp.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=76 & tf_hotspot_df$tf_counts<=125]
this_plot_3 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_y_continuous(expand=c(0,0), limits=c(0,250000)) +
			#scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="red", hjust=-0.5, size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_3
dev.off()


output_file_name <- paste0("~/Dropbox", "/", "tf_hotspot_prom_enh_assoc_dist_125_max_merged_100bp.pdf")       
pdf(output_file_name, width=6, height=6)

tf_hotspot_50 = tf_hotspot_df[tf_hotspot_df$tf_counts>=126]
this_plot_4 <- ggplot(tf_hotspot_50, aes(x=as.factor(tf_hotspot_50$tf_counts), y=total_sites, fill=total_sites)) + 
			geom_bar(stat="identity") + coord_flip() +
			scale_y_continuous(expand=c(0,0), limits=c(0,250000)) +
			#scale_fill_gradient("Count", low = "green", high = "red") + 
			geom_text(aes(label=total_sites), color="red", hjust=-1, size=2.5) +
			ylab("HepG2 sites count") +
			xlab("TF counts")
this_plot_4
dev.off()

