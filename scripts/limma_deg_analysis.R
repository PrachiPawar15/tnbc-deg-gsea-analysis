############################################################
# Transcriptomic DEG Analysis Workflow — TNBC Case Study
# Data source: GEO (GSE206912)
# Author: Prachi Pawar
# Environment: R / RStudio
# Description: End-to-end workflow including preprocessing,
# PCA, clustering, Limma-based DEA, and GSEA.
############################################################

# Step 1: Download and load data from GEO
library(GEOquery)
geo_accession <- "GSE206912"
gse_list <- getGEO(geo_accession, GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse_list[[1]]
#Getting Expression Data gene_exp<- exprs(gse) gene_exp<- data.frame(gene_exp) #Explore and View the Data
dim(gene_exp) #dim displays the dimensions of the data (number of rows and columns) colnames(gene_exp)#colnames shows the column names.
View(gene_exp)
g_expr<- t(gene_exp)  #transpose the data
# Extract Feature Data and Phenotype Data
f_data <- fData(gse) #Extracts feature data (f_data) from the GEO dataset. View(f_data)
p_data <- pData(gse) # Extracts phenotype data (p_data) from the GEO dataset. View(p_data)
#Replacing the probe ID's in the column names with Gene Symbols
colnames(g_expr)<- f_data$GENE_SYMBOL
# Finding and curating the indices of non-empty strings
keep_genes <- which(colnames(g_expr) != "")
# Subseting the dataframe based on non-empty indices
g_exprs <- g_expr[, keep_genes]
#Standard Normalization of Data using Z-Score Standardization
exp_std<- (g_exprs - mean(g_exprs))/ sd(g_exprs)
#density plot of Z-score Standardized data
density_plot<- plot(density(exp_std, main="Density plot of Z-score standardized data", xlab="Z value", ylab="Density"))


#Step 2: Unsupervised learning Using K-Means clustering and Elbow Method #Loading required libraries
library(cluster) library(factoextra) library(ggfortify) library(edgeR) library(ggplot2)
# Visualize Elbow Plot:
#Visualizing the elbow plot to help determine an appropriate number of clusters.
fviz_nbclust(exp_std, kmeans, method = "wss")
#performing k-means clustering with 2 centers and 25 different initial configurations.
km <- kmeans(exp_std, centers = 2, nstart = 25) table(km$cluster)
#View and Assign Clusters
View(km) #Views the k-means clustering results
final_data <- cbind(exp_std, cluster = km$cluster) #adds the cluster assignment to the data
# Combine the expression data with the cluster assignments
final_data_std <- data.frame(exp_std, cluster = km$cluster)
# View the final data frame
View(final_data_std)


#Step 3: Principle Componant Analysis representing different clusters
gene_data_pca <- final_data_std[, -c(1, ncol(final_data_std))]
# Performing PCA
pca_result <- prcomp(gene_data_pca) # Printing summary of PCA summary(pca_result) 
# Accessing PC scores
pc_scores <- pca_result$x
# Adding cluster information to PC scores
pc_scores_with_cluster <- cbind(pc_scores, Cluster = final_data_std$cluster)
# Creating a data frame for plotting
plot_data <- as.data.frame(pc_scores_with_cluster)
# Defining colors for each cluster
cluster_colors <- c("purple", "orange")
# Plotting PCA results using ggplot2 with different colors for each cluster ggplot(plot_data, aes(x = PC1, y = PC2, color = as.factor(Cluster))) + geom_point(size = 3) +
labs(title = "PCA of Gene Expression Data", x = "Principal Component 1",
y = "Principal Component 2") + scale_color_manual(values = cluster_colors,
labels = c("Cluster1","Cluster2")) +# Set custom colors labs(color ="Clusters")
theme_minimal()

#Step 4: Performing Differential Expression Analysis #Loading necessary libraries
library(limma)
# Specifying design matrix for differential expression analysis
design_matrix <- model.matrix(~ 0 + factor(final_data_std$cluster)) #design matrix colnames(design_matrix)=c("CL1","CL2")
View(design_matrix)
final_data_std<-final_data_std[,(1:ncol(final_data_std)-1)]
final_data_std<-data.frame(t(final_data_std)) #convert final_data to dataframe
#Fitting the linear model to gene expression data using design matrix in lmFit function
fit <- lmFit(final_data_std, design_matrix)
#Making contrast with design matrix
con<-makeContrasts(CL1-CL2, levels = design_matrix)
#Fitting the contrast
fit_contrasts <- contrasts.fit(fit, con)
#Empirical Bayes moderation is applied using eBayes
fit_ebayes <- eBayes(fit_contrasts) View(fit_ebayes)
#Extracts the differential expression results for each cluster using topTable.
#coef specifies the contrast of interest
results_std <- topTable(fit_ebayes, number = Inf) View(results_std)
# Volcano Plot for visualizing differentially expressed genes
ggplot(data	=	results_std,	aes(x	=	results_std[,	"logFC"],	y	=	-log10(results_std[, "P.Value"])))+
geom_vline(xintercept = c(-0.6, 0.6), col = "black", linetype = 'dashed') + geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + geom_point() +
theme_set(theme_classic(base_size = 12) + theme(
axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1), color = 'black'),
axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1), color = 'black'),
plot.title = element_text(hjust = 0.5)
))
# Adding 'diffexpressed' column to results_std
results_std$diffexpressed <- "NOT SIGNIFICANT" results_std$diffexpressed[results_std$logFC > 0.6 & results_std$P.Value < 0.05] <- "UP" results_std$diffexpressed[results_std$logFC < -0.6 & results_std$P.Value < 0.05] <- "DOWN"
# Updated ggplot with the correct dataframe (results_std) and labels
ggplot(data = results_std, aes(x = results_std[, "logFC"], y = -log10(results_std[, "P.Value"]), col = diffexpressed)) +
geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') + geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + geom_point(size = 2) +
scale_color_manual(values = c("green", "grey", "blue"), 
labels = c("Downregulated", "Non-significant", "Upregulated")) + coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) +
labs(color = 'Expression ', x = expression("Log Fold Change"), y = expression("- log"[10]*"P-value")) +
scale_x_continuous(breaks = seq(-10, 10, 2)) View(results_std)
# Filtering results to get more significant genes #For upregulated genes
filtered_res_up <- results_std[abs(results_std$logFC) > 1.25 & (results_std$adj.P.Val) < 1e- 10, ]
#For downregulated genes
filtered_res_down<- results_std[abs(results_std$logFC) > -2 & (results_std$adj.P.Val) < 1e- 15, ]


# Step 5: Performing GO enrichment analysis for upregulated genes #Loading required libraries
library(clusterProfiler) library(org.Hs.eg.db) #Filtered results
filtered_res_up <- results_std[abs(results_std$logFC) > 1.25 & (results_std$adj.P.Val) < 1e- 10, ]
filtered_res_down <- results_std[abs(results_std$logFC) > -2 & (results_std$adj.P.Val) < 1e- 15, ]
# Function to perform enrichment analysis and save results perform_enrichment <- function(gene_identifiers, ontology) { enrich_result <- enrichGO(gene = gene_identifiers,
OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = ontology, pvalueCutoff = 0.05,
qvalueCutoff = 0.2, readable = TRUE)
return(enrich_result)
}
#Performing enrichment analysis for upregulated genes gene_identifiers_up <- rownames(filtered_res_up) UP_BP<-perform_enrichment(gene_identifiers_up, "BP") #dotplot for upregulated genes
en_dotplot_ex_u<-dotplot(UP_BP, showCategory = 5) en_dotplot_ex_u <- dotplot(
UP_BP,
showCategory = 5,
title = "Enrichment Dotplot for Upregulated genes", color = "p.adjust"
)
print(en_dotplot_ex_u)
# Performing enrichment analysis for downregulated genes gene_identifiers_down <- rownames(filtered_res_down) DOWN_BP<-perform_enrichment(gene_identifiers_down, "BP") #dotplot for downregulated genes
en_dotplot_ex_d<-dotplot(UP_BP, showCategory = 5) en_dotplot_ex_d <- dotplot(
DOWN_BP,
showCategory = 5,
title = "Enrichment Dotplot for Downregulated genes", color = "p.adjust"
)
#Save all the files
write.csv(geodf,"GEO_data.csv") 
write.csv(exp_norm_df,"Normalized Expression Data.csv") 
write.csv(K_means_clustering2,"K_means Clustering Data.csv")
write.csv(cluster_data,"Cluster Data.csv") 
write.csv(plot_data,"PCA Plot data.csv") 
write.csv(DEA_results,"DEA_Results.csv") 
write.csv(filtered_results_up, "filtered_results_up.csv") 
write.csv(filtered_results_down, "filtered_results_down.csv")


#Step 6: PPI analysis done using STRING database


#Step 7: The node interactions are analyzed and the Degree centrality calculation to find Hub genes
#Importing required libraries
import pandas as pd import networkx as nx
# Read interaction data of upregulated genes
interaction_data_upregulated	=
pd.read_csv('results/string_interactions_upreg_ug.tsv', sep='\t')
# Create an empty graph
G_upregulated = nx.Graph()
# Add edges to the graph based on the interaction data for downregulated genes
for index, row in interaction_data_upregulated.iterrows(): gene_a = row['#node1']
gene_b = row['node2'] G_upregulated.add_edge(gene_a, gene_b)
# Calculate node degrees for downregulated genes node_degrees_upregulated = dict(G_upregulated.degree()) # Identify Hub Genes for Downregulated Genes
# Sort the genes by degree and select the top hub genes
sorted_genes_upregulated = sorted(node_degrees_upregulated.items(), key=lambda x: x[1], reverse=True)
top_hub_genes_upregulated = sorted_genes_upregulated[:10] # Assuming you want top 10 hub genes
# Print or save the results
print("Top 10 Hub Genes for Upregulated Genes:") for gene, degree in top_hub_genes_upregulated:
print(f"Gene: {gene}, Degree: {degree}")
# Read interaction data of downregulated genes
interaction_data_downregulated	=
pd.read_csv('results/string_interactions_downwreg_ug.tsv', sep='\t')
#Create an empty graph
G_downregulated = nx.Graph()
#Add edges to the graph based on the interaction data for downregulated genes
for index, row in interaction_data_downregulated.iterrows(): gene_a = row['#node1']
gene_b = row['node2'] G_downregulated.add_edge(gene_a, gene_b)
# Calculate node degrees for downregulated genes node_degrees_downregulated = dict(G_downregulated.degree()) # Identify Hub Genes for Downregulated Genes
# Sort the genes by degree and select the top hub genes
sorted_genes_downregulated = sorted(node_degrees_downregulated.items(), key=lambda x: x[1], reverse=True)
top_hub_genes_downregulated = sorted_genes_downregulated[:10]	# Assuming you want top 10 hub genes
# Print or save the results
print("Top 10 Hub Genes for Downregulated Genes:") for gene, degree in top_hub_genes_downregulated:
print(f"Gene: {gene}, Degree: {degree}")

