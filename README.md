# A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles

# Gene Expression Analysis and Cell Line Prediction

This project involves the analysis of gene expression data and the prediction of cell lines based on this data. The data is initially processed and analyzed in R, and then further analyzed and modeled in Python.

## Workflow

1. **Data Preparation**: The gene expression data is loaded from multiple files and combined into a single data frame. Metadata about the samples is retrieved from the GEO database.

2. **Differential Expression Analysis**: The `limma` package in R is used to perform a differential expression analysis. The genes with |logFC| ≥ 1 and pvalues
< 0.05 are considered differentially expressed genes (DEGs) in the dataset. The genes that are common between the datasets are taken for the Pathway analysis. The results are visualized using `ggplot2`.
<img width="767" alt="Screenshot 2024-03-26 at 4 27 44 PM" src="https://github.com/Anube9/A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles/assets/112353734/322c6929-26de-4d78-8f90-565055832321"><br>
<img width="748" alt="Screenshot 2024-03-26 at 4 28 32 PM" src="https://github.com/Anube9/A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles/assets/112353734/328c2165-44f4-4fa1-af13-cf32d81e5b3d"><br>

4. **Functional Enrichment and Pathway Analysis & Construction of Protein-Protein Interaction Network** : Functional enrichment analysis of the DEGs was implemented
using DAVID (Database for Annotation, Visualization, and Integrated Discovery) along with pathway analysis. GO enrichment analysis of the DEGs was performed.
Using STRING a protein-protein interaction network of the DEGs is constructed with a confidence of 0.7.<br>
![image](https://github.com/Anube9/A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles/assets/112353734/6ae6f2df-8bad-40f5-b34d-332aad3a3f16)<br>
Network of the cell_line_1_D<br>
![image](https://github.com/Anube9/A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles/assets/112353734/9d65f063-065b-48eb-88e0-1d7305cf4270)<br>
Network of the cell_line_1_D_R<br>
![image](https://github.com/Anube9/A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles/assets/112353734/b0ad08d2-2f45-4746-80d4-82124328da0d)<br>
Network of the cell_line_2_D<br>
![image](https://github.com/Anube9/A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles/assets/112353734/7e359d42-e3dc-43e0-8613-90a5aadd44d6)<br>
Network of the cell_line_2_D_R<br>
In this study, two datasets from the cell line MDAMB-453 were analyzed.In the first dataset,genes with nodal degree >5 were identified as hub genes in the network, including MYC, ALB, EGR1,HSPG2, and CDH2.The second dataset from MDA-MB453 identified EGFR, ALB, EDNI, HSPG2, EDNI, and TFBG1 as hub genes based on nodal degree. In the ACC-422 cell line treated with drug, radiation and drug (Fig 7and Fig 8) were included in the analysis. In the ACC-422first dataset, FOS, ESR1, EGR1, FBN1, CD4, and CDH2 were identified as hub genes based on nodal degree. In the second dataset of ACC-422, TGFB1, IGF1, CD4, ESR1,ALB, and AR were found to be the hub genes based on nodal degree.<br>
6. **Machine Learning Prediction**: The gene expression data and cell line labels are imported into Python. A decision tree classifier is trained on a subset of the data, and then used to predict the cell line for the remaining samples.<br>
 **Performance Evaluation**: The performance of the classifier is evaluated by generating a classification report and a confusion matrix, and by plotting a ROC curve. The decision tree is also visualized.<br>
![image](https://github.com/Anube9/A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles/assets/112353734/d68e7efd-0080-401f-a649-1a0ae2a6ce38)<br>
The above ROC curve shows the performance of the decision tree model with 88% of accuracy.<br>
![image](https://github.com/Anube9/A-Machine-Learning-Approach-to-Breast-cancer-Cell-Line-Identification-using-Gene-Expression-Profiles/assets/112353734/3dc4ca80-4dd5-48bc-95d9-aab6527fd531)<br>
The above figure shows the decision tree of the model showing the number of samples, value of prediction and gini which is the entropy of the dataset which summed up to be zero.<br>
**Conclusion:** 
This project demonstrates the use of bioinformatics and machine learning techniques to analyze gene expression data and predict cell lines. It provides a workflow that can be adapted for similar analyses.
