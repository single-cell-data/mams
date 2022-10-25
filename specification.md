
# Introduction
The past two decades have seen a rapid expansion of high throughput genomic and imaging technologies that have revolutionized the ability of researchers to capture the molecular and histological characteristics of biological samples. For example, assays such as single-cell RNA-seq can capture the states of individual cells within heterogeneous and complex tissues. Several major consortia have been funded to utilize single cell assays to create cellular atlases of healthy and disease tissues. These groups are generating large datasets that contain multi-modal single-cell data collected from longitudinally- and spatially-related biological specimens from different organs across different conditions. The majority of datasets are designed to answer a specific set of questions within a particular biological or clinical context. However, as more data becomes available, the ability to combine and integrate datasets from different settings is becoming increasingly apparent. Two major roadblocks to combining and integrating datasets are that 1) the data is stored in a wide variety of file formats or programming language-specific libraries, classes, or data structures and 2) the metadata about the matrix and the corresponding analysis that produced or utilized the matrix is not well standardized.

While a wide range of experimental protocols and platforms are available to generate molecular and histological data, an important commonality across these technologies is that they often produce a matrix of features that are measured in a set of observations. These feature and observation matrices (FOMs) are foundational for storing raw data from molecular assays (e.g. raw counts) and derived data from down-stream analytical tools (e.g. normalized matrix, reduced dimensional matrix). A variety of file formats are used to store FOMs on file systems in different representations. For example, Tab Separated Value (tsv/txt) files can be used to store dense matrices while Market Exchange (.mtx) files can be used to efficiently store sparse matrices. Although platform-independent, these formats do not readily capture relationships between matrices and annotations for features and observations. Several packages also exist that can capture relationships between matrices and annotations include AnnData in Python, the Seurat object in R, and the SingleCellExperiment package in R/BioConductor. In contrast to the simple flat file formats, these objects can capture complex relationships between FOMs as well as annotation data produced during the analysis such as quality control (QC metrics) and cluster labels. However, each package may store different sets of data or label the same type of data differently. Even if different datasets are stored the same format or type of object, harmonization of datasets still may require substantial manual curation before they can be combined and integrated. 

Lastly, a major goal of high-quality analysis workflows is to promote reproducibility by storing information related to provenance such as software version, function calls, and selected parameters that produced the matrix or annotation. However, there is a high degree of variability in which different analytical tools and software packages capture this information. Even if the data was produced by a workflow captured in a Docker container or versioned in GitHub, this information will often be lost when converting the data between formats or transferring between tools. Thus, there is a need to develop metadata standards for FOMs related to provenance to ensure this information can be readily captured and maintained throughout the dataset-specific analysis and during integration with other datasets.

In order to facilitate sharing and interoperability across consortiums and technologies and to promote reproducibility related to provenance, a detailed metadata schema describing the characteristics of FOMs can be used to serve as a standard for the community. Therefore, we developed the Matrix and Analysis Metadata Standards (MAMS) to capture the relevant information about the data matrices and annotations that are produced during common and complex analysis workflows for single-cell data. MAMS defines fields that describe what type of data is contained within a matrix, relationships between matrices, and provenance related to the tool or algorithm that created the matrix. These standards can serve as a roadmap for tool developers and data curators to ensure that their systems have the capability to store and retrieve relevant information needed for integration. All of the metadata fields are independent of the platform, programming language, and specific tool and thus can be used to support efforts to harmonize data across consortiums. 

# Quick reference table

| Field | Description | Class |
| -------------- | -------------------------- | --- |
| id | Denotes the unique id of the matrix, annotation data frame, or graph and should be unique within the set of objects within the same class of the dataset. | General |
| dataset_id | All matrices and data frames within this group should have observations and features that belong to a superset of observations and features that encompass an entire dataset. | General |
| class | Denotes the class of the matrix, array, data frame. | General |
| data_type | Explicitly describes the type of data stored in the FOM | FOM |
| representation | Preferred representation of the matrix. | FOM |
| obs_unit | Biological unit of the observations | FOM |
| processing | Used to describe the nature of the data contained within the matrix. | FOM |
| analyte | Used to describe the biological analytes being quantified in the matrix. | FOM |
| modality | Describes the modality of the matrix. | FOM | 
| obs_subset | Describes the subset of observations that are present in the FOM. | FOM |
| feature_subset | Describes the subset of features that are present in the FOM | FOM |
| parent_id | Denotes the id(s) of the parent matrices that were used to produce the matrix | FOM |
| parent_relationship | Denotes the type of relationship with the parent matrix or matrices | FOM |
| edge_metric | Name of the distance or similarity metric used to create the edges between features. | ONG |
| metric_type | “distance” indicates that smaller values denote more relatedness between features (e.g. euclidean distance) while “similarity” indicates that larger values denote more relatedness between observations (e.g. Pearson correlation) | ONG |
| edge_metric | Name of the distance or similarity metric used to create the edges between features. | FNG |
| metric_type | “distance” indicates that smaller values denote more relatedness between features (e.g. euclidean distance) while “similarity” indicates that larger values denote more relatedness between observations (e.g. Pearson correlation) | FNG |
| record_id | Unique id to denote a combination of entries for record_package_name, record_package_version, record_function_name, and record_function_parameters | 
| record_package_name | Name of the package, tool, or software that ran the algorithm to produce the matrix, annotation, or graph. | REC |
| record_package_version | Version of the package, tool, or software that ran the algorithm to produce the matrix, annotation, or graph. | REC |
| record_function_name | Name of the function or mathematical operation used to produce the matrix, annotation, or graph. | REC |
| record_function_parameters | Key/value pairs describing the primary parameters and their values used in the function call. | REC |
| record_workflow_link | Public link to workflow that ran the tool. | REC |


# General fields

These are fields that apply to all matrix, array and data frame related classes including `FOM`, `VAR`, `OBS`, `FID`, and `OID`. 

**Field:** id  
**Value:** Character string  
**Description:** Denotes the unique id of the matrix, annotation data frame, or graph and should be unique within the set of objects within the same class of the dataset. This ID could be a randomized unique ID such as a UUID or could be a unique combination of other fields from the FOM schema. For example, the IDs for FOMs could be a combination of `modality` and `processing` fields with the option of including algorithm_name. In this scenario, the ids “rna.raw”, “rna.normalize”, “rna.scale” could be used describe the data matrices and “rna.reduction.pca” and “rna.embedding.umap” could be used to describe the reduced dimensional objects.  

**Field:** dataset_id   
**Value:**  Character string  
**Description:** All matrices and data frames within this group should have observations and features that belong to a superset of observations and features that encompass an entire dataset.  Also applies to the `REC` class. 

**Field:** class  
**Value:** Character string  
**Description:** Denotes the class of the matrix, array, data frame. One of `FOM`, `VAR`, `OBS`, `FID`, or `OID`. Each class will have a correpsonding set of metadata fields associated with it.


# Feature Obsersvation Matrix class (FOM) 
A feature and observation matrix (FOM) is a data matrix that contains measurements of molecular features in biological entities. Examples of features include genes, genomic regions or peaks, transcripts, proteins, antibodies derived tags, signal intensities, cell type counts. Examples of observations include cells, cell pools, beads, spots, subcellular regions, and regions of interest (ROIs). Measurements may include transcript counts, protein abundances, signal intensities and velocity estimates. 

## Matrix description fields

**Field:** data_type  
**Value:** Character string    
**Description:** Explicitly describes the type of data stored in the FOM (e.g. int, int64, double, enum/categorical, etc). 

**Field:** representation  
**Value:** Character string    
**Description:** Preferred representation of the matrix.  

**Field:** representation_description  
**Value:** Character string    
**Description:** A brief description of the `representation` tag 

**Suggested categories for `representation` and `representation_description`:**

| representation | representation_description | Example |
| -------------- | ---------------------------| ----------------|
| sparse | The matrix contains zeros for the majority of the measurements. |
| dense | The matrix contains non-zeros for the majority of the measurements. |

**Field:** obs_unit  
**Value:** Character string    
**Description:** Biological unit of the observations  

**Field:** obs_unit_description  
**Value:** Character string    
**Description:** A brief description of the `obs_unit` tag 

**Suggested categories for `obs_unit` and `obs_unit_description`:**

| obs_unit | obs_unit_description | Notes/Examples |
| -------------- | ---------------------------| ----------------|
| bulk | Features are quantified for a collection of cells (e.g. bulk) such as tissue or culture | Bulk RNA-seq or ATAC-seq |
| cell | Features are quantified for individual cells | Single-cell RNA-seq |
| nucleus | Features are quantified for individual nuclei | Single-nucleus RNA-seq (snuc-seq) or single-nucleus ATAC-seq |


**Field:** processing  
**Value:** Character string    
**Description:** Used to describe the nature of the data contained within the matrix. This field will help distinguish between the matrices that are produced at different stages of an analysis workflow. This field should not be used to infer a history of previous steps in the analysis workflow. So users should not have to use a particular sequence of tags and should not assume that each matrix has gone through a set of prior tags. Rather, the history of previous steps should be captured in the provenance (LOG) fields.   

**Field:** processing_description  
**Value:** Character string    
**Description:** A brief description of the `processing` tag.   

**Suggested categories for `processing` and `processing_description`:**

| processing | processing_description | Notes/Examples |
| ---------- | -----------------------| ------- | 
| raw | Original measurements have not been altered |
| counts | Raw data for assays that produce integer-like data such as scRNA-seq. Child of `raw`.|  |
| intensities | Raw data for assays that produce continuous data. Child of `raw`. | |
| decontaminated | Measurements have been corrected for background signal.  | Ambient RNA removal in single-cell RNA-seq. |
| lograw | The log of the raw data |
| logcounts | The log of the raw counts | |
| logintensities | The log of the raw intensity values. | |
| corrected | Measurements have been corrected for observation-level covariates | |
| normalized | Data that has been normalized for differences in overall signal abundance between observations. | Correcting for total number of UMIs or reads in each cell. |
| lognormalized | Data that has been log transformed after normalizing for differences in overall signal abundance between observations. Child of `normalized`.| |
| centered | Data with features have been made to center around a standard quantity such as the mean or median. | Mean-centered data |
| scaled | Data with features have been centered around a standard quantity and standardized to have similar variances or ranges. | Z-scored data |
| reduction | A matrix containing a data dimensionality reduction generally useful for input into tools for downstream analysis such as clustering or 2D-embedding. | PCA, ICA, Autoencoders |
| embedding | A matrix containing a low dimensional embedding (usually 2D or 3D) generally used for visualization. Child of `reduction` | UMAP, tSNE | 

**Field:** analyte    
**Value:** Character string or list  
**Description:** Used to describe the biological analytes being quantified in the matrix. If more than one analyte was used to generate the values in a multi-modal analysis, then a list can be used to capture multiple analytes (e.g. a 2-D embedding generated from combined graph of RNA and protein data). 

**Field:** analyte_description  
**Value:** Character string    
**Description:** A brief description of the `analyte` field.   

**Suggested categories for `analyte` and `analyte_description`:**

| analyte | analyte_description | Notes/Examples |
| ------- | ------------------- | --------------- |
| rna | Used for technologies that measure RNA expression levels. | This should generally be used for assays listed under “RNA assay” (EFO_0001457) from the OLS. |
| dna | Used for technologies that measure features of DNA. | This should generally be used for assays listed under “DNA assay” (EFO_0001456) from the OLS. |
| chromatin | Used for technologies that measure open chromatin regions of DNA. |
| protein | Used for technologies that measure protein expression levels. | This should generally be used for assays listed under “protein assay” (EFO_0001458) from the OLS. Examples include CITE-seq, Total-seq, CODEX, MIBI. |
| morphology | Used for morphological measurements often derived from imaging technologies (e.g. cell size or shape). |
| lipid | Used for technologies that measure lipid levels. |
| metabolite | Used for technologies that measure metabolite levels. |
| methylation | Used for technologies that measure methtylation levels. |

**Field:** modality  
**Value:** Character string    
**Description:** Describes the modality of the matrix. If features or observations are of mixed modalities, then `feature_modality` in the `var` class or `observation_modality` in the FOM class should be used, respectively. This field may often be the same as another field or a combination of other fields such as `analyte` or `species`.

## Subset fields

**Field:** obs_subset  
**Value:** Character string    
**Description:** Describes the subset of observations that are present in the FOM. This field can be used to denote different subsets of observations that may be required at different stages of analysis after removing poor quality observations. Several suggested cateogories are provided for common quality control steps. This term can also be used describe a matrix containing a subset of observations based on biological characteristics. For example, after initial clustering and cell type identification, a subset of cells belonging to a particular cell type (e.g. T-cells) may be isolated and re-clustered to better characterize transitionary cell states within the cell type. This subset may have its own normalization, dimensionality reductions, embeddings, etc. 

**Field:** obs_subset_description  
**Value:** Character string    
**Description:** A brief description of the `obs_subset` field.   

**Suggested categories for `obs_subset` and `obs_subset_description`:**

| obs_subset | obs_subset_description | Examples |
| ---------- | ---------------------- | --------------- |
| full | Observations have not been filtered or subsetted. |
| filtered | Observations that have enough signal above background. | Droplets (or cell barcodes) that have enough counts to be considered to be non-empty. Similar to the “filtered” matrix from CellRanger. |
| threshold | Observations that have a total signal above a certain threshold | Filtering to include cells with a total UMI or read count above a certain threshold across features. |
| detected | Observations that have minimum levels of detection across features. | Filtering include cells with at least 3 counts in at least 3 genes. |
| nonartifact | A general term to describe filtering that may occur due other quality control metrics. | Artifacts in single cell RNA-seq data include high contamination from ambient material, high mitochondrial percentage, or doublets/multiplets. |
| clean | An “analysis ready” set of observations. | Cells that have been filtered for total signal, detection across features, outliers, and any other artifacts deemed to be important to filter against for quality control purposes. |

**Field:** feature_subset  
**Value:** Character string    
**Description:** Describes the subset of features that are present in the FOM.   

**Field:** feature_subset_description   
**Value:** Character string    
**Description:** A brief description of the `feature_subset` field.   

**Suggested categories for `feature_subset` and `feature_subset_description`:**

| feature_subset | feature_subset_description | Examples |
| -------------- | -------------------------- | ------- |
| full | Features have not been filtered or subsetted. | |
| threshold | Features that have a total signal above a certain threshold | Only including features with a total UMI or read count above a certain threshold across observations. | 
| detected | Features that have a minimum level of detection across observations. Only including features with at least 3 counts in at least 3 observations. |
| variable | Features that have minimum level of variability across all cells | The top 2,000 most variable features across observations |

## Provenance fields
These fields are used to link the FOM to a provenance record and describe the relationship between the FOM and its parent matrices (if any). 

**Field:** record_id  
**Value:** Character string or list  
**Description:** ID used to link FOMs to a specific record or set of records that produced the FOM. This should match an ID in the `REC` class.  

**Field:** parent_id  
**Value:** Character string or list    
**Description:** Denotes the id(s) of the parent matrices that were used to produce the matrix.

**Field:** parent_relationship  
**Value:** Character string   
**Description:** Denotes the type of relationship with the parent matrix or matrices. 

**Field:** parent_relationship_description  
**Value:** Character string   
**Description:** A brief description of the `parent_relationship_description` field.   

**Suggested categories for `parent_relationship` and `parent_relationship_description`:**

| parent_relationship | parent_relationship_description | Examples |
| -------------- | -------------------------- | ------- |
| transformation | Values have been modified but feature and observation dimensions are the same. | Normalization and Log2 transformation |
| subset | Values have been not been modified but a subset of feature and/or observation dimensions have been selected.  | Selecting all observations with a minimun number of counts |
| concatenation | Matrices that have been aggregated.  | Concatenating matrices of observations from two samples with the same features |
| reduction | Features have been reduced to a lower number of dimensions. | PCA, tSNE, UMAP |
| factorization | Decomposition of a matrix into two matrices of features and factors as well as factors and observations. | NMF, LDA |
| aggregation | Aggregating values for groups of features or observations.  | Taking the average of each feature within each group of observations |


# Observation ID class (OID)
An observation_id is character vector or combination of character vectors used to denote the unique ID of each observation. The number of elements in the vector(s) should be the same length as the number of observations in the FOM. If multiple vectors, then the combination of elements across the vectors for each combination should be unique. The number of elements in the vector(s) should be the same length as the number of observations in the FOM. Often compound IDs are created by concatenating multiple types of IDs together when matrices from different sources need to be combined. For example, a cell barcode may uniquely define a cell within a given sample, but the same cell barcode may be used across samples. To combine cell matrices from different samples, a sample ID can be concatenated with the cell barcode to make each observation ID unique across a group of samples. If the `OID` consists of multiple vectors, then the combination of elements at each position across the vectors should be unique. If the OID is a single character vector that is a compound ID, then a character delimiter should be used to separate the fields within a string and denoted with the optional delim field. For example, if “Sample_A” has a cell with a barcode of “ACGT” and the delimiter is chosen to be “.”, then the compound ID would be “Sample_A.ACGT”.

**Field:** oid_header  
**Value:** Character string  
**Description:** Describes headers contained with the observation IDs. For example if the matrix contains RNA expression data for cells, then this field can be labeled “cell_id”. If the OID is a compound ID with multiple vectors, then this should be a vector of the same length, essentially representing the headers of the OID matrix. If the OID is a single-string compound ID, then this field should denote each component of the ID separated by the delim character. For example, if the sample_id and cell_id are present in the id, then this field should be “sample_id.cell_id”.

**Field:** oid_header_delim  
**Value:** Character string  
**Description:** A character string denoting the delimiter that can be used to separate compound observation IDs. The recommended delimiter is a period (e.g. “.”). If this field is not included, then it will be assumed that the OID is not a compound ID. 

# Feature ID class (FID)
A feature_id is a character vector or combination of character vectors used to denote the unique ID of each feature. The number of elements in the vector(s) should be the same length as the number of features in the FOM. If multiple vectors, then the combination of elements across the vectors should be unique. 

**Field:** fid_header  
**Value:** Character string  
**Description:** Describes headers contained with the feature IDs. For example if the matrix contains RNA expression data for cells, then this field can denote type of ID used (e.g. ENSEMBL). If the FID is a compound ID with multiple vectors, then this should be a vector of the same length, essentially representing the headers of the FID matrix. If the FID is a single-string compound ID, then this field should denote each component of the ID separated by the delim character. For example, if the ENSENBL ID and GENE SYMBOL are concatenated in the FID with the period used as a delimiter, then this field should be “ENSEMBL_ID.SYMBOL”.

**Field:** fid_header_delim  
**Value:** Character string  
**Description:** A character string denoting the delimiter that can be used to separate compound feature IDs. The recommended delimiter is a period (e.g. “.”). If this field is not included, then it will be assumed that the FID is not a compound ID. 


# Observation Annotation class (OBS)
`OBS` objects are matrices or data frames with the same number of observations as its corresponding FOM. It is used to store annotations and observation-level metadata such as quality control metrics (e.g. total counts), sample demographics (e.g. age, disease status), and analysis results (e.g. cluster labels, trajectory scores).

**Field:** record_id  
**Value:** Character string or list  
**Description:** ID used to link an OBS annotation to a specific record or set of records. This should match an ID in the REC class.  


# Feature Annotation class (VAR)
`VAR` obects are matricies or data frames with the same number of features as its corresponding FOM. It is used to store annotations and feature-level metadata such as IDs (e.g. Ensembl, Gene Symbol), reference information (e.g. chromosome coordinates, gene biotype), and analysis results (e.g. variability metrics, cluster labels).

**Field:** record_id  
**Value:** Character string or list  
**Description:** ID used to link an FEA annotation to a specific record or set of records. This should match an ID in the REC class.  

# Observation Neighborhood Graph class (ONG)
Observation Neighborhood Graphs (i.e. adjacency matrices) can be used to store the correlation, similarity, or distance between pairs of observations. These measurements are often used by graph-based clustering and visualization tools.

**Field:** edge_metric  
**Value:** Character string  
**Description:** Name of the distance or similarity metric used to create the edges between observations.  

**Field:** metric_type  
**Value:** Character string. One of “distance” or “similarity”.  
**Description:** “distance” indicates that smaller values denote more relatedness between observations (e.g. euclidean distance) while “similarity” indicates that larger values denote more relatedness between observations (e.g. Pearson correlation). 

# Feature Neighborhood Graph class (FNG)
Feature neighborhood graphs (i.e. adjacency matrices) can be used to store the correlation, similarity, or distance between pairs of features. These measurements are often used by graph-based clustering and visualization tools.

**Field:** edge_metric  
**Value:** Character string  
**Description:** Name of the distance or similarity metric used to create the edges between features.  

**Field:** metric_type  
**Value:** Character string. One of “distance” or “similarity”.  
**Description:** “distance” indicates that smaller values denote more relatedness between features (e.g. euclidean distance) while “similarity” indicates that larger values denote more relatedness between observations (e.g. Pearson correlation). 


# Provenance record class (REC)
These fields are used to capture provenance about the tool, software package, version, functions, and parameters used to create the various FOMs, annotations, and graphs. 

**Field:** record_id  
**Value:** Character string  
**Description:** Unique id to denote a combination of entries for `record_package_name`, `record_package_version`, `record_function_name`, and `record_function_parameters` 

**Field:** record_package_name  
**Value:** Character string  
**Description:** Name of the package, tool, or software that ran the algorithm to produce the matrix, annotation, or graph.

**Field:** record_package_version  
**Value:** Character string  
**Description:** Version of the package, tool, or software that ran the algorithm to produce the matrix, annotation, or graph.

**Field:** record_function_name  
**Value:** Character string    
**Description:** Name of the function or mathematical operation used to produce the matrix, annotation, or graph.

**Field:** record_function_parameters  
**Value:** Character string or list   
**Description:** Key/value pairs describing the primary parameters and their values used in the function call.

**Field:** record_workflow_link  
**Value:** Character string  
**Description:** Public link to workflow that ran the tool. For example, links to workflow scripts such as CWL, WDL, Nextflow, etc. stored in a public repository such as GitHub, DockerHub, etc. 

**Field:** record_runtime_start  
**Value:** Character string in ISO 8601 format  
**Description:** The start time of the algorithm or operation.  

**Field:** record_runtime_end  
**Value:** Character string in ISO 8601 format  
**Description:** The finishing time of the algorithm or operation.  

**Field:** record_runtime_duration  
**Value:** Character string  
**Description:** The total duration of the algorithm or operation.  


