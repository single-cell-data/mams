# Overview

We adopt key-value pair list structure to store the metadata information of a data set. A list object in R, a YAML file or a JSON file works for this purpose. The basic structure consists of multiple data sets at the first level. Within each data set, we set fields for `FOM`, `OBS`, `VAR` and etc. Then we got all the FOMs in each data set embedded in the `FOM` field of each data set and so on so forth. Finally, under the element for each matrix, the fields to describe the metadata will be listed.. These fields are descibed in detail in the [Specifications](https://github.com/mvfki/mams/blob/main/specification.md). 

# Extra fields

## Linking fields
These fields are used to link the FOM to its corresponding ID and annotation objects. 

**Field:** oid  
**Value:** Character string  
**Description:** `id` of `OID` object that contains the observation IDs for this FOM.

**Field:** fid  
**Value:** Character string  
**Description:** `id` of `FID` object that contains the feature IDs for this FOM.

**Field:** obs  
**Value:** Character string or list  
**Description:** `id`(s) of `OBS` object(s) that contain annotations for the observation in this FOM.

**Field:** var  
**Value:** Character string or list  
**Description:** `id`(s) of `VAR` object(s) that contain annotations for the features in this FOM.


# Example

Here we take an example of a YAML file to how the list is structured and how elements match to the specifications. Comments on the right of some lines match to the explanation following this example. 

```{YAML}
pbmc8k:                     # [1]
  FOM:                      # [2]
    fom1:                   # [3]
      filepath: pbmc8k_seurat_raw.rds
      accessor: GetAssayData(object = pbmc8k_seurat_raw, slot = "counts", assay =
        "RNA")
      oid: oid1
      fid: fid1
      obs: obs1
      var: var1
      data_type: int
      representation: sparse
      obs_unit: cell
      processing: counts
      feature_subset: full
      modality: rna
      analyte: rna
      obs_subset: full
      record_id: CellRanger.count
  OID:                      # [2]
    oid1:                   # [3]
      filepath: pbmc8k_seurat_raw.rds
      accessor: '`colnames(pbmc8k_seurat_raw)`'
    oid2:
      filepath: pbmc8k_seurat_nonempty.rds
      accessor: '`colnames(pbmc8k_seurat_nonempty)`'
  FID:
    fid1:
      filename: pbmc8k_seurat_raw.rds
      accessor: '`rownames(Assays(pbmc8k_seurat_raw, slot = "RNA))`'
    fid2:
      filename: pbmc8k_seurat_raw.rds
      accessor: '`rownames(Assays(pbmc8k_seurat_raw, slot = "ADT))`'
  OBS:
    obs1:
      filepath: pbmc8k_seurat_raw.rds
      accessor: '`pbmc8k_seurat_raw[[]]`'
    obs2:
      filepath: pbmc8k_seurat_nonempty.rds
      accessor: '`pbmc8k_seurat_nonempty[[]]`'
  VAR:
    var1:
      filepath: pbmc8k_seurat_raw.rds
      accessor: '`pbmc8k_seurat_raw[["RNA"]][[]]`'
      feature_modality: rna
    var2:
      filepath: pbmc8k_seurat_raw.rds
      accessor: '`pbmc8k_seurat_raw[["ADT"]][[]]`'
      feature_modality: protein
  ONG:
    ong1:
      filepath: pbmc8k_seurat_filtered.rds
      accessor: '`Graphs(pbmc8k_seurat_filtered, "RNA_nn")`'
      edge_metric: euclidean
      metric_type: distance
      parent_id: fom11
      record_id: FindNeighbors.RNA.pca
  REC:
    CellRanger.count:
      record_package_name: CellRanger
      record_package_version: 3.0.0
      record_function_name: count
    NormalizeData.RNA:
      record_package_name: Seurat
      record_package_version: 4.1.0
      record_function_name: NormalizeData
      record_function_parameters: '[assay=RNA,normalization.method=LogNormalize,scale.factor=10000,margin=1,verbose=TRUE]'
```

#### Example Explanation

1. `pbmc8k`: The `dataset_id` of this data set, and all metadata for this data set should fall inside this field. The exact key name of this field can be arbitrarily set. 
2. `FOM`, `OID`, etc.: The fields under the `dataset_id` the groups all FOMs, all OIDs, etc. together, respectively. These key names should be fixed. 
3. `fom1`, `oid1`, `oid2`, etc.: The elements under the categories described in 2. and the keys stand for unique identifiers of a FOM, OID and etc. (i.e. `id` for each element). These names can be arbitrarily decided. 
