<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/limiric.png" alt="limiric_logo" height="150" />
</p>

## Contents
[Description](#description) | [Installation](#installation) | [Quickstart](#quickstart)  |  [Basic Usage](#basic-usage)  |  [Output explained](#output-explained)  |  [More Information](#more-information) | [Input file format](#input-file-format) | [Downstream](#downstream)


## Description

Single cell RNA sequencing quality control package for sample-specific damaged cell detection through low dimension mitochondrial and ribosomal cluster selection.

## Installation
### Prerequisites
```limiric``` requires the following packages to be installed in your ```R``` environment  
* ```cowplot``` Wilke, 2024. [GitHub repo](https://github.com/wilkelab/cowplot)
* ```DropletQC``` Muskovic, 2024. [GitHub repo](https://github.com/powellgenomicslab/DropletQC)
* ```dplyr``` Wickham _et al_, 2023. [GitHub repo](https://github.com/tidyverse/dplyr)
* ```ggplot2``` Wickham, 2016. [GitHub repo](https://github.com/tidyverse/ggplot2)
* ```Matrix``` Bates _et al_, 2024. [Github repo](https://github.com/cran/Matrix)
* ```png``` Urbanek, 2022. [GitHub](https://github.com/cran/png)
* ```Seurat``` Satija and Farrell _et al_, 2015. [Website](https://satijalab.org/seurat/), [GitHub repo](https://github.com/satijalab/seurat)
* ```SoupX``` Young and Behjati, 2020. [GitHub repo](https://github.com/constantAmateur/SoupX)
<br>

If not already, you can install them together

```R

packages <- c("cowplot", "dplyr", "ggplot2", "Matrix", "png", "Seurat", "SoupX")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }

```
<br>

But ```DropletQC``` must be installed as follows 
```R
devtools::install_github("powellgenomicslab/DropletQC")
```

<br>

Be sure to load these into your environment before continuing with ```limiric``` installation 
```R
packages <- c("cowplot", "dplyr", "DropletQC", "ggplot2", "Matrix", "png", "Seurat", "SoupX")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg)
  }
}

```
> Please check the associated documentation if problems occur in the installation of any of the prerequisite packages


<br><br>

### Package installation
After all the prerequisites are installed, you can install the latest development version of ```limiric``` 

```R
devtools::install_github("AlicenJoyHenning/limiric")
```
<br>
<br>

## _Quickstart_

#### Basic Usage


1. Detect damaged cells in a sample using the filtered alignment files (```barcodes.tsv.gz```, ```features.tsv.gz```, ```counts.mtx.gz```). All you need to input 
is your sample name, the directory where your alignment files are stored, and a directory where the ```limiric``` output can be created.  
```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```  
> **NB** Please keep the following in mind:
> * Your files must be zipped, see [Input file format](#input-file-format) for more
> * The above usage assumes the sample is of human origin (See parameter: ```Organism```)
> * When compiling the data, all genes will be retained (See parameter: ```MinCells```)
> * Red blood cells will automatically be removed from the sample before damaged cell detection (See parameter: ```FilterRBC```)
> * Your output ```Seurat``` object will be filtered, containing only undamaged cells (See parameter: ```FilterOutput```)
>

<br>
<br>

2. Alternatively, you can use a ```Seurat``` object as input
```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    SeuratInput  = seurat_object,
    OutputPath   = "/home/user/alignment/limiric/"
)
```

<br>
<br>
  
3. If you have more than one sample, you may want to define a sample list
instead of running each individually.
You can do this through the ```sample_list``` parameter.
```R
sample_list <- list(

    list(ProjectName  = "SRR1234567",
         FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
         OutputPath   = "/home/user/alignment/limiric/"),
    
    list(ProjectName  = "SRR1234568",
         FilteredPath = "/home/user/alignment/SRR1234568/filtered/",
         OutputPath   = "/home/user/alignment/limiric/"),
    
    list(ProjectName  = "SRR1234569",
         FilteredPath = "/home/user/alignment/SRR1234569/filtered/",
         OutputPath   = "/home/user/alignment/limiric/")
)

GSE1234567 <- limiric(sample_list = sample_list)

```
> **NB** Please note that you must define your sample list with each parameter specified, i.e. ```limiric```
> won't know what to do if your list looks like this :
> ```R
> sample_list <- list(
>
>    list("SRR1234567",
>         "/home/user/alignment/SRR1234567/filtered/",
>         "/home/user/alignment/limiric/"),
>    list("SRR1234568",
>         "/home/user/alignment/SRR1234568/filtered/",
>         "/home/user/alignment/limiric/"),
>    list("SRR1234569",
>         "/home/user/alignment/SRR1234569/filtered/",
>         "/home/user/alignment/limiric/")
>  )
>``` 
> 
<br>

## Output explained   

Before interpretting ```limiric```'s outputs, it may be helpful to understand a bit of what's going on in the background.
<br><br>
As you know, the magic of scRNA-seq revolves around the identification of distinct cell populations. In majority of scRNA-seq analysis workflows, distinct populations are identified according to the
clusters that form when dimensionality reduction is performed on the most variable genes in a sample. For instance, the typical ```Seurat``` workflow uses the top two or three thousand variable genes.  

In our case, instead of trying to isolate distinct cell-type associated populations, we are only interested in two broad populations: _damaged_ and _healthy_ cells.
In general, most of what we know about how these two populations differ can be summarized by looking at the expression of mitochondrial and ribosomal genes. Thus, if we perform dimensionality reduction and clustering 
using _only_ these genes, we expect to see a division of cells into two clusters associated with the cell's status as either _damaged_ or _healthy_.  

This is the basis of the ```low dimension mitochondrial & ribosomal clustering```, aka ```limiric```, algorithm. 

### Output directory 

The ```limiric``` output will be created in the ```OutputPath```directory with the following structure  

```
OutputPath/
|
├── RBCQC
|
├── CellQC
|
└── Filtered
```

<br>

### **_RBCQC_** 

Before damaged cells can be identified, ```limiric```first removes red blood cells, or cells that are highly contaminated with haemoglobin, from your data. This is done under the assumption that these cells will not be informative to your study. If this assumption should not be true, you can avoid this filtering using ```FilterRBC = FALSE```.
<br>

The output here comes in the form of a scatter plot where the red blood cells are coloured in blue. The percentage of the total cells that were removed is also given in the plot. 

<br><br>
<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/RBCQC.png" height="500" />
</p>

> **NB** Given many scRNA-seq protocols, such as the ```10X Genomics``` protocol, advise for globin treatment, we hope for the contamination percentage to be as low as possible.
<br>

### **_CellQC_** 

This directory contains core diagnostic plots for the ```limiric``` damaged cell detection algorithm. As shown below, you will see four tSNE plots.
<br><br>
In this example dataset, you can see the cells divide into two distinct clusters; but which contain _damaged_ cells and which contain _healthy_ cells? From literature we know that damaged cells have high mitochondrial expression and low ribosomal expression. This is because damaged cells, such as those undergoing apoptosis, are characterised by compromised cell membranes. In this case, free cytoplasmic RNA- like ribosomal RNA- is more likely to have escaped before being captured and sequenced, meaning its expression will be low. While the RNA enclosed within the mitochondria, which has its own membrane, is retained, meaning its expression will be high. But many other factors likely contribute to this phenomenon, including impaired ribosomal biogenesis in damaged cells, leading to reduced ribosomal expression. 
<br><br>
<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/CellQC.png" height="500" />
</p>
<br><br>

Now answering the question of which cluster is which is easy with the **_CellQC_** tSNEs where the mitochondrial and ribosomal gene expression of each cell is made visible.  In the example below, the ```Mitochondrial gene expression``` tSNE shows that cells in the smaller cluster express mitochondrial genes at a very high level, suggesting they are likely damaged. The same conclusion can be reached looking at the ```Ribosomal gene expression``` tSNE where cells in the smaller cluster expresses ribosomal genes at a very low level.  


This can be further verified by looking at the ```Complexity score``` tSNE. This score measures the total number of mitochondrial and ribosomal genes expressed in a cell. Cells that express a large number of these genes are more likely to be metabolically active, undamaged cells while those that only express a few are likely not. In the  ```Complexity score``` plot, the cluster suspected to be damaged contains cells with very low complexities verifying by another metric that the cells are likely damaged. 
<br>
Together, these metrics are used by the ```limiric``` algorithm to annotate the _damaged_ cells, as shown in the fourth and final tSNE. So in summary, by looking at the ```limiric``` tSNE plots, we can immediately tell that we have a small cluster of cells that are likley damaged and should be removed from the sample. This makes sense as, unless something went drastically wrong during sample preparations, you expect there to be far more healthy than damaged cells in your scRNA-seq data.
<br>
<br>



### **_Filtered_** 
<br>

The main output of the ```limiric``` function comes in the form of a ```barcodes.csv``` containing annotations for each cell barcode of the input data. This can easily be incorporated into an existing scRNA-seq analysis workflow for filtering (see [example](#downstream)).

#### ProjectName_barcodes.csv  

|     Barcode      | Limiric |
|:---------------:|:-------:|
| AAACCTGAGATAGGAG|  cell   |
| AAACCTGAGCTATGCT|  cell   |
| AAACCTGAGCTGTTCA| damaged |
| AAACCTGCACATTAGC| damaged |
| ... | ... |


Additionally, ```limiric``` will output a ```Seurat``` object with ready-filtered barcodes for seamless integration at the start of a ```Seurat``` pre-processing workflow. 



## More Information
#### _Slightly-less_ Basic Usage

1. Perform ambient RNA correction
It is highly advised to correct for ambient RNA in your scRNAseq data. If you haven't already performed some kind of correction,
the ```SoupX``` parameter allows you to do so. This will occur before anything else in the ```limiric``` workflow.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    SoupX        = TRUE,
    RawPath      = "/home/user/alignment/SRR1234567/raw/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```


2. Isolate immune cells
<br>

If you have a sample where only the immune cells are of interest, include the following ```IsolateCD45```
parameter. This will isolate the immune cells present in the sample, then identify damaged cells.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    IsolateCD45  = TRUE,
    OutputPath   = "/home/user/alignment/limiric/"
)
```

<br>

> **NB** This will change your output directory structure by adding a new ```IMCQC``` layer
>
> ```
> OutputPath/
> ├── CellQC
> ├── IMCQC
> ├── RBCQC
> └── Filtered
> ```

### **_IMCQC_**
This output, like RBCQC, contains a scatter plot showing the removed cell in blue and the retained cells in grey.

<br>

<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/IMCQC.png" height="500" />
</p>

<br>
<br>

3. Combine ```limiric``` annotations with ```DropletQC```
<br>

Detect damaged cells and compare results with those from ```DropletQC```.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    DropletQC    = TRUE,
    VelocytoPath = "/home/user/alignment/velocyto/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```
> **NB** This will change your output directory structure by adding a new ```DropletQC``` layer
>
> ```
> OutputPath/
> ├── CellQC
> ├── DropletQC
> ├── RBCQC
> └── Filtered
> ```
>
### **_DropletQC_**
This will output a scatter plot and tSNE showing the cells annotated as _damaged_ by both ```limiric``` and ```DropletQC```.

<br>

<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/DropletQC.png" height="500" />
</p>

<br>

4. Combine previous conditions
<br>
Perform ambient RNA correction with ```SoupX```, filter red blood cells, isolate immune cells, detect damaged cells, and compare against ```DropletQC```.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    SoupX        = TRUE,
    RawPath      = "/home/user/alignment/SRR1234567/raw/",
    DropletQC    = TRUE,
    IsolateCD45  = TRUE,
    VelocytoPath = "/home/user/alignment/velocyto/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```
> **NB** This will mean your output directory looks like this
>
> ```
> OutputPath/
> ├── CellQC
> ├── DropletQC
> ├── IMCQC
> ├── RBCQC
> └── Filtered
> ```

<br>
<br>


5. Process multiple samples with the same conditions as in **Example 4**

```R
sample_list <- list(
    list(ProjectName  = "SRR1234567",
         FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
         SoupX        = TRUE,
         RawPath      = "/home/user/alignment/SRR1234567/raw/",
         DropletQC    = TRUE,
         IsolateCD45  = TRUE,
         VelocytoPath = "/home/user/alignment/velocyto/",
         OutputPath   = "/home/user/alignment/limiric/"),
    
    list(ProjectName  = "SRR1234568",
         FilteredPath = "/home/user/alignment/SRR1234568/filtered/",
         SoupX        = TRUE,
         RawPath      = "/home/user/alignment/SRR1234568/raw/",
         DropletQC    = TRUE,
         IsolateCD45  = TRUE,
         VelocytoPath = "/home/user/alignment/velocyto/",
         OutputPath   = "/home/user/alignment/limiric/"),
    
    list(ProjectName  = "SRR1234569",
         FilteredPath = "/home/user/alignment/SRR1234569/filtered/",
         SoupX        = TRUE,
         RawPath      = "/home/user/alignment/SRR1234569/raw/",
         DropletQC    = TRUE,
         IsolateCD45  = TRUE,
         VelocytoPath = "/home/user/alignment/velocyto/",
         OutputPath   = "/home/user/alignment/limiric/")
)

GSE1234567 <- limiric(sample_list = sample_list)
```

## Input file format

If your alignment output files are not zipped (end with a ```.gz``` extension), you will need to find a way to do this 
before using ```limiric```. If you have a ```Linux``` or ```mac``` machine you can open the terminal in the directory where your files are stored and use ```gzip```

<br>

```bash
cd path/to/files
gzip * 
```
> This assumes that your files are the only items inside the directory. As with the standard output of many alignment algorithms such as ```STARsolo``` and ```CellRanger```, each sample should
> be stored in its own directory with the following file naming convention :
>
> path/
> |
> ├── matrix.mtx
> |
> ├── barcodes.tsv
> |
> └── features.tsv

<br>

If you have a ```Windows``` machine, Lord be with you. No, there are many ways to get around it with the simplest probably being to install ```Windows Subsystem for Linux```  [🐧](https://learn.microsoft.com/en-us/windows/wsl/install) that creates a new terminal environment for you with ```Linux``` capabilites. From there, you can do the same as above.

## Downstream

### Adding barcode annotations to pre-existing Seurat object 
```R
# Add output to pre-existing Seurat object
limiric_annotations <- read.csv2("path/to/ProjectName_barcodes.csv",
                      sep = ",",
                      col.names = c("barcodes", "limiric"))

seurat@meta.data$limiric <- limiric_annotations$limiric[match(rownames(seurat@meta.data), limiric_annotations$barcodes)]

# Using pre-existing dimensionality reductions, visualise the limiric annotations
limiric_visual <- DimPlot(seurat, group.by = limiric) 

# Filter cells marked as damaged by limiric
seurat_filtered <- subset(seurat, limiric == "cell")
seurat_filtered$limiric <- NULL # once used, remove the column 
```
