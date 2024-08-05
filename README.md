<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/limiric.png" alt="limiric_logo" height="150" />
</p>

## Contents
[Description](#description)
[Installation](#installation)
[Quickstart](#quickstart)
[Basic Usage](#basic-usage)
[Output explained](#output-explained)
[More Information](#more-information)



## Description

Single cell RNA sequencing quality control package for sample-specific damaged cell detection through low dimension mitochondrial and ribosomal cluster selection.

## Installation
### Prerequisites
The ```limiric``` package requires the following packages to be installed in your ```R``` environment  
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

required_packages <- c("cowplot", "dplyr", "ggplot2", "Matrix", "png", "Seurat", "SoupX")

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
<i>Please check the associated documentation if problems occur in the installation of any of the prerequisite packages.</i>
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


1. Detect damaged cells in a sample using the filtered alignment files (```barcodes.tsv.gz```, ```features.tsv.gz```, ```counts.mtx.gz```) as input. All you need to specify 
is a character string related to your sample name, the directory where your alignment files are stored, and a directory where the ```limiric``` output can be created.  
```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```  
> **NB** Please keep the following in mind: 
> * The above usage assumes the sample is of human origin (See parameter: ```Organism```)
> * When compiling the data, all genes will be retained (See parameter: ```MinCells```)
> * Red blood cells will automatically be removed from the sample before damaged cell detection (See parameter: ```FilterRBC```)
> * Your output ```Seurat``` object will be filtered, containing only undamaged cells (See parameter: ```FilterOutput```)
>
> 

<br>
<br>

2. Alternatively, you can use a ```Seurat``` object as input and your output will be identical to the previous example.
```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    SeuratInput  = seurat_object,
    OutputPath   = "/home/user/alignment/limiric/"
)
```

<br>
<br>
  
3. If you have more than one sample, it may be easier for you to define a sample list
instead of running each individually. You can do this through the ```sample_list``` parameter.
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
>    list("SRR1234567","/home/user/alignment/SRR1234567/filtered/","/home/user/alignment/limiric/"),
>    list("SRR1234568","/home/user/alignment/SRR1234568/filtered/","/home/user/alignment/limiric/"),
>    list("SRR1234569", "/home/user/alignment/SRR1234569/filtered/","/home/user/alignment/limiric/")
>  )
>``` 
> Please try to avoid giving ```limiric``` an existential crisis
> 
<br>

## Output explained   

Before interpretting ```limiric```'s outputs, it may be helpful to understand a bit of what is going on in the background.
<br><br>
As you know, the magic of scRNA-seq revolves around the identification of distinct cell populations. As you also know, distinct populations are identified according to the
clusters that form when dimensionality reduction is performed on the most variable genes in a sample- for example, the typical ```Seurat``` workflow uses the top two or three thousand variable genes.  

In our case, instead of trying to isolate _distinct_ populations, such lymphoid and myeloid cells in immune-related datasets, we are only interested in two populations: _damagaged_ cells and _healthy_ cells.
In general, most of what we know about how these two populations differ can be summarized by looking at the expression of mitochondrial and ribosomal genes. Thus, if we perform dimensionality reduction and clustering 
using _only_ these genes, we expect to see a division of cells into two clusters associated with the cell's status as either _damaged_ or _healthy_.  
This is the basis of the ```low dimension mitochondrial & ribosomal clustering```, aka ```limiric```, algorithm. 

Now lets look at the specifics of the ```limiric``` output that will be created in the directory you specified with the following core structure

```
OutputPath/
├── CellQC
├── RBCQC
└── Filtered
```

<br>

### **_RBCQC_** 

Before damaged cells can be identified, it makes sense to remove red blood cells or cells that are highly contaminated with haemoglobin since they will not be informative to your study, 
unless red blood cells would be important to your study in which case you can avoid this using ```FilterRBC = FALSE```.
<br>

The output here comes in the form of a scatter plot where the red blood cells are coloured in blue. The percentage of the total cells that were removed is also given in the plot. 

<br><br>
<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/RBCQC.png" height="500" />
</p>

> **NB** Given many scRNA-seq protocols, such as the ```10X Genomics``` protocol, advise for globin treatment to prevent this, we hope for the contamination percentage to be as low as possible.
<br>

### **_CellQC_** 

This directory contains core diagnostic plots for the ```limiric``` damaged cell detection algorithm. As shown below, you will see four tSNE dimensionally reduced representations of your cells.
<br><br>
In this example dataset, _damaged_ and _healthy_ cells divide into two distinct clusters. But which is which? From literature we know that _damaged_ cells have **high** mitochondrial expression and **low** ribosomal expression.
This is most likely because _damaged_ cells, such as those undergoing apoptosis, are characterised by compromised cell membranes where free cytoplasmic RNA, such as ribosomal RNA, is more likely to have escaped before being captured and sequenced while the RNA enclosed within the mitochondria is retained. Ribosomal biogenesis is also said to be impaired in _damaged_ cells, leading to reduced ribosomal expression. 

<br>

From the example below, we can see the smaller cluster is more blue in the ```Mitochondrial gene expression``` tSNE, related to **high** mitochondrial expression, 
as well as more grey in the ```Ribosomal gene expression``` tSNE, related to **low** ribosomal expression. This tells us that the smaller cluster represents the _damaged_ cells- as we would hope given they will be removed!  

<br>

This can be verified by looking at the ```Complexity score``` tSNE where cells that show expression of a higher number of genes present in the mitochondrial and ribosomal gene set are more likely to be metabolically active, _healthy_ cells. As you can in this plot, the _healthy_ cells are a darker blue, related to high complexity, while the _damaged_ cells are a light grey, related to low complexity. Together, these metrics are used by the ```limiric``` algorithm to annotate the _damaged_ cells, as shown in the fourth and final tSNE.
<br>
<br>
<br>
<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/CellQC.png" height="500" />
</p>
<br>
<br>



### **_Filtered_** 
<br>

This directory houses the filtered ```Seurat``` object that can be read into your ```R``` environment again for further use
```R
ProjectName <- readRDS("OutputPath/Filtered/ProjectName.rds")
````

<br>

It also includes a ```.csv``` containing the ```limiric``` annotation of all the cell barcdoes present in the input alignment file that can be re-added to the meta.data slot of your ```Seurat``` object

```R

limiric_annotations <- read.csv2(nf_csv,
                      sep = ",",
                      col.names = c("barcodes", "limiric"))

seurat@meta.data$limiric <- limiric_annotations$limiric[match(rownames(seurat@meta.data), limiric_annotations$barcodes)]

```
<br>

## More Information
#### _Slightly-less_ Basic Usage

1. Perform Ambient RNA Correction Prior to Damaged Cell Detection
Detect damaged cells after performing ambient RNA correction with ```SoupX```

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    SoupX        = TRUE,
    RawPath      = "/home/user/alignment/SRR1234567/raw/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```


2. Isolate Immune Cells and Identify Damaged Cells
First, isolate the immune cells present in the sample, then identify damaged cells.

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

### IMCQC
This output, like RBCQC, contains a scatter plot showing the removed cell in blue and the retained cells in grey.

<br>

<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/IMCQC.png" height="500" />
</p>

<br>
<br>

3. Combine ```limiric``` annotations with ```DropletQC```
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
> **NB** This will change your output directory structure by adding a new ```IMCQC``` layer
>
> ```
> OutputPath/
> ├── CellQC
> ├── DropletQC
> ├── RBCQC
> └── Filtered
> ```
>
### DropletQC
This will output a scatter plot and tSNE showing the cells annotated as _damaged_ by both ```limiric``` and ```DropletQC```.

<br>

<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/DropletQC.png" height="500" />
</p>

<br>

4. Combine previous conditions
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


