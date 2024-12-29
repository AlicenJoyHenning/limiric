
<br>

<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/inst/extdata/limiric_light.png" alt="limiric_logo" height="140" width="370">
</p>



## Contents
[Description](#description) | [Installation](#installation) | [Basic usage](#quickstart)  | [Extended usage](#extended-usage) | [Output explained](#output-explained)  |  [Input file format](#input-file-format) | [Downstream](#downstream)


## Description

<br>
Single-cell RNA sequencing (scRNA-seq) is a well-established technique in the era of next-generation sequencing. From identifying novel cell types to characterizing treatment responses, its widespread applications are undeniably valuable to the field of molecular and cell biology. The reliability of these applications, however, depend entirely on the quality of upstream pre-processing, with cell-level filtering being an important component.
<br><br>

For droplet-based protocols, _low quality_ cells are those that originate from droplets that contain more than one cell (doublet), no cells (empty droplet), or a damaged cell. While there are many scRNA-seq quality control tools available to identify doublets and empty droplets, few are specialised in identifing damaged cells. Damaged cell detection is more often achieved by setting thresholds for metrics such as average mitochondrial gene expression or feature counts per droplets. These thresholds, even when dynamically calculated, vary across cell types, tissues, treatment conditions, and species, with no established ground truth to define them precisely in any context. As a result, applying them in isolation can easily tip the delicate balance of stringency, excluding many true cells from downstream analysis, and leniency, letting many contaminating damaged cells remain. But searching for true damaged droplets in a sample-specific manner is tedious and not always intuitive, leading to user-defined filtering that lacks reproducibility. 

<br>

```limiric``` is an R package for scRNA-seq cell-level quality control that automates the sample-specific detection of damaged cells in one function. It operates on the biological principle that damaged and healthy cells can be differentiated by the complexity of their mitochondrial and ribosomal gene expression, and that this complexity can be effectively resolved through clustering in lower-dimensional space.  Beyond predicting damaged cells, ```limiric``` estimates and corrects for red blood cell contamination— an uncommon but significant artifact of single-cell isolation protocols for which no automated tool currently exists. If relevant to a user's investigation, ```limiric``` can also be used isolate immune (CD45⁺) populations, a step performed before removing damaged cells.

<br>

Built around the community standard ```Seurat``` suite, ```limiric``` is designed to integrate seamlessly into a user's pre-existing scRNA-seq analysis workflow in R. To streamline this integration, ```limiric``` offers ambient RNA correction (```SoupX```) functionality, a community-recommended first step of pre-processing. However, the main output of ```limiric``` is a standard, two-column comma separated value file (cell_identifier, limiric_annotation) that can easily be used in any analysis workflow. 

<br>

<br>

## Installation
### Prerequisites
```limiric``` makes use of the following packages
* ```cowplot``` Wilke, 2024. [GitHub repo](https://github.com/wilkelab/cowplot)
* ```devtools``` Wickham and Hester, 2022. [GitHub repo](https://github.com/r-lib/devtools)
* ```dplyr``` Wickham _et al_, 2023. [GitHub repo](https://github.com/tidyverse/dplyr)
* ```ggplot2``` Wickham, 2016. [GitHub repo](https://github.com/tidyverse/ggplot2)
* ```Matrix``` Bates _et al_, 2024. [Github repo](https://github.com/cran/Matrix)
* ```png``` Urbanek, 2022. [GitHub](https://github.com/cran/png)
* ```Seurat``` Satija and Farrell _et al_, 2015. [Website](https://satijalab.org/seurat/), [GitHub repo](https://github.com/satijalab/seurat)
* ```SoupX``` Young and Behjati, 2020. [GitHub repo](https://github.com/constantAmateur/SoupX)
<br>

Of these, the following need to be made available in your ```R``` environment before ```limiric``` can be installed. If not already, you can install them together with the following code 

```R

packages <- c("cowplot", "devtools", "dplyr", "ggplot2", "Matrix", "png", "Seurat", "SoupX")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }


```

> Please check the associated documentation if problems occur in the installation of any of the prerequisite packages


<br><br>

### Package installation
After all the prerequisites are installed, you can install the latest development version of ```limiric``` from ```GitHub``` using the ```devtools``` package:

```R
devtools::install_github("AlicenJoyHenning/limiric", build_vignettes = TRUE)
```

<br>


### Verify installation 
To ensure ```limiric``` has correctly installed, run the following to see if you can view the package vignette and the function help page. 

```R
library(limiric)
help(package = "limiric")
?limiric()
```
<br>

Additionally, you can perform a test run using a small example dataset stored in the package called ```test_data```. If you need assistance setting this up, please see [Test run](#test-run) for more details. 

<br>
<br>
<br>

## _Quickstart_

### Test run 

Access the available test data (```test_data```) and 
run the ```limiric``` function with the following parameters (ensures test_run is as quick as possible). 

```R

# Load example Seurat object from the limiric package
data("test_data", package = "limiric")

test <- limiric(
     project_name = "test_run",
     filter_rbc   = FALSE,
     seurat_input = test_data,
     output_path  = tempdir()
   )
```

<br>

If anything other than the output below is displayed, please check the installation of prerequisites or open an ```issue``` on this ```GitHub``` page.
```R
# Beginning limiric analysis for test_run ...
# ✔ Seurat object created
# ✔ limiric damaged cell predictions
# ✔ limiric analysis complete.

```
> **Note:** This shouldn't take more than 10 seconds to run but it will vary depending on your machine.


<br> 

#### Basic Usage


1. Detect damaged cells in a sample using the filtered alignment files (```barcodes.tsv.gz```, ```features.tsv.gz```, ```counts.mtx.gz```). All you need to input 
is your sample name, the directory where your alignment files are stored, and a directory where the ```limiric``` output can be created.  
```R
SRR1234567 <- limiric(
    project_name  = "SRR1234567",
    filtered_path = "/home/user/alignment/SRR1234567/filtered/",
    output_path   = "/home/user/alignment/limiric/"
)
```  
> **NB** Please keep the following in mind:
> * Your files must be zipped, see [Input file format](#input-file-format) for more
> * The above usage assumes the sample is of human origin (See parameter: ```organism```)
> * When compiling the data, all genes will be retained (See parameter: ```min_cells```)
> * Red blood cells will automatically be removed from the sample before damaged cell detection (See parameter: ```filter_rbc```)
> * Your output ```Seurat``` object will be filtered, containing only undamaged cells (See parameter: ```filter_output```)
>

<br>
<br>

2. Alternatively, you can use a ```Seurat``` object as input
```R
SRR1234567 <- limiric(
    project_name  = "SRR1234567",
    seurat_input  = seurat_object,
    output_path   = "/home/user/alignment/limiric/"
)
```

<br>
<br>
  
3. If you have more than one sample, you may want to define a sample list
instead of running each individually.
You can do this through the ```sample_list``` parameter.
```R
sample_list <- list(

    list(project_name  = "SRR1234567",
         filtered_path = "/home/user/alignment/SRR1234567/filtered/",
         output_path   = "/home/user/alignment/limiric/"),
    
    list(project_name  = "SRR1234568",
         filtered_path = "/home/user/alignment/SRR1234568/filtered/",
         output_path   = "/home/user/alignment/limiric/"),
    
    list(project_name  = "SRR1234569",
         filtered_path = "/home/user/alignment/SRR1234569/filtered/",
         output_path   = "/home/user/alignment/limiric/")
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


<br><br>


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

The ```limiric``` output will be created in the ```output_path```directory with the following structure  

```
output_path/
|
├── RBCQC
|
├── CellQC
|
└── Filtered
```

<br>

### **_RBCQC_** 

Before damaged cells can be identified, ```limiric``` removes red blood cells, or cells that are highly contaminated with haemoglobin, from your data. This is done under the assumption that red blood cells will not be informative to your study. If this assumption should not be true, you can avoid this filtering using ```filter_rbc = FALSE```.
> **NB** Given many scRNA-seq protocols, such as the ```10X Genomics``` protocol, advise for globin treatment, we hope for the contamination percentage to be as low as possible.

<br>
<table>
  <tr>
    <td>
      The output here comes in the form of a scatter plot where the red blood cells are coloured in blue. The percentage of the total cells that were removed is also given in the plot. 
    </td>
    <td>
      <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/inst/extdata/RBCQC.png" alt="Scatter plot" style="float: right; margin-left: 200px;" width="500">
    </td>
  </tr>
</table>

<br>
<br>

### **_CellQC_** 

This directory contains core diagnostic plots for the ```limiric``` damaged cell detection algorithm. As shown below, you will see four tSNE plots. In this example dataset, you can see the cells divide into two distinct clusters; but which contain _damaged_ cells and which contain _healthy_ cells? 

<br>

From literature we know that damaged cells have high mitochondrial expression and low ribosomal expression. This is largely because damaged cells, such as those undergoing apoptosis, are characterised by compromised cell membranes. In this case, free cytoplasmic RNA- like ribosomal RNA- is more likely to have escaped before being captured and sequenced, meaning its expression will be low. While the RNA enclosed within the mitochondria, which has its own membrane, is retained, meaning its expression will be high. But many other factors likely contribute to this phenomenon, including impaired ribosomal biogenesis in damaged cells. 

<br>

Now answering the question of which cluster is which is easy with the <code>limiric</code> tSNEs where the mitochondrial and ribosomal gene expression of each cell is made visible. 

<br>

<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/inst/extdata/CellQC.png" style="float: right; margin-left: 200px;" width="500">
</p>

<br>

- A `Mitochondrial gene expression` tSNE shows that cells in the smaller cluster express mitochondrial genes at a very high level, suggesting they are likely damaged.

<br>

- The same conclusion can be reached looking at the `Ribosomal gene expression` tSNE where cells in the smaller cluster express ribosomal genes at a very low level.

<br>

- The `Complexity score` measures the total number of mitochondrial and ribosomal genes expressed in a cell. Cells that express a large number of these genes are more likely to be metabolically active, undamaged cells while those that only express a few are likely not. In the `Complexity score` tSNE, the smaller cluster contains cells with very low complexities. This verifies by another metric that these cells are likely damaged.

<br><br>

Together, these metrics are used by the ```limiric``` algorithm to annotate the _damaged_ cells, as shown in the fourth and final tSNE. In summary, by looking at the ```limiric``` tSNE plots, we can immediately tell that we have a small cluster of cells that are likley damaged and should be removed from the sample. This makes sense as, unless something went drastically wrong during sample preparations, you expect there to be far more healthy than damaged cells in your scRNA-seq data.

<br><br>


### **_Filtered_** 
<br>

The main output of the ```limiric``` function comes in the form of a ```barcodes.csv``` containing annotations for each cell barcode of the input data. This can easily be incorporated into an existing scRNA-seq analysis workflow for filtering (see [example](#downstream)).

#### project_name_barcodes.csv  

|     Barcode      | Limiric |
|:---------------:|:-------:|
| AAACCTGAGATAGGAG|  cell   |
| AAACCTGAGCTATGCT|  cell   |
| AAACCTGAGCTGTTCA| damaged |
| AAACCTGCACATTAGC| damaged |
| ... | ... |

<br>
<br>

Additionally, ```limiric``` will output a ```Seurat``` object with ready-filtered barcodes for seamless integration at the start of a ```Seurat``` pre-processing workflow. 

<br>
<br>

## Extended usage
#### _Slightly-less_ Basic Usage

##### 1. Perform ambient RNA correction
Standard good practice highly advises to correct for ambient RNA in your scRNAseq data. If you haven't already performed some kind of correction,
the ```SoupX``` parameter allows you to do so. This will occur before anything else in the ```limiric``` workflow.

```R
SRR1234567 <- limiric(
    project_name  = "SRR1234567",
    filtered_path = "/home/user/alignment/SRR1234567/filtered/",
    soupx        = TRUE,
    raw_path      = "/home/user/alignment/SRR1234567/raw/",
    output_path   = "/home/user/alignment/limiric/"
)
```

<br>

##### 2. Isolate immune cells
If you have a sample where only the immune cells are of interest, include the following ```isolate_cd45```
parameter. This will isolate the immune cells present in the sample, then identify damaged cells.

```R
SRR1234567 <- limiric(
    project_name  = "SRR1234567",
    filtered_path = "/home/user/alignment/SRR1234567/filtered/",
    isolate_cd45  = TRUE,
    output_path   = "/home/user/alignment/limiric/"
)
```

<br>

> **NB** This will change your output directory structure by adding a new ```IMCQC``` layer
>
> ```
> output_path/
> ├── CellQC
> |
> ├── IMCQC
> |
> ├── RBCQC
> |
> └── Filtered
> ```

<table>
  <tr>
    <td>
      This output, like RBCQC, contains a scatter plot showing the removed cell in blue and the retained cells in grey. 
    </td>
    <td>
      <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/inst/extdata/IMCQC.png" alt="Scatter plot" style="float: right; margin-left: 200px;" width="500">
    </td>
  </tr>
</table>

<br>
<br>

##### 3. Combine previous condition
Perform ambient RNA correction with ```SoupX```, filter red blood cells, isolate immune cells and detect damaged cells.

```R
SRR1234567 <- limiric(
    project_name  = "SRR1234567",
    filtered_path = "/home/user/alignment/SRR1234567/filtered/",
    soupx         = TRUE,
    raw_path      = "/home/user/alignment/SRR1234567/raw/",
    isolate_cd45  = TRUE,
    velocyto_path = "/home/user/alignment/velocyto/",
    output_path   = "/home/user/alignment/limiric/"
)
```
> **NB** This will mean your output directory looks like this
>
> ```
> output_path/
>
> ├── CellQC
> |
> ├── IMCQC
> | 
> ├── RBCQC
> | 
> └── Filtered
> ```

<br>
<br>


##### 4. Process multiple samples with the same conditions as above

```R
sample_list <- list(
    list(project_name  = "SRR1234567",
         filtered_path = "/home/user/alignment/SRR1234567/filtered/",
         soupx         = TRUE,
         raw_path      = "/home/user/alignment/SRR1234567/raw/",
         isolate_cd45  = TRUE,
         output_path   = "/home/user/alignment/limiric/"),
    
    list(project_name  = "SRR1234568",
         filtered_path = "/home/user/alignment/SRR1234568/filtered/",
         soupx         = TRUE,
         raw_path      = "/home/user/alignment/SRR1234568/raw/",
         isolate_cd45  = TRUE,
         output_path   = "/home/user/alignment/limiric/"),
    
    list(project_name  = "SRR1234569",
         filtered_path = "/home/user/alignment/SRR1234569/filtered/",
         soupx         = TRUE,
         raw_path      = "/home/user/alignment/SRR1234569/raw/",
         isolate_cd45  = TRUE,
         output_path   = "/home/user/alignment/limiric/")
)

GSE1234567 <- limiric(sample_list = sample_list)
```

<br>
<br>



## Input file format

If your alignment output files are not zipped (end with a ```.gz``` extension), you will need to find a way to do this 
before using ```limiric```. If you have a ```Linux``` or ```mac``` machine you can open the terminal in the directory where your files are stored and use ```gzip```

<br>

```bash
cd path/to/directory
gzip * 
```
> This assumes that your files are the only items inside the directory. As with the standard output of many alignment algorithms such as ```STARsolo``` and ```CellRanger```, each sample should
> be stored in its own directory with the following file naming convention :
> ```
> path/
> |
> ├── matrix.mtx
> |
> ├── barcodes.tsv
> |
> └── features.tsv
> ```

<br>

If you have a ```Windows``` machine, the simplest solution is to install ```Windows Subsystem for Linux```  [🐧](https://learn.microsoft.com/en-us/windows/wsl/install) which creates a new terminal environment for you with ```Linux``` capabilites. From there, you can do the same as above. Note that ```WSL``` uses a different path structure to access file systems where ```Windows``` directories are mounted under ```/mnt/``` 

```bash
cd /mnt/c/Users/path/to/directory
gzip *

```

<br>
<br>



## Downstream

### Adding barcode annotations to pre-existing Seurat object 
```R
# Add output to pre-existing Seurat object
limiric_annotations <- read.csv2("path/to/project_name_barcodes.csv",
                      sep = ",",
                      col.names = c("barcodes", "limiric"))

seurat@meta.data$limiric <- limiric_annotations$limiric[match(rownames(seurat@meta.data), limiric_annotations$barcodes)]

# Using pre-existing dimensionality reductions, visualise the limiric annotations
limiric_visual <- DimPlot(seurat, group.by = limiric) 

# Filter cells marked as damaged by limiric
seurat_filtered <- subset(seurat, limiric == "cell")
seurat_filtered$limiric <- NULL # once used, remove the column 
```

