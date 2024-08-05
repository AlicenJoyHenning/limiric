<p align="center">
  <img src="https://github.com/AlicenJoyHenning/limiric/blob/master/images/limiric.png" alt="limiric_logo" height="180" />
</p>
<hr style="border: 1px solid black; width: 50%;"/>

## _Description_

Single cell RNA sequencing quality control package for sample-specific damaged cell detection through low dimension mitochondrial and ribosomal cluster selection.

## Installation 

Install latest development version by running the following 

```R
devtools::install_github("AlicenJoyHenning/limiric")
```

## Quickstart 

#### 1. Basic Usage
Detect damaged cells in a single sample using the filtered alignment files.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```



#### 2. Perform Ambient RNA Correction Prior to Damaged Cell Detection
Detect damaged cells after performing ambient RNA correction.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    SoupX        = TRUE,
    RawPath      = "/home/user/alignment/SRR1234567/raw/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```



#### 3. Isolate Immune Cells and Identify Damaged Cells
First, isolate the immune cells present in the sample, then identify damaged cells.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    IsolateCD45  = TRUE,
    OutputPath   = "/home/user/alignment/limiric/"
)
```



#### 4. Combine limiric Damaged Cell Annotations with DropletQC
Detect damaged cells and compare results with those from DropletQC.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    FilteredPath = "/home/user/alignment/SRR1234567/filtered/",
    DropletQC    = TRUE,
    VelocytoPath = "/home/user/alignment/velocyto/",
    OutputPath   = "/home/user/alignment/limiric/"
)
```



#### 5. Combine All Previous Conditions
Perform ambient RNA correction with SoupX, filter red blood cells, isolate immune cells, detect damaged cells, and compare against DropletQC.

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



#### 6. Process Multiple Samples with the Same Conditions as in 5
Define a list of samples and run the function with this list.

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



#### 7. Use a Seurat Object as Input
Alternatively, use a Seurat object as input.

```R
SRR1234567 <- limiric(
    ProjectName  = "SRR1234567",
    SeuratInput  = seurat_object,
    OutputPath   = "/home/user/alignment/limiric/"
)
```

## Output examples 

The pacakge will 
