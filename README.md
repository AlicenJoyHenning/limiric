---
editor_options: 
  markdown: 
    wrap: 72
---

![DamageDetective](docs/damagedetective.svg)

## Contents

[Description](#description) \| [Installation](#installation) \| [Basic
usage](#quickstart) \| [Extended usage](#extended-usage) \| [Output
explained](#output-explained) \| [Input file format](#input-file-format)
\| [Downstream](#downstream)

## Description {#description}

<br> Single-cell RNA sequencing (scRNA-seq) is a well-established
technique in the era of next-generation sequencing. From identifying
novel cell types to characterizing treatment responses, its widespread
applications are undeniably valuable to the field of molecular and cell
biology. The reliability of these applications, however, depend entirely
on the quality of upstream pre-processing, with cell-level filtering
being an important component. <br><br>

For droplet-based protocols, *low quality* cells are those that
originate from droplets that contain more than one cell (doublet), no
cells (empty droplet), or a damaged cell. While there are many scRNA-seq
quality control tools available to identify doublets and empty droplets,
few are specialised in identifing damaged cells. Damaged cell detection
is more often achieved by setting thresholds for metrics such as average
mitochondrial gene expression or UMI counts per droplets. These
thresholds, even when dynamically calculated, vary across cell types,
tissues, treatment conditions, and species, with no established ground
truth to define them precisely in any context. As a result, applying
them in isolation can easily tip the delicate balance of stringency,
excluding many true cells from downstream analysis, and leniency,
letting many contaminating damaged cells remain. But searching for true
damaged droplets in a sample-specific manner is tedious and not always
intuitive, leading to user-defined filtering that lacks reproducibility.

<br>

`DamageDetective` is an R package for scRNA-seq cell-level quality control that
automates the sample-specific detection of damaged cells in one
function. It operates on the biological principle that damaged and
healthy cells can be differentiated by the complexity of their
mitochondrial and ribosomal gene expression, and that this complexity
can be effectively resolved through clustering in lower-dimensional
space. Beyond predicting damaged cells, `DamageDetective` estimates and corrects
for red blood cell contamination— an uncommon but significant artifact
of single-cell isolation protocols for which no automated tool currently
exists. If relevant to a user's investigation, `DamageDetective` can also be
used isolate immune (CD45⁺) populations, a step performed before
removing damaged cells.

<br>

Built around the community standard `Seurat` suite, `DamageDetective` is
designed to integrate seamlessly into a user's pre-existing scRNA-seq
analysis workflow in R. To streamline this integration, `DamageDetective` offers
ambient RNA correction (`SoupX`) functionality, a community-recommended
first step of pre-processing. However, the main output of `DamageDetective` is a
standard, two-column comma separated value file (cell_identifier,
DamageDetective_annotation) that can easily be used in any analysis workflow.

<br>

<br>

## Installation {#installation}

### Prerequisites

`DamageDetective` makes use of the following packages \* `cowplot` Wilke, 2024.
[GitHub repo](https://github.com/wilkelab/cowplot) \* `devtools` Wickham
and Hester, 2022. [GitHub repo](https://github.com/r-lib/devtools) \*
`dplyr` Wickham *et al*, 2023. [GitHub
repo](https://github.com/tidyverse/dplyr) \* `ggplot2` Wickham, 2016.
[GitHub repo](https://github.com/tidyverse/ggplot2) \* `Matrix` Bates
*et al*, 2024. [Github repo](https://github.com/cran/Matrix) \* `png`
Urbanek, 2022. [GitHub](https://github.com/cran/png) \* `Seurat` Satija
and Farrell *et al*, 2015. [Website](https://satijalab.org/seurat/),
[GitHub repo](https://github.com/satijalab/seurat) \* `SoupX` Young and
Behjati, 2020. [GitHub repo](https://github.com/constantAmateur/SoupX)
<br>

Of these, the following need to be made available in your `R`
environment before `DamageDetective` can be installed. If not already, you can
install them together with the following code

``` r

packages <- c("cowplot", "devtools", "dplyr", "ggplot2", "Matrix", "png", "Seurat", "SoupX")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
```

> Please check the associated documentation if problems occur in the
> installation of any of the prerequisite packages

<br><br>

### Package installation

After all the prerequisites are installed, you can install the latest
development version of `DamageDetective` from `GitHub` using the `devtools`
package:

``` r
devtools::install_github("AlicenJoyHenning/DamageDetective", build_vignettes = TRUE)
```

<br>

### Verify installation

To ensure `DamageDetective` has correctly installed, run the following to see if
you can view the package vignette and the function help page.

``` r
library(DamageDetective)
help(package = "DamageDetective")
?DamageDetective()
```

<br>

Additionally, you can perform a test run using a small example dataset
stored in the package called `test_data`. If you need assistance setting
this up, please see [Test run](#test-run) for more details.

<br> <br> <br>

## *Quickstart* {#quickstart}

### Test run {#test-run}

Access the available test data (`test_data`) and run the `DamageDetective`
function with the following parameters (ensures test_run is as quick as
possible).

``` r

# Load example Seurat object from the DamageDetective package
data("test_data", package = "DamageDetective")

test <- DamageDetective(
     project_name = "test_run",
     filter_rbc   = FALSE,
     seurat_input = test_data,
     output_path  = tempdir()
   )
```

<br>

If anything other than the output below is displayed, please check the
installation of prerequisites or open an `issue` on this `GitHub` page.

``` r
# Beginning DamageDetective analysis for test_run ...
# ✔ Seurat object created
# ✔ DamageDetective damaged cell predictions
# ✔ DamageDetective analysis complete.
```

> **Note:** This shouldn't take more than 10 seconds to run but it will
> vary depending on your machine.

<br>

#### Basic Usage

1.  Detect damaged cells in a sample using the filtered alignment files
    (`barcodes.tsv.gz`, `features.tsv.gz`, `counts.mtx.gz`). All you
    need to input is your sample name, the directory where your
    alignment files are stored, and a directory where the `DamageDetective`
    output can be created.\

``` r
SRR1234567 <- DamageDetective(
    project_name  = "SRR1234567",
    filtered_path = "/home/user/alignment/SRR1234567/filtered/",
    output_path   = "/home/user/alignment/DamageDetective/"
)
```

> **NB** Please keep the following in mind: \* Your files must be
> zipped, see [Input file format](#input-file-format) for more \* The
> above usage assumes the sample is of human origin (See parameter:
> `organism`) \* When compiling the data, all genes will be retained
> (See parameter: `min_cells`) \* Red blood cells will automatically be
> removed from the sample before damaged cell detection (See parameter:
> `filter_rbc`) \* Your output `Seurat` object will be filtered,
> containing only undamaged cells (See parameter: `filter_output`)

<br> <br>

2.  Alternatively, you can use a `Seurat` object as input

``` r
SRR1234567 <- DamageDetective(
    project_name  = "SRR1234567",
    seurat_input  = seurat_object,
    output_path   = "/home/user/alignment/DamageDetective/"
)
```

<br> <br>

3.  If you have more than one sample, you may want to define a sample
    list instead of running each individually. You can do this through
    the `sample_list` parameter.

``` r
sample_list <- list(

    list(project_name  = "SRR1234567",
         filtered_path = "/home/user/alignment/SRR1234567/filtered/",
         output_path   = "/home/user/alignment/DamageDetective/"),
    
    list(project_name  = "SRR1234568",
         filtered_path = "/home/user/alignment/SRR1234568/filtered/",
         output_path   = "/home/user/alignment/DamageDetective/"),
    
    list(project_name  = "SRR1234569",
         filtered_path = "/home/user/alignment/SRR1234569/filtered/",
         output_path   = "/home/user/alignment/DamageDetective/")
)

GSE1234567 <- DamageDetective(sample_list = sample_list)
```

> **NB** Please note that you must define your sample list with each
> parameter specified, i.e. `DamageDetective` won't know what to do if your list
> looks like this :
>
> ``` r
> sample_list <- list(
>
>    list("SRR1234567",
>         "/home/user/alignment/SRR1234567/filtered/",
>         "/home/user/alignment/DamageDetective/"),
>    list("SRR1234568",
>         "/home/user/alignment/SRR1234568/filtered/",
>         "/home/user/alignment/DamageDetective/"),
>    list("SRR1234569",
>         "/home/user/alignment/SRR1234569/filtered/",
>         "/home/user/alignment/DamageDetective/")
>  )
> ```
>
> <br>

<br><br>

## Output explained {#output-explained}

Before interpretting `DamageDetective`'s outputs, it may be helpful to
understand a bit of what's going on in the background. <br><br> As you
know, the magic of scRNA-seq revolves around the identification of
distinct cell populations. In majority of scRNA-seq analysis workflows,
distinct populations are identified according to the clusters that form
when dimensionality reduction is performed on the most variable genes in
a sample. For instance, the typical `Seurat` workflow uses the top two
or three thousand variable genes.

In our case, instead of trying to isolate distinct cell-type associated
populations, we are only interested in two broad populations: *damaged*
and *healthy* cells. In general, most of what we know about how these
two populations differ can be summarized by looking at the expression of
mitochondrial and ribosomal genes. Thus, if we perform dimensionality
reduction and clustering using *only* these genes, we expect to see a
division of cells into two clusters associated with the cell's status as
either *damaged* or *healthy*.

This is the basis of the
`low dimension mitochondrial & ribosomal clustering`, aka `DamageDetective`,
algorithm.

### Output directory

The `DamageDetective` output will be created in the `output_path`directory with
the following structure

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

### ***RBCQC***

Before damaged cells can be identified, `DamageDetective` removes red blood
cells, or cells that are highly contaminated with haemoglobin, from your
data. This is done under the assumption that red blood cells will not be
informative to your study. If this assumption should not be true, you
can avoid this filtering using `filter_rbc = FALSE`. \> **NB** Given
many scRNA-seq protocols, such as the `10X Genomics` protocol, advise
for globin treatment, we hope for the contamination percentage to be as
low as possible.

<br>
31b8e172-b470-440e-83d8-e6b185028602:dAB5AHAAZQA6AE8AQQBCAGwAQQBHAFkAQQBOAFEAQgBoAEEARABjAEEATgB3AEEAeQBBAEMAMABBAFoAQQBCAGsAQQBEAFkAQQBNAHcAQQB0AEEARABRAEEATgBnAEEAeQBBAEQASQBBAEwAUQBBADQAQQBEAFEAQQBZAGcAQgBtAEEAQwAwAEEAWQBRAEIAbQBBAEQARQBBAE8AUQBBADUAQQBEAFUAQQBZAFEAQQB4AEEARwBJAEEATQBnAEIAaQBBAEQAawBBAAoAcABvAHMAaQB0AGkAbwBuADoATQBRAEEAdwBBAEQAYwBBAE0AdwBBADQAQQBBAD0APQAKAHAAcgBlAGYAaQB4ADoACgBzAG8AdQByAGMAZQA6AFAAQQBCADAAQQBHAEUAQQBZAGcAQgBzAEEARwBVAEEAUABnAEEASwBBAEMAQQBBAEkAQQBBADgAQQBIAFEAQQBjAGcAQQArAEEAQQBvAEEASQBBAEEAZwBBAEMAQQBBAEkAQQBBADgAQQBIAFEAQQBaAEEAQQArAEEAQQBvAEEASQBBAEEAZwBBAEMAQQBBAEkAQQBBAGcAQQBDAEEAQQBWAEEAQgBvAEEARwBVAEEASQBBAEIAdgBBAEgAVQBBAGQAQQBCAHcAQQBIAFUAQQBkAEEAQQBnAEEARwBnAEEAWgBRAEIAeQBBAEcAVQBBAEkAQQBCAGoAQQBHADgAQQBiAFEAQgBsAEEASABNAEEASQBBAEIAcABBAEcANABBAEkAQQBCADAAQQBHAGcAQQBaAFEAQQBnAEEARwBZAEEAYgB3AEIAeQBBAEcAMABBAEkAQQBCAHYAQQBHAFkAQQBJAEEAQgBoAEEAQwBBAEEAYwB3AEIAagBBAEcARQBBAGQAQQBCADAAQQBHAFUAQQBjAGcAQQBnAEEASABBAEEAYgBBAEIAdgBBAEgAUQBBAEkAQQBCADMAQQBHAGcAQQBaAFEAQgB5AEEARwBVAEEASQBBAEIAMABBAEcAZwBBAFoAUQBBAGcAQQBIAEkAQQBaAFEAQgBrAEEAQwBBAEEAWQBnAEIAcwBBAEcAOABBAGIAdwBCAGsAQQBDAEEAQQBZAHcAQgBsAEEARwB3AEEAYgBBAEIAegBBAEMAQQBBAFkAUQBCAHkAQQBHAFUAQQBJAEEAQgBqAEEARwA4AEEAYgBBAEIAdgBBAEgAVQBBAGMAZwBCAGwAQQBHAFEAQQBJAEEAQgBwAEEARwA0AEEASQBBAEIAaQBBAEcAdwBBAGQAUQBCAGwAQQBDADQAQQBJAEEAQgBVAEEARwBnAEEAWgBRAEEAZwBBAEgAQQBBAFoAUQBCAHkAQQBHAE0AQQBaAFEAQgB1AEEASABRAEEAWQBRAEIAbgBBAEcAVQBBAEkAQQBCAHYAQQBHAFkAQQBJAEEAQgAwAEEARwBnAEEAWgBRAEEAZwBBAEgAUQBBAGIAdwBCADAAQQBHAEUAQQBiAEEAQQBnAEEARwBNAEEAWgBRAEIAcwBBAEcAdwBBAGMAdwBBAGcAQQBIAFEAQQBhAEEAQgBoAEEASABRAEEASQBBAEIAMwBBAEcAVQBBAGMAZwBCAGwAQQBDAEEAQQBjAGcAQgBsAEEARwAwAEEAYgB3AEIAMgBBAEcAVQBBAFoAQQBBAGcAQQBHAGsAQQBjAHcAQQBnAEEARwBFAEEAYgBBAEIAegBBAEcAOABBAEkAQQBCAG4AQQBHAGsAQQBkAGcAQgBsAEEARwA0AEEASQBBAEIAcABBAEcANABBAEkAQQBCADAAQQBHAGcAQQBaAFEAQQBnAEEASABBAEEAYgBBAEIAdgBBAEgAUQBBAEwAZwBBAGcAQQBBAG8AQQBJAEEAQQBnAEEAQwBBAEEASQBBAEEAOABBAEMAOABBAGQAQQBCAGsAQQBEADQAQQBDAGcAQQBnAEEAQwBBAEEASQBBAEEAZwBBAEQAdwBBAGQAQQBCAGsAQQBEADQAQQBDAGcAQQBnAEEAQwBBAEEASQBBAEEAZwBBAEMAQQBBAEkAQQBBADgAQQBHAGsAQQBiAFEAQgBuAEEAQwBBAEEAYwB3AEIAeQBBAEcATQBBAFAAUQBBAGkAQQBHAGcAQQBkAEEAQgAwAEEASABBAEEAYwB3AEEANgBBAEMAOABBAEwAdwBCAG4AQQBHAGsAQQBkAEEAQgBvAEEASABVAEEAWQBnAEEAdQBBAEcATQBBAGIAdwBCAHQAQQBDADgAQQBRAFEAQgBzAEEARwBrAEEAWQB3AEIAbABBAEcANABBAFMAZwBCAHYAQQBIAGsAQQBTAEEAQgBsAEEARwA0AEEAYgBnAEIAcABBAEcANABBAFoAdwBBAHYAQQBHAHcAQQBhAFEAQgB0AEEARwBrAEEAYwBnAEIAcABBAEcATQBBAEwAdwBCAGkAQQBHAHcAQQBiAHcAQgBpAEEAQwA4AEEAYgBRAEIAaABBAEgATQBBAGQAQQBCAGwAQQBIAEkAQQBMAHcAQgBwAEEARwA0AEEAYwB3AEIAMABBAEMAOABBAFoAUQBCADQAQQBIAFEAQQBaAEEAQgBoAEEASABRAEEAWQBRAEEAdgBBAEYASQBBAFEAZwBCAEQAQQBGAEUAQQBRAHcAQQB1AEEASABBAEEAYgBnAEIAbgBBAEMASQBBAEkAQQBCAGgAQQBHAHcAQQBkAEEAQQA5AEEAQwBJAEEAVQB3AEIAagBBAEcARQBBAGQAQQBCADAAQQBHAFUAQQBjAGcAQQBnAEEASABBAEEAYgBBAEIAdgBBAEgAUQBBAEkAZwBBAGcAQQBIAE0AQQBkAEEAQgA1AEEARwB3AEEAWgBRAEEAOQBBAEMASQBBAFoAZwBCAHMAQQBHADgAQQBZAFEAQgAwAEEARABvAEEASQBBAEIAeQBBAEcAawBBAFoAdwBCAG8AQQBIAFEAQQBPAHcAQQBnAEEARwAwAEEAWQBRAEIAeQBBAEcAYwBBAGEAUQBCAHUAQQBDADAAQQBiAEEAQgBsAEEARwBZAEEAZABBAEEANgBBAEMAQQBBAE0AZwBBAHcAQQBEAEEAQQBjAEEAQgA0AEEARABzAEEASQBnAEEAZwBBAEgAYwBBAGEAUQBCAGsAQQBIAFEAQQBhAEEAQQA5AEEAQwBJAEEATgBRAEEAdwBBAEQAQQBBAEkAZwBBACsAQQBBAG8AQQBJAEEAQQBnAEEAQwBBAEEASQBBAEEAOABBAEMAOABBAGQAQQBCAGsAQQBEADQAQQBDAGcAQQBnAEEAQwBBAEEAUABBAEEAdgBBAEgAUQBBAGMAZwBBACsAQQBBAG8AQQBQAEEAQQB2AEEASABRAEEAWQBRAEIAaQBBAEcAdwBBAFoAUQBBACsAQQBBAD0APQAKAHMAdQBmAGYAaQB4ADoA:31b8e172-b470-440e-83d8-e6b185028602

<br> <br>

### ***CellQC***

This directory contains core diagnostic plots for the `DamageDetective` damaged
cell detection algorithm. As shown below, you will see four tSNE plots.
In this example dataset, you can see the cells divide into two distinct
clusters; but which contain *damaged* cells and which contain *healthy*
cells?

<br>

From literature we know that damaged cells have high mitochondrial
expression and low ribosomal expression. This is largely because damaged
cells, such as those undergoing apoptosis, are characterised by
compromised cell membranes. In this case, free cytoplasmic RNA- like
ribosomal RNA- is more likely to have escaped before being captured and
sequenced, meaning its expression will be low. While the RNA enclosed
within the mitochondria, which has its own membrane, is retained,
meaning its expression will be high. But many other factors likely
contribute to this phenomenon, including impaired ribosomal biogenesis
in damaged cells.

<br>

Now answering the question of which cluster is which is easy with the
<code>DamageDetective</code> tSNEs where the mitochondrial and ribosomal gene
expression of each cell is made visible.

<br>

<p align="center">

<img src="https://github.com/AlicenJoyHenning/DamageDetective/blob/master/inst/extdata/CellQC.png" style="float: right; margin-left: 200px;" width="500"/>

</p>

<br>

-   A `Mitochondrial gene expression` tSNE shows that cells in the
    smaller cluster express mitochondrial genes at a very high level,
    suggesting they are likely damaged.

<br>

-   The same conclusion can be reached looking at the
    `Ribosomal gene expression` tSNE where cells in the smaller cluster
    express ribosomal genes at a very low level.

<br>

-   The `Complexity score` measures the total number of mitochondrial
    and ribosomal genes expressed in a cell. Cells that express a large
    number of these genes are more likely to be metabolically active,
    undamaged cells while those that only express a few are likely not.
    In the `Complexity score` tSNE, the smaller cluster contains cells
    with very low complexities. This verifies by another metric that
    these cells are likely damaged.

<br><br>

Together, these metrics are used by the `DamageDetective` algorithm to annotate
the *damaged* cells, as shown in the fourth and final tSNE. In summary,
by looking at the `DamageDetective` tSNE plots, we can immediately tell that we
have a small cluster of cells that are likley damaged and should be
removed from the sample. This makes sense as, unless something went
drastically wrong during sample preparations, you expect there to be far
more healthy than damaged cells in your scRNA-seq data.

<br><br>

### ***Filtered***

<br>

The main output of the `DamageDetective` function comes in the form of a
`barcodes.csv` containing annotations for each cell barcode of the input
data. This can easily be incorporated into an existing scRNA-seq
analysis workflow for filtering (see [example](#downstream)).

#### project_name_barcodes.csv

|     Barcode      | DamageDetective |
|:----------------:|:-------:|
| AAACCTGAGATAGGAG |  cell   |
| AAACCTGAGCTATGCT |  cell   |
| AAACCTGAGCTGTTCA | damaged |
| AAACCTGCACATTAGC | damaged |
|       ...        |   ...   |

<br> <br>

Additionally, `DamageDetective` will output a `Seurat` object with
ready-filtered barcodes for seamless integration at the start of a
`Seurat` pre-processing workflow.

<br> <br>

## Extended usage {#extended-usage}

#### *Slightly-less* Basic Usage

##### 1. Perform ambient RNA correction

Standard good practice highly advises to correct for ambient RNA in your
scRNAseq data. If you haven't already performed some kind of correction,
the `SoupX` parameter allows you to do so. This will occur before
anything else in the `DamageDetective` workflow.

``` r
SRR1234567 <- DamageDetective(
    project_name  = "SRR1234567",
    filtered_path = "/home/user/alignment/SRR1234567/filtered/",
    soupx        = TRUE,
    raw_path      = "/home/user/alignment/SRR1234567/raw/",
    output_path   = "/home/user/alignment/DamageDetective/"
)
```

<br>

##### 2. Isolate immune cells

If you have a sample where only the immune cells are of interest,
include the following `isolate_cd45` parameter. This will isolate the
immune cells present in the sample, then identify damaged cells.

``` r
SRR1234567 <- DamageDetective(
    project_name  = "SRR1234567",
    filtered_path = "/home/user/alignment/SRR1234567/filtered/",
    isolate_cd45  = TRUE,
    output_path   = "/home/user/alignment/DamageDetective/"
)
```

<br>

> **NB** This will change your output directory structure by adding a
> new `IMCQC` layer
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

+---------------------------------+---------------------------------+
| This output, like RBCQC,        | ![Scatter                       |
| contains a scatter plot showing | p                               |
| the removed cell in blue and    | lot](htt%20ps://github.com/Alic |
| the retained cells in grey.     | enJoyHenning%20/DamageDetective/blob/ma |
|                                 | ster/inst/extdat%20a/IMCQC.png) |
+---------------------------------+---------------------------------+

<br> <br>

##### 3. Combine previous condition

Perform ambient RNA correction with `SoupX`, filter red blood cells,
isolate immune cells and detect damaged cells.

``` r
SRR1234567 <- DamageDetective(
    project_name  = "SRR1234567",
    filtered_path = "/home/user/alignment/SRR1234567/filtered/",
    soupx         = TRUE,
    raw_path      = "/home/user/alignment/SRR1234567/raw/",
    isolate_cd45  = TRUE,
    velocyto_path = "/home/user/alignment/velocyto/",
    output_path   = "/home/user/alignment/DamageDetective/"
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

<br> <br>

##### 4. Process multiple samples with the same conditions as above

``` r
sample_list <- list(
    list(project_name  = "SRR1234567",
         filtered_path = "/home/user/alignment/SRR1234567/filtered/",
         soupx         = TRUE,
         raw_path      = "/home/user/alignment/SRR1234567/raw/",
         isolate_cd45  = TRUE,
         output_path   = "/home/user/alignment/DamageDetective/"),
    
    list(project_name  = "SRR1234568",
         filtered_path = "/home/user/alignment/SRR1234568/filtered/",
         soupx         = TRUE,
         raw_path      = "/home/user/alignment/SRR1234568/raw/",
         isolate_cd45  = TRUE,
         output_path   = "/home/user/alignment/DamageDetective/"),
    
    list(project_name  = "SRR1234569",
         filtered_path = "/home/user/alignment/SRR1234569/filtered/",
         soupx         = TRUE,
         raw_path      = "/home/user/alignment/SRR1234569/raw/",
         isolate_cd45  = TRUE,
         output_path   = "/home/user/alignment/DamageDetective/")
)

GSE1234567 <- DamageDetective(sample_list = sample_list)
```

<br> <br>

## Input file format {#input-file-format}

If your alignment output files are not zipped (end with a `.gz`
extension), you will need to find a way to do this before using
`DamageDetective`. If you have a `Linux` or `mac` machine you can open the
terminal in the directory where your files are stored and use `gzip`

<br>

``` bash
cd path/to/directory
gzip * 
```

> This assumes that your files are the only items inside the directory.
> As with the standard output of many alignment algorithms such as
> `STARsolo` and `CellRanger`, each sample should be stored in its own
> directory with the following file naming convention :
>
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

If you have a `Windows` machine, the simplest solution is to install
`Windows Subsystem for Linux`
[🐧](https://learn.microsoft.com/en-us/windows/wsl/install) which
creates a new terminal environment for you with `Linux` capabilites.
From there, you can do the same as above. Note that `WSL` uses a
different path structure to access file systems where `Windows`
directories are mounted under `/mnt/`

``` bash
cd /mnt/c/Users/path/to/directory
gzip *
```

<br> <br>

## Downstream {#downstream}

### Adding barcode annotations to pre-existing Seurat object

``` r
# Add output to pre-existing Seurat object
DamageDetective_annotations <- read.csv2("path/to/project_name_barcodes.csv",
                      sep = ",",
                      col.names = c("barcodes", "DamageDetective"))

seurat@meta.data$DamageDetective <- DamageDetective_annotations$DamageDetective[match(rownames(seurat@meta.data), DamageDetective_annotations$barcodes)]

# Using pre-existing dimensionality reductions, visualise the DamageDetective annotations
DamageDetective_visual <- DimPlot(seurat, group.by = DamageDetective) 

# Filter cells marked as damaged by DamageDetective
seurat_filtered <- subset(seurat, DamageDetective == "cell")
seurat_filtered$DamageDetective <- NULL # once used, remove the column 
```
