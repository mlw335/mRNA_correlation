Gene Correlation Explorer

An interactive Shiny application for exploring gene–gene correlation structure in Escherichia coli using a large RNA-seq compendium.

This tool enables users to query a gene of interest, identify highly correlated partners, and visualize correlation neighborhoods as an interactive heatmap with downloadable results.

------------------------------------------------------------

Data Source

This application relies on data from:

Tjaden, B. (2023).
Escherichia coli transcriptome assembly from a compendium of RNA-seq data sets.
RNA Biology, 20(1), 77–84.
https://doi.org/10.1080/15476286.2023.2189331

This tool would not be possible without the above work.

------------------------------------------------------------

System Requirements

R version 4.2 or higher

RStudio (strongly recommended)

Internet connection (first run only, for package installation)

Supported operating systems:
- Windows
- macOS
- Linux

------------------------------------------------------------

Installation and Setup

1. Download the application

Download or clone this repository to your local machine.

2. Download required data files

Due to size constraints, large data files are not included in the repository.

Download the required files from the following Box folder:

Temporary data download (Box):
https://cornell.box.com/s/xex2urj4ot6zmworcaq63vn4tbnp5qw9

After downloading, place the files into the data directory so that it contains:

data/
  correlation_matrix_all_genes.rds
  Supplementary_Table_2.txt

Note: These data are provided temporarily via Box and will be archived on a permanent repository (for example, Zenodo) upon publication.

3. Open the project

Double-click the file:
GeneCorrelationExplorer.Rproj

This ensures the correct working directory.

4. Launch the application

Open any script in RStudio and click "Run App".

The graphical user interface will automatically open in your web browser.

------------------------------------------------------------

Using the Application

1. Enter a gene of interest (for example, xxxX or xxx)

2. Adjust filtering parameters:
   - Correlation threshold
   - Minimum overlap (n)
   - Number of top correlated genes
   - Optional: additional genes of interest

3. Click "Run Analysis"

------------------------------------------------------------

Output

The application generates:

1. An interactive correlation heatmap
   - Zoomable and pannable
   - Downloadable as SVG

2. A table of correlated genes
   - Filterable and sortable
   - Downloadable as CSV

These outputs can be saved for downstream analysis or publication.

------------------------------------------------------------

Notes

Required R packages are installed automatically on first run if missing.

Large correlation matrices are loaded locally; sufficient RAM is recommended.

For best performance, run the app from the project root via the RStudio project file.

------------------------------------------------------------

Questions and Issues

For questions, issues, or suggestions, please contact the authors or open an issue in the repository.
