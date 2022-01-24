### __Introduction__

TeaProt is an online Shiny tool that integrates upstream transcription factor enrichment analysis with downstream pathway analysis through an easy-to-use interactive interface. TeaProt maps user’s omics data with online databases to provide a collection of annotations on drug-gene interactions, subcellular localizations, phenotypic functions, gene-disease associations and enzyme-gene interactions, usefull for further analyses. Users can combine TeaProt and urPTMdb for a novel and easy-to-use online proteomics/transcriptomics analysis pipeline featuring novel underrepresented genesets to allow the discovery of downstream cellular processes, upstream transcriptional regulation and classes of PTMs potentially regulated by a users’ intervention. 

<br>
<br>

---

<br>

### __Tutorial__
#### 1. Uploading your data
* Convert the file to the right format
  + accepted formats include '.csv', '.txt', '.xls', '.xlsx'
* Make sure your file contains the following types of columns:
  + Identifiers (gene names/ UniProt ID/ ENSEMBL ID)
  + P-values
  + Fold change values (log2)
* Click __"Download demo data"__ for clarity
* Once the above is checked, press __"Browse"__ to upload your data

#### 2. Preparing for analysis
* Select the identifier column from the drop-down box
* Select the p-value column from the drop-down box
* Select the fold change column from the drop-down box
* Choose the type of species of which your data is sourced from
* Choose a p-value cut off as a determinant for significance
* Choose a (log2) fold change cutoff as a determinant for significance

#### 3. Start the analysis
* Press __"Start"__ to initiate analysis

#### 4. View analysis
* Press __"Analysis"__ on the sidebar to view the results and annotated datasets

<br>
<br>

---

<br>

### __Resources__
TeaProt relies on information from several databases to allow the annotation of user-uploaded data.
* BRENDA
* DisGeNet
* IMPC
* DGIdb
* Human cell atlass
* CHEA3
* MSigDB
* CoffeeProt

