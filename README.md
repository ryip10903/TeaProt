### __Introduction__

TeaProt is an online Shiny tool that integrates upstream transcription factor enrichment analysis with downstream pathway analysis through an easy-to-use interactive interface. TeaProt maps user’s omics data with online databases to provide a collection of annotations on drug-gene interactions, subcellular localizations, phenotypic functions, gene-disease associations and enzyme-gene interactions, usefull for further analyses. Users can combine TeaProt and urPTMdb for a novel and easy-to-use online proteomics/transcriptomics analysis pipeline featuring novel underrepresented genesets to allow the discovery of downstream cellular processes, upstream transcriptional regulation and classes of PTMs potentially regulated by a users’ intervention. 

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

---

<br>

### __Browser compatibility__

<style>
.basic-styling td,
.basic-styling th {
  border: 1px solid #555;
  padding: 1rem;
}
</style>

<div class="ox-hugo-table basic-styling">
<div></div>
<div class="table-caption">
  <span class="table-number"></span>
</div>

|OS     |version           |   Chrome   |  Firefox |Microsoft Edge|  Safari  |
|-------|:----------------:|:----------:|:--------:|:------------:|:--------:|
|Linux  |Ubuntu 20.04.1 LTS|87.0.4280.88|78.0.1    |n/a           |n/a       |
|MacOS  |10.13.6           |87.0.4280.67|83.0      |n/a           |13.1.2    |
|Windows|10                |87.0.4280.88|83.0      |87.0.664.55   |n/a       |

</div>

---

<br>

### __Contact__
For technical support, please email support@coffeeprot.com. To contact the Parker lab, please contact ben.parker@unimelb.edu.au.

---

<br>

### __Citation__
<p align="justify">Citation details will be added soon.</p>

---

<br>

### __Acknowledgements__
<p align="justify"> This research was supported by use of the Nectar Research Cloud and by the University of Melbourne Research Platform Services. The Nectar Research Cloud is a collaborative Australian research platform supported by the National Collaborative Research Infrastructure Strategy. This work was funded by an Australian National Health and Medical Research Council Ideas Grant (APP1184363) and The University of Melbourne Driving Research Momentum program. </p>

Human Protein Atlas subcellular localization data was obtained from http://www.proteinatlas.org and has previously been described in  <a href="http://dx.doi.org/10.1126%2Fscience.aal3321">Thul PJ et al., A subcellular map of the human proteome. Science. (2017)</a>. Drug-gene interaction data was obtained from DGIdb (https://www.dgidb.org/downloads). Genotype-phenotype associations were downloaded from the International Mouse Phenotyping Consortium (IMPC, www.mousephenotype.org). Enzymatic annotations were retrieved from BRENDA (https://www.brenda-enzymes.org/). Disease-Gene annotations were retrieved from DisGeNet (https://www.disgenet.org/). Transcription factor data was downloaded from CHEA3 (maayanlab.cloud/chea3/). The DNA vector image used in the TeaProt banner on the Welcome page was obtained from Vecteezy (<a href="https://www.vecteezy.com/vector-art/1270772-human-dna-design">Human dna design  Vectors by Vecteezy</a>).

<br>
