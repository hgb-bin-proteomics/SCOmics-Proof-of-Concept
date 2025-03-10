# SCOmics-Proof-of-Concept

A proof-of-concept implementation for automated analysis of (single cell) multi-omic proteomics
and transcriptomics data based on the publication by [Fulcher et al.](https://www.nature.com/articles/s41467-024-54099-z).

Implemented in [Shiny](https://shiny.posit.co/) using [R 4.4.2](https://www.r-project.org/).

Example data can be found here: [Cajun-data/nanoSPLITS](https://github.com/Cajun-data/nanoSPLITS/tree/main/Pooled_C10Cells)

## Running the App

- Install the [R programming language](https://www.r-project.org/).
- Install the [required packages](https://github.com/hgb-bin-proteomics/SCOmics-Proof-of-Concept/blob/2272596bf5ad40033a048189261c7729eafa5aab/scomics/app.R#L1-L15).
- In an R shell run:
  ```R
  shiny::runApp('scomics')
  ```

## Screenshots

### Proteomics

- **Screenshot 1**
  ![Screenshot1](screenshots/1.png)
- **Screenshot 2**
  ![Screenshot2](screenshots/2.png)

### Transcriptomics

- **Screenshot 3**
  ![Screenshot3](screenshots/3.png)

### Clustering

- **Screenshot 4**
  ![Screenshot4](screenshots/4.png)

## Contact

- [micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)
