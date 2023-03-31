# RxGRID for CRISPR-guided drug repurposing

RxGRID (Rapid proXimity Guidance for Repurposing Investigational Drugs) is a method for prioritizing candidates for drug repurposing based on network proximity to genetic hits from high-throughput screens. For additional documentation, please refer to our preprint, available [here](https://www.medrxiv.org/content/10.1101/2022.07.02.22277181v1.full). 

![Alt text](graphical_overview.png?raw=true "RxGRID Overview")

## Dependencies
RxGRID is built in Python 3 and requires NetworkX, NumPy, pandas, Matplotlib, and seaborn, in addition to the standard Python library. It was tested using Python 3.7, NetworkX 2.2, NumPy 1.21.6, pandas 1.2.3, Matplotlib 3.5.3, and seaborn 0.12.2. We provide an optional conda environment containing all necessary dependencies, which can be activated by running the following:
```
$ conda env create -f RxGRID_environment.yml
$ conda activate RxGRID
```

## Installation

Clone the repository as follows:

```
$ git clone https://github.com/henrycousins/RxGRID
$ cd RxGRID
```

## Simple usage

RxGRID can be run through the command line as follows:

```
python RxGRID.py --rnk_folder RNK_FOLDER --output_folder OUTPUT_FOLDER
```
where ```RNK_FOLDER``` is the name of the folder containing the collection of screens (as rnk files) to be analyzed and ```OUTPUT_FOLDER``` is the name of the existing folder to which the results should be saved. For instance, the following code performs RxGRID on the set of screens contained in ```/screens```, saving the results to ```/results```:

<!-- the path to the .rnk file containing the gene list, ```GMT_FILE``` is the path to the .gmt file containing the gene sets, and ```OUTPUT_FOLDER``` is name of the directory where results should be stored. (Please note that ```OUTPUT_FOLDER``` is created automatically.) For instance, the following code performs GSPA on the GSE4183 gene list using KEGG pathway gene sets, saving the results to ```/outputs```: -->

```
python RxGRID.py --rnk_folder screens --output_folder results
```

## Advanced usage

Additional parameters can be defined by the user.
```
usage: RxGRID.py [-h] [--gene_hit_threshold GENE_HIT_THRESHOLD]
                 [--drug_hit_threshold DRUG_HIT_THRESHOLD]
                 [--centrality_metric CENTRALITY_METRIC]
                 [--mutual_hit_freq MUTUAL_HIT_FREQ]
                 [--drug_exclusions_path DRUG_EXCLUSIONS_PATH]
                 [--rnk_folder RNK_FOLDER] [--output_folder OUTPUT_FOLDER]
                 [--save_plots SAVE_PLOTS]

optional arguments:
  -h, --help            show this help message and exit
  --gene_hit_threshold GENE_HIT_THRESHOLD
                        Significance cutoff for genes (default=0.05)
  --drug_hit_threshold DRUG_HIT_THRESHOLD
                        Significance cutoff for drugs (default=0.01)
  --centrality_metric CENTRALITY_METRIC
                        Primary centrality metric for ranking drugs (degree
                        ["ndc"] or betweenness ["btw"]; default="ndc"
  --mutual_hit_freq MUTUAL_HIT_FREQ
                        Fraction of screens in which drug must be hit to meet
                        overall significance (inputs >=1 are interpreted as
                        counts, rather than fractions; default = 1/3)
  --drug_exclusions_path DRUG_EXCLUSIONS_PATH
                        CSV file containing user-specified list of drugs to
                        exclude from downstream analysis
                        (default="data/drug_exclusions.csv")
  --rnk_folder RNK_FOLDER
                        Name of folder containing the screens, as rnk files
                        (default="screens")
  --output_folder OUTPUT_FOLDER
                        Name of folder to which results should be saved
                        (default="results")
  --save_plots SAVE_PLOTS
                        Whether to save descriptive plots (default=False)
```


## Detailed file descriptions

File | Type | Extension | Columns | Comments
--- | --- | --- | --- | ---
```RNK_FILE``` | Input | .rnk | [gene ID, gene score] | [Description](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29)
```<jobID>_summary``` | Output | .csv | [screenwise hit count, hit status {per-screen}, percentile rank {mean, std, per-screen}, centrality score {mean, std, per-screen}] | Created automatically
<!-- 

```
.
└── gspa
    ├── LICENSE           # LICENSE
    ├── README.md         # README
    ├── embeddings        # contains gene embeddings and IDs as pickled arrays
    ├── environments      # optional dependency environments
    ├── functions.py      # helper functions
    ├── gene_sets         # contains example gene sets
    ├── gspa.py           # main
    ├── outputs           # example directory for results
    ├── overview.png      # overview graphic
    ├── rnk_files         # contains example rnk files
    └── data              # supplementary data
        ├── ppi           # related to ppi network and embeddings
        ├── sarscov2      # related to sars-cov-2 screens
        └── set_matching  # related to comparing similar gene sets
``` -->

## License

This software is available under an MIT License.

## Sustainability

We anticipate releasing core updates on a 6-month basis for a minimum of 2 years and will respond to any issues regarding the codebase, which is open-source, as they arise. Supporting funding for RxGRID is provided by the National Institutes of Health (GM007365, GM102365, LM013337, and HG011316), by the Donald and Delia Baxter Foundation, by the National Science Foundation (1953686 and 1953415), and by the Knight-Hennessy Scholarships.