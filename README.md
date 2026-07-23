# MosaicViz

## Introduction 
MosaicViz is a python script that quantifies nucleotide read counts at target genomic positions using bam-readcount. It calculates nucleotides percentages across family trios and
generates an Excel output file with color-coded conditional formatting for rapid visual screening of **potential mosaic variants**.

It leverages [bam-readcount](https://github.com/genome/bam-readcount?tab=readme-ov-file) v1.0.1 to determine the count of each nucleotide in each variant position.

The implemented script reprocess [bam-readcount](https://github.com/genome/bam-readcount?tab=readme-ov-file) output, retaining only the depth, ref and alt alleles count data. Additionally, it incorporates supplementary information, including pedigree details, unique identifiers, and the gene symbol associated with each variant. 

The Excel file is composed of four sheets dedicated to the inspection of **SNVs**, **deletions**, **insertions**, and **delins**.

## Project organization 
```
MosaicViz/
├── input/
│   ├── testbams/            # BAM and BAI alignment files directory 
│   ├── hg19_simple_no_chr.fasta
│   ├── hg19_simple_no_chr.fasta.fai
│   └── info.csv            # Sample metadata and target regions
├── output/                 
│   └── results.xlsx        # Generated Excel results and log files
├── src/
│   ├── modules/
│   │   ├── __init__.py
│   │   ├── bamreadcount_exe.py
│   │   ├── indels_analysis.py
│   │   ├── preprocessing.py
│   │   └── snv_analysis.py
│   ├── config.py
│   └── MosaicViz.py 
├── environment.yml         
├── LICENSE
└── README.md               
```

## Conda env 
```
conda env create --file environment.yml
conda activate mosaicviz_env
```

## Input 
- Directory with **BAM** files of the trios being studied
- **CSV** file that details the trio identifiers and pedigree, along with the coordinates of the variants 

CSV example for SNV, deletion, insertion, and a delins respectively: 
```
sampleID,projectID,sampleType,geneSymbol,coordinates,BAMname
2076R,p649,I,DYNC1H1,14:102461022:T/C,OPTI111_recal
2077R,p649,F,DYNC1H1,14:102461022:T/C,OPTI112_recal
2078R,p649,M,DYNC1H1,14:102461022:T/C,OPTI113_recal
19N3158,p705,I,DYNC1H1,14:102452883:TTTACCCGT/-,OPTI135_recal
19N3159,p705,F,DYNC1H1,14:102452883:TTTACCCGT/-,OPTI136_recal
19N3160,p705,M,DYNC1H1,14:102452883:TTTACCCGT/-,OPTI137_recal
19N2563,P133,I,TSC1,9:135787826:-/A,OPTIIH_recal
19N2564,P133,F,TSC1,9:135787826:-/A,OPTIFH_recal
19N2565,P133,M,TSC1,9:135787826:-/A,OPTIMH_recal
20N0349,p714,I,COL4A1,13:110826811:CCTG/GCAC,OPTI141_recal
20N0350,p714,F,COL4A1,13:110826811:CCTG/GCAC,OPTI142_recal
20N0351,p714,M,COL4A1,13:110826811:CCTG/GCAC,OPTI143_recal
```
Column description
| Column | Description |
| ------| -----------|
| sampleID   | unique ID associated to each member of the family trio |
| projectID | unique ID which associated each sample to a specific trio research project |
| sampleType    | "I" for proband, "F" for father, "M" for mother |
| geneSymbol | symbol of the gene where the variant is detected |
| coordinates | genomic coordinates pinpointing the location of the variant along with the reference and alternative alleles |
| BAMname| BAM name of the trio member |

## Script execution: 

```
python3 MosaicViz.py -r PATH_REFERENCE -c PATH_CSV -d PATH_DIR_BAMS -o PATH_OUT_DIR
```                                         

```
python3 MosaicViz.py --help
Arguments:
  -h, --help              
  -v, --version            
  -b , --bam-readcount     path to bam-readcount executable file (default: bam-readcount)
  -r , --reference         path to reference sequence in fasta format
  -c , --csv               csv file with informations about trio identifiers, pedigree along with variant coordinates, and bam files name (default: None)
  -d , --bams-directory    path to the directory containing all the BAMs to be investigated (default: None)
  -o , --output-dir        path to output dir (default: None)
  -m , --mapq              minimum mapping quality of reads used for counting (default: 0)
  -e , --extension         variant upstream and downstream analysis extension value for INDEL inspection (default: 1)
  -l , --lower-threshold   lower threshold value to highlight proband/parent alternative allele percentage associated to potential mosaicism (default: 1)
  -u , --upper-threshold   upper threshold value to highlight proband/parent alternative allele percentage associated to potential mosaicism (default: 40)
```   

## Output

The SNV/DELINS-specific output includes the count and percentage abundance of each base type (A, C, G, T, N) at the variant position.
In the specific case of DELINS, we display only the aboundance related to the given position, since we assume is sufficient for mosaicism detection. 
In addition, the background noise is calculated to discriminate between false positive and potential cases of mosaicism. We defined the background noise as the maximum percentage among bases that were neither reference nor alternative within all trio's members.

In the INDEL-specific output are present 2 additional columns: the first displays the count of the identified target INDEL at the given position, while the second its aboundance percentage. The deletion analysis sheet displays the counts in a column marked “DEL:-deletion”. Similarly, for insertions, the data appears under “INS:+insertion”.

The four different Excel sheets are generated only when all variants type are available in the CSV file. 
At the following link you can find an example: [MosaicVizOutput](https://github.com/AureliaIuliana/MosDetection/blob/main/outputExample.xlsx) 

The parameters `--lower-threshold` and `--upper-threshold` highlight in red all variant percentages that fall between the lower and upper thresholds, establishing the bounds for mosaicism. 
Reference allele percentages (greater than the calculated trio background noise) and variant percentages outside the above thresholds are highlighted in yellow.
In the INDEL analysis sheet, the starting position of the INDEL variant is highlithed in order to differentiate between upstream and downstream positions.

##### _SNV_ _analysis_

In the following example it is possible to detect 3 SNV mosaicisms: two proband mosaicisms (projectID p703 and Pr1) and a father mosiacism (projectID p736).
The analysis of project p736 suggests that G/A is a inherited heterozygosity, distinguishing it from _de_ _novo_ variant.
<img width="1251" alt="image" src="https://github.com/user-attachments/assets/1664cfae-49bd-4077-b37f-3edd58f24a01">

##### _INDEL_ _analysis_

`--extension` parameter permits to inpect variant flanking regions. As anticipated, the starting position of the INDEL variant is highlighted to distinguished it from upstream and downstream positions, facilitating data visualization. 
In this insertion example the percentage of A-base-insertion associated to the mother is not yellow-highlithed since is lesser than background noise, while the percentage of A-base-insertion of the proband is highlighted and indicates that 48,6 % of the total reads count contains the analyzed insertion. 
![image](https://github.com/user-attachments/assets/94f53af2-3abe-454a-9a02-f8a150006d26)

In this deletion example is possibile to visualize that 57% of the reads are reference at the given position, while 42,7% are alternative. 26061 reads cointains the analyzed deletion. 
Since is a 9-base-deletion, and `baseAnalysisExtensionValue=1` we are able to  visualize 1-base upstream and 9-base downstream the given variant position.  
![image](https://github.com/user-attachments/assets/a6e07fb1-b797-44e8-9ea9-f0e88b03cab9)





