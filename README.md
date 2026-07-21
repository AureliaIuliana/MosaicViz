# MosaicViz

## Introduction 
MosaicViz is a python script implemented to rapidly detect **potential mosaicism in trios**, distinguishing inherited heterozygosity from de novo variants, based on read depth count and percentage of each nucleotide at a specific variant position.

We utilized [bam-readcount](https://github.com/genome/bam-readcount?tab=readme-ov-file) v1.0.1 to determine the count of each nucleotide in each variant position.

The implemented script reprocess [bam-readcount](https://github.com/genome/bam-readcount?tab=readme-ov-file) output, retaining only the depth, ref and alt alleles count data. Additionally, it incorporates supplementary information, including pedigree details, unique identifiers, and the gene symbol associated with each variant. This approach allows for faster detection, saving time compared to manual count check in IGV.

The script execution produces an Excel file composed of four sheets, dedicated to the analysis of **SNVs**, **deletions**, **insertions**, and **delins** being explored, respectively. 
## Project organization 
```
MosaicViz/
├── input/
│   ├── allbams/            # BAM and BAI alignment files directory 
│   ├── hg19_simple_no_chr.fasta
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
│   └── MosDetection.py 
├── environment.yml         
├── LICENSE
└── README.md               
```

## Installation
```
conda env create --file environment.yml
conda activate mosaicviz_env
```

## Input
The input data comprises **BAM** files of the trios being studied, a **CSV** file that details the trio identifiers and pedigree, along with the coordinates of the variants. 

This is the ordered content of the CSV file: 

| Column | Description |
| ------| -----------|
| sampleID   | unique ID associated to each member of the family trio |
| projectID | unique ID which associated each sample to a specific trio research project |
| sampleType    | "I" for proband, "F" for father, "M" for mother |
| geneSymbol | symbol of the gene where the variant is detected |
| coordinates | genomic coordinates pinpointing the location of the variant along with the reference and alternative alleles |
| BAMname| BAM name of the trio member |

CSV example for a SNV, a deletion, a insertion, and a delins respectively: 
```
sampleID,projectID,sampleType,geneSymbol,coordinates,BAMname
2076R,p649,I,DYNC1H1,14:102461022:T/C,OPTI111
2077R,p649,F,DYNC1H1,14:102461022:T/C,OPTI112
2078R,p649,M,DYNC1H1,14:102461022:T/C,OPTI113
19N3158,p705,I,DYNC1H1,14:102452883:TTTACCCGT/-,OPTI135
19N3159,p705,F,DYNC1H1,14:102452883:TTTACCCGT/-,OPTI136
19N3160,p705,M,DYNC1H1,14:102452883:TTTACCCGT/-,OPTI137
19N2563,P133,I,TSC1,9:135787826:-/A,OPTIIH
19N2564,P133,F,TSC1,9:135787826:-/A,OPTIFH
19N2565,P133,M,TSC1,9:135787826:-/A,OPTIMH
20N0349,p714,I,COL4A1,13:110826811:CCTG/GCAC,OPTI141
20N0350,p714,F,COL4A1,13:110826811:CCTG/GCAC,OPTI142
20N0351,p714,M,COL4A1,13:110826811:CCTG/GCAC,OPTI143
```

To streamline access and analysis, BAM files are stored into a single directory. 

## Analysis 

#### Script execution: 

```
python3 MosaicViz.py bam-readcount referenceSequence MAPQvalue CSV BAMsDirectory baseAnalysisExtensionValue lowerThresholdValue upperThresholdValue
```                                         
#### Parameters:
- **`bam-readcount`**: path to bam-readcount executable file 
- **`referenceSequence`**: path to reference sequence in fasta format 
- **`MAPQvalue`**: minimum mapping quality of reads used for counting
- **`CSV`**: CSV file with informations about trio identifiers and pedigree along with variant coordinates and bam files name
- **`BAMsDirectory`**: path to the directory containing all the BAMs to be investigated
- **`baseAnalysisExtensionValue`**: variant upstream and downstream analysis extension value (for INDEL analysis)
- **`lowerThresholdValue`**: lower threshold value to highlight proband/parent alternative allele percentage associated to potential mosaicism
- **`upperThresholdValue`**: upper threshold value to highlight proband/parent alternative allele percentage associated to potential mosaicism

Each BAM file is analyzed with _bam-readcount_ providing as parameters: 
```
bam-readcount -f referenceSequence -q MAPQvalue BAMfile genomicCoordinate
```
- **`referenceSequence`**: path to reference sequence in fasta format 
- **`MAPQvalue`**: minimum mapping quality of reads used for counting
- **`BAMfile`**: the BAM file name associated to each individual
- **`genomicCoordinate`** the genomic coordinate of the variant
  

## Output

The SNV/DELINS-specific output includes the count and percentage abundance of each base type (A, C, G, T, N) at the variant position.
In the specific case of DELINS, we display only the aboundance related to the given position, since we assume is sufficient for mosaicism detection. 
In addition, the background noise is calculated to discriminate between false positive and potential cases of mosaicism. We defined the background noise as the maximum percentage among bases that were neither reference nor alternative within all trio's members.

In the INDEL-specific output are present 2 additional columns: the first displays the count of the identified target INDEL at the given position, while the second its aboundance percentage. The deletion analysis sheet displays the counts in a column marked “DEL:-deletion”. Similarly, for insertions, the data appears under “INS:+insertion”.

The four different Excel sheets are generated only when all variants type are available in the CSV file. 
At the following link you can find an example: [MosaicVizOutput](https://github.com/AureliaIuliana/MosDetection/blob/main/outputExample.xlsx) 

The parameters `lowerThresholdValue` and `upperThresholdValue` highlight in red all variant percentages that fall between the lower and upper thresholds, establishing the bounds for mosaicism. 
Reference allele percentages (greater than the calculated trio background noise) and variant percentages outside the above thresholds are highlighted in yellow.
In the INDEL analysis sheet, the position of the INDEL variant is highlithed in order to differentiate between upstream and downstream positions.

#### Example
We set `lowerThresholdValue=1`,`upperThresholdValue=40` and `baseAnalysisExtensionValue=1`. 

##### _SNV_ _analysis_

In the following example it is possible to detect 3 SNV mosaicisms: two proband mosaicisms (projectID p703 and Pr1) and a father mosiacism (projectID p736).
The analysis of project p736 suggests that G/A is a inherited heterozygosity, distinguishing it from _de_ _novo_ variant.
<img width="1251" alt="image" src="https://github.com/user-attachments/assets/1664cfae-49bd-4077-b37f-3edd58f24a01">

##### _INDEL_ _analysis_

`baseAnalysisExtensionValue` parameter permits to inpect variant flanking regions. As anticipated, the position of the INDEL variant is highlighted to distinguished it from upstream and downstream positions, facilitating data visualization. The percentage of A-base-insertion associated to the mother is not yellow-highlithed since is lesser than background noise, while the percentage of A-base-insertion of the proband is highlighted and indicates that 48,6 % of the total reads count contains the analyzed insertion. 
![image](https://github.com/user-attachments/assets/94f53af2-3abe-454a-9a02-f8a150006d26)

In this deletion example is possibile to visualize that 57% of the reads are reference at the given position, while 42,7% are alternative. 26061 reads cointains the analyzed deletion. 
Since is a 9-base-deletion, and `baseAnalysisExtensionValue=1` we are able to  visualize 1-base upstream and 9-base downstream the given variant position.  
![image](https://github.com/user-attachments/assets/a6e07fb1-b797-44e8-9ea9-f0e88b03cab9)





