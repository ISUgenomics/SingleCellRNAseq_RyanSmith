# SCRNAseq runthrough


### Data
```
/work/gif/remkv6/Smith_Ryan/01_Monocle

#transposed Severo and Smith datasets in excel to obtain these files

SeveroHemocyteRNA.txt
SmithHemocyteRNA.txt
SmithHemocyteRNATransposePower2.txt #transposed, and inversed the log2 normalization
SeveroHemocyteRNAOrig.txt


#concatenate datasets, add missing genes to each dataset to create expression matrices
awk 'NR>3' SeveroHemocyteRNAOrig.txt |cut -f -26 |cat - <(awk 'NR>5{print $1,"0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"} ' SmithHemocyteRNATransposePower2.txt |tr " " "\t" ) <(awk  '{print $1,"0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"}' AllGenesSmithSevero.txt |tr " " "\t")|sort -u -k 1,1 >SeveroAllGenes.txt


awk 'NR>4' SmithHemocyteRNATransposePower2.txt |cat - <(awk 'NR>3{print $1,"0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"} ' SeveroHemocyteRNAOrig.txt |tr " " "\t" ) <(awk  '{print $1,"0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"}' AllGenesSmithSevero.txt |tr " " "\t")|sort -u -k 1,1 >SmithAllGenes.txt



paste SeveroAllGenes.txt SmithAllGenes.txt |cut -f 1-26,28- >tailer.txt
paste  <(awk 'NR==4' SeveroHemocyteRNAOrig.txt|cut -f -26) <(awk 'NR==1' SmithHemocyteRNATransposePower2.txt )   >header.txt
cat header.txt tailer.txt >AllGenesSmithSevero.txt

cat <(awk 'NR==1' SmithHemocyteRNATransposePower2.txt ) SmithAllGenes.txt >SmithExpressionMatrix
cat <(awk 'NR==4' SeveroHemocyteRNAOrig.txt|cut -f -26) SeveroAllGenes.txt >SeveroExpressionMatrix

```

### Run Monocle for Smith data only, expected result 8 clusters?
```

#Create gene metadata from gff and missing genes added from studies
wget ftp://ftp.ensemblgenomes.org/pub/release-46/metazoa/gff3/anopheles_gambiae/Anopheles_gambiae.AgamP4.46.gff3.gz

less -S Anopheles_gambiae.AgamP4.46.gff3 |awk -F"\t" '$3=="gene"{print $9}' |sed 's/ID=gene://g' |sed 's/;/\t/1' >gene_metadata
cat gene_metadata <(awk  'NR>1{print $1}' AllGenesSmithSevero.txt) |sort -u -k1,1 >AllGeneMetadata

#cell metadata
vi cell_metadata #copied from excel
cat cell_metadata <(awk 'NR==4' SeveroHemocyteRNAOrig.txt|cut -f -26  |tr "\t" "\n" |awk '{print $1,"naive"}' |tr " " "\t") >Cell_metadata_Combine

 less cell_metadata |awk 'NR<241' >Smithcell_metadata
```

### Smith only experiment

```
SmithAllGenes.txt
cell_metadata
AllGeneMetadata

awk 'NR>1' SmithExpressionMatrix >SmithExpressionMatrix.noheaders

# Load the data
library(monocle3)
library(dplyr)
expression_matrix<- as.matrix("SmithExpressionMatrixSparseHeadersIntact",header = T,quote = "",row.names = 1)
cell_metadata<- as.matrix("RedoCellMetadata",header = T,quote = "",row.names = 1)
gene_metadata <- as.matrix("AllGeneMetadata",header = T,quote = "",row.names = 1)

cds <- new_cell_data_set(expression_matrix,
cell_metadata = cell_metadata,
gene_metadata = gene_metadata)
#Error: row.names of cell_metadata must be equal to colnames of expression_data
#In addition: Warning message:
#In storage.mode(from) <- "double" : NAs introduced by coercion

#have tried removing the header and first column to make them equal, changed the first column of SmithExpressionMatrix from "cell" to "gene".  Nothing has gotten past this error

Perhaps it needs to be space delimited?
#downloaded rdf files from their tutorial, but not sure how to open.
#They showt that you can import cell ranger objects, but the format there is not sensical to me and I do not have the reads to repeat.

page(cell_metadata, method = "print")
 #################################################################################
 plate cao_cluster            cao_cell_type
cele-001-001.CATGACTCAA   001          20     Unclassified neurons
cele-001-001.AAGACGGCCA   001           6                 Germline
cele-001-001.GCCAACGCCA   001          13 Intestinal/rectal muscle
cele-001-001.ATAGGAGTAC   001          27        Vulval precursors
cele-001-001.CTCGTCTAGG   001           2             Coelomocytes
cele-001-001.AAGTTGCCAT   001           6                 Germline
cele-001-001.GTTGAAGGAT   001           9                     <NA>
cele-001-001.AACTACGGCT   001          19 Ciliated sensory neurons
cele-001-001.CATTCTAGAA   001          10                Failed QC
cele-001-001.AATACTCTTC   001           7               Seam cells
cele-001-001.GAGGCTTATT   001          19 Ciliated sensory neurons
cele-001-001.ATGCCGGACG   001          25      Non-seam hypodermis
cele-001-001.AGCGGTAACG   001          11     Pharyngeal epithelia
###################################################################################


 page(expression_matrix, method = "print")
###################################################################################
20271 x 42035 sparse Matrix of class "dgCMatrix"

WBGene00000001 .  .  . .  .  . .  . .  . . .  . . . .  . .  .  .  .  .  .  .  .  .  .  .  .  .
WBGene00000002 .  .  . .  .  . .  . .  . . .  . . . .  . .  .  .  .  .  .  .  .  .  .  .  .  .
WBGene00000003 .  .  . .  .  . .  . .  . . .  . . . .  . .  .  .  .  .  .  .  .  .  .  .  2  .
WBGene00000004 .  .  . .  .  . .  . .  . . .  . . . .  . .  .  .  .  .  .  .  .  .  .  .  .  .
WBGene00000005 .  .  . .  .  . .  . .  . . .  . . . .  . .  .  .  .  .  .  .  .  .  .  .  .  .
WBGene00000006 .  .  . .  .  . 1  . .  . . .  . . . .  . .  .  .  .  .  .  .  2  .  .  .  .  .
WBGene00000007 .  .  . .  .  . .  . .  . . .  . . . .  . .  .  .  .  .  .  .  .  .  .  .  .  .
WBGene00000008 .  .  . .  .  . .  . .  . . .  . . . .  . .  .  .  .  .  .  .  .  .  .  .  .  .
WBGene00000009 .  .  . .  .  . .  . .  . . .  . . . .  . .  .  .  .  .  .  .  .  .  .  .  .  .
###################################################################################

 page(gene_annotation, method = "print")
###################################################################################
gene_short_name
WBGene00000001           aap-1
WBGene00000002           aat-1
WBGene00000003           aat-2
WBGene00000004           aat-3
WBGene00000005           aat-4
WBGene00000006           aat-5
WBGene00000007           aat-6
WBGene00000008           aat-7
WBGene00000009           aat-8
WBGene00000010           aat-9
############################################################################

> dim(cell_metadata)
[1] 42035     4
> dim(expression_matrix)
[1] 20271 42035
> dim(gene_annotation)
[1] 20271     1


#So I need to format the


pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSSM_gene_annotation)
HSMM <- new_cell_data_set(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
cell_metadata = cell_metadata,
gene_metadata = gene_annotation)
```

### Smith +Severo experiment
```
AllGenesSmithSevero.txt
Cell_metadata_Combine
AllGeneMetadata


```
