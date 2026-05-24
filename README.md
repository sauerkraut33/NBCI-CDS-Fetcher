# NCBI-CDS-FETCHER

A semi-automated Python program for retrieving coding sequences (CDS) from NCBI for a target gene across a list of species.

The program is designed to reduce repetitive manual searching in NCBI while still allowing manual review of transcript variants, CDS length, protein/product annotation, and alignment quality.

## What the Program Currently Does

The current version:

1. Reads a species list from an Excel file.
2. Prompts the user to enter one or more target gene, such as `CALM1` or `CALM1,CNGA1,AIPL1`.
3. Searches NCBI Nucleotide for matching RefSeq mRNA records for each species.
4. Keeps RefSeq mRNA records with:
    - `NM_` accessions: RefSeq curated mRNA
    - `XM_` accessions: RefSeq predicted mRNA
    Prioritizes NM_ records when both NM_ and XM_ records are found for the same gene/species search.
5. Fetches GenBank records for the candidate accessions.
6. Extracts CDS features matching the target gene.
7. Calculates CDS length from the extracted CDS sequence.
8. Creates a gene-specific output folder.
9. Writes all retrieved CDS sequences into a combined `.txt` file.
10. Generates a tracking sheet with accession numbers, sequence types, CDS lengths, notes, and NCBI record titles.
    Records rejected candidates in the tracking sheet if an NCBI search result is retrieved but the CDS gene annotation does not exactly match the target gene.

## Input File

The program expects a species list spreadsheet in the `input/` directory:
```
input/Additional_sampling_species.xlsx
```
The spreadsheet should contain these columns:
```
Taxonomic group
Species (Common Name)
Species (Scientific Name)
```
Example:
```
Monotremes | Platypus | Ornithorhynchus anatinus
Rodentia | Mouse | Mus musculus
Primate | Human | Homo sapiens
```
## How to Run

Install required packages:
```
pip install biopython pandas openpyxl
```
Run the script:
```
python fetch_cds.py
```
When prompted, enter the target gene name:
```
Enter target gene(s), separated by commas, e.g. CALM1,CNGA1,AIPL1: CALM1,CNGA1,AIPL1
```

## Current Output

For a target gene such as CALM1, the program creates a folder:
```
output/CALM1/
output/CNGA1/
output/AIPL1/
```
Inside a folder, e.g folder CALM1 it creates:
```
CALM1.txt
CALM1_tracking_sheet.csv
```
The .txt file contains the retrieved CDS sequences in FASTA-style format, with the species common name added at the beginning of each sequence header before the NCBI-style CDS annotation.
It is saved as .txt first so the user can manually review and remove unwanted sequences before converting it to .fas or .fasta for alignment in MEGA.

The tracking sheet contains:
```
Taxonomic group
Species (Common Name)
Species (Scientific Name)
Accession Numbers
Sequence Type
CDS Length (bp)
Notes
Title
```

## Current Selection Logic

At this stage, the program retrieves matching NM_ and XM_ RefSeq mRNA candidates that contain a CDS feature matching the target gene.

If both NM_ and XM_ records are found for the same gene/species search, the current version prioritizes NM_ curated RefSeq records for the main output. If no NM_ records are found, XM_ predicted RefSeq records are used.

If an NCBI search result is retrieved but its CDS gene annotation does not exactly match the target gene, the sequence is not added to the combined .txt file. The tracking sheet records this as a rejected candidate and lists the CDS gene annotation that was found.

The program does not yet automatically decide which transcript variant should be kept as the final selected sequence.

This is intentional. The combined sequence file is meant to be manually reviewed before alignment.

Example Output Structure
```
output/
  CALM1/
    CALM1.txt
    CALM1_tracking_sheet.csv
```
Later, after manual review, the .txt file can be saved or renamed as:
```
CALM1.fas
```
and opened in MEGA for MUSCLE codon alignment.

## Planned Features
### 1. Manual deletion and transcript variant management

I plan to add a function for manually removing unwanted transcript variants from the main sequence file.

The goal is:

Keep one selected representative CDS per species in the main gene .txt file.
Move removed transcript variants into a separate file or folder.
Mark selected and removed sequences in the tracking sheet.

Possible future output structure:
```
output/
  CALM1/
    CALM1.txt
    CALM1_tracking_sheet.csv
    CALM1_extra_transcript_variants.txt
```
This would allow extra transcript variants to be preserved in case they are needed later, while keeping the main sequence file cleaner for alignment.

### 2. Protein/product annotation checking

I plan to add an additional user input for the expected protein or product keyword related to the target gene.

For example, for CALM1, the expected protein/product keyword could be:
```
calmodulin
```
The program could then check whether the CDS product annotation contains the expected keyword, such as:
```
calmodulin-1
calmodulin-1 isoform 2
neo-calmodulin
```
This would help flag records where the gene name matches but the protein/product annotation may need manual review.

## Notes

This program is not intended to fully replace manual sequence review.

The main goal is to automate repetitive NCBI searching and CDS extraction while keeping the final transcript selection, sequence cleanup, and alignment quality control transparent and manually reviewable.
