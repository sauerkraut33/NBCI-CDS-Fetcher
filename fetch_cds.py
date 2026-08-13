from Bio import Entrez, SeqIO
import time
from pathlib import Path
import pandas as pd


Entrez.email = "runchen.yang@mail.utoronto.ca"


def read_species_list(file_path: Path) -> list[dict]:
    """Read the species spreadsheet and return species records."""
    df = pd.read_excel(file_path)

    print("Columns found:")
    print(df.columns.tolist())

    species_records = []

    for _, row in df.iterrows():
        species_records.append({
            "taxonomic_group": row["Taxonomic group"],
            "common_name": row["Species (Common Name)"],
            "scientific_name": row["Species (Scientific Name)"],
        })

    return species_records

def build_candidates_from_nucleotide_ids(ids: list[str], gene: str, species: str) -> list[dict]:
    """Build RefSeq mRNA candidate records from NCBI Nucleotide IDs."""
    if not ids:
        return []

    with Entrez.esummary(db="nucleotide", id=",".join(ids)) as handle:
        summaries = Entrez.read(handle)

    candidates = []

    for item in summaries:
        accession = str(item.get("AccessionVersion", ""))
        title = str(item.get("Title", ""))
        length = int(item.get("Length", 0))

        # Keep only RefSeq mRNA records for now
        if not accession.startswith(("NM_", "XM_")):
            continue

        if accession.startswith("NM_"):
            record_type = "RefSeq curated mRNA"
        else:
            record_type = "RefSeq predicted mRNA"

        candidates.append({
            "gene": gene,
            "species": species,
            "accession": accession,
            "title": title,
            "mRNA_length": length,
            "record_type": record_type,
            "notes": "",
        })

    time.sleep(0.4)
    return candidates


def search_ncbi_candidates(gene: str, species: str, max_results: int = 20) -> list[dict]:
    """Search NCBI Nucleotide for mRNA records matching one gene and one species."""
    query = (
        f'{gene}[Gene Name] AND "{species}"[Organism] '
        f'AND biomol_mrna[PROP]'
    )

    print(f"Searching: {query}")

    with Entrez.esearch(db="nucleotide", term=query, retmax=max_results) as handle:
        search_record = Entrez.read(handle)

    return build_candidates_from_nucleotide_ids(
        ids=[str(nucleotide_id) for nucleotide_id in search_record["IdList"]],
        gene=gene,
        species=species,
    )


def find_ncbi_gene_matches(gene: str, species: str, max_results: int = 20) -> list[dict]:
    """
    Find species-specific NCBI Gene records whose official symbol or alias
    exactly matches the requested gene.
    """
    query = f'{gene}[All Fields] AND "{species}"[Organism]'
    print(f"Searching NCBI Gene for an alternative gene name: {query}")

    with Entrez.esearch(db="gene", term=query, retmax=max_results) as handle:
        search_record = Entrez.read(handle)

    gene_ids = [str(gene_id) for gene_id in search_record["IdList"]]
    if not gene_ids:
        return []

    with Entrez.esummary(db="gene", id=",".join(gene_ids)) as handle:
        summary_record = Entrez.read(handle)

    summaries = summary_record["DocumentSummarySet"]["DocumentSummary"]
    requested_upper = gene.upper()
    exact_matches = []

    for item in summaries:
        official_symbol = str(item.get("Name", "")).strip()
        aliases = [
            alias.strip()
            for alias in str(item.get("OtherAliases", "")).split(",")
            if alias.strip()
        ]
        organism = str(item.get("Organism", {}).get("ScientificName", "")).strip()

        searchable_names = [official_symbol, *aliases]
        name_matches = any(
            name.upper() == requested_upper
            for name in searchable_names
        )

        if organism.lower() != species.lower() or not name_matches:
            continue

        exact_matches.append({
            "gene_id": str(item.attributes.get("uid", "")),
            "official_symbol": official_symbol,
            "aliases": aliases,
            "description": str(item.get("Description", "")).strip(),
        })

    time.sleep(0.4)
    return exact_matches


def get_gene_linked_candidates(
        gene_id: str,
        official_symbol: str,
        species: str,
) -> list[dict]:
    """Retrieve RefSeq RNA records linked directly to one NCBI GeneID."""
    print(f"Checking RefSeq RNA records linked to NCBI GeneID {gene_id}")

    with Entrez.elink(
            dbfrom="gene",
            db="nucleotide",
            id=gene_id,
            linkname="gene_nuccore_refseqrna",
    ) as handle:
        link_records = Entrez.read(handle)

    nucleotide_ids = []
    for link_set in link_records:
        for link_database in link_set.get("LinkSetDb", []):
            for link in link_database.get("Link", []):
                nucleotide_ids.append(str(link["Id"]))

    time.sleep(0.4)
    return build_candidates_from_nucleotide_ids(
        ids=nucleotide_ids,
        gene=official_symbol,
        species=species,
    )


def get_alternative_gene_candidates(
        requested_gene: str,
        species: str,
        max_results: int = 20,
) -> tuple[list[dict], str, str]:
    """
    Resolve a requested gene through NCBI Gene, then combine candidates found
    by the official-symbol search and the GeneID link.
    """
    gene_matches = find_ncbi_gene_matches(
        gene=requested_gene,
        species=species,
        max_results=max_results,
    )

    alternative_matches = [
        match for match in gene_matches
        if match["official_symbol"].upper() != requested_gene.upper()
    ]

    if not alternative_matches:
        return [], "", ""

    alternative_symbols = {
        match["official_symbol"].upper()
        for match in alternative_matches
    }
    if len(alternative_symbols) > 1:
        symbols = ", ".join(
            sorted(match["official_symbol"] for match in alternative_matches)
        )
        return [], "", (
            "multiple possible alternative gene names were found "
            f"({symbols}); manual review required"
        )

    gene_match = alternative_matches[0]
    official_symbol = gene_match["official_symbol"]

    symbol_candidates = search_ncbi_candidates(
        gene=official_symbol,
        species=species,
        max_results=max_results,
    )
    linked_candidates = get_gene_linked_candidates(
        gene_id=gene_match["gene_id"],
        official_symbol=official_symbol,
        species=species,
    )

    candidates_by_accession = {}
    for candidate in [*symbol_candidates, *linked_candidates]:
        candidates_by_accession[candidate["accession"]] = candidate

    note = (
        f"found using the alternative gene name {official_symbol} "
        f"provided by NCBI Gene"
    )
    return list(candidates_by_accession.values()), official_symbol, note

def prioritize_nm_over_xm(candidates: list[dict]) -> list[dict]:
    """
    If curated NM_ records are found for a species/gene,
    keep only NM_ candidates for the main output.
    If no NM_ records are found, keep XM_ candidates.
    """
    nm_candidates = [
        candidate for candidate in candidates
        if candidate["accession"].startswith("NM_")
    ]

    if nm_candidates:
        return nm_candidates

    return candidates

def fetch_genbank(accession: str):
    """Fetch one NCBI nucleotide record in GenBank format."""
    print(f"Fetching GenBank record: {accession}")

    with Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text",
    ) as handle:
        record = SeqIO.read(handle, "genbank")

    time.sleep(0.4)
    return record

def extract_cds(record, expected_gene: str) -> tuple[list[dict], list[str]]:
    """
    Extract CDS features from a GenBank record that match the expected gene.
    Also return gene names found in CDS features for clearer tracking notes.
    """
    cds_list = []
    found_cds_genes = []

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene_name = feature.qualifiers.get("gene", [""])[0]
        product = feature.qualifiers.get("product", [""])[0]
        protein_id = feature.qualifiers.get("protein_id", [""])[0]
        db_xrefs = feature.qualifiers.get("db_xref", [])
        is_low_quality_protein = "LOW QUALITY PROTEIN" in product.upper()

        if gene_name and gene_name not in found_cds_genes:
            found_cds_genes.append(gene_name)

        gene_matches = gene_name.upper() == expected_gene.upper()

        cds_seq = feature.extract(record.seq)
        location = format_cds_location(feature)

        if gene_matches:
            cds_list.append({
                "accession": record.id,
                "gene": gene_name,
                "product": product,
                "protein_id": protein_id,
                "db_xrefs": db_xrefs,
                "location": location,
                "cds_length": len(cds_seq),
                "sequence": str(cds_seq),
                "is_low_quality_protein": is_low_quality_protein,
            })

    return cds_list, found_cds_genes

def build_ncbi_style_header(common_name: str, cds: dict, variant_note: str = "") -> str:
    """Build a FASTA header similar to NCBI's downloaded CDS FASTA header."""
    clean_common_name = clean_name(common_name)

    # Add transcript variant to the common name if present.
    # Example: Mouse_Transcript_Variant_2
    if variant_note:
        clean_variant = clean_name(variant_note)
        display_name = f"{clean_common_name}_{clean_variant}"
    else:
        display_name = clean_common_name

    accession = cds["accession"]
    protein_id = cds["protein_id"]
    gene = cds["gene"]
    product = cds["product"]
    location = cds["location"]

    # db_xref part
    db_xref_text = ""
    if cds["db_xrefs"]:
        db_xref_text = f" [db_xref={','.join(cds['db_xrefs'])}]"

    header = (
        f">{display_name} "
        f"lcl|{accession}_cds_{protein_id}_1 "
        f"[gene={gene}]"
        f"{db_xref_text} "
        f"[protein={product}] "
        f"[protein_id={protein_id}] "
        f"[location={location}] "
        f"[gbkey=CDS]"
    )

    return header

def clean_name(name: str) -> str:
    """Make a name safe for FASTA headers."""
    return (
        str(name)
        .strip()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
        .replace(",", "")
    )


def get_variant_note(title: str) -> str:
    """Extract transcript variant information from the NCBI title if present."""
    title_lower = title.lower()

    if "transcript variant" not in title_lower:
        return ""

    # Example:
    # "Mus musculus calmodulin 1 (Calm1), transcript variant 2, mRNA"
    parts = title.split(",")

    for part in parts:
        if "transcript variant" in part.lower():
            return part.strip()

    return "transcript variant mentioned"

def add_note(existing_notes: str, new_note: str) -> str:
    """Add a note using ';' if there is already an existing note."""
    if not new_note:
        return existing_notes

    if not existing_notes:
        return new_note

    return existing_notes + "; " + new_note

def format_cds_location(feature) -> str:
    """Convert a CDS feature location to a 1-based NCBI-like location string."""
    location = feature.location

    # Simple location, e.g. 195..644
    if not hasattr(location, "parts"):
        start = int(location.start) + 1
        end = int(location.end)
        return f"{start}..{end}"

    # Compound location, e.g. join(100..200,300..400)
    parts = []
    for part in location.parts:
        start = int(part.start) + 1
        end = int(part.end)
        parts.append(f"{start}..{end}")

    return "join(" + ",".join(parts) + ")"

def write_combined_sequence_txt(fasta_records: list[dict], output_path: Path) -> None:
    """Write all selected CDS sequences into one FASTA file."""
    with open(output_path, "w") as fasta_file:
        for record in fasta_records:
            fasta_file.write(record["header"] + "\n")

            sequence = record["sequence"]
            for i in range(0, len(sequence), 70):
                fasta_file.write(sequence[i:i + 70] + "\n")

            fasta_file.write("\n")

def get_target_genes() -> list[str]:
    """Prompt for genes entered with commas, on separate lines, or both."""
    print(
        "Enter target gene(s), separated by commas or one gene per line.\n"
        "Press Enter on a blank line when finished."
    )

    gene_lines = []
    while True:
        line = input()
        if not line.strip():
            break
        gene_lines.append(line)

    gene_input = "\n".join(gene_lines)
    genes = [
        gene.strip().upper()
        for line in gene_input.splitlines()
        for gene in line.split(",")
        if gene.strip()
    ]

    if not genes:
        raise ValueError("No target genes entered.")

    return genes


def process_candidates(
        candidates: list[dict],
        expected_gene: str,
        taxonomic_group: str,
        common_name: str,
        scientific_name: str,
        recovery_note: str = "",
) -> tuple[list[dict], list[dict]]:
    """
    Fetch and validate candidates for one gene/species search.
    Return tracking rows and FASTA records without writing global output yet.
    """
    tracking_rows = []
    fasta_records = []

    if not candidates:
        return tracking_rows, fasta_records

    candidates = prioritize_nm_over_xm(candidates)

    for candidate in candidates:
        accession = candidate["accession"]
        title = candidate["title"]
        sequence_type = candidate["record_type"] + " CDS"

        try:
            record = fetch_genbank(accession)
            cds_entries, found_cds_genes = extract_cds(record, expected_gene)

        except Exception as error:
            tracking_rows.append({
                "Taxonomic group": taxonomic_group,
                "Species (Common Name)": common_name,
                "Species (Scientific Name)": scientific_name,
                "Accession Numbers": accession,
                "Sequence Type": sequence_type,
                "CDS Length (bp)": "-",
                "Notes": f"error fetching/extracting CDS: {error}",
                "Title": title,
            })
            continue

        if not cds_entries:
            found_gene_text = (
                ", ".join(found_cds_genes)
                if found_cds_genes
                else "no CDS gene annotation found"
            )

            tracking_rows.append({
                "Taxonomic group": taxonomic_group,
                "Species (Common Name)": common_name,
                "Species (Scientific Name)": scientific_name,
                "Accession Numbers": accession,
                "Sequence Type": sequence_type,
                "CDS Length (bp)": "-",
                "Notes": (
                    "candidate found by NCBI search, but CDS gene annotation "
                    "does not match target gene; "
                    f"found CDS gene(s): {found_gene_text}"
                ),
                "Title": title,
            })
            continue

        for cds in cds_entries:
            variant_note = get_variant_note(title)
            notes = variant_note

            if cds["is_low_quality_protein"]:
                notes = add_note(notes, "LOW QUALITY PROTEIN")

            notes = add_note(notes, recovery_note)

            fasta_header = build_ncbi_style_header(
                common_name=common_name,
                cds=cds,
                variant_note=variant_note,
            )

            fasta_records.append({
                "header": fasta_header,
                "sequence": cds["sequence"],
            })

            tracking_rows.append({
                "Taxonomic group": taxonomic_group,
                "Species (Common Name)": common_name,
                "Species (Scientific Name)": scientific_name,
                "Accession Numbers": accession,
                "Sequence Type": sequence_type,
                "CDS Length (bp)": cds["cds_length"],
                "Notes": notes,
                "Title": title,
            })

    return tracking_rows, fasta_records


def main():
    target_genes = get_target_genes()

    input_dir = Path("input")
    base_output_dir = Path("output")

    species_file = input_dir / "Additional_sampling_species.xlsx"
    species_records = read_species_list(species_file)

    for target_gene in target_genes:
        print(f"\n==============================")
        print(f"Starting gene: {target_gene}")
        print(f"==============================")

        gene_output_dir = base_output_dir / target_gene
        gene_output_dir.mkdir(parents=True, exist_ok=True)

        tracking_rows = []
        fasta_records = []

        for species_record in species_records:
            taxonomic_group = species_record["taxonomic_group"]
            common_name = species_record["common_name"]
            scientific_name = species_record["scientific_name"]

            print(f"\nSearching {target_gene} in {common_name} ({scientific_name})")

            candidates = search_ncbi_candidates(
                gene=target_gene,
                species=scientific_name,
                max_results=20,
            )

            primary_rows, primary_fasta_records = process_candidates(
                candidates=candidates,
                expected_gene=target_gene,
                taxonomic_group=taxonomic_group,
                common_name=common_name,
                scientific_name=scientific_name,
            )

            if primary_fasta_records:
                tracking_rows.extend(primary_rows)
                fasta_records.extend(primary_fasta_records)
                continue

            print(
                "No CDS sequence was added by the original search. "
                "Trying NCBI Gene symbol resolution and GeneID links."
            )

            try:
                fallback_candidates, resolved_gene, fallback_note = (
                    get_alternative_gene_candidates(
                        requested_gene=target_gene,
                        species=scientific_name,
                        max_results=20,
                    )
                )
            except Exception as error:
                fallback_candidates = []
                resolved_gene = ""
                fallback_note = f"error checking NCBI Gene: {error}"

            fallback_rows, fallback_fasta_records = process_candidates(
                candidates=fallback_candidates,
                expected_gene=resolved_gene,
                taxonomic_group=taxonomic_group,
                common_name=common_name,
                scientific_name=scientific_name,
                recovery_note=fallback_note,
            )

            if fallback_fasta_records:
                tracking_rows.extend(fallback_rows)
                fasta_records.extend(fallback_fasta_records)
                continue

            if fallback_note:
                failure_note = fallback_note
            elif primary_rows:
                failure_note = primary_rows[0]["Notes"]
            else:
                failure_note = (
                    "no NM_/XM_ RefSeq mRNA candidates found; "
                    "no alternative gene name found through NCBI Gene"
                )

            tracking_rows.append({
                "Taxonomic group": taxonomic_group,
                "Species (Common Name)": common_name,
                "Species (Scientific Name)": scientific_name,
                "Accession Numbers": "-",
                "Sequence Type": "-",
                "CDS Length (bp)": "-",
                "Notes": failure_note,
                "Title": "-",
            })

        columns = [
            "Taxonomic group",
            "Species (Common Name)",
            "Species (Scientific Name)",
            "Accession Numbers",
            "Sequence Type",
            "CDS Length (bp)",
            "Notes",
            "Title",
        ]

        tracking_df = pd.DataFrame(tracking_rows, columns=columns)

        tracking_output_path = gene_output_dir / f"{target_gene}_tracking_sheet.csv"
        sequence_output_path = gene_output_dir / f"{target_gene}.txt"

        tracking_df.to_csv(tracking_output_path, index=False)
        write_combined_sequence_txt(fasta_records, sequence_output_path)

        print(f"\nDone with {target_gene}.")
        print(f"Tracking sheet saved to: {tracking_output_path}")
        print(f"Combined sequence text file saved to: {sequence_output_path}")

    print("\nAll requested genes complete.")

if __name__ == "__main__":
    main()
