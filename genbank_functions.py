from Bio import SeqIO, BiopythonParserWarning
from collections import defaultdict
import urllib.request
import gzip
import warnings
from typing import Dict, Optional


def get_genbank_file(ref_seq: str, assembly: str) -> str:
    """
    Download a compressed GenBank file from the NCBI FTP server and extract it.

    Args:
        ref_seq (str): NCBI annotated reference sequence identifier (e.g., 'GCF_000001405').
        assembly (str): Assembly version (e.g., '39').

    Returns:
        str: Path to the unzipped extracted GenBank file.

    Raises:
        ValueError: if FTP request to NCBI fails with given parameters.
    """

    ftp_header = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
    digits = ref_seq[4:7] + "/" + ref_seq[7:10] + "/" + ref_seq[10:13] + "/"
    file_name = ref_seq + "_" + assembly + "/" + ref_seq + "_" + assembly + "_genomic"
    file_ext = ".gbff.gz"

    ftp_url = ftp_header + digits + file_name + file_ext
    try:
        gb_zip, _ = urllib.request.urlretrieve(ftp_url)
    except urllib.error.HTTPError as e:
        raise ValueError("Unable to retrieve genome with RefSeq " + ref_seq + " and assembly version " + assembly)
    with open(ref_seq + "_" + assembly + ".gb", "wb") as f_unzip:
        with gzip.open(gb_zip, "rb") as f_zip:        
            f_unzip.write(f_zip.read())
    return ref_seq + "_" + assembly + ".gb"

def extract_genes_info(genbank_file: str) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Extract gene information (name, NCBI gene ID, UTR, CDS) for every gene in a GenBank file.

    Args:
        genbank_file (str): Path to the GenBank file to parse.

    Returns:
        Dict[str, Dict[str, Optional[str]]]: Dictionary with locus tags as keys and another dictionary as values.
            Each inner dictionary contains:
                - "gene_name" (Optional[str]): The name of the gene.
                - "gene_id" (Optional[str]): The NCBI GeneID extracted from db_xref.
                - "UTR" (Optional[str]): The untranslated region sequence.
                - "CDS" (Optional[str]): The coding sequence.
    """

    gene_dict = defaultdict(dict)  # Dictionary to store gene info
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                gene_name = feature.qualifiers.get("gene", [None])[0]
                db_xref = feature.qualifiers.get("db_xref", [None])
                gene_id = None
                for info in db_xref:
                    if "GeneID" in info:
                        gene_id = info.replace("GeneID:", "")
                        break

                # CDS information
                cds_feature = None
                for cds in record.features:
                    if cds.type == "CDS" and cds.qualifiers.get("locus_tag") == [locus_tag]:
                        cds_feature = cds
                        break

                if cds_feature:
                    start, end = cds_feature.location.start, cds_feature.location.end
                    strand = cds_feature.location.strand
                    if strand == 1:  # Forward strand
                        utr_start = max(0, start - 50)
                        utr_seq = record.seq[utr_start:start]
                    else:  # Reverse strand, we need to reverse complement
                        utr_start = end
                        utr_seq = record.seq[utr_start:utr_start + 50].reverse_complement()

                    cds_seq = cds_feature.extract(record.seq)
                    # Save the gene information in the dictionary
                    gene_dict[locus_tag] = {
                        "gene_name": gene_name,
                        "gene_id": gene_id,
                        "UTR": utr_seq,
                        "CDS": cds_seq
                    }
    return gene_dict

def gene_id_to_seq(gene_id: str, genes_info: dict) -> str:
    """
    Given an NCBI gene ID and a `genes_info` dictionary, return the sequence for that gene.

    Args:
        gene_id (str): NCBI gene ID, e.g. "69578261" for EcoRI in Bifidobacterium longum.
        genes_info(dict): dictionary mapping locus tags to a dictionary that contains a key-value
                          pair with key `gene_id` and value of NCBI Gene ID as a string.

    Returns:
        str: DNA sequence for the requested gene, empty string if gene is not found.

    Raises:
        ValueError: if gene with provided gene_id is not located in genes_info.
    """

    seq = ""
    for locus in genes_info.keys():
        if genes_info[locus]["gene_id"] == gene_id:
            seq = genes_info[locus]["CDS"]
            break
    if seq == "":
        raise ValueError("Gene with gene_id " + gene_id + " not located in genome")
    return seq

def plasmid_gb_to_seq(plasmid_fp: str) -> str:
    """
    Given a filepath to a GenBank file containing a plasmid sequence, return the sequence.

    Args:
        plasmid_fp (str): filepath for GenBank file with plasmid sequence.

    Returns:
        str: full DNA sequence of plasmid from GenBank file.
    """

    # Dangerous to suppress warnings but NCBI has many GenBank files with too-long locus names
    warnings.simplefilter("ignore", BiopythonParserWarning)
    for record in SeqIO.parse(plasmid_fp, "genbank"):
        return(record.seq)

if __name__ == "__main__":

    gb_file = get_genbank_file("GCF_001890905.1", "Aspac1")
    print("Extracted genbank file named " + gb_file)

    genes_info = extract_genes_info(gb_file)
    for idx, (locus_tag, gene_info) in enumerate(genes_info.items()):
        print(f"Locus Tag: {locus_tag}")
        print(f"Gene Name: {gene_info['gene_name']}")
        print(f"UTR: {gene_info['UTR']}")
        print(f"CDS: {gene_info['CDS']}")
        print(f"Gene ID: {gene_info['gene_id']}")
        print("-" * 40)

        if idx == 5:  # Stop after printing 5 entries
            break

    expected_seq = "ATGGTCCGGGGTATTATCTTCCAGCTCGTATTCGTGGGACTGGTCGTCGACTTCGTGCTGCGCCTGTCCAAAAGGGGCCGGCCGCAGGCCCTGTGGTCCAACCGGCCCCTGCTTTTGCTTGGTGGCGCGACCGCGCTCTCGCTGCTACTGATCTATATTCGGAGTGTGTACCGGACGATTGAGCTGCTACACGGATGGACGAGCTCGACGATGCATAACGAGATGTTGCTGATCGGGCTGGACGGCGCGATTATGGTCCCCGCCTTTTCGGTTTATAACCTGCTTCATCCGGGGTACCTGTTGCCGAAGGTACAGCGGGAGGTCGGGTACCTTGACGCTCGGGGGTTGCAGATGGAAATGGAGGCGGTTAAGGATGGTCGGTACCAGATAGTCGAGGTGGAGGGTGGAAAGGATGACTCTGTTACGTTGCTGGACAGGAGCTGA"
    ret_seq = gene_id_to_seq("30974164", genes_info)
    print("\nRetrieved sequence for ID 30974164 matches expected? " + str(ret_seq == expected_seq))

    plasmid_seq = plasmid_gb_to_seq("seq_files/pFC330_333.gb")
    print("\nStart of pFC330_333 plasmid sequence: " + plasmid_seq[:40] + "...")
