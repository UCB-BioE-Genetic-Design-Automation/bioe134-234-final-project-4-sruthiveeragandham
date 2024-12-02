from genbank_functions import get_genbank_file, extract_genes_info, gene_id_to_seq
from cf_functions import generate_aspergillus_ko_cf, cf_shorthand
from crispr_functions import find_protospacer

if __name__ == "__main__":
    genome_name = input("Enter NCBI Ref Seq for desired genome (press enter to use default): ")
    if genome_name != "":
        assembly_name = input("Enter assembly name for desired genome: ")
    else:
        genome_name = "GCF_001890905.1"
        assembly_name = "Aspac1"
        gene_id = "30974164"

    gb_file = get_genbank_file(genome_name, assembly_name)
    print("\nExtracted genbank file named " + gb_file)
    genes_info = extract_genes_info(gb_file)
    print("First 5 genes in genome GenBank file:")
    for idx, (locus_tag, gene_info) in enumerate(genes_info.items()):
        print(f"Locus Tag: {locus_tag}")
        print(f"Gene Name: {gene_info['gene_name']}")
        print(f"UTR: {gene_info['UTR']}")
        print(f"CDS: {gene_info['CDS']}")
        print(f"Gene ID: {gene_info['gene_id']}")
        print("-" * 40)

        if idx == 5:  # Stop after printing 5 entries
            break

    gene_id = input("Enter Gene ID for gene to knock out: ")
    ko_seq = gene_id_to_seq(gene_id, genes_info)

    protospacer = find_protospacer(ko_seq)
    print("\nProtospacer sequence is: " + str(protospacer))
    ko_cf = generate_aspergillus_ko_cf(protospacer)
    print("\nGenerated the following construction file to knock out GeneID:" + gene_id)
    print(cf_shorthand(ko_cf))
