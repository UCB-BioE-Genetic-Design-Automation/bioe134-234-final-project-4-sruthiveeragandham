import os
import sys
sys.path.append(os.getcwd() + "/PyDNA_CF_Simulator/")

from pydna_cf_simulator.construction_file import PCR, Digest, Ligate, ConstructionFile
from pydna_cf_simulator.polynucleotide import Polynucleotide, dsDNA, oligo, plasmid
from genbank_functions import plasmid_gb_to_seq
from crispr_functions import find_protospacer

def cf_shorthand(cf: ConstructionFile) -> str:
    """
    Given a ConstructionFile, returns string representation of the steps and sequences in
    standard construction file shorthand.

    Args:
        cf (ConstructionFile): construction file to convert to shorthand.

    Returns:
        str: string representation of construction file.

    Raises:
        ValueError: if the ConstructionFile input contains any unrecognized steps.
    """
    
    shorthand = ""
    shorthand += "Sequences:\n"
    for name, poly in cf.sequences:
        seq = poly.sequence[:40] + "..." if len(poly.sequence) > 40 else poly.sequence
        shorthand += f"{name: <15}{seq}\n"
    shorthand += "\nSteps:\n"
    for step in cf.steps:
        if step.operation == "PCR":
            shorthand += f"pcr {step.forward_oligo} {step.reverse_oligo} {step.template} {step.output}\n"
        elif step.operation == "Digest":
            shorthand += f"digest {step.dna} {str(step.enzymes)} {step.fragSelect} {step.output}\n"
        elif step.operation == "Ligate":
            shorthand += f"ligate {str(step.dnas)} {step.output}\n"
        elif step.operation == "GoldenGate":
            shorthand += f"golden_gate {str(step.dnas)} {step.enzyme} {step.output}\n"
        elif step.operation == "Gibson":
            shorthand += f"gibson {str(step.dnas)} {step.output}\n"
        elif step.operation == "Transform":
            shorthand += f"transform {step.dna} {step.strain} {str(step.antibiotics)} {step.temperature} {step.output}\n"
        else:
            raise ValueError("Invalid Step: " + str(step.operation))
    return shorthand

def generate_aspergillus_ko_cf(protospacer: str) -> ConstructionFile: 
    """
    Given a protospacer sequence, writes a ConstructionFile to create a plasmid that can be used
    to knock out a gene with a region that sufficiently matches the protospacer sequence, using
    pre-existing plasmids pFC334 and pFC330-333.

    Args:
        protospacer (str): protospacer sequence to insert before sgRNA scaffold.

    Returns:
        ConstructionFile: construction file with steps to produce knockout plasmid.
    """

    # List of (name: str, sequence: Polynucleotide) tuples representing sequences
    sequences = []
    # Add template plasmid sequences
    pFC334_seq = plasmid_gb_to_seq("seq_files/pFC334.gb")
    pFC330_333_seq = plasmid_gb_to_seq("seq_files/pFC330_333.gb")
    sequences.append(("pFC334", plasmid(pFC334_seq)))
    sequences.append(("pFC330-333", plasmid(pFC330_333_seq)))

    # List of Step objects
    steps = []

    # STEP 1A/1B: PCR with pFC334 template for fragments with correct protospacer
    # Fragment 1 forward primer matches start of gdpA promoter with a uracil residue
    frag1_fwd = oligo("GGGTTTAAU")
    sequences.append(("frag1-fwd", frag1_fwd))
    # Fragment 1 reverse primer matches HH ribozyme sequence up to overlap with Fragment 2 forward primer
    frag1_rev = oligo("AGGCGGGACTACTCAGGCACTCCTGCTTTGCUCATTCGA")
    sequences.append(("frag1-rev", frag1_rev))
    # Fragment 2 forward primer matches end of HH ribozyme sequence + protospacer
    frag2_fwd = oligo("AGTAAGCUCGTC" + protospacer)
    sequences.append(("frag2-fwd", frag2_fwd))
    # Fragment 2 reverse primer matches end of trpC terminator with a uracil residue
    frag2_rev = oligo("UAATTCTGG")
    sequences.append(("frag2-rev", frag2_rev))
    # Add one PCR step for each fragment
    steps.append(PCR("frag1-fwd", "frag1-rev", "pFC334", "frag1"))
    steps.append(PCR("frag2-fwd", "frag2-rev", "pFC334", "frag2"))

    # STEP 2A/2B: USER digest of PCR fragments with UDG and Endonuclease VIII
    # Removes uracil residues and creates sticky ends
    steps.append(Digest("frag1", ["UDG", "Nei"], 1, "frag1-sticky"))
    steps.append(Digest("frag2", ["UDG", "Nei"], 1, "frag2-sticky"))

    # STEP 3 (independent of STEP 1 and STEP 2): PacI/Nt.BbvCI digest of pFC330-333
    steps.append(Digest("pFC330-333", ["PacI", "Nt.BbvCI"], 1, "pFC330-dig")) 

    # STEP 4: USER fusion of digested pFC330-333 with sticky ends + PCR fragments
    steps.append(Ligate(["frag1-sticky", "frag2-sticky", "pFC330-dig"], "ko-plasmid"))

    return ConstructionFile(steps, sequences)

if __name__ == "__main__":
    target = 'CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC'
    protospacer = find_protospacer(target)
    ko_cf = generate_aspergillus_ko_cf(protospacer)
    print(cf_shorthand(ko_cf))
