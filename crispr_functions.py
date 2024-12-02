def find_PAM(target: str) -> int:
    """
    Finds the position of the first occurrence of the 'GG' PAM sequence in a DNA target sequence.

    Args:
        target (str): The DNA sequence to be searched. The input string must consist of valid DNA
                      nucleotides ('A', 'T', 'C', 'G').

    Returns:
        int: The index of the first occurrence of 'GG' in the target sequence after the 23rd position.

    Raises:
        ValueError: If the target sequence is invalid or if no 'GG' PAM sequence is found.
    """

    valid_nucleotides = {'A', 'T', 'C', 'G'}
    if any(char not in valid_nucleotides for char in target):
        raise ValueError("DNA sequence contains invalid characters. Allowed characters: A, T, C, G.")

    if not len(target) > 26:
        raise ValueError("Target sequence is not long enough to locate valid PAM")

    pam_index = target.find("GG", 23)

    if pam_index == -1:
      raise ValueError("No PAM sequence (GG) located in target string")

    return pam_index

def find_protospacer(target: str) -> str:
    """
    Identifies the protospacer sequence within a DNA target sequence based on the location of the PAM sequence.

    The protospacer is a 20 base pair sequence that immediately precedes the PAM sequence ('NGG'). This function
    locates the PAM sequence using the `find_PAM` function and then extracts the 20 bp protospacer sequence
    located just before it.

    Args:
        target (str): The DNA sequence to be searched. The input string must consist of valid DNA
                      nucleotides ('A', 'T', 'C', 'G').

    Returns:
        str: The 20 base pair protospacer sequence found immediately before the PAM sequence.

    Raises:
        ValueError: If the target sequence is invalid or if no valid PAM sequence is found.
    """

    # Find the index of the PAM sequence using the find_PAM function
    pam_index = find_PAM(target)
    protospacer = target[(pam_index - 21):(pam_index - 1)]

    return protospacer

if __name__ == "__main__":

    try:
        # Test out `find_PAM`
        print("Found PAM at location: " + str(find_PAM('CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC')))

        # Test out `find_protospacer`
        target = 'CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC'
        print("Found protospacer at location: " + str(find_protospacer(target)))
        print("Length of protospacer should be 20, is " + str(len(find_protospacer(target))))
    except Exception as e: 
        print(f"Error: {str(e)}")
