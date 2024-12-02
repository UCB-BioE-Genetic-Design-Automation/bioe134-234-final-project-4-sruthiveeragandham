import pytest
import re
from genbank_functions import get_genbank_file

def test_invalid_refseq():
    with pytest.raises(ValueError):
        get_genbank_file("Invalid", "Invalid")

def test_valid_refseq():
    gb_filename = get_genbank_file("GCF_016861735.1", "AchevalieriM1_assembly01")
    print(gb_filename)
    assert re.findall(r"\.gb", gb_filename) is not None
