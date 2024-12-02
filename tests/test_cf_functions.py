import pytest
from cf_functions import cf_shorthand

import os
import sys
sys.path.append(os.getcwd() + "/PyDNA_CF_Simulator/")

from pydna_cf_simulator.construction_file import Step, ConstructionFile

def test_invalid_step():
    invalid_step = Step("SpillReagent", "mess")
    invalid_cf = ConstructionFile([invalid_step], [])
    with pytest.raises(ValueError):
        cf_shorthand(invalid_cf)
