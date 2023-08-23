# Project configuration, particularly for logging.

import logging

from ._version import __version__
from .seqcol import SeqColClient
from .seqcol import (
    parse_fasta,
    explain_flag,
    fasta_to_seqcol,
    build_sorted_name_length_pairs,
)

from .const import *
from .utilities import *


__classes__ = ["SeqColClient"]
__all__ = __classes__ + ["build_sorted_name_length_pairs"]
