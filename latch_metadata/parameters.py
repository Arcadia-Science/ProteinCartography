from dataclasses import dataclass
import typing
import typing_extensions

from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import SnakemakeParameter, SnakemakeFileParameter, SnakemakeFileMetadata
from latch.types.file import LatchFile
from latch.types.directory import LatchDir


# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters, file_metadata

generated_parameters = {
    "mode": SnakemakeParameter(
        display_name="Mode",
        type=str,
        default="search",
    ),
    "analysis_name": SnakemakeParameter(
        display_name="Analysis Name",
        type=str,
        default="actin-demo-small",
    ),
    "input_dir": SnakemakeParameter(
        display_name="Input Dir",
        type=LatchDir,
    ),
    "output_dir": SnakemakeParameter(
        display_name="Output Dir",
        type=LatchDir,
    ),
    "max_blast_hits": SnakemakeParameter(
        display_name="Max Blast Hits",
        type=int,
        default=10,
    ),
    "max_foldseek_hits": SnakemakeParameter(
        display_name="Max Foldseek Hits",
        type=int,
        default=10,
    ),
    "max_structures": SnakemakeParameter(
        display_name="Max Structures",
        type=int,
        default=10,
    ),
    "plotting_modes": SnakemakeParameter(
        display_name="Plotting Modes",
        type=typing.List[str],
        default=["pca_umap"],
    ),
}

file_metadata = {
    "input_dir": SnakemakeFileMetadata(
        path="latch://36681.account/2024-proteincartography-on-latch-demo-1733788167.0832448/demo/search-mode/input/",
        config=True,
    ),
    "output_dir": SnakemakeFileMetadata(
        path="latch://36681.account/2024-proteincartography-on-latch-demo-1733788167.0832448/demo/search-mode/output/",
        config=True,
    ),
}
