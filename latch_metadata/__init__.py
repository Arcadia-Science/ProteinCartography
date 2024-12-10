from latch.types.metadata import SnakemakeMetadata, LatchAuthor, EnvironmentConfig
from latch.types.directory import LatchDir

from .parameters import generated_parameters, file_metadata

SnakemakeMetadata(
    output_dir=LatchDir("latch:///your_output_directory"),
    display_name="Your Workflow Name",
    author=LatchAuthor(
        name="Your Name",
    ),
    env_config=EnvironmentConfig(
        use_conda=False,
        use_container=False,
    ),
    cores=4,
    # Add more parameters
    parameters=generated_parameters,
    file_metadata=file_metadata,

)
