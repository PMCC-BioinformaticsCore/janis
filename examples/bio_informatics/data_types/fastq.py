from typing import Dict

from pipeline_definition.types.common_data_types import File
from pipeline_definition.types.data_types import DataType


class FastQ(File):
    @staticmethod
    def name():
        return "FASTQ"

    @staticmethod
    def doc():
        return "FASTQ files are text files containing sequence data with quality score, there are different types" \
               "with no standard: https://www.drive5.com/usearch/manual/fastq_files.html"