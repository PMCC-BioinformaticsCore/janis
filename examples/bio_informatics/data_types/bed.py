from pipeline_definition.types.common_data_types import File
from pipeline_definition.types.data_types import DataType


class Bed(File):
    @staticmethod
    def name():
        return "bed"
