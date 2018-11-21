from pipeline_definition.types.common_data_types import File


class RefFasta(File):

    @staticmethod
    def name():
        return "FastaRef"
