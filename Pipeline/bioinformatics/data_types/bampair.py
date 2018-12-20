from Pipeline import File
from Pipeline.bioinformatics.data_types.bam import Bam


class BamPair(File):

    @staticmethod
    def name():
        return "BamPair"

    @staticmethod
    def secondary_files():
        return ["^.bai"]

    @staticmethod
    def doc():
        return "A Bam and bai as the secondary"

