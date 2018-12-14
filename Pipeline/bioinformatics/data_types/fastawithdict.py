from Pipeline.bioinformatics.data_types.fasta import Fasta


class FastaWithDict(Fasta):

    @staticmethod
    def name():
        return "FastaWithDict"

    @staticmethod
    def secondary_files():
        return [*super().secondary_files(), "^.dict"]
