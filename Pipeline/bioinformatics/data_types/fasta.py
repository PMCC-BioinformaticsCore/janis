from Pipeline import File


class Fasta(File):

    @staticmethod
    def name():
        return "Fasta"

    @staticmethod
    def secondary_files():
        return [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]
