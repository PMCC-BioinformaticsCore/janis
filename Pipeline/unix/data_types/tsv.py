from Pipeline import DataType, File


class Tsv(File):
    @staticmethod
    def name():
        return "tsv"

    def doc(self):
        return "A tab separated file"