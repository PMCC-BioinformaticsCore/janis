from Pipeline import DataType, File


class Tsv(File):
    @staticmethod
    def name():
        return "tsv"

    @staticmethod
    def doc():
        return "A tab separated file"