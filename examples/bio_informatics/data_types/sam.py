from Pipeline import File


class Sam(File):

    @staticmethod
    def name():
        return "SAM"

    @staticmethod
    def doc():
        return "Tab-delimited text file that contains sequence alignment data"
