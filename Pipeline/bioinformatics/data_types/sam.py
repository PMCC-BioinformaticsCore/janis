from Pipeline import File


class Sam(File):

    @staticmethod
    def name():
        return "SAM"

    def doc(self):
        return "Tab-delimited text file that contains sequence alignment data"
