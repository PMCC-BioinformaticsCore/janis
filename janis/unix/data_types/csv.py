from janis import File


class Csv(File):
    @staticmethod
    def name():
        return "csv"

    def doc(self):
        return "A comma separated file"
