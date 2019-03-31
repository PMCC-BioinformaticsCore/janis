from janis.types.common_data_types import File


class ZipFile(File):
    @staticmethod
    def name():
        return "Zip"

    def doc(self):
        return "A zip archive, ending with .zip"
