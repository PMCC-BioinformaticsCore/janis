from janis import File


class JsonFile(File):
    @staticmethod
    def name():
        return "jsonFile"

    def doc(self):
        return "A JSON file file"
