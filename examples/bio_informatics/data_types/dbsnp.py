from Pipeline import File


class Dbsnp(File):
    @staticmethod
    def name():
        return "DBSNP"

    @staticmethod
    def doc():
        pass

    # @classmethod
    # def schema(cls) -> Dict:
    #     return {
    #         "path": {"type": "string", "required": True}
    #     }
