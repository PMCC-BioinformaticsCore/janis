from Pipeline import File


class VcfIdx(File):

    @staticmethod
    def name():
        return "VCFIDX"

    @staticmethod
    def secondary_files():
        return [".idx"]