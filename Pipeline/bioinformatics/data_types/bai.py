from Pipeline import File


class Bai(File):

    @staticmethod
    def name():
        return "BAI"

    def doc(self):
        return "Indexed BAM file (https://www.biostars.org/p/15847/), http://software.broadinstitute.org/software/igv/bam"
