from Pipeline import File


class FastQ(File):
    @staticmethod
    def name():
        return "FASTQ"

    @staticmethod
    def doc():
        return "FASTQ files are text files containing sequence data with quality score, there are different types" \
               "with no standard: https://www.drive5.com/usearch/manual/fastq_files.html"