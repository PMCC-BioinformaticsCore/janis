from abc import ABC

from Pipeline import File


class Vcf(File):
    @staticmethod
    def name():
        return "VCF"

    def doc(self):
        return """
    Variant Call Format:

    The Variant Call Format (VCF) specifies the format of a text file 
    used in bioinformatics for storing gene sequence variations. 

    Documentation: https://samtools.github.io/hts-specs/VCFv4.3.pdf
    """.strip()


class VcfIdx(Vcf):

    @staticmethod
    def name():
        return "VCFIDX"

    @staticmethod
    def secondary_files():
        return [".idx"]


class CompressedVcf(File):
    @staticmethod
    def name():
        return "compressed-vcf-gz"

    def doc(self):
        return ".vcf.gz"


class TabixIdx(CompressedVcf):
    @staticmethod
    def name():
        return "vcf-gz-tbi"

    @staticmethod
    def secondary_files():
        return [".tbi"]

    def doc(self):
        return ".vcf.gz with .vcf.gz.tbi file"


class GVCF(Vcf):
    @staticmethod
    def name():
        return "gVCF"

    def doc(self):
        return """
Section 5.5: Representing unspecified alleles and REF only blocks (gVCF)
Documentation: https://samtools.github.io/hts-specs/VCFv4.3.pdf
        """