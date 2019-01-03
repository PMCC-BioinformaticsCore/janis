from abc import ABC

from Pipeline import File


class VcfBase(File, ABC):
    @staticmethod
    def name():
        return "VCF"

    @staticmethod
    def doc():
        return """
    Variant Call Format:

    The Variant Call Format (VCF) specifies the format of a text file 
    used in bioinformatics for storing gene sequence variations. 

    Documentation: https://samtools.github.io/hts-specs/VCFv4.3.pdf
    """.strip()


class VcfIdx(VcfBase):

    @staticmethod
    def name():
        return "VCFIDX"

    @staticmethod
    def secondary_files():
        return [".idx"]


class TabixIdx(VcfBase):
    @staticmethod
    def name():
        return "vcf-gz-tbi"

    @staticmethod
    def secondary_files():
        return [".tbi"]

    @staticmethod
    def doc():
        return ".vcf.gz with .vcf.gz.tbi file"


class GVCF(VcfBase):
    @staticmethod
    def name():
        return "gVCF"

    @staticmethod
    def doc():
        return """
Section 5.5: Representing unspecified alleles and REF only blocks (gVCF)
Documentation: https://samtools.github.io/hts-specs/VCFv4.3.pdf
        """