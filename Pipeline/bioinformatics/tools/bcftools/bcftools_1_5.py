from Pipeline.bioinformatics.tools.bcftools.bcftoolsbase import BcfToolsBase


class BcfTools_1_5(BcfToolsBase):
    @staticmethod
    def docker():
        # Todo: Create a docker container with 1_9
        return "biocontainers/bcftools:v1.5_cv2"