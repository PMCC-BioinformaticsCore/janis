from Pipeline.bioinformatics.tools.bcftools.bcftoolsbase import BcfToolsBase


class BcfTools_1_5(BcfToolsBase):
    @staticmethod
    def docker():
        # Todo: Create a docker container with 1_9
        raise Exception("A docker container that contains bcftools v1.9 must be created in order to use this tool")
        # return "biocontainers/bcftools:v1.5_cv2"
