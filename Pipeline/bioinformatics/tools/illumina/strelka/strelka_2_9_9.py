import yaml

from Pipeline.bioinformatics.tools.illumina.strelka.base import StrelkaBase


class Strelka_2_9_9(StrelkaBase):
    @staticmethod
    def docker():
        return "illusional/strelka"


if __name__ == "__main__":
    print(Strelka_2_9_9().help())
    cwl = Strelka_2_9_9().cwl(with_docker=True)
    print("\n\n" + yaml.dump(cwl))
