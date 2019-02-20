from setuptools import setup, find_packages
from janis import __version__
from janis.utils.janisconstants import GITHUB_URL

DESCRIPTION = "Contains classes and helpers to build a workflow, and provide options to convert to CWL / WDL"

######## SHOULDN'T NEED EDITS BELOW THIS LINE ########

with open("./README.md") as readme:
    long_description = readme.read()

setup(
    name="janis pipelines",
    version=__version__,
    description=DESCRIPTION,
    url=GITHUB_URL,
    author="Michael Franklin, Evan Thomas, Mohammad Bhuyan",
    author_email="michael.franklin@petermac.org",
    license="GNU",
    packages=["janis"] + ["janis." + p for p in sorted(find_packages('./janis'))],
    install_requires=["networkx>=2.1", "ruamel.yaml >= 0.12.4, <= 0.15.77", "six", "illusional.wdlgen",
                      "illusional.cwlgen>=0.0.5"],
    zip_safe=False,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    extras_require={
        'bioinformatics': ['janis-pipelines.bioinformatics']
    }
)
