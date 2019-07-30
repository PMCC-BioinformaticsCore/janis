from setuptools import setup, find_packages

# Version information is found in the __init__ file of `janis/`
DESCRIPTION = "Contains classes and helpers to build a workflow, and provide options to convert to CWL / WDL"

######## SHOULDN'T NEED EDITS BELOW THIS LINE ########

with open("./README.md") as readme:
    long_description = readme.read()

vsn = {}
with open("./janis/__meta__.py") as fp:
    exec(fp.read(), vsn)
__version__ = vsn["__version__"]
githuburl = vsn["GITHUB_URL"]


setup(
    name="janis pipelines",
    version=__version__,
    description=DESCRIPTION,
    url=githuburl,
    author="Michael Franklin, Richard Lupat",
    author_email="michael.franklin@petermac.org",
    license="GNU",
    keywords=["pipelines", "bioinformatics"],
    packages=["janis"]
    + ["janis." + p for p in sorted(find_packages("./janis"))]
    + ["toolbuilder"],
    install_requires=["janis-pipelines.core>=0.5.0", "janis-pipelines.unix"],
    extras_require={
        "bioinformatics": "janis-pipelines.bioinformatics>=0.5.0",
        "runner": "janis-pipelines.runner>=0.5.0",
    },
    entry_points={"console_scripts": ["janisbuilder=toolbuilder.main:process_args"]},
    zip_safe=False,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
