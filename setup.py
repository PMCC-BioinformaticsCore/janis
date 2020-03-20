from setuptools import setup, find_packages

# Version information is found in the __init__ file of `janis/`
DESCRIPTION = "Contains classes and helpers to build a workflow, and provide options to convert to CWL / WDL"

JANIS_CORE_VERSION = "v0.9.9"
JANIS_ASSISTANT_VERSION = "v0.9.11"
JANIS_UNIX_VERSION = "v0.9.0"
JANIS_BIOINFORMATICS_VERSION = "v0.9.6"
JANIS_PIPELINES_VERSION = "v0.9.3"
JANIS_TEMPLATES_VERSION = "v0.9.4"


######## SHOULDN'T NEED EDITS BELOW THIS LINE ########

min_core_version = f"janis-pipelines.core>=" + JANIS_CORE_VERSION
min_assistant_version = f"janis-pipelines.runner>=" + JANIS_ASSISTANT_VERSION
min_unix_version = f"janis-pipelines.unix>=" + JANIS_UNIX_VERSION
min_bioinf_version = f"janis-pipelines.bioinformatics>=" + JANIS_BIOINFORMATICS_VERSION
min_pipes_version = f"janis-pipelines.pipelines>=" + JANIS_PIPELINES_VERSION
min_templs_version = f"janis-pipelines.templates>=" + JANIS_TEMPLATES_VERSION

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
    install_requires=[
        min_core_version,
        min_assistant_version,
        min_unix_version,
        min_bioinf_version,
        min_pipes_version,
        min_templs_version,
    ],
    extras_require={"bioinformatics": [min_bioinf_version, min_pipes_version]},
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
