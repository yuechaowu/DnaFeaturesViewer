from setuptools import setup, find_packages

version = {}
with open("dna_features_viewer/version.py") as fp:
    exec(fp.read(), version)


setup(
    name="dna_features_viewer",
    version=version["__version__"],
    author="Zulko",
    description="Plot features from DNA sequences (e.g. Genbank) with Python",
    long_description=open("pypi-readme.rst").read(),
    url="https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer",
    license="MIT",
    keywords="DNA Sequence Feature Genbank Biopython Matplotlib",
    packages=find_packages(exclude="docs"),
    install_requires=["matplotlib>=3", "Biopython", "packaging"],
)
