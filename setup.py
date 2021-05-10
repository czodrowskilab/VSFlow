from setuptools import setup

with open("README.md", "r") as readme:
    long_description = readme.read()

setup(
    name="vsflow",
    version="1.0.1",
    description="Virtual Screening Workflow: substructure-, fingerprint- and shape-based screening",
    author="Sascha Jung",
    author_email="sascha.jung@tu-dortmund.de",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/czodrowskilab/VSFlow",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    keywords="virtual screening, substructure, fingerprints, shape, ligand, cheminformatics",
    packages=["vslib"],
    package_data={"vslib": ["resources/BaseFeatures.fdef", "resources/MinimalFeatures.fdef", "resources/DejaVuSansMono.ttf"]},
    scripts=["vsflow"],
    python_requires=">=3.7.*"
)
