from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="autocbei",
    version="1.8",
    packages=find_packages(exclude=["cbeiUnitTests"]),
    description = '"autoCBEI.py" can calculate the potential CBEI loci of the target CDSs and perform statistics and plots.',
    author = "Haopeng Yu",
    author_email = "atlasbioin4@gmail.com",
    url = "https://github.com/atlasbioinfo/CRISPR-CBEI/tree/master/autoCBEI",
    python_requires='>=3.7.0',
    license='Apache-2.0',
    entry_points = {
        'console_scripts' : [
            'autocbei = autocbei:mainCBEI'
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    install_requires = ["biopython","matplotlib"],
    package_data = {
        '':['*.fa']
    }
)
