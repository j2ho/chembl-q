# chembl-curator/
# Complete package structure for ChEMBL bioactivity curation

# ============================================================================
# setup.py
# ============================================================================

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="chembl-curator",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A Python package for curating ChEMBL bioactivity data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/username/chembl-curator",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "chembl-curator=chembl_curator.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "chembl_curator": ["data/*.txt", "config/*.yaml"],
    },
)

