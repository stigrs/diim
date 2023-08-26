# Copyright (c) 2023 Stig Rune Sellevag
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as fh:
    requirements = fh.readlines()

setuptools.setup(
    name="pydiim",
    version="0.2.0",
    author="Stig Rune Sellevag",
    author_email="stig-rune.sellevag@ffi.no",
    license="MIT License",
    description="Dynamic Inoperability Input-Output Model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="git@github.com:stigrs/diim.git",
    packages=setuptools.find_packages(),
    install_requires=[req for req in requirements if req[:2] != "# "],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
)
