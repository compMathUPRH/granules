import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="granules",
    version="1.0.0",
    author="Lemuel I. Rivera Cantú",
    author_email="lemuel.rivera6@upr.edu",
    description="testing package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/compMathUPRH/granules.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['pylint>=','GraphVIZ>='],
)
