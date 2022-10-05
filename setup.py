import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bmdms-np", # Replace with your own username
    version="0.0.1",
    author="Sangwon Lee",
    author_email="sw.lee@yonsei.ac.kr",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: AGPL-3.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)