from setuptools import setup, find_packages

setup(
    name="acruisepy",
    version="0.1.0",
    description="Provides tools for detecting peaks within ACRUISE measurements",
    url="https://github.com/wacl-york/acruise-peakid/acruisepy",
    author="Stuart Lacy",
    author_email="stuart.lacy@york.ac.uk",
    license="MIT",
    install_requires=[
        "pandas",
        "matplotlib",
        "numpy",
    ],
    packages=find_packages(include=["acruisepy", "acruisepy.*"]),
    zip_safe=False,
)
