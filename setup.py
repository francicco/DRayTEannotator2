from setuptools import setup, find_packages

setup(
    name="drayte",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "matplotlib",
    ],
    entry_points={
        "console_scripts": [
            "drayte=drayte.pipeline.runner:main",
            "drayte-summary=drayte.reporting.summary_cli:main",
        ],
    },
)
