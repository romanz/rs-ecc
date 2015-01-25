#!/usr/bin/env python
from setuptools import setup
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):

    def finalize_options(self):
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import sys
        import pytest
        sys.exit(pytest.main(['tests']))

setup(
    name="ecc",
    version="0.1",
    description="Error Correction Code Library",
    author="Roman Zeyde",
    author_email="roman.zeyde@gmail.com",
    license="MIT",
    url="http://github.com/romanz/ecc",
    packages=['ecc'],
    tests_require=['pytest'],
    cmdclass={'test': PyTest},
    platforms=['POSIX'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Information Technology",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 2.7",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: System :: Networking",
        "Topic :: Communications",
    ],
)
