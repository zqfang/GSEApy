#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# All project metadata is in pyproject.toml [project].
# This file exists solely for the Rust extension (setuptools_rust).
import os
import pathlib

from setuptools import setup
from setuptools.command.sdist import sdist as SdistCommand
from setuptools_rust import Binding, RustExtension


class CargoModifiedSdist(SdistCommand):
    def make_release_tree(self, base_dir, files):
        """Include Rust sources in the sdist."""
        files.append("Cargo.toml")
        files += [str(f) for f in pathlib.Path("src").glob("**/*.rs") if f.is_file()]
        super().make_release_tree(base_dir, files)


setup(
    rust_extensions=[
        RustExtension(
            "gseapy.gse",
            "Cargo.toml",
            debug="DEBUG" in os.environ,
            binding=Binding.PyO3,
        )
    ],
    cmdclass={"sdist": CargoModifiedSdist},
    zip_safe=False,
)
