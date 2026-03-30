#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Metadata is declared in pyproject.toml [project].
# This file only configures the Rust extension and sdist command.
import os
import pathlib

from setuptools import setup
from setuptools.command.sdist import sdist as SdistCommand
from setuptools_rust import Binding, RustExtension


class CargoModifiedSdist(SdistCommand):
    def make_release_tree(self, base_dir, files):
        """Stages the files to be included in archives"""
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
