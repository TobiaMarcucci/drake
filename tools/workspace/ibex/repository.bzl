# -*- mode: python -*-

load("@drake//tools/workspace:github.bzl", "github_archive")

def ibex_repository(
        name,
        mirrors = None):
    github_archive(
        name = name,
        repository = "ibex-team/ibex-lib",
#        commit = "ibex-2.8.9",
#        sha256 = "fee448b3fa3929a50d36231ff2f14e5480a0b82506594861536e3905801a6571",  # noqa
        commit = "ibex-2.7.4",
        sha256 = "2b32a1e51766476c9baab4017001e3e9ce2b6b102e2ac7492002305483607036",  # noqa
        build_file = "@drake//tools/workspace/ibex:package.BUILD.bazel",
        mirrors = mirrors,
    )
