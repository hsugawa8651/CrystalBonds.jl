# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara

using Documenter
using CrystalBonds

makedocs(;
    sitename = "CrystalBonds.jl",
    modules = [CrystalBonds],
    pages = [
        "Home" => "index.md",
        "Quick Start" => "guide.md",
        "Theory" => "theory.md",
        "Examples" => "examples.md",
        "Generating _hr.dat" => "generating_hr_dat.md",
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    repo = "https://github.com/hsugawa8651/CrystalBonds.jl",
    warnonly = false,
)

deploydocs(;
    repo = "github.com/hsugawa8651/CrystalBonds.jl",
    devbranch = "main",
)
