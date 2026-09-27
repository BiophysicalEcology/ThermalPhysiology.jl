using Documenter
using DocumenterVitepress
using ThermalPhysiology
using CairoMakie
using Unitful

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

# Helpers for the figures, loaded in the examples with `using Main.FigureHelpers`
include("figure_helpers.jl")

makedocs(
    modules = [ThermalPhysiology],
    sitename = "ThermalPhysiology.jl",
    authors = "Michael Kearney et al.",
    clean = true,
    doctest = false,
    checkdocs = :exports,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/ThermalPhysiology.jl", # this must be the full URL!
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/BiophysicalEcology/ThermalPhysiology.jl",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
