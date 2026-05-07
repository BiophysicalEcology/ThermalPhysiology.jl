using Documenter
using ThermalPhysiology
using Unitful

DocMeta.setdocmeta!(ThermalPhysiology, :DocTestSetup,
    :(using ThermalPhysiology, Unitful); recursive=true)

makedocs(
    sitename = "ThermalPhysiology.jl",
    authors  = "Michael Kearney",
    modules  = [ThermalPhysiology],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://BiophysicalEcology.github.io/ThermalPhysiology.jl/stable/",
    ),
    pages = [
        "Overview"  => "index.md",
        "Models"    => "models.md",
        "Fitting"   => "fitting.md",
        "API"       => "api.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(
    repo   = "github.com/BiophysicalEcology/ThermalPhysiology.jl.git",
    target = "build",
    push_preview = true,
)
