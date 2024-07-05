using RayTracing
using Documenter

DocMeta.setdocmeta!(RayTracing, :DocTestSetup, :(using RayTracing); recursive=true)

makedocs(;
    modules=[RayTracing],
    authors="Jessica Kent, Ted Marku",
    sitename="RayTracing.jl",
    format=Documenter.HTML(;
        canonical="https://jessica-kent.github.io/RayTracing.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jessica-kent/RayTracing.jl",
    devbranch="main",
)
