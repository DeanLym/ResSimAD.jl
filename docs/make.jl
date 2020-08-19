using Documenter, ResSimAD

makedocs(
    sitename="ResSimAD.jl",
    format = Documenter.HTML(
        assets=["assets/favicon.ico", asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css)],
    ),
    pages = [
        "Home" => "index.md",
        "Basic workflow" => "workflow.md",
        "Advanced workflow" => "workflow2.md",
        "Precompile" => "precompile.md",
        "Benchmark" => "benchmark.md",
        "Python usage" => "python.md",
        "Example models" => "examples.md",
        "Input options" => "input.md",
        "API functions" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/DeanLym/ResSimAD.jl",
)
