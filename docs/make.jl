using Documenter, ResSimAD

makedocs(
    sitename="ResSimAD.jl",
    format = Documenter.HTML(
        assets=["assets/favicon.ico", asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css)],
    ),
    pages = [
        "Home" => "index.md",
        "Basic workflow" => "workflow.md",
        "Example Models" => "examples.md",
        "Input Options" => "input.md",
        "Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/DeanLym/ResSimAD.jl",
)
