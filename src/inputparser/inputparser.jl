module InputParser

using DelimitedFiles
using DataFrames
using CSV

using Memento
const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)

using ..Grid:compute_cell_depth
using ..Facility:get_ctrl_mode, get_limit

include("utils.jl")
include("gridparser.jl")
include("rockparser.jl")
include("fluidparser.jl")
include("facilityparser.jl")
include("schedulerparser.jl")
include("nsolverparser.jl")
include("lsolverparser.jl")

function parse_input(options::Dict)
    info(LOGGER, "Parsing input options")
    keywords = keys(options)
    # Parse grid
    grid_opt = parse_grid(options, keywords)
    nc = grid_opt["nc"]
    # Parse rock
    rock_opt = parse_rock(options, keywords, nc)
    # Parse fluid
    fluid_opt = parse_fluid(options, keywords, nc)
    # Parse well
    facility_opt = parse_facility(options, keywords)
    # Parse scheduler
    scheduler_opt = parse_scheduler(options, keywords)
    # Parse nsolver
    nsolver_opt = parse_nsolver(options, keywords)
    # Parse lsolver
    lsolver_opt = parse_lsolver(options, keywords)
    parsed_options = Dict([
        ("grid_opt", grid_opt),
        ("rock_opt", rock_opt),
        ("fluid_opt", fluid_opt),
        ("facility_opt", facility_opt),
        ("scheduler_opt", scheduler_opt),
        ("nsolver_opt", nsolver_opt),
        ("lsolver_opt", lsolver_opt),
    ])
    return parsed_options
end

end
