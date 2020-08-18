function default_log_format()
    fmt = DefaultFormatter("[{level}]: {msg}")
    opt = Dict{Symbol, Any}(
        :colors => Dict{AbstractString, Symbol}(
            "trace" => :normal,
            "debug" => :blue,
            "info" => :normal,
            "notice" => :cyan,
            "warn" => :magenta,
            "error" => :red,
            "critical" => :yellow,
            "alert" => :white,
            "emergency" => :red,
        )
    )
    return fmt, opt
end

function __init__()
    logger = getlogger("root")
    fmt, opt = default_log_format()
    logger.handlers["console"] = DefaultHandler(stdout, fmt, opt)
end

"""
    add_log_file(fn)

Direct log to file `fn`

# Examples
```jldoctest
julia> using ResSimAD:add_log_file

julia> add_log_file("sim.log");

```

"""
function add_log_file(fn)
    fmt, opt = default_log_format()
    push!(getlogger("root"), DefaultHandler(fn, fmt, opt))
end

"""
    silence()

Do not show any message except for "error" messages. See `ResSimAD.verbose()` for more information.

# Examples
```jldoctest
julia> using ResSimAD:silence

julia> silence();

```

"""
function silence()
    setlevel!(getlogger("root"), "error"; recursive=true)
end

"""
    debug()

Show "debug" messages. See `ResSimAD.verbose()` for more information.

# Examples
```jldoctest
julia> using ResSimAD:debug

julia> debug();

```

"""
function debug()
    setlevel!(getlogger("root"), "debug"; recursive=true)
end

"""
    verbose(level::String)

Set verbose level. Show messages with higher level than the set `level`.
See Memento.jl for more information.

# Arguments
- `level::String`: verbose levels (low to high)
    - `"debug"`: verbose message used for debugging
    - `"info"`: general information
    - `"notice"`: important events that are part of normal execution
    - `"warn"`: warning that may cause the program to fail
    - `"error"`: error that causes the program to terminate

# Examples
```jldoctest
julia> using ResSimAD:verbose

julia> verbose("debug");

```

"""
function verbose(level::String)
    setlevel!(getlogger("root"), level; recursive=true);
end
