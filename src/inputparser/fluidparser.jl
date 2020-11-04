
function parse_fluid(options, keywords, nc)
    fluid_opt = Dict()
    fluid_opt["type"] = options["fluid"]
    if fluid_opt["type"] == "OW"
        return parse_ow_fluid(options, keywords, nc, fluid_opt)
    else
        error(LOGGER, "Unsupported fluid type $(fluid_opt["type"])")
    end
end


function parse_ow_fluid(options, keywords, nc, fluid_opt)
    # Phase density
    v = ("ρo", "ρw")
    values = (49.1, 64.79)
    for (p, val) in zip(v, values)
        if p ∈ keywords
            fluid_opt[p] = options[p]
        else
            notice(LOGGER, "Keyword \"$p\" missing (default value $val will be used)")
            fluid_opt[p] = val
        end
    end
    pvct_keys = ("pref", "bref",  "c", "μref", "cμ")
    # PVTO
    v = ("PVDO", "PVCDO")
    p = check_complementary_keywords(v, keywords)
    val = options[p]
    fluid_opt["PVTO_TYPE"] = p
    if p == "PVDO"
        check_keyword_type(p, options, (String, DataFrame))
        if isa(val, String)
            info(LOGGER, "Reading $p table from $val")
            fluid_opt["PVTO"] = read_pvt_table(val, p)
        else
            fluid_opt["PVTO"] = val
        end
    else  #"PVCDO"
        check_keyword_type(p, options, (String, DataFrame, Dict))
        if isa(val, String)
            info(LOGGER, "Reading $p table from $val")
            fluid_opt["PVTO"] = read_pvtc_table(val, p)
        elseif isa(val, Dict)
            #
            for v in pvct_keys
                if !(v ∈ keys(val))
                    error(LOGGER, "Missing key $v in the PVCDO input dictionary" )
                end
            end
            fluid_opt["PVTO"] = val
        else
            for v in pvct_keys
                if !(v ∈ names(val))
                    error(LOGGER, "Missing column $v in the PVCDO input dataframe" )
                end
            end
            fluid_opt["PVTO"] = val
        end
    end
    # PVTW
    p = "PVTW"
    check_required_keyword(p, keywords)
    check_keyword_type(p, options, (String, DataFrame, Dict))
    val = options[p]
    if isa(val, String)
        info(LOGGER, "Reading $p table from $val")
        fluid_opt["PVTW"] = read_pvtc_table(val, p)
    elseif isa(val, Dict)
        for v in pvct_keys
            if !(v ∈ keys(val))
                error(LOGGER, "Missing key $v in the PVTW input dictionary" )
            end
        end
        fluid_opt["PVTW"] = val
    else
        for v in pvct_keys
            if !(v ∈ names(val))
                error(LOGGER, "Missing column $v in the PVTW input dataframe" )
            end
        end
        fluid_opt["PVTW"] = val
    end
    swof_keys = ("sw", "krw", "kro", "pcw")
    swof_corey_keys = ("swi", "sor", "aw", "ao", "krw0", "kro0")
    # Relative permeability
    v = ("SWOF", "SWOFCorey")
    p = check_complementary_keywords(v, keywords)
    fluid_opt["SWOF_TYPE"] = p
    val = options[p]
    if p == "SWOF"
        check_keyword_type(p, options, (String, DataFrame))
        if isa(val, String)
            info(LOGGER, "Reading $p table from $val")
            fluid_opt["SWOF"] = read_swof_table(val)
        else # DataFrame
            for v in swof_keys
                if !(v ∈ names(val))
                    error(LOGGER, "Missing column $v in the SWOF input dataframe" )
                end
            end
            fluid_opt["SWOF"] = val
        end
    else
        check_keyword_type(p, options, (DataFrame, Dict))
        if isa(val, DataFrame)
            for v in swof_corey_keys
                if !(v ∈ names(val))
                    error(LOGGER, "Missing column $v in the SWOFCorey input dataframe" )
                end
            end
            fluid_opt["SWOF"] = val
        else # Dict
            for v in swof_corey_keys
                if !(v ∈ keys(val))
                    error(LOGGER, "Missing column $v in the SWOFCorey input dictionary" )
                end
            end
            fluid_opt["SWOF"] = val
        end
    end
    # Initial solution
    p = "sw"
    check_required_keyword(p, keywords)
    fluid_opt["sw"] = parse_vector_keyword(p, options, nc)

    v = ("equil", "po")
    p = check_complementary_keywords(v, keywords)
    fluid_opt["po_type"] = p
    if p == "po"
        fluid_opt[p] = parse_vector_keyword(p, options, nc)
    else
        check_keyword_type(p, options, (Tuple{Float64, Float64},))
        fluid_opt[p] = options[p]
    end
    return fluid_opt
end



function read_swof_table(fn::String)
    df = read_table_file(fn, "SWOF")
    rename!(df, [:sw, :krw, :kro, :pcw])
    # Add flat extrapolation
    insert!.(eachcol(df), 1, [-1.e30, df[1, :krw], df[1, :kro], df[1, :pcw]])
    insert!.(eachcol(df), size(df)[1]+1, [1.e30, df[end,:krw], df[end,:kro], df[end,:pcw]])
    #
    return df
end

function read_pvtc_table(fn::String, keyword::String)
    df = read_table_file(fn, keyword)
    vecs = [:pref, :bref, :c, :μref, :cμ]
    rename!(df, vecs)
    return df
end

function read_pvt_table(fn::String, keyword::String)
    df = read_table_file(fn, keyword)
    rename!(df, [:p, :b, :μ])
    # Add flat extrapolation
    insert!.(eachcol(df), 1, [-1.e30, df[1,:b], df[1,:μ]])
    insert!.(eachcol(df), size(df)[1]+1, [1.e30, df[end,:b], df[end,:μ]])
    return df
end

