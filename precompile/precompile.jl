using PackageCompiler

dir = @__DIR__

if !(@isdefined replace_default)
    replace_default = 1
end
if !(@isdefined sysimage_path)
    sysimage_path = "ResSimADSysImg.so"
end

if replace_default == 0
    println("Precompiling ResSimAD.jl.")
    println("Custom system image will be created at $sysimage_path.")
    create_sysimage(:ResSimAD; sysimage_path=sysimage_path,
                    precompile_execution_file=joinpath(dir, "precompile_execution_file.jl"))
else
    println("Precompiling ResSimAD.jl.")
    println("Default system image will be replaced.")
    create_sysimage(:ResSimAD; replace_default=true,
                    precompile_execution_file=joinpath(dir, "precompile_execution_file.jl"))
end
