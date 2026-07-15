module Plotting

"""
    getext()

Return the currently loaded plotting backend extension. Tries PythonPlotExt, PyPlotExt,
and MakieExt in order and errors if no backend is loaded.
"""
function getext()
    pkg = parentmodule(@__MODULE__)
    ext = Base.get_extension(pkg, :PythonPlotExt)
    !isnothing(ext) && return ext
    ext = Base.get_extension(pkg, :PyPlotExt)
    !isnothing(ext) && return ext
    ext = Base.get_extension(pkg, :MakieExt)
    !isnothing(ext) && return ext
    error("No plotting backend loaded. Please load one of: PythonPlot, PyPlot, GLMakie, CairoMakie, or WGLMakie.")
end

"""
    getext_makie()

Return the Makie backend extension, erroring if it is not loaded. Used to gate features
which are only available with Makie (the interactive GUIs and 3D plots).
"""
function getext_makie()
    pkg = parentmodule(@__MODULE__)
    ext = Base.get_extension(pkg, :MakieExt)
    !isnothing(ext) && return ext
    error("This feature is only available with a Makie backend. Please load one of: GLMakie (interactive), WGLMakie, or CairoMakie (static plots only).")
end

end
