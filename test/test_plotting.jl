import Test: @test, @testset, @test_throws
import HiSol
import HiSol: WindowConstraint

wc800 = WindowConstraint(800e-9, :SiO2; thickness=1e-3, LIDT=2000)

@testset "no backend loaded" begin
    # these tests are only meaningful if no plotting backend has been loaded yet
    # in this session
    if isnothing(Base.get_extension(HiSol, :MakieExt)) &&
       isnothing(Base.get_extension(HiSol, :PyPlotExt)) &&
       isnothing(Base.get_extension(HiSol, :PythonPlotExt))
        @test_throws "No plotting backend loaded" HiSol.Plotting.getext()
        @test_throws "only available with a Makie backend" HiSol.Plotting.getext_makie()
        @test_throws "No plotting backend loaded" HiSol.design_space_a_energy(
            160e-9, :He, 800e-9, 10e-15, 5;
            input_constraint=wc800, output_constraint=wc800, Nplot=32)
    end
end

import CairoMakie

@testset "CairoMakie static plots" begin
    figs, params, a, energy, ratios = HiSol.design_space_a_energy(
        160e-9, :He, 800e-9, 10e-15, 5;
        input_constraint=wc800, output_constraint=wc800, Nplot=64)
    @test length(figs) == 2
    @test figs[1] isa CairoMakie.Makie.Figure
    @test figs[2] isa CairoMakie.Makie.Figure
    p = params(125e-6, 200e-6)
    @test isapprox(p.flength, 2.8737700334469096; rtol=1e-6)

    fig, paramsf = HiSol.Design.aplot_energy_maxlength(
        200e-9, :HeJ, 800e-9, 10e-15, 100e-6, 2;
        input_constraint=wc800, output_constraint=wc800, Na=32)
    @test fig isa CairoMakie.Makie.Figure

    wc1030 = WindowConstraint(1030e-9, :MgF2; thickness=1e-3)
    fig, f = HiSol.Compressor.plot_optimise(220e-15, 50e-15, :Kr, 1030e-9, 200e-6, 1.5;
                                            input_constraint=wc1030,
                                            output_constraint=wc1030, Na=64)
    @test fig isa CairoMakie.Makie.Figure
end

@testset "3D design space" begin
    fig, data = HiSol.design_space_a_energy_3D(
        range(140e-9, 220e-9, 6), :He, 800e-9, 10e-15, 5;
        input_constraint=wc800, output_constraint=wc800, Nplot=(24, 24, 6))
    @test fig isa CairoMakie.Makie.Figure
    @test size(data.margin) == (24, 24, 6)
    @test data.thirdaxis == :λ_target
    @test length(data.third) == 6
    @test any(data.margin .< 1) # design space is not empty
    @test any(data.margin .> 1)

    fig, data = HiSol.design_space_a_energy_3D(
        160e-9, :He, 800e-9, range(8e-15, 12e-15, 4), 5;
        input_constraint=wc800, output_constraint=wc800, Nplot=(16, 16, 4))
    @test data.thirdaxis == :τfwhm
    @test size(data.margin) == (16, 16, 4)
    @test any(data.margin .< 1)

    # exactly one of λ_target and τfwhm must be a range
    @test_throws ArgumentError HiSol.design_space_a_energy_3D(
        range(140e-9, 220e-9, 4), :He, 800e-9, range(8e-15, 12e-15, 4), 5;
        input_constraint=wc800, output_constraint=wc800)
    @test_throws ArgumentError HiSol.design_space_a_energy_3D(
        160e-9, :He, 800e-9, 10e-15, 5;
        input_constraint=wc800, output_constraint=wc800)
end

@testset "window thickness plots" begin
    fig = HiSol.Focusing.plot_window_thickness_variable(
        125e-6, 2, 200e-6, 10e-15, 800e-9, 800e-9)
    @test fig isa CairoMakie.Makie.Figure
    fig = HiSol.Focusing.plot_window_thickness_fixed(
        125e-6, 2, 20e9, 800e-9, 800e-9, 10e-3)
    @test fig isa CairoMakie.Makie.Figure
end
