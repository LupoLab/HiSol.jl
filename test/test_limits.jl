import Test: @test, @testset
using HISOL
import Luna: Tools
import Roots: find_zero

@testset "barrier suppression intensity" begin
    @test isapprox(HISOL.Limits.barrier_suppression_intensity(:Xe), 8.66320e17, rtol=1e-4)
    @test isapprox(HISOL.Limits.barrier_suppression_intensity(:Kr), 1.53336e18, rtol=1e-4)
    @test isapprox(HISOL.Limits.barrier_suppression_intensity(:Ar), 2.46849e18, rtol=1e-4)
    @test isapprox(HISOL.Limits.barrier_suppression_intensity(:Ne), 8.65198e18, rtol=1e-4)
    @test isapprox(HISOL.Limits.barrier_suppression_intensity(:He), 1.46225e19, rtol=1e-4)
end

##
@testset "maximum soliton order - ionisation" begin
    a = 350e-6
    τ = [5e-15, 10e-15, 30e-15]
    gas = [:Xe, :Kr, :Ar, :Ne, :HeJ]
    λ0 = [400e-9, 800e-9, 1200e-9, 2000e-9]
    λzdfrac = [1/2, 2/3, 3/4]
    @testset "$gi - λ0 $(1e9λ0i) nm - λzd $(1e9λ0i*frac) - $(1e15τi) fs" for gi in gas, λ0i in λ0, frac in λzdfrac, τi in τ
        λzdi = λ0i*frac
        S_ion = 10
        Isupp = HISOL.Limits.barrier_suppression_intensity(gi)
        # find energy where peak intensity reaches Isupp/S_ion purely using functions from Luna
        pressure = Tools.pressureZDW(a, gi, λzdi)
        EmaxLuna = find_zero(500e-6) do E
            p = Tools.capillary_params(E, τi, λ0i, a, gi; P=pressure)
            intensity = p.I0
            p.I0 - Isupp/S_ion
        end
        # get max soliton order from max energy
        p = Tools.capillary_params(EmaxLuna, τi, λ0i, a, gi; P=pressure)
        NmaxLuna = p.N
        # compare to HISOL.jl
        Nmax = HISOL.Limits.Nmax_ion(λzdi, gi, λ0i, τi; S_ion)
        @test isapprox(NmaxLuna, Nmax, rtol=1e-3)
    end
end

##
@testset "maximum soliton order - self-focusing" begin
    a = 350e-6
    τ = [5e-15, 10e-15, 30e-15]
    # gas = [:Xe, :Kr, :Ar, :Ne, :HeJ]
    gas = [:HeJ]
    λ0 = [400e-9, 800e-9, 1200e-9, 2000e-9]
    λzdfrac = [1/2, 2/3, 3/4]
    @testset "$gi - λ0 $(1e9λ0i) nm - λzd $(1e9λ0i*frac) - $(1e15τi) fs" for gi in gas, λ0i in λ0, frac in λzdfrac, τi in τ
        λzdi = λ0i*frac
        S_sf = 5
        # find energy where peak power reaches Pcrit/S_sf using functions purely from Luna
        pressure = Tools.pressureZDW(a, gi, λzdi)
        EmaxLuna = find_zero(500e-6) do E
            p = Tools.capillary_params(E, τi, λ0i, a, gi; P=pressure)
            p.P0 - p.Pcr/S_sf
        end
        # get max soliton order from max energy
        p = Tools.capillary_params(EmaxLuna, τi, λ0i, a, gi; P=pressure)
        NmaxLuna = p.N
        # compare to HISOL.jl
        Nmax = HISOL.Limits.Nmax_sf(λzdi, gi, λ0i, τi; S_sf)
        @test isapprox(NmaxLuna, Nmax, rtol=1e-3)
    end
end