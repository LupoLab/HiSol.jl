import Test: @test, @testset
using HiSol
import Luna.PhysData: density

@testset "dispersion" begin
    a = [100e-6, 200e-6, 400e-6]
    gas = [:Xe, :Kr, :Ar, :Ne, :HeJ]
    λ0 = [400e-9, 800e-9, 1200e-9, 2000e-9]
    λzdfrac = [1/2, 2/3, 3/4]
    @testset "$gi - λ0 $(1e9λ0i) nm - λzd $(1e9λ0i*frac) - $(1e6ai) μm" for gi in gas, λ0i in λ0, frac in λzdfrac, ai in a
        λzdi = λ0i*frac
        β2δ = HiSol.HCF.δ(gi, λ0i, λzdi)/ai^2
        density = HiSol.HCF.ZDW_density(λzdi, ai, gi)
        β2Δ = HiSol.HCF.Δ(gi, λ0i, density*ai^2)/ai^2
        pressure = HiSol.HCF.ZDW_pressure(λzdi, ai, gi)
        β2n = HiSol.HCF.dispersion(ai, gi, pressure, λ0i)
        @test isapprox(β2δ, β2n; rtol=1e-5)
        @test isapprox(β2Δ, β2n; rtol=1e-5)
    end
end

@testset "nonlinear coeff" begin
    λ0 = 800e-9
    a = (50:50:300) * 1e-6
    gas = [:HeJ, :Ne, :Ar, :Kr]
    λzd = (400:50:800) * 1e-9
    m = 1:4
    @testset "$(ai*1e6) μm - $gi - ZDW $(λi*1e9) nm - $(m)th mode" for ai in a, gi in gas, λi in λzd, m in m
        pressure = HiSol.HCF.ZDW_pressure(λi, ai, gi; m)
        ρ = density(gi, pressure) 
        γ1 = HiSol.HCF.γ(ai, gi, pressure, λ0; m)
        γ2 = HiSol.HCF.γ(ai, gi, λ0; λzd=λi, m)
        γ3 = HiSol.HCF.γ(ai, gi, λ0; ρasq=ρ*ai^2, m)
        @test isapprox(γ1, γ2; rtol=1e-6)
    end
end