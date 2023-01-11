import Test: @test, @testset
using HiSol

@testset "dispersion" begin
    a = [100e-6, 200e-6, 400e-6]
    gas = [:Xe, :Kr, :Ar, :Ne, :HeJ]
    λ0 = [400e-9, 800e-9, 1200e-9, 2000e-9]
    λzdfrac = [1/2, 2/3, 3/4]
    @testset "$gi - λ0 $(1e9λ0i) nm - λzd $(1e9λ0i*frac) - $(1e6ai) μm" for gi in gas, λ0i in λ0, frac in λzdfrac, ai in a
        λzdi = λ0i*frac
        β2δ = HiSol.HCF.δ(gi, λ0i, λzdi)/ai^2
        pressure = HiSol.HCF.ZDW_pressure(λzdi, ai, gi)
        β2n = HiSol.HCF.dispersion(ai, gi, pressure, λ0i)
        @test isapprox(β2δ, β2n; rtol=1e-5)
    end
end