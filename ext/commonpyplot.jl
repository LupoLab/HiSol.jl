# Common implementation for PyPlotExt and PythonPlotExt.
# Expects the including module to define:
#   plt - the plotting module (PyPlot.plt or PythonPlot.pyplot)
#
# Each render function takes a NamedTuple of prepared plot data from the corresponding
# function in HiSol.Design, HiSol.Compressor or HiSol.Focusing.

import Printf: @sprintf

function cmap_colours(num, cmap="viridis"; cmin=0, cmax=0.8)
    cm = getproperty(plt.cm, Symbol(cmap))
    [cm(v) for v in range(cmin, cmax; length=num)]
end

function design_space_a_energy(d)
    patch = plt.Polygon(1e6d.cropped.vertices; closed=false,
                        facecolor="0.5", edgecolor="none", alpha=0.3)

    fig = plt.figure()
    fig.set_size_inches(12, 3.5)
    plt.subplot(1, 4, 1)
    plt.pcolormesh(1e6d.a, 1e6d.energy, d.ratios.loss; cmap="Spectral_r", rasterized=true)
    plt.clim(0, 2)
    plt.contour(1e6d.a, 1e6d.energy, Float64.(d.idcs.all); levels=[0.5], colors="k")
    plt.axhline(1e6d.full.loss; color="0.4")
    plt.gca().add_patch(patch)
    plt.ylabel("Energy (μJ)")
    plt.xlabel("Core radius (μm)")
    plt.title("Loss")
    plt.ylim(extrema(1e6d.energy))
    plt.xlim(extrema(1e6d.a))
    plt.subplot(1, 4, 2)
    plt.pcolormesh(1e6d.a, 1e6d.energy, d.ratios.fiss; cmap="Spectral_r", rasterized=true)
    plt.clim(0, 2)
    plt.plot(d.full.length*1e6, d.energyb*1e6, color="0.4")
    plt.gca().set_yticklabels([])
    plt.xlabel("Core radius (μm)")
    plt.title("Fission length")
    plt.ylim(extrema(1e6d.energy))
    plt.xlim(extrema(1e6d.a))
    plt.subplot(1, 4, 3)
    plt.pcolormesh(1e6d.a, 1e6d.energy, d.ratios.Nmin; cmap="Spectral_r", rasterized=true)
    plt.clim(0, 2)
    plt.plot(1e6d.ab, 1e6d.full.Nmin, color="0.4")
    plt.xlabel("Core radius (μm)")
    plt.gca().set_yticklabels([])
    plt.title("Minimum soliton order")
    plt.ylim(extrema(1e6d.energy))
    plt.xlim(extrema(1e6d.a))
    plt.subplot(1, 4, 4)
    plt.pcolormesh(1e6d.a, 1e6d.energy, d.ratios.Nmax; cmap="Spectral_r", rasterized=true)
    plt.clim(0, 2)
    plt.plot(1e6d.ab, 1e6d.full.Nmax, color="0.4")
    plt.xlabel("Core radius (μm)")
    plt.gca().set_yticklabels([])
    plt.title("Maximum soliton order")
    plt.ylim(extrema(1e6d.energy))
    plt.xlim(extrema(1e6d.a))

    fig.tight_layout()

    fig2 = plt.figure()
    fig2.set_size_inches(8, 3.5)
    gs = fig2.add_gridspec(1, 2)
    gss = gs[1].subgridspec(1, 2; width_ratios=(1, 0.05), wspace=0.05)
    ax = fig2.add_subplot(gss[1])
    img = ax.pcolormesh(1e6d.a, 1e6d.energy, d.Lfiss_ok, rasterized=true)
    ax.set_xlabel("Core radius (μm)")
    ax.set_ylabel("Energy (μJ)")
    ax.set_title("Fission length")
    cax = fig2.add_subplot(gss[2])
    fig2.colorbar(img; cax, label="Fission length (m)")
    gss = gs[2].subgridspec(1, 2; width_ratios=(1, 0.05), wspace=0.05)
    ax = fig2.add_subplot(gss[1])
    img = ax.pcolormesh(1e6d.a, 1e6d.energy, d.N_ok, rasterized=true)
    ax.set_xlabel("Core radius (μm)")
    ax.set_ylabel("Energy (μJ)")
    ax.set_title("Soliton order")
    cax = fig2.add_subplot(gss[2])
    fig2.colorbar(img; cax, label="Soliton order")
    fig2.tight_layout()

    fig, fig2
end

function aplot_energy_maxlength(d)
    a, ratios, idcs, agood = d.a, d.ratios, d.idcs, d.agood
    cols = cmap_colours(4)
    fig = plt.figure()
    l, = plt.plot(1e6a, ratios.loss; label="\$L_l / 1.5L_f\$", c=cols[1])
    plt.fill_between(1e6a[idcs.loss], ones(count(idcs.loss)), ratios.loss[idcs.loss]; color=l.get_color(), alpha=0.25)
    l, = plt.plot(1e6a, ratios.fiss; label="\$L_{max} / 1.5L_f\$", c=cols[2])
    plt.fill_between(1e6a[idcs.fiss], ones(count(idcs.fiss)), ratios.fiss[idcs.fiss]; color=l.get_color(), alpha=0.25)
    l, = plt.plot(1e6a, ratios.Nmin; label="\$N / N_{min}\$", c=cols[3])
    plt.fill_between(1e6a[idcs.Nmin], ones(count(idcs.Nmin)), ratios.Nmin[idcs.Nmin]; color=l.get_color(), alpha=0.25)
    l, = plt.plot(1e6a, ratios.Nmax; label="\$N_{max} / N\$", c=cols[4])
    plt.fill_between(1e6a[idcs.Nmax], ones(count(idcs.Nmax)), ratios.Nmax[idcs.Nmax]; color=l.get_color(), alpha=0.25)
    plt.axhline(1; color="0.5")
    if length(agood) > 1
        plt.fill_between(1e6*agood, 0, 1, color="r", alpha=0.2)
    end
    plt.legend()
    plt.ylim(0, 1.5*ratios.loss[1])
    plt.xlim(1e6.*extrema(a))
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Ratio")
    fig
end

function plot_optimise(d)
    a = d.a
    fig = plt.figure()
    fig.set_size_inches(15, 4)
    plt.subplot(1, 4, 1)
    plt.plot(1e6a, d.broadfac)
    if ~isnothing(d.opt)
        plt.plot(1e6d.opt.a, d.factor, "o"; color="k", label=@sprintf("%.1f μm", 1e6d.opt.a))
    end
    if ~isnothing(d.dot)
        plt.plot(1e6d.dot.a, d.dot.broadening_factor, "o"; color="b", label=@sprintf("%.1f μm", 1e6d.dot.a))
    end
    plt.axhline(d.factor; linestyle="--", color="k", label="Required")
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.legend(;frameon=false, fontsize=10)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Broadening factor")

    plt.subplot(1, 4, 2)
    plt.plot(1e6a, d.flength; label="HCF")
    plt.plot(1e6a, d.Ltot, "--"; label="Total")
    if ~isnothing(d.opt)
        plt.plot(1e6d.opt.a, d.opt.flength, "o"; color="k", label=@sprintf("%.2f m", d.opt.flength))
        plt.plot(1e6d.opt.a, d.opt.Ltot, "o"; fillstyle="none", color="k", label=@sprintf("%.2f m", d.opt.Ltot))
    end
    if ~isnothing(d.dot)
        plt.plot(1e6d.dot.a, d.dot.flength, "o"; color="b", label=@sprintf("%.2f m", d.dot.flength))
        plt.plot(1e6d.dot.a, d.dot.Ltot, "o"; fillstyle="none", color="b", label=@sprintf("%.2f m", d.dot.Ltot))
    end
    plt.legend(;frameon=false, fontsize=10)
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Length (m)")

    plt.subplot(1, 4, 3)
    plt.plot(1e6a, 100*d.transmission)
    if ~isnothing(d.opt)
        plt.plot(1e6d.opt.a, 100*d.opt.transmission, "o"; color="k", label=@sprintf("%.1f %%", 100*d.opt.transmission))
        plt.legend(;frameon=false, fontsize=10)
    end
    if ~isnothing(d.dot)
        plt.plot(1e6d.dot.a, 100*d.dot.transmission, "o"; color="b", label=@sprintf("%.1f %%", 100*d.dot.transmission))
        plt.legend(;frameon=false, fontsize=10)
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(ymin=0)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Transmission (%)")

    plt.subplot(1, 4, 4)
    plt.plot(1e6a, d.intensity*1e-4)
    plt.axhline(d.Isupp*1e-4/d.S_ion; linestyle="--", color="r", label="Limit")
    if ~isnothing(d.opt)
        plt.plot(1e6d.opt.a, d.opt.intensity*1e-4, "o"; color="k", label=@sprintf("%.2e W/cm\$^{-2}\$", d.opt.intensity*1e-4))
    end
    if ~isnothing(d.dot)
        plt.plot(1e6d.dot.a, d.dot.intensity*1e-4, "o"; color="b", label=@sprintf("%.2e W/cm\$^{-2}\$", d.dot.intensity*1e-4))
    end
    plt.xlim(1e6.*extrema(a))
    plt.ylim(0, 2*d.Isupp*1e-4/d.S_ion)
    plt.legend(;frameon=false, fontsize=10)
    plt.xlabel("Core radius (μm)")
    plt.ylabel("Intensity (W/cm\$^{-2}\$)")

    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.suptitle(@sprintf("Input peak power: %.2f GW | Pressure: %.2f bar", d.P0*1e-9, d.pressure))

    fig
end

function plot_window_thickness_variable(d)
    fig = plt.figure()
    plt.plot(d.distance*1e2, d.tNL*1e3; label="Nonlinear limit")
    plt.plot(d.distance*1e2, d.tP*1e3; label="Pressure limit")
    plt.plot(d.dOpt*1e2, d.tOpt*1e3, "k.";
             label=@sprintf("%.2f mm thickness, %.2f cm away, %.2f mm aperture",
                            d.tOpt*1e3, d.dOpt*1e2, 1e3d.apOpt))
    plt.axvline(d.dOpt*1e2; linestyle="--", color="0.5")
    if ~isnothing(d.LIDT_distance)
        plt.axvline(d.LIDT_distance*1e2; linestyle="--", color="r")
    end
    plt.xlabel("Distance (cm)")
    plt.ylabel("Window thickness (mm)")
    plt.legend()
    rax = plt.gca().twinx()
    rax.plot(d.distance*1e2, d.aperture*1e3, "k--")
    rax.set_ylabel("Aperture radius (mm)")
    fig
end

function plot_window_thickness_fixed(d)
    fig = plt.figure()
    plt.plot(d.distance*1e2, d.tNL*1e3; label="Nonlinear limit")
    plt.axhline(d.tP*1e3; color="C1", label="Pressure limit")
    if ~isnothing(d.crossing)
        plt.plot(d.crossing.distance*1e2, d.crossing.thickness*1e3, "k.";
                 label=@sprintf("%.2f mm thickness, %.2f cm away",
                                d.crossing.thickness*1e3, d.crossing.distance*1e2))
    end
    plt.xlabel("Distance (cm)")
    plt.ylabel("Window thickness (mm)")
    plt.legend()
    fig
end

# Makie-only features
const MAKIEONLY = "is only available with a Makie backend. Please load one of: GLMakie (interactive), WGLMakie, or CairoMakie (static plots only)."

design_space_3D(args...; kwargs...) = error("The 3D design-space plot "*MAKIEONLY)
interactive_design_space(args...; kwargs...) = error("The interactive design-space GUI "*MAKIEONLY)
interactive_design_space_3D(args...; kwargs...) = error("The interactive 3D design-space GUI "*MAKIEONLY)
