module MakieExt
import Makie
import Printf: @sprintf
import HiSol
import HiSol: Design
import HiSol.Focusing: WindowConstraint, DamageConstraint, FixedConstraint, NoConstraint

const RATIO_CMAP = Makie.Reverse(:Spectral)

function newfig(; size=(800, 600))
    fig = Makie.Figure(; size)
    if nameof(Makie.current_backend()) in (:GLMakie, :WGLMakie)
        Makie.DataInspector(fig)
    end
    fig
end

function display_fig(fig)
    try
        backend = Makie.current_backend()
        if nameof(backend) == :GLMakie && isdefined(backend, :Screen)
            display(backend.Screen(), fig)
        else
            display(fig)
        end
    catch
        # no display available (e.g. headless CairoMakie) -- the figure is still returned
    end
    fig
end

function cmap_colours(num, cmap=:viridis; cmin=0, cmax=0.8)
    cg = Makie.cgrad(cmap)
    [cg[v] for v in range(cmin, cmax; length=num)]
end

#===================================================================================
Static 2D plots
===================================================================================#

# draw one criteria-ratio panel: heatmap of the ratio plus overlays
function ratio_panel!(ax, a, energy, ratio)
    # rasterize: embed the heatmap as a bitmap in vector (SVG/PDF) output
    hm = Makie.heatmap!(ax, 1e6a, 1e6energy, permutedims(ratio);
                        colormap=RATIO_CMAP, colorrange=(0, 2), rasterize=4)
    Makie.xlims!(ax, extrema(1e6a)...)
    Makie.ylims!(ax, extrema(1e6energy)...)
    hm
end

function design_space_panels!(gl, d)
    axs = [Makie.Axis(gl[1, ii]; xlabel="Core radius (μm)", title) for
           (ii, title) in enumerate(("Loss", "Fission length",
                                     "Minimum soliton order", "Maximum soliton order"))]
    axs[1].ylabel = "Energy (μJ)"
    for ax in axs[2:end]
        ax.yticklabelsvisible = false
    end
    Makie.linkaxes!(axs...)

    hm = ratio_panel!(axs[1], d.a, d.energy, d.ratios.loss)
    Makie.contour!(axs[1], 1e6d.a, 1e6d.energy, permutedims(Float64.(d.idcs.all));
                   levels=[0.5], color=:black)
    Makie.hlines!(axs[1], 1e6d.full.loss; color=(:black, 0.4))
    if size(d.cropped.vertices, 1) > 2
        Makie.poly!(axs[1], Makie.Point2f.(1e6d.cropped.vertices[:, 1],
                                           1e6d.cropped.vertices[:, 2]);
                    color=(:gray, 0.3), strokewidth=0)
    end

    ratio_panel!(axs[2], d.a, d.energy, d.ratios.fiss)
    Makie.lines!(axs[2], 1e6d.full.length, 1e6d.energyb; color=(:black, 0.4))

    ratio_panel!(axs[3], d.a, d.energy, d.ratios.Nmin)
    Makie.lines!(axs[3], 1e6d.ab, 1e6d.full.Nmin; color=(:black, 0.4))

    ratio_panel!(axs[4], d.a, d.energy, d.ratios.Nmax)
    Makie.lines!(axs[4], 1e6d.ab, 1e6d.full.Nmax; color=(:black, 0.4))

    Makie.Colorbar(gl[1, 5], hm; label="Criterion ratio")
    axs
end

function design_space_a_energy(d)
    fig = newfig(size=(1200, 350))
    design_space_panels!(fig.layout, d)

    fig2 = newfig(size=(800, 350))
    ax = Makie.Axis(fig2[1, 1]; xlabel="Core radius (μm)", ylabel="Energy (μJ)",
                    title="Fission length")
    hm = Makie.heatmap!(ax, 1e6d.a, 1e6d.energy, permutedims(d.Lfiss_ok); rasterize=4)
    Makie.Colorbar(fig2[1, 2], hm; label="Fission length (m)")
    ax2 = Makie.Axis(fig2[1, 3]; xlabel="Core radius (μm)", ylabel="Energy (μJ)",
                     title="Soliton order")
    hm2 = Makie.heatmap!(ax2, 1e6d.a, 1e6d.energy, permutedims(d.N_ok); rasterize=4)
    Makie.Colorbar(fig2[1, 4], hm2; label="Soliton order")

    display_fig(fig)
    display_fig(fig2)
    fig, fig2
end

function aplot_energy_maxlength(d)
    a, ratios, idcs, agood = d.a, d.ratios, d.idcs, d.agood
    cols = cmap_colours(4)
    labels = ("L_loss / 1.5 L_fiss", "L_max / 1.5 L_fiss", "N / N_min", "N_max / N")
    fig = newfig(size=(700, 450))
    ax = Makie.Axis(fig[1, 1]; xlabel="Core radius (μm)", ylabel="Ratio")
    for (ratio, idx, col, label) in zip(ratios, idcs, cols, labels)
        Makie.lines!(ax, 1e6a, ratio; color=col, label)
        if count(idx) > 1
            Makie.band!(ax, 1e6a[idx], ones(count(idx)), ratio[idx];
                        color=(col, 0.25))
        end
    end
    Makie.hlines!(ax, 1; color=(:black, 0.5))
    if length(agood) > 1
        Makie.band!(ax, 1e6agood, zeros(length(agood)), ones(length(agood));
                    color=(:red, 0.2))
    end
    Makie.axislegend(ax; framevisible=false)
    Makie.ylims!(ax, 0, 1.5*ratios.loss[1])
    Makie.xlims!(ax, 1e6 .* extrema(a)...)
    display_fig(fig)
    fig
end

function plot_optimise(d)
    a = d.a
    fig = newfig(size=(1500, 400))
    Makie.Label(fig[0, 1:4],
                @sprintf("Input peak power: %.2f GW | Pressure: %.2f bar", d.P0*1e-9, d.pressure);
                fontsize=16)

    ax1 = Makie.Axis(fig[1, 1]; xlabel="Core radius (μm)", ylabel="Broadening factor")
    Makie.lines!(ax1, 1e6a, d.broadfac)
    Makie.hlines!(ax1, d.factor; linestyle=:dash, color=:black, label="Required")
    if ~isnothing(d.opt)
        Makie.scatter!(ax1, 1e6d.opt.a, d.factor; color=:black,
                       label=@sprintf("%.1f μm", 1e6d.opt.a))
    end
    if ~isnothing(d.dot)
        Makie.scatter!(ax1, 1e6d.dot.a, d.dot.broadening_factor; color=:blue,
                       label=@sprintf("%.1f μm", 1e6d.dot.a))
    end
    Makie.axislegend(ax1; framevisible=false)
    Makie.ylims!(ax1; low=0)

    ax2 = Makie.Axis(fig[1, 2]; xlabel="Core radius (μm)", ylabel="Length (m)")
    Makie.lines!(ax2, 1e6a, d.flength; label="HCF")
    Makie.lines!(ax2, 1e6a, d.Ltot; linestyle=:dash, label="Total")
    if ~isnothing(d.opt)
        Makie.scatter!(ax2, 1e6d.opt.a, d.opt.flength; color=:black,
                       label=@sprintf("%.2f m", d.opt.flength))
        Makie.scatter!(ax2, 1e6d.opt.a, d.opt.Ltot; color=:black, marker=:circle,
                       strokecolor=:black, strokewidth=1, alpha=0.5,
                       label=@sprintf("%.2f m", d.opt.Ltot))
    end
    if ~isnothing(d.dot)
        Makie.scatter!(ax2, 1e6d.dot.a, d.dot.flength; color=:blue,
                       label=@sprintf("%.2f m", d.dot.flength))
        Makie.scatter!(ax2, 1e6d.dot.a, d.dot.Ltot; color=:blue, marker=:circle,
                       strokecolor=:blue, strokewidth=1, alpha=0.5,
                       label=@sprintf("%.2f m", d.dot.Ltot))
    end
    Makie.axislegend(ax2; framevisible=false)
    Makie.ylims!(ax2; low=0)

    ax3 = Makie.Axis(fig[1, 3]; xlabel="Core radius (μm)", ylabel="Transmission (%)")
    Makie.lines!(ax3, 1e6a, 100*d.transmission)
    if ~isnothing(d.opt)
        Makie.scatter!(ax3, 1e6d.opt.a, 100*d.opt.transmission; color=:black,
                       label=@sprintf("%.1f %%", 100*d.opt.transmission))
    end
    if ~isnothing(d.dot)
        Makie.scatter!(ax3, 1e6d.dot.a, 100*d.dot.transmission; color=:blue,
                       label=@sprintf("%.1f %%", 100*d.dot.transmission))
    end
    (~isnothing(d.opt) || ~isnothing(d.dot)) && Makie.axislegend(ax3; framevisible=false)
    Makie.ylims!(ax3; low=0)

    ax4 = Makie.Axis(fig[1, 4]; xlabel="Core radius (μm)", ylabel="Intensity (W/cm²)")
    Makie.lines!(ax4, 1e6a, d.intensity*1e-4)
    Makie.hlines!(ax4, d.Isupp*1e-4/d.S_ion; linestyle=:dash, color=:red, label="Limit")
    if ~isnothing(d.opt)
        Makie.scatter!(ax4, 1e6d.opt.a, d.opt.intensity*1e-4; color=:black,
                       label=@sprintf("%.2e W/cm²", d.opt.intensity*1e-4))
    end
    if ~isnothing(d.dot)
        Makie.scatter!(ax4, 1e6d.dot.a, d.dot.intensity*1e-4; color=:blue,
                       label=@sprintf("%.2e W/cm²", d.dot.intensity*1e-4))
    end
    Makie.axislegend(ax4; framevisible=false)
    Makie.ylims!(ax4, 0, 2*d.Isupp*1e-4/d.S_ion)

    for ax in (ax1, ax2, ax3, ax4)
        Makie.xlims!(ax, 1e6 .* extrema(a)...)
    end

    display_fig(fig)
    fig
end

function plot_window_thickness_variable(d)
    fig = newfig(size=(700, 450))
    ax = Makie.Axis(fig[1, 1]; xlabel="Distance (cm)", ylabel="Window thickness (mm)")
    Makie.lines!(ax, d.distance*1e2, d.tNL*1e3; label="Nonlinear limit")
    Makie.lines!(ax, d.distance*1e2, d.tP*1e3; label="Pressure limit")
    Makie.scatter!(ax, d.dOpt*1e2, d.tOpt*1e3; color=:black,
                   label=@sprintf("%.2f mm thickness, %.2f cm away, %.2f mm aperture",
                                  d.tOpt*1e3, d.dOpt*1e2, 1e3d.apOpt))
    Makie.vlines!(ax, d.dOpt*1e2; linestyle=:dash, color=(:black, 0.5))
    if ~isnothing(d.LIDT_distance)
        Makie.vlines!(ax, d.LIDT_distance*1e2; linestyle=:dash, color=:red)
    end
    Makie.axislegend(ax; framevisible=false)

    rax = Makie.Axis(fig[1, 1]; yaxisposition=:right, ylabel="Aperture radius (mm)")
    Makie.hidespines!(rax)
    Makie.hidexdecorations!(rax)
    Makie.lines!(rax, d.distance*1e2, d.aperture*1e3; color=:black, linestyle=:dash)
    Makie.linkxaxes!(ax, rax)

    display_fig(fig)
    fig
end

function plot_window_thickness_fixed(d)
    fig = newfig(size=(700, 450))
    ax = Makie.Axis(fig[1, 1]; xlabel="Distance (cm)", ylabel="Window thickness (mm)")
    Makie.lines!(ax, d.distance*1e2, d.tNL*1e3; label="Nonlinear limit")
    Makie.hlines!(ax, d.tP*1e3; color=:orange, label="Pressure limit")
    if ~isnothing(d.crossing)
        Makie.scatter!(ax, d.crossing.distance*1e2, d.crossing.thickness*1e3; color=:black,
                       label=@sprintf("%.2f mm thickness, %.2f cm away",
                                      d.crossing.thickness*1e3, d.crossing.distance*1e2))
    end
    Makie.axislegend(ax; framevisible=false)
    display_fig(fig)
    fig
end

#===================================================================================
3D design space
===================================================================================#

third_scale_label(thirdaxis) = thirdaxis == :λ_target ?
    (1e9, "RDW wavelength (nm)") : (1e15, "Pulse duration (fs)")

# margin can contain Inf (e.g. when the HCF does not fit at all); clamp for rendering
clean_margin(margin) = clamp.(replace(margin, NaN => 10.0, Inf => 10.0), 0.0, 10.0)

function isosurface!(ax, d)
    zscale, _ = third_scale_label(d.thirdaxis)
    vol = clean_margin(d.margin)
    if nameof(Makie.current_backend()) in (:GLMakie, :WGLMakie)
        # 3D contour (VolumeLike) takes interval endpoints, not vectors, for the axes
        Makie.contour!(ax, extrema(1e6d.a), extrema(1e6d.energy), extrema(zscale*d.third),
                       vol; levels=[1.0], alpha=0.5, colormap=[:dodgerblue, :dodgerblue])
    else
        # CairoMakie cannot render volume plots: show the design-space region as a
        # point cloud instead
        idcs = findall(<(1), vol)
        step = max(1, length(idcs) ÷ 20_000)
        idcs = idcs[1:step:end]
        pts = [Makie.Point3f(1e6d.a[I[1]], 1e6d.energy[I[2]], zscale*d.third[I[3]])
               for I in idcs]
        Makie.scatter!(ax, pts; color=(:dodgerblue, 0.15), markersize=4)
    end
end

function make_axis3(gridpos, thirdaxis)
    _, zlabel = third_scale_label(thirdaxis)
    Makie.Axis3(gridpos; xlabel="Core radius (μm)", ylabel="Energy (μJ)", zlabel,
                perspectiveness=0.2, protrusions=(60, 60, 30, 30))
end

function design_space_3D(d)
    fig = newfig(size=(800, 700))
    ax = make_axis3(fig[1, 1], d.thirdaxis)
    isosurface!(ax, d)
    display_fig(fig)
    fig
end

#===================================================================================
GUI helpers
===================================================================================#

const GASES = ["He", "Ne", "Ar", "Kr", "Xe"]
const MATERIALS = ["SiO2", "MgF2"]

# extract seed values for the constraint controls from a constraint object
function seedvals(c::WindowConstraint)
    (type="Window", thickness=isnothing(c.thickness) ? 0.0 : 1e3*c.thickness,
     material=c.n2 isa Symbol ? string(c.n2) : "SiO2", Bmax=c.Bmax,
     LIDT=isnothing(c.LIDT) ? 0.0 : float(c.LIDT), S_fluence=c.S_fluence, distance=0.5)
end
function seedvals(c::DamageConstraint)
    (type="Mirror", thickness=0.0, material="SiO2", Bmax=0.2,
     LIDT=c.LIDT, S_fluence=c.S_fluence, distance=0.5)
end
function seedvals(c::FixedConstraint)
    (type="Fixed", thickness=0.0, material="SiO2", Bmax=0.2,
     LIDT=0.0, S_fluence=5.0, distance=c.distance)
end
function seedvals(c::NoConstraint)
    (type="None", thickness=0.0, material="SiO2", Bmax=0.2,
     LIDT=0.0, S_fluence=5.0, distance=0.5)
end
seedvals(c) = error("Constraints of type $(typeof(c)) cannot be used in the GUI. "*
                    "Use WindowConstraint, DamageConstraint, FixedConstraint or NoConstraint.")

# returns a getter for the current value and the observable to watch for changes
function numbox!(gl, row, label, startvalue)
    Makie.Label(gl[row, 1], label; halign=:left, fontsize=12)
    tb = Makie.Textbox(gl[row, 2]; stored_string=string(startvalue),
                       validator=Float64, width=70, height=22, fontsize=12)
    (() -> parse(Float64, tb.stored_string[])), tb.stored_string
end

fmtval(v) = v isa Integer ? string(v) : string(round(v; digits=2))

# label + full-width slider + editable value textbox in columns 1-3 of gl.
# Moving the slider updates the textbox; committing a value in the textbox
# (enter or defocus) moves the slider.
function labelled_slider!(gl, row, label, range, startvalue)
    Makie.Label(gl[row, 1], label; halign=:left, fontsize=12)
    sl = Makie.Slider(gl[row, 2]; range, startvalue)
    tb = Makie.Textbox(gl[row, 3]; stored_string=fmtval(sl.value[]),
                       validator=Float64, width=60, height=22, fontsize=12)
    Makie.on(sl.value) do v
        tb.displayed_string[] = fmtval(v)
    end
    Makie.on(tb.stored_string) do s
        v = tryparse(Float64, s)
        if !isnothing(v)
            Makie.set_close_to!(sl, v)
            # the typed value snaps to the slider grid; always show the snapped value
            tb.displayed_string[] = fmtval(sl.value[])
        end
    end
    sl
end

# build the control panel for one length constraint; returns a closure which constructs
# the constraint from the current control values plus a vector of observables which fire
# when any of the controls change
function constraint_controls!(gl, title, seed)
    v = seedvals(seed)
    Makie.Label(gl[1, 1:2], title; font=:bold, halign=:left)
    Makie.Label(gl[2, 1], "Type"; halign=:left, fontsize=12)
    typemenu = Makie.Menu(gl[2, 2]; options=["Window", "Mirror", "Fixed", "None"],
                          default=v.type, fontsize=12)
    thickness, thickness_obs = numbox!(gl, 3, "Thickness (mm, 0=auto)", v.thickness)
    Makie.Label(gl[4, 1], "Material"; halign=:left, fontsize=12)
    matmenu = Makie.Menu(gl[4, 2]; options=MATERIALS, default=v.material, fontsize=12)
    Bmax, Bmax_obs = numbox!(gl, 5, "Bmax", v.Bmax)
    LIDT, LIDT_obs = numbox!(gl, 6, "LIDT (J/m², 0=off)", v.LIDT)
    S_fluence, S_fluence_obs = numbox!(gl, 7, "S_fluence", v.S_fluence)
    distance, distance_obs = numbox!(gl, 8, "Distance (m)", v.distance)
    Makie.rowgap!(gl, 4)

    watch = Any[typemenu.selection, matmenu.selection, thickness_obs, Bmax_obs,
                LIDT_obs, S_fluence_obs, distance_obs]

    function makecon(λ0)
        type = typemenu.selection[]
        if type == "Window"
            t = thickness()
            l = LIDT()
            WindowConstraint(λ0, Symbol(matmenu.selection[]);
                             thickness=(t == 0 ? nothing : 1e-3*t),
                             Bmax=Bmax(),
                             LIDT=(l == 0 ? nothing : l),
                             S_fluence=S_fluence())
        elseif type == "Mirror"
            DamageConstraint(λ0, LIDT(); S_fluence=S_fluence())
        elseif type == "Fixed"
            FixedConstraint(distance())
        else
            NoConstraint()
        end
    end

    makecon, watch
end

shorterror(e) = first(sprint(showerror, e), 200)

const LOGO = Base.RefValue{Any}(nothing)

function logo_image()
    if isnothing(LOGO[])
        LOGO[] = Makie.FileIO.load(joinpath(dirname(@__DIR__), "assets", "logo.png"))
    end
    LOGO[]
end

# place the LUPO logo in a decoration-free axis at the given grid position
function logo!(gl, pos)
    img = logo_image()
    h = 40
    w = round(Int, h*size(img, 2)/size(img, 1))
    ax = Makie.Axis(gl[pos...]; aspect=Makie.DataAspect(), height=Makie.Fixed(h),
                    width=Makie.Fixed(w), halign=:left, tellwidth=false)
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    Makie.image!(ax, rotr90(img))
    ax
end

const STALE_TEXT = "Parameters changed — press Compute to update"

# make status show a warning whenever any of the watched controls change
function watchstale!(status, watch)
    for obs in watch
        Makie.on(obs) do _
            status.text = STALE_TEXT
            status.color = :chocolate
        end
    end
end

#===================================================================================
Interactive 2D design space
===================================================================================#

function interactive_design_space(; λ_target, gas, λ0, τfwhm, maxlength,
                                    S_sf, S_ion, S_fiss, S_loss,
                                    input_constraint, output_constraint, Nplot)
    fig = newfig(size=(1450, 850))
    controls = fig[1, 1] = Makie.GridLayout(; valign=:top)
    plots = fig[1, 2] = Makie.GridLayout()
    Makie.colsize!(fig.layout, 1, Makie.Fixed(460))

    # ---- controls: label | slider | value in columns 1-3 ----
    logo!(controls, (1, 1:3))
    Makie.Label(controls[2, 1], "Gas"; halign=:left, fontsize=12)
    gasmenu = Makie.Menu(controls[2, 2:3]; options=GASES, default=string(gas))

    sliders = [
        labelled_slider!(controls, 3, "λ_target (nm)", 60:1:1000, round(1e9λ_target)),
        labelled_slider!(controls, 4, "λ0 (nm)", 200:5:2500, round(1e9λ0)),
        labelled_slider!(controls, 5, "τfwhm (fs)", 3:0.5:500, 1e15τfwhm),
        labelled_slider!(controls, 6, "Length (m)", 0.5:0.1:20, maxlength),
        labelled_slider!(controls, 7, "S_sf", 1:0.5:20, S_sf),
        labelled_slider!(controls, 8, "S_ion", 1:0.5:40, S_ion),
        labelled_slider!(controls, 9, "S_fiss", 1:0.1:5, S_fiss),
        labelled_slider!(controls, 10, "S_loss", 0.5:0.1:5, S_loss),
    ]
    # Auto(false) makes the slider column take up all leftover width
    Makie.colsize!(controls, 2, Makie.Auto(false))
    Makie.colsize!(controls, 3, Makie.Fixed(60))
    Makie.rowgap!(controls, 4)

    congl = controls[11, 1:3] = Makie.GridLayout()
    incon_gl = congl[1, 1] = Makie.GridLayout(; valign=:top)
    makecon_in, watch_in = constraint_controls!(incon_gl, "Input constraint", input_constraint)
    outcon_gl = congl[1, 2] = Makie.GridLayout(; valign=:top)
    makecon_out, watch_out = constraint_controls!(outcon_gl, "Output constraint", output_constraint)

    Makie.Label(controls[12, 1], "Grid points"; halign=:left, fontsize=12)
    nplot_options = sort(unique([128, 256, 512, Nplot]))
    # an integer `default` is interpreted as an index into the options
    nplotmenu = Makie.Menu(controls[12, 2:3]; options=nplot_options,
                           default=findfirst(==(Nplot), nplot_options))
    Makie.Label(controls[13, 1], "Boundary lines"; halign=:left, fontsize=12)
    boundstoggle = Makie.Toggle(controls[13, 2]; active=true, halign=:left)

    computebutton = Makie.Button(controls[14, 1:3]; label="Compute", fontsize=16)
    status = Makie.Label(controls[15, 1:3], "Ready"; halign=:left, fontsize=12)
    readout = Makie.Label(controls[16, 1:3], "Click in a map to inspect a design point";
                          halign=:left, fontsize=12, justification=:left)

    watchstale!(status, vcat(Any[gasmenu.selection, nplotmenu.selection,
                                 boundstoggle.active],
                             Any[sl.value for sl in sliders],
                             watch_in, watch_out))

    # ---- plot axes: 2x2 criteria maps | ratio colorbar | Lfiss and N maps + colorbars ----
    titles = ("Loss", "Fission length", "Minimum soliton order", "Maximum soliton order")
    axs = [Makie.Axis(plots[divrem(ii-1, 2) .+ (1, 1)...];
                      xlabel="Core radius (μm)", ylabel="Energy (μJ)", title=titles[ii])
           for ii in 1:4]
    Makie.Colorbar(plots[1:2, 3]; colormap=RATIO_CMAP, colorrange=(0, 2),
                   label="Criterion ratio")
    axLf = Makie.Axis(plots[1, 4]; xlabel="Core radius (μm)", ylabel="Energy (μJ)",
                      title="Fission length (m)")
    axN = Makie.Axis(plots[2, 4]; xlabel="Core radius (μm)", ylabel="Energy (μJ)",
                     title="Soliton order")
    Makie.linkaxes!(axs..., axLf, axN)
    Lflims = Makie.Observable((0.0, 1.0))
    Nlims = Makie.Observable((0.0, 1.0))
    Makie.Colorbar(plots[1, 5]; colormap=:viridis, limits=Lflims, label="Fission length (m)")
    Makie.Colorbar(plots[2, 5]; colormap=:viridis, limits=Nlims, label="Soliton order")

    paramsf = Base.RefValue{Any}(nothing) # (a, energy) -> spec at the last-computed settings

    function recompute!()
        status.text = "Computing..."
        status.color = :black
        try
            gasnow = Symbol(gasmenu.selection[])
            λt = sliders[1].value[]*1e-9
            λ0now = sliders[2].value[]*1e-9
            τnow = sliders[3].value[]*1e-15
            mlnow = sliders[4].value[]
            S_sfnow, S_ionnow, S_fissnow, S_lossnow = (s.value[] for s in sliders[5:8])
            incon = makecon_in(λ0now)
            outcon = makecon_out(λ0now)
            Nnow = nplotmenu.selection[]

            a, energy, ratios, paramst, idcs, params = Design.maxlength_limitratios(
                λt, gasnow, λ0now, τnow, mlnow;
                input_constraint=incon, output_constraint=outcon,
                S_sf=S_sfnow, S_ion=S_ionnow, S_fiss=S_fissnow, S_loss=S_lossnow,
                Nplot=Nnow)

            bdata = if boundstoggle.active[]
                try
                    ab, energyb, full, cropped = Design.boundaries(
                        λt, gasnow, λ0now, τnow, mlnow;
                        input_constraint=incon, output_constraint=outcon,
                        S_sf=S_sfnow, S_ion=S_ionnow, S_fiss=S_fissnow,
                        Nplot=Nnow)
                    (;ab, energyb, full, cropped)
                catch
                    nothing
                end
            else
                nothing
            end

            aμ, eμ = 1e6a, 1e6energy
            for (ax, ratio) in zip(axs, ratios)
                Makie.empty!(ax)
                Makie.heatmap!(ax, aμ, eμ, permutedims(ratio);
                               colormap=RATIO_CMAP, colorrange=(0, 2))
            end
            Makie.contour!(axs[1], aμ, eμ, permutedims(Float64.(idcs.all));
                           levels=[0.5], color=:black)
            if ~isnothing(bdata)
                Makie.hlines!(axs[1], 1e6bdata.full.loss; color=(:black, 0.4))
                if size(bdata.cropped.vertices, 1) > 2
                    Makie.poly!(axs[1], Makie.Point2f.(1e6bdata.cropped.vertices[:, 1],
                                                       1e6bdata.cropped.vertices[:, 2]);
                                color=(:gray, 0.3), strokewidth=0)
                end
                Makie.lines!(axs[2], 1e6bdata.full.length, 1e6bdata.energyb; color=(:black, 0.4))
                Makie.lines!(axs[3], 1e6bdata.ab, 1e6bdata.full.Nmin; color=(:black, 0.4))
                Makie.lines!(axs[4], 1e6bdata.ab, 1e6bdata.full.Nmax; color=(:black, 0.4))
            end

            Lfiss_ok = copy(paramst.Lfiss)
            Lfiss_ok[.~idcs.all] .= NaN
            N_ok = copy(paramst.N)
            N_ok[.~idcs.all] .= NaN
            anygood = any(idcs.all)
            Lflims[] = anygood ?
                (0.0, maximum(x for x in Lfiss_ok if ~isnan(x))) : (0.0, 1.0)
            Nlims[] = anygood ?
                extrema(x for x in N_ok if ~isnan(x)) : (0.0, 1.0)
            Makie.empty!(axLf)
            Makie.heatmap!(axLf, aμ, eμ, permutedims(Lfiss_ok);
                           colormap=:viridis, colorrange=Lflims[])
            Makie.empty!(axN)
            Makie.heatmap!(axN, aμ, eμ, permutedims(N_ok);
                           colormap=:viridis, colorrange=Nlims[])

            Makie.xlims!(axs[1], extrema(aμ)...)
            Makie.ylims!(axs[1], extrema(eμ)...)

            paramsf[] = (ai, ei) -> params(ai, ei, τnow)
            status.text = anygood ? "Done" : "Done (design space is empty)"
            status.color = :black
        catch e
            status.text = "Error: "*shorterror(e)
            status.color = :red
        end
        nothing
    end

    Makie.on(computebutton.clicks) do _
        recompute!()
    end

    # click in any map to inspect the full system spec at that point
    Makie.on(Makie.events(fig).mousebutton) do ev
        if ev.button == Makie.Mouse.left && ev.action == Makie.Mouse.press
            isnothing(paramsf[]) && return Makie.Consume(false)
            for ax in (axs..., axLf, axN)
                if Makie.is_mouseinside(ax)
                    pos = Makie.mouseposition(ax)
                    try
                        p = paramsf[](pos[1]*1e-6, pos[2]*1e-6)
                        readout.text = @sprintf(
                            "a = %.1f μm, %.1f μJ:\npressure = %.3f bar\nN = %.2f (%.2f–%.2f)\nL_fiss = %.3f m\nL_loss = %.3f m\nmax. HCF length = %.3f m",
                            pos[1], pos[2], p.pressure, p.N, p.Nmin, p.Nmax,
                            p.Lfiss, p.Lloss, p.flength)
                    catch e
                        readout.text = "Error: "*shorterror(e)
                    end
                    break
                end
            end
        end
        return Makie.Consume(false)
    end

    recompute!()
    display_fig(fig)
    fig
end

#===================================================================================
Interactive 3D design space
===================================================================================#

function interactive_design_space_3D(; λ_target, gas, λ0, τfwhm, maxlength,
                                       S_sf, S_ion, S_fiss, S_loss,
                                       input_constraint, output_constraint, Nplot)
    sweep_λt = λ_target isa AbstractVector
    thirdaxis = sweep_λt ? :λ_target : :τfwhm
    zscale, zlabel = third_scale_label(thirdaxis)

    fig = newfig(size=(1450, 850))
    controls = fig[1, 1] = Makie.GridLayout(; valign=:top)
    plots = fig[1, 2] = Makie.GridLayout()
    Makie.colsize!(fig.layout, 1, Makie.Fixed(460))

    # ---- controls (the swept parameter has no slider) ----
    logo!(controls, (1, 1:3))
    Makie.Label(controls[2, 1], "Gas"; halign=:left, fontsize=12)
    gasmenu = Makie.Menu(controls[2, 2:3]; options=GASES, default=string(gas))

    sliders = [labelled_slider!(controls, 3, "λ0 (nm)", 200:5:2500, round(1e9λ0))]
    if sweep_λt
        push!(sliders, labelled_slider!(controls, 4, "τfwhm (fs)", 3:0.5:500, 1e15τfwhm))
    else
        push!(sliders, labelled_slider!(controls, 4, "λ_target (nm)", 60:1:1000,
                                        round(1e9λ_target)))
    end
    push!(sliders, labelled_slider!(controls, 5, "Length (m)", 0.5:0.1:20, float(maxlength)))
    push!(sliders, labelled_slider!(controls, 6, "S_sf", 1:0.5:20, float(S_sf)))
    push!(sliders, labelled_slider!(controls, 7, "S_ion", 1:0.5:40, float(S_ion)))
    push!(sliders, labelled_slider!(controls, 8, "S_fiss", 1:0.1:5, float(S_fiss)))
    push!(sliders, labelled_slider!(controls, 9, "S_loss", 0.5:0.1:5, float(S_loss)))
    # Auto(false) makes the slider column take up all leftover width
    Makie.colsize!(controls, 2, Makie.Auto(false))
    Makie.colsize!(controls, 3, Makie.Fixed(60))
    Makie.rowgap!(controls, 4)

    congl = controls[10, 1:3] = Makie.GridLayout()
    incon_gl = congl[1, 1] = Makie.GridLayout(; valign=:top)
    makecon_in, watch_in = constraint_controls!(incon_gl, "Input constraint", input_constraint)
    outcon_gl = congl[1, 2] = Makie.GridLayout(; valign=:top)
    makecon_out, watch_out = constraint_controls!(outcon_gl, "Output constraint", output_constraint)

    computebutton = Makie.Button(controls[11, 1:3]; label="Compute", fontsize=16)
    status = Makie.Label(controls[12, 1:3], "Ready"; halign=:left, fontsize=12)

    # the slice slider is not watched: moving it does not invalidate the computation
    watchstale!(status, vcat(Any[gasmenu.selection],
                             Any[sl.value for sl in sliders],
                             watch_in, watch_out))

    # ---- plots: rotatable isosurface + 2D slice along the third axis ----
    ax3 = make_axis3(plots[1, 1], thirdaxis)
    axsl = Makie.Axis(plots[1, 2]; xlabel="Core radius (μm)", ylabel="Energy (μJ)")
    Makie.Colorbar(plots[1, 3]; colormap=RATIO_CMAP, colorrange=(0, 2),
                   label="Worst criterion ratio")
    slicegl = plots[2, 1:3] = Makie.GridLayout()
    Makie.Label(slicegl[1, 1], "Slice"; halign=:right)
    sliceslider = Makie.Slider(slicegl[1, 2]; range=0:0.001:1, startvalue=0.5)
    slicelabel = Makie.Label(slicegl[1, 3], ""; halign=:left)

    data = Base.RefValue{Any}(nothing)

    sliceindex() = isnothing(data[]) ? 1 :
        clamp(round(Int, sliceslider.value[]*(length(data[].third)-1)) + 1, 1,
              length(data[].third))

    function updateslice!()
        isnothing(data[]) && return
        d = data[]
        k = sliceindex()
        Makie.empty!(axsl)
        Makie.heatmap!(axsl, 1e6d.a, 1e6d.energy, clean_margin(d.margin[:, :, k]);
                       colormap=RATIO_CMAP, colorrange=(0, 2))
        Makie.contour!(axsl, 1e6d.a, 1e6d.energy, Float64.(d.margin[:, :, k] .< 1);
                       levels=[0.5], color=:black)
        slicelabel.text = @sprintf("%s = %.1f", zlabel, zscale*d.third[k])
        axsl.title = @sprintf("%s = %.1f", zlabel, zscale*d.third[k])
        nothing
    end

    function recompute!()
        status.text = "Computing..."
        status.color = :black
        try
            gasnow = Symbol(gasmenu.selection[])
            λ0now = sliders[1].value[]*1e-9
            if sweep_λt
                λt = λ_target
                τnow = sliders[2].value[]*1e-15
            else
                λt = sliders[2].value[]*1e-9
                τnow = τfwhm
            end
            mlnow = sliders[3].value[]
            S_sfnow, S_ionnow, S_fissnow, S_lossnow = (s.value[] for s in sliders[4:7])
            incon = makecon_in(λ0now)
            outcon = makecon_out(λ0now)

            d = Design.design_space_3D_data(λt, gasnow, λ0now, τnow, mlnow;
                                            input_constraint=incon, output_constraint=outcon,
                                            S_sf=S_sfnow, S_ion=S_ionnow,
                                            S_fiss=S_fissnow, S_loss=S_lossnow,
                                            Nplot)
            data[] = d
            Makie.empty!(ax3)
            isosurface!(ax3, d)
            updateslice!()
            status.text = any(d.margin .< 1) ? "Done" : "Done (design space is empty)"
            status.color = :black
        catch e
            status.text = "Error: "*shorterror(e)
            status.color = :red
        end
        nothing
    end

    Makie.on(computebutton.clicks) do _
        recompute!()
    end
    Makie.on(sliceslider.value) do _
        updateslice!()
    end

    recompute!()
    display_fig(fig)
    fig
end

end # module
