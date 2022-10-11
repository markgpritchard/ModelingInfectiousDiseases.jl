
module Chapter2Additions

# Additional functions that were not included in the programmes that accompany the book

using CairoMakie, ..MID_23, ..MID_24
import ..MID_23: plot_sir23_vals, plot_sir23!
import ..MID_24: plot_sir24_vals, plot_sir24!

export plot_sirs23_24

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Results of programmes 2.3 and 2.4 side-by-side 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This function is not included in the programmes that accompany the book, but 
# Figure 2.8 does compare the density-dependent and frequency-dependent results. 
# Here, we plot the results of the two programmes side-by-side using the same
# parameters for each.

function plot_sirs23_24(; beta = 520 / 365, gamma = 1 / 7, mu = 1 / (70 * 365), nu = 1 / (70 * 365), 
        rho = .5, X0 = .2, Y0 = 1e-6, N0 = 1, duration = 100 * 365)

    sol23 = run_sir23(; beta, gamma, mu, nu, rho, X0, Y0, N0, duration, saveat = 1)
    sol24 = run_sir24(; beta, gamma, mu, nu, rho, X0, Y0, N0, duration, saveat = 1)

    fig = Figure()
    ax1 = Axis(fig[1, 1]); ax2 = Axis(fig[2, 1]); ax3 = Axis(fig[3, 1])
    ax4 = Axis(fig[1, 2]); ax5 = Axis(fig[2, 2]); ax6 = Axis(fig[3, 2])

    xs23, X23, Y23, Z23, N23 = plot_sir23_vals(sol23)
    xs24, X24, Y24, Z24, N24 = plot_sir24_vals(sol24)

    plot_sir23!(ax1, ax2, ax3, xs23, X23, Y23, Z23, N23)
    plot_sir24!(ax4, ax5, ax6, xs24, X24, Y24, Z24, N24)
    fig[3, 3] = Legend(fig, ax3)

    linkyaxes!(ax1, ax4); linkyaxes!(ax2, ax5); linkyaxes!(ax3, ax6)
    for ax ∈ [ax4, ax5, ax6] hideydecorations!(ax; grid = false) end 

    Label(fig[0, :], "p2.3.jl and p2.4.jl: SIR models with infection-induced mortality")
    ax1.title = "Density-dependent transmission"; ax1.titlefont = "Makie"
    ax4.title = "Frequency-dependent transmission"; ax4.titlefont = "Makie"
    resize_to_layout!(fig)

    return fig 
end 

end # module Chapter2Additions
