#=

This example shows how to ...

To run the example, add the following packages to your environment:
Eegle, Plots, Plots.PlotMeasures, LinearAlgebra, Statistics, DSP, DataFrames 

=#


############################
# 0. SETUP AND DEPENDENCIES
############################

push!(LOAD_PATH, joinpath(@__DIR__))
using PeakFrequency
using Eegle, Plots, Plots.PlotMeasures, LinearAlgebra, Statistics, DSP, DataFrames 


plotSpectra(S, sr, wl, DC, taper, flabels, func, smoother, plotArgs) =
  plot(Spectra(S, sr, wl, DC, taper, flabels, func, smoother); plotArgs...)


############################
# 1. DATA INPUT
############################

dir = joinpath(abspath(homedir(), "...", "..."), "Example Data")  #Insert your own path



files=getFilesInDir(dir, ext=(".txt",))
session= "T1"
condition= "oc"
ph = "pre"

files = filter(file -> contains(file, session) && contains(file, condition) && contains(file, ph), files)
𝐗=readASCII(files)


const leads = ["FP1","FP2","FZ","C5", "CZ", "C6","O1", "O2"]

## general parameters: sr=sampling rate
# maxf= max freq. for analysis; k=#files; n=#electrodes;
sr, maxf, k, n = 128, 32, length(𝐗), length(leads)

# spectra plot arguments to make nice plots (tuple)
spPlotArgs=(label=["FP1" "FP2" "FZ" "C5" "CZ" "C6" "O1" "O2"],
            fmax = maxf,
             left_margin = 2mm,
             bottom_margin = 2mm,
             ytitle = "Amplitude (\\muV)",
             xtickfont = font(11, "Times"),
             ytickfont = font(11, "Times"),
             legend = true) 


############################
# 2. AMPLITUDE SPECTRA (R, N)
############################

# parameters for finding peak alpha frequencies
#wl=1024; #2048; # FFT window length
wl=128;
DC=false
maxb=f2b(32, 128, 128; DC=DC); # FFT bin corresponding to maxf
#maxb=f2b(maxf, sr, wl; DC=DC); # FFT bin corresponding to maxf
slide=wl÷2; # Welch sliding in samples: FFT resolution=sr/wl
tape=slepians(sr, wl, golden);
# utility for 'plotSaveSpectra' function
spArgs=(sr, wl, DC, string(tape.kind),
        fdf(sr, wl; DC=DC)[1:maxb], #flabels
        sqrt, noSmoother) # func and smoothing

# compute ampitude spectra for all files (takes some time)
plan=Planner(plan_exhaustive, 10.0, wl) #pre-compute a planner to go fast

# Raw amplitude spectra
𝐒=spectra([X for X∈𝐗], sr, wl;
   tapering=tape, slide=slide, planner=plan,func=sqrt);

for i in 1:length(𝐗)
p = plotSpectra(real(𝐒[i].y), spArgs..., spPlotArgs) 
display(p)
end

############################
# 3. INDIVIDUAL ALPHA FREQUENCY
############################


OK, PAF, CoG, OKvec,reports = peakFrequency(
    𝐒,
    sr,
    wl,
    5:15,
    20, 2, 2;
    maxType = :max,
    channels = [7,8]);

df = reports_to_table(reports)
display(df)


######## PLOT #########


for i in eachindex(𝐒)
    p = plotSpectra(real(𝐒[i].y), spArgs..., spPlotArgs)

    if !ismissing(PAF[i])
        vline!(p, [PAF[i]],
            label = "PAF",
            linewidth = 2)
    end

    if !ismissing(CoG[i])
        vline!(p, [CoG[i]],
            label = "CoG",
            linestyle = :dash,
            linewidth = 2)
    end

    display(p)
end


for i in eachindex(𝐒)

    meanSpec = vec(mean(real(𝐒[i].y), dims=2))

    p = plot(
        fdf(sr, wl; DC=DC)[1:maxb],
        meanSpec[1:maxb],
        title = "Mean spectrum across channels | subject $i",
        xlabel = "Frequency (Hz)",
        ylabel = "Amplitude (\\muV)",
        linewidth = 2,
        label = "Mean spectrum"
    )

    if !ismissing(PAF[i])
        vline!(p, [PAF[i]], label="PAF", linewidth=2)
    end

    if !ismissing(CoG[i])
        vline!(p, [CoG[i]], label="CoG", linestyle=:dash, linewidth=2)
    end

    display(p)
end