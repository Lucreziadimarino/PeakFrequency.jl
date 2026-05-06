module PeakFrequency


using FourierAnalysis, DataFrames

export  remove_aperiodic,
        find_peak,
        peakFrequency,
        reports_to_table

# -------------------------------------------------------------------------
# remove_aperiodic
# -------------------------------------------------------------------------
#
# Given a single power spectrum `spec`, sampling rate `sr`, window length `t`,
# and a frequency interval `frange` (Hz), this function estimates and removes
# the aperiodic 1/f spectral component within the selected frequency range.
#
# The function preserves the full spectral vector, but only the bins
# corresponding to `frange` are modeled and corrected.
#
# -------------------------------------------------------------------------
# MODELING STEPS
# -------------------------------------------------------------------------
#
# 1) Frequency Restriction
#    The frequency bins corresponding to `frange` are identified through `f2b`.
#
# 2) Log-Log Linear Modeling
#    The selected spectral segment is modeled in the log-log domain as:
#       log P(f) = β0 + β1 log(f) + ε(f)
#    where β0 is the intercept and β1 is the spectral slope.
#
# 3) Least-Squares Fit
#    The coefficients β = [β0, β1] are estimated by least-squares regression.
#
# 4) Aperiodic Component Removal
#    The fitted background is subtracted from the log spectrum and the
#    corrected spectrum is reconstructed in the original domain as:
#       P_corr(f) = exp(log P(f) - (β0 + β1 log(f)))
#
# 5) Full-Spectrum Output
#    The original full spectrum is preserved, but the bins inside `frange`
#    are replaced by their corrected version.
#
# -------------------------------------------------------------------------
# RETURN
# -------------------------------------------------------------------------
#
# The function returns:
#
#   spec_detr : full spectrum with the aperiodic component removed
#               within the selected range
#   β         : vector of regression coefficients [β0, β1]
#
# -------------------------------------------------------------------------
function remove_aperiodic(
    spec::AbstractVector{<:Real},
    sr::Real,
    t::Integer,
    frange
)
    a1 = f2b(first(frange), sr, t)
    a2 = f2b(last(frange), sr, t)

    f = fdf(sr, t)
    fsub = f[a1:a2]

    logSpec = log.(spec[a1:a2])
    X = [ones(length(fsub)) log.(fsub)]
    β = X \ logSpec
    bg = X * β

    spec_detr = copy(Float64.(spec))
    spec_detr[a1:a2] .= exp.(logSpec .- bg)

    return spec_detr, β
end


# -------------------------------------------------------------------------
# remove_aperiodic
# -------------------------------------------------------------------------
#
# Matrix method of `remove_aperiodic`.
#
# Given a matrix of spectra `S`, where each column corresponds to one
# spectrum, this function applies the same aperiodic-component removal
# procedure column by column.
#
# -------------------------------------------------------------------------
# RETURN
# -------------------------------------------------------------------------
#
# The function returns:
#
#   Sout   : matrix of corrected spectra, one column per input spectrum
#   betas  : matrix of regression coefficients (2 × nSpectra), where each
#            column contains [β0, β1] for the corresponding input spectrum
#
# This method makes the function compatible with generic matrix inputs,
# including views and subsets of larger spectral matrices.
# -------------------------------------------------------------------------
function remove_aperiodic(
    S::AbstractMatrix{<:Real},
    sr::Real,
    t::Integer,
    frange
)
    Sout = Matrix{Float64}(undef, size(S))
    betas = Matrix{Float64}(undef, 2, size(S, 2))

    for j in 1:size(S, 2)
        Sout[:, j], β = remove_aperiodic(view(S, :, j), sr, t, frange)
        betas[:, j] = β
    end

    return Sout, betas
end


# -------------------------------------------------------------------------
# find_peak
# -------------------------------------------------------------------------
#
# Given a single spectrum `spec`, sampling rate `sr`, window length `t`,
# a frequency interval `frange`, the maximum number of smoothing iterations
# `maxiter`, and the maximum allowed number of minima `maxmins`
# and maxima `maxmaxs`,
# this function detects a spectral peak within the selected range.
#
# This function is intentionally band-independent and can therefore be used
# as a peak detector for alpha, theta, beta, or any other spectral
# interval of interest.
#
# -------------------------------------------------------------------------
# PROCESSING STEPS
# -------------------------------------------------------------------------
#
# 1) Frequency Restriction
#    The bins corresponding to `frange` are identified through `f2b`.
#
# 2) Iterative Smoothing
#    A Blackman smoothing filter is applied iteratively until:
#       0 < number of minima ≤ maxmins
#       0 < number of maxima ≤ maxmaxs
#    or until `maxiter` is reached.
#
# 3) Local Extrema Detection
#    A bin j is defined as local maximum (or minimum) if it satisfies
#    dominance conditions over ±1, ±2 and ±3 neighboring bins.
#
# 4) Elimination of Weak Maxima
#    Maxima whose amplitude relative to the closest minimum does not exceed
#    `eliminateTh` are removed.
#
# 5) Final Peak Selection
#    If one or more maxima survive, the final peak is selected according to:
#       :max   -> highest remaining maximum
#       :first -> first remaining maximum
#       :last  -> last remaining maximum
#
# -------------------------------------------------------------------------
# RETURN
# -------------------------------------------------------------------------
#
# The function returns a NamedTuple containing:
#
#   ok        : boolean flag indicating whether a valid peak was found
#   reason    : diagnostic message
#   peakHz    : detected peak frequency in Hz
#   peakBin   : corresponding peak bin
#   mins      : vector of detected minima bins
#   maxs      : vector of detected maxima bins
#   nsmooth   : number of smoothing iterations performed
#   spectrum  : final spectrum used for peak detection
#
# -------------------------------------------------------------------------
function find_peak(
    spec::AbstractVector{<:Real},
    sr::Real,
    t::Integer,
    frange,
    maxiter::Integer,
    maxmins::Integer,
    maxmaxs::Integer;
    maxType::Symbol = :max,
    eliminateTh::Real = 1.1,
    DC::Bool = false
)

    
    a1 = f2b(first(frange), sr, t; DC=DC)
    a2 = f2b(last(frange), sr, t; DC=DC)
   

    S_ = copy(Float64.(spec))

    mins = Int[]
    maxs = Int[]
    nsmooth = 0
    foundOK = false
    reason = "OK"

    while true
        nsmooth += 1
        empty!(mins)
        empty!(maxs)

        for j = (a1+3):(a2-3)
            if S_[j-1] >= S_[j] <= S_[j+1] &&
               S_[j-2] >= S_[j] <= S_[j+2] &&
               S_[j-3] >= S_[j] <= S_[j+3]
                push!(mins, j)
            end

            if S_[j-1] <= S_[j] >= S_[j+1] &&
               S_[j-2] <= S_[j] >= S_[j+2] &&
               S_[j-3] <= S_[j] >= S_[j+3]
                push!(maxs, j)
            end
        end

        foundOK = (0 < length(mins) <= maxmins &&
                   0 < length(maxs) <= maxmaxs)

        if foundOK || nsmooth > maxiter
            break
        else
            S_ = smooth(blackmanSmoother, S_)
        end
    end

    if isempty(maxs)
        return (
            ok = false,
            reason = "No maxima found in the selected range",
            peakHz = missing,
            peakBin = missing,
            mins = copy(mins),
            maxs = copy(maxs),
            nsmooth = nsmooth,
            spectrum = S_
        )
    end

    # Remove weak maxima
    if length(maxs) > 1 && !isempty(mins)
        j = 1
        while j <= length(maxs)
            _, pos = findmin(abs.(mins .- maxs[j]))
            if S_[maxs[j]] / S_[mins[pos]] < eliminateTh
                deleteat!(maxs, j)
            else
                j += 1
            end
        end
    end

    if isempty(maxs)
        return (
            ok = false,
            reason = "All maxima eliminated by threshold",
            peakHz = missing,
            peakBin = missing,
            mins = copy(mins),
            maxs = copy(maxs),
            nsmooth = nsmooth,
            spectrum = S_
        )
    end

    # Final peak selection
    peakBin =
        if maxType == :max
            maxs[argmax(S_[maxs])]
        elseif maxType == :first
            first(maxs)
        elseif maxType == :last
            last(maxs)
        else
            error("Unsupported maxType: $maxType")
        end

    peakHz = b2f(peakBin, sr, t; DC=DC)

    return (
        ok = true,
        reason = "OK",
        peakHz = peakHz,
        peakBin = peakBin,
        mins = copy(mins),
        maxs = copy(maxs),
        nsmooth = nsmooth,
        spectrum = S_
    )
end


# -------------------------------------------------------------------------
# peakFrequency
# -------------------------------------------------------------------------
#
# Given k spectral objects in `𝐒` (one per subject), sampling rate `sr`,
# window length `t`, a computational frequency range `frange` (Hz),
# the maximum number of smoothing iterations `maxiter`,
# and the maximum allowed number of minima `maxmins`
# and maxima `maxmaxs`,
# this function estimates the Individual Peak Frequency 
# for each subject.
#
# In contrast with the original implementation, the present version:
#   (i) operates at the channel level rather than on a channel-averaged spectrum,
#   (ii) removes the aperiodic 1/f spectral background before peak detection,
#   (iii) estimates both Peak Frequency (PF) and Center of Gravity (CoG),
#   (iv) allows restricting the analysis to a selected subset of channels
#        (e.g., occipital channels only),
#   (v) relies on a band-independent peak detector (`find_peak`).
#
# -------------------------------------------------------------------------
# PROCESSING STEPS
# -------------------------------------------------------------------------
#
# 1) Frequency Restriction
#    The computational range `frange` is converted to frequency bins through `f2b`.
#
# 2) Channel Selection
#    If `channels` is provided, only the selected channels are processed.
#    Otherwise, all channels are included in the analysis.
#
# 3) Aperiodic Component Removal
#    For each selected channel, the 1/f background is estimated and removed
#    by calling `remove_aperiodic`.
#
# 4) Peak Detection
#    The corrected spectrum of each channel is passed to `find_peak`,
#    which performs smoothing, local-extrema detection, weak-peak elimination,
#    and final peak selection.
#
# 5) Channel-wise PF Estimation
#    For each valid channel:
#       PF_ch = detected peak frequency returned by `find_peak`
#
# 6) Channel-wise CoG Estimation
#    CoG is computed as:
#       CoG_ch = Σ(f_i * P_corr(f_i)) / Σ(P_corr(f_i))
#
# 7) Subject-level Aggregation
#    Subject-level PF and CoG are computed as the mean of the channel-wise
#    estimates across valid channels.
#
# 8) Minimum Channel Constraint
#    A subject estimate is considered valid only if the number of valid channels
#    is at least `minValidChannels`.
#
# -------------------------------------------------------------------------
# RETURN
# -------------------------------------------------------------------------
#
# The function returns:
#
#   OK        : global boolean, true if all subjects have a valid PF estimate
#   PF       : vector of subject-level Peak Frequency values (Hz)
#   CoG       : vector of subject-level Center of Gravity values (Hz)
#   OKvec     : boolean vector indicating whether each subject has a valid estimate
#   reports   : vector of NamedTuples, one per subject, containing:
#
#       PF             : subject-level Peak Frequency
#       CoG             : subject-level Center of Gravity
#       ok              : validity flag for the subject
#       reason          : diagnostic message
#       f               : frequency vector in the selected computational range
#       RawSpec         : raw spectra matrix (frequency × channels)
#       DetrendedSpec   : 1/f-corrected spectra matrix (frequency × channels)
#       validChannels   : indices of channels contributing to the final estimate
#       betas           : matrix of 1/f regression coefficients [β0, β1]
#       sr              : sampling rate
#       wl              : window length
#       DC              : DC correction flag
#
# -------------------------------------------------------------------------
function peakFrequency(𝐒::AbstractVector, sr, t, frange, maxiter, maxmins, maxmaxs;
                                 maxType::Symbol = :max,
                                 eliminateTh::Real = 1.1,
                                 minValidChannels::Int = 1,
                                 channels::AbstractVector, 
                                 DC::Bool = false)

    
    f = fdf(sr, t; DC=DC)
    a1 = f2b(first(frange), sr, t; DC=DC)
    a2 = f2b(last(frange), sr, t; DC=DC)

    k = length(𝐒)

    OKvec = falses(k)
    PF   = Vector{Union{Missing,Float64}}(undef, k)
    CoG   = Vector{Union{Missing,Float64}}(undef, k)
    reports = Vector{NamedTuple}(undef, k)

    for i = 1:k

        subj = 𝐒[i]
        nchan = size(subj.y, 2)

        peakList = Float64[]
        cogList  = Float64[]
        rawMat   = fill(NaN, a2-a1+1, nchan)
        detrMat  = fill(NaN, a2-a1+1, nchan)
        betaMat  = fill(NaN, 2, nchan)
        validChannelIdx = Int[]

        ok = true
        reason = "OK"

        use_channels = isnothing(channels) ? collect(1:nchan) : channels

        # Remove aperiodic component on the full spectrum only for selected channels
        Ysel = Matrix(subj.y[:, use_channels])
        Ysel_detr, betas_sel = remove_aperiodic(Ysel, sr, t, frange)

        for (loc, ch) in enumerate(use_channels)

            # Full raw and full detrended spectrum for this channel
            spec_raw  = subj.y[:, ch]
            spec_detr = Ysel_detr[:, loc]

            # Save only the computational range for plotting/reporting
            rawMat[:, ch]  .= spec_raw[a1:a2]
            detrMat[:, ch] .= spec_detr[a1:a2]
            betaMat[:, ch] .= betas_sel[:, loc]

            pk = find_peak(
                spec_detr,
                sr,
                t,
                frange,
                maxiter,
                maxmins,
                maxmaxs;
                maxType = maxType,
                eliminateTh = eliminateTh,
                DC = DC
            )

            if !pk.ok
                continue
            end

            # Channel-wise PF
            push!(peakList, pk.peakHz)
            push!(validChannelIdx, ch)

            # Channel-wise CoG computed on the final spectrum used for peak detection
            α1 = f2b(first(frange), sr, t; DC=DC)
            α2 = f2b(last(frange), sr, t; DC=DC)

            Sused = pk.spectrum
            num = sum(f[j] * Sused[j] for j in α1:α2)
            den = sum(Sused[j] for j in α1:α2)

            if den > 0
                push!(cogList, num / den)
            end
        end

        # Subject-level aggregation
        if length(peakList) >= minValidChannels
            PF[i] = mean(peakList)
        else
            PF[i] = missing
            ok = false
            reason = "Insufficient valid channels (PF)"
        end

        if length(cogList) >= minValidChannels
            CoG[i] = mean(cogList)
        else
            CoG[i] = missing
            ok = false
            reason = "Insufficient valid channels (CoG)"
        end

        OKvec[i] = ok

        reports[i] = (
            PF = PF[i],
            CoG = CoG[i],
            ok  = ok,
            reason = reason,
            f   = f[a1:a2],
            RawSpec = rawMat,
            DetrendedSpec = detrMat,
            validChannels = validChannelIdx,
            betas = betaMat,
            sr = sr,
            wl = t,
            DC = DC,
            taper = nothing,
            func = identity,
            smoother = nothing
        )
    end

    OK = all(OKvec)

    return OK, PF, CoG, OKvec, reports
end


# -------------------------------------------------------------------------
# reports_to_table
# -------------------------------------------------------------------------
#
# Given the vector of subject-level reports returned by `peakFrequency,
# this function converts the main estimation outputs into a tabular format.
#
# -------------------------------------------------------------------------
# OUTPUT TABLE
# -------------------------------------------------------------------------
#
# The resulting DataFrame contains one row per subject and the following columns:
#
#   OK        : boolean flag indicating whether the subject-level estimate is valid
#   PF_Hz    : subject-level Peak Frequency (Hz)
#   CoG_Hz    : subject-level Center of Gravity (Hz)
#   Valid     : indices of channels contributing to the estimate
#
# This function is intended for quick inspection, reporting, and downstream analysis.
# -------------------------------------------------------------------------

function reports_to_table(reports)
    DataFrame(
        OK     = getfield.(reports, :ok),
        PF_Hz = getfield.(reports, :PF),
        CoG_Hz = getfield.(reports, :CoG))
end

end # module