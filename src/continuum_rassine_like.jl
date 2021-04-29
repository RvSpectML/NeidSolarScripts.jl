module Continuum

 using Statistics
 using DSP
 using DataInterpolations
 using SortFilters, RollingFunctions, NaNMath
 #https://github.com/sairus7/SortFilters.jl

# constants
speed_of_light_mks = 3e8 # m/s
fwhm_sol = 7.9e3 # m/s

function replace_nans!(f::AV2; half_width::Integer = 6 ) where {  T2<:Real, AV2<:AbstractVector{T2} }
  idx_nan = findall(isnan.(f))
  for i in idx_nan
  end


end

function smooth(f::AV2; half_width::Integer = 6 ) where {  T2<:Real, AV2<:AbstractVector{T2} }
  w=DSP.hanning(1+2*half_width)
  if 1+2*half_width<3
    return f
  elseif length(f)<1+2*half_width
    return f
  else
    return y_smooth = conv(f, w/sum(w))[1+half_width:end-half_width]
  end
end

function calc_rolling_max(f::AV2; width::Integer = 14 ) where {  T2<:Real, AV2<:AbstractVector{T2} }
 shift = floor(Int64,width//2)
 z = similar(f)
 #z[(shift):(length(f)-shift)] .= rollmax(f,width)
 z[(shift):(length(f)-shift)] .= rolling(NaNMath.maximum,f,width)
 z[1:shift-1] .= z[shift]
 z[(length(f)-shift+1):end] .= z[length(f)-shift]
 return z
end

function find_local_maxima(f::AV2; half_width::Integer = 7 ) where {  T2<:Real, AV2<:AbstractVector{T2} }
  rollingmax = calc_rolling_max(f,width=half_width*2)
  findall(f .== rollingmax)
end

function calc_rolling_median(λ::AV1, f::AV2; width::Integer = 4, verbose::Bool = false ) where { T1<:Real, T2<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2} }
  @assert length(λ) == length(f)
  @assert width >= 2
  @assert width%2 == 0
  shift = floor(Int64,width//2)
  f_filtered = similar(f)
  #if any(isnan.(f))
  #   f_filtered .= running(NaNMath.median,f,width)
  #else
     f_filtered[(1+shift):(length(f)-shift)] .= movsort(f, width, 0.5)
     f_filtered[1:shift] .= f_filtered[1+shift]
     f_filtered[(length(f)-shift+1):end] .= f_filtered[length(f)-shift]
  #end
  return f_filtered
end


function calc_rolling_median_and_max_to_clip(λ::AV1, f::AV2; width::Integer = 8, cr_pad_width::Integer = 1,
          num_clippings::Integer = 3, verbose::Bool = false ) where { T1<:Real, T2<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2} }
  @assert length(λ) == length(f)
  @assert width >= 4
  @assert width%2 == 0
  @assert 0<= cr_pad_width <= 7
  @assert 1 <= num_clippings <= 10
  shift = floor(Int64,width//2)
  high_quantile = (length(λ) > 200) ? 0.99 : 0.95
  p = (0.25, 0.5, 0.75, high_quantile)
  function iqr(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} }
      (mq[i][3]-mq[i][1])    # standard iqr
  end
  function iqr_upper(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} }
      2*(mq[i][3]-mq[i][2])  # double Q3-Q2
  end
  q75(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} } = mq[i][3]
  qhi(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} } = mq[i][4]

  output = similar(λ)
  mask = trues(length(λ))
  f_median_filtered = similar(f)
  f_clean = deepcopy(f)

  if verbose flush(stdout)  end
  for i in 1:num_clippings
     moving_quartiles = movsort(f_clean, width, p)
     f_median_filtered[(1+shift):(length(f)-shift)] .= map(i->moving_quartiles[i][2], (1+2*shift):length(moving_quartiles) )
     f_median_filtered[1:shift] .= f_median_filtered[1+shift]
     f_median_filtered[(length(f)-shift+1):end] .= f_median_filtered[length(f)-shift]
     max_iqr = map(i->q75(moving_quartiles,i)+1.5*iqr_upper(moving_quartiles,i),1:length(λ))
     idx = (1+2*shift):length(moving_quartiles)
     output[(1+shift):(length(f)-shift)] .= max_iqr[idx]
     #=
     # Rassine uses a second pass with simple high quantile with a wider window (to avoid being affected by lines).
     # For now we use the same window
     width_pass2 = 100  # TODO: Automate for orders.  FOr now arbtrary rather than 100A, while testing on small chunks
     shift_pass2 = floor(Int64,width_pass2//2)
     max_99 = movsort(f, width_pass2, 0.99)
     output[(1+shift):(length(f)-shift)] .= max.(output[idx]),max_99[idx])
     =#
     #output[(1+shift):(length(f)-shift)] .= max.(output[idx],map(i->qhi(moving_quartiles,i),idx))
     output[(1+shift):(length(f)-shift)] .= map(i->max(output[i],qhi(moving_quartiles,i)),idx)

     # Fill in edges
     output[1:shift] .= output[1+shift]
     output[(length(f)-shift+1):end] .= output[length(f)-shift]

     # Identify & replace high outliers, plus neighboring pixels
     for i in 1:length(mask)
       if f[i] <= output[i] continue end
       for j in max(1,i-cr_pad_width):min(i+cr_pad_width,length(mask))
          mask[j] = false
       end
     end
     #mask .&= (f.<=output)
     if verbose println("# pass ", i, ", removing ", sum((f.>output)), " hi pixels for a total of ", sum(.!mask), ".") end
     f_clean[.!mask] .= f_median_filtered[.!mask]
  end
  return (f_filtered=f_median_filtered, f_threshold=output, f_clean=f_clean)
end

function longest_cluster_same_sign(y::AV) where { T<:Real, AV<:AbstractVector{T} }
  longest = 0
  current = 1
  this_start = 1
  start = 1
  stop = 1
  sign_last = sign(first(y))
  for i in 2:length(y)
    sign_this = sign(y[i])
    if sign_this == sign_last
      current += 1
    else
      if current > longest
        longest = current
        start = this_start
        stop = i-1
      end
      current = 1
      this_start = i
    end
    sign_last = sign_this
  end
  if current > longest
    longest = current
    start = this_start
    stop = length(y)
  end
  return (len=longest, start=start, stop=stop)
end

function calc_rollingpin_radius(λ::AV1, f::AV2; fwhm::Real = fwhm_sol,
            S1_width_factor::Real = 40, S2_width_factor::Real = 10,
            min_R_factor::Real = 100, ν::Real = 1, verbose::Bool = false ) where { T1<:Real, AV1<:AbstractVector{T1}, T2<:Real, AV2<:AbstractVector{T2}  }
   λlo = minimum(λ)
   S1_width = S1_width_factor*fwhm/speed_of_light_mks*λlo
   S2_width = S2_width_factor*S1_width
   S1_width = 2*floor(Int64,S1_width/2)
   S2_width = 2*floor(Int64,S2_width/2)
   f_max1 = calc_rolling_max(f,width=S1_width)
   f_max2 = calc_rolling_max(f,width=S2_width)
   (longest_len, longest_start, longest_stop) = longest_cluster_same_sign(f_max2.-f_max1)
   longest_Δλoverλ = (λ[longest_stop]-λ[longest_start]) / ((λ[longest_start]+λ[longest_stop])/2)
   #penalty = (f_max2.-f_max1)./f_max2
   #conversion_fwhm_sig = (10*min_λ/(sqrt(8*log(2))*speed_of_light_kmps))
   #min_radius = fwhm * conversion_fwhm_sig # par_R
   #radius = min_radius .* λ./min_λ
   rollingpin_R_min = min_R_factor*λlo*(fwhm/speed_of_light_mks)/log(8*log(2))
   #rollingpin_R_max = max_R_factor*rollingpin_R_min  # Arbitrary, not same as Rassine's auto setting, see https://github.com/MichaelCretignier/Rassine_public/blob/master/Rassine.py#L718
   rollingpin_R_max = max(λlo*longest_Δλoverλ,rollingpin_R_min*2)
   #rollingpin_R_min = 2*floor(Int64,rollingpin_R_min/2)
   #rollingpin_R_max = 2*floor(Int64,rollingpin_R_max/2)
   if verbose
     println("longest_len = ", longest_len)
     println("longest_Δλoverλ = ", longest_Δλoverλ)
     println("rollingpin_R_min = ", rollingpin_R_min)
     println("rollingpin_R_max = ", rollingpin_R_max)
   end
   rollingpin_radius = (λ./λlo).*(rollingpin_R_min .+ (rollingpin_R_max-rollingpin_R_min).*((f_max2.-f_max1)./f_max2).^ν)
   #plot!(lambda,f_max1,color=:red)
   #plot!(lambda,f_max2,color=:magenta)
   #plot!(lambda,penalty)
   return rollingpin_radius
end

function calc_continuum_anchors(λ::AV1, f::AV2; radius::AV3,
            stretch_factor::Real = 1.0, fwhm::Real = 7.3 #= km/s=#,
                        verbose::Bool = false ) where { T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3} }
  @assert length(λ) == length(f)
  #speed_of_light_kmps = 3e5
  #par_stretching = 2# 0.5# 0.5 #11.385
  extrema_λ = extrema(λ)
  min_λ = first(extrema_λ)
  extrema_f = extrema(f)
  normalization = (extrema_f[2]-extrema_f[1])/(extrema_λ[2]-extrema_λ[1])
  normalization *= stretch_factor
  #= v1
  waves = λ.-λ'
  mask = sign.(waves).>0
  #distance = sign.(waves).*sqrt.(waves.^2 .+((f.-f')./normalization).^2)
  distance[mask] = sqrt.(waves[mask].^2 .+((f.-f')[mask]./normalization).^2)
  =#
  #= v2
  distance = zeros(length(λ),length(λ))
  mask = sign.(λ.-λ').>0
  distance[mask] .= sqrt.((λ.-λ')[mask].^2 .+((f.-f')[mask]./normalization).^2)
  =#
  #= v3
  distance = zeros(length(λ),length(λ))
  for i in 1:length(λ)
    for j in 1:length(λ)
      if λ[i]<=λ[j] continue end
      distance[i,j] = sqrt((λ[i]-λ[j])^2 +((f[i]-f[j])/normalization)^2)
    end
  end
  #return distance
  =#
  # v4
  function calc_dist(i::Integer, j::Integer)  # TODO OPT: eliminate sqrt
    if λ[i]<=λ[j]
      return 0.0
    else
      sqrt((λ[i]-λ[j])^2 +((f[i]-f[j])/normalization)^2)
    end
  end
  numero = range(1,stop=length(λ))
  mask = falses(length(λ))
  dist = zeros(length(λ))
  #conversion_fwhm_sig = (10*min_λ/(sqrt(8*log(2))*speed_of_light_kmps))
  #min_radius = fwhm * conversion_fwhm_sig # par_R
  #println("# old min_radius = ", min_radius)
  if verbose println("# passed radius = ", extrema(radius)) end
  #radius = min_radius .* λ./min_λ
  keep = Vector{Int64}()
  sizehint!(keep,max(32,min(256,2^floor(Int64,log2(length(λ))-2))))
  push!(keep,1)
  j = 1
  while length(λ)-j>2
    par_R = radius[j]
    if j==1 par_R /= 1.5 end
    map!(i->calc_dist(i,j),dist,1:length(λ))
    #mask .= (0.0 .< distance[:,j].<2.0*par_R)
    mask .= (0.0 .< dist.<2.0*par_R)
    while sum(mask)==0
      par_R *= 1.5
      #mask .= (0.0 .< distance[:,j].<2.0*par_R)
      mask .= (0.0 .< dist.<2.0*par_R)
    end
    p1 = [ λ[j] f[j]/normalization ]
    p2 = hcat(λ[mask], f[mask]./normalization)
    delta = p2 .- p1
    #c = sqrt.(delta[:,1].^2 .+delta[:,2].^2)
    c = sqrt.(vec(sum(delta.^2,dims=2)))
    harg = (par_R.^2 .-0.25.*c.^2)
    #=
    if any(harg.<0)
      println("# j = ", j, " dist = ", distance[j:j+10,j])
      println("# R^2 = ", par_R^2)
      println("# c^2/4 = ", c.^2/4)
      println("# harg = ", harg)
    end
    h = zeros(length(harg))
    hmask = harg.>0
    h[hmask] .= sqrt.(harg[hmask])
    =#
    h = sqrt.(harg)
    cx = p1[1] .+ 0.5.*delta[:,1] .-h./c.*delta[:,2]
    cy = p1[2] .+ 0.5.*delta[:,2] .+h./c.*delta[:,1]
    #return (cx,cy)
    #cond1 = (cy.-p1[2]).>=0
    #theta = zeros(length(cy))
    mintheta = Inf
    argmintheta = 0
    for i in 1:length(cy)
      if cy[i]-p1[2] >0
        thetai = -acos((cx[i]-p1[1])/par_R)+π
      else
        thetai = -asin((cy[i]-p1[2])/par_R)+π
      end
      if thetai < mintheta
        mintheta = thetai
        argmintheta = i
      end
    end
    #theta[cond1] .= -acos.((cx[cond1].-p1[1])./par_R).+π
    #theta[.!cond1] .= -asin.((cy[.!cond1].-p1[2])./par_R).+π
    #j2 = argmin(theta)
    j2 = argmintheta
    j = numero[mask][j2]
    push!(keep,j)
  end
  if verbose println("# using ", length(keep), " anchor points") end
  return keep
end

function calc_continuum_from_anchors( λ::AV1, f::AV2, anchors::AV3; λout::AV4 = λ,
              verbose::Bool = false ) where {
                T1<:Real, T2<:Real, T3<:Integer, T4<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3}, AV4<:AbstractVector{T4}  }
  calc_continuum_from_anchors_cubic( λ, f, anchors, λout=λout,verbose=verbose)
end

function calc_continuum_from_anchors_linear( λ::AV1, f::AV2, anchors::AV3; λout::AV4 = λ,
              verbose::Bool = false ) where {
                T1<:Real, T2<:Real, T3<:Integer, T4<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3}, AV4<:AbstractVector{T4}  }
  #interp_cubic = CubicSpline(f[anchors],λ[anchors])
  interp_linear = LinearInterpolation(f[anchors],λ[anchors])
  function extrap(x::Real)
    if x<first(interp_linear.t)
      return first(interp_linear.u)
    elseif x>last(interp_linear.t)
      return last(interp_linear.u)
    else
      return interp_linear(x)
    end
  end
  return extrap.(λout)
end

function calc_continuum_from_anchors_cubic( λ::AV1, f::AV2, anchors::AV3; λout::AV4 = λ,
              verbose::Bool = false ) where {
                T1<:Real, T2<:Real, T3<:Integer, T4<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3}, AV4<:AbstractVector{T4}  }
  interp_cubic = CubicSpline(f[anchors],λ[anchors])
  interp_cubic.(λout)
  #=
  interp_linear = LinearInterpolation(f[anchors],λ[anchors])
  function extrap(x::Real)
    if x<first(interp_linear.t)
      return first(interp_linear.u)
    elseif x>last(interp_linear.t)
      return last(interp_linear.u)
    else
      return interp_linear(x)
    end
  end
  return extrap.(λout)
  =#
end

function replace_edge_anchor_vals!(f::AV1, n::Integer = 3 ) where { T1<:Real, AV1<:AbstractVector{T1} }
  @assert length(f) >= 2*n+1
  f[1:(n+1)] .= f[n+1]
  f[end-n:end] .= f[length(f)-n-1]
  return f
end

function find_clean_anchors_by_slope(anchors::AV1, f::AV2; threshold::Real = 0.995, verbose::Bool = false ) where { T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3} }
  nanchors = length(anchors)
  @assert nanchors == length(f)
  @assert nanchors >= 8
  Δabsslope = zeros(length(anchors))
  for i in 2:nanchors-1
    slope_hi = (f[i+1]-f[i])/(anchors[i+1]-anchors[i])
    slope_lo = (f[i]-f[i-1])/(anchors[i]-anchors[i-1])
    Δabsslope[i] = abs(slope_hi-slope_lo)
  end
  threshold = quantile(Δabsslope[2:end-1],threshold)
  mask = Δabsslope.<threshold
  return mask
end

function filter_anchor_groups(anchors::AV1; min_num::Integer = 3, verbose::Bool = false ) where { T1<:Real, AV1<:AbstractVector{T1} }
  Δpix = anchors[2:end] .- anchors[1:end-1]
  nanchors = length(anchors)
  @assert nanchors == length(f)
  @assert nanchors >= 8
  Δabsslope = zeros(length(anchors))
  for i in 2:nanchors-1
    slope_hi = (f[i+1]-f[i])/(anchors[i+1]-anchors[i])
    slope_lo = (f[i]-f[i-1])/(anchors[i]-anchors[i-1])
    Δabsslope[i] = abs(slope_hi-slope_lo)
  end
  threshold = quantile(Δabsslope[2:end-1],threshold)
  mask = Δabsslope.<threshold
  return mask
end


function merge_nearby_anchors(anchors::AV1, λ::AV2; threshold::Real = 1, verbose::Bool = false ) where { T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3} }
  Δλ = λ[anchors[2:end]] .- λ[anchors[1:end-1]]
  close_anchor_pair_mask = Δλ .<= threshold
  nclose = sum(close_anchor_pair_mask)
  if verbose println("# Found ", nclose, " close anchor pairs out of ", length(anchors), ".")  end
  if nclose == 0
    return anchors
  else
    println("# Anchors = ", anchors)
    println("Close anchor pairs: ", findall(close_anchor_pair_mask), " => ", anchors[findall(close_anchor_pair_mask)])
  end
  old_anchor_mask = falses(length(anchors))
  if first(Δλ) > threshold
    old_anchor_mask[1] = true
  end
  for i in 2:length(Δλ)
    if (Δλ[i-1] > threshold) && (Δλ[i] > threshold)
      old_anchor_mask[i] = true
    end
  end
  if last(Δλ) > threshold
    old_anchor_mask[length(anchors)] = true
  end
  old_anchors = anchors[old_anchor_mask]
  println("# Old anchors = ", old_anchors)
  close_anchor_idx = (1:(length(anchors)-1))[close_anchor_pair_mask]
  new_anchors = zeros(Int64,nclose)
  resize!(new_anchors,0)
  idx_start = 1
  for i in 2:nclose
    if close_anchor_idx[i]-close_anchor_idx[i-1] == 1
      continue
    else
      idx_stop = i-1

      if idx_start == idx_stop
        merged_anchor_idx = idx_start
      elseif (idx_stop-idx_start)%2 == 0
        merged_anchor_idx = idx_start + convert(Int64,(idx_stop-idx_start)//2)
      else
        merged_anchor_idx = idx_start + convert(Int64,(idx_stop-idx_start+1)//2)
      end
      #=
      num_in_group = idx_stop-idx_start+1
      min_equidist = Inf
      merged_anchor_idx = 0
      for j in 0:(num_in_group-1)
        λ_anchor_prev = λ[anchors[max(close_anchor_idx[idx_start]-1,1)]]
        λ_anchor_next = λ[anchors[min(close_anchor_idx[idx_stop]+1,length(anchors))]]
        λthis = λ[anchors[close_anchor_idx[idx_start+j]]]
        equidist = abs((λthis-λ_anchor_prev)-(λ_anchor_next-λthis))
        #equidist = abs(log(λthis^2/(λ_anchor_prev*λ_anchor_next)))
        if equidist<min_equidist
          min_equidist = equidist
          merged_anchor_idx = idx_start+j
        end
      end
      =#
      merged_anchor = anchors[close_anchor_idx[merged_anchor_idx]]
      push!(new_anchors,merged_anchor)
      idx_start = i
    end
  end
  if idx_start == nclose == 1
    push!(new_anchors,anchors[close_anchor_idx[idx_start]])
  elseif idx_start != nclose
    idx_stop = nclose
    if idx_start == idx_stop
      merged_anchor_idx = idx_start
    elseif (idx_stop-idx_start)%2 == 0
      merged_anchor_idx = idx_start + convert(Int64,(idx_stop-idx_start)//2)
    else
      merged_anchor_idx = idx_start + convert(Int64,(idx_stop-idx_start+1)//2)
    end
    #=
    num_in_group = idx_stop-idx_start+1
    min_equidist = Inf
    merged_anchor_idx = 0
    for j in 0:(num_in_group-1)
      λ_anchor_prev = λ[anchors[max(close_anchor_idx[idx_start]-1,1)]]
      λ_anchor_next = λ[anchors[min(close_anchor_idx[idx_stop]+1,length(anchors))]]
      λthis = λ[anchors[close_anchor_idx[idx_start+j]]]
      equidist = abs(log(λthis^2/(λ_anchor_prev*λ_anchor_next)))
      if equidist<min_equidist
        min_equidist = equidist
        merged_anchor_idx = idx_start+j
      end
    end
    =#
    merged_anchor = anchors[close_anchor_idx[merged_anchor_idx]]
    push!(new_anchors,merged_anchor)
  end
  if verbose println("# Keeping ", length(new_anchors), " new anchors, ", length(old_anchors), " old anchors.")  end
  anchors_out = mergesorted(old_anchors,new_anchors)
  return anchors_out
end

function calc_continuum(λ::AV1, f_obs::AV2; λout::AV3 = λ, fwhm::Real = fwhm_sol, ν::Real =1.0,
  stretch_factor::Real = 1.0, merging_threshold::Real = 0.04, verbose::Bool = false ) where { T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3} }
 clip_width_A = 1.0
 clip_width_pix = 2*floor(Int64,clip_width_A/(λ[2]-λ[1])/2)
 f_smooth = Continuum.smooth(f_obs, half_width=6)
 idx_local_maxima = find_local_maxima(f_smooth, half_width=7)
 #(f_median_filtered, f_clip_threshold, f_clean) = calc_rolling_median_and_max_to_clip(λ,f_obs, width=clip_width_pix, verbose=verbose)
 #rollingpin_radius = calc_rollingpin_radius(λ, f_median_filtered, fwhm=fwhm)
 rollingpin_radius = calc_rollingpin_radius(λ, f_smooth, fwhm=fwhm, verbose=false, ν=ν)

 anch_orig = calc_continuum_anchors(λ[idx_local_maxima],f_smooth[idx_local_maxima],radius=rollingpin_radius[idx_local_maxima], stretch_factor=stretch_factor, verbose=verbose)
 anch_orig = idx_local_maxima[anch_orig]
 anchor_vals = f_smooth[anch_orig]
 replace_edge_anchor_vals!(anchor_vals)
 anch_mask = find_clean_anchors_by_slope(anch_orig,anchor_vals, verbose=verbose)
 #anchors_merged = anch_orig[anch_mask]
 # #=
 if merging_threshold > 0
     anchors_merged = merge_nearby_anchors(anch_orig[anch_mask],λ, threshold=merging_threshold, verbose=verbose)
 else
    anchors_merged = anch_orig[anch_mask]
 end
 # =#
 #f_median_filtered = calc_rolling_median(λ,f_obs, width=4)
 if length(anchors_merged)>10
   continuum = calc_continuum_from_anchors_cubic(λ,f_smooth,anchors_merged, λout=λout, verbose=verbose)
 else
   continuum = calc_continuum_from_anchors_linear(λ,f_smooth,anchors_merged, λout=λout, verbose=verbose)
 end
 return (anchors=anchors_merged, continuum=continuum, f_filtered=f_smooth)
end


# Code for mergesorted from
# https://discourse.julialang.org/t/looking-for-merge-algorithm/47473/2
import Base.Order: Ordering, Forward, ord, lt

indices(x,i)=eachindex(x)

# merge sorted vectors vl and vr into v
# from indices lo to hi in v
function mergesorted!(v::AbstractVector,
                      lo::Int, hi::Int,
                      vl::AbstractVector,
                      lol::Int, hil::Int,
                      vr::AbstractVector,
                      lor::Int, hir::Int,
                      order::Ordering)
    c = lol
    p = lor
    nl = hil
    nr = hir
    i = lo
    @inbounds while c <= nl && p <= nr && i <= hi
        if lt(order, vr[p], vl[c])
            v[i] = vr[p]
            p = p+1
            i = i+1
        else
            v[i] = vl[c]
            c = c+1
            i = i+1
        end
    end
    @inbounds while p <= nr && i <= hi
        v[i] = vr[p]
        i = i+1
        p = p+1
    end
    @inbounds while c <= nl && i <= hi
        v[i] = vl[c]
        i = i+1
        c = c+1
    end
    v
end

function mergesorted!(v::AbstractVector, vl::AbstractVector,
                      vr::AbstractVector, order::Ordering)
    inds = indices(v,1)
    indsl = indices(vl,1)
    indsr = indices(vr,1)
    mergesorted!(v,first(inds),last(inds),vl,first(indsl),last(indsl),
                 vr,first(indsr),last(indsr),order)
end

"""
    mergesorted!(v, vl, vr; lt=isless, by=identity, rev::Bool=false, order::Ordering=Forward)
Merge sorted vectors `vl` and `vr`, overwriting vector `v`. Assumes
that `vl` and `vr` are sorted and does not check whether `vl` or `vr`
are sorted. You could used `issorted` to check if `vl` and `vr` are sorted.
If length of `v` is less than the sum of the lengths of
`vl` and `vr`, this simply stops when all indices in `v` are filled.
The `by` keyword lets you provide a function that will be applied to
each element before comparison; the `lt` keyword allows providing a
custom "less than" function; use `rev=true` to reverse the sorting
order. These options are independent and can be used together in all
possible combinations: if both `by` and `lt` are specified, the `lt`
function is applied to the result of the `by` function; `rev=true`
reverses whatever ordering specified via the `by` and `lt` keywords.
"""
function mergesorted!(v::AbstractVector,
                      vl::AbstractVector,
                      vr::AbstractVector;
                      lt=isless,
                      by=identity,
                      rev::Bool=false,
                      order::Ordering=Forward)
    ordr = ord(lt,by,rev,order)
    mergesorted!(v, vl, vr, ordr)
end

"""
    mergesorted(vl, vr; lt=isless, by=identity, rev::Bool=false, order::Ordering=Forward)
Merge sorted vectors `vl` and `vr`. See [`mergesorted!`](@ref).
"""
function mergesorted(vl::AbstractVector,
                      vr::AbstractVector;
                      lt=isless,
                      by=identity,
                      rev::Bool=false,
                      order::Ordering=Forward)
    v = similar(promote_type(typeof(vl),typeof(vr)), length(vl)+length(vr))
    ordr = ord(lt,by,rev,order)
    mergesorted!(v, vl, vr, ordr)
end

end # module Continuum
