module Continuum

 using Statistics
 using DataInterpolations
 using SortFilters, RollingFunctions, NaNMath
 #https://github.com/sairus7/SortFilters.jl

# constants
speed_of_light_mks = 3e8 # m/s
fwhm_sol = 7.3e3 # m/s

 function calc_rolling_max(f::AV2; width::Integer = 10 ) where {  T2<:Real, AV2<:AbstractVector{T2} }
 shift = floor(Int64,width//2)
 z = similar(f)
 #z[(shift):(length(f)-shift)] .= rollmax(f,width)
 z[(shift):(length(f)-shift)] .= rolling(NaNMath.maximum,f,width)
 z[1:shift-1] .= z[shift]
 z[(length(f)-shift+1):end] .= z[length(f)-shift]
 return z
end

function calc_rolling_median_and_max_to_clip(λ::AV1, f::AV2; width::Integer = 10, cr_pad_width::Integer = 1,
          verbose::Bool = false ) where { T1<:Real, T2<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2} }
  @assert length(λ) == length(f)
  @assert width >= 4
  @assert width%2 == 0
  @assert 0<= cr_pad_width <= 7
  shift = floor(Int64,width//2)
  p = (0.25, 0.5, 0.75, 0.95)
  function iqr(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} }
      (mq[i][3]-mq[i][1])    # standard iqr
  end
  function iqr_upper(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} }
      # moving_quartiles[i][3]-moving_quartiles[i][1]    # standard iqr
      2*(mq[i][3]-mq[i][2])  # double Q3-Q2
  end
  #q2(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} } = mq[i][2]
  q75(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} } = mq[i][3]
  q95(mq::AVNT, i::Integer) where { T<:Real, NT<:NTuple{4,T}, AVNT<:AbstractVector{NT} } = mq[i][4]

  output = similar(λ)
  mask = trues(length(λ))
  f_median_filtered = similar(f)
  f_clean = deepcopy(f)
  num_clippings = 5
  flush(stdout)
  for i in 1:num_clippings
     moving_quartiles = movsort(f_clean, width, p)
     f_median_filtered[(1+shift):(length(f)-shift)] .= map(i->moving_quartiles[i][2], (1+2*shift):length(moving_quartiles) )
     f_median_filtered[1:shift] .= f_median_filtered[1+shift]
     f_median_filtered[(length(f)-shift+1):end] .= f_median_filtered[length(f)-shift]
     max_iqr = map(i->q75(moving_quartiles,i)+1.5*iqr_upper(moving_quartiles,i),1:length(λ))
     idx = (1+2*shift):length(moving_quartiles)
     output[(1+shift):(length(f)-shift)] .= max_iqr[idx]
     #=
     # Pass two with simple high quantile with wide window (to avoid being affected by lines)
     width = 100  # TODO: Automate for orders.  FOr now arbtrary rather than 100A, while testing on small chunks
     shift = floor(Int64,width//2)
     max_99 = movsort(f, width, 0.98)  # TODO Replace with 0.99 after expand window
     output[(1+shift):(length(f)-shift)] .= max.(output[idx]),max_99[idx])
     =#
     output[(1+shift):(length(f)-shift)] .= max.(output[idx],map(i->q95(moving_quartiles,i),idx))

     # Fill in edges
     output[1:shift] .= output[1+shift]
     output[(length(f)-shift+1):end] .= output[length(f)-shift]
     # Identify & replace high outliers
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
  return (f_median_filtered, output, f_clean)
end


function calc_rollingpin_radius(λ::AV1, f::AV2; fwhm::Real = fwhm_sol,
            S1_width_factor::Real = 40, S2_width_factor::Real = 10, min_R_factor::Real = 10, max_R_factor::Real = 15, ν::Real = 1 ) where { T1<:Real, AV1<:AbstractVector{T1}, T2<:Real, AV2<:AbstractVector{T2}  }
   λlo = minimum(λ)
   S1_width = S1_width_factor*fwhm/speed_of_light_mks*λlo
   S2_width = S2_width_factor*S1_width
   S1_width = 2*floor(Int64,S1_width/2)
   S2_width = 2*floor(Int64,S2_width/2)
   f_max1 = calc_rolling_max(f,width=S1_width)
   f_max2 = calc_rolling_max(f,width=S2_width)
   #penalty = (f_max2.-f_max1)./f_max2
   #conversion_fwhm_sig = (10*min_λ/(sqrt(8*log(2))*speed_of_light_kmps))
   #min_radius = fwhm * conversion_fwhm_sig # par_R
   #radius = min_radius .* λ./min_λ
   rollingpin_R_min = min_R_factor*λlo*(fwhm/speed_of_light_mks)/log(8*log(2))
   rollingpin_R_max = max_R_factor*rollingpin_R_min  # Arbitrary, not same as Rassine's auto setting, see https://github.com/MichaelCretignier/Rassine_public/blob/master/Rassine.py#L718
   #rollingpin_R_min = 2*floor(Int64,rollingpin_R_min/2)
   #rollingpin_R_max = 2*floor(Int64,rollingpin_R_max/2)
   #rollingpin_ν = 1.0  # TODO Experiment for NEID data
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
  function calc_dist(i::Integer, j::Integer)  # TODO OPT: reliminate sqrt
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
  #=
  interp_cubic = CubicSpline(f[keep],λ[keep])
  interp_linear = LinearInterpolation(f[keep],λ[keep])
  function extrap(x::Real)
    if x<first(interp_linear.t)
      return first(interp_linear.u)
    elseif x>last(interp_linear.t)
      return last(interp_linear.u)
    else
      return interp_linear(x)
    end
  end
  =#
  #  return (anchors = keep, continuum=interp_cubic.(λ))
  return keep
  #return (anchors = keep, continuum=extrap.(λ))
end

function calc_continuum_from_anchors( λ::AV1, f::AV2, anchors::AV3,
              verbose::Bool = false ) where {
                T1<:Real, T2<:Real, T3<:Integer, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3} }
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
  #  return (anchors = keep, continuum=interp_cubic.(λ))
  return extrap.(λ)
end


function replace_edge_anchor_vals!(f::AV1, n::Integer = 3 ) where { T1<:Real, AV1<:AbstractVector{T1} }
  @assert length(f) >= 2*n+1
  f[1:(n+1)] .= f[n+1]
  f[end-n:end] .= f[length(f)-n-1]
  return f
end

function find_clean_anchors_by_slope(anchors::AV1, f::AV2; threshold::Real = 0.95, verbose::Bool = false ) where { T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3} }
  nanchors = length(anchors)
  @assert nanchors == length(f)
  @assert nanchors >= 8
  Δabsslope = zeros(length(anchors))
  for i in 2:nanchors-1
    slope_hi = (f[i+1]-f[i])/(anchors[i+1]-anchors[i])
    slope_lo = (f[i]-f[i-1])/(anchors[i]-anchors[i-1])
    Δabsslope[i] = abs(slope_hi)-abs(slope_lo)
  end
  threshold = quantile(Δabsslope[2:end-1],0.995)  # TODO: Update to 99.5% once have large enough dataset
  mask = Δabsslope.<threshold
  return mask
end

function merge_nearby_anchors(anchors::AV1, λ::AV2; threshold::Real = 1, verbose::Bool = false ) where { T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3} }
  Δλ = λ[anchors[2:end]] .- λ[anchors[1:end-1]]
  close_anchor_idx = (1:(length(anchors)-1))[Δλ .<= threshold]
  nclose = length(close_anchor_idx)
  if verbose println("# Found ", nclose, " close anchors.")  end
  if nclose == 0
    return anchors
  end
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
      merged_anchor = anchors[close_anchor_idx[merged_anchor_idx]]
      push!(new_anchors,merged_anchor)
      idx_start = i
    end
  end
  if idx_start != nclose
    idx_stop = nclose
    #println("# start = ", idx_start, " stop = ", idx_stop, " anchors_idx = ", close_anchor_idx[idx_start:idx_stop])
    if idx_start == idx_stop
      merged_anchor_idx = idx_start
      #merged_anchor = close_anchor_idx[idx_start]
    elseif (idx_stop-idx_stop)%2 == 0
      merged_anchor_idx = idx_start + floor(Int64,(idx_stop-idx_start)//2)
      #merged_anchor = anchors[close_anchor_idx[merged_anchor_idx]]
      #println("# 1 start = ", idx_start, " stop = ", idx_stop, " mid = ", merged_anchor_idx)
      #println("# 2 start = ", close_anchor_idx[idx_start], " stop = ", close_anchor_idx[idx_stop], " mid = ",  close_anchor_idx[merged_anchor_idx] )
      #println("# 3 start = ", λ[close_anchor_idx[idx_start]], " stop = ", λ[close_anchor_idx[idx_stop]], " mid = ",  λ[close_anchor_idx[merged_anchor_idx]] )
    else
      merged_anchor_idx = idx_start + floor(Int64,(idx_stop-idx_start+1)//2)
      #merged_anchor = anchors[close_anchor_idx[merged_anchor_idx]]
      #println("# 1 start = ", idx_start, " stop = ", idx_stop, " mid = ", merged_anchor_idx)
      #println("# 2 start = ", close_anchor_idx[idx_start], " stop = ", close_anchor_idx[idx_stop], " mid = ",  close_anchor_idx[merged_anchor_idx] )
      #println("# 3 start = ", λ[close_anchor_idx[idx_start]], " stop = ", λ[close_anchor_idx[idx_stop]], " mid = ",  λ[close_anchor_idx[merged_anchor_idx]] )
    end
    merged_anchor = anchors[close_anchor_idx[merged_anchor_idx]]
    push!(new_anchors,merged_anchor)
  end
  return new_anchors
end

function calc_continuum(λ::AV1, f_obs::AV2; fwhm::Real = fwhm_sol ) where { T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3} }
 clip_width_A = 1.0
 clip_width_pix = 2*floor(Int64,clip_width_A/(λ[2]-λ[1])/2)
 (f_median_filtered, f_clip_threshold, f_clean) = calc_rolling_median_and_max_to_clip(λ,f_obs, width=clip_width_pix)
 rollingpin_radius = calc_rollingpin_radius(λ, f_median_filtered, fwhm=fwhm)
 anch_orig = calc_continuum_anchors(λ,f_median_filtered,radius=rollingpin_radius, stretch_factor=4.0)
 anchor_vals = f_median_filtered[anch_orig]
 replace_edge_anchor_vals!(anchor_vals)
 anch_mask = find_clean_anchors_by_slope(anch_orig,anchor_vals)
 anchors_merged = merge_nearby_anchors(anch_orig[anch_mask],λ, threshold=0.04)
 continuum = calc_continuum_from_anchors(λ,f_median_filtered,anchors_merged)
 return (anchors_merged, continuum, f_median_filtered)
end

end # module Continuum