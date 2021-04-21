module Continuum

using RCall
using RvSpectMLBase
using EchelleInstruments

function calc_continuum_model_order(spectrum::AbstractSpectra2D; order_idx::Integer )
    possible_pix = get_inst_module(get_inst(spectrum)).get_pixel_range(get_inst(spectrum),order_idx)
    bad_pix = get_inst_module(get_inst(spectrum)).bad_col_ranges(get_inst(spectrum),order_idx)
    pix_rng = EchelleInstruments.calc_complement_index_ranges(possible_pix,bad_pix)
    pix = mapreduce(p->collect(p),vcat,pix_rng)
    afs_inputs = zeros(Float64,length(pix),2)
    afs_inputs[:,1] .= spectrum.λ[pix,order_idx]
    afs_inputs[:,2] .= spectrum.flux[pix,order_idx]
    @assert !any(isnan.(afs_inputs))
    #=
    wv = mapreduce(p->spec.λ[p,order_idx],vcat,pix_rng)
    @assert !any(isnan.(wv))
    inten = mapreduce(p->convert(Vector{Float64},spec.flux[p,order_idx]),vcat,pix_rng)
    @assert !any(isnan.(inten))
    afs_inputs = hcat(wv,inten)
    =#
    #df = DataFrame("wv"=>wv,"intes"=>inten)
    afs_src = joinpath("src","AFS.R")
    R"source($afs_src)"
    afs_output_R = R"AFS($afs_inputs,0.95,0.25)"
    afs_output = rcopy(afs_output_R) 
    continuum = zeros(eltype(spectrum.flux),size(spectrum.flux[:,order_idx]))
    continuum = fill(NaN,size(spectrum.flux[:,order_idx]))
    continuum[pix] .= afs_output
    return continuum
end

function calc_continuum_model(spectrum::AbstractSpectra2D )
    vec_of_orders = map(ord->calc_continuum_model_order(spectrum,order_idx=ord), min_order(get_inst(spectrum)):max_order(get_inst(spectrum)) )
    output = fill(NaN, size(spectrum.flux))
    for (i,ord) in enumerate(min_order(get_inst(spectrum)):max_order(get_inst(spectrum)))
        output[:,ord] .= vec_of_orders[i]
    end
    return output
end

end
