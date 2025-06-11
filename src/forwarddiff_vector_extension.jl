# ==================================================================================================
# vectorized derivative
# ==================================================================================================
function ForwardDiff.DerivativeConfig(f::F,
                          y::Vector{Vector{Y}},
                          x::X,
                          tag::T=ForwardDiff.Tag(f, X)) where {F,X <: Real,Y <: Real,T}
    duals = [similar(y_, ForwardDiff.Dual{T,Y,1}) for y_ in y]
    return ForwardDiff.DerivativeConfig{T,typeof(duals)}(duals)
end

@inline function ForwardDiff.derivative(f::F, x::Vector{R}) where {F, R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    one_ = one(eltype(x))
    zero_ = zero(one_)

    duals = fill(ForwardDiff.Dual{T}(zero_, one_), length(x))
    duals .+= x

    res = f(duals)
    return map(r -> ForwardDiff.extract_derivative(T, r), res)
end

@inline function ForwardDiff.derivative!(result::Vector, #{R},
    f::F, x::Vector{R}) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    one_ = one(eltype(x))
    zero_ = zero(one_)

    duals = fill(ForwardDiff.Dual{T}(zero_, one_), length(x))
    duals .+= x

    ydual = f(duals)
    return ForwardDiff.extract_derivative!(T, result, ydual)
end

# ==================================================================================================
# seeding for gradient and jacobian
# ==================================================================================================
function ForwardDiff.seed!(duals::Vector{Vector{ForwardDiff.Dual{T,V,N}}}, x,
               seeds::NTuple{N,ForwardDiff.Partials{N,V}}) where {T,V,N}
    for dual_inds in 1:N
        duals[dual_inds] .= ForwardDiff.Dual{T,V,N}.(view(x, dual_inds)[1], Ref(seeds[dual_inds])) # TODO: quite expensive
    end

    return duals
end

# ==================================================================================================
# vectorized gradient
# ==================================================================================================
function ForwardDiff.GradientConfig(f::F,
                        x::Vector{Vector{V}},
                        ::ForwardDiff.Chunk{N}=ForwardDiff.Chunk{length(x)}(),
                        ::T=ForwardDiff.Tag(f, V)) where {F,V,N,T}
    seeds = ForwardDiff.construct_seeds(ForwardDiff.Partials{N,V})
    duals = [similar(x_, ForwardDiff.Dual{T,V,N}) for x_ in x] # TODO: slow?
    return ForwardDiff.GradientConfig{T,V,N,typeof(duals)}(seeds, duals)
end

function ForwardDiff.gradient(f, x::Vector{<:Vector{<:Number}}, cfg::ForwardDiff.GradientConfig{T}=ForwardDiff.GradientConfig(f, x), ::Val{CHK}=Val{true}()) where {T,CHK}
    CHK && ForwardDiff.checktag(T, f, x[1])
    if ForwardDiff.chunksize(cfg) == length(x)
        return ForwardDiff.vector_mode_gradient(f, x, cfg)
    else
        return ForwardDiff.chunk_mode_gradient(f, x, cfg)
    end
end

function ForwardDiff.vector_mode_gradient(f::F, x::Vector{<:Vector{<:Number}}, cfg::ForwardDiff.GradientConfig{T}) where {T,F}
    ydual = ForwardDiff.vector_mode_dual_eval!(f, cfg, x)
    grads = [collect(y.partials) for y in ydual]

    return grads
end

# ==================================================================================================
# vectorized jacobian
# ==================================================================================================
function ForwardDiff.JacobianConfig(f::F,
                        x::Vector{Vector{V}},
                        ::ForwardDiff.Chunk{N}=ForwardDiff.Chunk{length(x)}(),
                        ::T=ForwardDiff.Tag(f, V)) where {F,V,N,T}
    seeds = ForwardDiff.construct_seeds(ForwardDiff.Partials{N,V})
    duals = [similar(x_, ForwardDiff.Dual{T,V,N}) for x_ in x]
    return ForwardDiff.JacobianConfig{T,V,N,typeof(duals)}(seeds, duals)
end

function ForwardDiff.jacobian(f, x::Vector{<:Vector{<:Number}}, cfg::ForwardDiff.JacobianConfig{T}=ForwardDiff.JacobianConfig(f, x), ::Val{CHK}=Val{true}()) where {T,CHK}
    CHK && ForwardDiff.checktag(T, f, x[1])
    if ForwardDiff.chunksize(cfg) == length(x)
        return ForwardDiff.vector_mode_jacobian(f, x, cfg)
    else
        return ForwardDiff.chunk_mode_jacobian(f, x, cfg)
    end
end

function ForwardDiff.vector_mode_jacobian(f::F, x::Vector{<:Vector{<:Number}}, cfg::ForwardDiff.JacobianConfig{T}) where {T,F}
    N = ForwardDiff.chunksize(cfg)
    ydual = ForwardDiff.vector_mode_dual_eval!(f, cfg, x)
    result = [similar(ydual[1], valtype(eltype(ydual[1])), length(ydual[1]), N) for _ in 1:length(ydual)]
    for i in 1:length(result)   ForwardDiff.extract_jacobian!(T, result[i], ydual[i], N)   end
    for i in 1:length(result)   ForwardDiff.extract_value!(T, result[i], ydual[i])   end
    return result
end
