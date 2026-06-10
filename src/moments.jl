using Base.Iterators

"""
    trace_product(o1, o2; scale=0)
"""
function trace_product(
    o1::Operator, o2::Operator; 
    scale=0)
    if length(o1.strings) < length(
        o2.strings)
        return trace_product(
            o2, o1; scale)
    end

    checklength(o1, o2)
    N = qubitlength(o1)
    tr = zero(scalartype(o1))

    if length(o1.strings) != length(
        o1.coeffs)
        throw(DimensionMismatch(
            "Length mismatch in o1"))
    end
    if length(o2.strings) != length(
        o2.coeffs)
        throw(DimensionMismatch(
            "Length mismatch in o2"))
    end

    d = emptydict(o2)
    @inbounds for i in eachindex(
        o2.strings)
        insert!(
            d, o2.strings[i], 
            o2.coeffs[i])
    end

    @inbounds for i in eachindex(
        o1.strings)
        p1 = o1.strings[i]
        c1 = o1.coeffs[i]
        c2 = get(d, p1, nothing)
        isnothing(c2) && continue
        p, k = prod(p1, p1)
        tr += c1 * c2 * k
    end

    (scale == 0) && (scale = 2.0^N)
    return tr * scale
end

function trace_product(
    o1::Operator{<:PauliStringTS}, 
    o2::Operator{<:PauliStringTS}; 
    scale=0)
    checklength(o1, o2)
    Ls = qubitsize(o1)
    Ps = periodicflags(o1)
    tr = zero(scalartype(o1))

    d = emptydict(o2)
    for (p2, c2) in zip(
        o2.strings, o2.coeffs)
        insert!(d, p2, c2)
    end

    for (p1, c1) in zip(
        o1.strings, o1.coeffs)
        c2 = get(d, p1, nothing)
        isnothing(c2) && continue
        rep1 = representative(p1)
        p, k = prod(rep1, rep1)
        f = c1 * c2 * k
        for s in all_shifts(Ls, Ps)
            shifted = shift(
                rep1, Ls, Ps, s)
            if shifted == rep1
                tr += f
            end
        end
    end
    if iszero(scale)
        scale = 2.0^Base.prod(Ls)
    end
    num_translations = Base.prod(
        L for (L, p) in zip(Ls, Ps) if p)
    return tr * scale * num_translations
end

Base.@deprecate oppow(
    o::AbstractOperator, k::Int) o^k

Base.:^(
    o::AbstractOperator, 
    k::Int) = Base.power_by_squaring(o, k)

"""
    trace_product(A, k, B, l; scale=0)
"""
function trace_product(
    A::AbstractOperator, k::Int, 
    B::AbstractOperator, l::Int; 
    scale=0)
    @assert typeof(A) == typeof(B)
    m = div(k + l, 2)
    n = k + l - m
    if k < m
        C = A^k * B^(m - k)
        D = B^n
    elseif k > m
        C = A^m
        D = A^(k - m) * B^l
    else
        C = A^k
        D = B^l
    end
    return trace_product(
        C, D; scale=scale)
end

"""
    trace_product(A; scale=0)
"""
function trace_product(
    A::Operator; scale=0)
    c = get_coeffs(A)
    N = qubitlength(A)
    val = iszero(scale) ? 2.0^N : scale
    return sum(c.^2) * val
end

trace_product(
    A::Operator{<:PauliStringTS}; 
    scale=0) = trace_product(
        A, A; scale=scale)

"""
    trace_product_z(o1, o2; scale=0)
"""
function trace_product_z(
    o1::AbstractOperator, 
    o2::AbstractOperator; scale=0)
    if iszero(scale)
        scale = 2.0^qubitlength(o1)
    end
    tr = zero(scalartype(o1))

    for i in eachindex(o1.strings)
        p1 = o1.strings[i]
        c1 = o1.coeffs[i]
        for j in eachindex(o2.strings)
            p2 = o2.strings[j]
            c2 = o2.coeffs[j]

            p, k = prod(p1, p2)
            if xcount(p) == ycount(p) == 0
                tr += c1 * c2 * k
            end
        end
    end

    return tr * scale
end

"""
    moments(H, kmax; start=1, scale=0)
"""
function moments(
    H::AbstractOperator, kmax::Int; 
    start=1, scale=0)
    return [trace_product(
        H, k; scale=scale) for 
        k in start:kmax]
end

function trace_product(
    o::Operator, p::PauliString; 
    scale=0)
    checklength(o, p)
    c = get_coeff(o, p)
    N = qubitlength(o)
    (scale == 0) && (scale = 2.0^N)
    return c * scale
end

trace_product(
    p::PauliString, o::Operator; 
    scale=0) = trace_product(
        o, p; scale=scale)

function trace_product(
    o1::Operator{<:PauliStringTS}, 
    o2::PauliStringTS; scale=0)
    checklength(o1, o2)
    Ls = qubitsize(o1)
    Ps = periodicflags(o1)
    tr = zero(scalartype(o1))
    i = findfirst(==(o2), o1.strings)
    isnothing(i) && return tr
    rep1 = representative(o2)
    p, k = prod(rep1, rep1)
    c1 = o1.coeffs[i]
    c2 = (1im)^ycount(o2)
    f = c1 * c2 * k
    for s in all_shifts(Ls, Ps)
        shifted = shift(
            rep1, Ls, Ps, s)
        if shifted == rep1
            tr += f
        end
    end
    if iszero(scale)
        scale = 2.0^Base.prod(Ls)
    end
    num_translations = Base.prod(
        L for (L, p) in zip(Ls, Ps) if p)
    return tr * scale * num_translations
end

trace_product(
    p::PauliStringTS, 
    o::Operator{<:PauliStringTS}; 
    scale=0) = trace_product(
        o, p; scale=scale)

function trace_product(
    s1::P, s2::P; 
    scale=0) where {P<:PauliString}
    N = qubitlength(s1)
    if s1 == s2
        (iszero(scale)) && (scale = 2.0^N)
        return scale
    else
        return 0
    end
end

function trace_product(
    s1::P, s2::P; 
    scale=0) where {P<:PauliStringTS}
    N = qubitlength(s1)
    if s1 == s2
        (scale == 0) && (scale = 2.0^N)
        Ls = qubitsize(s1)
        Ps = periodicflags(s1)
        num_translations = Base.prod(
            L for (L, p) in zip(
                Ls, Ps) if p)
        return scale * num_translations
    else
        return 0
    end
end

# === OUR OPTIMIZED FAST-PATH ===

# Check for anticommuting core
function is_anticommuting_core(
    A::AbstractOperator)
    if typeof(A) <: Operator{<:PauliStringTS}
        return false
    end
    if !(typeof(A) <: Operator)
        return false
    end

    strs = A.strings
    len = length(strs)
    if len < 2 || len > 3
        return false
    end
    
    if len == 2
        return commutator(
            strs[1], strs[2])[2] != 0
    else
        return commutator(
            strs[1], strs[2])[2] != 0 &&
               commutator(
            strs[1], strs[3])[2] != 0 &&
               commutator(
            strs[2], strs[3])[2] != 0
    end
end

# Analytical lookup
function analytic_trace(
    A::AbstractOperator, k::Int)
    if k % 2 != 0
        return zero(scalartype(A))
    end
    
    c = get_coeffs(A)
    base_val = sum(c.^2)
    half_k = div(k, 2)
    
    if half_k % 2 == 0
        return base_val^(half_k)
    else
        return -base_val^(half_k)
    end
end

# Main optimized trace function
function trace_product(
    A::AbstractOperator, 
    k::Int; scale=0)
    N = qubitlength(A)
    current_scale = iszero(scale) ? 
                    2.0^N : scale
    
    if is_anticommuting_core(A)
        return analytic_trace(
            A, k) * current_scale
    end
    
    m = div(k, 2)
    n = k - m
    C = A^m
    if k % 2 == 0
        return trace_product(
            C; scale=scale)
    end
    D = A^n
    return trace_product(
        C, D; scale=scale)
end
EOF
