# This is a simple set of datastructures and functions to work with fractional differential equations
# as sums of differential terms of arbitrary fractional order.

# A Kabla, 2021, 2024.



# Term of a differential equation, expressed as:
#     coef * d^order / dt^order
# Each parameter could have a set value, or be a symbol or expression to evaluate later.

struct DETerm
    coef::Union{Real, Expr, Symbol}
    order::Union{Real, Expr, Symbol}
end


# Struct to hold a list of terms for the left and right handside of the equation

struct DiffEqu
    leftvar::Symbol     # e.g. ϵ for rheology
    rightvar::Symbol    # and σ on the other side
    leftde::Set{DETerm}
    rightde::Set{DETerm}
end


#  Basic constructor of an empty equation

function DiffEqu(leftvar::Symbol, rightvar::Symbol)
    return(DiffEqu(leftvar, rightvar, Set{DETerm}(), Set{DETerm}()))
end

# To add terms to the diff equation

function push!(de::DiffEqu, var::Symbol, t::Tuple{Union{Real, Expr, Symbol},Union{Real, Expr, Symbol}})
    if var==de.leftvar
        Base.push!(de.leftde, DETerm(t...))
    elseif var==de.rightvar
        Base.push!(de.rightde, DETerm(t...))
    end
    return(nothing)
end

function push!(de::DiffEqu, var::Symbol, tt::Tuple{Vararg{Tuple}})
    for t in tt
        push!(de,var,t)
    end
    return(nothing)
end


#  Contructor for full fractional diff equation

function DiffEqu(;kwargs...)
    vars = keys(kwargs)
    de = DiffEqu(vars...)#  , values(kwargs)...)
    push!(de, vars[1], kwargs[vars[1]])
    push!(de, vars[2], kwargs[vars[2]])
    return(de)
end







# Returns a suitable derivative expression

function diffexpr(var::Symbol, order)
    if order==0
        return(var)
    elseif order==1  #(integer)
        dvar=Symbol("d"*string(var)*"dt")
        return(dvar)
    elseif typeof(order)<:Real && floor(order)==order  #(integer value)
        dvar=Symbol("d$order"*string(var)*"dt$order")
        return(:($dvar,))
    else
        return(:(dᵅdtᵅ($var, $order)))   #  May need time data too.
    end
end

# Returns a laplace/frequency expression

function laplaceexpr(order, var=:s)
    if order==0
        return(1)
    elseif order==1  #(integer)
        return(var)
    else 
        return(:($var^$order))
    end
end




#
#   Placeholder function - does nothing at the moment but will eventually return a derivative. May need time data too.
#
function dᵅ_dtᵅ(vat, order)
    return(order)
end


# Tools to operate on expressions

function addexpr(ex, tex)
    return( :($ex + $tex) )
end

function multexpr(ex, tex)
    if ex==1
        return(tex)
    elseif tex==1
        return(ex)
    else
        return( :($ex * $tex) )
    end
end


#
# Builing the differential equation
#

function builddeexpr(var::Symbol, de::Set{DETerm})
    ex = :()
    if length(de)==1
        t = first(de)
        ex = multexpr(t.coef, diffexpr(var, t.order))
    else
        s = copy(de)
        t = pop!(s)
        ex = multexpr(t.coef, diffexpr(var, t.order))
        for t in s
            ex = addexpr(ex, multexpr(t.coef, diffexpr(var, t.order)) )
        end
    end
    return(ex)
end

function builddeexpr(de::DiffEqu)
    leftex = builddeexpr(de.leftvar, de.leftde)
    rightex = builddeexpr(de.rightvar, de.rightde)
    return( :($leftex - $rightex) )
end



#
# Building the laplace / harmonic responses
#

function buildlaplaceexpr(de::Set{DETerm}, var=:s)
    ex = :()
    if length(de)==1
        t = first(de)
        ex = multexpr(t.coef, laplaceexpr(t.order, var))
    else
        s = copy(de)
        t = pop!(s)
        ex = multexpr(t.coef, laplaceexpr(t.order, var))
        for t in s
            ex = addexpr(ex, multexpr(t.coef, laplaceexpr(t.order, var)) )
        end
    end
    return(ex)
end


function builddynamicexpr(de::DiffEqu)
    leftex = buildlaplaceexpr(de.leftde, :(im * ω))
    rightex = buildlaplaceexpr(de.rightde, :(im * ω))
    return( :($leftex / $rightex) )
end


function buildrelaxationexpr(de::DiffEqu)
    leftex = buildlaplaceexpr(de.leftde)
    rightex = buildlaplaceexpr(de.rightde)
    return( :($leftex / (s * $rightex) ) )
end


function buildcreepexpr(de::DiffEqu)
    leftex = buildlaplaceexpr(de.leftde)
    rightex = buildlaplaceexpr(de.rightde)
    return( :($rightex / (s * $leftex) ) )
end





#
#   Some demo code
#

println("=====================")
println("Maxwell model")
println("=====================")


demax = DiffEqu(ϵ = (:η,1), σ = ((:k,0), (:η,1)) )

println("Differential equation (cost function)")
println(builddeexpr(demax))
println("Relaxation modulus")
println(buildrelaxationexpr(demax))
println("Creep compliance")
println(buildcreepexpr(demax))
println("Dynamic modulus")
println(builddynamicexpr(demax))


println("=====================")
println("Fract Zener model")
println("=====================")


defzener = DiffEqu( ϵ = ( (:cₐ,:α), (:cᵧ,:γ), (:(cₐ*cᵧ/cᵦ),:(α+γ-β)) ), 
                    σ = ( (1,0), (:(cₐ/cᵦ),:(α-β)))  )


println("Differential equation (cost function)")
println(builddeexpr(defzener))
println("Relaxation modulus")
println(buildrelaxationexpr(defzener))
println("Creep compliance")
println(buildcreepexpr(defzener))
println("Dynamic modulus")
println(builddynamicexpr(defzener))


