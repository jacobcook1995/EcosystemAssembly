# A script to run Lyapunov analysis of our inhibition model.
using Assembly
using JLD
using SymPy

# function to set up and run a lyapunov analysis of a simple case
function lya()
    # only a small number of strains and metabolites for this analysis
    N = 4
    M = 4
    O = 4
    mR = 2.0
    sdR = 0.0
    mq = 1.0
    sdq = 0.1
    mK  = 0.1
    sdK = 0.01
    mk = 10.0
    sdk = 1.0
    # Now make the parameter set, can just use the normal function for now
    ps = initialise(N,M,O,mR,sdR,mq,sdq,mK,sdK,mk,sdk)
    # save this parameter set
    jldopen("Temp/Paras/psL.jld","w") do file
        write(file,"ps",ps)
    end
    # Preallocate Jacobian
    J = Array{Sym,2}(undef,ps.N+ps.M,ps.N+ps.M)
    # Find Jacobian using function
    @time J1 = Jacobian(ps,J)
    @time J2 = Jacobian_test(ps)
    # Now want to define partciluar concentrtion and population values to use
    N1 = 10.0
    N2 = 27.0
    N3 = 12.0
    N4 = 123.1
    M1 = 3.4
    M2 = 7.8
    M3 = 3.4
    M4 = 2.1
    for i = 1:ps.N+ps.M
        for j = 1:ps.N+ps.M
            J1[i,j] = subs(J1[i,j],"N1"=>N1,"N2"=>N2,"M1"=>M1,"M2"=>M2)
            J1[i,j] = subs(J1[i,j],"N3"=>N3,"N4"=>N4,"M3"=>M3,"M4"=>M4)
            J2[i,j] = subs(J2[i,j],"N1"=>N1,"N2"=>N2,"M1"=>M1,"M2"=>M2)
            J2[i,j] = subs(J2[i,j],"N3"=>N3,"N4"=>N4,"M3"=>M3,"M4"=>M4)
        end
    end
    println(J1)
    println(J2)
    return(nothing)
end

@time lya()
