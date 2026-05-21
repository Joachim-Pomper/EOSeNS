module EOSeNS

include("utils/utils.jl")
include("equation_of_state/neos.jl")
include("tov_equation/tov.jl")

using .Utils
using .NEOS
using .TOV

export Utils, NEOS, TOV

for name in names(NEOS; all=false)
    if Base.isexported(NEOS, name)
        @eval export $name
    end
end

for name in names(TOV; all=false)
    if Base.isexported(TOV, name)
        @eval export $name
    end
end

end
