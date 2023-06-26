@recipe function plot(r::ResultGap; x=0:0.01:2)
    seriestype := :path
    line := 2
    xlabel --> "q [GeV]"
    ylabel --> "Dynamical Mass [GeV]"
    label --> ""
    q, mass = x, r.dynamicmass.(x)
end
