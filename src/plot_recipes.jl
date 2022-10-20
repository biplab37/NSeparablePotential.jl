@recipe function plot(r::ResultGap)
    seriestype := :path
    line := 2
    xlabel -> "T"
    ylabel -> "massgap"
    Tlist, massgap = r.Tlist, r.massgap
end