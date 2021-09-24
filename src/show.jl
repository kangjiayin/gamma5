using Plots
for i=1:32
    if i==1
        scatter(P2,F1k[1,:],xscale=:log10)
    else
        scatter!(P2,F1k[i,:],xscale=:log10)
    end
end
