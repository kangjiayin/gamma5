using Plots
test=Array{Float64}(undef,kstep,1)
for i=1:kstep
    test[i]=F1[i,1]
end
    
# using JSON3
scatter(meshk,test)
# open("/Users/kjy/Desktop/program/julia/Gamma5/F1.json","w") do f
#     JSON3.write(f, F1)
#     println(f)
# end