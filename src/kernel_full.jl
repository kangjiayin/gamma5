using Plots
kernel=[kernel11 kernel12 kernel13 kernel14;kernel21 kernel22 kernel23 kernel24;kernel31 kernel32 kernel33 kernel34;kernel41 kernel42 kernel43 kernel44;]
#δ=Diagonal(ones(4*dim))
kernel_eva=kernel-I
right=[ones(dim) ;zeros(3*dim)]


solution=-kernel_eva\right
F1=Array{Float64}(undef, kstep, zstep)
F2=Array{Float64}(undef, kstep, zstep)
F3=Array{Float64}(undef, kstep, zstep)
F4=Array{Float64}(undef, kstep, zstep)
FF=Array{Array}(undef,4,1)
FF[1]=F1
FF[2]=F2
FF[3]=F3
FF[4]=F4

Threads.@threads for i=1:4*dim
    f=Int((i-1)÷dim+1)
    k=getk((i-1)%dim+1)
    z=getz((i-1)%dim+1)
    FF[f][k,z]=solution[i]
end

F1k=zeros(kstep)
for i=1:kstep
    for j=1:zstep
        F1k[i]+=F1[i,j]*weightz[j]
    end
end
# plot(meshk,F1k,scale=:log10)