ker=kernel11-I
right=(-1).*ones(dim)
root=ker\right
E1=Array{Float64}(undef,kstep,zstep)
Threads.@threads for i=1:dim
    # f=Int((i-1)÷dim+1)
    kk=getk((i-1)%dim+1)
    zk=getz((i-1)%dim+1)
    E1[kk,zk]=root[i]
end
print(E1)