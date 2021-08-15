
# a=ones(100)
# b=0.
# @btime for i=1:10
#     global a,b
#     aa=a
#     bbb=aa[i]
#     b=bbb+b
# end
function GaussChebyshevIntegral2(f,n)
    x,w=gausschebyshev(n::Int64,2)
    out=dot(w,[f(x[i]) for i=1:n])
    out::Float64
end
