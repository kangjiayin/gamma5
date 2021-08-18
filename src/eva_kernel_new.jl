# using ProfileVega
# include("./integral.jl")

# function upgrateconst(i,j)
#     QPlus2=qPlus2[j]
#     QSub2=qSubt2[j]
#     Kdotp=kdotp[i]
#     Pdotq=pdotq[j]
#     innerA1=A1[j]
#     innerB1=B1[j]
#     innerA2=A2[j]
#     innerB2=B2[j]
#     weightzk=weightz[i]
#     weightzq=weightz[j]
#     K2=k[i]
#     Q2=k[j]
#     Zk=z[i]
#     Zq=z[j]
# end

@time Threads.@threads for i = 1:dim
    for j=1:dim
        K2=k[i]::Float64
        Q2=k[j]::Float64
        Zk=z[i]::Float64
        Zq=z[j]::Float64
        # QPlus2=qPlus2[j]
        # QSub2=qSubt2[j]
        Kdotp=kdotp[i]::Float64
        Pdotq=pdotq[j]::Float64
        innerA1=A1[j]::Float64
        innerB1=B1[j]::Float64
        innerA2=A2[j]::Float64
        innerB2=B2[j]::Float64
        #weightzk=weightz[i]
        Weightzq=weightz[Int(getz(Int(j)))]::Float64
        Weightq=weightk[Int(getz(Int(j)))]::Float64
        #这里的a,b,c是纯纯的过程参数，可以不用管
        a=sqrt(K2*Q2)::Float64
        b=Zk*Zq::Float64
        c=sqrt((1-Zk^2)*(1-Zq^2))::Float64
        b=a*b::Float64
        c=a*c::Float64
        Kdotq(y)=b+c*y::Float64
        d=K2+Q2::Float64
        Ksubq2(y)=d-2*Kdotq(y)::Float64
        allkindsofweight=Weightzq*Weightq*Q2*1/(16*pi^3)::Float64
        kernel11[i,j]=GaussChebyshevIntegral64(y->
        -allkindsofweight*
        D(Ksubq2(y))/branch[Int(j)]*
        (-4*innerB1*innerB2 + innerA1*innerA2*(P2 - 4*Q2))
        )::Float64

        # kernel12[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # (-2*P2*(innerA2*innerB1 + innerA1*innerB2)+4*(innerA2*innerB1 - innerA1*innerB2)*Pdotq)
        # )::Float64

        # kernel13[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # (-2*Kdotp*((innerA2*innerB1 + innerA1*innerB2)*Kdotp + 2*(-innerA2*innerB1 + innerA1*innerB2)*Kdotq(y)))
        # )::Float64

        # kernel14[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # (4*innerA1*innerA2*(P2*Kdotq(y) - Kdotp*Pdotq))
        # )::Float64

        # kernel21[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((2*(K2^2*(P2*(innerA2*innerB1 + innerA1*innerB2) + 2*(-innerA2 *innerB1 + innerA1* innerB2)* Pdotq) + K2 *(-((innerA2* innerB1 + innerA1* innerB2) *Kdotp^2) + 2*(innerA2* innerB1 + innerA1* innerB2)* Pdotq^2 + 2*Kdotq(y) *(-P2 *(innerA2* innerB1 + innerA1 *innerB2) + 4 *(innerA2* innerB1 - innerA1* innerB2) *Pdotq) + Q2* (P2 *(innerA2* innerB1 + innerA1* innerB2) + 6 *(-innerA2* innerB1 + innerA1 *innerB2)* Pdotq) + 2*Kdotp*((innerA2* innerB1 - innerA1* innerB2) *Kdotq(y) - (innerA2 *innerB1 + innerA1 *innerB2)* Pdotq)) + Kdotp *((innerA2* innerB1 + innerA1* innerB2)* Kdotp *(-Q2 + 4*Kdotq(y)) - 2*Kdotq(y)*(3 *(-innerA2 *innerB1 + innerA1* innerB2) *Q2 + 4* (innerA2* innerB1 - innerA1 *innerB2) *Kdotq(y) + (innerA2 *innerB1 + innerA1* innerB2) *Pdotq))))/(3 *(-P2*K2 + Kdotp^2)*(K2 + Q2 - 2*Kdotq(y))))
        # )::Float64

        # kernel22[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((Kdotp *((-4 *innerB1* innerB2 + innerA1* innerA2* (P2 + 4* Q2))*Kdotp *(-Q2 + 4*Kdotq(y)) - 2 *Kdotq(y)* (-4* innerB1* innerB2 + innerA1* innerA2* (P2 - 8 *Q2 + 16 *Kdotq(y)))* Pdotq) + K2^2 *(P2^2 *innerA1* innerA2 + P2 *(-4 *innerB1 *innerB2 + 4* innerA1* innerA2* Q2) - 8* innerA1* innerA2 *Pdotq^2) + K2 *(4 *P2 *innerA1 *innerA2 *Q2^2 + (4 *innerB1 *innerB2 - innerA1* innerA2 *(P2 + 4*Q2)) *Kdotp^2 - 2* P2* (P2* innerA1* innerA2 - 4 *innerB1* innerB2)* Kdotq(y) + 2 *Kdotp *(4* innerB1* innerB2 + innerA1* innerA2*(-P2 - 4* Q2 + 4* Kdotq(y)))* Pdotq + 2* (-4* innerB1* innerB2 + innerA1* innerA2* (P2 + 16 *Kdotq(y)))* Pdotq^2 + Q2 *(P2 *(-4* innerB1* innerB2 + innerA1 *innerA2 *(P2 - 8* Kdotq(y))) - 16 *innerA1* innerA2* Pdotq^2)))/(3 *(-P2* K2 + Kdotp^2)* (K2 + Q2 - 2*Kdotq(y))))
        # )::Float64

        # kernel23[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((2 *Kdotp *(K2^2 *(P2 *innerA1 *innerA2 *Kdotp + (4* innerB1 *innerB2 + innerA1* innerA2 *(P2 - 4* Q2 - 4 *Kdotq(y)))* Pdotq) + Kdotp* ((4* innerB1 *innerB2 + innerA1 *innerA2* (P2 + 8*Q2 - 16 *Kdotq(y))) *Kdotq(y)^2 + innerA1* innerA2 *Kdotp^2 *(-Q2 + 4* Kdotq(y)) - 2 *innerA1 *innerA2* Kdotp* Kdotq(y)* Pdotq) + K2 *(-innerA1* innerA2 *Kdotp^3 - 2 *innerA1* innerA2 *Kdotp^2 *Pdotq + Kdotq(y) *(-4* innerB1* innerB2 + innerA1 *innerA2 *(-P2 - 8* Q2 + 16 *Kdotq(y)))* Pdotq + Kdotp *((-4 *innerB1 *innerB2 + innerA1* innerA2* (-3* P2 + 4* Q2))* Kdotq(y) + 4* innerA1* innerA2* Kdotq(y)^2 + innerA1* innerA2 *(P2 *Q2 + 2 *Pdotq^2)))))/(3*(-P2* K2 + Kdotp^2)* (K2 + Q2 - 2 *Kdotq(y))))
        # )::Float64

        # kernel24[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((2 *(K2^2 *(P2* ((-innerA2* innerB1 + innerA1 *innerB2)* Kdotp + 2 *(innerA2* innerB1 + innerA1* innerB2) *Kdotq(y)) + 2 *P2* (-innerA2* innerB1 + innerA1* innerB2) *Pdotq + 4 *(innerA2 *innerB1 + innerA1* innerB2) *Pdotq^2) + Kdotp *(2* P2 *(-innerA2* innerB1 + innerA1 *innerB2)* Kdotq(y)^2 - (innerA2* innerB1 - innerA1* innerB2)* Kdotp^2 *(-Q2 + 4* Kdotq(y)) + 2* Kdotp* Kdotq(y)* (-((innerA2* innerB1 + innerA1* innerB2)* Q2) + 4 *(innerA2* innerB1 + innerA1 *innerB2) *Kdotq(y) + (innerA2 *innerB1 - innerA1* innerB2)* Pdotq)) + K2* ((innerA2* innerB1 - innerA1* innerB2) *Kdotp^3 + 2 *P2 *Kdotq(y) *((innerA2* innerB1 + innerA1 *innerB2)* Q2 - 2* (innerA2* innerB1 + innerA1 *innerB2)* Kdotq(y) + (innerA2 *innerB1 - innerA1 *innerB2)* Pdotq) - 2* Kdotp^2 *((innerA2 *innerB1 + innerA1 *innerB2)* Kdotq(y) + (-innerA2* innerB1 + innerA1* innerB2)* Pdotq) + Kdotp *(4 *Kdotq(y)* (P2 *(innerA2* innerB1 - innerA1 *innerB2) - 2* (innerA2* innerB1 + innerA1 *innerB2) *Pdotq) - (innerA2* innerB1 -innerA1 *innerB2) *(P2 *Q2 + 2* Pdotq^2)))))/(3* (-P2 *K2 + Kdotp^2) *(K2 + Q2 - 2* Kdotq(y))))
        # )::Float64

        # kernel31[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # (-((4* ((innerA2 *innerB1 + innerA1 *innerB2)* Kdotp^3 - 2* Kdotp^2 *((-innerA2* innerB1 + innerA1* innerB2) *Q2 + (innerA2 *innerB1 - innerA1* innerB2)* Kdotq(y) + (innerA2* innerB1 + innerA1 *innerB2)* Pdotq) + Kdotp *(-K2 *(P2 *(innerA2* innerB1 + innerA1 *innerB2) + (innerA2* innerB1 - innerA1* innerB2)* Pdotq) + Kdotq(y)* (P2 *(innerA2* innerB1 + innerA1 *innerB2) + 4* (innerA2* innerB1 - innerA1* innerB2) *Pdotq) + Pdotq* (3 *(-innerA2* innerB1 + innerA1* innerB2) *Q2 + (innerA2 *innerB1 + innerA1* innerB2)* Pdotq)) + P2* (K2 *(2 *(-innerA2* innerB1 + innerA1 *innerB2) *Q2 + 3 *(innerA2* innerB1 - innerA1 *innerB2)* Kdotq(y) + (innerA2* innerB1 + innerA1* innerB2) *Pdotq) - Kdotq(y) *(3 *(-innerA2* innerB1 + innerA1 *innerB2) *Q2 + 4 *(innerA2* innerB1 - innerA1* innerB2) *Kdotq(y) + (innerA2 *innerB1 + innerA1 *innerB2) *Pdotq))))/(3 *Kdotp *(-P2* K2 + Kdotp^2)* (K2 + Q2 - 2* Kdotq(y)))))
        # )::Float64

        # kernel32[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((2* ((4* innerB1* innerB2 - innerA1* innerA2 *(P2 + 4* Q2)) *Kdotp^3 + 2* Kdotp^2 *(-4 *innerB1 *innerB2 + innerA1 *innerA2 *(P2 + 4 *Kdotq(y))) *Pdotq + P2 *(K2* (4* innerB1* innerB2 - innerA1* innerA2* (P2 - 4 *Q2 + 12* Kdotq(y))) + Kdotq(y) *(-4 *innerB1* innerB2 + innerA1* innerA2* (P2 - 8* Q2 + 16* Kdotq(y)))) *Pdotq + Kdotp* ((4 *innerB1* innerB2 - innerA1* innerA2* (P2 - 8* Q2))* Pdotq^2 + Kdotq(y) *(-P2^2 *innerA1 *innerA2 + 4 *P2 *(innerB1* innerB2 - innerA1 *innerA2* Q2) - 16* innerA1* innerA2 *Pdotq^2) + K2* (P2^2* innerA1 *innerA2 + P2* (-4* innerB1 *innerB2 + 4* innerA1* innerA2* Q2) + 4* innerA1* innerA2* Pdotq^2))))/(3* Kdotp* (-P2 *K2 + Kdotp^2)* (K2 + Q2 - 2* Kdotq(y))))
        # )::Float64

        # kernel33[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((-3 *P2 *K2^2 *(4* innerB1* innerB2 + innerA1* innerA2* (P2 - 4 *Q2)) - 4 *innerA1* innerA2* Kdotp^4 + 2 *P2 *Kdotq(y)^2 *(-4 *innerB1 *innerB2 + innerA1* innerA2* (-P2 - 8* Q2 + 16 *Kdotq(y))) + 8* innerA1* innerA2* Kdotp^3 *Pdotq + 2* Kdotp* (4 *innerB1 *innerB2 + innerA1* innerA2* (3 *P2 + 8* Q2 - 16 *Kdotq(y))) *Kdotq(y) *Pdotq + K2* ((12* innerB1* innerB2 + innerA1* innerA2* (7* P2 - 12 *Q2)) *Kdotp^2 + P2 *(4 *innerA1 *innerA2 *Q2^2 + 6 *(4 *innerB1* innerB2 + innerA1* innerA2* (P2 - 4* Kdotq(y))) *Kdotq(y) - Q2* (4* innerB1 *innerB2 + innerA1 *innerA2* (P2 + 8* Kdotq(y)))) + 2* Kdotp* (-4* innerB1* innerB2 + innerA1* innerA2 *(-3 *P2 + 4 *Q2 + 4* Kdotq(y))) *Pdotq) + Kdotp^2* ((P2* innerA1* innerA2 + 4* innerB1* innerB2)* Q2 - 4* innerA1 *innerA2* Q2^2 - 8* (2 *innerB1 *innerB2 + innerA1* innerA2* (P2 - 2* Kdotq(y)))* Kdotq(y) - 4* innerA1* innerA2* Pdotq^2))/(3* (-P2 *K2 + Kdotp^2) *(K2 + Q2 - 2* Kdotq(y))))
        # )::Float64

        # kernel34[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((2 *(2 *(innerA2 *innerB1 - innerA1* innerB2)* Kdotp^4 - 4 *Kdotp^3 *((innerA2* innerB1 + innerA1* innerB2)* Kdotq(y) + (innerA2* innerB1 - innerA1* innerB2)* Pdotq) + 4* Kdotp* (P2 *(innerA2 *innerB1 + innerA1* innerB2) *(K2 - Kdotq(y))* Kdotq(y) + P2 *(innerA2* innerB1 - innerA1* innerB2)* (K2 - Kdotq(y))* Pdotq - (innerA2 *innerB1 + innerA1 *innerB2)* K2* Pdotq^2) + Kdotp^2 *(Q2 *(P2 *(-innerA2* innerB1 + innerA1* innerB2) + 2 *(innerA2 *innerB1 + innerA1* innerB2)* Pdotq) + K2* (P2 *(-5 *innerA2* innerB1 + 5 *innerA1 *innerB2) + 6 *(innerA2* innerB1 + innerA1* innerB2) *Pdotq) + 2* (innerA2* innerB1 - innerA1* innerB2) *(3 *P2 *Kdotq(y) + Pdotq^2)) + P2 *(2 *P2 *(innerA2* innerB1 - innerA1* innerB2)* Kdotq(y)^2 + 3* K2^2* (P2* (innerA2* innerB1 - innerA1 *innerB2) - 2 *(innerA2* innerB1 + innerA1* innerB2) *Pdotq) + K2 *(Q2* (P2 *(innerA2* innerB1 - innerA1* innerB2) - 2* (innerA2* innerB1 + innerA1* innerB2) *Pdotq) + Kdotq(y)* (P2 *(-6* innerA2* innerB1 + 6 *innerA1* innerB2) + 8* (innerA2 *innerB1 + innerA1* innerB2)* Pdotq)))))/(3 *Kdotp *(-P2* K2 + Kdotp^2)* (K2 + Q2 - 2 *Kdotq(y))))
        # )::Float64

        # kernel41[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((4 *innerA1 *innerA2* (-2 *Q2* Kdotp^2 - P2 *Q2* Kdotq(y) + Kdotp *(Q2 + 2 *Kdotq(y)) *Pdotq + K2* (-P2 *(-2* Q2 + Kdotq(y)) + (Kdotp - 2 *Pdotq)* Pdotq)))/(3 *(-P2* K2 + Kdotp^2) *(K2 + Q2 - 2 *Kdotq(y))))
        # )::Float64

        # kernel42[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((4 *(innerA2* innerB1 + innerA1 *innerB2)* (2 *Q2* Kdotp^2 + P2* Q2 *Kdotq(y) - Kdotp* (Q2 + 2 *Kdotq(y))* Pdotq + K2* (P2 *(-2* Q2 + Kdotq(y)) + Pdotq* (-Kdotp + 2 *Pdotq))))/(3* (-P2* K2 + Kdotp^2)* (K2 + Q2 - 2 *Kdotq(y))))
        # )::Float64

        # kernel43[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # (-((2* Kdotp* (2* P2 *(-innerA2 *innerB1 + innerA1* innerB2)* Kdotq(y)^2 - (innerA2 *innerB1 - innerA1* innerB2) *Kdotp^2* (Q2 + 2 *Kdotq(y)) + 2 *Kdotp *Kdotq(y) *(-((innerA2 *innerB1 + innerA1* innerB2)* Q2) + 2 *(innerA2* innerB1 + innerA1* innerB2) *Kdotq(y) + 2* (innerA2* innerB1 - innerA1* innerB2)* Pdotq) + K2^2* (P2 *(-innerA2* innerB1 + innerA1* innerB2) + 2* (innerA2 *innerB1 + innerA1* innerB2) *Pdotq) + K2* ((innerA2* innerB1 - innerA1* innerB2) *Kdotp^2 - 2 *(innerA2 *innerB1 + innerA1* innerB2)* Kdotp *Kdotq(y) + 2* (-innerA2 *innerB1 + innerA1* innerB2)* Pdotq^2 + 2 *Kdotq(y)* (P2 *(innerA2* innerB1 - innerA1* innerB2) - 2* (innerA2* innerB1 + innerA1 *innerB2) *Pdotq) + Q2* (P2* (innerA2* innerB1 - innerA1* innerB2) + 2 *(innerA2* innerB1 + innerA1 *innerB2)* Pdotq))))/(3 *(-P2 *K2 +Kdotp^2) *(K2 + Q2 - 2* Kdotq(y)))))
        # )::Float64

        # kernel44[i,j]=GaussChebyshevIntegral64(y->
        # -allkindsofweight*
        # D(Ksubq2(y))/branch[j]*
        # ((4* innerA1* innerA2* Q2^2* Kdotp^2 + Q2* Kdotp^2* (4 *innerB1* innerB2 + innerA1* innerA2* (P2 - 8* Kdotq(y))) + 2 *(P2* innerA1* innerA2 + 4* innerB1 *innerB2)* Kdotq(y)* (Kdotp^2 + P2* Kdotq(y) - 2 *Kdotp *Pdotq) + K2^2 *(P2^2 *innerA1* innerA2 + 4* P2* (innerB1* innerB2 + innerA1* innerA2* Q2) - 8 *innerA1* innerA2* Pdotq^2) - K2* (4 *P2 *innerA1* innerA2* Q2^2 + (4* innerB1 *innerB2 + innerA1* innerA2* (P2 + 4* Q2))* Kdotp^2 + P2* Q2* (4 *innerB1* innerB2 + innerA1* innerA2* (P2 - 8* Kdotq(y))) + 2* P2 *Kdotq(y) *(4 *innerB1* innerB2 + innerA1* innerA2* (P2 + 4* Kdotq(y))) - 16* innerA1* innerA2* Kdotp* Kdotq(y) *Pdotq - 2 *(P2 *innerA1* innerA2 + 4 *innerB1 *innerB2)* Pdotq^2))/(3 *(-P2* K2 + Kdotp^2)* (K2 + Q2 - 2 *Kdotq(y))))
        # )::Float64
    end
end


