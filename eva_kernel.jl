# using ProfileVega
# include("./integral.jl")
##计算kernel
kernel11=[GaussChebyshevIntegral(
                y->(
                    -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
                    *(A1(j)*A2(j)*(P2-4*k(j))-4*B1(j)*B2(j))
                    )
,ystep) for i=1:dim,j=1:dim];

# result = Array{Float64}(undef, dim, dim)
# @time Threads.@threads for i = 1:dim
#     for j = 1:dim
#         ft(y) = -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))
#             *(A1(j)*A2(j)*(P2-4*k(j))-4*B1(j)*B2(j))
#         result[i, j] = GaussChebyshevIntegral(ft, ystep)
#     end
# end

kernel12=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *(-2*P2*(A2(j)*B1(j) + A1(j)*B2(j))+4*(A2(j)*B1(j) - A1(j)*B2(j))*pdotq(j))
        )
,ystep) for i=1:dim,j=1:dim];


kernel13=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *(-2*kdotp(i)*((A2(j)*B1(j) + A1(j)*B2(j))*kdotp(i) + 2*(-A2(j)*B1(j) + A1(j)*B2(j))*kdotq(i, j, y)))
        )
,ystep) for i=1:dim,j=1:dim]; 

kernel14=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *(4*A1(j)*A2(j)*(P2*kdotq(i, j, y) - kdotp(i)*pdotq(j)))
        )
,ystep) for i=1:dim,j=1:dim];

kernel21=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((2*(k(i)^2*(P2*(A2(j)*B1(j) + A1(j)*B2(j)) + 2*(-A2(j) *B1(j) + A1(j)* B2(j))* pdotq(j)) + k(i) *(-((A2(j)* B1(j) + A1(j)* B2(j)) *kdotp(i)^2) + 2*(A2(j)* B1(j) + A1(j)* B2(j))* pdotq(j)^2 + 2*kdotq(i, j,y) *(-P2 *(A2(j)* B1(j) + A1(j) *B2(j)) + 4 *(A2(j)* B1(j) - A1(j)* B2(j)) *pdotq(j)) + k(j)* (P2 *(A2(j)* B1(j) + A1(j)* B2(j)) + 6 *(-A2(j)* B1(j) + A1(j) *B2(j))* pdotq(j)) + 2*kdotp(i)*((A2(j)* B1(j) - A1(j)* B2(j)) *kdotq(i, j, y) - (A2(j) *B1(j) + A1(j) *B2(j))* pdotq(j))) + kdotp(i) *((A2(j)* B1(j) + A1(j)* B2(j))* kdotp(i) *(-k(j) + 4*kdotq(i, j, y)) - 2*kdotq(i, j, y)*(3 *(-A2(j) *B1(j) + A1(j)* B2(j)) *k(j) + 4* (A2(j)* B1(j) - A1(j) *B2(j)) *kdotq(i, j, y) + (A2(j) *B1(j) + A1(j)* B2(j)) *pdotq(j)))))/(3*(-P2*k(i) + kdotp(i)^2)*(k(i) + k(j) - 2*kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel22=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((kdotp(i) *((-4 *B1(j)* B2(j) + A1(j)* A2(j)* (P2 + 4* k(j)))*kdotp(i) *(-k(j) + 4*kdotq(i, j, y)) - 2 *kdotq(i, j, y)* (-4* B1(j)* B2(j) + A1(j)* A2(j)* (P2 - 8 *k(j) + 16 *kdotq(i, j, y)))* pdotq(j)) + k(i)^2 *(P2^2 *A1(j)* A2(j) + P2 *(-4 *B1(j) *B2(j) + 4* A1(j)* A2(j)* k(j)) - 8* A1(j)* A2(j) *pdotq(j)^2) + k(i) *(4 *P2 *A1(j) *A2(j) *k(j)^2 + (4 *B1(j) *B2(j) - A1(j)* A2(j) *(P2 + 4*k(j))) *kdotp(i)^2 - 2* P2* (P2* A1(j)* A2(j) - 4 *B1(j)* B2(j))* kdotq(i, j, y) + 2 *kdotp(i) *(4* B1(j)* B2(j) + A1(j)* A2(j)*(-P2 - 4* k(j) + 4* kdotq(i, j, y)))* pdotq(j) + 2* (-4* B1(j)* B2(j) + A1(j)* A2(j)* (P2 + 16 *kdotq(i, j, y)))* pdotq(j)^2 + k(j) *(P2 *(-4* B1(j)* B2(j) + A1(j) *A2(j) *(P2 - 8* kdotq(i, j, y))) - 16 *A1(j)* A2(j)* pdotq(j)^2)))/(3 *(-P2* k(i) + kdotp(i)^2)* (k(i) + k(j) - 2*kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel23=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((2 *kdotp(i) *(k(i)^2 *(P2 *A1(j) *A2(j) *kdotp(i) + (4* B1(j) *B2(j) + A1(j)* A2(j) *(P2 - 4* k(j) - 4 *kdotq(i, j, y)))* pdotq(j)) + kdotp(i)* ((4* B1(j) *B2(j) + A1(j) *A2(j)* (P2 + 8*k(j) - 16 *kdotq(i, j, y))) *kdotq(i, j, y)^2 + A1(j)* A2(j) *kdotp(i)^2 *(-k(j) + 4* kdotq(i, j, y)) - 2 *A1(j) *A2(j)* kdotp(i)* kdotq(i, j, y)* pdotq(j)) + k(i) *(-A1(j)* A2(j) *kdotp(i)^3 - 2 *A1(j)* A2(j) *kdotp(i)^2 *pdotq(j) + kdotq(i, j, y) *(-4* B1(j)* B2(j) + A1(j) *A2(j) *(-P2 - 8* k(j) + 16 *kdotq(i, j, y)))* pdotq(j) + kdotp(i) *((-4 *B1(j) *B2(j) + A1(j)* A2(j)* (-3* P2 + 4* k(j)))* kdotq(i, j, y) + 4* A1(j)* A2(j)* kdotq(i, j, y)^2 + A1(j)* A2(j) *(P2 *k(j) + 2 *pdotq(j)^2)))))/(3*(-P2* k(i) + kdotp(i)^2)* (k(i) + k(j) - 2 *kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim]

kernel24=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((2 *(k(i)^2 *(P2* ((-A2(j)* B1(j) + A1(j) *B2(j))* kdotp(i) + 2 *(A2(j)* B1(j) + A1(j)* B2(j)) *kdotq(i, j, y)) + 2 *P2* (-A2(j)* B1(j) + A1(j)* B2(j)) *pdotq(j) + 4 *(A2(j) *B1(j) + A1(j)* B2(j)) *pdotq(j)^2) + kdotp(i) *(2* P2 *(-A2(j)* B1(j) + A1(j) *B2(j))* kdotq(i, j, y)^2 - (A2(j)* B1(j) - A1(j)* B2(j))* kdotp(i)^2 *(-k(j) + 4* kdotq(i, j, y)) + 2* kdotp(i)* kdotq(i, j, y)* (-((A2(j)* B1(j) + A1(j)* B2(j))* k(j)) + 4 *(A2(j)* B1(j) + A1(j) *B2(j)) *kdotq(i, j, y) + (A2(j) *B1(j) - A1(j)* B2(j))* pdotq(j))) + k(i)* ((A2(j)* B1(j) - A1(j)* B2(j)) *kdotp(i)^3 + 2 *P2 *kdotq(i, j, y) *((A2(j)* B1(j) + A1(j) *B2(j))* k(j) - 2* (A2(j)* B1(j) + A1(j) *B2(j))* kdotq(i, j, y) + (A2(j) *B1(j) - A1(j) *B2(j))* pdotq(j)) - 2* kdotp(i)^2 *((A2(j) *B1(j) + A1(j) *B2(j))* kdotq(i, j, y) + (-A2(j)* B1(j) + A1(j)* B2(j))* pdotq(j)) + kdotp(i) *(4 *kdotq(i, j, y)* (P2 *(A2(j)* B1(j) - A1(j) *B2(j)) - 2* (A2(j)* B1(j) + A1(j) *B2(j)) *pdotq(j)) - (A2(j)* B1(j) -A1(j) *B2(j)) *(P2 *k(j) + 2* pdotq(j)^2)))))/(3* (-P2 *k(i) + kdotp(i)^2) *(k(i) + k(j) - 2* kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel31=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *(-((4* ((A2(j) *B1(j) + A1(j) *B2(j))* kdotp(i)^3 - 2* kdotp(i)^2 *((-A2(j)* B1(j) + A1(j)* B2(j)) *k(j) + (A2(j) *B1(j) - A1(j)* B2(j))* kdotq(i, j, y) + (A2(j)* B1(j) + A1(j) *B2(j))* pdotq(j)) + kdotp(i) *(-k(i) *(P2 *(A2(j)* B1(j) + A1(j) *B2(j)) + (A2(j)* B1(j) - A1(j)* B2(j))* pdotq(j)) + kdotq(i, j, y)* (P2 *(A2(j)* B1(j) + A1(j) *B2(j)) + 4* (A2(j)* B1(j) - A1(j)* B2(j)) *pdotq(j)) + pdotq(j)* (3 *(-A2(j)* B1(j) + A1(j)* B2(j)) *k(j) + (A2(j) *B1(j) + A1(j)* B2(j))* pdotq(j))) + P2* (k(i) *(2 *(-A2(j)* B1(j) + A1(j) *B2(j)) *k(j) + 3 *(A2(j)* B1(j) - A1(j) *B2(j))* kdotq(i, j, y) + (A2(j)* B1(j) + A1(j)* B2(j)) *pdotq(j)) - kdotq(i, j, y) *(3 *(-A2(j)* B1(j) + A1(j) *B2(j)) *k(j) + 4 *(A2(j)* B1(j) - A1(j)* B2(j)) *kdotq(i, j, y) + (A2(j) *B1(j) + A1(j) *B2(j)) *pdotq(j)))))/(3 *kdotp(i) *(-P2* k(i) + kdotp(i)^2)* (k(i) + k(j) - 2* kdotq(i, j, y)))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel32=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((2* ((4* B1(j)* B2(j) - A1(j)* A2(j) *(P2 + 4* k(j))) *kdotp(i)^3 + 2* kdotp(i)^2 *(-4 *B1(j) *B2(j) + A1(j) *A2(j) *(P2 + 4 *kdotq(i, j, y))) *pdotq(j) + P2 *(k(i)* (4* B1(j)* B2(j) - A1(j)* A2(j)* (P2 - 4 *k(j) + 12* kdotq(i, j, y))) + kdotq(i, j, y) *(-4 *B1(j)* B2(j) + A1(j)* A2(j)* (P2 - 8* k(j) + 16* kdotq(i, j, y)))) *pdotq(j) + kdotp(i)* ((4 *B1(j)* B2(j) - A1(j)* A2(j)* (P2 - 8* k(j)))* pdotq(j)^2 + kdotq(i, j, y) *(-P2^2 *A1(j) *A2(j) + 4 *P2 *(B1(j)* B2(j) - A1(j) *A2(j)* k(j)) - 16* A1(j)* A2(j) *pdotq(j)^2) + k(i)* (P2^2* A1(j) *A2(j) + P2* (-4* B1(j) *B2(j) + 4* A1(j)* A2(j)* k(j)) + 4* A1(j)* A2(j)* pdotq(j)^2))))/(3* kdotp(i)* (-P2 *k(i) + kdotp(i)^2)* (k(i) + k(j) - 2* kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel33=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((-3 *P2 *k(i)^2 *(4* B1(j)* B2(j) + A1(j)* A2(j)* (P2 - 4 *k(j))) - 4 *A1(j)* A2(j)* kdotp(i)^4 + 2 *P2 *kdotq(i, j, y)^2 *(-4 *B1(j) *B2(j) + A1(j)* A2(j)* (-P2 - 8* k(j) + 16 *kdotq(i, j, y))) + 8* A1(j)* A2(j)* kdotp(i)^3 *pdotq(j) + 2* kdotp(i)* (4 *B1(j) *B2(j) + A1(j)* A2(j)* (3 *P2 + 8* k(j) - 16 *kdotq(i, j, y))) *kdotq(i, j, y) *pdotq(j) + k(i)* ((12* B1(j)* B2(j) + A1(j)* A2(j)* (7* P2 - 12 *k(j))) *kdotp(i)^2 + P2 *(4 *A1(j) *A2(j) *k(j)^2 + 6 *(4 *B1(j)* B2(j) + A1(j)* A2(j)* (P2 - 4* kdotq(i, j, y))) *kdotq(i, j, y) - k(j)* (4* B1(j) *B2(j) + A1(j) *A2(j)* (P2 + 8* kdotq(i, j, y)))) + 2* kdotp(i)* (-4* B1(j)* B2(j) + A1(j)* A2(j) *(-3 *P2 + 4 *k(j) + 4* kdotq(i, j, y))) *pdotq(j)) + kdotp(i)^2* ((P2* A1(j)* A2(j) + 4* B1(j)* B2(j))* k(j) - 4* A1(j) *A2(j)* k(j)^2 - 8* (2 *B1(j) *B2(j) + A1(j)* A2(j)* (P2 - 2* kdotq(i, j, y)))* kdotq(i,j, y) - 4* A1(j)* A2(j)* pdotq(j)^2))/(3* (-P2 *k(i) + kdotp(i)^2) *(k(i) + k(j) - 2* kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel34=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((2 *(2 *(A2(j) *B1(j) - A1(j)* B2(j))* kdotp(i)^4 - 4 *kdotp(i)^3 *((A2(j)* B1(j) + A1(j)* B2(j))* kdotq(i, j, y) + (A2(j)* B1(j) - A1(j)* B2(j))* pdotq(j)) + 4* kdotp(i)* (P2 *(A2(j) *B1(j) + A1(j)* B2(j)) *(k(i) - kdotq(i, j, y))* kdotq(i, j, y) + P2 *(A2(j)* B1(j) - A1(j)* B2(j))* (k(i) - kdotq(i, j, y))* pdotq(j) - (A2(j) *B1(j) + A1(j) *B2(j))* k(i)* pdotq(j)^2) + kdotp(i)^2 *(k(j) *(P2 *(-A2(j)* B1(j) + A1(j)* B2(j)) + 2 *(A2(j) *B1(j) + A1(j)* B2(j))* pdotq(j)) + k(i)* (P2 *(-5 *A2(j)* B1(j) + 5 *A1(j) *B2(j)) + 6 *(A2(j)* B1(j) + A1(j)* B2(j)) *pdotq(j)) + 2* (A2(j)* B1(j) - A1(j)* B2(j)) *(3 *P2 *kdotq(i, j, y) + pdotq(j)^2)) + P2 *(2 *P2 *(A2(j)* B1(j) - A1(j)* B2(j))* kdotq(i, j, y)^2 + 3* k(i)^2* (P2* (A2(j)* B1(j) - A1(j) *B2(j)) - 2 *(A2(j)* B1(j) + A1(j)* B2(j)) *pdotq(j)) + k(i) *(k(j)* (P2 *(A2(j)* B1(j) - A1(j)* B2(j)) - 2* (A2(j)* B1(j) + A1(j)* B2(j)) *pdotq(j)) + kdotq(i, j, y)* (P2 *(-6* A2(j)* B1(j) + 6 *A1(j)* B2(j)) + 8* (A2(j) *B1(j) + A1(j)* B2(j))* pdotq(j))))))/(3 *kdotp(i) *(-P2* k(i) + kdotp(i)^2)* (k(i) + k(j) - 2 *kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel41=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((4 *A1(j) *A2(j)* (-2 *k(j)* kdotp(i)^2 - P2 *k(j)* kdotq(i, j, y) + kdotp(i) *(k(j) + 2 *kdotq(i, j, y)) *pdotq(j) + k(i)* (-P2 *(-2* k(j) + kdotq(i, j, y)) + (kdotp(i) - 2 *pdotq(j))* pdotq(j))))/(3 *(-P2* k(i) + kdotp(i)^2) *(k(i) + k(j) - 2 *kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel42=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((4 *(A2(j)* B1(j) + A1(j) *B2(j))* (2 *k(j)* kdotp(i)^2 + P2* k(j) *kdotq(i, j, y) - kdotp(i)* (k(j) + 2 *kdotq(i, j, y))* pdotq(j) + k(i)* (P2 *(-2* k(j) + kdotq(i, j, y)) + pdotq(j)* (-kdotp(i) + 2 *pdotq(j)))))/(3* (-P2* k(i) + kdotp(i)^2)* (k(i) + k(j) - 2 *kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel43=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *(-((2* kdotp(i)* (2* P2 *(-A2(j) *B1(j) + A1(j)* B2(j))* kdotq(i, j, y)^2 - (A2(j) *B1(j) - A1(j)* B2(j)) *kdotp(i)^2* (k(j) + 2 *kdotq(i, j, y)) + 2 *kdotp(i) *kdotq(i, j, y) *(-((A2(j) *B1(j) + A1(j)* B2(j))* k(j)) + 2 *(A2(j)* B1(j) + A1(j)* B2(j)) *kdotq(i, j, y) + 2* (A2(j)* B1(j) - A1(j)* B2(j))* pdotq(j)) + k(i)^2* (P2 *(-A2(j)* B1(j) + A1(j)* B2(j)) + 2* (A2(j) *B1(j) + A1(j)* B2(j)) *pdotq(j)) + k(i)* ((A2(j)* B1(j) - A1(j)* B2(j)) *kdotp(i)^2 - 2 *(A2(j) *B1(j) + A1(j)* B2(j))* kdotp(i) *kdotq(i, j, y) + 2* (-A2(j) *B1(j) + A1(j)* B2(j))* pdotq(j)^2 + 2 *kdotq(i, j, y)* (P2 *(A2(j)* B1(j) - A1(j)* B2(j)) - 2* (A2(j)* B1(j) + A1(j) *B2(j)) *pdotq(j)) + k(j)* (P2* (A2(j)* B1(j) - A1(j)* B2(j)) + 2 *(A2(j)* B1(j) + A1(j) *B2(j))* pdotq(j)))))/(3 *(-P2 *k(i) +kdotp(i)^2) *(k(i) + k(j) - 2* kdotq(i, j, y)))))
        )
,ystep) for i=1:dim,j=1:dim];

kernel44=[GaussChebyshevIntegral(
    y->(
        -D(ksubq2(i,j,y))/(branch(qPlus2(j))*branch(qSubt2(j)))*sqrt(1-y^2)*sqrt(1-z(j)^2)
        *((4* A1(j)* A2(j)* k(j)^2* kdotp(i)^2 + k(j)* kdotp(i)^2* (4 *B1(j)* B2(j) + A1(j)* A2(j)* (P2 - 8* kdotq(i, j, y))) + 2 *(P2* A1(j)* A2(j) + 4* B1(j) *B2(j))* kdotq(i, j, y)* (kdotp(i)^2 + P2* kdotq(i, j, y) - 2 *kdotp(i) *pdotq(j)) + k(i)^2 *(P2^2 *A1(j)* A2(j) + 4* P2* (B1(j)* B2(j) + A1(j)* A2(j)* k(j)) - 8 *A1(j)* A2(j)* pdotq(j)^2) - k(i)* (4 *P2 *A1(j)* A2(j)* k(j)^2 + (4* B1(j) *B2(j) + A1(j)* A2(j)* (P2 + 4* k(j)))* kdotp(i)^2 + P2* k(j)* (4 *B1(j)* B2(j) + A1(j)* A2(j)* (P2 - 8* kdotq(i, j, y))) + 2* P2 *kdotq(i, j, y) *(4 *B1(j)* B2(j) + A1(j)* A2(j)* (P2 + 4* kdotq(i, j, y))) - 16* A1(j)* A2(j)* kdotp(i)* kdotq(i, j, y) *pdotq(j) - 2 *(P2 *A1(j)* A2(j) + 4 *B1(j) *B2(j))* pdotq(j)^2))/(3 *(-P2* k(i) + kdotp(i)^2)* (k(i) + k(j) - 2 *kdotq(i, j, y))))
        )
,ystep) for i=1:dim,j=1:dim];
