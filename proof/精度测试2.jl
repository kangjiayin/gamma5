kkk=[0.00010255209872165668, 0.00011417051335584808, 0.00013834102133968853, 0.00018212757130146563, 0.0002598456057921237, 0.000400432708855992, 0.0006638421378970346, 0.0011783424541740922, 0.002227526771379863, 0.004457968249600969, 0.009384267183563286, 0.020634544929970737, 0.047046312702388855, 0.11036987507699117, 0.264309122997184, 0.6408691073856274, 1.560381033311806, 3.783448670482154, 9.0604433438239, 21.255650922653995, 48.462420828459635, 106.56133083588223, 224.31743431315593, 448.92838678681306, 848.6497252625143, 1506.3822299197059, 2497.298492066069, 3848.4391412029468, 5490.656866800006, 7228.513931125004, 8758.828970867351, 9751.14124884138]
div1=Array{Float64}(undef, 32, 32)
div=[0 for i=1:32]
Threads.@threads for k in kkk
    global div1
    include("/Users/kjy/Desktop/program/julia/Gamma5/proof/精度测试.jl")
    div[k,:]=div
end   #for 