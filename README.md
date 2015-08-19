# Impala-avx2-parquet-scanner

    optimize parquet scanner with avx2
    
        1. implement fixed length encoding with avx2
    
        2. evaluate predicate (comparison, in, between) on encoded data(encoding format is 1).
    
                         base      avx2-optimized
        TPCH Q6          3.5175       2.4525
        TPCH Q12         11.9525      8.14
        TPCH Q14         10.855       5.4075
        TPCH Q13         84.5875      58.445
    
        selectivity      base      avx2-optimized
        1%               14.49        6.365
        10%              14.67        7.295
        20%              15.93        8.37
        30%              15.7         9.27
        50%              16.0625      10.7275
        80%              16.5525      12.755
        100%             16.0675      13.375



