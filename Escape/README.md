estimate distance is the main function
    pseudocode:
        if p1 and p2 are in L2+, then run L2_to_L2 <br />
        run Ln_to_L0 for p1
        run Ln_to_L0 for p2
        if connected through L1, return
        run L0_to_L0 (or L0_to_L0_optimized if PLL)
        return

L2_to_L2 explores the L2+ neighborhoods of P1 and P2
    pseudocode
        BiBFS both meighborhoods
        if in same neighborhood, return distance
        else, return L2 points they are connected to and the distance to L2

Ln_to_L0
    pseudocode:
        if p in L0, return
        if p in L1
            if p1 and p2 are connected through L1, return distance and early connection flag
            run L1_to_L0
            return
        if p in L2+
            if L2+ neighborhood not yet traversed, run Ln_to_L2
            collect neighbors of the L2 point in L1
            if p1 and p2 are connected through L1, return distance and early connection flag
            run L1_to_L0
            return

Ln_to_L2
    BiBFS until we reach L2

L0_to_L0
    BiBFS between source and destination nodes until a path is found