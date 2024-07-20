#!/bin/bash

# exe/createCGraph livejournal
exe/createCGraph pokec

# exe/generateInput livejournal
exe/generateInput pokec

# exe/sample BiBFS livejournal 0 10000
exe/sample BiBFS pokec 0 10000

# exe/createCGraph gbbs-soc-livejournal-sanitized_spanner_FGVbase
# cp inputs/livejournal_input.txt inputs/gbbs-soc-livejournal-sanitized_spanner_FGVbase_input.txt
# exe/sample BiBFS gbbs-soc-livejournal-sanitized_spanner_FGVbase 0 10000

# exe/createCGraph gbbs-soc-livejournal-sanitized_spanner_MPVXbase
# cp inputs/livejournal_input.txt inputs/gbbs-soc-livejournal-sanitized_spanner_MPVXbase_input.txt
# exe/sample BiBFS gbbs-soc-livejournal_sanitized_spanner_MPVXbase 0 10000

# exe/createCGraph gbbs-soc-livejournal-sanitized_spanner_MPVXcompact
# cp inputs/livejournal_input.txt inputs/gbbs-soc-livejournal-sanitized_spanner_MPVXcompact_input.txt
# exe/sample BiBFS gbbs-soc-livejournal-sanitized_spanner_MPVXcompact 0 10000

exe/createCGraph gbbs-pokec-sanitized_spanner_FGVbase
cp inputs/pokec_input.txt inputs/gbbs-pokec-sanitized_spanner_FGVbase_input.txt
exe/sample BiBFS gbbs-pokec-sanitized_spanner_FGVbase 0 10000

exe/createCGraph gbbs-pokec-sanitized_spanner_MPVXbase
cp inputs/pokec_input.txt inputs/gbbs-pokec-sanitized_spanner_MPVXbase_input.txt
exe/sample BiBFS gbbs-pokec-sanitized_spanner_MPVXbase 0 10000

exe/createCGraph gbbs-pokec-sanitized_spanner_MPVXcompact
cp inputs/pokec_input.txt inputs/gbbs-pokec-sanitized_spanner_MPVXcompact_input.txt
exe/sample BiBFS gbbs-pokec-sanitized_spanner_MPVXcompact 0 10000

# exe/createCGraph dblp-5-random
# cp inputs/dblp_input.txt inputs/dblp-5-random_input.txt
# exe/sample BiBFS dblp-5-random 0 10000

# exe/createCGraph dblp-8-random-wtd
# cp inputs/dblp_input.txt inputs/dblp-8-random-wtd_input.txt
# exe/sample BiBFS dblp-8-random-wtd 0 10000

# exe/createCGraph dblp-8-random
# cp inputs/dblp_input.txt inputs/dblp-8-random_input.txt
# exe/sample BiBFS dblp-8-random 0 10000

# exe/generateL0seed livejournal seed
exe/generateL0seed pokec seed

# exe/generateSparsifiedCOO livejournal_seed_10-0
# exe/createCGraph livejournal_seed_10-0_sparsified
# cp inputs/livejournal_input.txt inputs/livejournal_seed_10-0_sparsified_input.txt
# exe/sample BiBFS livejournal_seed_10-0_sparsified 0 10000


# exe/generateSparsifiedCOO livejournal_seed_8-0
# exe/createCGraph livejournal_seed_8-0_sparsified
# cp inputs/livejournal_input.txt inputs/livejournal_seed_8-0_sparsified_input.txt
# exe/sample BiBFS livejournal_seed_8-0_sparsified 0 10000

# exe/generateSparsifiedCOO livejournal_seed_6-0
# exe/createCGraph livejournal_seed_6-0_sparsified
# cp inputs/livejournal_input.txt inputs/livejournal_seed_6-0_sparsified_input.txt
# exe/sample BiBFS livejournal_seed_6-0_sparsified 0 10000

exe/generateSparsifiedCOO pokec_seed_10-0
exe/createCGraph pokec_seed_10-0_sparsified
cp inputs/pokec_input.txt inputs/pokec_seed_10-0_sparsified_input.txt
exe/sample BiBFS pokec_seed_10-0_sparsified 0 10000


exe/generateSparsifiedCOO pokec_seed_8-0
exe/createCGraph pokec_seed_8-0_sparsified
cp inputs/pokec_input.txt inputs/pokec_seed_8-0_sparsified_input.txt
exe/sample BiBFS pokec_seed_8-0_sparsified 0 10000

exe/generateSparsifiedCOO pokec_seed_6-0
exe/createCGraph pokec_seed_6-0_sparsified
cp inputs/pokec_input.txt inputs/pokec_seed_6-0_sparsified_input.txt
exe/sample BiBFS pokec_seed_6-0_sparsified 0 10000
