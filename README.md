# Efficient_path_planning

**Overview:**
This repository implements a path-planning and task assignment algorithm for very large multi-robot systems that should fulfill a global Boolean specification.

There are 2 main programs as it follows:

1) main_TAPF_chantry.m is for solving the classical **TAPF** problem (number of regions to be reached is equal with the number of robots; for this simulations ht_chantry map from [1] is used; start-goal configurations are randomly selected from the benchmark map);

Note: if after solving the first LP problem, the cell capacity s*, is greater than 1, in order to ensure collision avoidance the second LP problem will be solved using s* intermediary markings and cell capacity will be fixed to 1.

2) main_TAPF_aisle_boolean.m is for **global Boolean-based goal** (extended version of TAPF problem) and a fixed team of 100 robots. The Boolean specification has multiple disjunctions per term. The simulations are performed using a warehouse environment from [1], where the narrow corridors (of unitary width) are contributing to the problem complexity by requiring multiple intermediary markings.

Note: if after solving the first LP problem, the cell capacity s* is greater than 1, in order to ensure collision avoidance the second problem (MILP) will be solved using s* intermediary markings and cell capacity will be fixed to 1.

Note for users:
1) Users can choose if plots with trajectories will be displayed (by modifying the plot_animation flag).
2) Solving the ILP formulations of the problems might take a while for large teams of mobile agents. 

All optimization problems are solved using intlinprog solver (for the relaxed LPs problems, intcon option is [ ], meaning that no variable is imposed to be integer);

**Note:** All simulations were performed on a computer equipped with an AMD Ryzen 9 9950X 16 CPU and 64 GB RAM. Different perfomances may be obtained with different specifications.

[1] Stern, R., Sturtevant, N., Felner, A., Koenig, S., Ma, H., Walker, T., Li, J., Atzmon, D., Cohen, L., Kumar, T.K. and Barták, R., 2019. Multi-agent pathfinding: Definitions, variants, and benchmarks. In Proceedings of the International Symposium on Combinatorial Search (Vol. 10, No. 1, pp. 151-158).
