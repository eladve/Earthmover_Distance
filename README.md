# Python code for Earth Mover Distance between Multidimensional Pointsets 

The file EMD.py implements a nice little function computing earthmover distance (EMD) between two multidimensional point sets efficiently, using network flow algorithms. This is the right way to compute EMD between two sets of points -- a lot of libraries implement it very inefficiently.

Sample use: `resulting_flow, EMD_value = compute_optimal_transport(S, T)`
 
Then `EMD_value` will be the EMD between the two point-sets (i.e. between the two distributions that are uniform on those point-sets). And `resulting_flow` is the actual optimal flow. The sets `S` and `T` need to be lists of numpy vectors. Each list should consist of unique vectors. If you want to include a vector more than once, use the optional weight parameter.

**More info:**

This function computes the  earthmover distance (EMD) between two high-dimensional pointsets, which we will call Source (S) and Target (T). (But recall the EMD is symmetric so those names are just for convenience.)

The EMD between two point-sets is defined as the EMD between the distributions defined to be uniform on the support on those two pointsets.

The cardinalities of the two point-sets do not have to be equal. If they were equal there would be a slightly more efficient algorithm based on maximum matching. But since the point sets might be of different cardinalitie, we will use an approach based on min cost flow. In any case, we will absolutely not use linear programming or the simplex algorithm or any kind of optimization of that sort -- that is wasteful and terrible. Instead we use network flow algorithms like God and Karp intended.

The EMD relies on an underlying distance metric that says how close two points are to one another. We set this by default to be L2, but it can be changed to any other metric. (It needs to be symmetric.) In particular you might want to choose some normalized-L2.

The code below is pretty ugly. It should be using numpy much more to do the operations using numpy functions rather than loops. Also there is a potential improvement in running time, remarked inside the code.

This is barely tested -- use at your own risk.
