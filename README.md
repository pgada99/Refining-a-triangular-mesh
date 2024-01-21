# Refining-a-triangular-mesh
The project we have decided to develop involves writing a series of functions capable of refining a two-dimensional domain, which is provided to us as material.
Our task is to further cut the domain, and for this purpose, we have written a series of functions with different functionalities, such as:

Marking the bases of triangles marked for refinement.
Refining the triangular mesh in order to obtain a conforming mesh. To achieve this result, once we cut the side of a triangle along its midpoint, we deactivate it, iterating this process until there are no more triangle sides to cut (this operation must be carried out for each triangle in the mesh). At the same time, the adjacency information for the sides of the triangulation is also updated.
Updating the adjacency information for the triangles.
Uniformly refining the triangulation (each triangle is divided into 4 similar triangles by joining the midpoints of the sides of the parent triangle).
