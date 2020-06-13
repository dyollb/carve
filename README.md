# carve
Carve is a fast, robust constructive solid geometry library

Carve provides code to perform boolean operations on polyhedral & triangle meshes.
It is relatively fast and robust. A nice feature of carve is its generic framework to 
collect different types of results from the intersection test using custom collectors 
(e.g. union, intersection, imprint, etc.).
So-called hooks can be registered. Hooks are callbacks that receive the original face, and the new (cut) face, allowing user
defined operations like interpolation from input to output.

Carve was originally developed by [Tobias Sargeant](https://github.com/folded/carve). The original repository is not maintained.
This is an attempt to modernize the library, and eventually improve performance and robustness.
