## VTK-m Volume Renderer

This repository contains a volume rendering algorithm I wrote and implemented in the VTK-m framework (http://m.vtk.org/index.php/Main_Page).

The VTK-m framework promotes effective use of high levels of concurrency through functional programming abstractions. Algorithms can be written once and easily target heterogenous platforms.

Here I used VTK-m to target CUDA and Threading Building Blocks (TBB) implementations of my volume renderer, and benchmarked them against a serial implementation.
