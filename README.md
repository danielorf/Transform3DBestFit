Transform3D class used to calculate best fit 3D transform between two sets of cartesian (X,Y,Z) points.
4x4 Transform matrix, 6DoF vector, Siemens NC style 6DoF vector, 4x4 Rotation matrix, 4x4 Translation matrix, fit error RMS and Rotation matrix determinant.
This project relies on the Math.Net Numerics library.


    Usage:          The Transform3D class takes in two sets of 3D points and calculates the best-fit transform and fit error RMS.  
                    Transform3D works with any number of points so long as the number and order of points in both sets match.
                    Incoming point arrays or matrices must be formatted such that each row represents a point with column0 as "X", 
                    column1 as "Y" and column2 as "Z".
                
    Terminology:    "actuals" refers to the starting points and "nominals" refers to ending points.  The calculated transform is 
                    the best fit for actuals->nominals transform
 
    References:     Details of SVD and rigid transform calculation can be found here:
                        http://graphics.stanford.edu/~smr/ICP/comparison/eggert_comparison_mva97.pdf
                        http://igl.ethz.ch/projects/ARAP/svd_rot.pdf
                        http://en.wikipedia.org/wiki/Wahba%27s_problem
                        http://en.wikipedia.org/wiki/Kabsch_algorithm
                        http://nghiaho.com/?page_id=671
 
    Math.NET License Info:  
                   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
                   documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
                   the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and 
                   to permit persons to whom the Software is furnished to do so, subject to the following conditions:

                   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
                   THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
                   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
                   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.