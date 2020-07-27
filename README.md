# Effect of Mixed Precision Computing on H-Matrix Vector Multiplication in BEM Analysis

Full paper at https://arxiv.org/abs/1911.00093, authored by Rise Ooi, Takeshi Iwashita, Takeshi Fukaya, Akihiro Ida, Rio Yokota. This is an accepted manuscript to International Conference on High Performance Computing in Asia-Pacific Region (HPCAsia2020), January 15--17, 2020, Fukuoka, Japan.

## Abstract

Hierarchical Matrix (H-matrix) is an approximation technique which splits a target dense matrix into multiple submatrices, and where a selected portion of submatrices are low-rank approximated. The technique substantially reduces both time and space complexity of dense matrix vector multiplication, and hence has been applied to numerous practical problems.
In this paper, we aim to accelerate the H-matrix vector multiplication by introducing mixed precision computing, where we employ both binary64 (FP64) and binary32 (FP32) arithmetic operations. We propose three methods to introduce mixed precision computing to H-matrix vector multiplication, and then evaluate them in a boundary element method (BEM) analysis. The numerical tests examine the effects of mixed precision computing, particularly on the required simulation time and rate of convergence of the iterative (BiCG-STAB) linear solver. We confirm the effectiveness of the proposed methods.

## First author's private comment

Personally, I think there is a huge potential of going into this direction as this essentially shows the effectiveness of applying mixed precision to solvers that can solve a wide variety of problems. I think the impact of pursuing the ideas I wrote in the Future Works section would be fruitful. Although I am no longer a full time academic and cannot update this repository as often as I used to, I try to do research whenever possible. So, if you have any questions, please feel free to contact me at the email written in the PDF.
