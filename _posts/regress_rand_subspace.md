---
layout: post
title: "Lower Bound for Linear Regression between Random Subspaces"
date:   2025-05-20
author: Satwik
# image:  /blog_assets/thumbs/welcome.png   # ← optional; delete this line if you have no thumbnail
usemathjax: true                     # ← only if the post contains LaTeX
---


**[Hanson-Wright inequality {\cite[Thm. 6.2.2]{Vershynin_2018}}]**
Let $A\in\R^{n\times n}$ be arbitrary matrix and let 
$X=(X_1,\dots,X_n)^{\top}\in\R^{n}$
have independent, mean-zero, sub-Gaussian coordinates.
Then, for every $t\ge0$,
$$
\Pr\Bigl\{\,
      \bigl|\,X^{\top}AX - \E\bigl[X^{\top}AX\bigr]\bigr|
       \ge t
   \Bigr\}
 \le 
2\exp\Bigl[
   -\,c 
   \min\Bigl\{
        \tfrac{t^{2}}{K^{4}\,\|A\|_{F}^{2}},
         \tfrac{t}{K^{2}\,\|A\|_{2}}
   \Bigr\}
\Bigr],
$$
where 
$
K=\max_{i}\|X_i\|_{\psi_2}
$
and $c>0$ is an absolute universal constant.






Let $\mA,\mB\in\R^{m\times k}$ have i.i.d.\ $\cN(0,1)$ entries, and assume $k\le m/2$.  Draw 
\[
\vx \sim \cN(\vzero,\mI_m),
\]
independently of $\mA,\mB$.  We will show that there exist absolute constants $c_1,c_2>0$ such that
\[
\Pr\Bigl\{\inf_{\mC\in\R^{k\times k}}
    \E_{\vx}\bigl\|\vx^{\top}\mA\mC-\vx^{\top}\mB\bigr\|_2^2
   < c_1\,m\,k\Bigr\}
\le 2\exp\bigl(-c_2\,m\,k\bigr).
\]


\noindent See that, for any fixed $\mC\in\R^{k\times k}$,
\[
\E_{\vx}\bigl\|\vx^{\top}(\mA\mC-\mB)\bigr\|_2^2
=\tr\bigl((\mA\mC-\mB)^{\top}(\mA\mC-\mB)\E_{\vx}[\vx \vx^{\top}]\bigr) =\tr\bigl((\mA\mC-\mB)^{\top}(\mA\mC-\mB)\bigr)
=\|\mA\mC-\mB\|_F^2.
\]
Hence we have
\[
\inf_{\mC}\E_{\vx}\|\vx^{\top}\mA\mC-\vx^{\top}\mB\|_2^2
 = 
\inf_{\mC}\|\mA\mC-\mB\|_F^2.
\]

\medskip

\noindent Our goal now is to derive a lower bound for $\inf_{\mC}\|\mA\mC-\mB\|_F^2$. We will use a standard fact that random Gaussian matrices have full column rank almost surely (with probability $1$). Thus, the matrix $A$ has full column rank and the concatenated matrix $[A,B] \in \br^{m \times 2k}$ has full rank as well implying that the columns of $A$ and $B$ are linearly independent of each other. 

The function $f(\mC) =\|\mA\mC-\mB\|_F^2$ has a unique minimizer $\mC^{*}= \mA^{\dagger}\mB$ given by the Moore-Penrose pseudoinverse $\mA^\dagger = (\mA^{\top}\mA)^{-1}\mA^{\top}$ which exists since $\mA$ has full column rank. The orthogonal
projector onto its column space is given by
\[
\mP_{\mA} = \mA\,\mA^\dagger,
\quad \text{and let }
\mR = \mI_m - \mP_{\mA}.
\]

\noindent By expansion with $\mC=\mA^\dagger\mB$, see that
\[
\inf_{\mC}\|\mA\mC-\mB\|_F^2 
 =  \|\mA(\mA^{\top}\mA)^{-1}\mA^{\top}\mB-\mB\|_F^2
 = \|(\mI-\mP_{\mA})\,\mB\|_F^2
 = \sum_{j=1}^k\bigl\|\mR\,\vb_j\bigr\|_2^2,
\]
where $\vb_j\in\R^m$ is the $j$th column of $\mB$.

\medskip


\noindent Conditional on $\mA$, the matrices $\mP_{\mA}, \mR$ are fixed symmetric matrices of rank $k$ and $m-k$ respectively. Additionally, 
\[
\mP_{\mA}^2 = \mA(\mA^{\top}\mA)^{-1}\mA^{\top}\mA(\mA^{\top}\mA)^{-1}\mA^{\top} = \mA(\mA^{\top}\mA)^{-1}\mA^{\top} = \mP_{\mA}.
\]
This implies that all the eigenvalues of $\mP_{\mA}$ and consequently $\mR$ are either $0$ or $1$. Since $\mR$ has rank $m - k$,
\[
\|\mR\|_2\le1,
\qquad
\tr(\mR)=\|\mR\|_F^2=m-k.
\]

Each column $\vb_j\sim N(\vzero,\mI_m)$ is independent of $\mR$, and
\[
\E\|\mR\,\vb_j\|_2^2=\tr(R)=m-k. 
\]
 By the Hanson–Wright inequality
(Thm~\ref{thm:hansonwright}) , for any $t>0$,
\[
\Pr\Bigl\{\bigl|\|\mR\,\vb_j\|_2^2 - (m-k)\bigr|\ge t
       \Bigm| \mA\Bigr\}
 \le 
2\exp\Bigl(-c\,
  \min\bigl\{\tfrac{t^2}{\|\mR\|_F^2},\,\tfrac{t}{\|\mR\|_2}\bigr\}\Bigr)
 \le 
2\exp\Bigl(-c\min\{t^2/(m-k),\,t\}\Bigr).
\]
Setting $t=(m-k)/2$ and using $m-k\ge m/2$, we get
\[
\Pr\bigl\{\|\mR\,\vb_j\|_2^2 < \tfrac12(m-k) \big| \mA\bigr\}
 \le 
2\exp\bigl(-c_3\,m\bigr).
\]
Because the $b_j$’s are independent, a union bound over $j=1,\dots,k$ yields
\[
\Pr\Bigl\{\exists\,j:\,\|\mR\,\vb_j\|_2^2<\tfrac12(m-k)\Bigr\}
 \le 
k\cdot2e^{-c_3m}
 \le 
2e^{-c_4\,m},
\]
where in the last inequality we used $k\le m/2$ and $n e^{-cn} = e^{-O(n)}$.


Thus, with probability $\ge 1-2e^{-cm}$, each
$\|\mR\,\vb_j\|_2^2\ge\tfrac12(m-k)$, so
\[
\inf_{\mC}\|\mA\mC-\mB\|_F^2
=\sum_{j=1}^k\|\mR\,\vb_j\|_2^2
 \ge 
k\cdot\tfrac12(m-k)
 \ge 
\tfrac14\,m\,k.
\]

