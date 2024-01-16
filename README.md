
# Wireless Reliability Fairness Optimization
We present an algorithm for wireless reliability fairness optimization that optimizes the min-max outage probability first published in [INFOCOM 2011](http://www.cs.cityu.edu.hk/~cheewtan/Tan-INFOCOM2011.pdf) (See longer version in [IEEE/ACM Transactions on Networking](file:///Users/siyac/Downloads/TanTon2015.pdf)). An overview survey is available in [Wireless Network Optimization by Perron-Frobenius Theory]{http://www.nowpublishers.com/article/Details/NET-048}.


## The Problem Statement
An outage event occurs at the Ith receiver when the received SINR falls below a given reliability threshold, i.e., $SINR_I(\pmb{p})<\beta_I$ for $I=1, \ldots, L$. So we are interested in minimizing the worst-case outage probability to ensure reliability fairness, which is formulated as follows :


$$
\begin{cases}
a_1x+b_1y+c_1z=d_1\\
a_2x+b_2y+c_2z=d_2\\
a_3x+b_3y+c_3z=d_3\\
\end{cases}
$$



$\text{minimize} \max_l P(SINR_l(\pmb{p})<\beta_l)$

$subject  to: \pmb{p} \in \mathcal{P}$

$variables:  p.$

where $SINR_l(\pmb{p}=R_{ll} G_{ll} p_l /(\sum_{j \neq l} R_{lj} G_{lj} p_j+n_{l})$for all $l$ where $R_{lj}, \forall l, j$ are random variables that model fading, and models general power constraint set, e.g., a single total power constraint.


## Analytical Solution
Under the Rayleigh fading model, the above nonconvex stochastic program can be simplified because the outage probability (please see Kandukuri and Boyd TWC 2002 for more details) of the Ith receiver can be given analytically by :
$$
P\left({SINR}_l(\mathbf{p})<\beta_l\right)=1-e^{\frac{-y \beta}{p}} \prod_{j \neq l}\left(1+\frac{\beta_l F_{l j} p_j}{p_l}\right)^{-1}
$$
where
$$
\begin{aligned}
& F_{l j}= \begin{cases}0, & l=j \\
G_{\mid j} / G_{\|}, & I \neq j .\end{cases} \\
& \mathbf{v}=\left(\frac{n_l}{G_{\|}}, \cdots, \frac{n_L}{G_{L L}}\right) .
\end{aligned}
$$

Next, we give an analytical solution by applying nonnegative matrix theory and nonlinear Perron-Frobenius theory. For illustration, consider the single total power constraint, then the optimal value and solution of (*) are, respectively, given as follows:
$$
\begin{gathered}
1-\mathrm{e}^{-\rho\left(\mathbf{B}\left(\mathbf{p}^*\right)+\frac{1}{\mathrm{p}} \mathbf{v} \mathbf{1}^{\top}\right)} \\
\mathbf{p}^*=\mathbf{x}\left(\mathbf{B}\left(\mathbf{p}^*\right)+\frac{1}{\mathrm{p}} \mathbf{v} \mathbf{1}^{\top}\right)
\end{gathered}
$$


where $\mathbf{x}(\cdot)$ is the right eigenvector corresponding to the Perron-Frobenius eigenvalue $\rho(\cdot)$, and we define
$$
B_{l j}= \begin{cases}0, & l=j \\ \frac{p_l}{p_j} \log \left(1+\frac{\beta_l F_{l j} p_j}{p_l}\right), & l \neq j .\end{cases}
$$

Observe that the spectrum of $\mathbf{B}$ and its rank-one perturbation capture the optimality entirely. For details of the proof and general idea, please refer to [INFOCOM 2011](http://www.cs.cityu.edu.hk/~cheewtan/Tan-INFOCOM2011.pdf). Interestingly, this nonlinear Perron-Frobenius theory approach solves an open problem in [Kandukuri and Boyd TWC 2002](http://www.stanford.edu/~boyd/papers/outage.html) for the interference-limited special case.

## The Algorithm
Using the nonlinear Perron-Frobenius theory, an optimal algorithm is given below to solve the stochastic program (for details: see INFOCOM 2011):
1) Update Power $\mathbf{p}(\mathrm{k}+1)$ :
$$
p_1(k+1)=-\log P\left({SINR}_1(\mathbf{p}(k))>\beta_l\right) p_1(k) \quad \forall I .
$$
2) Nomalize Power $\mathbf{p}(k+1)$ :
$$
\begin{aligned}
& \mathbf{P}(\mathrm{k}+1) \leftarrow \frac{\mathbf{p}(\mathrm{k}+1) \cdot \mathrm{p}}{\mathbf{1}^{\top} \mathbf{p}(\mathrm{k}+1)} \quad \text { if } \quad=\left\{\mathbf{p} \mid \mathbf{1}^{\top} \mathbf{p} \leqslant p\right\} . \\
& \mathbf{P}(\mathrm{k}+1) \leftarrow \frac{\mathbf{p}(\mathrm{k}+1) \cdot \mathrm{p}}{\max _{\mathrm{j}} \mathrm{p}_{\mathrm{j}}(\mathrm{k}+1)} \quad \text { if } \quad=\left\{\mathbf{p}\left|\mathrm{p}_{\mathrm{I}} \leqslant \mathrm{p} \quad \forall\right|\right\} .
\end{aligned}
$$

## The MATLAB Code
Below is an example for the stochastic problem with a single total power constraint:\\

%======================

G = [3.1929    0.1360    0.2379  0.3;
    0.0702    2.8835    0.2436   0.3;
    0.1702    0.8835    2.4436   0.3;
    0.0693    0.0924    0.3060   2.3];
    
n = [0.05;0.05;0.05;0.05];

beta = [1;2;2;2];

pmax = 4;

[p,power_evolution]=worst_outage_prob_min(G,n,beta,pmax)

plot(power_evolution,'-*','linewidth',1.5)

legend('User 1','User 2','User 3','User 4');

%======================


# Max-min Weighted SINR Optimization : Analytical solution and Algorithm
The analytical solution and a fast algorithm for the max-min weighted SINR optimization presented here has its roots in a series of work published in a 2013 IEEE/ACM Transactions on Networking paper and a 2011 IEEE Transactions on Signal Processing paper. An overview survey is available in Wireless Network Optimization by Perron-Frobenius Theory.

The Problem Statement
Maximizing the minimum weighted signal-to-interference-and-noise radio (SINR) under the total power constraint is formulated as follows :

maximize $\min _l \frac{{SINR}_l(\mathbf{p})}{\beta_{\mid}}$

subject to $\mathbf{1}^{\top} \mathbf{p} \leqslant P, \mathbf{p} \geqslant \mathbf{0}$

variables: p.


where ${SINR}_{\mathbf{I}}(\mathbf{p})=G_{\|} p_l /\left(\Sigma_{j \neq l} G_{\mid j} p_j+n_l\right)$ for all $I$, and $\boldsymbol{\beta}=\left(\beta_1, \ldots, \beta_L\right)^{\top} \geqslant 0$ is a given weight vector to reflect priority among users (larger weight means higher priority). A total power budget is given by $P$.

## Analytical Solution
Let us define the following nonnegative matrix :
$$
\mathbf{B}=\mathbf{F}+(1 / P) \mathbf{1 1}^{\top}
$$
and denote
$$
\begin{aligned}
& \mathbf{v}=\left(\frac{1}{G_{11}}, \ldots, \frac{1}{G_{L L}}\right) \\
& F_{l j}= \begin{cases}0, & l=j \\
G_{l j}, & I \neq j .\end{cases}
\end{aligned}
$$

The optimal value and solution of $(*)$ are given, respectively, by
$$
\gamma^*=\frac{1}{\rho({diag}(\boldsymbol{\beta} \cdot \mathbf{v}) \mathbf{B})}
$$
and
$$
\mathbf{P}^*=\left(P / \mathbf{1}^{\top} \mathbf{x}({diag}(\boldsymbol{\beta} \circ \mathbf{v}) \mathbf{B})\right) \mathbf{x}({diag}(\boldsymbol{\beta} \circ \mathbf{v}) \mathbf{B}),
$$
where $\circ$ denotes Schur product and $\mathbf{x}(\cdot)$ denotes the right eigenvector corresponding to the Perron-Frobenius eigenvalue $\rho(\cdot)$.
## A Short Proof Using the Classical Linear Perron-Frobenius Theorem
It can be shown that solving the optimization problem in $(*)$ is equivalent to solving the following fixed-point equation:
$$
\frac{1}{\mathbf{Y}^*} \mathbf{p}^*={diag}(\boldsymbol{\beta} \circ \mathbf{v})\left(\mathbf{F} \mathbf{p}^*+\mathbf{1}\right), \quad \mathbf{1}^{\top} \mathbf{p}^*=\mathrm{P}
$$

A Short Proof Using the Classical Linear Perron-Frobenius Theorem
It can be shown that solving the optimization problem in $(*)$ is equivalent to solving the following fixed-point equation:
$$
\frac{1}{\bigvee^*} \mathbf{p}^*={diag}(\boldsymbol{\beta} \circ \mathbf{v})\left(\mathbf{F} \mathbf{p}^*+\mathbf{1}\right), \quad \mathbf{1}^{\top} \mathbf{p}^*=\mathrm{P}
$$

Now, observe that :
$$
\frac{1}{\gamma^*} \mathbf{p}^*={diag}(\boldsymbol{\beta} \circ \mathbf{v})\left(\mathbf{F}+\frac{1}{\mathrm{P}} \mathbf{1 1}^{\top}\right) \mathbf{p}^*
$$

therefore (∗) can be solved analytically as an eigenvalue problem by the classical linear Perron-Frobenius theorem.

For a different proof, e.g., the nonlinear Perron−Frobenius theory
, please see IEEE/ACM Transactions on Networking in 2013 and IEEE Transactions on Signal Processing in 2011.
