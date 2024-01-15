# Wireless-Network-Optimization


## The Problem Statement
An outage event occurs at the Ith receiver when the received SINR falls below a given reliability threshold, i.e., $\operatorname{SINR}_I(\mathbf{p})<\beta_I$ for $I=1, \ldots, L$. So we are interested to minimize the worst-case outage probability to ensure reliability fairness, which is formulated as follows :
$\operatorname{minimize} \max _1 P\left(\operatorname{SINR}_1(\mathbf{p})<\beta_l\right)$
subject to $p \in$
variables: p.
where $\operatorname{SINR}_{\mid}(\mathbf{p})=R_{\|} G_{\|} p_l /\left(\sum_{j \neq l} R_{\mid j} G_{\mid j} p_j+n_{\mid}\right)$for all I where $R_{\mid j}, \forall l, j$ are random variables that model fading, and models general power constraint set, e.g., a single total power constraint.
