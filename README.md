# Wireless-Network-Optimization

\begin{equation} \begin{array}[c]{rl} \mbox{maximize} & \displaystyle\min_{l}\frac{{\mathbf{SINR}}_l(\mathbf{p})}{\beta_l}\\ \mbox{subject to} & \mathbf{1}^\top\mathbf{p}\le\bar{P},\mathbf{p}\geq\mathbf{0}\\ \mbox{variables:} & \mathbf{p}. \end{array} \tag{*} \end{equation}
where $\mathbf{SINR_l(\mathbf{p})}=G_{ll}p_l/(\sum_{j\ne l} G_{lj}p_j+n_l)$ for all $l$, and $\boldsymbol{\beta}=(\beta_{1},\dots,\beta_{L})^\top \ge 0$ is a given weight vector to reflect priority among users (larger weight means higher priority). A total power budget is given by $\bar{P}$.
