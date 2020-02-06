## Calculation of risk

For more detailied explaination of calculation of risk refer to : https://arxiv.org/pdf/1606.00931.pdf


In linear cox regression :
$$ risk(i) = \beta^Tx_i $$

In NN cox regression:
$$ r_i = risk(i) = NN(x_i) $$
And let $$p_i  = exp(r_i)$$
The log likelihood of risk is defined as :

The patients are sorted based on descending events
\begin{align}
log \prod_{i=n:event=death}^1 \frac{p_i}{\sum_{j\ge i} p_j} \\
= \sum_{i=n}^1 log(p_i) - log(cumsum(p_i)) \\
= \sum_{i=n}^1 r_i - log(cumsum(p_i)) \\
= \sum_{i=n:all~events}^1 ( r_i - log(cumsum(p_i))) * Event(i) \\
\end{align}

The above assume $$Event(i) =1 $$ is death and  0 is censor.

## Calculation of uncertainity 

Similar to:  https://arxiv.org/abs/1705.07115
Here we make intuitive argument that 

Because the log likelihood is very similar form to softmax function, we use same constant i.e. (refer to equation 10)
\begin{align}
\Lambda(W, \sigma) = 1/\sigma^2 \Lambda(W) + log \sigma
\end{align}
