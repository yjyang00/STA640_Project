# STA-640
640 Project

Authors: Yicheng Shen, Yanjiao Yang, Huiying Lin

Do fixed/random effects account for unmeasured confounding? 

<img width="715" alt="Screenshot 2023-04-19 at 9 58 08 PM" src="https://user-images.githubusercontent.com/67173948/233238251-43b5dc0b-5c09-47bf-966e-ed30e645d135.png">


A common saying in econometrics is that fixed/random effects in analyzing panel data absorb unmeasured confounding. In other words, with repeated measurements of the same unit, a model with unit-specific fixed/random effect can bypass unmeasured confounding. Some discussion can be found in Angrist and Pischke (2009, Mostly Harmless Econometrics, Chapter 5). A recent bold example is Dee et al. (2023, Nature Communications). But this claim is a mystery, why and in what sense? Hazlett and Wainstein (2020, Political Analysis) has some useful discussion.


$$
\begin{align*}
E\left(\mathrm{Y}_{i t} \mid A_i, \mathrm{X}_{i t}, t, \mathrm{D}_{i t}\right)& =\alpha+\lambda_t+\rho \mathrm{D}_{i t}+A_i^{\prime} \gamma+\mathrm{X}_{i t} \delta \\
Y_{it} &= \alpha_i+\lambda_t+\rho \mathrm{D}_{i t}+A_i^{\prime} \gamma+\mathrm{X}_{i t} \delta + \epsilon_{it} 
\end{align*}
$$


<img width="744" alt="Screenshot 2023-04-19 at 9 54 11 PM" src="https://user-images.githubusercontent.com/67173948/233237710-33e78c3d-2184-4224-847b-c75a2dd563ac.png">



