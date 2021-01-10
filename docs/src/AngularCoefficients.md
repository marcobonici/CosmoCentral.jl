# Angular Coefficients
Central quantities are the Angular Coefficients ``C_\ell``. Actually we
implement only the Limber approximation to evaluate the ``C_\ell``,
[according to](https://arxiv.org/abs/1910.09273):

```math
C_{i j}^{AB}(\ell)=\frac{c}{H_0} \int_{z_{\min }}^{z_{\max }} \mathrm{d} z \frac{W_{i}^{A}(z) W_{j}^{B}(z)}{E(z) r^{2}(z)} P_{\delta \delta}\left(\frac{\ell+1 / 2}{r(z)}, z\right)
```

```@docs
CosmoCentral.AngularCoefficientsStruct
CosmoCentral.ComputeAngularCoefficients
```
