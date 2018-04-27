#  Modular arbitrary-order ocean-atmosphere model: The Tangent Linear and Adjoint model #

## Description : ##

 The Tangent Linear and Adjoint model model are implemented in the same way as the nonlinear model, with a tensor storing the different terms. The Tangent Linear (TL) tensor \f$\mathcal{T}_{i,j,k}^{TD}\f$ is defined as:

\f[ \mathcal{T}_{i,j,k}^{TL} = \mathcal{T}_{i,k,j} + \mathcal{T}_{i,j,k} \f]

while the Adjoint (AD) tensor \f$\mathcal{T}_{i,j,k}^{AD}\f$ is defined as:

\f[ \mathcal{T}_{i,j,k}^{AD} = \mathcal{T}_{j,k,i} + \mathcal{T}_{j,i,k} . \f]

where \f$ \mathcal{T}_{i,j,k}\f$ is the tensor of the nonlinear model.

These two tensors are used to compute the trajectories of the models, with the equations

\f[  \frac{d\delta y_i}{dt} = \sum_{j=1}^{ndim}\sum_{k=0}^{ndim} \, \mathcal{T}_{i,j,k}^{TL} \, y^{\ast}_k \; \delta y_j . \f]

\f[   -\frac{d\delta y_i}{dt} = \sum_{j=1}^{ndim} \sum_{k=0}^{ndim} \, \mathcal{T}_{i,j,k}^{AD} \, y^{\ast}_k \; \delta y_j . \f]

where \f$\boldsymbol{y}^{\ast}\f$ is the point where the Tangent model is defined (with \f$y_0^{\ast}=1\f$).

## Implementation : ##

The two tensors are implemented in the module tl_ad_tensor and must be initialized (after calling params::init_params and aotensor_def::aotensor) by calling tl_ad_tensor::init_tltensor() and tl_ad_tensor::init_adtensor(). The tendencies are then given by the routine tl(t,ystar,deltay,buf) and ad(t,ystar,deltay,buf). An integrator with the Heun method is available in the module rk2_tl_ad_integrator and a fourth-order Runge-Kutta integrator in rk4_tl_ad_integrator. An example on how to use it can be found in the test file test_tl_ad.f90
