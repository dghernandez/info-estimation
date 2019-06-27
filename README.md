# info-estimation
Estimating the mutual information between discrete variables with limited samples. This is the code for the paper https://www.mdpi.com/1099-4300/21/6/623.


![](https://github.com/dghernandez/info-estimation/blob/master/scheme2.jpg)


The different functions are embedded in a Mathematica Package :"InfoHDP.m". This file needs to be placed in your local Package folder. Once this is done, you can load the Package to any notebook using:
```
<<InfoHDP.m`
```
The application of our method in a simple example can be found in the notebook example0.nb.

# Core calculations and conditions
The main parts of the code are rather simple and they can be implemented in any programming language: (1) a maximization of the marginal log-likelihood (or posterior) over the hyperparameter beta, and (2) the evaluation of the posterior information in such beta. 

First, we need to obtain the ML (or MAP) estimate for the hyperparameter beta (from Eq. (14) in the paper, or Eq. (21) for the symmetric binary case). I recommend to do this maximization with the argument log(beta), and explore the interval log(0.01)<log(beta)<log(100), which is usually enough. Only terms with coincidences in the large entropy variable (n_x>1) need to be considered, as the others add a constant term in beta. In an undersampled regime, there would be many repeated terms and they can be grouped together for a more efficient evaluation (using multiplicities). Secondly, we evaluate the posterior information (see Eq. (16) in the paper, or Eq. (20) for the symmetric binary case) in the beta found previously. In this evaluation, all occupied states (n_x>0) need to be included.

Our method needs coincidences on the large entropy variable, which starts to happen when N> exp(H(X)/2). If there are no coincidences, then the marginal likelihood is flat on beta. If there are few coincidences and no prior on beta is used, the maximum may be attain for beta tending to zero or infinity. For such cases the posterior information is still well-defined and takes the values of H(Y) or zero, respectively.
