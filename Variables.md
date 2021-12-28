# ParMooN SPADE

# Generation of Monte Carlo realizations 
## Variables 
### From Database (changes to ```Database.h```, ```Database.C``` and ```ReadParam.C```)
```N_Realisations```  - Number of Realizations, read from ReadIn file \
``` LengthScale```  - Length Scale for correlation function, read from ReadIn File \
``` Eigen Precent``` - Cut-off percentage for singular values, read from ReadIn file \
``` E ``` - Scaling factor for standard deviation (bump)function\
``` disp ``` - Displacement of center of standard deviation (bump) function 
``` power ``` - Power of standard deviation (bump) function (flatness of bump) \

### Declared in code
```N_U``` - Velocity degree of freedom 
```N_P``` - Pressure degree of freedom 


### Calculated 
``` C ``` - Covariance Matrix, size ```N_U x N_U```, stored in **row major format** \
``` C1 ``` - Correlation Matrix, size ```N_U x N_U```, stored in **row major format** \
``` C = U x S x Vt``` is the singular value decomposition (LAPACK, row major form) of covariance matrix with,\ 

```S``` consists of ```N_U``` singular values \
```U``` is the left singular matrix, dimension ```N_U x N_U```, in row major form & \
```Vt``` is the right singular matrix, dimension ```N_U x N_U```, in row major form 

Considering only the ```modDim``` singular values which retain ```SVPercent``` or more of the energy, we get the truncated singular matrix ```Ut``` with dimension ```N_U x modDim```

```Z``` has dimension ```modDim x N_Realisations```, each row ```k``` has a normal distribution with standard deviation ```S[k]```



```RealizationVector``` is product of ```Ut``` and ```Z``` with dimension ```N_U x N_Realisations```, i.e. each column is a realization for ```N_U``` degress of freedom

# DO Initialization 
## Variables 
### From Database (changes to ```Database.h```, ```Database.C``` and ```ReadParam.C```)
```SVPercent``` - Cut-off percentage for singular values of perturbation matrix ```PerturbationVector```, read from ReadIn file 

### Calculated
```MeanVector``` = Mean of all realizations for a given DOF, calculated from ```RealizationVector```. Since ```RealizationVector``` is in row major format, mean for ith DOF is mean of ith block of length ```N_Realisations```.\
 ```MeanVector``` is stored in __row major format__. Size ```N_U x 1```

```PerturbationVector``` = Matrix of perturbation i.e. mean removed from realization matrix. Subtract ```MeanVector``` from each column of ```RealizationVector```.
```PerturbationVector``` is stored in __row major format__. Size ```N_U x N_Realisations```

