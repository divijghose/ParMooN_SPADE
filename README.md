# ParMooN SPADE

# Generation of Monte Carlo realizations 
## Variables 
### From Database (changes to ```Database.h```)
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

