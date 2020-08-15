## @package NereidesExceptions
# Module for raising custom exceptions for using with Nereides
# Least Squares Spectral Element flow solver.

## Class for mesh generation erros, inherits from python Exception class
class MeshGenerationError:
    
    ## The constructor fpr the MeshGenerationError object
    # \param[in] value The content that defines the kind of error, which will
    #                  be given to the user.
    # \param[out] self.parameter Where the value of the exception will be 
    #                            stored for the user to access it.
    def __init__(self, value):
        self.parameter = value
        
    def __str__(self):
        return repr(self.parameter)


## Class for Gauss-Quadrature generation errors, inherits from python Exception class
class GaussQuadratureError:
    
    ## The constructor fpr the GaussQuadError object
    # \param[in] value The content that defines the kind of error, which will
    #                  be given to the user.
    # \param[out] self.parameter Where the value of the exception will be 
    #                            stored for the user to access it.
    def __init__(self, value):
        self.parameter = value
        
    def __str__(self):
        return repr(self.parameter)
    
## Class for GlobalNumbering generation errors, inherits from python Exception class
class GlobalNumberingError:
    
    ## The constructor fpr the GlobalNumberingError object
    # \param[in] value The content that defines the kind of error, which will
    #                  be given to the user.
    # \param[out] self.parameter Where the value of the exception will be 
    #                            stored for the user to access it.
    def __init__(self, value):
        self.parameter = value
        
    def __str__(self):
        return repr(self.parameter)