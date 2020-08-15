import numpy


## \class variable
# \brief Contains all data relative to the definition of variables and operations
# that can be performed on variables. These include, variable name (ex.: pressure); 
# symbol used to define the variable (ex.: p); order (ex.: -1) this order is
# relative to the main order defined, that is, in our example it is -1 that means
# that if we are using main order p=10, then this variable will be computed
# using order p_pressure = p+(-1) = 10+(-1) = 9; q_type (quadrature type), defines
# which type of collocation points should be used for this variable (ex.: gauss).
#
class variable:
    # variable class contains the following variables:
    
    ## \var name 
    # \brief The full name of the variable.
    #
    # Can be a string of any length.\n\n
    # \b TYPE: \c string
    
    ## \var symbol 
    # \brief Small group of letters, in general only one, to use to identify the
    # variable.
    #
    # Can be a string of any length, but should be as small as possible.\n\n
    # \b TYPE: \c string
    
    ## \var order 
    # \brief The relative order to be used to compute the variable. If main
    # order is p the order used for this variable will be p+order.\n\n
    # \b TYPE: \c int16
    
    ## \var q_type 
    # \brief The quadrature type used to compute this variable. 
    # Can be any one of settings.Q_TYPES .\n\n
    # \b TYPE: \c string
    
    ## \brief The constructor of the variable object, generates the %variable.
    #
    # Depending on the definitiion of the variable
    def __init__(name, notation, order, q_type):
        print o