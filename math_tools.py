import numpy
import scipy
import sympy
from nereides import _math_tools

class poly2d:
    def __init__(self, **parameters):
        # parameters can be:
        # coefficients: a 2d matrix with the coefficients of the 2d
        #               polynomial, lower indices correspond to higherpowers
        #               first dimension, rows, is y variable, second dimension,
        #               columns, is x variable
        #
        # polyX: a poly1d polynomial which is the x variable
        #
        # polyY: a poly1d polynomial which is the y variable
        #
        # the 2d polynomial will be the product of the two.
        
        # part for parameters with coefficients
        if parameters.has_key("coeffs"):
            coeffs = parameters["coeffs"]
        
        # part for parameters with polynomial in x and polynomial in y
        if (parameters.has_key("polyX") and parameters.has_key("polyY")):
            if isinstance(parameters["polyX"], scipy.poly1d):
                polyX = parameters["polyX"]
            else:
                polyX = scipy.poly1d(parameters["polyX"])
                
            if isinstance(parameters["polyY"], scipy.poly1d):
                polyY = parameters["polyY"]
            else:
                polyY = scipy.poly1d(parameters["polyY"])
            
            # multiply both polynomials to get a 2d polynomial
            coeffY,coeffX =  numpy.meshgrid(polyY.coeffs,
                                            polyX.coeffs)
                                            
            coeffs = coeffY * coeffX
            
        # Remove the extra zeros
        if coeffs.shape != (1,1):
            # Remove first the 0 coefficients on the higher order of the variables
            # first on y, the columns
            columns_to_remove = 0
            for column in coeffs.T:
                if (column == 0.).all():
                    columns_to_remove += 1
                else:
                    # get out of the for loop as soon as one column is not 0.
                    break
            
            # Now the x, the rows
            rows_to_remove = 0
            for row in coeffs:
                if (row == 0.).all():
                    rows_to_remove += 1
                else:
                    # get out of the for loop as soon as one row is not 0.
                    break
        
            # remove the first rows and columns that are zeros
            coeffs = coeffs[rows_to_remove:, columns_to_remove:]
            if coeffs.shape == (0,):
                coeffs = numpy.array([[0]])
        
        # store the coefficients
        self.__dict__['coeffs'] = coeffs
        
        # compute the order and store it
        self.__dict__['order'] = (coeffs.shape[0]-1,coeffs.shape[1]-1)
        
        # define the variables
        self.__dict__['variable'] = ('x', 'y')
        
    def __array__(self, t=None):
        if t:
            return numpy.asarray(self.coeffs, t)
        else:
            return numpy.asarray(self.coeffs)
        
    def __repr__(self):
        vals = repr(self.coeffs)
        vals = vals[6:-1]
        return "poly2d(%s)" % vals
    
    def __len__(self):
        return self.order
    
    def __call__(self, x, y, dx_dy=None):
        return poly2deval(self, x, y, dx_dy)
    
    def __getitem__(self, val):
        ind = [self.order[0] - val[0], self.order[1] - val[1]]
        if (ind[0]<0) or (ind[1]<0):
            return 0
        if (val[0]<0) or (val[1]<0):
            return 0
        return self.coeffs[ind[0],ind[1]]
        
    def deriv(self, dx_dy):
        return poly2dderiv(self, dx_dy)
        
## \class TransfiniteMapping_2d 
# \brief Implements a transfinite mapping for transforming local coordinates
#  (u,v) onto global coordinates (x,y).
#
# In order to have higher accuracy in nonpolygonal domains, the elements used
# for the partitioning of the domain should have curved edges. In order to
# circumvent the problems that arise with curved parametric elements, such as
# isoparametric, one implements a transfinite mapping. See paper from Gordon
# and Hall (http://www.springerlink.com/content/p5517206551228j6/fulltext.pdf)
# for more detailed information. \n
#
# This class implements, given the definition of the edges of an element, using
# straigh lines (implemented), circles (to implement) or NURBS (to implement),
# compute the transfinite mapping T from the local coordinates (u,v) onto the
# global coordinates:
# 
# \f[
#     \mathbf{T}:(u,v)\rightarrow (x,y)
# \f]
#
# with:
#
# \f[
#    \mathbf{T}(u,v) = \left[
#                         \begin{array}{c}
#                            X(u,v) \\
#                            Y(u,v)
#                         \end{array}
#                      \right]
# \f]
#
# For the mapping one uses the vector valued bilinearly blended map given 
# for an arbitrary local domain \f$ S = [u_{s},u_{e}]\times[v_{s},v_{e}]\f$ by:
#
# \f[
#    \mathbf{T}(u,v) = \left[
#                         \begin{array}{c}
#                            X(u,v) \\
#                            Y(u,v)
#                         \end{array}
#                      \right] = \left(1-\frac{u-u_{s}}{u_{e}-u_{s}}\right)\mathbf{\Gamma}(u_{s},v) +
#                                \frac{u-u_{s}}{u_{e}-u_{s}}\mathbf{\Gamma}(u_{e},v) + 
#                                \left(1-\frac{v-v_{s}}{v_{e}-v_{s}}\right)\mathbf{\Gamma}(u,v_{s}) +
#                                \frac{v-v_{s}}{v_{e}-v_{s}}\mathbf{\Gamma}(u,v_{e}) - 
#                                \left(1-\frac{u-u_{s}}{u_{e}-u_{s}}\right)\left(1-\frac{v-v_{s}}{v_{e}-v_{s}}\right)\mathbf{\Gamma}(u_{s},v_{s}) -
#                                \left(1-\frac{u-u_{s}}{u_{e}-u_{s}}\right)\frac{v-v_{s}}{v_{e}-v_{s}}\mathbf{\Gamma}(u_{s},v_{e}) -
#                                \left(1-\frac{v-v_{s}}{v_{e}-v_{s}}\right)\frac{u-u_{s}}{u_{e}-u_{s}}\mathbf{\Gamma}(u_{e},v_{s}) - 
#                                \frac{v-v_{s}}{v_{e}-v_{s}}\frac{u-u_{s}}{u_{e}-u_{s}}\mathbf{\Gamma}(u_{e},v_{e})
# \f]
#
# where \f$ \mathbf{\Gamma} \f$ is the parametrization of the boundary of the
# element and \f$ \mathbf{\Gamma}(u,{v_{s}}) \f$ is the parametrization of
# the bottom edge, \f$ \mathbf{\Gamma}(u_{e},{v}) \f$ is the parametrization of
# the right edge, \f$ \mathbf{\Gamma}(u,{v_{e}}) \f$ is the parametrization of
# the top edge and \f$ \mathbf{\Gamma}(u_{s},{v}) \f$ is the parametrization of
# the left edge.
#
# For this particular case our local space is: \f$ S = [-1,1]\times[-1,1]\f$
# hence the transfinite mapping resumes to:
#
# \f[
#    \mathbf{T}(u,v) = \left[
#                         \begin{array}{c}
#                            X(u,v) \\
#                            Y(u,v)
#                         \end{array}
#                      \right] = \left(1-\frac{u+1}{2}\right) \mathbf{\Gamma}(-1,v) +
#                                \frac{u+1}{2}\mathbf{\Gamma}(1,v) + 
#                                \left(1-\frac{v+1}{2}\right) \mathbf{\Gamma}(u,-1) +
#                                \frac{v+1}{2}\mathbf{\Gamma}(u,1) - 
#                                \left(1-\frac{u+1}{2}\right) \left(1-\frac{v+1}{2}\right) \mathbf{\Gamma}(-1,-1) -
#                                \left(1-\frac{u+1}{2}\right) \frac{v+1}{2}\mathbf{\Gamma}(-1,1) -
#                                \left(1-\frac{v+1}{2}\right) \frac{u+1}{2}\mathbf{\Gamma}(1,-1) - 
#                                \frac{v+1}{2}\frac{u+1}{2}\mathbf{\Gamma}(1,1)
# \f]
#
# For usage in the integration and presentation of results, one needs:
#    \li the mapping : \f$ X(u,v)\f$ and \f$Y(u,v)\f$.
#    \li the partial derivatives of the mapping : \f$\frac{\partial X(u,v)}{\partial u}\f$,
#        \f$\frac{\partial X(u,v)}{\partial v}\f$, \f$\frac{\partial Y(u,v)}{\partial u}\f$ and
#        \f$\frac{\partial Y(u,v)}{\partial v}\f$.
#    \li the Jacobian determinant of the mapping : \f$\det [J(u,v)] = 
#       \frac{\partial X(u,v)}{\partial u}\frac{\partial Y(u,v)}{\partial v} -
#       \frac{\partial Y(u,v)}{\partial u}\frac{\partial X(u,v)}{\partial v}\f$
#
# \n Hence, when instantiated, this class will create all this functions, as
# regular python functions that can be evaluated in numpy.array arrays, for
# fast evaluation. \n
#
# Since one will use straight lines, circular arcs or NURBS for the
# specification of the boundaries, and correspondingly the element's edges,
# one needs to compute the derivatives for all this types of curves. This would
# be troublesome and tedious. A different way is followed here. The usage of
# sympy is made. Sympy is a python library for symbolic mathematics, allowing
# us to easily define our functions and then compute their derivatives, multiply
# the derivatives with in order to generate the Jacobian and simplifying them
# in order to have fast evaluatable functions. For more information on
# sympy look at:
# \li main page : http://code.google.com/p/sympy/
# \li documentation page : http://docs.sympy.org/
# \li new wiki page : http://wiki.sympy.org/wiki/Main_Page
# \li interactive web version : http://live.sympy.org/
# 
# \todo <b> FOR NOW ONLY STRAIGHT LINES ARE IMPLEMENTED, NEXT WORK IS EXTENDING TO
#    CIRCLE ARCS AND NURBS, WORK NEEDS TO BE DONE IN ORDER TO GET THIS INFO
#    FROM THE MESH GENERATOR. </b>
#
class TransfiniteMapping_2d:
    ## \var element
    # \brief Holds the information on the type of element.
    #
    #    \li 3 : triangular straight element
    #    \li 4 : quadrilateral straight element
    #
    
    ## \var X
    # \brief Y evaluatable function for the transfinite mapping.
    #
    # It is a 2 variable inputs (u and v) function of the type sympy.lambdify
    # that takes as inputs numpy.array arrays and returns the values of the
    # function X of the transfinite mapping at all the input pairs. It returns an object equal to the
    # ones given as input. That is, if u and v are numpy.array of type
    # float64 of size 2 x 3, the returned array will be a numpy.array of type
    # float64 of size 2 x 3. 
    #
    
    ## \var Y
    # \brief Y evaluatable function for the transfinite mapping.
    #
    # It is a 2 variable inputs (u and v) function of the type sympy.lambdify
    # that takes as inputs numpy.array arrays and returns the values of the
    # function Y of the transfinite mapping at all the input pairs. It returns an object equal to the
    # ones given as input. That is, if u and v are numpy.array of type
    # float64 of size 2 x 3, the returned array will be a numpy.array of type
    # float64 of size 2 x 3. 
    #
    
    ## \var X_func
    # \brief X analytical function for the transfinite mapping.
    #
    # It is an object of the type sympy.function that represents the X(u,v) function
    # of the transfinite mapping. This one cannot be evaluated
    # at a numpy.array. Instead one can see how the function for the mapping
    # is defined analytically. This is also used for the differentiation.
    #
    
    ## \var Y_func
    # \brief Y analytical function for the transfinite mapping.
    #
    # It is an object of the type sympy.function that represents the Y(u,v) function
    # of the transfinite mapping. This one cannot be evaluated
    # at a numpy.array. Instead one can see how the function for the mapping
    # is defined analytically. This is also used for the differentiation.
    #
    
    ## \var dX_du
    # \brief \f$\frac{\partial X}{\partial u}\f$ evaluatable function for the transfinite mapping.
    #
    # It is a 2 variable inputs (u and v) function of the type sympy.lambdify
    # that takes as inputs numpy.array arrays and returns the values of the
    # function \f$\frac{partial X}{\partial u}\f$ of the transfinite mapping at all the input pairs. It returns an object equal to the
    # ones given as input. That is, if u and v are numpy.array of type
    # float64 of size 2 x 3, the returned array will be a numpy.array of type
    # float64 of size 2 x 3. 
    #
    
    ## \var dX_dv
    # \brief \f$\frac{\partial X}{\partial v}\f$ evaluatable function for the transfinite mapping.
    #
    # It is a 2 variable inputs (u and v) function of the type sympy.lambdify
    # that takes as inputs numpy.array arrays and returns the values of the
    # function \f$\frac{partial X}{\partial v}\f$ of the transfinite mapping at all the input pairs. It returns an object equal to the
    # ones given as input. That is, if u and v are numpy.array of type
    # float64 of size 2 x 3, the returned array will be a numpy.array of type
    # float64 of size 2 x 3.
    #
    
    ## \var dY_du
    # \brief \f$\frac{\partial Y}{\partial u}\f$ evaluatable function for the transfinite mapping.
    #
    # It is a 2 variable inputs (u and v) function of the type sympy.lambdify
    # that takes as inputs numpy.array arrays and returns the values of the
    # function \f$\frac{partial Y}{\partial u}\f$ of the transfinite mapping at all the input pairs. It returns an object equal to the
    # ones given as input. That is, if u and v are numpy.array of type
    # float64 of size 2 x 3, the returned array will be a numpy.array of type
    # float64 of size 2 x 3.
    #
    
    ## \var dY_dv
    # \brief \f$\frac{\partial Y}{\partial v}\f$ evaluatable function for the transfinite mapping.
    #
    # It is a 2 variable inputs (u and v) function of the type sympy.lambdify
    # that takes as inputs numpy.array arrays and returns the values of the
    # function \f$\frac{partial Y}{\partial v}\f$ of the transfinite mapping at all the input pairs. It returns an object equal to the
    # ones given as input. That is, if u and v are numpy.array of type
    # float64 of size 2 x 3, the returned array will be a numpy.array of type
    # float64 of size 2 x 3.
    #
    
    ## \var J
    # \brief Jacobian determinant evaluatable function for the transfinite mapping.
    #
    # It is a 2 variable inputs (u and v) function of the type sympy.lambdify
    # that takes as inputs numpy.array arrays and returns the values of the
    # function \f$\det [J(u,v)] = 
    #       \frac{\partial X(u,v)}{\partial u}\frac{\partial Y(u,v)}{\partial v} -
    #       \frac{\partial Y(u,v)}{\partial u}\frac{\partial X(u,v)}{\partial v}\f$ of the transfinite mapping at all the input pairs. It returns an object equal to the
    # ones given as input. That is, if u and v are numpy.array of type
    # float64 of size 2 x 3, the returned array will be a numpy.array of type
    # float64 of size 2 x 3.
    #
    
    ## \brief Constructor for the TransfiniteMapping_2d object, generated the
    # transfinite mapping on an element given the element's definition.
    #
    # Upton instantiation, the object generates:
    #    \li the mapping : \f$ X(u,v)\f$ and \f$Y(u,v)\f$.
    #    \li the partial derivatives of the mapping : \f$\frac{\partial X(u,v)}{\partial u}\f$,
    #        \f$\frac{\partial X(u,v)}{\partial v}\f$, \f$\frac{\partial Y(u,v)}{\partial u}\f$ and
    #        \f$\frac{\partial Y(u,v)}{\partial v}\f$.
    #    \li the Jacobian determinant of the mapping : \f$\det [J(u,v)] = 
    #       \frac{\partial X(u,v)}{\partial u}\frac{\partial Y(u,v)}{\partial v} -
    #       \frac{\partial Y(u,v)}{\partial u}\frac{\partial X(u,v)}{\partial v}\f$
    #
    ## \param[in] nodes
    # It is a vector of type numpy.array, that contains the list of nodes
    # that define the element. The nodes should be given counterclockwise
    # following the convention used on the mesh generation of elements2d.nodes
    # for each element.
    #
    ## \param[out] X_func
    # See #X_func
    #
    ## \param[out] Y_func
    # See #Y_func
    #
    ## \param[out] X
    # See #X
    #
    ## \param[out] Y
    # See #Y
    #
    ## \param[out] dX_du
    # See #dX_du
    #
    ## \param[out] dX_dv
    # See #dX_dv
    ## \param[out] dY_du
    # See #dY_du
    #
    ## \param[out] dY_dv
    # See #dY_dv
    #
    ## \param[out] J
    # See #J
    #
    ## \param[out] element
    # See #element
    
    def __init__(self, nodes):
        # nodes is an array with the x and y coordinates of the nodes that
        # define the element. Each row contains the x and y coordinates of each
        # node. Hence nodes, has as many rows as the number of nodes that
        # define the element, i.e., for triangular elements 3 nodes (rows) and
        # for quadrilateral elements 4 nodes (rows). The number of columns is
        # always 2, since we are in 2d.
        
        try:
            # compute the type of element: 3- triangle, 4- quadrilateral
            self.element = nodes.shape[0]
        except AttributeError:
            print "nodes must be a numpy.array of N_nodes x 2 values"
            del(self)
            return
        
        # if all the nodes are  straight lines just compute directly the
        # transfinite mapping
        
        # first define the independent variables u and v as sympy.variables
        # this is needed to define the functions that define the mappings
        u,v = sympy.symbols("uv")
        # start with the function for X(u,v)
        self.X_func = (nodes[3,0] - nodes[0,0])*(v + 1.0)*0.5 +\
                                  (nodes[1,0] - nodes[0,0])*(u + 1.0)*(1.0 - v)*0.25 +\
                                  (nodes[2,0] - nodes[3,0])*(u + 1.0)*(v + 1.0)*0.25 +\
                                  nodes[0,0]
        # now with the function for Y(u,v)
        self.Y_func = (nodes[3,1] - nodes[0,1])*(v + 1.0)*0.5 +\
                                  (nodes[1,1] - nodes[0,1])*(u + 1.0)*(1.0 - v)*0.25 +\
                                  (nodes[2,1] - nodes[3,1])*(u + 1.0)*(v + 1.0)*0.25 +\
                                  nodes[0,1]
                                
        # The sympy functions just generated cannot be evaluated, they are
        # just useful to compute derivatives and other things.
        # For rapid evaluation of these functions in points (u,v), one needs
        # to create what is called a lambdified function evaluatable at 
        # numpy.array objects, which is what we want.
        
        # the X evaluatable mapping
        self.X = sympy.lambdify((u,v), self.X_func, numpy)
        # the Y evaluatable mapping
        self.Y = sympy.lambdify((u,v), self.Y_func, numpy)
        
        # now, since we will mostly use the first derivatives of X and Y, i.e.,
        # dX/du, dX/dv, dY/du and dY/dv, one generates the lamdified versions
        # of these functions right away
        
        # store the mapping derivatives
        dX_du = self.X_func.expand().diff(u)
        dX_dv = self.X_func.expand().diff(v)
        dY_du = self.Y_func.expand().diff(u)
        dY_dv = self.Y_func.expand().diff(v)
        
        # the X evaluatable derivatives
        self.dX_du = sympy.lambdify((u,v), dX_du, numpy)
        self.dX_dv = sympy.lambdify((u,v), dX_dv, numpy)
        # the Y evaluatable derivatives
        self.dY_du = sympy.lambdify((u,v), dY_du, numpy)
        self.dY_dv = sympy.lambdify((u,v), dY_dv, numpy)
        
        # compute the Jacobian determinant, J(u,v), as an evaluatable function
        # recall that J(u,v) = dX/du * dY/dv - dY/du * dX/dv.
        # sympy.expand is called since in this way the resulting expression
        # will be simpler to evaluate.
        self.J = sympy.lambdify((u,v), \
                             sympy.expand(dX_du*dY_dv - dY_du*dX_dv), numpy)
        
        
class basis_function_2d:
    # function that defines a basis function and its derivatives
    
    def __init__(self, coeffs, d_coeffs):
        # coeffs has on the first row the coefficients for the x part,
        # that is phi(x), and the second row has the coefficients for the
        # y part, that is psi_(y).
        # d_coeffs is the same but for the  derivatives in 
        # x and y, respectively.
        
        self.__dict__["coeffs"] = numpy.array(coeffs, dtype="float64")
        self.__dict__["d_coeffs"] = numpy.array(d_coeffs, dtype="float64")
        
        # define the order of the basis function
        self.__dict__["order"] = [self.coeffs[0].shape[0], self.coeffs[1].shape[0]]
        
    def __len__(self):
        # return the order of the basis function
        return self.order
    
    def __call__(self, u, v, du_dv = None):
        p = self.order
        
        # make sure that u and v are of type float64 and ordered in 
        # fortran fashion
        u = numpy.array(u, dtype="float64", order="fortran")
        v = numpy.array(v, dtype="float64", order="fortran")
        
        # allocate memory space for the result array
        result = numpy.zeros_like(u)
    
        # check if du_dv was given
        if du_dv != None:
            # check the size of du_dv, it must be vector of length 2
            try:
                if len(du_dv) != 2:
                    print "du_dy should be a vector of size 2: [du dv]"
                    return
            except TypeError:
                print "du_dv must be a vector of size 2: [du dv]"
                return
            
            # try to convert both values to int
            try:
                du_dv[0] = int(du_dv[0])
                du_dv[1] = int(du_dv[1])
            except ValueError:
                print "du_dv must be a vector with two integers"
                return
                
            # check if they are bigger or equal to zero
            if (du_dv[0]<0) or (du_dv[1]<0):
                print "du_dv must be a vector with two positive values"
        
        if (du_dv == None) or (du_dv == [0, 0]):
            result = _math_tools.polyeval_2d(u, v, self.coeffs[0], self.coeffs[1], result)
            return result
        elif du_dv == [1, 0]:
            result = _math_tools.polyeval_2d(u, v, self.d_coeffs[0], self.coeffs[1], result)
            return result
        elif du_dv == [0, 1]:
            result = _math_tools.polyeval_2d(u, v, self.coeffs[0], self.d_coeffs[1], result)
            return result
        else:
            # compute the coefficients of the derivatives
            du_coeffs = scipy.polyder(self.coeffs[0], m=du_dv[0])
            dv_coeffs = scipy.polyder(self.coeffs[1], m=du_dv[1])
            result = _math_tools.polyeval_2d(u, v, du_coeffs, dv_coeffs, result)
            return result
        
def poly2deval_fortran(x, y, polyXY, dx_dy=None):
    # make sure that x, y and polyXY are float64 and ordered in fortran mode
    x = numpy.array(x, dtype="float64", order="fortran")
    y = numpy.array(y, dtype="float64", order="fortran")
    polyXY = numpy.array(polyXY, dtype="float64", order="fortran")
    
    # check if dx_dy was given
    if dx_dy != None:
        # check the size of dx_dy, it must be vector of length 2
        try:
            if len(dx_dy) != 2:
                print "dx_dy should be a vector of size 2: [dx dy]"
                return
        except TypeError:
            print "dx_dy must be a vector of size 2: [dx dy]"
            return
        
        # try to convert both values to int
        try:
            dx_dv[0] = int(dx_dy[0])
            dx_dy[1] = int(dx_dy[1])
        except ValueError:
            print "dx_dy must be a vector with two integers"
            return
            
        # check if they are bigger or equal to zero
        if (dx_dy[0]<0) or (dx_dy[1]<0):
            print "dx_dy must be a vector with two positive values"
    
    # allocate memory space for the result array
    result = numpy.zeros_like(x)
    
    if (dx_dy == None) or (dx_dy == [0, 0]):
        result = _math_tools.polyeval_2d(x, y, polyXY, result)
        return result
    else:
         # compute the coefficients of the derivative polynomial in x
        polyXY[0,:] = scipy.polyder(polyXY[0],m=dx_dy[0])
        # compute the coefficients of the derivative polynomial in y
        polyXY[1,:] = scipy.polyder(polyXY[1],m=dx_dy[1])
        # evaluate the 2d polynomial in x and y
        result = _math_tools.polyeval_2d(x, y, polyXY, result)
        return resultcipy.polyder(polyXY[1],m=dx_dy[1])
        # evaluate the 2d polynomial in x and y
        result = _math_tools.polyeval_2d(x, y, polyXY, result)
        return result
    
def poly2deval(coeffs_or_poly2d, x, y, dx_dy=None):    
    x = numpy.array(x, dtype="float64")
    y = numpy.array(y, dtype="float64")
    
    # check if x and y have the same dimensions
    if x.shape != y.shape:
        print "x and y must have the same shape"
        return
    
    # test if coeffs_or_poly2d is of type poly2d
    if isinstance(coeffs_or_poly2d, poly2d):
        mypoly = coeffs_or_poly2d
    else:
        # if not, transform it in a 2d polynomial
        mypoly = poly2d(coeffs=coeffs_or_poly2d)
        
    # allocate memory space for result of evaluation at x and y points
    if x.shape != ():
        if x.shape[0] == x.size:
            result_middle = numpy.zeros([mypoly.order[0]+1, x.shape[0]], dtype="float64")
            result = numpy.zeros_like(x)
        
        else:
            result_middle = numpy.zeros([mypoly.order[0]+1, x.shape[0], x.shape[1]], dtype="float64")
            result = numpy.zeros_like(x)
    else:
        result_middle = numpy.zeros(mypoly.order[0]+1, dtype="float64")
        result = 0.
    
    # check if 2d polynomial is to be computed at the points but
    # differentiating first is requested
    if (dx_dy == None) or (dx_dy == [0, 0]):
        # if no differentiation is requested
        coeffs = mypoly.coeffs
        order = mypoly.order
    else:
        # if differentiation is requested, first differentiate, then
        # store the coefficients and order, in order to evaluate the
        # differentiated polynomial at the requested points
        mypoly2dDeriv = mypoly.deriv(dx_dy)
        coeffs = mypoly2dDeriv.coeffs
        order = mypoly2dDeriv.order
        
    # compute int x_order * y_order multiplications instead of
    # x_order! * y_order! multiplications
    
    # for that computes first the polynomial at y values
    for i in range(0, order[0]+1):
        for j in range(0, order[1]+1):
            result_middle[i] = result_middle[i] * y + coeffs[i,j]
    
    # for that computes the polynomial at x values
    for i in range(0, order[0]+1):
        result = result * x + result_middle[i]
    
    return result

def poly2dderiv(coeffs_or_poly2d, dx_dy):
    
    # check the size of dx_dy, it must be vector of length 2
    try:
        if len(dx_dy) != 2:
            print "dx_dy should be a vector of size 2: [dx dy]"
            return
    except TypeError:
        print "dx_dy must be a vector of size 2: [dx dy]"
        return
    
    # try to convert both values to int
    try:
        dx_dy[0] = int(dx_dy[0])
        dx_dy[1] = int(dx_dy[1])
    except ValueError:
        print "dx_dy must be a vector with two integers"
        return
        
    # check if they are bigger or equal to zero
    if (dx_dy[0]<0) or (dx_dy[1]<0):
        print "dx_dy must be a vector with two positive values"
    
    # test if coeffs_or_poly2d is of type poly2d
    if isinstance(coeffs_or_poly2d, poly2d):
        mypoly = coeffs_or_poly2d
    else:
        mypoly = poly2d(coeffs=coeffs_or_poly2d)
        
    # check if the order of the derivative is higher than the order of the
    # 2d polynomial
    if (dx_dy[0]>mypoly.order[0]) or (dx_dy[1]>mypoly.order[1]):
        return poly2d(coeffs=numpy.array([[0]]))
        
    # compute the derivatives
    
    coeffs = mypoly.coeffs.copy()
    
    # y derivative
    for derivY in range(0, dx_dy[1]):
        coeffs = (coeffs * numpy.arange(coeffs.shape[1]-1,-1,-1))
        coeffs = coeffs[:,:-1]
        
    # x derivative
    for derivX in range(0, dx_dy[0]):
        coeffs = ((coeffs.T) * numpy.arange(coeffs.shape[0]-1,-1,-1)).T
        coeffs = coeffs[:-1,:]
        
    return poly2d(coeffs=coeffs)    
    
    
    