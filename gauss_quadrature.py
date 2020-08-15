## \package nereides.gauss_quadrature
# \brief Gauss Lobatto related classes and methods
#

# import NereidesException for handling with custom exceptions
import NereidesExceptions

# import math_tools for working with 2d polynomials
import math_tools
import _math_tools

# import definition variable
import settings

# from utilities import memoize decorator
from utilities import memoize

# import numpy for using optimized array objects
import numpy
import scipy
from scipy import special
import pylab

        
# class that holds the coefficients for the Lagrange polynomial of order p
# it is based on numpy.poly1d hence it has methods such as deriv, which
# returns the coefficients of the polynomial that represents the derivative
# of our polynomial
class lagrange_p:
    def __init__(self, p, q_type="lobatto"):
        # check if the quadrature type is one of the available ones
        available_q_types = settings.Q_TYPES
        if q_type not in available_q_types:
            raise NereidesExceptions.GaussQuadratureError("Polynomial type not\
                                     available: " + str(q_type) + " should be one of ",\
                                     str(available_q_types), "\n")
            return
        
        # store the order of the Lagrange polynomial
        self.order = p
        
        if q_type == "lobatto":
            # set the type of Lagrange polynomial
            self.q_type = "lobatto"
            # generate the Lobatto collocated Lagrange polynomial
            self.__compute_coefficients()
            
        if q_type == "gauss":
            # set the type of Lagrange polynomial
            self.q_type = "gauss"
            # generate the Gauss collocated Lagrange polynomial
            self.__compute_coefficients(q_type="gauss")
        
    def __compute_coefficients(self, q_type="lobatto"):
        # make local copy of polynomial order
        p = self.order
        
        # compute the lagrange polynomial coefficients
        
        # first compute the collocation points
        myCollocation_points = collocation_points(self.order, q_type=q_type)
        
        # first allocate memory space for the matrix
        self.coefficients = numpy.zeros((p+1,p+1), dtype="float64")
        # compute the coefficients for each lagrange polynomial and store them
        # as a row of lagrange_p

        # collocation points index list
        collocation_points_index_list = numpy.arange(p+1)

        for row,lagrange_point in enumerate(myCollocation_points):
            # create a vector that selects all the collocation points except the one
            # for which the lagrange polynomial evaluates to 1, that is the
            # lagrange_point
            lagrange_zeros = myCollocation_points[collocation_points_index_list != row]
            self.coefficients[row,:] = scipy.poly(lagrange_zeros) / \
                                       (lagrange_point - lagrange_zeros).prod()
    
    def deriv(self,m=1):
        # return the mth derivative of the lagrange_polynomials
        p = self.order
        # first allocate memory space for the matrix
        dlagrange_p = numpy.zeros((p+1,p+1-m), dtype="float64")

        # compute the coefficients of the derivative polynomial
        for row,lagrange in enumerate(self.coefficients):
            dlagrange_p[row,:] = scipy.polyder(lagrange, m=m)
        
        # return the derivative
        return dlagrange_p
    
    def evaluate_at(self, x):
        p = self.order
        # allocate memory space for the matrix with the values at the collocation
        # points, for all derivatives of the Lagrange polynomials
        lagrange_p_at_x = numpy.zeros([p+1,x.size], dtype="float64")

        # return an array whose rows are the value of the corresponding lagrange
        # polynomial computed at points given by x
        for row, lagrange in enumerate(self.coefficients):
            lagrange_p_at_x[row] = scipy.polyval(lagrange, x)
        
        # return the values
        return lagrange_p_at_x
    
    def evaluate_deriv_at(self, x, m=1):
        p = self.order
        # allocate memory space for the matrix with the values at the collocation
        # points, for all derivatives of the Lagrange polynomials
        dlagrange_p_at_x = numpy.zeros([p+1,x.size], dtype="float64")

        # first compute the coefficients of the polynomial that represents the
        # derivative of the lagrange polynomial
        dlagrange_p = self.deriv(m=m)
        
        # return an array whose rows are the value of the derivative of the
        #  corresponding lagrange polynomial computed at points given by x        
        for row, dlagrange in enumerate(dlagrange_p):
            dlagrange_p_at_x[row] = scipy.polyval(dlagrange, x)
        
        # return the values
        return dlagrange_p_at_x
    
    def inner_product(self, derivative_n, derivative_m, p_int):
        # compute the array of inner products
        #     <d_n(phi_n), d_m(phi_m)>
        # where d_n represents the nth order derivative and d_m represents
        # the mth order derivative.
        # returns an array with the all possible combinations of inner products
        # between the lagrange polynomials, with derivatives computed as
        # given at the input parameters. p_int is the integration order to be
        # used.
        
        # compute the integration collocation points
        int_collocation_points = collocation_points(p_int, q_type=self.q_type)
        
        # compute the integration weights
        gl_weights = gauss_weights(p_int, q_type=self.q_type)
        
        # compute the values of of the lagrange polynomials at the collocation
        # points
        n_phis = self.evaluate_deriv_at(int_collocation_points, m=derivative_n)
        m_phis = self.evaluate_deriv_at(int_collocation_points, m=derivative_m)
        
        # allocate memory space for the array of inner products
        phis_integrated = numpy.zeros([self.order+1,self.order+1],\
                                      dtype="float64")
                                           
        # no scalling is made here
        for n,n_phi in enumerate(n_phis):
            for m,m_phi in enumerate(m_phis):
                phis_integrated[n,m] = (n_phi*m_phi*gl_weights).sum()
                
        return phis_integrated
    
    
# class that holds the coefficients for the 2d Lagrange polynomial of order p
# it is based on lagrange_p hence it has methods such as deriv, which
# returns the coefficients of the polynomial that represents the derivative
# of our polynomial, but adapted for the 2d case.
class lagrange_p_2d:
    def __init__(self, p, q_type=("lobatto","lobatto")):
        # store the order of the lagrange polynomial, it is the same on both
        # directions, it is like this, since there will be no difference
        # between both directions, since, in general, there is no way of
        # doing this, because the elements, in general, are not always lined
        # up with the axis.
        
        self.__dict__["order"] = p
        
        try:
            p[0]
        except TypeError:
            self.__dict__["order"] = [p,p]
            
        self.__dict__["q_type"] = q_type
        
        # generate the lagrange polynomials
        self.__compute_coefficients()

    def __compute_coefficients(self):
        # make local copy of polynomial order
        p = self.order
        q_type = self.q_type
        
        # compute the 1d lagrange polynomial coefficients, which make up the
        # 2d, since varphi_ij(u,v) = phi_i(u) * psi_j(v)
        lagrange_1d_u = lagrange_p(p[0],q_type=q_type[0])
        
        if (p[0] == p[1]) and (q_type[0] == q_type[1]):
            lagrange_1d_v = lagrange_1d_u
        else:
            lagrange_1d_v = lagrange_p(p[1], q_type[1])
        
        # store the coefficients locally
        coeffs_u = lagrange_1d_u.coefficients
        
        if (p[0] == p[1]) and (q_type[0] == q_type[1]):
            coeffs_v = coeffs_u
        else:
            coeffs_v = lagrange_1d_v.coefficients
        
        # since the first derivative is much likely to be used very often
        # compute it right away to reduce the overhead of computing it everytime
        
        # allocate memory space
        d_coeffs_u = numpy.zeros([p[0]+1, p[0]], dtype="float64")
        if (p[0] != p[1]) or (q_type[0] != q_type[1]):
            d_coeffs_v = numpy.zeros([p[1]+1, p[1]], dtype="float64")
        
        # loop over polynomials of the basis function and compute derivatives
        for i in range(0, p[0]+1):
            d_coeffs_u[i,:] = scipy.polyder(coeffs_u[i], m=1)
            
        if (p[0] != p[1]) or (q_type[0] != q_type[1]): 
            for i in range(0, p[1]+1):
                d_coeffs_v[i,:] = scipy.polyder(coeffs_v[i], m=1)
        else:
            d_coeffs_v = d_coeffs_u
        
        # the 2d lagrange polynomials are stored in an array
        
        # initialize self.basis as a tuple, for later filling
        self.__dict__["basis_functions"] = []
        
        # check the 0 padding in the coefficients in order to have the coefficient
        # vectors all with the same size
        if p[0] > p[1]:
            zero_padding = [0, p[0]-p[1]]
        elif p[0] < p[1]:
            zero_padding = [p[1]-p[0], 0]
        else:
            zero_padding = [0, 0]
        
        # loop over all combinations of basis functions and create them
        for phi_i in range(0, p[0]+1):
            # append one more row to basis functions
            self.basis_functions.append([])
            for psi_j in range(0, p[1]+1):
                # append one more basis function
                self.basis_functions[phi_i].append(\
                     math_tools.basis_function_2d(coeffs=numpy.array([numpy.concatenate((numpy.zeros(zero_padding[0]), coeffs_u[phi_i])), numpy.concatenate((numpy.zeros(zero_padding[1]), coeffs_v[psi_j]))]),\
                                                  d_coeffs=numpy.array([numpy.concatenate((numpy.zeros(zero_padding[0]), d_coeffs_u[phi_i])), numpy.concatenate((numpy.zeros(zero_padding[1]), d_coeffs_v[psi_j]))])))
        
        # store coeffs and d_coeffs    
        self.__dict__["_coeffs_u"] = numpy.asfortranarray(coeffs_u)
        if (p[0] == p[1]) and (q_type[0] == q_type[1]):
            self.__dict__["_coeffs_v"] = self.__dict__["_coeffs_u"]
        else:
            self.__dict__["_coeffs_v"] = numpy.asfortranarray(coeffs_v)
        self.__dict__["_d_coeffs_u"] = numpy.asfortranarray(d_coeffs_u)
        if (p[0] == p[1]) and (q_type[0] == q_type[1]):
            self.__dict__["_d_coeffs_v"] = self.__dict__["_d_coeffs_u"]
        else:
            self.__dict__["_d_coeffs_v"] = numpy.asfortranarray(d_coeffs_v)
        
    def __len__(self):
        return self.order    
        
    def __call__(self, u, v, du_dv=None):
        p = self.order
        
        u = numpy.array(u, dtype="float64", order="fortran")
        v = numpy.array(v, dtype="float64", order="fortran")
        
        # allocate memory space for result of evaluation at u and v points
        if u.shape != ():
            if u.shape[0] == u.size:
                result = numpy.zeros([p[0]+1,p[1]+1,u.shape[0]], dtype="float64", order="fortran")
            
            else:
                result = numpy.zeros([p[0]+1,p[1]+1,u.shape[0],u.shape[1]], dtype="float64", order="fortran")
        else:
            result = numpy.zeros([p[0]+1,p[1]+1], dtype="float64", order="fortran")
            
        # loop over all basis functions and evaluate them at x and y
        # this code was the older one, it is not much slower, an order of 2 times
        # but it is always good to speed things up hence the newer code below
##        for phi_i in range(0, p+1):
##            for psi_j in range(0, p+1):
##                result[phi_i, psi_j] = self[phi_i,psi_j](u, v, du_dv)
        
        # this is the faster code, less elegant but faster
        if (du_dv == None) or (du_dv == [0, 0]):
            result = _math_tools.polyeval_all_2d(u, v, self._coeffs_u, self._coeffs_v, result)
            return result
        elif du_dv == [1, 0]:
            result = _math_tools.polyeval_all_2d(u, v, self._d_coeffs_u, self._coeffs_v, result)
            return result
        elif du_dv == [0, 1]:
            result = _math_tools.polyeval_all_2d(u, v, self._coeffs_u, self._d_coeffs_v, result)
            return result
        else:
            print "unavailable derivative"
    
    def __getitem__(self, i):
        return self.basis_functions[i[0]][i[1]]
    
    def all_inner_products(self, transfinite, p_int):
        # Computes the inner product between the two functions: f_0 and f_1, that is,
        # computes the following integral:
        #
        #    \int\int\_{\Omega_{e}}frac{\partial^{n^{x}_{0}}}{\partial x^{n^{x}_{0}}}\frac{\partial^{n^{y}_{0}}}{\partial y^{n^{y}_{0}}}(f_{0}(x,y))\cdot \frac{\partial^{n^{x}_{1}}}{\partial x^{n^{x}_{1}}}\frac{\partial^{n^{y}_{1}}}{\partial y^{n^{y}_{1}}}(f_{1}(x,y))\mathrm{d}x\mathrm{d}y
        #
        # Where n^{x}_{i} and n^{y}_{i} are defined in dx_dy, and state the order
        # of the derivative present in f_i in the integral.
        #
        # Notice that the integral is performed in \Omega_{e}, the element, in
        # global coordinates.
        #
        # The integral is not an analytical one, instead it is a performed a numerical
        # integration using a Gauss-Lobatto quadrature or Radau quadrature in 2d.
        # see for the 1d case:
        #     http://mathworld.wolfram.com/LobattoQuadrature.html
        #     Atkinson, An Introduction to Numerical Analysis
        #
        # Hence, a transformation of coordinates is made,
        # using transfinite, in order to pass the integral to \Omega', the (u,v)
        # general coordinates. So that the Gauss Quadrature is easily computed
        #
        # transfinite :: is a valid TransfiniteMapping_2d instance
        #
        # p_int :: since Gauss-Lobatto quadrature is used for the integration then
        #          p_int specifies the order of the quadrature used.
        
        p = self.order
        q_type = self.q_type
        
        # compute the 1d collocation coordinates based on the integration order
        #collocation_coordinates = collocation_points(p_int)
        collocation_coordinates_u = collocation_points(p[0],q_type=q_type[0])
        collocation_coordinates_v = collocation_points(p[1],q_type=q_type[1])
        
        # allocate memory space for the collocation coordinates
        # and make sure they are of type float64 and have fortran ordering
        u = numpy.zeros([p[0]+1, p[1]+1], dtype="float64", order="fortran")
        v = numpy.zeros([p[0]+1, p[1]+1], dtype="float64", order="fortran")
        
        # compute the u,v coordinates of the collocation points based on the 1d
        # collocation_coordinates
        v[:,:], u[:,:] = numpy.meshgrid(collocation_coordinates_v, collocation_coordinates_u)
        
        # compute the basis functions:
        # varphi_ij
        # dx_varphi_ij
        # dy_varphi_ij
        # which will be needed to compute all the inner product of the basis
        # functions
        varphi_ij = self(u, v, du_dv=[0,0])
        dv_varphi_ij = self(u, v, du_dv=[0,1])
        du_varphi_ij = self(u, v, du_dv=[1,0])
        
        # compute the transfinite mappings
        dv_T_Y = numpy.asfortranarray(transfinite.dY_dv(u, v))
        if dv_T_Y.shape==(1,):
            dv_T_Y = numpy.asfortranarray(numpy.ones_like(u)*dv_T_Y)
        du_T_Y = numpy.asfortranarray(transfinite.dY_du(u, v))
        if du_T_Y.shape==(1,):
            du_T_Y = numpy.asfortranarray(numpy.ones_like(u)*du_T_Y)
        dv_T_X = numpy.asfortranarray(transfinite.dX_dv(u, v))
        if dv_T_X.shape==(1,):
            dv_T_X = numpy.asfortranarray(numpy.ones_like(u)*dv_T_X)
        du_T_X = numpy.asfortranarray(transfinite.dX_du(u, v))
        if du_T_X.shape==(1,):
            du_T_X = numpy.asfortranarray(numpy.ones_like(u)*du_T_X)
        
        # compute J for all possible combinations of derivatives:
        # 0 derivatives
        # 1 derivative
        # 2 derivatives
        J_0 = numpy.asfortranarray(numpy.abs(transfinite.J(u,v)))
        if J_0.shape==(1,):
            J_0 = numpy.asfortranarray(numpy.ones_like(u)*J_0)
        J_1 = numpy.asfortranarray(numpy.sign(transfinite.J(u,v)))
        if J_1.shape==(1,):
            J_1 = numpy.asfortranarray(numpy.ones_like(u)*J_1)
        J_2 = numpy.asfortranarray(1.0/numpy.abs(transfinite.J(u,v)))
        if J_2.shape==(1,):
            J_2 = numpy.asfortranarray(numpy.ones_like(u)*J_2)
        
        # compute the weights for the Gauss-Lobatto quadrature
        weights_u = gauss_weights(p[0], q_type=q_type[0])
        weights_v = gauss_weights(p[1], q_type=q_type[1])
        
        # compute the inner product
        # the 6 comes from the fact that we wish to compute the following six
        # inner products: <phi_ij, phi_nm>, <phi_ij, du_phi_nm>, <phi_ij, dv_phi_nm>,
        # <du_phi_ij, du_phi_nm>, <dv_phi_ij, dv_phi_nm> and <du_phi_ij, dv_phi_nm>
        result = numpy.zeros([6, p[0]+1, p[1]+1, p[0]+1, p[1]+1], dtype="float64", order="fortran")
        result = _math_tools.all_basis_inner_products(varphi_ij, du_varphi_ij,\
                                                      dv_varphi_ij, du_T_X, dv_T_X, \
                                                      du_T_Y, dv_T_Y, J_0, J_1, J_2,\
                                                      weights_u, weights_v, result)
        return result
    
    def all_other_inner_products(self, otherLagrange_p_2d, du_dv, transfinite, p_int):
        # Computes the inner product between the all the basis functions of the
        # lagrange_p_2d object witht all the basis functions of the lagrange_p_2d
        # object given as input (otherLagrange_p_2d). That is computes the following integral:
        #
        #    \int\int\_{\Omega_{e}}frac{\partial^{n^{x}_{0}}}{\partial x^{n^{x}_{0}}}\frac{\partial^{n^{y}_{0}}}{\partial y^{n^{y}_{0}}}(f_{0}(x,y))\cdot \frac{\partial^{n^{x}_{1}}}{\partial x^{n^{x}_{1}}}\frac{\partial^{n^{y}_{1}}}{\partial y^{n^{y}_{1}}}(f_{1}(x,y))\mathrm{d}x\mathrm{d}y
        #
        # Where n^{x}_{i} and n^{y}_{i} are defined in dx_dy, and state the order
        # of the derivative present in f_i in the integral.
        #
        # Notice that the integral is performed in \Omega_{e}, the element, in
        # global coordinates.
        #
        # The integral is not an analytical one, instead it is a performed a numerical
        # integration using a Gauss-Lobatto quadrature or Radau quadrature in 2d.
        # see for the 1d case:
        #     http://mathworld.wolfram.com/LobattoQuadrature.html
        #     Atkinson, An Introduction to Numerical Analysis
        #
        # Hence, a transformation of coordinates is made,
        # using transfinite, in order to pass the integral to \Omega', the (u,v)
        # general coordinates. So that the Gauss Quadrature is easily computed
        #
        # otherLagrange_p_2d :: is a valid lagrange_p_2d object
        #
        # du_dv :: is a 2x2 matrix specifying which inner product to make. Examples:
        #          du_dv = [[0,0],[0,0]] the inner products are:
        #                  <phi_ij, psi_nm>
        #
        #          du_dv = [[1,0],[0,0]] the inner products are:
        #                  <du_phi_ij, psi_nm>
        #
        #          du_dv = [[0,1],[1,0]] the inner products are:
        #                  <dv_phi_ij, du_psi_nm>
        #
        # transfinite :: is a valid TransfiniteMapping_2d instance
        #
        # p_int :: since Gauss-Lobatto quadrature is used for the integration then
        #          p_int specifies the order of the quadrature used.
        
        # the order of the self Lagrange polynomial
        myp = self.order
        myq_type = self.q_type
        
        du_dv = numpy.array(du_dv)
        
        # the order of the other Lagrange polynomial
        otherp = otherLagrange_p_2d.order
        
        # compute the 1d collocation coordinates based on the integration order
        collocation_coordinates_u = collocation_points(myp[0],q_type=myq_type[0])
        collocation_coordinates_v = collocation_points(myp[1],q_type=myq_type[1])
        
        # allocate memory space for the collocation coordinates
        # and make sure they are of type float64 and have fortran ordering
        u = numpy.zeros([myp[0]+1, myp[1]+1], dtype="float64", order="fortran")
        v = numpy.zeros([myp[0]+1, myp[1]+1], dtype="float64", order="fortran")
        
        # compute the u,v coordinates of the collocation points based on the 1d
        # collocation_coordinates
        v[:,:], u[:,:] = numpy.meshgrid(collocation_coordinates_v, collocation_coordinates_u)
        
        # compute the basis functions:
        # varphi_ij
        # dx_varphi_ij
        # dy_varphi_ij
        # which will be needed to compute all the inner product of the basis
        # functions
        # the choice of which to compute is based upon the du_dv input
        varphi_ij = self(u, v, du_dv=[0,0])
        du_varphi_ij = self(u, v, du_dv=[1,0])
        dv_varphi_ij = self(u, v, du_dv=[0,1])
            
        # compute the basis functions:
        # other_varphi_ij
        # other_dx_varphi_ij
        # other_dy_varphi_ij
        # other_which will be needed to compute all the inner product of the basis
        # functions
        # the choice of which to compute is based upon the du_dv input
        other_varphi_ij = otherLagrange_p_2d(u, v, du_dv=[0,0])
        du_other_varphi_ij = otherLagrange_p_2d(u, v, du_dv=[1,0])
        dv_other_varphi_ij = otherLagrange_p_2d(u, v, du_dv=[0,1])

        
        # compute the transfinite mappings
        dv_T_Y = numpy.asfortranarray(transfinite.dY_dv(u, v))
        if dv_T_Y.shape==(1,):
            dv_T_Y = numpy.asfortranarray(numpy.ones_like(u)*dv_T_Y)
        du_T_Y = numpy.asfortranarray(transfinite.dY_du(u, v))
        if du_T_Y.shape==(1,):
            du_T_Y = numpy.asfortranarray(numpy.ones_like(u)*du_T_Y)
        dv_T_X = numpy.asfortranarray(transfinite.dX_dv(u, v))
        if dv_T_X.shape==(1,):
            dv_T_X = numpy.asfortranarray(numpy.ones_like(u)*dv_T_X)
        du_T_X = numpy.asfortranarray(transfinite.dX_du(u, v))
        if du_T_X.shape==(1,):
            du_T_X = numpy.asfortranarray(numpy.ones_like(u)*du_T_X)
        
        # compute J for all possible combinations of derivatives:
        # 0 derivatives
        # 1 derivative
        # 2 derivatives
        if du_dv.sum() == 0:
            Jacobian = numpy.asfortranarray(numpy.abs(transfinite.J(u,v)))
            if Jacobian.shape==(1,):
                Jacobian = numpy.asfortranarray(numpy.ones_like(u)*Jacobian)
        elif du_dv.sum() == 1:
            Jacobian = numpy.asfortranarray(numpy.sign(transfinite.J(u,v)))
            if Jacobian.shape==(1,):
                Jacobian = numpy.asfortranarray(numpy.ones_like(u)*Jacobian)
        elif du_dv.sum() == 2:
            Jacobian = numpy.asfortranarray(1.0/numpy.abs(transfinite.J(u,v)))
            if Jacobian.shape==(1,):
                Jacobian = numpy.asfortranarray(numpy.ones_like(u)*Jacobian)
        else:
            raise(ValueError("Too many derivatives"))

        # compute the type of du_dv
        # example of results: 
        #       du_dv = [[0,0],[0,0]] --> 0
        #       du_dv = [[1,0],[0,0]] --> 1
        #       du_dv = [[0,1],[0,0]] --> 2
        #       du_dv = [[0,1],[0,1]] --> 10
        type_du_dv = int(du_dv[0,0] + du_dv[0,1]*2 + du_dv[1,0]*4 + du_dv[1,1]*8)
        
        # compute the weights for the Gauss-Lobatto or Gauss-Gauss quadrature
        weights_u = gauss_weights(myp[0], q_type=myq_type[0])
        weights_v = gauss_weights(myp[1], q_type=myq_type[1])
        
        # compute the inner product
        # the 6 comes from the fact that we wish to compute the following six
        # inner products: <phi_ij, phi_nm>, <phi_ij, du_phi_nm>, <phi_ij, dv_phi_nm>,
        # <du_phi_ij, du_phi_nm>, <dv_phi_ij, dv_phi_nm> and <du_phi_ij, dv_phi_nm>
        result = numpy.zeros([myp[0]+1, myp[1]+1, otherp[0]+1, otherp[1]+1], dtype="float64", order="fortran")
        result = _math_tools.all_other_basis_inner_products(varphi_ij, du_varphi_ij,\
                                                      dv_varphi_ij, \
                                                      other_varphi_ij, du_other_varphi_ij,\
                                                      dv_other_varphi_ij, \
                                                      du_T_X, dv_T_X, du_T_Y, dv_T_Y, \
                                                      Jacobian,\
                                                      weights_u, weights_v, type_du_dv, result)
        return result    
    
    def all_function_inner_products(self, f, transfinite, p_int, coeffs, derivs):
        # method that computes all the combinations of inner products between
        # the function f and the basis functions, returning a matrix with
        # dimensions (p+1)x(p+1) with the inner products of f with the
        # basis function varphi_ij
        p = self.order
        
        coeffs = numpy.array(coeffs, dtype="float64", order="fortran")
        derivs = numpy.array(derivs, dtype="int16", order="fortran")
        
        # compute the 1d collocation coordinates based on the integration order
        collocation_coordinates = collocation_points(p_int)
        
        # allocate memory space for the collocation coordinates
        # and make sure they are of type float64 and have fortran ordering
        u = numpy.zeros([p_int+1, p_int+1], dtype="float64", order="fortran")
        v = numpy.zeros([p_int+1, p_int+1], dtype="float64", order="fortran")
        
        # compute the u,v coordinates of the collocation points based on the 1d
        # collocation_coordinates
        v[:,:], u[:,:] = numpy.meshgrid(collocation_coordinates, collocation_coordinates)
        
        # compute the function f at the collocation points
        # allocate memory space
        evaluated_f = numpy.zeros_like(u)
        evaluated_f[:,:] = f(transfinite.X(u,v), transfinite.Y(u,v))
        
        # compute the basis functions:
        # varphi_ij
        # dx_varphi_ij
        # dy_varphi_ij
        # which will be needed to compute all the inner product of the basis
        # functions
        varphi_ij = self(u, v, du_dv=[0,0])
        du_varphi_ij = self(u, v, du_dv=[1,0])
        dv_varphi_ij = self(u, v, du_dv=[0,1])
        
        # compute the transfinite mappings
        dv_T_Y = numpy.asfortranarray(transfinite.dY_dv(u, v))
        if dv_T_Y.shape==(1,):
            dv_T_Y = numpy.asfortranarray(numpy.ones_like(u)*dv_T_Y)
        du_T_Y = numpy.asfortranarray(transfinite.dY_du(u, v))
        if du_T_Y.shape==(1,):
            du_T_Y = numpy.asfortranarray(numpy.ones_like(u)*du_T_Y)
        dv_T_X = numpy.asfortranarray(transfinite.dX_dv(u, v))
        if dv_T_X.shape==(1,):
            dv_T_X = numpy.asfortranarray(numpy.ones_like(u)*dv_T_X)
        du_T_X = numpy.asfortranarray(transfinite.dX_du(u, v))
        if du_T_X.shape==(1,):
            du_T_X = numpy.asfortranarray(numpy.ones_like(u)*du_T_X)
        
        # compute J for all possible combinations of derivatives:
        # 0 derivatives
        # 1 derivative
        J_0 = numpy.asfortranarray(numpy.abs(transfinite.J(u,v)))
        if J_0.shape==(1,):
            J_0 = numpy.asfortranarray(numpy.ones_like(u)*J_0)
        J_1 = numpy.asfortranarray(numpy.sign(transfinite.J(u,v)))
        if J_1.shape==(1,):
            J_1 = numpy.asfortranarray(numpy.ones_like(u)*J_1)
        
        # compute the weights for the Gauss-Lobatto quadrature
        weights = gauss_weights(p_int)
        # compute the inner product
        result = numpy.zeros([p[0]+1, p[1]+1], dtype="float64", order="fortran")
        result = _math_tools.all_function_basis_inner_products(evaluated_f, varphi_ij, du_varphi_ij,\
                                                               dv_varphi_ij, du_T_X, dv_T_X, \
                                                               du_T_Y, dv_T_Y, J_0, J_1, \
                                                               weights, coeffs, derivs, result)
        return result
    
    def all_function_inner_products_v2(self, f, transfinite, p_int, coeffs, du_dv):
        # method that computes all the combinations of inner products between
        # the function f and the basis functions, returning a matrix with
        # dimensions (p+1)x(p+1) with the inner products of f with the
        # basis function varphi_ij
        p = self.order
        q_type= self.q_type
        
        coeffs = numpy.array(coeffs, dtype="float64", order="fortran")
        du_dv = numpy.array(du_dv, dtype="int16", order="fortran")
        
        # compute the 1d collocation coordinates based on the integration order
        collocation_coordinates_x = collocation_points(p_int, q_type=q_type[0])
        collocation_coordinates_y = collocation_points(p_int, q_type=q_type[1])
        
        # allocate memory space for the collocation coordinates
        # and make sure they are of type float64 and have fortran ordering
        u = numpy.zeros([p_int+1, p_int+1], dtype="float64", order="fortran")
        v = numpy.zeros([p_int+1, p_int+1], dtype="float64", order="fortran")
        
        # compute the u,v coordinates of the collocation points based on the 1d
        # collocation_coordinates
        v[:,:], u[:,:] = numpy.meshgrid(collocation_coordinates_x, collocation_coordinates_y)
        
        # compute the function f at the collocation points
        # allocate memory space
        evaluated_f = numpy.zeros_like(u)
        evaluated_f[:,:] = f(transfinite.X(u,v), transfinite.Y(u,v))
        
        # compute the basis functions:
        # varphi_ij
        # dx_varphi_ij
        # dy_varphi_ij
        # which will be needed to compute all the inner product of the basis
        # functions
        varphi_ij = self(u, v, du_dv=[0,0])
        du_varphi_ij = self(u, v, du_dv=[1,0])
        dv_varphi_ij = self(u, v, du_dv=[0,1])
        
        # compute the transfinite mappings
        dv_T_Y = numpy.asfortranarray(transfinite.dY_dv(u, v))
        if dv_T_Y.shape==(1,):
            dv_T_Y = numpy.asfortranarray(numpy.ones_like(u)*dv_T_Y)
        du_T_Y = numpy.asfortranarray(transfinite.dY_du(u, v))
        if du_T_Y.shape==(1,):
            du_T_Y = numpy.asfortranarray(numpy.ones_like(u)*du_T_Y)
        dv_T_X = numpy.asfortranarray(transfinite.dX_dv(u, v))
        if dv_T_X.shape==(1,):
            dv_T_X = numpy.asfortranarray(numpy.ones_like(u)*dv_T_X)
        du_T_X = numpy.asfortranarray(transfinite.dX_du(u, v))
        if du_T_X.shape==(1,):
            du_T_X = numpy.asfortranarray(numpy.ones_like(u)*du_T_X)
        
        # compute J for all possible combinations of derivatives:
        # 0 derivatives
        # 1 derivative
        J_0 = numpy.asfortranarray(numpy.abs(transfinite.J(u,v)))
        if J_0.shape==(1,):
            J_0 = numpy.asfortranarray(numpy.ones_like(u)*J_0)
        J_1 = numpy.asfortranarray(numpy.sign(transfinite.J(u,v)))
        if J_1.shape==(1,):
            J_1 = numpy.asfortranarray(numpy.ones_like(u)*J_1)
        
        # compute the weights for the Gauss-Lobatto quadrature
        weights_x = gauss_weights(p_int, q_type=q_type[0])
        weights_y = gauss_weights(p_int, q_type=q_type[1])
        
        # compute the inner product
        result = numpy.zeros([p[0]+1, p[1]+1], dtype="float64", order="fortran")
        result = _math_tools.all_function_basis_inner_products_v2(evaluated_f, varphi_ij, du_varphi_ij,\
                                                               dv_varphi_ij, du_T_X, dv_T_X, \
                                                               du_T_Y, dv_T_Y, J_0, J_1, \
                                                               weights_x, weights_y, coeffs, du_dv, result)
        return result
    
@memoize    
def collocation_points(p, q_type="lobatto"):
    # computes the collocation points for different types of quadratures. The
    # implemented quadrature types (q_type) are:
    #     - lobatto
    #     - gauss
    
    # LOBATTO
    # The collocation points are p+1 points, where the endpoints are -1 and 1 
    # and the inner points are the roots of the derivative of the Legendre 
    # polynomial of order p, dlegendre_p. The roots of dlegendre_p
    # are the same as the roots of the Jacobi polynomial of order p-1 with 
    # alpha=1 and beta=1, jacobi_{p-1}^{11}. Hence one computes these roots
    # using the built in function.
    
    # GAUSS
    # The collocation points are p+1 points, contrasting with LOBATTO quadrature
    # there are no fixed end points. The inner points are the roots
    # of the Legendre polynomial of order p_int, legendre_p_int. Since the roots
    # of the Legendre polynomial of order p_int are the same as the ones of
    # the Jacobi polynomial of order p_int with alpha=0 and beta=0 
    # (since they are identical), one 
    # computes these roots using the built in function.
    
    # check if the quadrature type is one of the available ones
    available_q_types = settings.Q_TYPES
    if q_type not in available_q_types:
        raise NereidesExceptions.GaussQuadratureError("Quadrature type not\
                                 available: " + str(q_type) + " should be one of ",\
                                 str(available_q_types), "\n")
        return
    
    # LOBATTO
    if q_type == "lobatto":
        # Root finding of derivative of Legendre polynomial of order p
        
        # allocate memory space of the roots
        collocation_points = numpy.zeros(p+1, dtype="float64")
        
        # [0] at the end is used because special.orthogonal.j_roots returns two
        # vectors, the first, the index 0, is the roots, the second, index 1, is
        # the weights for the Gauss-Jacobi quadrature, which we do not use
        if p !=1:
            collocation_points[1:-1] = (special.orthogonal.j_roots(p-1,1,1)[0]).real
        collocation_points[0] = -1.0
        collocation_points[-1] = 1.0

        # The collocation points just computed are the local collocation points,
        # relative to the element. One also needs, for future computations, the global
        # coordinates of the collocation points. 
        
    # GAUSS
    if q_type == "gauss":
        # Root finding of Legendre polynomial of order p+1
        
        # allocate memory space of the roots
        collocation_points = numpy.zeros(p+1, dtype="float64")
        
        # [0] at the end is used because special.orthogonal.j_roots returns two
        # vectors, the first, the index 0, is the roots, the second, index 1, is
        # the weights for the Gauss-Jacobi quadrature, which we do not use
        collocation_points = special.orthogonal.j_roots(p+1,0,0)[0]

        # The collocation points just computed are the local collocation points,
        # relative to the element. One also needs, for future computations, the global
        # coordinates of the collocation points.
    
    return collocation_points

@memoize
def gauss_weights(p_int, q_type="lobatto"):
    # Computes the Gauss quadrature weights in order to perform integration. The
    # implemented quadrature types (q_type) are:
    #     - lobatto
    #     - gauss
    
    # LOBATTO
    # The collocation points are p_int+1 points, where the endpoints are -1 and 1 
    # and the inner points are the roots of the derivative of the Legendre 
    # polynomial of order p_int, dlegendre_p_int. The roots of dlegendre_p_int 
    # are the same as the roots of the Jacobi polynomial of order p_int-1 with 
    # alpha=1 and beta=1, jacobi_{p_int-1}^{11}. Hence one computes these roots
    # using the built in function.
    
    # GAUSS
    # The collocation points are p_int+1 points, contrasting with LOBATTO quadrature
    # there are no fixed end points. The inner points are the roots
    # of the Legendre polynomial of order p_int, legendre_p_int. Since the roots
    # of the Legendre polynomial of order p_int are the same as the ones of
    # the Jacobi polynomial of order p_int with alpha=0 and beta=0 
    # (since they are identical), one 
    # computes these roots using the built in function.
    
    # check if the quadrature type is one of the available ones
    available_q_types = settings.Q_TYPES
    if q_type not in available_q_types:
        raise NereidesExceptions.GaussQuadratureError("Quadrature type not\
                                 available: " + str(q_type) + " should be one of ",\
                                 str(available_q_types), "\n")
        return
    
    # LOBATTO
    if q_type == "lobatto":
        # Root finding of derivative of Legendre polynomial of order p_int

        # [0] at the end is used because special.orthogonal.j_roots returns two
        # vectors, the first, the index 0, is the roots, the second, index 1, is
        # the weights for the Gauss-Jacobi quadrature, which we do not use
            
        if p_int!=1:
            dlegendre_p_int_roots = special.orthogonal.j_roots(p_int-1,1,1)[0]

            # create the collocation points vector, which is simply the endpoints -1 and 1
            # together with the dlegendre_p_int_roots
            int_collocation_points = numpy.concatenate(([-1], dlegendre_p_int_roots, [1]))
        else:
            int_collocation_points = numpy.array([-1,1])

        # The collocation points just computed are the local collocation points,
        # relative to the element. One also needs, for future computations, the global
        # coordinates of the collocation points.
        
        # compute the coefficients for the Legendre polynomial of order p_int
        legendre_p_int = special.legendre(p_int)

        # compute the weights for the Gauss-Lobatto quadrature using the expression
        #     w_i = 2 / ( p_int . (p_int + 1) . [P_p_int(x_i)]^2 )
        #     x_i = the roots of derivative of Legendre polynomial of order p_int-1
        gq_weights = 2.0 / ( p_int*(p_int+1) * scipy.polyval(legendre_p_int, \
                     int_collocation_points)**2 ).real
    
    # GAUSS
    if q_type == "gauss":
##        # Root finding of Legendre polynomial of order p_int
##
##        # [0] at the end is used because special.orthogonal.j_roots returns two
##        # vectors, the first, the index 0, is the roots, the second, index 1, is
##        # the weights for the Gauss-Jacobi quadrature, which we do not use
##        int_collocation_points = special.orthogonal.j_roots(p_int+1,0,0)[0]
##
##        # The collocation points just computed are the local collocation points,
##        # relative to the element. One also needs, for future computations, the global
##        # coordinates of the collocation points.
##        
##        # compute the coefficients for the Legendre polynomial of order p_int
##        legendre_p_int = special.legendre(p_int)
##
##        # compute the weights for the Gauss-Lobatto quadrature using the expression
##        #     w_i = 2 / ( p_int . (p_int + 1) . [P_p_int(x_i)]^2 )
##        #     x_i = the roots of derivative of Legendre polynomial of order p_int-1
##        gq_weights = 2.0 * (1.0 - int_collocation_points**2) /\
##                     ((scipy.polyval(legendre_p_int, int_collocation_points)**2) *\
##                      (p_int+1)**2)\
                     
        gq_weights = special.orthogonal.j_roots(p_int+1,0,0)[1]
                
    return gq_weights

###############################################################################
## DEPRECATED
###############################################################################
##
##def function_inner_product(function, lagrange_deriv, a, p, p_int, mesh):
##    # compute the inner product over each element of function with the lagrange
##    # polynomials.
##    # It computes the following vector (using the global numbering):
##    #
##    #      SUM_k a^k<function, d^k_1(phi_n)*d^k_2(phi_m)>
##    #
##    # lagrange_deriv is an array whose rows are indexed with index k of the 
##    # previous expression and whose columns contain the order of the derivative
##    # to use in the k term, for each lagrange polynomial. a is a vector with
##    # the coefficients of the a^k.
##    # Example:
##    #
##    # consider that lagrange_deriv is:
##    #     1  0
##    #     0  2
##    #     3  4
##    #
##    # and a is
##    #     5
##    #     6
##    #     7
##    #
##    # then this function computes
##    #
##    #     5 * <function, d_1(phi_n)*phi_m> + 6 * <function, phi_n*d_2(phi_m)> +
##    #     7 * <function, d_3(phi_n)*d_4(phi_m)>
##    #
##    # for each element and for each n and m
##    #
##    # function should be a function evaluatable in 2 dimensions, x and y.
##    
##    # compute lagrange polynomials
##    mylagrange = lagrange_p(p)
##    
##    # compute the gauss lobatto weights
##    gl_weights = gauss_weights(p_int)
##    
##    # compute integration collocation points
##    int_collocation_points = collocation_points(p_int)
##    
##    # with them, compute arrays of x and y coordinates
##    int_collocation_y, int_collocation_x = numpy.meshgrid(int_collocation_points,\
##                                                  int_collocation_points)
##    
##    # make a local copy to the number of main nodes
##    N_nodes = mesh.N_nodes
##    
##    # make a local copy to the number of elements
##    N_elements = mesh.N_elements
##    
##    # make a local copy to the number of edges
##    N_edges = mesh.N_edges
##    
##    # compute the total number of nodes, assuming that all elements have the
##    # same spectral order p (order of the polynomial approximation)
##    N_total_nodes = N_edges*(p-1) + N_nodes + ((p-1)**2)*N_elements
##    
##    # allocate memory space for global integral vector
##    F = numpy.zeros(N_total_nodes, dtype="float64")
##    
##    # go along the elements and perform the integration on each element
##    for element, nodes in enumerate(mesh.elements.nodes):
##        # allocate memory space for local integral matrix
##        F_e = numpy.zeros([p+1,p+1], dtype="float64")
##        
##        for k,deriv in enumerate(lagrange_deriv):
##            # lower_left_corner point coordinates
##            lower_left_corner = mesh.nodes.coordinates[nodes[0]]
##            # upper right corner point coordinates
##            upper_right_corner = mesh.nodes.coordinates[nodes[2]]
##            
##            # compute the x size of the node
##            element_x_size = upper_right_corner[0] - \
##                             lower_left_corner[0]
##            # compute the y size of the node
##            element_y_size = upper_right_corner[1] - \
##                             lower_left_corner[1]
##            
##            # evaluate the x lagrange polynomials at the collocation points
##            # with the derivative given by lagrange_deriv
##            phi_x_at_points = mylagrange.evaluate_deriv_at(int_collocation_points, m=deriv[0])
##            # evaluate the y lagrange polynomials at the collocation points
##            # with the derivative given by lagrange_deriv
##            psi_y_at_points = mylagrange.evaluate_deriv_at(int_collocation_points, m=deriv[1])
##            
##            # evaluate function at the  collocation points
##            # but remap the coordinates to the element coordinates
##            function_at_collocation_points = function(0.5*(int_collocation_x \
##                                           + 1.0)*element_x_size + lower_left_corner[0],\
##                                           0.5*(int_collocation_y + 1.0)*\
##                                           element_y_size + lower_left_corner[1])
##                                        
##            # compute the inner product
##            for n, phi_n in enumerate(phi_x_at_points):
##                for m, psi_m in enumerate(psi_y_at_points):
##                    temp_integral = numpy.zeros_like(function_at_collocation_points)
##                    # first multiply each row of function_at_collocation_points
##                    # by each element of gl_weights and each value of phi_n
##                    temp_integral = function_at_collocation_points * psi_m *\
##                                    gl_weights
##                    # now multiply each column of temp_integral by each value
##                    # of psi_m and each value of gl_weights
##                    for column in range(0,p_int+1):
##                        temp_integral[:,column] = temp_integral[:,column]*\
##                                                  phi_n * gl_weights
##                    
##                    # the integral is now obtained just summing all the
##                    # element of the temp_integral array
##                    integral = temp_integral.flatten().sum()
##                    
##                    # the integral must be rescaled
##                    # each derivative must be multiplied by 2/(b-a)
##                    # and the hole integral by (b-a)/2
##                    # where (b-a) is the size of the element, and this should
##                    # be made for each of the directions x and y
##                    scalling = ((2.0/element_x_size)**deriv[0]) *\
##                               ((2.0/element_y_size)**deriv[1]) *\
##                               (element_y_size/2.0)*(element_x_size/2.0)
##                    F_e[n,m] = F_e[n,m] + a[k]*integral*scalling
##        
##        # add F_e to the F vector, with the gathering matrix
##        for n in range(0,p+1):
##            for m in range(0,p+1):
##                F[mesh.global_numbering[element,n,m]] = F[mesh.global_numbering[element,n,m]] + F_e[n,m]
##    
##    return F    
        
        
def inner_product(f_0, f_1, transfinite, p_int, u_v, dx_dy):
    # Computes the inner product between the two functions: f_0 and f_1, that is,
    # computes the following integral:
    #
    #    \int\int\_{\Omega_{e}}frac{\partial^{n^{x}_{0}}}{\partial x^{n^{x}_{0}}}\frac{\partial^{n^{y}_{0}}}{\partial y^{n^{y}_{0}}}(f_{0}(x,y))\cdot \frac{\partial^{n^{x}_{1}}}{\partial x^{n^{x}_{1}}}\frac{\partial^{n^{y}_{1}}}{\partial y^{n^{y}_{1}}}(f_{1}(x,y))\mathrm{d}x\mathrm{d}y
    #
    # Where n^{x}_{i} and n^{y}_{i} are defined in dx_dy, and state the order
    # of the derivative present in f_i in the integral.
    #
    # Notice that the integral is performed in \Omega_{e}, the element, in
    # global coordinates.
    #
    # The integral is not an analytical one, instead it is a performed a numerical
    # integration using a Gauss-Lobatto quadrature or Radau quadrature in 2d.
    # see for the 1d case:
    #     http://mathworld.wolfram.com/LobattoQuadrature.html
    #     Atkinson, An Introduction to Numerical Analysis
    #
    # Hence, a transformation of coordinates is made,
    # using transfinite, in order to pass the integral to \Omega', the (u,v)
    # general coordinates. So that the Gauss Quadrature is easily computed.
    #
    # f_0 :: is a valid python fuction that takes two inputs, u and v or 
    #        x and y and returns the function f_0 evaluated at (u,v) or (x,y)
    #        the output of f_0 has the same shape as the inputs u and v and
    #        x and y.
    #
    # f_1 :: is a valid python fuction that takes two inputs, u and v or 
    #        x and y and returns the function f_1 evaluated at (u,v) or (x,y)
    #        the output of f_0 has the same shape as the inputs u and v and
    #        x and y.
    #
    # transfinite :: is a valid TransfiniteMapping_2d instance
    #
    # p_int :: since Gauss-Lobatto quadrature is used for the integration then
    #          p_int specifies the order of the quadrature used.
    #
    # u_v :: is a 2d vector of boolean values, that specifies if the functions
    #        f_0 and f_1 are defined as functions of (u,v) or as fuctions of 
    #        (x,y). Example, typically the basis functions varphi_ij(u,v) are
    #        defined in (u,v) and the forcing functions F(x,y) are specified
    #        in (x,y).
    #
    # dx_dy :: specifies the order of the derivative in x and y of the functions 
    #          whose inner product one wants to compute. This does not compute
    #          the derivative it just serves as information to the routine. The
    #          derivative of the functions should be given. Notice that this
    #          parameters will only have effect for functions defined in (u,v)
    #          not for functions defined in (x,y). For functions defined
    
    # compute the 1d collocation coordinates based on the integration order
    collocation_coordinates = collocation_points(p_int)
    
    # allocate memory space for the collocation coordinates
    # and make sure they are of type float64 and have fortran ordering
    u = numpy.zeros([p_int+1, p_int+1], dtype="float64", order="fortran")
    v = numpy.zeros([p_int+1, p_int+1], dtype="float64", order="fortran")
    
    # compute the u,v coordinates of the collocation points based on the 1d
    # collocation_coordinates
    v[:,:], u[:,:] = numpy.meshgrid(collocation_coordinates, collocation_coordinates)
    
    if u_v[0] == True:
        if dx_dy[0] == [0,0]:
            evaluated_f_0 = numpy.asfortranarray(f_0(u, v))
        elif dx_dy[0] == [1,0]:
            du_f_0 = f_0(u, v, du_dv=[1,0])
            dv_f_0 = f_0(u, v, du_dv=[0,1])
            dv_T_Y = transfinite.dY_dv(u, v)
            du_T_Y = transfinite.dY_du(u, v)
            evaluated_f_0 = numpy.asfortranarray(dv_T_Y*du_f_0 - du_T_Y*dv_f_0)
            
        else:
            du_f_0 = f_0(u, v, du_dv=[1,0])
            dv_f_0 = f_0(u, v, du_dv=[0,1])
            dv_T_X = transfinite.dX_dv(u, v)
            du_T_X = transfinite.dX_du(u, v)
            
            evaluated_f_0 = numpy.asfortranarray(du_T_X*dv_f_0 - dv_T_X*du_f_0)
    else:
        evaluated_f_0 = numpy.asfortranarray(f_0(transfinite.X(u,v), transfinite.Y(u,v)))
        
    if u_v[1] == True:
        if dx_dy[1] == [0,0]:
            evaluated_f_1 = numpy.asfortranarray(f_1(u, v))
        elif dx_dy[1] == [1,0]:
            du_f_1 = f_1(u, v, du_dv=[1,0])
            dv_f_1 = f_1(u, v, du_dv=[0,1])
            dv_T_Y = transfinite.dY_dv(u, v)
            du_T_Y = transfinite.dY_du(u, v)
            
            evaluated_f_1 = numpy.asfortranarray(dv_T_Y*du_f_1 - du_T_Y*dv_f_1)
            
        else:
            du_f_1 = f_1(u, v, du_dv=[1,0])
            dv_f_1 = f_1(u, v, du_dv=[0,1])
            dv_T_X = transfinite.dX_dv(u, v)
            du_T_X = transfinite.dX_du(u, v)
            
            evaluated_f_1 = numpy.asfortranarray(du_T_X*dv_f_1 - dv_T_X*du_f_1)
    
    else:
        evaluated_f_1 = numpy.asfortranarray(f_1(transfinite.X(u,v), transfinite.Y(u,v)))
    
    # compute the power of J to be multiplied
    J_count = (numpy.array(u_v) * numpy.array(dx_dy).T).sum()
    
    # compute J depending on J_count
    if J_count == 0:
        J = numpy.asfortranarray(numpy.abs(transfinite.J(u,v)))
    elif J_count == 1:
        J = numpy.asfortranarray(numpy.sign(transfinite.J(u,v)))
    else:
        J = numpy.asfortranarray(1.0/numpy.abs(transfinite.J(u,v)))
    
    weights = gauss_weights(p_int)
    
    # compute the inner product
    result = 0.
    result = _math_tools.inner_product(evaluated_f_0, evaluated_f_1, J, weights, result)
    return result

def compute_element_analytical_residual(element_solution, analytical_solution,\
                                        transfinite, myLagrange2d, p_int):
    p = myLagrange2d.order
    
    # given the vector with the solution on the element and the
    # function with the analytical solution on the element, compute the residual
    # with gauss-lobatto quadrature of order p_int, taking into account the
    # transfinite mapping
    
    # compute the 1d collocation coordinates based on the integration order
    collocation_coordinates = collocation_points(p_int)
    
    # allocate memory space for the collocation coordinates
    # and make sure they are of type float64 and have fortran ordering
    u = numpy.zeros([p_int+1, p_int+1], dtype="float64", order="fortran")
    v = numpy.zeros([p_int+1, p_int+1], dtype="float64", order="fortran")
    
    # compute the u,v coordinates of the collocation points based on the 1d
    # collocation_coordinates
    v[:,:], u[:,:] = numpy.meshgrid(collocation_coordinates, collocation_coordinates)

    # allocate memory space for the evaluated residual at the integration
    # collocation points
    evaluated_residual = numpy.zeros_like(u)
    
    # compute the numerical solution at the integration collocation points
    for px in range (0,p[0]+1):
        for py in range (0,p[1]+1):
            evaluated_residual = evaluated_residual + myLagrange2d[px,py](u, v) * element_solution[px,py]
    
    # compute the analytical solution at the integration collocation points
    # and subtract from the evaluated_residual, that now is just the numerical
    # solution
    evaluated_residual = numpy.asfortranarray(evaluated_residual - analytical_solution(transfinite.X(u,v),transfinite.Y(u,v)))
    
    # compute the Jacobian of the transformation of coordinates
    J = numpy.asfortranarray(numpy.abs(transfinite.J(u,v)))
    
    # compute the integration weights for Gauss-Lobatto quadrature
    weights = gauss_weights(p_int)
    
    # compute the residual integral
    result = 0.
    result = _math_tools.inner_product(evaluated_residual, evaluated_residual, J, weights, result)
    return result
    