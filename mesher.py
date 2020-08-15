## \package nereides.mesher
# \brief Mesh related classes and methods.
#
# \detailed Package for generating meshes on a given domain, for using with Nereides
# Least Squares Spectral Element flow solver.
# 
# Currently the following types of meshes are available:
#     - Structured 2d square mesh

# import NereidesException for handling with custom exceptions
import NereidesExceptions
import gauss_quadrature
# fortran routines so that things run faster
import _mesher
import math_tools
import settings

# import numpy for using optimized array objects
import numpy
import pylab

# import pysparse for using optimized sparse arrays
import pysparse

## Types of available 2d element types.
VALID_ELEMENT_TYPES_2d = ("rectangle","quad")

## Types of available mesh sources.
# \li \b generate_rectangle_2d : generates a 2d rectangular %mesh.
# \li \b read_geompack++_2d : reads %mesh data from data files generated from geompack++ mesh generator.
VALID_MESH_SOURCES = ("generate_rectangle_2d", "read_geompack++_2d")

## \class mesh 
# \brief Contains all data that defines the %mesh and all methods
# used for the construction of the various types of meshes.
#
# Two different methods for generating a mesh may be used:
#    \li \c generate_rectangle_2d : for generating a rectangular %mesh with
#                                   rectangular angles.
#    \li \c read_geompack++_2d : for reading %mesh data information from
#                                geompack++ (http://members.shaw.ca/bjoe) %mesh
#                                generator.
#
# The available names of the methods for generating a mesh are specified in 
# #VALID_MESH_SOURCES.
#
# This class stores and computes, mainly, information regarding:
#    \li nodes : see nodes2d
#    \li elements : see elements2d
#    \li edges : see edges2d
#    \li global numbering : see mesh.generate_global_numbering and 
#                           mesh.global_numbering
#
class mesh:
    # mesh class contains the following variables:
    
    ## \var dimension 
    # \brief The number of dimensions of the computational domain.
    #
    # Can take any value 
    # included in #VALID_SPACIAL_DIMENSIONS.\n\n
    # \b TYPE: \c int
    
    ## \var element_type
    # \brief The type of element used to construct the %mesh.
    # 
    # Can take any value
    # included in #VALID_ELEMENT_TYPES_2d.\n\n
    # \b TYPE: \c str
    
    ## \var N_elements
    # \brief The number of elements of the %mesh.
    # 
    # \b TYPE: \c int
    
    ## \var N_nodes
    # \brief The number of main nodes of the %mesh.
    #
    # \b TYPE: \c int
    
    ## \var N_edges
    # \brief The number of edges of the %mesh.
    #
    # \b TYPE: \c int
    
    ## \var nodes
    # \brief Object that contains data for the nodes of the %mesh.
    #
    # Instance of object of type nodes2d containing all the information related
    # to the main nodes, such as coordinates. See nodes2d for more information.
    # \n\n
    # \b TYPE: \c nodes2d
    
    ## \var edges
    # \brief Object that contains data for the edges of the %mesh.
    #
    # Instance of object of type edges2d containing all the information related
    # to the main edges of the %mesh, such as the nodes that define them. 
    # See edges2d for more information.
    # \n\n
    # \b TYPE: \c edges2d
    
    ## \var elements
    # \brief Object that contains data for the elements of the %mesh.
    #
    # Instance of object of type elements2d containing all the information
    # related to the elements of the %mesh, such as the nodes that define
    # them. See elements2d for more information.
    # \n\n
    # \b TYPE: \c elements2d
    
    ## \var N_total_nodes
    # \brief The total number of nodes, main and inner nodes.
    #
    # \b TYPE: \c int32
    
    ## \var global_numbering
    # \brief Array that contains the global number of each node, main and inner.
    #
    # It is an array that whose indices identify the node: [e,i,j], where e is
    # the element number, i is the x local index and j is the y local index of
    # the node. The value of the array at those indices is thus the global
    # index of the node.\n\n
    #
    # \b TYPE: 3d array of \c int32
    
    ## \var boundary_nodes_indices
    # \brief Array that contains the list with the indices of nodes, main and 
    # inner, that lie in a boundary with prescribed boundary conditions.
    #
    # The list contains in each row the data for each node at the
    # boundary. The row contains the index of the element, the x and y indices
    # of the node, thus uniquely identifying the node, the last column contains
    # the boundary condition identifier, specifying the boundary condition
    # prescribed on that node. This data, jointly with
    # the global numbering completely identifies the node.\n\n
    #
    # \b TYPE: array of \c int32
    
    ## \var boundary_nodes_coordinates
    # \brief Array that contains the list of coordinates of the nodes,
    # main and inner, that lie in a boundary with prescribed boundary 
    # conditions.
    #
    # The nodes here are listed in the same order as in 
    # #boundary_nodes_indices, specifying the x and y coordinates of the nodes.
    # \n\n
    # 
    # \b TYPE: array of \c float64
    
    ## \var x_spacing
    # \brief The spacing of the %mesh in the \e x direction. \b Only in rectangle
    # meshes.
    #
    # \b TYPE: \c float64
    
    ## \var y_spacing
    # \brief The spacing of the %mesh in the \e y direction. \b Only in rectangle
    # meshes.
    #
    # \b TYPE: \c float64
    
    ## \brief The constructor of the mesh object, generates the %mesh 
    # requested for the computational domain defined.
    #
    # Depending on the mesh_source input, a mesh will be generated according
    # to the requested algorithm.
    #    
    ## \param[in] mesh_source
    # Specifies the source for construction of mesh data, can take the values
    # present in #VALID_MESH_SOURCES.\n
    # Notice 
    # that a check is made to confirm that this parameter is present in
    # #VALID_MESH_SOURCES, if not the 
    # object is destroyed and an exception of the type 
    # NereidesExceptions.MeshGenerationError is raised. \n
    # \b TYPE: \c str \n
    #
    ## \param[in] **parameters
    # Consists of the list of paramerters dependent on the mesh_source selected.
    # It is a dictionary. The parameters depending on the mesh_source are:\n
    # \b generate_rectangle_2d \n
    #    \li \c mesh_size --> see \ref mesh_mesh_size.
    #    \li \c lower_left_corner --> see \ref mesh_lower_left_corner.
    #    \li \c upper_right_corner --> see \ref mesh_upper_right_corner.
    #    \li \c boundary_condition_segments --> see \ref mesh_boundary_condition_segments.
    # \n
    #
    # \b read_geompack++_2d \n
    #    \li regions_filename --> the file format used in this file can be checked in
    #                             http://members.shaw.ca/bjoe/techrep/zcs08-01.pdf
    #    \li curves_filename --> the file format used in this file can be checked in
    #                            http://members.shaw.ca/bjoe/techrep/zcs08-01.pdf
    #    \li mesh_filename --> the file format used in this file can be checked in
    #                          http://members.shaw.ca/bjoe/techrep/zcs08-01.pdf
    # \n
    #
    # \b TYPE: \c dictionary \n
    #
    ## \param[out] self.dimension 
    # See #dimension \n
    #
    ## \param[out] self.element_type 
    # See #element_type \n
    #
    ## \param[out] self.N_elements
    # See #N_elements \n
    #
    ## \param[out] self.N_nodes
    # See #N_nodes \n
    #
    ## \param[out] self.N_edges
    # See #N_edges \n
    #
    ## \param[out] self.nodes
    # See #nodes \n
    #
    ## \param[out] self.edges
    # See #edges \n
    #
    ## \param[out] self.elements
    # See #elements \n
    #
    ## \param[out] self.elements.boundary_edges
    # see #elements
    #
    ## \param[out] self.elements.on_boundary
    # see #elements
    #    
    def __init__(self, mesh_source, **parameters):
        
        if mesh_source not in VALID_MESH_SOURCES:
            # delete the object since mesh_source is an invalid value
            del(self)
            # raise an Exception of the type 
            # NereidesExceptions.MeshGenerationError
            raise NereidesExceptions.MeshGenerationError(\
                    "Mesh source not valid: " + str(mesh_source) + ". " + \
                    "Must be one of " + str(VALID_MESH_SOURCES))
            return
        
        # check if the parameters are compatible with mesh_source
        
        # parameters for generate_rectangle_2d
        if mesh_source == "generate_rectangle_2d":
            # list of needed parameters for generate_rectangle_2d source
            needed_parameters = ["mesh_size",\
                                 "lower_left_corner",\
                                 "upper_right_corner",\
                                 "boundary_condition_segments"]
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return
                    
        # parameters for read_geompack++_2d
        if mesh_source == "read_geompack++_2d":
            # list of needed parameters for  read_geompack++_2d
            needed_parameters = ["regions_filename",\
                                 "curves_filename",\
                                 "mesh_filename"]
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return        
        
        # depending on mesh_source, call the correct mesh generation function
        # which can be generate_rectangle_mesh, or  read_geompack_2d (which
        # just reads the data output from geompack++ mesh generator)
        if (mesh_source == "generate_rectangle_2d"):
            # output info to user
            print "Generating rectangle 2d mesh..."
            # set mesh dimension
            self.dimension = 2
            # set element type valid values are the ones in VALID_ELEMENT_TYPES_2d
            self.element_type = "rectangle"
            # generate 2d rectangle mesh using the input parameters
            self.__generate_rectangle_mesh(parameters["mesh_size"],\
                                         parameters["lower_left_corner"],\
                                         parameters["upper_right_corner"],\
                                         parameters["boundary_condition_segments"])
                                        
        elif (mesh_source == "read_geompack++_2d"):
            # output info to user
            print "Reading geompack++ mesh data..."
            # set mesh_dimension
            self.dimension = 2
            # set element type valid values are the ones in VALID_ELEMENT_TYPES_2d
            self.element_type = "quad"
            # generate 2d quad mesh using data present in output files of
            # geompack++
            self.__read_geompack_2d(parameters["regions_filename"],\
                                  parameters["curves_filename"],\
                                  parameters["mesh_filename"])
                                
    ## \brief Generates a rectangular %mesh.
    #
    # \anchor mesher_generate_rectangle_mesh_brief <b>Brief description</b>
    #
    # Generates a %mesh made of rectangles on the specified domain and defines
    # the boundary conditions on the segments of the boundary with prescribed
    # boundary conditions. The rectangular domain is specified with two points
    # one at the lower left corner and another at the upper right corner. 
    # Notice that this domain has boundaries parallel to the \e x and \e y axis.
    # \n
    #
    ## \anchor mesher_generate_rectangle_mesh_check <b>Parameter check</b>
    #
    ## \anchor mesher_generate_rectangle_mesh_check_mesh_size \em mesh_size
    #
    # \c mesh_size is a vector of type \c int of dimension 2, with the
    # number of elements of the computational domain in the \e x and \e y 
    # directions. A check is made to the data given by the user in order to 
    # confirm it is a vector with this shape. Additionaly the a check is made
    # to insure that the values are positive. If not, an \c Exception of the type 
    # NereidesExceptions.MeshGenerationError is raised, with all changes to the
    # object not implemented. A \c cast to \c int is made, hence if the values
    # are floating point numbers, they will be truncated and no error appears,
    # unless the obtained \c int values are invalid.
    #
    ## \anchor mesher_generate_rectangle_mesh_check_lower_left_corner \em lower_left_corner
    #
    # \c lower_left_corner is a vector of \c float of dimension 2, with the
    # \e x and \y coordinates. A check is made to the data given by the user in
    # order to confirm it is of this shape, if not, an \c Exception of the type
    # NereidesExceptions.MeshGenerationError is raised, with all changes to the
    # object not implemented. The data inputed is \c cast into a \c numpy.array
    # of \c float64
    #
    ## \anchor mesher_generate_rectangle_mesh_check_upper_right_corner \em upper_right_corner
    #
    # The same as for \c lower_left_corner.
    #
    ## \anchor mesher_generate_rectangle_mesh_check_boundary_condition_segments \em boundary_condition_segments
    #
    # \c boundary_condition_segments is an array of \c int of dimension \e n x 
    # 4 where \e n is the number of boundaries with prescribed boundary
    # conditions and the dimension 4 contains the number that idenfitifies the
    # boundary (0 to 3),
    # the start and stop indices of the elements of the corresponding
    # boundary, where the boundary condition is defined and to the boundary 
    # condition 
    # identification.
    # A check is made to confirm:
    #     \li that the array has dimensions \e n x 4, a cast to\c int is made
    #     \li that the first element of each row has values in 0,1,2 or 3.
    #     \li that the last element of each row has values in the range 0 to 
    #         \f$ N_{bc}-1 \f$, where \f$ N_{bc} \f$ is the number of different
    #         boundary conditions specified in \c boundary_conditions.
    #
    # If not, an \c 
    # Exception of the type NereidesExceptions.MeshGenerationError is raised,
    # with all changes to the object not implemented.
    #
    ## \anchor mesher_generate_rectangle_mesh_check_boundary_conditions \em boundary_conditions
    #
    # \c boundary_conditions is a vector of \c float whose dimension is the
    # number of different boundary conditions specified. This is a list of
    # valid python functions dependent on \e x and \e y. Here the check is just
    # to verify that this function names correpond to valid functions.
    # If not, an \c 
    # Exception of the type NereidesExceptions.MeshGenerationError is raised,
    # with all changes to the object not implemented.
    #
    ## \anchor mesher_generate_rectangle_mesh_algorithm \b Algorithm
    #
    # The output of this function is a mesh object. In order to obtain all the
    # data contained on such an object this function peforms the following:
    # \li Compute the \e x and \e y spacing.
    # \li Compute the main nodes, constructing a nodes2d object. This
    #     automatically numbers the main nodes of the %mesh, from left to right
    #     and from bottom to top.
    # \li Compute the edges, constructing an edges2d object.
    # \li Compute the elements, constructing an elements2d object.
    # \li Compute, for the boundary elements, which edges they contain.
    #
    ## \anchor mesher_generate_rectangle_mesh_algorithm_compute_delx_dely <em>Compute x and y spacing</em>
    #
    # The \e x and \e y lengths of the computational domain, \f$ \Delta x_{i} 
    # \f$, (\e i = 1,2, corresponds to \e x and \e y, respectively) are 
    # computed using the \e x and \e y coordinates of the \c upper_right_corner
    # and \c lower_left_corner points, in the following way:
    # \f[
    #     \Delta x_{i} = P^{ruc}_{i} - P^{llc}_{i}
    # \f]
    #
    # The \e x and \e y spacings (\f$ \delta x \f$ and \f$ \delta y \f$) are
    # computed in the following way:
    # \f[
    #     \delta x_{i} = \frac{\Delta x_{i}}{N^{e}_{x_{i}}}
    # \f]
    #    
    ## \anchor mesher_generate_rectangle_mesh_algorithm_compute_main_nodes <em>Compute the main nodes coordinates</em>
    #
    # The main nodes coordinates are generated by the creation of a nodes2d
    # object with the generate parameter set to true. This object also stores
    # all data relative to the main nodes of the %mesh. Look at nodes2d for a
    # complete explanation on the generation of the main nodes coordinates and
    # structure of the data.
    #
    ## \anchor mesher_generate_rectangle_mesh_algorithm_compute_edges <em>Compute the edges of the %mesh</em>
    #
    # The %mesh main edges are generated by the creation of an edges2d object
    # with the generate parameter set to true. This object also stores all data
    # relative to the main edges of the %mesh. Look at edges2d for a complete
    # explanation on the generation of the main edges and structure of the
    # data.
    #
    ## \anchor mesher_generate_rectangle_mesh_algorithm_compute_elements <em>Compute the elements of the %mesh</em>
    # The %mesh elements are generated by the creation of an elements2d object
    # with the generate parameter set to true. This object also stores all data
    # relative to the elements of the %mesh. Look at elements2d for a complete
    # explanation on the generation of the elements and structure of the data.
    #
    ## \anchor mesher_generate_rectangle_mesh_algorithm_compute_elements <em>Compute edges of boundary elements</em>
    #
    # Create a list with all boundary elements and  which edges they have.
    #
    ## \anchor mesher_generate_rectangle_mesh_definitions \b Definitions
    #
    #    \li \f$ \mathbf{x}^{n} \f$ represents the vector with the coordinates
    #        of a node.
    #    \li \f$ x_{i} \f$ with \f$ i=1,2 \f$ represents the \e x and \e y 
    #         coordinates, respectively.
    #    \li \f$ \Delta x_{i} \f$ represents the \e x and \e y lengths of the
    #         computational domain.
    #    \li \f$ P^{ruc} \f$ represents the Right Upper Corner point.
    #    \li \f$ P^{llc} \f$ represents the Lower Left Corner point.
    #    \li \f$ \delta x_{i} \f$ represents the spacing in the \f$ x_{i} \f$
    #        direction.
    #    \li \f$ N^{e}_{x_{i}} \f$ represents the number of elements in the 
    #        \f$ x_{i} \f$ direction. That is, the \c mesh_size elements.
    #
    ## \anchor mesh_generate_rectangle_mesh_parameters
    #
    ## \param[in] mesh_size  \anchor mesh_mesh_size
    # A two dimensional vector, where the first
    # value specifies the number of elements in the \e x
    # direction and the second value specifies the
    # number of elements in the \e y direction. \n
    # \b TYPE: vector of \c int of size 2.
    #
    ## \param[in] lower_left_corner \anchor mesh_lower_left_corner
    # A two dimensional vector with the \e x and \e y
    # coordinates of the lower left corner point
    # that defines the domain.\n
    # \b TYPE: vector of \c float of size 2.
    #
    ## \param[in] upper_right_corner \anchor mesh_upper_right_corner
    # A two dimensional vector with the x and y
    # coordinates of the upper right corner 
    # point that defines the domain.\n
    # \b TYPE: vector of \c float of size 2.
    #
    ## \param[in] boundary_condition_segments \anchor mesh_boundary_condition_segments
    # Defines the boundary segments where boundary conditions are prescribed.
    # It can specify any arbitrary number of different boundary conditions.
    # It is an array where each row contains 4 values:
    #     \li id number of the face where the boundary condition segment is
    #         prescribed, takes values from 0 to 3, where 0 is the bottom, 1
    #         left, 2 top and 3 right.
    #     \li start node of boundary condition segment
    #     \li end node of boundary condition segment
    #     \li boundary condition id, a number that specifies which boundary
    #         condition is to be used, it relates to the boundary_conditions
    #         vector defined next. Takes values from 0 to size of 
    #         boundary_conditions - 1.
    #
    # \b TYPE: array of \c int of size \e n x 4, where \e n is the number of
    # boundary_segments.
    #
    ## \param[out] self.N_elements
    # see #N_elements
    #
    ## \param[out] self.N_nodes
    # see #N_nodes
    #
    ## \param[out] self.N_edges
    # see #N_edges
    #
    ## \param[out] self.boundary_conditions
    # see #boundary_conditions
    #
    ## \param[out] self.x_spacing
    # see #x_spacing
    #
    ## \param[out] self.y_spacing
    # see #y_spacing
    #
    ## \param[out] self.nodes
    # see #nodes
    #
    ## \param[out] self.edges
    # see #edges
    #
    ## \param[out] self.elements
    # see #elements
    #
    ## \param[out] self.elements.boundary_edges
    # see #elements
    #
    ## \param[out] self.elements.on_boundary
    # see #elements
    #
    def __generate_rectangle_mesh(self, mesh_size, \
            lower_left_corner, upper_right_corner, \
            boundary_condition_segments=None):
        
        # check input parameter validity
        
        # check mesh_size
        # must be a vector of dimension 2 of type int
        # convert to numpy.array and then check the shape
        mesh_size = numpy.array(mesh_size, dtype="int")
        if mesh_size.shape != (2,):
            # since it does not have the correct shape
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError("Mesh size shape not \
                valid:" + str(mesh_size.shape) + ", must be (2,)")
            return
        
        # check if the values are negative
        if mesh_size.min() <= 0:
            # since there are negative values
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError("Mesh size vector \
                cannot have negative values or zero:" + str(mesh_size))
            return
        
        # check lower_left_corner
        # must be a vector of dimension 2 of type float64
        # convert to numpy array of float64 and check the shape
        lower_left_corner = numpy.array(lower_left_corner, dtype="float64")
        if lower_left_corner.shape != (2,):
            # since it does not have the correct shape
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError("lower_left_corner \
                shape not valid:" + lower_left_corner.shape + ", must be (2,)")
            return
        
        # check upper_right_corner
        # must be a vector of dimension 2 of type \c float64
        # convert to numpy array of \c float64 and check the shape
        upper_right_corner = numpy.array(upper_right_corner, dtype="float64")
        if upper_right_corner.shape != (2,):
            # since it does not have the correct shape
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError("upper_right_corner \
                shape not valid:" + upper_right_corner.shape + ", must be (2,)")
            return
        
        # check boundary_condition_segments
        # convert to numpy.array check that has dimension n x 4 where n is
        # the dimension of the boundary_conditions vector
        boundary_condition_segments = numpy.array(boundary_condition_segments,\
                dtype="int64")
        if boundary_condition_segments.shape[1] != 4:
            # since it does not have the correct shape
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError(\
                    "boundary_condition_segments shape not valid:" + \
                    str(boundary_condition_segments.shape) + ", must be (%d,4)" % \
                    boundary_conditions.shape[0])
            return
        
        if boundary_condition_segments[:,-1].min() < 0:
            # since it has invalid values
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError(\
                    "boundary_condition_segments boundary condition \
                    identifiers negative:" + \
                    str(boundary_condition_segments[:,-1].min()) + ", must be >= 0")
            return
        
        # check if the boundary identifiers, the first column are in the range
        # 0 to 3
        if boundary_condition_segments[:,0].max() > 3:
            # since it has invalid values
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError(\
                    "boundary_condition_segments boundary \
                    identifiers too big:" + \
                    str(boundary_condition_segments[:,0].max())\
                    + ", must be < 4")
            return
        if boundary_condition_segments[:,0].min() < 0:
            # since it has invalid values
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError(\
                    "boundary_condition_segments boundary   \
                    identifiers negative:" + \
                    boundary_condition_segments[:,0].min() + ", must be > 0")
            return
        
        # compute the number of elements, number of nodes and number of edges
        N_elements = mesh_size[0]*mesh_size[1]
        self.N_elements = N_elements
        N_nodes = (mesh_size[0]+1)*(mesh_size[1]+1)
        self.N_nodes = N_nodes
        self.N_edges = N_elements + N_nodes - 1
        
        # compute the x and y spacing
        
        # first compute the x and y lengths of the computational domain
        x_length = upper_right_corner[0] - lower_left_corner[0]
        y_length = upper_right_corner[1] - lower_left_corner[1]
        
        # compute the x and y spacing
        x_spacing = x_length / mesh_size[0]
        y_spacing = y_length / mesh_size[1]
        
        self.x_spacing = x_spacing
        self.y_spacing = y_spacing
        
        # compute nodes
        print "computing nodes..."
        self.nodes = nodes2d("generate_rectangle_2d",\
                             spacings=numpy.array([x_spacing,y_spacing]),\
                             mesh_size=mesh_size,\
                             lower_left_corner=lower_left_corner,\
                             upper_right_corner=upper_right_corner,\
                             boundary_condition_segments=boundary_condition_segments)
                                            
        # compute the edges
        print "computing edges..."
        self.edges = edges2d("generate_rectangle_2d", N_nodes=N_nodes, N_elements=N_elements, mesh_size=mesh_size,\
                                            boundary_condition_segments=boundary_condition_segments,\
                                            element_type=self.element_type)
                                            
        # compute the elements
        print "computing elements..."
        self.elements = elements2d("generate_rectangle_2d", N_elements=N_elements,\
                              mesh_size=mesh_size,\
                              element_type=self.element_type)
                            
        # determine in which elements the nodes are
        #print "finding nodes in elements..."
        #self.__find_nodes()
        
        # determine which edges are in the boundary element
        print "finding boundary edges in boundary elements..."
        self.__generate_boundary_element_edges()
        
    ## \brief Reads %mesh data from a geompack++ generated %mesh.
    #
    # \anchor mesher_read_geompack_2d_brief <b>Brief description</b>
    #
    # Reads %mesh information from a geompack generated %mesh, comprising:
    #    \li nodes : x and y coordinates of nodes of the %mesh.
    #    \li elements : list of nodes that specify each element.
    #    \li edges : list of edges definitions, with start and end point.
    #    \li boundary conditions : identification of the curves that make up
    #                              the boundary of the domain.
    #
    # These informations are specified, as can be seen in 
    # http://members.shaw.ca/bjoe/techrep/zcs08-01.pdf, in three different
    # files:
    #    \li regions file : This file contains information on the definition
    #                       of the regions of the computational domain. Contains
    #                       the nodes that define the geometry of the boundaries.
    #                       And also a list with the nodes that define a boundary.
    #    \li curves file : This file contains information on the definition of
    #                      all the curves that make up the computational domain.
    #                      This curves use the nodes specified in the regions file.
    #                      Three types of curves can be used: straight lines,
    #                      circles and NURBS.
    #    \li %mesh file : This file contains information on the %mesh. Mainly
    #                     a list of the nodes that define the %mesh, a list of
    #                     the nodes that make up each element, a list of the 
    #                     curves to which every created node belongs and a list 
    #                     of the curves to which each edge belongs, in order to
    #                     recover the boundary conditions and the correct shape
    #                     of the edge.
    # \n
    #
    ## \anchor mesher_read_geompack_2d_check <b>Parameter check</b>
    #
    # A check is made to insure that the filenames for the regions, curves and
    # %mesh files gven are strings. If not a NereidesExceptions.MeshGenerationError
    # exception is raised, specifying the filename that is incorrect.
    #
    ## \anchor mesher_read_geompack_2d_algorithm <b>Algorithm</b>
    #
    # The information that comes as output from geompack++ is not as separated
    # as, for example, in Triangle from Shewchuk, see 
    # http://www.cs.cmu.edu/~quake/triangle.html. In triangle each file contains
    # the required information for each part of the data, that is, there is a
    # nodes file, an elements file, an edges file, a neighbors file, etc.. They
    # are dependent on each other but all the information for, say, the nodes
    # is all on the nodes file. \n
    # 
    # On geompack++ that is not the case. The information is spread on all
    # files, being necessary to read  information for all of them to be able to
    # reconstruct information for the nodes, for example. The algorithm used
    # to read the data from the three files is highly dependent on the file 
    # formats used by geompack for each file. Here we will not specify the file
    # formats very extensively. We will try to restrict to the minimum possible
    # needed to correctly understand the code. For more information on the file
    # formats one recomends the reading of 
    # http://members.shaw.ca/bjoe/techrep/zcs08-01.pdf . \n
    #
    # Shortly, the user to generate a %mesh in geompack++ gives the program a list
    # of nodes, specifying the type of each node, they can be varied but basically
    # they can be nodes that will appear in the %mesh as as they are given by 
    # the user, that is, keeping the position, or can appear in the %mesh but at
    # a different position along the line they are used to define or may not
    # appear at all, if they are used only as auxiliary points to build the
    # lines. Also the user pecifies regions of the domain, defining the nodes
    # that lie on the boundary of each region. Also the type of region is 
    # specified, being it for instance a hole in the domain. \n
    #
    # Additionaly the user defines how certain nodes are connected, specifying
    # the shape of the boundaries and the boundary condition at that curve.
    # Three types of curves (connections between nodes) can be used:
    # straight lines, circles and NURBS (http://en.wikipedia.org/wiki/NURBS).
    # For now, this code only recognizes straight lines. \n
    #
    # With this information and another file, not necessary for this code, but
    # that parametrizes the %mesh, indicating, for instance, the number of
    # elements or the minimum or maximum size of the elements, geompack++
    # outputs a %mesh file. This file contains a list of all the nodes present
    # at the %mesh, that is the x and y coordinates of the nodes and the code
    # of the node. This code is of the form: 20*xtraptr + vertyp, where xtraptr
    # is the index of the extra data list that corresponds to this node (a node
    # might not have any extra data, hence xtraptr = 0) and vertyp is the type
    # of node, specifying how the node was created, if it is a node that came
    # from a node in the regions file that could not be moved, if it is a node
    # from the regions file that was moved or if it is a node that was created.
    # More subtleties exist, but for that we recommend reading the pdf file
    # mentioned above.\n
    #
    # The extra information for the nodes specifies to which curve the node
    # belongs. The nodes that do not have extra information, means that they
    # where specified at the regions file, hence the data from that file can be
    # used for the boundary markers.\n
    #
    # After this data, comes the elements data that specify for each element
    # the index of the nodes of the %mesh that define the element. Following it,
    # a list of the code that specifies a code for each edge of each element.
    # This code is of the following form: +-(20*indcur + edgtyp), where indcur
    # is the index of the curve to which the edge belongs and edgtyp is the
    # type of edge. Edgtyp is not relevant for now. Indcur is used to define
    # the boundary condition markers of the edges of the %mesh.
    #
    # The following lines, try to specify how the extration of this data is made,
    # in order to obtain:
    #    \li \c coordinates : a list of the x and y coordinates of the nodes of
    #                         the %mesh.
    #    \li \c nodes_bc_markers : a list of the boundary condition markers of
    #                              the nodes of the %mesh.
    #    \li \c element_nodes : a list of the indices of the nodes that define
    #                           the element of the %mesh.
    #    \li \c edges_nodes : a list of the nodes that define the edges of the
    #                         %mesh.
    #    \li \c edges_bc_markers : a list of the boundary condition markers of
    #                              the edges of the %mesh.
    # \n
    #
    # The strategy followed is the following:\n
    # 
    # From the regions file read only the data specifying the type of node,
    # filling the first row of a N_nodes_region x 2 array, called regions_nodes.
    # After this regions file is closed since no more data from it is needed. \n
    #
    # From the curves file read for each curve the type of curve, the boundary
    # condition marker (store in curves_bc_markers) and the list of nodes that
    # specify the curve (nodes_on_curve). Using the curves_bc_markers and
    # nodes_on_curve one can specify which node of the regions file belongs to
    # which curve. \b Notice that if a node belongs to two curves the bc_marker
    # it will have wis the onde from the last curve read, since the data is
    # overwritten. The curves file is closed since no more data is needed. \n
    #
    # From the mesh file fill in coordinates array with the x and y coordinates
    # of the nodes of the %mesh and the nodes_xtra_info_indices, specifying
    # the index of extra info data that corresponds to each node. \b Notice
    # that not all nodes have extra info. The onde defined on the regions file
    # do not need extra info. For these points -1 is assigned as the index.\n
    #
    # Next fill in the nodes_bc_markers. \b Notice that the first nodes of the
    # %mesh are the onde fo the regions file of type 2, exactly in the same 
    # order. Hence, for these the second column of the regions_nodes array is
    # used to fill in the nodes_bc_markers of the first nodes of the %mesh.\n
    #
    # After this, extra info contains data for the nodes at the boundary. Store
    # this data in nodes_xtra_info. After that,
    # loop over nodes_xtra_info_indices data and use it as an index for
    # nodes_xtra_info to fill the nodes_bc_markers of the nodes at the boundary
    # that are not regions nodes. All the other will have a boundary condition
    # marker equal to zero.\n
    #
    # For getting the edges information, one reads the data from the %mesh file
    # that states the curve to which each edge of the element belongs to. In
    # order to overcome the overlapping of edge information for two neighboring
    # elements and avoiding a search to see if the edge was already introduced,
    # one identifies each edge by the indices of the start and end nodes and
    # the boundary condition marker of the edge are stored in a matrix, whose
    # indices are the start and end indices of the edge. This marix is sparse,
    # hence a sparse matrix is used to reduce the amount of memory used. Also
    # this matrix is symmetrical, hence both order of the start and end edges is
    # used, so that each edge is always defined by (start_edge, end_edge) and
    # (end_edge, start_edge).\n
    #
    # To get back the edges_bc_markers one just takes the non-zero data of the
    # sparse matrix using the data() method. The start and end nodes of the
    # edges can be recovered using the keys() method. Afterwards, one just takes
    # the cases where start_edge < end_edge, so that one does not have the edges
    # repeated. The bc_markers of the edges that are -1 (not on a boundary) are
    # converted to 0. They are initialized as -1 since if they where 0 they
    # would be lost in the sparse matrix. Negative, because there is no
    # possibility of getting a negative boundary condition marker for a curve.\n
    #
    # Finally, with:
    #    \li \c coordinates
    #    \li \c nodes_bc_markers
    #    \li \c element_nodes
    #    \li \c edges_nodes
    #    \li \c edges_bc_markers
    # \n
    #
    # One just generates the nodes (nodes2d), elements (elements2d) and edges
    # (edges2d) objects of the %mesh and compute the boundary_element_edges
    # (#__generate_boundary_element_edges).
    #
    ## \anchor mesher_read_geompack_2d_parameters <b>Parameters</b>
    #
    ## \param[in] regions_filename
    # The filename of a correct geompack++ regions file, see 
    # http://members.shaw.ca/bjoe/techrep/zcs08-01.pdf for more information on
    # the file format.
    #
    ## \param[in] curves_filename
    # The filename of a correct geompack++ cruves file, see 
    # http://members.shaw.ca/bjoe/techrep/zcs08-01.pdf for more information on
    # the file format.
    #
    ## \param[in] mesh_filename
    # The filename of a correct geompack++ %mesh file, see 
    # http://members.shaw.ca/bjoe/techrep/zcs08-01.pdf for more information on
    # the file format.
    #
    ## \param[out] self.nodes
    # See #nodes.
    #
    ## \param[out] self.elements
    # See #elements.
    #
    ## \param[out] self.edges
    # See #edges.
    #
    ## \param[out] self.edges
    # See #edges.
    #
    ## \param[out] self.elements.boundary_edges
    # see #elements
    #
    ## \param[out] self.elements.on_boundary
    # see #elements
    #
    def __read_geompack_2d(self, regions_filename, curves_filename, mesh_filename):
        
        # check if regions filename is a proper filename and file exists
        if (type(regions_filename) != str):
            # since it is not a string
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError("Regions filename is \
                not valid:" + str(type(regions_filename)) + ", must be a string")
            return
        
        # check if curves filename is a proper filename and file exists
        if (type(curves_filename) != str):
            # since it is not a string
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError("Curves filename is \
                not valid:" + str(type(curves_filename)) + ", must be a string")
            return
        
        # check if mesh filename is a proper filename and file exists
        if (type(mesh_filename) != str):
            # since it is not a string
            # raise an exception of the type 
            # NereidesExceptions.MeshGenerationError
            del(self)
            raise NereidesExceptions.MeshGenerationError("Mesh filename is \
                not valid:" + str(type(mesh_filename)) + ", must be a string")
            return
        
        # read nodes data from geompack++ mesh generator
        # open regions file
        regions_file = open(regions_filename, "r")
        # read the number of nodes
        N_nodes_regions = int(regions_file.readline().strip())
        # allocate memory space for the nodes of the regions
        # here we allocate memory space for an array with two columns and
        # N_nodes_region rows. The two columns are for the type of node
        # (type 2 means it is a node to appear in the mesh, this ones are the
        # ones of importance, on what regards the boundary condition markers)
        # and for the bc_marker (boundary condition marker)
        regions_nodes = numpy.zeros([N_nodes_regions, 2], dtype="int32")
        # loop over the regions nodes and store the information regarding the
        # type of node
        for node in range(0, N_nodes_regions):
            regions_nodes[node,0] = int(regions_file.readline().strip().split()[2])
        
        # close regions file
        regions_file.close()
        
        # open curves file
        curves_file = open(curves_filename, "r")
        # read the number of curves
        N_curves = int(curves_file.readline().strip())
        # allocate memory space for curves boundary condition marker
        curves_bc_markers = numpy.zeros(N_curves, dtype="int16")
        # loop over the curves and fill in curves_bc_markers and the
        # bc_markers of the nodes that define the curves0
        for curve in range(0,N_curves):
            # ignore the empty line
            garbage = curves_file.readline()
            # read the curve type and the bc_marker of the curve
            curve_type, curves_bc_markers[curve] =  [int(value) for value in \
                                                   curves_file.readline().strip().split()]
            # read the nodes numbers that define the curve
            nodes_on_curve = [int(value) for value in \
                                                   curves_file.readline().strip().split()]
            # set the bc_markers of these nodes to the bc_marker of 
            # the curve
            # the -1 is put because nodes in geompack++ start on 1
            # and in python start on 0
            for node in nodes_on_curve:
                # 1 is added since internally 0 means that node is not at the boundary
                # and 1 means it is at the boundary but that no prescribed boundary
                # conditions exist, and >=2 means that the node has prescribed
                # boundary conditions, see nodes2d.bc_markers documentation for
                # more information. The user when defines the curves in geompack++
                # sets 0 for no prescribed boundary condition and 1,2, etc for the
                # bc_marker of the prescribed boundary condition at the curve, hence
                # the +1.
                regions_nodes[node-1,1] = curves_bc_markers[curve] + 1
                
        # close curves file
        curves_file.close()
        
        # drop all nodes that make up the curves that are not of type 2, that
        # is, nodes that appear on the mesh, all the other do not show up on
        # the mesh, hence no need to keep track of them
        regions_nodes = regions_nodes[regions_nodes[:,0]==2]
        
        # open mesh file
        mesh_file = open(mesh_filename, "r")
        # read the number of nodes
        N_nodes = int(mesh_file.readline().strip())
        self.N_nodes = N_nodes
        
        # allocate memory space for nodes coordinates
        coordinates = numpy.zeros([N_nodes,2], dtype="float64")
        # allocate memory space for nodes boundary condition markers
        nodes_bc_markers = numpy.zeros(N_nodes, dtype="int16")
        # nodes extra info indices
        nodes_extra_info_indices = numpy.zeros(N_nodes, dtype="int32")
        # loop over nodes and store the x and y coordinates
        for node in range(0, N_nodes):
            nodes_data = numpy.array([float(value) for value in mesh_file.readline().strip().split()])
            coordinates[node] = nodes_data[0:2]
            # the indices of the extra data for the nodes 
            # is given as 20*xtraptr + vertyp
            # we are only interested in xtraptr, hence we make an integer
            # division by 20 and subtract 1 since python indices start at 0 and
            # geompack++ indices start at 1. -1 in the index, means that the
            # node does not have extra info
            nodes_extra_info_indices[node] = (int(nodes_data[2])/20) - 1
        
        # the regions nodes of type 2, stored in regions_nodes are the first
        # nodes of the mesh, in the exact same order. Use the bc_markers of the
        # regions_nodes, to fill the bc_markers data of the nodes of the mesh
        for node, bc_marker in enumerate(regions_nodes[:,1]):
            nodes_bc_markers[node] = bc_marker
        
        # keep track of the last node whose bc_marker was updated
        # add 1 since we want to start on the on the following node, not again
        # on the last node
        last_node = node + 1
        
        # now that all curves' data were read and all the data relative to
        # the nodes that define the curve has been read, read the following
        # group of data from the mesh file, which contains data regarding
        # the nodes on the boundary that are not nodes present on the regions
        # file. That is, this nodes are nodes that were created during the
        # meshing process.
        
        # ignore the blank line
        garbage = mesh_file.readline()
        
        # read the number of nodes to which this data is related
        N_nodes_boundary = int(mesh_file.readline().strip())
        
        # read the data relative to each node and depending on the which
        # curve the node belongs, give the node the corresponding bc_marker
        # of the curve it belongs to
        
        # allocate memory space for nodes_xtra_info data
        nodes_xtra_info = numpy.zeros(N_nodes_boundary, dtype="int16")
        # loop over xtra infor of nodes
        for node in range(0, N_nodes_boundary):
            # the second value on each line is the only one that matters.
            # 1 is added since internally 0 means that node is not at the boundary
            # and 1 means it is at the boundary but that no prescribed boundary
            # conditions exist, and >=2 means that the node has prescribed
            # boundary conditions, see nodes2d.bc_markers documentation for
            # more information. The user when defines the curves in geompack++
            # sets 0 for no prescribed boundary condition and 1,2, etc for the
            # bc_marker of the prescribed boundary condition at the curve, hence
            # the +1.
            # The -1 is placed because the curves are numbered begining on 1
            # and python indexes the values on an array starting on 0, hence
            # curve 1 is on position 0 of curve_bc_markers.
            
            nodes_xtra_info[node] = curves_bc_markers[int(mesh_file.readline().strip().split()[1])-1] + 1
        
        # now loop over the nodes and store the info of the ones that have
        # extra data
        for node in range(0, N_nodes):
            if nodes_extra_info_indices[node] >= 0:
                nodes_bc_markers[node] = nodes_xtra_info[nodes_extra_info_indices[node]]
                
        # read the element data of mesh file
        # ignore empty line
        garbadge = mesh_file.readline()
        # read number of nodes per element and number of elements
        N_nodes_per_element, N_elements = [int(value) for value in mesh_file.readline().strip().split()]
        # set N_elements of mesh
        self.N_elements = N_elements
        # allocate memory space for element data
        element_nodes = numpy.zeros([N_elements, N_nodes_per_element], dtype=settings.index_INT_TYPE)
        # read nodes data and store in element_nodes
        for element in range(0, N_elements):
            element_nodes[element] = [int(value) for value in mesh_file.readline().strip().split()]
        
        # subtract 1 since geompack++ start numbering from 1 and python
        # start from 0
        element_nodes = element_nodes - 1
        
        # allocate memory space for sparse matrix that will store the bc_markers
        # of the edges and the edges information, that is, start and end nodes
        # the sparse matrix has 2*N_edges elements, since each edge is introduced
        # twice, one for node_0 - node_1 and another for node_1 - node_2
        # but in the end only one of the two will be used, specifically the one
        # where start_node = node_0 and end_node = node_1 if node_0 < node_1
        # or start_node = node_1 and end_node = node_0 if node_1 < node_0
        sparse_edges_bc_markers = pysparse.spmatrix.ll_mat(N_nodes, N_nodes, 2*(N_nodes+N_elements))
        
        # go along the elements and fill in the edges bc_markers
        #
        # Notice that for geompack++ edge i of an element is the one that is
        # defined by elements i and i+1.
        # Edges and introduced twice, for the normal order of the edges and for
        # the other order, that is, (i, i+1) and (i+1, i), so that we do not
        # missh anything. This is done in this way since 
        # pysparse.spmatrix.ll_mat_sym does not have implemented the keys()
        # method, hence one has to use a full sparse matrix, despite it is
        # symmetric.
        
        # ignore the empty line
        garbadge = mesh_file.readline()
        
        # read the edges data
        for element in range(0, N_elements):
            # read only the data that specifies to which line the edge belongs
            edges_data = [abs(int(value)) for value in mesh_file.readline().strip().split()[1:]]
            # data that specifies to which curve the node belongs is stored
            # in the following way:
            # +-(20*indcur + edgtyp) --> see geompack++ File formats for Regions
            #                            and Meshes for more information
            # we are only interested here in the curve index, that is, to
            # which curve does the edge belong. To get it we just need to divide
            # the edge_data by 20, and be carefull for both to be ints, in this
            # way we will get the integer division, which is what we want.
            
            # what is done below could have been done in a for loop, looping
            # over the edges. This would be more elegant, but slower, hence
            # since the number of edges is small, we just do the loop explicitly
            # by coding all the operations for each edge instead of the for loop.
            
            # edge 0
            
            # compute the curve to which the edge belongs
            # -1 is put because geompack++ start indexing curves at 1
            # and python indexes vectors starting at 0
            curve_index = edges_data[0]/20 - 1
            # check if curve_index >= 0
            # if not, it means that the edge is not at a boundary, hence just
            # set the bc_marker to 0, that is, do nothing, if true  set the 
            # bc_marker of the edge to the bc_marker of the curve
            # 1 is added since internally 0 means that node is not at the boundary
            # and 1 means it is at the boundary but that no prescribed boundary
            # conditions exist, and >=2 means that the node has prescribed
            # boundary conditions, see nodes2d.bc_markers documentation for
            # more information. The user when defines the curves in geompack++
            # sets 0 for no prescribed boundary condition and 1,2, etc for the
            # bc_marker of the prescribed boundary condition at the curve, hence
            # the +1.
            if curve_index >= 0:
                sparse_edges_bc_markers[element_nodes[element,0], element_nodes[element,1]] =  curves_bc_markers[curve_index] + 1
                sparse_edges_bc_markers[element_nodes[element,1], element_nodes[element,0]] =  curves_bc_markers[curve_index] + 1
            else:
                # assign -1 to the bc_marker of points that do not lie o a curve
                # this is done and not 0 since if 0 is assigned, later these
                # edgs will not be recoverable
                sparse_edges_bc_markers[element_nodes[element,0], element_nodes[element,1]] = -1
                sparse_edges_bc_markers[element_nodes[element,1], element_nodes[element,0]] = -1
            
            # edge 1
            
            # compute the curve to which the edge belongs
            # -1 is put because geompack++ start indexing curves at 1
            # and python indexes vectors starting at 0
            curve_index = edges_data[1]/20 - 1
            # check if curve_index >= 0
            # if not, it means that the edge is not at a boundary, hence just
            # set the bc_marker to 0, that is, do nothing, if true  set the 
            # bc_marker of the edge to the bc_marker of the curve
            # 1 is added since internally 0 means that node is not at the boundary
            # and 1 means it is at the boundary but that no prescribed boundary
            # conditions exist, and >=2 means that the node has prescribed
            # boundary conditions, see nodes2d.bc_markers documentation for
            # more information. The user when defines the curves in geompack++
            # sets 0 for no prescribed boundary condition and 1,2, etc for the
            # bc_marker of the prescribed boundary condition at the curve, hence
            # the +1.
            if curve_index >= 0:
                sparse_edges_bc_markers[element_nodes[element,1], element_nodes[element,2]] =  curves_bc_markers[curve_index] + 1
                sparse_edges_bc_markers[element_nodes[element,2], element_nodes[element,1]] =  curves_bc_markers[curve_index] + 1
            else:
                # assign -1 to the bc_marker of points that do not lie o a curve
                # this is done and not 0 since if 0 is assigned, later these
                # edgs will not be recoverable
                sparse_edges_bc_markers[element_nodes[element,1], element_nodes[element,2]] = -1
                sparse_edges_bc_markers[element_nodes[element,2], element_nodes[element,1]] = -1
            
            # edge 2
            
            # compute the curve to which the edge belongs
            # -1 is put because geompack++ start indexing curves at 1
            # and python indexes vectors starting at 0
            curve_index = edges_data[2]/20 - 1
            # check if curve_index >= 0
            # if not, it means that the edge is not at a boundary, hence just
            # set the bc_marker to 0, that is, do nothing, if true  set the 
            # bc_marker of the edge to the bc_marker of the curve
            # 1 is added since internally 0 means that node is not at the boundary
            # and 1 means it is at the boundary but that no prescribed boundary
            # conditions exist, and >=2 means that the node has prescribed
            # boundary conditions, see nodes2d.bc_markers documentation for
            # more information. The user when defines the curves in geompack++
            # sets 0 for no prescribed boundary condition and 1,2, etc for the
            # bc_marker of the prescribed boundary condition at the curve, hence
            # the +1.
            if curve_index >= 0:
                sparse_edges_bc_markers[element_nodes[element,2], element_nodes[element,3]] =  curves_bc_markers[curve_index] + 1
                sparse_edges_bc_markers[element_nodes[element,3], element_nodes[element,2]] =  curves_bc_markers[curve_index] + 1
            else:
                # assign -1 to the bc_marker of points that do not lie o a curve
                # this is done and not 0 since if 0 is assigned, later these
                # edgs will not be recoverable
                sparse_edges_bc_markers[element_nodes[element,2], element_nodes[element,3]] = -1
                sparse_edges_bc_markers[element_nodes[element,3], element_nodes[element,2]] = -1
            
            # edge 3
            
            # compute the curve to which the edge belongs
            # -1 is put because geompack++ start indexing curves at 1
            # and python indexes vectors starting at 0
            curve_index = edges_data[3]/20 - 1
            # check if curve_index >= 0
            # if not, it means that the edge is not at a boundary, hence just
            # set the bc_marker to 0, that is, do nothing, if true  set the 
            # bc_marker of the edge to the bc_marker of the curve
            # 1 is added since internally 0 means that node is not at the boundary
            # and 1 means it is at the boundary but that no prescribed boundary
            # conditions exist, and >=2 means that the node has prescribed
            # boundary conditions, see nodes2d.bc_markers documentation for
            # more information. The user when defines the curves in geompack++
            # sets 0 for no prescribed boundary condition and 1,2, etc for the
            # bc_marker of the prescribed boundary condition at the curve, hence
            # the +1.
            if curve_index >= 0:
                sparse_edges_bc_markers[element_nodes[element,3], element_nodes[element,0]] =  curves_bc_markers[curve_index] + 1
                sparse_edges_bc_markers[element_nodes[element,0], element_nodes[element,3]] =  curves_bc_markers[curve_index] + 1
            else:
                # assign -1 to the bc_marker of points that do not lie o a curve
                # this is done and not 0 since if 0 is assigned, later these
                # edgs will not be recoverable
                sparse_edges_bc_markers[element_nodes[element,3], element_nodes[element,0]] = -1
                sparse_edges_bc_markers[element_nodes[element,0], element_nodes[element,3]] = -1
                
        # now that all edges data has been read one just needs to organize it
        # organize edges_bc_markers
        # the edges' nodes are just the keys of the sparse_edges_bc_markers
        # sparse matrix and the edges_bc_makers are just the values of that matrix
        edges_nodes = numpy.array(sparse_edges_bc_markers.keys(), dtype="int32")
        edges_bc_markers = numpy.array(sparse_edges_bc_markers.values(), dtype="int32")
        # note that, as stated above, the edges are doubly stated, hence we
        # remove the onde that are repeated, and keep only the was whose
        # start_edge < end_edge
        edges_bc_markers = edges_bc_markers[edges_nodes[:,0] < edges_nodes[:,1]]
        edges_nodes = edges_nodes[edges_nodes[:,0] < edges_nodes[:,1]]
        
        # change the bc_markers that are -1 to 0
        edges_bc_markers[edges_bc_markers == -1] = 0
        
        # now compute the number of edges
        self.N_edges = edges_bc_markers.size
        
        # update data on nodes2d object of self    
        self.nodes = nodes2d("read_geompack++_2d",\
                             coordinates=coordinates,\
                             bc_markers=nodes_bc_markers)
        
        # update data on elements2d object of self
        self.elements = elements2d("read_geompack++_2d",\
                                   nodes=element_nodes,\
                                   edges_nodes=edges_nodes)
        
        # update data on edges2d object of self
        self.edges = edges2d("read_geompack++_2d",\
                             nodes=edges_nodes,\
                             bc_markers=edges_bc_markers)
        
        # determine which edges are in the boundary element
        print "finding boundary edges in boundary elements..."
        self.__generate_boundary_element_edges()
        
    ## \brief Generates the global numbering of the %mesh's nodes.
    #
    ## \anchor mesher_generate_global_numbering_brief <b>Brief description</b>
    #
    # This method, usable after the generation or loading of a %mesh, generates
    # the global numbering of all the nodes of the %mesh. It produces:
    #    \li global numbering to the nodes of the %mesh
    #    \li list of boundary nodes with prescribed boundary conditions
    #    \li list of coordinates of boundary nodes with prescribed boundary
    #        conditions
    #
    # Global node numbering is done, shortly, along the elements 
    # (0 \f$\rightarrow\f$ \c N_elements), and in the following order:
    #    \li main nodes at the boundary
    #    \li edges nodes at the boundary
    #    \li inner main nodes
    #    \li inner edges nodes
    #    \li inner nodes
    #
    # Nodes at boundaries with prescribed boundary conditions are put at the
    # end, after all the other nodes, in the following order:
    #    \li main nodes
    #    \li edges nodes
    #
    # #boundary_nodes_indices keeps the order in which the boundary nodes with
    # prescribed boundary conditions are numbered, hence the first numbered
    # boundary node is the first on the list and the last one is the last. The
    # same happens for #boundary_nodes_coordinates.
    #
    ## \anchor mesher_generate_global_numbering_check <b>Parameter check</b>
    #
    # No parameter check is performed.
    #
    ## \anchor mesher_generate_global_numbering_algorithm \b Algorithm
    #
    # As described above, the global numbering is to be performed in such a way
    # such that the nodes at boundaries with prescribed boundary conditions are
    # numbered last. Hence, shortly the algorithm used consists of:
    #    \li Number the main nodes at the boundary
    #    \li Number the nodes at boundary edges
    #    \li Number the inner main nodes
    #    \li Number the inner nodes at the edges
    #    \li Number the inner element nodes
    #
    # At first, we start with the following numbering of the elements and main
    # nodes:
    # \image html global_numbering_0.jpg
    # \image latex global_numbering_0-crop.eps "Starting numbering of elements and main nodes" width=10cm
    #
    ## \anchor mesher_generate_global_numbering_algorithm_main_nodes_boundary <b><em>Main nodes at the boundary</em></b>
    #
    # Go along the elements at the boundary.
    #
    # Go counter-clockwise along the nodes of the elements.
    #
    # Check if node is at the boundary (this is done, by checking if nodes.bc_markers
    # of the node is > 0, see nodes2d.bc_markers).
    #
    # Check if node was already numbered, if it was, use the numbering used
    # before, that is temporarily stored in node_global_numbering. If not,
    # continue.
    #
    # Check if the node does not have prescribed boundary conditions. if true
    # number the node, using the numbering that starts in 0. If not, continue.
    #
    # Number the node, using the prescribed boundary conditions numbering.
    # Compute the coordinates of the node and add them to the list of the
    # coordinates of boundary nodes with prescribed boundary conditions.
    #
    # Leading to:
    #
    # \image html global_numbering_1.jpg
    # \image latex global_numbering_1-crop.eps "Numbering of boundary main nodes" width=10cm
    #
    ## \anchor mesher_generate_global_numbering_algorithm_nodes_boundary_edges <b><em>Nodes at boundary edges</em></b>
    #
    # Go along the boundary elements.
    #
    # For each boundary element, go along the edges, counter-clockwise,
    # starting at edge number 3, that is the one that connects node 0 to node 1.
    #
    # Check if edge has prescribed boundary conditions (check for 
    # edges.bc_markers > 1, see edges2d.bc_markers). If true, compute the
    # coordinates of the nodes that lie on that edge and add them to the list
    # of nodes with prescribed boundary conditions, number them. If not,
    # continue.
    #
    # Number the nodes according to the non-prescribed boundary conditions
    # numbering.
    #
    # Leading to:
    #
    # \image html global_numbering_1_5.jpg
    # \image latex global_numbering_1_5-crop.eps "Numbering of boundary nodes at edges" width=10cm
    #
    ## \anchor mesher_generate_global_numbering_algorithm_inner_main_nodes <b><em>Inner main nodes</em></b>
    #
    # Go along all elements.
    #
    # Go along the nodes counter-clockwise.
    #
    # Check if the node is not at the boundary, if true, continue.
    #
    # Check if the node was already numbered, if true, give it its global
    # number, if not, continue.
    #
    # Give the node a new number, set the node as numbered and store its global
    # number in node_global_numbering, for temporary use.
    #
    # Leading to:
    #
    # \image html global_numbering_2.jpg
    # \image latex global_numbering_2-crop.eps "Numbering of inner main nodes" width=10cm
    #
    ## \anchor mesher_generate_global_numbering_algorithm_inner_nodes_edges <b><em>Inner nodes at the edges</em></b>
    #
    # Go along all the elements.
    #
    # Check if each of the neighbors of the element exists. If neighbor
    # exists, continue, if not it means it the corresponding edge is at the
    # boundary, hence those nodes were already numbered in previuous passes.
    #
    # Check if neighbor was already numbered. If true, use the numbering of
    # the edge of that neighbor element for the edge of the current element.
    # If not, continue.
    #
    # Number the nodes.
    #
    # After numbering all edges of the element, check the element as numbered.
    #
    # Leading to:
    #
    # \image html global_numbering_3.jpg
    # \image latex global_numbering_3-crop.eps "Numbering of inner nodes at the edges" width=10cm
    #
    ## \anchor mesher_generate_global_numbering_algorithm_inner_element_nodes <b><em>Inner element nodes</em></b>
    #
    # Go along all elements
    #
    # Number the nodes.
    #
    # Leading to:
    #
    # \image html global_numbering_4.jpg
    # \image latex global_numbering_4-crop.eps "Numbering of inner nodes" width=10cm
    #
    ## \anchor mesher_generate_global_numbering_definitions \b Definitions
    #
    ## \anchor mesher_generate_global_numbering_parameters
    ## \param[in] p
    # Defines the polynomial order inside the elements, the same for all
    # elements.\n
    # \b TYPE: \c int32
    #
    ## \param[out] self.global_numbering_zero_form
    # See #global_numbering_zero_form.
    #
    ## \param[out] self.boundary_nodes_indices
    # See #boundary_nodes_indices.
    #
    ## \param[out] self.boundary_nodes_coordinates
    # See #boundary_nodes_coordinates.
    #
    # \todo Still need to update the documentation on the numbering. Added
    # numbering of 1-forms and 2-forms degrees of freedom. The code is temporary
    # because what should be numbered are tangent degrees of freedom and now
    # only dx and dy degrees of freedom, which in the particular case of rectangular
    # meshes is the same.
    #
    def generate_global_numbering(self, p, zero_form=True,\
                                           one_form=True,\
                                           two_form=True):
        # function that, given a mesh object, return an array with the global
        # numbering of all the degrees of freedom of the computational domain, 
        # for the three types of variables that appear in the problem:
        #    0-forms
        #    1-forms (du part and dv part, in the local coordinates)
        #    2-forms
        # For the 0-forms they are the main nodes of the mesh,
        # which are the ones present in the mesh, and the spectral degrees of freedom
        # which are the ones that lie in the interior and at the edges of the
        # elements of the mesh and are originated from the spectral element method,
        # Lobatto locates in both u and v directions.
        #
        # As for the 1-forms they are Gauss located degrees of freedom in the
        # u direction and Lobatto located degrees of freedom in the v direction,
        # for the du component; for the dv component they are Lobatto located
        # degrees of freedom in the u direction and Gauss located degrees of
        # freedom for the v direction.
        #
        # As for the 2-forms they are Gauss located degrees of freedom in both
        # u and v directions.
        #
        # For the 0-forms the numbering occurs in the following way:
        # Simply put, the numbering first numbers the main degrees of freedom at the boundary
        # then the spectral nodes at the boundary. These nodes global numbers
        # will be the highest ones. After these, the main nodes are numbered
        # and then the spectral nodes at the edges which are not at the boundary
        # of the computational domain. Finally the inner nodes are numbered.
        
        
        
        # make a local copy to the number of main nodes
        N_nodes = self.N_nodes
        
        # make a local copy to the number of elements
        N_elements = self.N_elements
        
        # make a local copy to the number of edges
        N_edges = self.N_edges
        
        # compute the total number of nodes, assuming that all elements have the
        # same spectral order p (order of the polynomial approximation)
        N_total_nodes = N_edges*(p-1) + N_nodes + ((p-1)**2)*N_elements        
        self.N_total_nodes = N_total_nodes
        
        # compute the number of nodes with prescribed boundary conditions
        N_prescribed_boundary_nodes = self.nodes.bc_markers[\
                                      self.nodes.bc_markers>1].size
        
        # compute the number of boundary edges
        N_boundary_edges = self.edges.on_boundary.size
        N_prescribed_boundary_edges = self.edges.on_prescribed_boundary.size
        
        # make a reference to the nodes bc_markers array, just not to overload
        # the notation of the coding, since it is just a reference (like a
        # pointer in c) there is no memory overload,
        # that states if the node is at a boundary or not
        # 0 means that the node is at a boundary, value different from 0 means
        # that the value is at the boundary with the boundary condition
        # identified by the number
        nodes_bc_markers = self.nodes.bc_markers

        # create a vector that specifies if a given node was already numbered or
        # not, it is a boolean vector
        node_was_checked = numpy.zeros(N_nodes, dtype="bool")
        
        # create a vector that holds the global numbering of the main nodes
        node_global_numbering = numpy.zeros(N_total_nodes, dtype=settings.index_INT_TYPE)
        
        # set the node global numbering to -1, so that if problem occurs it is
        # easier to check it
        node_global_numbering -= 1
        
        
        # initialize current_number, which stores the current number to be given
        # to the new node found and the current_boundary_number whuch stores
        # the current number to be given to the new node found with prescribed
        # boundary conditions
        current_number = 0
        current_boundary_number = N_total_nodes - \
                                  N_prescribed_boundary_edges*(p-1) -\
                                  N_prescribed_boundary_nodes
        bc_counter = 0
                                
        # compute the total number of boundary nodes
        N_total_boundary_nodes = N_prescribed_boundary_edges*(p-1) + N_prescribed_boundary_nodes
        
        #----------------------------------------------------------------------
        # 0-form
        #----------------------------------------------------------------------
        
        # allocate memory space for the global numbering array, global_numbering_zero_form
        # it is a 3d array with dimensions N_elements x p x p, where p is the
        # order of the spectral element method
        self.global_numbering_zero_form = numpy.zeros([N_elements,p+1,p+1], dtype=settings.index_INT_TYPE)
        
        # set the global numbering to -1, so that if a problem occurs it is easier
        # to check it
        self.global_numbering_zero_form = self.global_numbering_zero_form - 1
        
        # define arrays that contain the parameters to be used when re-using inner
        # edge nodes that were already numbered
        # quadrilateral elements
        element_start_quad = numpy.array([[p,1],[p-1,p],[0,p-1],[1,0]])
        element_end_quad = numpy.array([[p+1,p],[0,p+1],[1,0],[p,1]])
        element_step_quad = numpy.array([[1,1],[-1,1],[1,-1],[1,1]])
        neighbor_start_quad = numpy.array([[p,p-1],[1,p],[0,1],[p-1,0]])
        neighbor_end_quad = numpy.array([[p+1,0],[p,p+1],[1,p],[0,1]])
        neighbor_step_quad = numpy.array([[1,-1],[1,1],[1,1],[-1,1]])
        # triangular elements
        element_start_tri = numpy.array([[p,1],[p,p-1],[0,1]])
        element_end_tri = numpy.array([[p+1,p],[p+1,0],[1,p]])
        element_step_tri = numpy.array([[1,1],[1,-1],[1,1]])
        neighbor_start_tri = numpy.array([[p,p-1],[0,1],[p-1,0]])
        neighbor_end_tri = numpy.array([[p+1,0],[1,p],[0,1]])
        neighbor_step_tri = numpy.array([[1,-1],[1,1],[-1,1]])
        
        # allocate memory space for array that will store the node
        # indices (element, and i and j and boundary condition marker)
        # of the nodes with prescribed boundary conditions
        self.boundary_nodes_indices = numpy.zeros([N_total_boundary_nodes,4],dtype=settings.index_INT_TYPE)
        
        # allocate memory space for the nodes coordinates
        self.boundary_nodes_coordinates =  numpy.zeros([N_total_boundary_nodes,2],dtype="float64")
        
        # compute the collocation points to use in the computation of the
        # boundary nodes coordinates
        collocation_points = gauss_quadrature.collocation_points(p)
        
        # number the main nodes at the boundary
        # go along all elements
        for element in numpy.arange(0, N_elements):
            # store the nodes of the element in a local variable
            element_nodes = self.elements.nodes[element] 
            
            # check if node 0 is at the boundary
            if nodes_bc_markers[element_nodes[0]] != 0:
                # if at the boundary check if it was numbered
                # already
                if node_was_checked[element_nodes[0]]:
                    # give the node its global number
                    self.global_numbering_zero_form[element,0,0] = node_global_numbering[\
                                                        element_nodes[0]]
                else:
                    if nodes_bc_markers[element_nodes[0]] == 1:
                        # node is just at the boundary
                        # number the node
                        self.global_numbering_zero_form[element,0,0] = current_number
                        # check the node as already numbered
                        node_was_checked[element_nodes[0]] = True
                        # update the global numbering to the node_global_numbering
                        # array
                        node_global_numbering[element_nodes[0]] = current_number
                        # update the current_number so that next node will have the 
                        # numbering correct
                        current_number += 1
                    else:
                        # store the node indices and bc_marker
                        self.boundary_nodes_indices[bc_counter] = numpy.array([element, 0, 0, self.nodes.bc_markers[element_nodes[0]]])
                        # store node coordinates
                        self.boundary_nodes_coordinates[bc_counter] = self.nodes.coordinates[element_nodes[0]]
                        bc_counter += 1
                        # the node has prescribed boundary conditions
                        # number the node
                        self.global_numbering_zero_form[element,0,0] = current_boundary_number
                        # check the node as already numbered
                        node_was_checked[element_nodes[0]] = True
                        # update the global numbering to the node_global_numbering
                        # array
                        node_global_numbering[element_nodes[0]] = current_boundary_number
                        # update the current_number so that next node will have the 
                        # numbering correct
                        current_boundary_number += 1
        
            # check if node 1 is at the boundary
            if nodes_bc_markers[element_nodes[1]] != 0:
                # if at the boundary check if it was numbered
                # already
                if node_was_checked[element_nodes[1]]:
                    # give the node its global number
                    self.global_numbering_zero_form[element,p,0] = node_global_numbering[\
                                                        element_nodes[1]]
                else:
                    if nodes_bc_markers[element_nodes[1]] == 1:
                        # node is just at the boundary
                        # number the node
                        self.global_numbering_zero_form[element,p,0] = current_number
                        # check the node as already numbered
                        node_was_checked[element_nodes[1]] = True
                        # update the global numbering to the node_global_numbering
                        # array
                        node_global_numbering[element_nodes[1]] = current_number
                        # update the current_number so that next node will have the 
                        # numbering correct
                        current_number += 1
                    else:
                        # store the node indices and bc_marker
                        self.boundary_nodes_indices[bc_counter] = numpy.array([element, p, 0, self.nodes.bc_markers[element_nodes[1]]])
                        # store node coordinates
                        self.boundary_nodes_coordinates[bc_counter] = self.nodes.coordinates[element_nodes[1]]
                        bc_counter += 1
                        # the node has prescribed boundary conditions
                        # number the node
                        self.global_numbering_zero_form[element,p,0] = current_boundary_number
                        # check the node as already numbered
                        node_was_checked[element_nodes[1]] = True
                        # update the global numbering to the node_global_numbering
                        # array
                        node_global_numbering[element_nodes[1]] = current_boundary_number
                        # update the current_number so that next node will have the 
                        # numbering correct
                        current_boundary_number += 1
            
            # check if node 2 is at the boundary
            if nodes_bc_markers[element_nodes[2]] != 0:
                # if at the boundary check if it was numbered
                # already
                if node_was_checked[element_nodes[2]]:
                    # give the node its global number
                    self.global_numbering_zero_form[element,p,p] = node_global_numbering[\
                                                        element_nodes[2]]
                else:
                    if nodes_bc_markers[element_nodes[2]] == 1:
                        # node is just at the boundary
                        # number the node
                        self.global_numbering_zero_form[element,p,p] = current_number
                        # check the node as already numbered
                        node_was_checked[element_nodes[2]] = True
                        # update the global numbering to the node_global_numbering
                        # array
                        node_global_numbering[element_nodes[2]] = current_number
                        # update the current_number so that next node will have the 
                        # numbering correct
                        current_number += 1
                    else:
                        # store the node indices and bc_marker
                        self.boundary_nodes_indices[bc_counter] = numpy.array([element, p, p, self.nodes.bc_markers[element_nodes[2]]])
                        # store node coordinates
                        self.boundary_nodes_coordinates[bc_counter] = self.nodes.coordinates[element_nodes[2]]
                        bc_counter += 1
                        # the node has prescribed boundary conditions
                        # number the node
                        self.global_numbering_zero_form[element,p,p] = current_boundary_number
                        # check the node as already numbered
                        node_was_checked[element_nodes[2]] = True
                        # update the global numbering to the node_global_numbering
                        # array
                        node_global_numbering[element_nodes[2]] = current_boundary_number
                        # update the current_number so that next node will have the 
                        # numbering correct
                        current_boundary_number += 1
            
            # check if node 3 is at the boundary
            if nodes_bc_markers[element_nodes[3]] != 0:
                # if at the boundary check if it was numbered
                # already
                if node_was_checked[element_nodes[3]]:
                    # give the node its global number
                    self.global_numbering_zero_form[element,0,p] = node_global_numbering[\
                                                        element_nodes[3]]
                else:
                    if nodes_bc_markers[element_nodes[3]] == 1:
                        # node is just at the boundary
                        # number the node
                        self.global_numbering_zero_form[element,0,p] = current_number
                        # check the node as already numbered
                        node_was_checked[element_nodes[3]] = True
                        # update the global numbering to the node_global_numbering
                        # array
                        node_global_numbering[element_nodes[3]] = current_number
                        # update the current_number so that next node will have the 
                        # numbering correct
                        current_number += 1
                    else:
                        # store the node indices and bc_marker
                        self.boundary_nodes_indices[bc_counter] = numpy.array([element, 0, p, self.nodes.bc_markers[element_nodes[3]]])
                        # store node coordinates
                        self.boundary_nodes_coordinates[bc_counter] = self.nodes.coordinates[element_nodes[3]]
                        bc_counter += 1
                        # the node has prescribed boundary conditions
                        # number the node
                        self.global_numbering_zero_form[element,0,p] = current_boundary_number
                        # check the node as already numbered
                        node_was_checked[element_nodes[3]] = True
                        # update the global numbering to the node_global_numbering
                        # array
                        node_global_numbering[element_nodes[3]] = current_boundary_number
                        # update the current_number so that next node will have the 
                        # numbering correct
                        current_boundary_number += 1
        
        # number the nodes at boundary edges
        # go along the elements at the boundary
        for boundary_element_index, element in enumerate(self.elements.on_boundary):
            # store the nodes of the element in a local variable
            element_nodes = self.elements.nodes[element]
            
            # generate the transfinite mapping fo the element
            myTransfinite = math_tools.TransfiniteMapping_2d(self.nodes.coordinates[element_nodes])
            
            # store the element boundary edges in a local variable
            # recall that -1 means that the edge is not at the boundary
            element_edges = self.elements.boundary_edges[boundary_element_index]
            
            # check if edge 3 is at the boundary
            if element_edges[3] != -1:
                # edge 3 is at the boundary
                # check if it has prescribed boundary conditions
                if self.edges.bc_markers[element_edges[3]] > 1:
                    for node in range(1,p):
                        # store the node indices and bc_markers
                        self.boundary_nodes_indices[bc_counter] = numpy.array([element, node, 0, self.edges.bc_markers[element_edges[3]]])
                        # compute node coordinates
                        x = myTransfinite.X(collocation_points[node],-1.0)
                        y = myTransfinite.Y(collocation_points[node],-1.0)
                        self.boundary_nodes_coordinates[bc_counter] = numpy.array([x,y])
                        bc_counter += 1
                    # number the nodes, there are p-1 nodes in each edge
                    # these are unique, since boundary edges are not shared
                    self.global_numbering_zero_form[element,1:p,0] = numpy.arange(\
                                                        current_boundary_number,\
                                                        current_boundary_number + p - 1,\
                                                        dtype="int32")
                    # update current_boundary_number
                    current_boundary_number = current_boundary_number + p - 1
                else:
                    # number the nodes according to the non prescribed boundary
                    # conditions numbering
                    self.global_numbering_zero_form[element,1:p,0] = numpy.arange(\
                                                        current_number,\
                                                        current_number + p - 1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p - 1
            
            # check if edge 0 is at the boundary
            if element_edges[0] != -1:
                # edge 0 is at the boundary
                # check if it has prescribed boundary conditions
                if self.edges.bc_markers[element_edges[0]] > 1:
                    for node in range(1,p):
                        # store the node indices and bc_markers
                        self.boundary_nodes_indices[bc_counter] = numpy.array([element, p, node, self.edges.bc_markers[element_edges[0]]])
                        # compute node coordinates
                        x = myTransfinite.X(1.0, collocation_points[node])
                        y = myTransfinite.Y(1.0, collocation_points[node])
                        self.boundary_nodes_coordinates[bc_counter] = numpy.array([x,y])
                        bc_counter += 1
                    # number the nodes, there are p-1 nodes in each edge
                    # these are unique, since boundary edges are not shared
                    self.global_numbering_zero_form[element,p,1:p] = numpy.arange(\
                                                        current_boundary_number,\
                                                        current_boundary_number + p - 1,\
                                                        dtype="int32")
                    # update current_boundary_number
                    current_boundary_number = current_boundary_number + p - 1
                else:
                    # number the nodes according to the non prescribed boundary
                    # conditions numbering
                    self.global_numbering_zero_form[element,p,1:p] = numpy.arange(\
                                                        current_number,\
                                                        current_number + p - 1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p - 1
            
            # check if edge 1 is at the boundary
            if element_edges[1] != -1:
                # edge 1 is at the boundary
                # check if it has prescribed boundary conditions
                if self.edges.bc_markers[element_edges[1]] > 1:
                    for node in range(1,p):
                        # store the node indices and bc_markers
                        self.boundary_nodes_indices[bc_counter] = numpy.array([element, node, p, self.edges.bc_markers[element_edges[1]]])
                        # compute node coordinates
                        x = myTransfinite.X(collocation_points[node], 1.0)
                        y = myTransfinite.Y(collocation_points[node], 1.0)
                        self.boundary_nodes_coordinates[bc_counter] = numpy.array([x,y])
                        bc_counter += 1
                    # number the nodes, there are p-1 nodes in each edge
                    # these are unique, since boundary edges are not shared
                    self.global_numbering_zero_form[element,1:p,p] = numpy.arange(\
                                                        current_boundary_number + p - 2,\
                                                        current_boundary_number - 1,\
                                                        -1,\
                                                        dtype="int32")
                    # update current_boundary_number
                    current_boundary_number = current_boundary_number + p - 1
                else:
                    # number the nodes according to the non prescribed boundary
                    # conditions numbering
                    self.global_numbering_zero_form[element,1:p,p] = numpy.arange(\
                                                        current_number + p - 2,\
                                                        current_number - 1,\
                                                        -1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p - 1
                    
            # check if edge 2 is at the boundary
            if element_edges[2] != -1:
                # edge 2 is at the boundary
                # check if it has prescribed boundary conditions
                if self.edges.bc_markers[element_edges[2]] > 1:
                    for node in range(1,p):
                        # store the node indices and bc_markers
                        self.boundary_nodes_indices[bc_counter] = numpy.array([element, 0, node, self.edges.bc_markers[element_edges[2]]])
                        # compute node coordinates
                        x = myTransfinite.X(-1, collocation_points[node])
                        y = myTransfinite.Y(-1, collocation_points[node])
                        self.boundary_nodes_coordinates[bc_counter] = numpy.array([x,y])
                        bc_counter += 1
                    # number the nodes, there are p-1 nodes in each edge
                    # these are unique, since boundary edges are not shared
                    self.global_numbering_zero_form[element,0,1:p] = numpy.arange(\
                                                        current_boundary_number + p - 2,\
                                                        current_boundary_number - 1,\
                                                        -1,\
                                                        dtype="int32")
                    # update current_boundary_number
                    current_boundary_number = current_boundary_number + p - 1
                else:
                    # number the nodes according to the non prescribed boundary
                    # conditions numbering
                    self.global_numbering_zero_form[element,0,1:p] = numpy.arange(\
                                                        current_number + p - 2,\
                                                        current_number - 1,\
                                                        -1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p - 1
        
        # number the inner main nodes
        # go along all the elements
        for element in numpy.arange(0, N_elements):
            # store the nodes of the element in a local variable
            element_nodes = self.elements.nodes[element]
            
            # check if node 0 is not at the boundary
            if nodes_bc_markers[element_nodes[0]] == 0:
                # if not at the boundary check if it was numbered
                # already
                if node_was_checked[element_nodes[0]]:
                    # give the node its global number
                    self.global_numbering_zero_form[element,0,0] = node_global_numbering[\
                                                        element_nodes[0]]
                else:
                    # number the node
                    self.global_numbering_zero_form[element,0,0] = current_number
                    # check the node as already numbered
                    node_was_checked[element_nodes[0]] = True
                    # update the global numbering to the node_global_numbering
                    # array
                    node_global_numbering[element_nodes[0]] = current_number
                    # update the current_number so that next node will have the 
                    # numbering correct
                    current_number += 1
                    
            # check if node 1 is not at the boundary
            if nodes_bc_markers[element_nodes[1]] == 0:
                # if not at the boundary check if it was numbered
                # already
                if node_was_checked[element_nodes[1]]:
                    # give the node its global number
                    self.global_numbering_zero_form[element,p,0] = node_global_numbering[\
                                                        element_nodes[1]]
                else:
                    # number the node
                    self.global_numbering_zero_form[element,p,0] = current_number
                    # check the node as already numbered
                    node_was_checked[element_nodes[1]] = True
                    # update the global numbering to the node_global_numbering
                    # array
                    node_global_numbering[element_nodes[1]] = current_number
                    # update the current_number so that next node will have the 
                    # numbering correct
                    current_number += 1
                    
            # check if node 2 is not at the boundary
            if nodes_bc_markers[element_nodes[2]] == 0:
                # if not at the boundary check if it was numbered
                # already
                if node_was_checked[element_nodes[2]]:
                    # give the node its global number
                    self.global_numbering_zero_form[element,p,p] = node_global_numbering[\
                                                        element_nodes[2]]
                else:
                    # number the node
                    self.global_numbering_zero_form[element,p,p] = current_number
                    # check the node as already numbered
                    node_was_checked[element_nodes[2]] = True
                    # update the global numbering to the node_global_numbering
                    # array
                    node_global_numbering[element_nodes[2]] = current_number
                    # update the current_number so that next node will have the 
                    # numbering correct
                    current_number += 1
                    
            # check if node 3 is not at the boundary
            if nodes_bc_markers[element_nodes[3]] == 0:
                # if not at the boundary check if it was numbered
                # already
                if node_was_checked[element_nodes[3]]:
                    # give the node its global number
                    self.global_numbering_zero_form[element,0,p] = node_global_numbering[\
                                                        element_nodes[3]]
                else:
                    # number the node
                    self.global_numbering_zero_form[element,0,p] = current_number
                    # check the node as already numbered
                    node_was_checked[element_nodes[3]] = True
                    # update the global numbering to the node_global_numbering
                    # array
                    node_global_numbering[element_nodes[3]] = current_number
                    # update the current_number so that next node will have the 
                    # numbering correct
                    current_number += 1
                    
        # number the inner nodes at the edges
        # create a vector that specifies if a given element was already numbered or
        # not, it is a boolean vector
        element_was_checked = numpy.zeros(N_elements, dtype="bool")
        
        # generate auxiliary neighbor numbering 0-1
        aux_neighbor = numpy.arange(0,4)
        
        # go along all the elements
        for element in numpy.arange(0, N_elements):
            # store the nodes of the element in a local variable
            element_nodes = self.elements.nodes[element]
            
            # store the element neighbors in a local variable
            # recall that -1 means that the element does not have that neighbor
            element_neighbors = self.elements.neighbors[element]
            
            # check if neighbor exists
            if element_neighbors[3] != -1:
                # check if neighbor 3 was already numbered
                if element_was_checked[element_neighbors[3]] == True:
                    # find which neighbor of the neighbor the element is
                    element_is_neighbor = aux_neighbor[self.elements.neighbors\
                                                [element_neighbors[3]] == element]
                                                
                    # give the nodes their global numbering
                    # generate i nodes indices of the element
                    element_i_indices = numpy.arange(element_start_quad[3,0],\
                                                     element_end_quad[3,0],\
                                                     element_step_quad[3,0])
                    
                    # generate j nodes indices of the element
                    element_j_indices = numpy.arange(element_start_quad[3,1],\
                                                     element_end_quad[3,1],\
                                                     element_step_quad[3,1])
                    
                    # generate i nodes indices of the neighbor
                    neighbor_i_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,0],\
                                                       neighbor_end_quad[element_is_neighbor,0],\
                                                       neighbor_step_quad[element_is_neighbor,0])
                    
                    # generate j nodes indices of the neighbor
                    neighbor_j_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,1],\
                                                       neighbor_end_quad[element_is_neighbor,1],\
                                                       neighbor_step_quad[element_is_neighbor,1])
                    
                    # actually give the numbering
                    self.global_numbering_zero_form[element,element_i_indices,element_j_indices]\
                                    = self.global_numbering_zero_form[element_neighbors[3],neighbor_i_indices,neighbor_j_indices]
                else:
                    # number the nodes
                    self.global_numbering_zero_form[element,1:p,0] = numpy.arange(\
                                                        current_number,\
                                                        current_number + p - 1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p - 1
            
            # check if neighbor exists
            if element_neighbors[0] != -1:
                # check if neighbor 0 was already numbered
                if element_was_checked[element_neighbors[0]] == True:
                    # find which neighbor of the neighbor the element is
                    element_is_neighbor = aux_neighbor[self.elements.neighbors\
                                                [element_neighbors[0]] == element]
                    # give the nodes their global numbering
                    # generate i nodes indices of the element
                    element_i_indices = numpy.arange(element_start_quad[0,0],\
                                                     element_end_quad[0,0],\
                                                     element_step_quad[0,0])
                    
                    # generate j nodes indices of the element
                    element_j_indices = numpy.arange(element_start_quad[0,1],\
                                                     element_end_quad[0,1],\
                                                     element_step_quad[0,1])
                    
                    # generate i nodes indices of the neighbor
                    neighbor_i_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,0],\
                                                       neighbor_end_quad[element_is_neighbor,0],\
                                                       neighbor_step_quad[element_is_neighbor,0])
                    
                    # generate j nodes indices of the neighbor
                    neighbor_j_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,1],\
                                                       neighbor_end_quad[element_is_neighbor,1],\
                                                       neighbor_step_quad[element_is_neighbor,1])
                    
                    # actually give the numbering
                    self.global_numbering_zero_form[element,element_i_indices,element_j_indices]\
                                    = self.global_numbering_zero_form[element_neighbors[0],neighbor_i_indices,neighbor_j_indices]
                else:
                    # number the nodes
                    self.global_numbering_zero_form[element,p,1:p] = numpy.arange(\
                                                        current_number,\
                                                        current_number + p - 1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p - 1
            
            # check if neighbor exists
            if element_neighbors[1] != -1:
                # check if neighbor 1 was already numbered
                if element_was_checked[element_neighbors[1]] == True:
                    # find which neighbor of the neighbor the element is
                    element_is_neighbor = aux_neighbor[self.elements.neighbors\
                                                [element_neighbors[1]] == element]
                    # give the nodes their global numbering
                    # generate i nodes indices of the element
                    element_i_indices = numpy.arange(element_start_quad[1,0],\
                                                     element_end_quad[1,0],\
                                                     element_step_quad[1,0])
                    
                    # generate j nodes indices of the element
                    element_j_indices = numpy.arange(element_start_quad[1,1],\
                                                     element_end_quad[1,1],\
                                                     element_step_quad[1,1])
                    
                    # generate i nodes indices of the neighbor
                    neighbor_i_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,0],\
                                                       neighbor_end_quad[element_is_neighbor,0],\
                                                       neighbor_step_quad[element_is_neighbor,0])
                    
                    # generate j nodes indices of the neighbor
                    neighbor_j_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,1],\
                                                       neighbor_end_quad[element_is_neighbor,1],\
                                                       neighbor_step_quad[element_is_neighbor,1])
                    
                    # actually give the numbering
                    self.global_numbering_zero_form[element,element_i_indices,element_j_indices]\
                                    = self.global_numbering_zero_form[element_neighbors[1],neighbor_i_indices,neighbor_j_indices]
                else:
                    # number the nodes
                    self.global_numbering_zero_form[element,1:p,p] = numpy.arange(\
                                                        current_number + p - 2,\
                                                        current_number - 1,\
                                                        -1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p - 1
                    
            # check if neighbor exists
            if element_neighbors[2] != -1:   
                # check if neighbor 2 was already numbered
                if element_was_checked[element_neighbors[2]] == True:
                    # find which neighbor of the neighbor the element is
                    element_is_neighbor = aux_neighbor[self.elements.neighbors\
                                                [element_neighbors[2]] == element]
                    # give the nodes their global numbering
                    # generate i nodes indices of the element
                    element_i_indices = numpy.arange(element_start_quad[2,0],\
                                                     element_end_quad[2,0],\
                                                     element_step_quad[2,0])
                    
                    # generate j nodes indices of the element
                    element_j_indices = numpy.arange(element_start_quad[2,1],\
                                                     element_end_quad[2,1],\
                                                     element_step_quad[2,1])
                    
                    # generate i nodes indices of the neighbor
                    neighbor_i_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,0],\
                                                       neighbor_end_quad[element_is_neighbor,0],\
                                                       neighbor_step_quad[element_is_neighbor,0])
                    
                    # generate j nodes indices of the neighbor
                    neighbor_j_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,1],\
                                                       neighbor_end_quad[element_is_neighbor,1],\
                                                       neighbor_step_quad[element_is_neighbor,1])
                    
                    # actually give the numbering
                    self.global_numbering_zero_form[element,element_i_indices,element_j_indices]\
                                    = self.global_numbering_zero_form[element_neighbors[2],neighbor_i_indices,neighbor_j_indices]
                else:
                    # number the nodes
                    self.global_numbering_zero_form[element,0,1:p] = numpy.arange(\
                                                        current_number + p - 2,\
                                                        current_number - 1,\
                                                        -1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p - 1
              
            # set the element as checked
            element_was_checked[element] = True
        
        # number inner element nodes
        # generate numbering array
        inner_nodes_numbers = numpy.arange(current_number,
                                                current_number + (p-1)**2).reshape([p-1,p-1], order="fortran")
        # go along all the elements
        for element in numpy.arange(0, N_elements):
            # number the nodes
            self.global_numbering_zero_form[element,1:p,1:p] = inner_nodes_numbers
            
            # update the inner nodes numbers
            inner_nodes_numbers += (p-1)**2
            
        #----------------------------------------------------------------------
        # 1-form
        #----------------------------------------------------------------------
        # \todo No boundary conditions are implemented for 1-forms, degrees of
        # freedom start numbering from edges, from element 0 to self.N_elements-1.
        # Implementation of boundary conditions is needed for Neumann boundary
        # conditions. Then the numbering shall follow the same procedure as for
        # 0-forms: boundary degrees of freedom are numbered last.
        
        p = p+1
        
        # du part
        
        # allocate memory space for the global numbering array, global_numbering_one_form_du
        # it is a 3d array with dimensions N_elements x p x p+1, where p is the
        # order of the spectral element method
        self.global_numbering_one_form_du = numpy.zeros([N_elements,p,p+1], dtype=settings.index_INT_TYPE)
        
        # set the global numbering to -1, so that if a problem occurs it is easier
        # to check it
        self.global_numbering_one_form_du = self.global_numbering_one_form_du - 1
        
        # number the degrees of freedom at the edges
        # create a vector that specifies if a given element was already numbered or
        # not, it is a boolean vector
        element_was_checked = numpy.zeros(N_elements, dtype="bool")
        
        # generate auxiliary neighbor numbering 0-1
        aux_neighbor = numpy.arange(0,4)
        
        # initialize current number of degrees of freedom
        current_number = 0
        
        # define arrays that contain the parameters to be used when re-using inner
        # edge degrees of freedom that were already numbered
        # quadrilateral elements
        element_start_quad = numpy.array([[p,1],[p-1,p],[0,p-1],[0,0]])
        element_end_quad = numpy.array([[p+1,p],[-1,p+1],[1,0],[p,1]])
        element_step_quad = numpy.array([[1,1],[-1,1],[1,-1],[1,1]])
        neighbor_start_quad = numpy.array([[p,p-1],[0,p],[0,1],[p-1,0]])
        neighbor_end_quad = numpy.array([[p+1,0],[p,p+1],[1,p],[-1,1]])
        neighbor_step_quad = numpy.array([[1,-1],[1,1],[1,1],[-1,1]])
        # triangular elements
##        element_start_tri = numpy.array([[p,1],[p,p-1],[0,1]])
##        element_end_tri = numpy.array([[p+1,p],[p+1,0],[1,p]])
##        element_step_tri = numpy.array([[1,1],[1,-1],[1,1]])
##        neighbor_start_tri = numpy.array([[p,p-1],[0,1],[p-1,0]])
##        neighbor_end_tri = numpy.array([[p+1,0],[1,p],[0,1]])
##        neighbor_step_tri = numpy.array([[1,-1],[1,1],[-1,1]])
        
        # go along all the elements
        for element in numpy.arange(0, N_elements):
            # store the nodes of the element in a local variable
            element_nodes = self.elements.nodes[element]
            
            # store the element neighbors in a local variable
            # recall that -1 means that the element does not have that neighbor
            element_neighbors = self.elements.neighbors[element]
            
            # since 1-forms du part only has degrees of freedom at edges 3 and 1
            # one only needs to check for neighbors 1 and 3 of the element and
            # number only degrees of freedom at this edges.
            
            # check if neighbor exists
            if element_neighbors[3] != -1:
                # check if neighbor 3 was already numbered
                if element_was_checked[element_neighbors[3]] == True:
                    # find which neighbor of the neighbor the element is
                    element_is_neighbor = aux_neighbor[self.elements.neighbors\
                                                [element_neighbors[3]] == element]
                                                
                    # give the nodes their global numbering
                    # generate i degrees of freedom indices of the element
                    element_i_indices = numpy.arange(element_start_quad[3,0],\
                                                     element_end_quad[3,0],\
                                                     element_step_quad[3,0])
                    
                    # generate j nodes indices of the element
                    element_j_indices = numpy.arange(element_start_quad[3,1],\
                                                     element_end_quad[3,1],\
                                                     element_step_quad[3,1])
                    
                    # generate i nodes indices of the neighbor
                    neighbor_i_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,0],\
                                                       neighbor_end_quad[element_is_neighbor,0],\
                                                       neighbor_step_quad[element_is_neighbor,0])
                    
                    # generate j nodes indices of the neighbor
                    neighbor_j_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,1],\
                                                       neighbor_end_quad[element_is_neighbor,1],\
                                                       neighbor_step_quad[element_is_neighbor,1])
                    
                    # actually give the numbering
                    self.global_numbering_one_form_du[element,element_i_indices,element_j_indices]\
                                    = self.global_numbering_one_form_du[element_neighbors[3],neighbor_i_indices,neighbor_j_indices]
                                    
                else:
                    # number the degrees of freedom
                    self.global_numbering_one_form_du[element,0:p,0] = numpy.arange(\
                                                        current_number,\
                                                        current_number + p,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p
            
            # if neighbor does not exist simply number the degrees of freedom
            else:
            # number the degrees of freedom
                self.global_numbering_one_form_du[element,0:p,0] = numpy.arange(\
                                                    current_number,\
                                                    current_number + p,\
                                                    dtype="int32")
                # update current_number
                current_number = current_number + p
                    
            # check if neighbor exists
            if element_neighbors[1] != -1:
                # check if neighbor 1 was already numbered
                if element_was_checked[element_neighbors[1]] == True:
                    # find which neighbor of the neighbor the element is
                    element_is_neighbor = aux_neighbor[self.elements.neighbors\
                                                [element_neighbors[1]] == element]
                    # give the nodes their global numbering
                    # generate i nodes indices of the element
                    element_i_indices = numpy.arange(element_start_quad[1,0],\
                                                     element_end_quad[1,0],\
                                                     element_step_quad[1,0])
                    
                    # generate j nodes indices of the element
                    element_j_indices = numpy.arange(element_start_quad[1,1],\
                                                     element_end_quad[1,1],\
                                                     element_step_quad[1,1])
                    
                    # generate i nodes indices of the neighbor
                    neighbor_i_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,0],\
                                                       neighbor_end_quad[element_is_neighbor,0],\
                                                       neighbor_step_quad[element_is_neighbor,0])
                    
                    # generate j nodes indices of the neighbor
                    neighbor_j_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,1],\
                                                       neighbor_end_quad[element_is_neighbor,1],\
                                                       neighbor_step_quad[element_is_neighbor,1])
                    
                    # actually give the numbering
                    self.global_numbering_one_form_du[element,element_i_indices,element_j_indices]\
                                    = self.global_numbering_one_form_du[element_neighbors[1],neighbor_i_indices,neighbor_j_indices]
                else:
                    # number the degrees of freedom
                    self.global_numbering_one_form_du[element,0:p,p] = numpy.arange(\
                                                        current_number + p - 1,\
                                                        current_number - 1,\
                                                        -1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p
            
            # if neighbor does not exist simply number the degrees of freedom
            else:
            # number the degrees of freedom
                self.global_numbering_one_form_du[element,0:p,p] = numpy.arange(\
                                                    current_number + p - 1,\
                                                    current_number - 1,\
                                                    -1,\
                                                    dtype="int32")
                # update current_number
                current_number = current_number + p
                  
            # set the element as checked
            element_was_checked[element] = True
        
        # number inner element degrees of freedom
        # generate numbering array
        if p > 1:
            inner_nodes_numbers = numpy.arange(current_number,\
                                               current_number + p*(p-1)).reshape([p,p-1], order="fortran")
        else:
            inner_nodes_numbers = numpy.arange(current_number,\
                                               current_number + (p-1)).flatten()
        
        # go along all the elements
        for element in numpy.arange(0, N_elements):
            if p > 1:
                # number the nodes
                self.global_numbering_one_form_du[element,0:p,1:p] = inner_nodes_numbers
                # update the inner nodes numbers
                inner_nodes_numbers += (p-1)*p
            else:
                # number the nodes
                self.global_numbering_one_form_du[element,0,1:p] = inner_nodes_numbers
                # update the inner nodes numbers
                inner_nodes_numbers += p
                
        # dv part
        
        # allocate memory space for the global numbering array, global_numbering_one_form_du
        # it is a 3d array with dimensions N_elements x p x p+1, where p is the
        # order of the spectral element method
        self.global_numbering_one_form_dv = numpy.zeros([N_elements,p+1,p], dtype=settings.index_INT_TYPE)
        
        # set the global numbering to -1, so that if a problem occurs it is easier
        # to check it
        self.global_numbering_one_form_dv = self.global_numbering_one_form_dv - 1
        
        # number the degrees of freedom at the edges
        # create a vector that specifies if a given element was already numbered or
        # not, it is a boolean vector
        element_was_checked = numpy.zeros(N_elements, dtype="bool")
        
        # generate auxiliary neighbor numbering 0-1
        aux_neighbor = numpy.arange(0,4)
        
        # initialize current number of degrees of freedom
        current_number = 0
        
        # define arrays that contain the parameters to be used when re-using inner
        # edge degrees of freedom that were already numbered
        # quadrilateral elements
        element_start_quad = numpy.array([[p,0],[p-1,p],[0,p-1],[0,0]])
        element_end_quad = numpy.array([[p+1,p],[-1,p+1],[1,-1],[p,1]])
        element_step_quad = numpy.array([[1,1],[-1,1],[1,-1],[1,1]])
        neighbor_start_quad = numpy.array([[p,p-1],[0,p],[0,0],[p-1,0]])
        neighbor_end_quad = numpy.array([[p+1,-1],[p,p+1],[1,p],[-1,1]])
        neighbor_step_quad = numpy.array([[1,-1],[1,1],[1,1],[-1,1]])
        # triangular elements
##        element_start_tri = numpy.array([[p,1],[p,p-1],[0,1]])
##        element_end_tri = numpy.array([[p+1,p],[p+1,0],[1,p]])
##        element_step_tri = numpy.array([[1,1],[1,-1],[1,1]])
##        neighbor_start_tri = numpy.array([[p,p-1],[0,1],[p-1,0]])
##        neighbor_end_tri = numpy.array([[p+1,0],[1,p],[0,1]])
##        neighbor_step_tri = numpy.array([[1,-1],[1,1],[-1,1]])
        
        # go along all the elements
        for element in numpy.arange(0, N_elements):
            # store the nodes of the element in a local variable
            element_nodes = self.elements.nodes[element]
            
            # store the element neighbors in a local variable
            # recall that -1 means that the element does not have that neighbor
            element_neighbors = self.elements.neighbors[element]
            
            # check if neighbor exists
            if element_neighbors[0] != -1:
                # check if neighbor 0 was already numbered
                if element_was_checked[element_neighbors[0]] == True:
                    # find which neighbor of the neighbor the element is
                    element_is_neighbor = aux_neighbor[self.elements.neighbors\
                                                [element_neighbors[0]] == element]
                    # give the nodes their global numbering
                    # generate i nodes indices of the element
                    element_i_indices = numpy.arange(element_start_quad[0,0],\
                                                     element_end_quad[0,0],\
                                                     element_step_quad[0,0])
                    
                    # generate j nodes indices of the element
                    element_j_indices = numpy.arange(element_start_quad[0,1],\
                                                     element_end_quad[0,1],\
                                                     element_step_quad[0,1])
                    
                    # generate i nodes indices of the neighbor
                    neighbor_i_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,0],\
                                                       neighbor_end_quad[element_is_neighbor,0],\
                                                       neighbor_step_quad[element_is_neighbor,0])
                    
                    # generate j nodes indices of the neighbor
                    neighbor_j_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,1],\
                                                       neighbor_end_quad[element_is_neighbor,1],\
                                                       neighbor_step_quad[element_is_neighbor,1])
                    
                    # actually give the numbering
                    self.global_numbering_one_form_dv[element,element_i_indices,element_j_indices]\
                                    = self.global_numbering_one_form_dv[element_neighbors[0],neighbor_i_indices,neighbor_j_indices]
                else:
                    # number the nodes
                    self.global_numbering_one_form_dv[element,p,0:p] = numpy.arange(\
                                                        current_number,\
                                                        current_number + p,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p
            
            # if neighbor does not exist simply number the degrees of freedom
            else:
            # number the degrees of freedom
                self.global_numbering_one_form_dv[element,p,0:p] = numpy.arange(\
                                                    current_number,\
                                                    current_number + p,\
                                                    dtype="int32")
                # update current_number
                current_number = current_number + p
                
            # check if neighbor exists
            if element_neighbors[2] != -1:   
                # check if neighbor 2 was already numbered
                if element_was_checked[element_neighbors[2]] == True:
                    # find which neighbor of the neighbor the element is
                    element_is_neighbor = aux_neighbor[self.elements.neighbors\
                                                [element_neighbors[2]] == element]
                    # give the nodes their global numbering
                    # generate i nodes indices of the element
                    element_i_indices = numpy.arange(element_start_quad[2,0],\
                                                     element_end_quad[2,0],\
                                                     element_step_quad[2,0])
                    
                    # generate j nodes indices of the element
                    element_j_indices = numpy.arange(element_start_quad[2,1],\
                                                     element_end_quad[2,1],\
                                                     element_step_quad[2,1])
                    
                    # generate i nodes indices of the neighbor
                    neighbor_i_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,0],\
                                                       neighbor_end_quad[element_is_neighbor,0],\
                                                       neighbor_step_quad[element_is_neighbor,0])
                    
                    # generate j nodes indices of the neighbor
                    neighbor_j_indices = numpy.arange(neighbor_start_quad[element_is_neighbor,1],\
                                                       neighbor_end_quad[element_is_neighbor,1],\
                                                       neighbor_step_quad[element_is_neighbor,1])
                    
                    # actually give the numbering
                    self.global_numbering_one_form_dv[element,element_i_indices,element_j_indices]\
                                    = self.global_numbering_one_form_dv[element_neighbors[2],neighbor_i_indices,neighbor_j_indices]
                else:
                    # number the nodes
                    self.global_numbering_one_form_dv[element,0,0:p] = numpy.arange(\
                                                        current_number + p - 1,\
                                                        current_number - 1,\
                                                        -1,\
                                                        dtype="int32")
                    # update current_number
                    current_number = current_number + p
                    
            else:
                # number the nodes
                self.global_numbering_one_form_dv[element,0,0:p] = numpy.arange(\
                                                    current_number + p - 1,\
                                                    current_number - 1,\
                                                    -1,\
                                                    dtype="int32")
                # update current_number
                current_number = current_number + p
              
            # set the element as checked
            element_was_checked[element] = True
        
        # number inner element degrees of freedom
        # generate numbering array
        if p > 1:
            inner_nodes_numbers = numpy.arange(current_number,\
                                               current_number + p*(p-1)).reshape([p-1,p], order="fortran")
        else:
            inner_nodes_numbers = numpy.arange(current_number,\
                                               current_number + (p-1)).flatten()
        
        # go along all the elements
        for element in numpy.arange(0, N_elements):
            if p > 1:
                # number the nodes
                self.global_numbering_one_form_dv[element,1:p,0:p] = inner_nodes_numbers
                # update the inner nodes numbers
                inner_nodes_numbers += (p-1)*p
            else:
                # number the nodes
                self.global_numbering_one_form_dv[element,1:p,0] = inner_nodes_numbers
                # update the inner nodes numbers
                inner_nodes_numbers += p
        
        #----------------------------------------------------------------------
        # 0-form
        #----------------------------------------------------------------------
        
        p = p - 1
        
        # allocate memory space for the global numbering array, global_numbering_one_form_du
        # it is a 3d array with dimensions N_elements x p x p+1, where p is the
        # order of the spectral element method
        self.global_numbering_two_form = numpy.zeros([N_elements,p,p], dtype=settings.index_INT_TYPE)
        
        # set the global numbering to -1, so that if a problem occurs it is easier
        # to check it
        self.global_numbering_two_form = self.global_numbering_two_form - 1
        
        # initialize current number of degrees of freedom
        current_number = 0
        
        # number inner element degrees of freedom
        # generate numbering array
        if p > 1:
            inner_nodes_numbers = numpy.arange(current_number,\
                                               current_number + (p**2)).reshape([p,p], order="fortran")
                                            
            # go along all the elements
            for element in numpy.arange(0, N_elements):
                # number the nodes
                self.global_numbering_two_form[element,0:p,0:p] = inner_nodes_numbers
                # update the inner nodes numbers
                inner_nodes_numbers += p**2
        
        else:
            inner_nodes_numbers = current_number
            
            # go along all the elements
            for element in numpy.arange(0, N_elements):
                # number the nodes
                self.global_numbering_two_form[element,0,0] = inner_nodes_numbers
                # update the inner nodes numbers
                inner_nodes_numbers += 1
            
    ## \brief Computes, for each node, the elements that contain a node.
    #
    ## \anchor mesher_generate_global_numbering_brief <b>Brief description</b>
    #
    # Takes the element object of self and updates the array nodes.where of node
    # object of self. This array contains, for each node, the list of elements
    # where the node is present. See nodes2d.where and nodes2d.update_where for
    # more information.
    #
    ## \anchor mesher_generate_global_numbering_parameters
    #
    ## \param[out] self.nodes.where
    # See #nodes and nodes2d.where.
    #    
    def __find_nodes(self):
        # takes the element object of self and updates in the node object
        # of self the array nodes.where which contains, for each node, the
        # elements that contain the node
        
        # get the maximum number of times a node appears on a different element
        node_max_repetition = numpy.bincount(self.elements.nodes.flatten()).max()
        
        # allocate memory space for a matrix that stores in which elements
        # a node appears, this is a N_nodes x node_max_repetition array
        # notice that the ordering used is fortran, this is just for memory
        # allocation, the user does not need to know this for accessing the
        # array. This has to be done since we will call a fortran function
        # to do the actual demanding computation, so that things go faster.
        # Recall that c stores arrays rows after rows and fortran columns after
        # columns, since python is made in c it stores arrays by default in c
        # fashion, hence the need to explicitly say the ordering here. If this
        # was not done, once the fortran funtion is called a copy of the array
        # is made, which creates a big memory overhead and also slows down the
        # code.
        # The last column holds the information that helps keep track of how
        # many elements were found that contain the current node, until the
        # moment, during the finding process, after the finding process
        # is finished, it hold the number of elements where the node was found.
        # Hence the +1 in the second dimension (node_max_repetition+1).
        # 1 is subtracted since -1 means that no more nodes where found. That
        # is, when going along a row, finding a -1 means that the node was not
        # found in no more elements.
        node_in_elements = numpy.zeros([self.N_nodes, node_max_repetition+1],\
                                       dtype='int32',\
                                       order="fortran")
        # this summing is done like this otherwise a copy of the array is made
        # and the fortran ordering is lost
        node_in_elements -= 1
                                   
        # set the last row to zero
        node_in_elements[:,-1] = 0
        
        # there is no need to convert to fortran ordering the self.element.nodes
        # array since when it was created it was already fortan, this was taken
        # care off
        node_in_elements = _mesher.get_elements(self.elements.nodes,\
                                                node_in_elements)
        
        # update nodes.where array
        self.nodes.update_where(node_in_elements)
        
    ## \brief Computes, for each element at the boundary, the edges they contain.
    #
    ## \anchor mesher_generate_global_numbering_brief <b>Brief description</b>
    #
    # Computes, for each element at the boundary, the edges they contain.
    #
    ## \anchor mesher_generate_global_numbering_parameters
    #
    ## \param[out] self.elements.boundary_edges
    # See #elements and elements2d.boundary_edges.
    #        
    def __generate_boundary_element_edges(self):
        # generate the edges information for the boundary elements
        # only the edges at the boundary

        # compute the number of elements at the boundary
        N_boundary_elements = self.elements.on_boundary.shape[0]
        
        # allocate memory space for boundary_element_edges array
        boundary_element_edges = numpy.zeros([N_boundary_elements, 4],\
                                             dtype="int32",\
                                             order="fortran")
        
        # set all values to -1, so that when -1 is found, it means that
        # that edge of the element is not at the boundary
        boundary_element_edges -= 1
                                    
        # create an array with the nodes of the boundary elements
        boundary_elements = numpy.array(self.elements.nodes[\
                            self.elements.on_boundary],\
                            copy=1,\
                            order="fortran")
        
        # create a copy with the nodes of the boundary eges
        boundary_edges = numpy.array(self.edges.nodes[\
                         self.edges.on_boundary],\
                         copy=1,\
                         order="fortran")
                        
        # compute boundary_element_edges
        boundary_element_edges = _mesher.get_edges_of_elements(boundary_elements,\
                                                      boundary_edges,\
                                                      boundary_element_edges)
        
        # update edges of boundary elements
        self.elements.update_boundary_edges(boundary_element_edges)
        
    ## \brief Plots the degrees of freedom of the problem, for the various
    # variables.
    #
    ## \anchor mesher_plot_nodes_brief <b>Brief description</b>
    #
    # Plots the degrees of freedom of the problem, for the various
    # variables. For each variable is used a different order for each variable,
    # in accordance with the space to which it belongs. It can be in 2D a 0-form, a
    # 1-form or a 2-form; in 3D it can be a 0-form, a 1-form, a 2-form or a 3.form.
    #
    # Typically the spaces are for 2D:
    # 0-forms:
    # \f[
    #     \Lambda^{0}_{p} = \mathrm{span}\left\{h^{p}_{m}(x)h^{p}_{n}(y)\right\}
    # \f]
    #
    # 1-forms:
    # \f[
    #     \Lambda^{1}_{p} = \mathrm{span}\left\{\tilde{h}^{p-1}_{m}(x)h^{p}_{n}(y)\right\}
    #                      \otimes\left\{h^{p}_{k}(x)\tilde{h}^{p-1}_{l}(y)\right\}
    # \f]
    #
    # 2-forms:
    # \f[
    #     \Lambda^{2}_{p} = \mathrm{span}\left\{\tilde{h}^{p-1}_{m}(x)\tilde{h}^{p-1}_{n}(y)\right\}
    # \f]
    #
    # Typically the spaces are for 3D:
    # 0-forms:
    # \f[
    #     \Lambda^{0}_{p} = \mathrm{span}\left\{h^{p}_{l}(x)h^{p}_{m}(y)h^{p}_{n}(z)\right\}
    # \f]
    #
    # 1-forms:
    # \f[
    #     \Lambda^{1}_{p} = \mathrm{span}\left\{\tilde{h}^{p-1}_{l}(x)h^{p}_{m}(y)h^{p}_{n}(z)\right\}
    #                      \otimes\left\{h^{p}_{l}(x)\tilde{h}^{p-1}_{m}(y)h^{p}_{n}(z)\right\}
    #                      \otimes\left\{h^{p}_{l}(x)h^{p}_{m}(y)\tilde{h}^{p-1}_{n}(z)\right\}
    # \f]
    #
    # 2-forms:
    # \f[
    #     \Lambda^{2}_{p} = \mathrm{span}\left\{h^{p}_{l}(x)\tilde{h}^{p-1}_{m}(y)\tilde{h}^{p-1}_{n}(z)\right\}
    #                      \otimes\left\{\tilde{h}^{p-1}_{l}(x)h^{p}_{m}(y)\tilde{h}^{p-1}_{n}(z)\right\}
    #                      \otimes\left\{\tilde{h}^{p-1}_{l}(x)\tilde{h}^{p-1}_{m}(y)h^{p}_{n}(z)\right\}
    # \f]
    #
    # 3-forms:
    # \f[
    #     \Lambda^{3}_{p} = \mathrm{span}\left\{\tilde{h}^{p-1}_{l}(x)\tilde{h}^{p-1}_{m}(y)\tilde{h}^{p-1}_{m}(y)\right\}
    # \f]
    #
    # Where we have used \f$ h^{p}(x)\f$ to designate the pth-order Lagrange
    # interpolant polynomial through the Lobatto points and \f$ \tilde{h}^{p}(x)\f$
    # for the pth-order Lagrange interpolant polynomial through the Gauss points.
    #
    # As an output, this function returns, according to the input parameters, one
    # plot for each approximated variable, showing or not the location of the
    # degrees of freedom and the global numbering of the degrees of freedom.
    #
    ## \anchor mesher_generate_global_numbering_parameters
    #
    ## \param[in] p
    # The order of the interpolant polynomial to use. This will be p for all
    # the ones using Lobatto points and p-1 por the ones using Gauss points. As
    # shown in the description of this function, \ref mesher_plot_nodes_brief "here".
    #
    # \b TYPE: \c INT.
    #
    ## \param[in] zero_form
    # A list of two booleans specifying for the 0-form variables if location of
    # degrees of freedom (first boolean) and/or numbering (second boolean)
    # are to be shown in plot.
    #
    # \b TYPE: (\c BOOL , \c BOOL ).
    #
    ## \param[in] one_form
    # A list of two booleans specifying for the 1-form variables if location of
    # degrees of freedom (first boolean) and/or numbering (second boolean)
    # are to be shown in plot.
    #
    # \b TYPE: (\c BOOL , \c BOOL ).
    #
    ## \param[in] two_form
    # A list of two booleans specifying for the 2-form variables if location of
    # degrees of freedom (first boolean) and/or numbering (second boolean)
    # are to be shown in plot.
    #
    # \b TYPE: (\c BOOL , \c BOOL ).
    #
    def plot_nodes(self, p, zero_form=(True,True),\
                         one_form=(True,True),\
                         two_form=(True,True)):
        # function that plots the mesh edges and the nodes, displaying also
        # the global numbering of the nodes
        
        # check if global numbering was already computed, if not, compute it
        try:
            self.global_numbering_zero_form
            self.global_numbering_one_form_du
        except NameError:
            self.generate_global_numbering(p)
        
        # initialize the plots
        zero_form_figure_number = 1
        zero_form_figure = pylab.figure(zero_form_figure_number)
        one_form_dx_figure_number = 2
        one_form_dx_figure = pylab.figure(one_form_dx_figure_number)
        one_form_dy_figure_number = 3
        one_form_dy_figure = pylab.figure(one_form_dy_figure_number)
        two_form_figure_number = 4
        two_form_figure = pylab.figure(two_form_figure_number)
        zero_form_figure.hold(True)
        one_form_dx_figure.hold(True)
        one_form_dy_figure.hold(True)
        two_form_figure.hold(True)
        
        # plot the mesh
        
        # start with the edges
        # loop over elements and plot the edges, some overlap because they are drawn
        # twice
        for element in range(0, self.N_elements):
            x_coordinates = numpy.concatenate((self.nodes.coordinates[self.elements.nodes[element]][:,0], numpy.array([self.nodes.coordinates[self.elements.nodes[element]][0,0]])))
            y_coordinates = numpy.concatenate((self.nodes.coordinates[self.elements.nodes[element]][:,1], numpy.array([self.nodes.coordinates[self.elements.nodes[element]][0,1]])))
            # plot the edges in zero form plot
            pylab.figure(zero_form_figure_number)
            pylab.plot(x_coordinates, y_coordinates,"k")
            # plot the edges in 1 form dx plot
            pylab.figure(one_form_dx_figure_number)
            pylab.plot(x_coordinates, y_coordinates,"k")
            # plot the edges in 1 form dy plot
            pylab.figure(one_form_dy_figure_number)
            pylab.plot(x_coordinates, y_coordinates,"k")
            # plot the edges in 2 form plot
            pylab.figure(two_form_figure_number)
            pylab.plot(x_coordinates, y_coordinates,"k")
            
        if zero_form[0] == True:
            # plot the degrees of freedom of zero_form
            pylab.figure(zero_form_figure_number)
            pylab.scatter(self.nodes.coordinates[:,0], self.nodes.coordinates[:,1])
        
            # define auxiliary list for node number inside element
            i_index = [0,-1,-1,0]
            j_index = [0,0,-1,-1]
            
        if zero_form[1] == True:
            # add the numbering to the degrees of freedom
            for element in range(0, self.N_elements):
                for n,node in enumerate(self.elements.nodes[element]):
                    pylab.text(self.nodes.coordinates[node,0],\
                               self.nodes.coordinates[node,1],\
                               str(self.global_numbering_zero_form[element,i_index[n],j_index[n]]))
                            
        # generate the inner degrees of freedom
        
        # generate the local coordinate degrees of freedom for order p (Lobatto points)
        zero_form_v_nodes,zero_form_u_nodes = numpy.meshgrid(gauss_quadrature.collocation_points(p), \
                                         gauss_quadrature.collocation_points(p))
        
        # generate the local coordinate degrees of freedom for 1 form dx of
        # order p (Gauss/Lobatto points)
        one_form_v_dx_nodes,one_form_u_dx_nodes = numpy.meshgrid(gauss_quadrature.collocation_points(p), \
                                         gauss_quadrature.collocation_points(p-1, q_type="gauss"))
        
        # generate the local coordinate degrees of freedom for 1 form dy of
        # order p (Lobatto/Gauss points)
        one_form_v_dy_nodes,one_form_u_dy_nodes = numpy.meshgrid(gauss_quadrature.collocation_points(p-1, q_type="gauss"), \
                                         gauss_quadrature.collocation_points(p))
                                        
        # generate the local coordinate degrees of freedom for 2 form of
        # order p (Gauss/Gauss points)
        two_form_v_nodes,two_form_u_nodes = numpy.meshgrid(gauss_quadrature.collocation_points(p-1, q_type="gauss"), \
                                         gauss_quadrature.collocation_points(p-1, q_type="gauss"))
        
        # loop over the elements and plot the inner nodes
        for element in range(0,self.N_elements):
            # compute the trasnfinite mapping
            myTransfinite = math_tools.TransfiniteMapping_2d(self.nodes.coordinates[self.elements.nodes[element]])
            
            # transform the local degrees of freedom to global
            # degrees of freedom using the transfinite mapping
            
            # 0-form
            zero_form_x_nodes = myTransfinite.X(zero_form_u_nodes, zero_form_v_nodes)
            zero_form_y_nodes = myTransfinite.Y(zero_form_u_nodes, zero_form_v_nodes)
            
            # 1-form dx
            one_form_x_dx_nodes = myTransfinite.X(one_form_u_dx_nodes, one_form_v_dx_nodes)
            one_form_y_dx_nodes = myTransfinite.Y(one_form_u_dx_nodes, one_form_v_dx_nodes)
            
            # 1-form dy
            one_form_x_dy_nodes = myTransfinite.X(one_form_u_dy_nodes, one_form_v_dy_nodes)
            one_form_y_dy_nodes = myTransfinite.Y(one_form_u_dy_nodes, one_form_v_dy_nodes)
            
            # 2-form
            two_form_x_nodes = myTransfinite.X(two_form_u_nodes, two_form_v_nodes)
            two_form_y_nodes = myTransfinite.Y(two_form_u_nodes, two_form_v_nodes)
            
            # plot the degrees of freedom
            if zero_form[0] == True:
                if p > 1:
                    # plot the inner degrees of freedom of the 0-form
                    pylab.figure(zero_form_figure_number)
                    pylab.scatter(zero_form_x_nodes[0,1:-1].flatten(),zero_form_y_nodes[0,1:-1].flatten(),color="b")
                    pylab.scatter(zero_form_x_nodes[-1,1:-1].flatten(),zero_form_y_nodes[-1,1:-1].flatten(),color="b")
                    pylab.scatter(zero_form_x_nodes[1:-1,:].flatten(),zero_form_y_nodes[1:-1,:].flatten(),color="b")
            
            if one_form[0] == True:
                # plot the inner degrees of freedom of the 1-form
                # dx part
                pylab.figure(one_form_dx_figure_number)
                pylab.scatter(one_form_x_dx_nodes.flatten(),one_form_y_dx_nodes.flatten(),color="g")
                # dy part
                pylab.figure(one_form_dy_figure_number)
                pylab.scatter(one_form_x_dy_nodes.flatten(),one_form_y_dy_nodes.flatten(),color="r")
            
            if two_form[0] == True:
                # plot the inner degrees of freedom of the 2-form
                pylab.figure(two_form_figure_number)
                pylab.scatter(two_form_x_nodes.flatten(),two_form_y_nodes.flatten(),color="g")
            
            # plot the numbering of the degrees of freedom
                
            if zero_form[1] == True:
                # add the numbering of the degrees of freedom of 0-form
                pylab.figure(zero_form_figure_number)
                for px in range(0,p+1):
                    for py in range(0,p+1):
                        if [px,py] not in [[0,0],[0,p],[p,p],[p,0]]:
                            pylab.text(zero_form_x_nodes[px,py], zero_form_y_nodes[px,py],\
                                       str(self.global_numbering_zero_form[element,px,py]))
                                    
            if one_form[1] == True:
                # add the numbering of the degrees of freedom of 1-form dx part
                pylab.figure(one_form_dx_figure_number)
                for px in range(0,p):
                    for py in range(0,p+1):
                        pylab.text(one_form_x_dx_nodes[px,py], one_form_y_dx_nodes[px,py],\
                                   str(self.global_numbering_one_form_du[element,px,py]))
                # add the numbering of the degrees of freedom of 1-form dy part
                pylab.figure(one_form_dy_figure_number)
                for px in range(0,p+1):
                    for py in range(0,p):
                        pylab.text(one_form_x_dy_nodes[px,py], one_form_y_dy_nodes[px,py],\
                                   str(self.global_numbering_one_form_dv[element,px,py]))
                                
            if two_form[1] == True:
                # add the numbering of the degrees of freedom of 2-form
                pylab.figure(two_form_figure_number)
                for px in range(0,p):
                    for py in range(0,p):
                        pylab.text(two_form_x_nodes[px,py], two_form_y_nodes[px,py],\
                                   str(self.global_numbering_two_form[element,px,py]))
            
        pylab.show()
        
        
        
## \class nodes2d
# \brief Contains all data that defines the %nodes of the %mesh and all methods
# used for the manipulation of the nodes.
#        
class nodes2d:
    # nodes2d class contains the following variables:
    
    ## \var coordinates
    # \brief The list of \em x and \em y coordinates of the main nodes of the
    # %mesh.
    # 
    # It is a 2d array that contains the \em x and \em y coordinates of the
    # main nodes of the mesh. Each row contains the coordinates of a node.
    # The row number correspondes to the node numbering given during the 
    # meshing process. See nodes2d.__compute_coordinates for a complete
    # description of how the coordinates for a rectangular mesh are created.
    # For other types of meshes, the coordinates are just obtained from
    # a file, or meshing lybrary.\n\n
    # \b TYPE: 2d array of \c float64
    
    ## \var bc_markers
    # \brief The list of the boundary condition markers of each main node,
    # identifying which prescribed boundary condition is present.
    #
    # It is a vector whose values have the boundary condition marker that
    # specifies which prescribed boundary condition is to be considered in each
    # main node.
    # 
    # The possible values are:
    #    \li 0 :: node not at the boundary.
    #    \li 1 :: node at a boundary without prescribed boundary conditions.
    #    \li >2 :: node at a boundary with prescribed boundary conditions. In
    #              this case the number obtained for the bc_marker is given by
    #              2 + bc_id, where bc_id is the index number of the boundary
    #              conditions vector.
    # \n\n
    # \b TYPE: vector of \c int32
    
    ## \var where
    # \brief The list of elements that contain each of the main nodes of the
    # %mesh.
    #
    # It is an array whose rows contain the element where the main node, whose
    # index is the row number, is present.
    # 
    # \b TYPE: array of \c int32
    
    ## \brief The constructor of the nodes2d object.
    #
    # For rectangular meshes the constructor calls private methods 
    # nodes2d.__compute_coordinates and nodes2d.__compute_bc_markers to
    # compute the nodes coordinates and the boundary condition markers.
    # For other types of meshes, the constructor just updates the nodes
    # coordinates data and bc_markers.
    #
    ## \param[in] mesh_source
    # Specifies the source for construction of mesh data, can take the values
    # present in #VALID_MESH_SOURCES.\n
    # Notice 
    # that a check is made to confirm that this parameter is present in
    # #VALID_MESH_SOURCES, if not the 
    # object is destroyed and an exception of the type 
    # NereidesExceptions.MeshGenerationError is raised. \n
    # \b TYPE: \c str \n
    #
    ## \param[in] **parameters
    # Consists of the list of paramerters dependent on the mesh_source selected.
    # It is a dictionary. The parameters depending on the mesh_source are:\n
    # \b generate_rectangle_2d \n
    #    \li \c spacings \anchor nodes2d_spacings --> It is a 2d vector whose 
    #                                                 values are the \em x and \em y spacing between
    #                                                 the nodes to construct.\n
    #                                                 \b TYPE: 2d vector of \c int
    #    \li \c mesh_size \anchor nodes2d_mesh_size --> It is a 2d vector whose values are the 
    #                                                   number of elements that the %mesh
    #                                                   has in the \em x and \em y directions.\n
    #                                                   \b TYPE: \c int32
    #    \li \c lower_left_corner  \anchor nodes2d_lower_left_corner --> It is a 2d vector that specifies
    #                                                                    the \em x and \em y coordinates of the
    #                                                                    lower left corner of the %mesh.\n
    #                                                                    \b TYPE: 2d vector of \c float64
    #    \li \c upper_right_corner \anchor nodes2d_upper_right_corner --> It is a 2d vector that specifies the \em x and \em y coordinates of the
    #                                                                     upper right corner of the %mesh.\n
    #                                                                     \b TYPE: 2d vector of \c float64
    #    \li \c boundary_condition_segments --> See \ref mesh_boundary_condition_segments "boundary condition_segments"
    # \n
    #
    # \b read_geompack++_2d \n
    #    \li coordinates --> See nodes2d.coordinates
    #    \li bc_markers --> See nodes2d.bc_markers
    # \n
    #
    # \b TYPE: \c dictionary \n
    #
    ## \param[out] self.coordinates
    # See #coordinates
    #
    ## \param[out] self.bc_markers
    # See #bc_markers
    #
    def __init__(self, mesh_source, **parameters):

        if mesh_source not in VALID_MESH_SOURCES:
            # delete the object since mesh_source is an invalid value
            del(self)
            # raise an Exception of the type 
            # NereidesExceptions.MeshGenerationError
            raise NereidesExceptions.MeshGenerationError(\
                    "Mesh source not valid: " + str(mesh_source) + ". " + \
                    "Must be one of " + str(VALID_MESH_SOURCES))
            return
        
        # check if the parameters are compatible with mesh_source
        
        # parameters for generate_rectangle_2d
        if mesh_source == "generate_rectangle_2d":
            # list of needed parameters for generate_rectangle_2d source
            needed_parameters = ["spacings",\
                                 "mesh_size",\
                                 "lower_left_corner",\
                                 "upper_right_corner",\
                                 "boundary_condition_segments"]
                                
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return
                    
        # parameters for read_geompack++_2d
        if mesh_source == "read_geompack++_2d":
            # list of needed parameters for  read_geompack++_2d
            needed_parameters = ["coordinates",\
                                 "bc_markers"]
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return        
        
        # depending on mesh_source, call the correct node generation function
        # which can be __compute_coordinates, or  other function (which
        # just reads the data given to the initialization of nodes2d)
        
        if (mesh_source == "generate_rectangle_2d"):
            # set mesh dimension
            # self.dimension = 2
            # set element type valid values are the ones in VALID_ELEMENT_TYPES_2d
            # self.element_type = "rectangle"
            # compute the nodes x and y coordinates
            self.__compute_coordinates(parameters["spacings"],\
                                       parameters["mesh_size"],\
                                       parameters["lower_left_corner"],\
                                       parameters["upper_right_corner"])
        
            # compute the nodes boundary condition markers
            self.__compute_bc_markers(parameters["mesh_size"],\
                                      parameters["boundary_condition_segments"])
                                        
        elif (mesh_source == "read_geompack++_2d"):
            # load data of nodes coordinates
            # and the nodes boundary condition markers
            self.coordinates = parameters["coordinates"]
            self.bc_markers = parameters["bc_markers"]
    
    ## \brief Generates the \em x and \em y coordinates of the main nodes of
    # the rectangular %mesh.
    #
    ## \anchor nodes2d_generate_global_numbering_brief <b>Brief description</b>
    # In the case of a rectangular mesh, the main nodes \em x and \em y
    # coordinates are computed and stored in nodes2d.coordinates.
    ## \anchor nodes2d_generate_global_numbering_algorithm \b Algorithm
    #
    # The nodes' coordinates \f$ \mathbf{x}^{n} =(x^{n}_{1},x^{n}_{2}) \f$ 
    # (with 
    # \e n in the range 0 to \f$ N^{\mathrm{nodes}}-1 \f$, the total number of 
    # nodes of the
    # %mesh) are generated in the following way:
    # \f[
    #     x^{n_{i}}_{i} = n_{i}\,\delta x_{i} + P^{llc}_{i}
    # \f]
    #    
    # This is implemented using the \c numpy.mgrid function to generate all the
    # indices \f$ n_{i} \f$ of each node of the %mesh and then multiplying the
    # resulting arrays by \f$ \delta x_{i} \f$ and adding \f$ P^{llc}_{i} \f$.
    #
    # This arrays, with \e x and \e y coordinates of the nodes are ordered in
    # the following way: along the rows, that is, keeping the row index
    # fixed and varying the column index, we move in the \e y direction and
    # along the columns, that is, keeping the column index fixed and varying
    # the row index, we move in the \e x direction.
    #
    ## \anchor nodes2d_generate_global_numbering_definitions \b Definitions
    #
    #    \li \f$ \mathbf{x}^{n} \f$ represents the vector with the coordinates
    #        of a node.
    #    \li \f$ x_{i} \f$ with \f$ i=1,2 \f$ represents the \e x and \e y 
    #         coordinates, respectively.
    #    \li \f$ \Delta x_{i} \f$ represents the \e x and \e y lengths of the
    #         computational domain.
    #    \li \f$ P^{ruc} \f$ represents the Right Upper Corner point.
    #    \li \f$ P^{llc} \f$ represents the Lower Left Corner point.
    #    \li \f$ \delta x_{i} \f$ represents the spacing in the \f$ x_{i} \f$
    #        direction.
    #    \li \f$ N^{e}_{x_{i}} \f$ represents the number of elements in the 
    #        \f$ x_{i} \f$ direction. That is, the \c mesh_size elements.
    #
    ## \anchor nodes2d_generate_global_numbering_parameters
    ## \param[in] spacings
    # See \ref nodes2d_spacings "here"
    ## \param[in] mesh_size
    # See \ref nodes2d_mesh_size "here"
    #
    ## \param[in] lower_left_corner
    # See \ref nodes2d_lower_left_corner "here"
    #
    ## \param[in] upper_right_corner
    # See \ref nodes2d_upper_right_corner "here"
    #
    ## \param[out] self.coordinates
    # See nodes2d.coordinates
    #
    def __compute_coordinates(self, spacings, mesh_size, lower_left_corner,\
                                                         upper_right_corner):
        x_spacing, y_spacing = spacings
        
        # start by generating two arrays of indices using numpy.meshgrid
        x_coordinates,y_coordinates = \
            numpy.meshgrid(numpy.arange(0, mesh_size[0]+1, dtype="float64"), \
                           numpy.arange(0, (mesh_size[1]+1), dtype="float64")) 
                                                      
        # now, just multiply these index matrices by the spacing and
        # adding the coordinates of the lower left corner point, obtaining
        # the coordinates of the main nodes
        # notice that using numpy to perform this multiplication and addition
        # one will have an overhead of memory, since a copy of the array is
        # made and passed to the fortran function which is then returned
        # this means that almost the double of the memory is used. Hence when
        # using very large arrays that are almost the size of the physical
        # memory, care should be taken!
        # that is using:
        #     x_coordinates = x_coordinates*x_spacing + lower_left_corner[0]
        #
        # instead one can use the two operations separately in this fashion
        #    x_coordinates *= x_spacing
        #    x_coordinates += lower_left_corner[0]
        #
        # this way there is no overhead of memory
        
        # compute in two steps:
        # x_coordinates = x_coordinates*x_spacing + lower_left_corner[0]
        x_coordinates *= numpy.float64(x_spacing)
        x_coordinates += numpy.float64(lower_left_corner[0])
        # compute in two steps:
        # y_coordinates = y_coordinates*y_spacing + lower_left_corner[1]
        y_coordinates *= numpy.float64(y_spacing)
        y_coordinates += numpy.float64(lower_left_corner[1])
        
        # number the main nodes and assembly array with x and y coordinates
        # of the nodes
        
        # we use numpy.reshape fuction to take the rows and
        # put them in a vector
        x_coordinates = numpy.array(numpy.reshape(x_coordinates,\
                                            (mesh_size[0]+1)*(mesh_size[1]+1)))
        y_coordinates = numpy.array(numpy.reshape(y_coordinates,\
                                            (mesh_size[0]+1)*(mesh_size[1]+1)))                                        
        
        # allocate memory space for the array that will store the x an y
        # coordinates of the nodes
        self.coordinates = numpy.zeros([(mesh_size[0]+1)*(mesh_size[1]+1),2])
        
        # store the x and y coordinates in the nodes_coordinates matrix and
        # then release memory from x_coordinates and y_coordinates, since
        # they are not needed anymore
        self.coordinates[:,0] = x_coordinates
        self.coordinates[:,1] = y_coordinates
        del(x_coordinates)
        del(y_coordinates)
    
    ## \brief Generates the boundary condition markers for the nodes of the 
    # %mesh.
    #
    ## \anchor nodes2d_generate_global_numbering_brief <b>Brief description</b>
    # The boundary condition markers are computed based upon the boundary
    # segments provided as input. Simply put, all the nodes inside a specific
    # boundary condition segment are labeled with the corresponding boundary
    # condition marker. Notice that the edges will have the same boundary
    # condition marker as the node at its end, that is the second node.
    #
    ## \anchor nodes2d_generate_global_numbering_algorithm \b Algorithm
    #
    # Simply loop over the boundary segments and use the start and end nodes
    # and the side specification and give to those nodes the boundary condition
    # markers specified.
    #
    ## \anchor nodes2d_generate_global_numbering_parameters
    ## \param[in] mesh_size
    # See \ref nodes2d_mesh_size "here"
    #
    ## \param[in] boundary_condition_segments
    # See \ref nodes2d_boundary_condition_segments "here"
    #
    ## \param[out] self.bc_markers
    # See nodes2d.bc_markers
    #    
    def __compute_bc_markers(self, mesh_size, boundary_condition_segments):    
        # Nodes boundary conditions marker
        # allocate memory space for the array that will store the boundary
        # condition marker
        self.bc_markers = numpy.zeros((mesh_size+1),dtype="int16")

        # set all boundary nodes bc_markers to 1, so that one knows
        # they are at the boundary
        # bottom boundary
        self.bc_markers[0:,0] = 1
        # right boundary
        self.bc_markers[mesh_size[0],0:] = 1
        # top boundary
        self.bc_markers[0:,mesh_size[1]] = 1
        # left boundary
        self.bc_markers[0,0:] = 1
        
        # take the information of boundary_condition_segments and 
        # and fill the bc_marker array
        for bc_segment in boundary_condition_segments:
            # the bottom boundary
            if bc_segment[0] == 0:
                self.bc_markers[bc_segment[1]:bc_segment[2],0] = bc_segment[3]+2
            # the right boundary
            if bc_segment[0] == 1:
                self.bc_markers[-1,bc_segment[1]:bc_segment[2]] = bc_segment[3]+2
            # the upper boundary
            if bc_segment[0] == 2:
                # remember that the number of nodes in the x direction is
                # mesh_size[0]+1 and the same for the y direction
                self.bc_markers[(mesh_size[0]+1-bc_segment[2]):\
                        (mesh_size[0]+1-bc_segment[1]),-1] = bc_segment[3]+2
            # the left boundary
            if bc_segment[0] == 3:
                # remember that the number of nodes in the x direction is
                # mesh_size[0]+1 and the same for the y direction
                self.bc_markers[0,(mesh_size[1]+1-bc_segment[2]):\
                        (mesh_size[1]+1-bc_segment[1])] = bc_segment[3]+2
        
        # convert the bc_marker from array to vector putting the columns
        # one next to the other in the same way as with x_coordinates
        # and y_coordinates, since this is the ordering of the main nodes
        # 0,1,2,N_x_elements, for the first line of nodes, the one on the bottom
        # N_x_elements+1, N_x_elements+2, ... for the second line of nodes
        self.bc_markers = numpy.array(numpy.reshape(self.bc_markers,\
                                            (mesh_size[0]+1)*(mesh_size[1]+1),\
                                            order="fortran"))   
        
    ## \brief Simply updates the value of nodes2d.where with the input values
    # given.
    #
    ## \anchor nodes2d_update_where_brief <b>Brief description</b>
    # It is just a function to make the update at nodes2d.where, instead
    # of making the user change the value himself.
    #
    ## \anchor nodes2d_update_where_parameters
    ## \param[in] node_in_elements
    # It is a 2d array where each row contains the index of the elements where
    # the nodes, whose index is the row number, can be found.\n
    # \b TYPE: 2d array of \c int32
    #
    ## \param[out] self.bc_markers
    # See nodes2d.where
    #     
    def update_where(self, node_in_elements):
        # updates the nodes.where attribute, that contains the list, for
        # each node, of the elements where it can be found
        self.where = node_in_elements
            

## \class edges2d
# \brief Contains all data that defines the %edges of the %mesh and all methods
# used for the manipulation of the edges.
#        
class edges2d:
    # edges2d class contains the following variables:
    
    ## \var on_boundary
    # \brief The list of the indices of all edges that lie on a boundary.
    # 
    # \b TYPE: vector of \c int32
    
    ## \var on_prescribed_boundary
    # \brief The list of the indices of all edges that lie on a boundary with
    # prescribed boundary conditions.
    # 
    # \b TYPE: vector of \c int32
    
    ## \var nodes
    # \brief The list of nodes that define all the edges.
    #
    # This array is of type \c int32 and has dimension
    # (N_nodes+N_elements-1) x 2
    # where (N_nodes+N_elements-1) is the number of edges and 2 accounts
    # for the two nodes that identify the edge.
    # 
    # \b TYPE: (N_nodes+N_elements-1) x 2 array of \c int32
    
    ## \var bc_markers
    # \brief The list of boundary condition markers of each edge.
    #
    # This vector is of type \c int32 and has dimension
    # (N_nodes+N_elements-1),
    # where (N_nodes+N_elements-1) is the number of edges.
    #
    # The possible values are:
    #    \li 0 :: node not at the boundary.
    #    \li 1 :: node at a boundary without prescribed boundary conditions.
    #    \li >2 :: node at a boundary with prescribed boundary conditions. In
    #              this case the number obtained for the bc_marker is given by
    #              2 + bc_id, where bc_id is the index number of the boundary
    #              conditions vector, mesh.boundary_conditions.
    # 
    # \b TYPE: (N_nodes+N_elements-1) vector of \c int32
    #
    
    ## \brief The constructor of the edges2d object.
    #
    # For rectangular meshes the constructor calls private methods 
    # edges2d.__generate_edges and edges2d.__generate_bc_markers to
    # compute the nodes that define the edges and the boundary condition
    # markers of the edges. For other types of meshes, the constructor just 
    # updates the nodes data and bc_markers.
    # Additionally it also computes the edges that lie at a boundary
    # and the edges that lie at a boundary with prescribed boundary conditions
    # edges2d.__generate_bc_markers.
    #
    ## \param[in] mesh_source
    # Specifies the source for construction of mesh data, can take the values
    # present in #VALID_MESH_SOURCES.\n
    # Notice 
    # that a check is made to confirm that this parameter is present in
    # #VALID_MESH_SOURCES, if not the 
    # object is destroyed and an exception of the type 
    # NereidesExceptions.MeshGenerationError is raised. \n
    # \b TYPE: \c str \n
    #
    ## \param[in] **parameters
    # Consists of the list of paramerters dependent on the mesh_source selected.
    # It is a dictionary. The parameters depending on the mesh_source are:\n
    # \b generate_rectangle_2d \n
    #    \li \c N_nodes \anchor edges2d_N_nodes --> Number that specifies the number of nodes present on the %mesh.
    #                                                \b TYPE: \c int32
    #    \li \c N_elements \anchor edges2d_N_elements --> Number that specifies the number of elements present on the %mesh.
    #                                                     \b TYPE: \c int32
    #    \li \c mesh_size  \anchor edges2d_mesh_size --> See \ref mesh_mesh_size "here".
    #    \li \c boundary_condition_segments \anchor edges2d_boundary_condition_segments --> See \ref mesh_boundary_condition_segments "here".
    #    \li \c element_type \anchor edges2d_element_type --> See \ref mesh.element_type
    # \n
    #
    # \b read_geompack++_2d \n
    #    \li nodes --> See edges2d.nodes
    #    \li bc_markers --> See edges2d.bc_markers
    # \n
    #
    # \b TYPE: \c dictionary \n
    #
    #
    ## \param[out] self.nodes
    # See #nodes
    #
    ## \param[out] self.bc_markers
    # See #bc_markers
    #
    ## \param[out] self.on_boundary
    # See edges2d.on_boundary.
    #
    ## \param[out] self.on_prescribed_boundary
    # See edges2d.on_prescribed_boundary
    #
    def __init__(self, mesh_source, **parameters):
        if mesh_source not in VALID_MESH_SOURCES:
            # delete the object since mesh_source is an invalid value
            del(self)
            # raise an Exception of the type 
            # NereidesExceptions.MeshGenerationError
            raise NereidesExceptions.MeshGenerationError(\
                    "Mesh source not valid: " + str(mesh_source) + ". " + \
                    "Must be one of " + str(VALID_MESH_SOURCES))
            return
        
        # check if the parameters are compatible with mesh_source
        
        # parameters for generate_rectangle_2d
        if mesh_source == "generate_rectangle_2d":
            # list of needed parameters for generate_rectangle_2d source
            needed_parameters = ["N_nodes", "N_elements", "mesh_size",\
                                 "boundary_condition_segments", "element_type"]
                                
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return
                    
        # parameters for read_geompack++_2d
        if mesh_source == "read_geompack++_2d":
            # list of needed parameters for  read_geompack++_2d
            needed_parameters = ["nodes", "bc_markers"]
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return
        
        # depending on mesh_source, call the correct element generation function
        # which can be __compute_coordinates, or  other function (which
        # just reads the data given to the initialization of elements2d)
        
        if (mesh_source == "generate_rectangle_2d"):
            element_type = parameters["element_type"]
            N_nodes = parameters["N_nodes"]
            N_elements = parameters["N_elements"]
            boundary_condition_segments = parameters["boundary_condition_segments"]
            mesh_size = parameters["mesh_size"]
            # fill in edges data
            if element_type == "rectangle":
                self.__generate_edges(N_nodes, N_elements, mesh_size,\
                                               boundary_condition_segments)
            
            # fill in edges boundary condition markers
            if element_type == "rectangle":
                self.__generate_bc_markers(N_nodes, N_elements,\
                                                    mesh_size,\
                                                    boundary_condition_segments)
            
            # order the edges data
            # putting node_0 < node_1
            self.nodes.sort(axis=1)
            
            # fill the first column, the one with the edge numbers, it is simply
            # going from 0 to the number_of_edges-1 (N_nodes+N_elements-1)
            # this is done only now so that there is no problem with the sort
            self.nodes[:,0] = numpy.arange(0, (N_nodes+N_elements-1), dtype="int32")
            
            # compute the boundary edges
            self.on_boundary = numpy.arange(0, N_elements+N_nodes-1)\
                               [self.bc_markers != 0]
                            
            # compute the edges with prescribed boundary conditions
            self.on_prescribed_boundary = numpy.arange(0, N_elements+N_nodes-1)\
                                          [self.bc_markers > 1]
        
        elif (mesh_source == "read_geompack++_2d"):
            # allocate memory space for self.nodes
            # it contains 3 columns since the first is the number of the edge
            # and the following two are the start and end nodes
            self.nodes = numpy.zeros([parameters["nodes"].shape[0],\
                                      parameters["nodes"].shape[1]+1], dtype="int32")
            # fill the second and third columns, the ones that contain the
            # start and end nodes
            self.nodes[:,1:] = parameters["nodes"]
            # fill the first column, the one with the edge numbers, it is simply
            # going from 0 to the number_of_edges-1 (N_nodes+N_elements-1)
            self.nodes[:,0] = numpy.arange(0, (parameters["nodes"].shape[0]), dtype="int32")
            self.bc_markers = parameters["bc_markers"]
            # compute the boundary edges
            self.on_boundary = numpy.arange(0, self.nodes.shape[0])[self.bc_markers != 0]
            # compute the edges with prescribed boundary conditions
            self.on_prescribed_boundary = numpy.arange(0, self.nodes.shape[0])[self.bc_markers > 1]
            
    ## \brief Generate the edges, that is, the list of nodes that define the
    # edges.
    #
    # First generate the boundary edges, starting with the left boundary, then
    # the right, then the upper and finally the lower boundary.
    #
    # Afterwards, the inner vertical edges are generated, followed by the inner
    # horizontal edges.
    #
    ## \param[in] N_nodes
    # Number that specifies the number of nodes present on the %mesh.
    # \b TYPE: \c int32
    #
    ## \param[in] N_elements
    # Number that specifies the number of elements present on the %mesh.
    # \b TYPE: \c int32
    #
    ## \param[in] mesh_size
    # See \ref mesh_mesh_size "here".
    #
    ## \param[in] boundary_condition_segments \anchor nodes2d_boundary_condition_segments
    # See \ref mesh_boundary_condition_segments "here".
    #
    ## \param[out] self.nodes
    # See edges2d.nodes
    #
    def __generate_edges(self, N_nodes, N_elements, mesh_size,\
                                        boundary_condition_segments): 
        # generate the edges array that contains all the edges that
        # make up the mesh, this array is of type int and has dimension
        # (N_nodes+N_elements-1) x 3
        # where (N_nodes+N_elements-1) is the number of edges and 3 accounts
        # for the two nodes that identify the edge and for the boundary 
        # condition marker
        
        # allocate memory space
        # the first column contains the edge number
        self.nodes = numpy.zeros([(N_nodes+N_elements-1),3],dtype="int32",\
                                 order="fortran")
        
        # first the boundary edges
        # lower boundary
        self.nodes[0:mesh_size[0],0] = numpy.arange(0,mesh_size[0])
        self.nodes[0:mesh_size[0],1] = numpy.arange(1,(mesh_size[0]+1))
        # right boundary
        self.nodes[mesh_size[0]:(mesh_size[0]+mesh_size[1]),0] = \
                numpy.arange(mesh_size[0],N_nodes-1,mesh_size[0]+1)
        self.nodes[mesh_size[0]:(mesh_size[0]+mesh_size[1]),1] = \
                numpy.arange((mesh_size[0]*2 + 1),N_nodes,mesh_size[0]+1)
        # upper boundary
        self.nodes[(mesh_size[0]+mesh_size[1]):(2*mesh_size[0]+mesh_size[1]),0] = \
                numpy.arange(N_nodes-1,N_nodes-mesh_size[0]-1,-1)
        self.nodes[(mesh_size[0]+mesh_size[1]):(2*mesh_size[0]+mesh_size[1]),1] = \
                numpy.arange(N_nodes-2,N_nodes-mesh_size[0]-2,-1)
        # lower boundary
        self.nodes[(2*mesh_size[0]+mesh_size[1]):(2*mesh_size[0]+2*mesh_size[1]),0] = \
                numpy.arange(N_nodes-1-mesh_size[0],0,-(mesh_size[0]+1))
        self.nodes[(2*mesh_size[0]+mesh_size[1]):(2*mesh_size[0]+2*mesh_size[1]),1] = \
                numpy.arange(N_nodes-2-2*mesh_size[0],-1,-(mesh_size[0]+1))
        
        # the vertical inner edges
        # generate the start nodes of the edges
        start_nodes = numpy.arange(1,1+(mesh_size[0]+1)*mesh_size[1],\
                                    mesh_size[0]+1)
        # generate the end nodes of the edges
        end_nodes = numpy.arange(2+mesh_size[0],\
                                1+mesh_size[0]+(mesh_size[0]+1)*mesh_size[1],\
                                mesh_size[0]+1)
        # fill the data for all the vertical edges, the following edges are
        # simply obtained by adding one to each start and end node of the
        # edges of the previous vertical line
        for v_edge in range(0,(mesh_size[0]-1)):
            self.nodes[(2*mesh_size[0]+2*mesh_size[1]+v_edge*mesh_size[1]):\
                  (2*mesh_size[0]+2*mesh_size[1]+(v_edge+1)*mesh_size[1]),0] = \
                  start_nodes + v_edge
            self.nodes[(2*mesh_size[0]+2*mesh_size[1]+v_edge*mesh_size[1]):\
                  (2*mesh_size[0]+2*mesh_size[1]+(v_edge+1)*mesh_size[1]),1] = \
                  end_nodes + v_edge
                
        # the horizontal inner edges
        # the vertical inner edges
        # generate the start nodes of the edges
        start_nodes = numpy.arange(mesh_size[0]+1,\
                                   2*mesh_size[0]+1)
        # generate the end nodes of the edges
        end_nodes = numpy.arange(mesh_size[0]+2,\
                                 2*(mesh_size[0]+1))
        # fill the data for all the horizontal edges, the following edges are
        # simply obtained by adding (1+mesh_size[0]) to each start and end node 
        # of the edges of the previous horizontal line
        # to simplify the numbering we introduce the start edge as the last
        # used edge of the vertical edges
        if mesh_size[0] == 1:
            last_edge = 2*mesh_size[1]+2
        else:
            last_edge = (2*mesh_size[0]+2*mesh_size[1]+(v_edge+1)*mesh_size[1])
        for h_edge in range(0,(mesh_size[1]-1)):
            self.nodes[(last_edge+h_edge*mesh_size[0]):\
                  (last_edge+(h_edge+1)*mesh_size[0]),0] = \
                  start_nodes + h_edge*(1+mesh_size[0])
            self.nodes[(last_edge+h_edge*mesh_size[0]):\
                  (last_edge+(h_edge+1)*mesh_size[0]),1] = \
                  end_nodes + h_edge*(1+mesh_size[0])
    
    ## \brief Generate the boundary condition markers of the edges.
    #
    # All inner edges have boundary marker value equal to zero, since they
    # are not on a boundary, the ones on the boundary (the first ones) have
    # boundary marker values given by the boundary_condition_segments
    # array. All edges inside the segments defined have the corresponding
    # boundary condition, all the others have boundary condition marker
    # value equal to 1, just to know that the value edge lies on an edge.
    # First set the boundary condition marker of all edges that lie in the 
    # boundary of the computational domain to 1.
    #
    ## \param[in] N_nodes
    # Number that specifies the number of nodes present on the %mesh.
    # \b TYPE: \c int32
    #
    ## \param[in] N_elements
    # Number that specifies the number of elements present on the %mesh.
    # \b TYPE: \c int32
    #
    ## \param[in] mesh_size
    # See \ref mesh_mesh_size "here".
    #
    ## \param[in] boundary_condition_segments \anchor nodes2d_boundary_condition_segments
    # See \ref mesh_boundary_condition_segments "here".
    #
    ## \param[out] self.nodes
    # See edges2d.bc_markers
    #
    def __generate_bc_markers(self, N_nodes,\
                                    N_elements,\
                                    mesh_size,\
                                    boundary_condition_segments):    
        # fill in the edges boundary marker part
        # all inner edges have boundary marker value equal to zero, since they
        # are not on a boundary, the ones on the boundary (the first ones) have
        # boundary marker values given by the boundary_condition_segments
        # array. All edges inside the segments defined have the corresponding
        # boundary condition, all the others have boundary condition marker
        # value equal to 1, just to know that the value edge lies on an edge.
        # First set the boundary condition marker of all edges that lie in the 
        # boundary of the computational domain to 1.
        #  first allocate memory space
        self.bc_markers = numpy.zeros([(N_nodes+N_elements-1)],dtype="int32")
        self.bc_markers[0:(2*mesh_size[0]+2*mesh_size[1])] = 1
        
        # take the information of boundary_condition_segments and 
        # and fill the boundary condition marker column of edges array
        for bc_segment in boundary_condition_segments:
            # the bottom boundary
            if bc_segment[0] == 0:
                self.bc_markers[bc_segment[1]:(bc_segment[2]-1)] = \
                                                                bc_segment[3]+2
            # the right boundary
            if bc_segment[0] == 1:
                self.bc_markers[(mesh_size[0]+bc_segment[1]):\
                      (mesh_size[0]+bc_segment[2]-1)] = bc_segment[3]+2
            # the upper boundary
            if bc_segment[0] == 2:
                self.bc_markers[(mesh_size[0]+mesh_size[1]+bc_segment[1]):\
                      (mesh_size[0]+mesh_size[1]+bc_segment[2]-1)] \
                      = bc_segment[3]+2
            # the left boundary
            if bc_segment[0] == 3:
                self.bc_markers[(2*mesh_size[0]+mesh_size[1]+bc_segment[1]):\
                      (2*mesh_size[0]+mesh_size[1]+bc_segment[2]-1)] \
                      = bc_segment[3]+2
                    
        
## \class elements2d
# \brief Contains all data that defines the %edges of the %mesh and all methods
# used for the manipulation of the edges.
#        
class elements2d:
    # elements2d class contains the following variables:
    
    ## \var on_boundary
    # \brief The list of the indices of all elements that lie on a boundary.
    # 
    # \b TYPE: vector of \c int32
    
    ## \var type
    # \brief The type of element used to construct the %mesh.
    # 
    # Can take any value
    # included in #VALID_ELEMENT_TYPES_2d.\n\n
    # \b TYPE: \c str
    
    ## \var nodes
    # \brief The list of nodes that define all the elements.
    #
    # This array is of type \c int32 and has dimension
    # N_elements x N_nodes_per_element. For the case of quadrilateral elements
    # N_nodes_per_element is 4, for triangular is 3, and so on. The nodes are
    # specified counterclockwise, starting at 0.
    # 
    # \b TYPE: N_elements x N_nodes_per_element array of \c int32
    
    ## \var neighbors
    # \brief The list of the neighbors of each element.
    #
    # It is an array of dimension N_elements x N_nodes_per_element.
    # The neighbors are specified in such a way as that in the case of 
    # triangular elements, the neighbor 0 is the one that is opposite to the
    # node 0, and so on. Hence for the case of quarilateral elements one has
    # that the neighbor 0 is the one that shares the edge 0, that is the one
    # defined by the nodes 1 and 2. Hence, in general, neighbor nh is the one
    # that shares the edge defined by the nodes nh+1 and nh+2. Notice that if
    # nh+2 bigger than the highest index of an element node, than nh+2 is 0.
    # 
    # \todo Make a figure with an element and its neighbors, specifying which
    # ones are neighbor 0, 1, 2 and 3.
    # 
    # \b TYPE: N_elements x N_nodes_per_element array of \c int32
    
    ## \var boundary_edges
    #
    # \brief The list of the edges of the elements that lie at the boundary.
    #
    # It is an array of dimension N_elements x N_nodes_per_element.
    # I
    # 
    # \b TYPE: N_elements x N_nodes_per_element array of \c int32
    
    ## \var mesh_source
    #
    # \brief Speficies the source of #mesh construction for the #mesh data used.
    #
    # Specifies the source for construction of #mesh data, can take the values
    # present in #VALID_MESH_SOURCES.\n
    # Notice 
    # that a check is made to confirm that this parameter is present in
    # #VALID_MESH_SOURCES, if not the 
    # object is destroyed and an exception of the type 
    # NereidesExceptions.MeshGenerationError is raised. \n
    # \b TYPE: \c str 
    
    ## \brief The constructor of the elements2d object.
    #
    # For rectangular meshes the constructor calls private method
    # elements2d.__generate_elements to
    # compute the nodes that define the elements. For other types of meshes,
    # the constructor just 
    # updates the nodes data.
    # Additionally it also computes the neighbors and for each element the edges
    # that lie at the boundary, edges2d.__generate_neighbors and
    # edges2d.update_boundary_edges.
    #
    ## \param[in] mesh_source
    # Specifies the source for construction of mesh data, can take the values
    # present in #VALID_MESH_SOURCES.\n
    # Notice 
    # that a check is made to confirm that this parameter is present in
    # #VALID_MESH_SOURCES, if not the 
    # object is destroyed and an exception of the type 
    # NereidesExceptions.MeshGenerationError is raised. \n
    # \b TYPE: \c str \n
    #
    ## \param[in] **parameters
    # Consists of the list of paramerters dependent on the mesh_source selected.
    # It is a dictionary. The parameters depending on the mesh_source are:\n
    # \b generate_rectangle_2d \n
    #    \li \c N_elements \anchor edges2d_N_elements --> Number that specifies the number of elements present on the %mesh.
    #                                                     \b TYPE: \c int32
    #    \li \c mesh_size  \anchor edges2d_mesh_size --> See \ref mesh_mesh_size "here".
    #    \li \c element_type \anchor edges2d_element_type --> See \ref mesh.element_type
    # \n
    #
    # \b read_geompack++_2d \n
    #    \li nodes --> See elements2d.nodes
    #    \li edges_nodes --> See edges2d.nodes
    # \n
    #
    # \b TYPE: \c dictionary \n
    #
    ## \param[out] self.nodes
    # See elements2d.nodes
    #
    ## \param[out] self.neighbors
    # See elements2d.neighbors
    #
    ## \param[out] self.on_boundary
    # See elements2d.on_boundary
    #
    ## \param[out] self.type
    # See elements2d.type
    #
    def __init__(self, mesh_source, **parameters):

        if mesh_source not in VALID_MESH_SOURCES:
            # delete the object since mesh_source is an invalid value
            del(self)
            # raise an Exception of the type 
            # NereidesExceptions.MeshGenerationError
            raise NereidesExceptions.MeshGenerationError(\
                    "Mesh source not valid: " + str(mesh_source) + ". " + \
                    "Must be one of " + str(VALID_MESH_SOURCES))
            return
        
        # check if the parameters are compatible with mesh_source
        
        # parameters for generate_rectangle_2d
        if mesh_source == "generate_rectangle_2d":
            # list of needed parameters for generate_rectangle_2d source
            needed_parameters = ["N_elements", "mesh_size", "element_type"]
                                
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return
                    
        # parameters for read_geompack++_2d
        if mesh_source == "read_geompack++_2d":
            # list of needed parameters for  read_geompack++_2d
            needed_parameters = ["nodes", "edges_nodes"]
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return
        
        # depending on mesh_source, call the correct element generation function
        # which can be __compute_coordinates, or  other function (which
        # just reads the data given to the initialization of elements2d)
        
        if (mesh_source == "generate_rectangle_2d"):
            mesh_size = parameters["mesh_size"]
            N_elements = parameters["N_elements"]
            element_type = parameters["element_type"]
            
            # define the type of elements
            self.type = element_type
            # define the generation method
            self.mesh_source = mesh_source
            # generate elements
            self.__generate_elements(N_elements, mesh_size)
            # generate neighbors data
            self.__generate_neighbors(N_elements=N_elements, mesh_size=mesh_size)
            # compute the boundary elements
            self.on_boundary = numpy.arange(0, N_elements)\
                               [(self.neighbors==-1).any(axis=1)]
                            
        elif (mesh_source == "read_geompack++_2d"):
            # define the type of elements
            self.type = "rectangle"
            # define the generation method
            self.mesh_source = mesh_source
            # update element nodes
            self.nodes = parameters["nodes"]
            # compute neighbors
            self.__generate_neighbors(edges_nodes=parameters["edges_nodes"])
            # compute the boundary elements
            self.on_boundary = numpy.arange(0, self.nodes.shape[0])\
                               [(self.neighbors==-1).any(axis=1)]
    
    ## \brief Simply updates the value of boundary_edges
    #
    ## \anchor elements2d_update_boundary_edges_brief <b>Brief description</b>
    # Update the boundary edges information, specifying, for each element, 
    # which edges are present and lie at the boundary.
    # It is an array whose dimensions are N_elements at boundary x 4.
    # On the rows one will have the edge number or -1 if the cooresponding
    # edge is not at the boundary. Using the usual opposite numbering
    # described in the edges class.
    #
    ## \anchor elementsed_update_boundary_edges_parameters
    ## \param[in] edges
    # See elements2d.boundary_edges
    #
    ## \param[out] self.boundary_edges
    # See elements2d.boundary_edges
    #    
    def update_boundary_edges(self, edges):
        # Update the boundary edges information, specifying, for each element, 
        # which edges are present and lie at the boundary.
        # It is an array whose dimensions are N_elements at boundary x 4.
        # On the rows one will have the edge number or -1 if the cooresponding
        # edge is not at the boundary. Using the usual opposite numbering
        # described in the edges class.
        self.boundary_edges = edges
    
    ## \brief Computes the elements.
    #
    ## \anchor elements2d_generate_elements_brief <b>Brief description</b>
    # Generate the elements array that contains all the information that
    # defines the elements, specifying, for each element, the nodes that
    # are contained in it. This is an array of int whose dimensions are
    # N_elements x 4.
    # 4 accounts for the 4 nodes that make up the element, recall that the
    # element is rectangular hence has 4 vertices.
    #
    # This method is only for generating rectangular meshes, hence the above
    # specification.
    #
    ## \anchor elementsed_update_boundary_edges_parameters
    ## \param[in] N_elements
    # Number that specifies the number of elements present on the %mesh.
    # \b TYPE: \c int32
    #
    ## \param[in] mesh_size
    # See \ref mesh_mesh_size "here".
    #
    ## \param[out] self.nodes
    # See elements2d.nodes
    #
    def __generate_elements(self, N_elements, mesh_size):    
        # Generate the elements array that contains all the information that
        # defines the elements, specifying, for each element, the nodes that
        # are contained in it. This is an array of int whose dimensions are
        # N_elements x 4.
        # 4 accounts for the 4 nodes that make up the element, recall that the
        # element is rectangular hence has 4 vertices.
        
        # allocate memory space
        self.nodes = numpy.zeros([N_elements,4], dtype="int32", order="fortran")
        
        # One must notice that the same node of consecutive elements is a
        # consecutive number, that is, if the node i (0 to 3) of an element,
        # whose number is e, has the number n then the node i of the 
        # consecutive element, numbered e+1, will be n+1.
        for row in range(0, mesh_size[1]):
            # node 0
            self.nodes[(row*mesh_size[0]):((row+1)*mesh_size[0]), 0] = \
                        numpy.arange(row*(mesh_size[0]+1),\
                                    ((row+1)*mesh_size[0]+row), dtype="int32")
            # node 1
            self.nodes[(row*mesh_size[0]):((row+1)*mesh_size[0]), 1] = \
                        numpy.arange(row*(mesh_size[0]+1)+1,\
                                    (row+1)*(mesh_size[0]+1), dtype="int32")
            # node 2
            self.nodes[(row*mesh_size[0]):((row+1)*mesh_size[0]), 2] = \
                        numpy.arange(row*(mesh_size[0]+1)+2+mesh_size[0],\
                                    (row+1)*(mesh_size[0]+1)+mesh_size[0]+1,\
                                    dtype="int32")
            # node 3
            self.nodes[(row*mesh_size[0]):((row+1)*mesh_size[0]), 3] = \
                        numpy.arange(row*(mesh_size[0]+1)+1+mesh_size[0],\
                                    (row+1)*(mesh_size[0]+1)+mesh_size[0],\
                                    dtype="int32")
    
    ## \brief Computes the neighbors of the elements.
    #
    ## \anchor elements2d_generate_neighbors_brief <b>Brief description</b>
    # Generate the neighbors array that contains all the information that
    # defines the neighbors of each element.
    #
    # -1 means that there is no neighbor on that position. See 
    # elements2d.neighbors for a specification on how the neighbors are
    # identified.
    #
    # This method is only for generating rectangular meshes, hence the above
    # specification.
    #
    ## \anchor elementsed_update_boundary_edges_parameters
    ## \param[in] N_elements
    # Number that specifies the number of elements present on the %mesh.
    # \b TYPE: \c int32
    #
    ## \param[in] mesh_size
    # See \ref mesh_mesh_size "here".
    #
    ## \param[out] self.neighbors
    # See elements2d.neighbors
    #    
    def __generate_neighbors(self, **parameters):
        # generate the neighbors information
        # -1 means that there is no neighbor

        # check if the parameters are compatible with mesh_source
        
        # parameters for generate_rectangle_2d
        if self.mesh_source == "generate_rectangle_2d":
            # list of needed parameters for generate_rectangle_2d source
            needed_parameters = ["N_elements", "mesh_size"]
                                
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return
                    
        # parameters for read_geompack++_2d
        if self.mesh_source == "read_geompack++_2d":
            # list of needed parameters for generate_rectangle_2d source
            needed_parameters = ["edges_nodes"]
                                
            # loop over all needed parameters and check if they exist
            for parameter in needed_parameters:
                if (parameters.has_key(parameter) == False):
                    # if not, delete the object
                    del(self)
                    # raise an Exception of type NereidesExceptions.MeshGenerationError
                    # stating which parameter is missing
                    raise NereidesExceptions.MeshGenerationError("Parameter for mesh_source \
                          is missing: " + str(parameter) + "\n")
                    return
        
        # different ways of generating the neighbors
        
        if self.mesh_source == "generate_rectangle_2d":
            N_elements = parameters["N_elements"]
            mesh_size = parameters["mesh_size"]
            # allocate memory space
            self.neighbors = numpy.zeros([N_elements, 4], dtype="int32")
            
            if (mesh_size[0] == 1) and (mesh_size[1] == 1):
                self.neighbors[0,:] = [-1, -1, -1, -1]
                
            elif (mesh_size[0] == 1) and (mesh_size[1] != 1):
                elements_list =  numpy.arange(0,N_elements)
                # lower element
                self.neighbors[0,:] = [-1, 1, -1, -1]
                # upper element
                self.neighbors[N_elements-1,:] = [-1, -1, -1, N_elements-2]
                # inner elements
                # the right neighbors
                self.neighbors[1:-1,0] = -1
                # the upper neighbors
                self.neighbors[1:-1,1] = elements_list[2:]
                # the left neighbors
                self.neighbors[1:-1,2] = -1
                # the lower neighbors
                self.neighbors[1:-1,3] = elements_list[0:-2]
                
            elif (mesh_size[0] != 1) and (mesh_size[1] == 1):
                elements_list =  numpy.arange(0,N_elements)
                # left element
                self.neighbors[0,:] = [1, -1, -1, -1]
                # right element
                self.neighbors[N_elements-1,:] = [-1, -1, N_elements-2, -1]
                # inner elements
                # the right neighbors
                self.neighbors[1:-1,0] = elements_list[2:]
                # the upper neighbors
                self.neighbors[1:-1,1] = -1
                # the left neighbors
                self.neighbors[1:-1,2] = elements_list[0:-2]
                # the lower neighbors
                self.neighbors[1:-1,3] = -1
                
            else:
                # first fill in the neighbors for the corner points
                # the lower left corner element (element 0) only has the upper and
                # right neighbors
                self.neighbors[0,:] = [1, mesh_size[0], -1, -1]
                # the lower right corner element (element mesh_size[0]-1) only has the
                # upper and the left neighbors
                self.neighbors[mesh_size[0]-1,:] = [-1, 2*mesh_size[0]-1,\
                                                         mesh_size[0]-2, -1]
                # the upper right corner element (element N_elements-1) only has the
                # left and lower neighbors
                self.neighbors[N_elements-1,:] = [-1, -1, N_elements-2,\
                                                           N_elements-mesh_size[0]-1]
                # the upper left corner elements (element N_elements-mesh_size[0]) only
                # has the  right and lower neighbors
                self.neighbors[N_elements-mesh_size[0],:] = \
                                                [N_elements-mesh_size[0]+1,\
                                                 -1, -1, N_elements-2*mesh_size[0]]
                
                # fill in information for the nodes at the left boundary
                elements_list =  numpy.arange(0,N_elements-mesh_size[0]+1,mesh_size[0])
                # neighbors on the right
                self.neighbors[elements_list[1:-1],0] = elements_list[1:-1]+1 
                # neighbors up
                self.neighbors[elements_list[1:-1],1] = elements_list[2:]
                # neighbors left
                self.neighbors[elements_list[1:-1],2] = -1
                # neighbors down
                self.neighbors[elements_list[1:-1],3] = elements_list[0:-2]
                
                # fill in information for the nodes at the right boundary
                elements_list =  numpy.arange(mesh_size[0]-1,N_elements,mesh_size[0])
                
                # neighbors on the right
                self.neighbors[elements_list[1:-1],0] = -1
                # neighbors up
                self.neighbors[elements_list[1:-1],1] = elements_list[2:]
                # neighbors left
                self.neighbors[elements_list[1:-1],2] = elements_list[1:-1]-1
                # neighbors down
                self.neighbors[elements_list[1:-1],3] = elements_list[0:-2]
                
                # fill in information for the nodes at the bottom boundary
                elements_list =  numpy.arange(0,mesh_size[0])
                # neighbors on the right
                self.neighbors[elements_list[1:-1],0] = elements_list[1:-1]+1
                # neighbors up
                self.neighbors[elements_list[1:-1],1] = elements_list[1:-1]+\
                                                         mesh_size[0]
                # neighbors left
                self.neighbors[elements_list[1:-1],2] = elements_list[1:-1]-1
                # neighbors down
                self.neighbors[elements_list[1:-1],3] = -1
                
                # fill in information for the nodes at the upper boundary
                elements_list =  numpy.arange(N_elements-mesh_size[0], N_elements)
                # neighbors on the right
                self.neighbors[elements_list[1:-1],0] = elements_list[1:-1]+1
                # neighbors up
                self.neighbors[elements_list[1:-1],1] = -1
                # neighbors left
                self.neighbors[elements_list[1:-1],2] = elements_list[1:-1]-1
                # neighbors down
                self.neighbors[elements_list[1:-1],3] = elements_list[1:-1]-\
                                                         mesh_size[0]
                                                        
                # fill in information for the inner nodes
                for start_element in range(1,mesh_size[0]-1):
                    elements_list = numpy.arange(start_element,\
                                            mesh_size[0]*(mesh_size[1]-1)+start_element+1,\
                                            mesh_size[0])
                    # neighbors on the right
                    self.neighbors[elements_list[1:-1],0] = elements_list[1:-1]+1
                    # neighbors up
                    self.neighbors[elements_list[1:-1],1] = elements_list[1:-1]+\
                                                             mesh_size[0]
                    # neighbors left
                    self.neighbors[elements_list[1:-1],2] = elements_list[1:-1]-1
                    # neighbors down
                    self.neighbors[elements_list[1:-1],3] = elements_list[1:-1]-\
                                                             mesh_size[0]                                        
            
        if self.mesh_source == "read_geompack++_2d":
            edges = parameters["edges_nodes"]
            # Each edge is shared by at most 2 elements
            # hence an array which has for each edge the elements where
            # it can be found has dimension N_edges x 2.
            # 
            # If the edge is present only in one element the other is -1.
            # Allocate memory space for array that stores the elemens where the edges
            # can be found
            where_edges = numpy.zeros([edges.shape[0],2], dtype="int32") - 1

            # allocate memory space for the array, of dimensions N_edges x 2, that stores
            # for each edge the position in the element, which corresponds to the neighbor
            # numbering
            edges_element_position = numpy.zeros([edges.shape[0],2], dtype="int32") - 1

            # allocate memory space for the vector that keeps track of the number of
            # elements found for each edge
            edges_index = numpy.zeros(edges.shape[0],dtype="int32")

            # store the edges in a sparse array
            # the i and j indices of the array corresponde to the start and end node of the
            # edge. The value of the array is the edge number. Hence the array is symmetric.
            # pysparse.spmatrix.ll_mat_sym writes only to the lower triangle part.
            look_up_edges = pysparse.spmatrix.ll_mat_sym(edges.max()+1, edges.shape[0])

            # fill in look_up_edges sparse matrix
            # go along the edges
            for edge,nodes in enumerate(edges):
                if nodes[0] >= nodes[1]:
                    look_up_edges[nodes[0],nodes[1]] = edge
                else:
                    look_up_edges[nodes[1],nodes[0]] = edge
            
            # go along the elements and determine in which element each edge exists
            for element,nodes in enumerate(self.nodes):
                # node 0 and node 1
                # the edge number
                edge = look_up_edges[nodes[0],nodes[1]]
                # update the where_edges array
                where_edges[edge, edges_index[edge]] = element
                # update the location of the edge in the element, to compute the neighbors
                edges_element_position[edge, edges_index[edge]] = 3
                # update the number of elements found for the edge
                edges_index[edge] += 1
                
                # node 1 and node 2
                # the edge number
                edge = look_up_edges[nodes[1],nodes[2]]
                # update the where_edges array
                where_edges[edge, edges_index[edge]] = element
                # update the location of the edge in the element, to compute the neighbors
                edges_element_position[edge, edges_index[edge]] = 0
                # update the number of elements found for the edge
                edges_index[edge] += 1
                
                # node 2 and node 3
                # the edge number
                edge = look_up_edges[nodes[2],nodes[3]]
                # update the where_edges array
                where_edges[edge, edges_index[edge]] = element
                # update the location of the edge in the element, to compute the neighbors
                edges_element_position[edge, edges_index[edge]] = 1
                # update the number of elements found for the edge
                edges_index[edge] += 1
                
                # node 3 and node 0
                # the edge number
                edge = look_up_edges[nodes[3],nodes[0]]
                # update the where_edges array
                where_edges[edge, edges_index[edge]] = element
                # update the location of the edge in the element, to compute the neighbors
                edges_element_position[edge, edges_index[edge]] = 2
                # update the number of elements found for the edge
                edges_index[edge] += 1
                
            # now build up the neighbors list for each element

            # neighbors are numbered in the following way:
            # neighbor 0: the one that shares nodes 1 and 2
            # neighbor 1: the one that shares nodes 2 and 3
            # neighbor 2: the one that shares nodes 3 and 0
            # neighbor 3: the one that shares nodes 0 and 1

            # this is done in this way to be compatible with neighbor numbering for triangles
            # where the neighbor is numbered according to the opposing node

            # allocate memory space for neighbors list -1 means that there is no neighbor
            self.neighbors = numpy.zeros([self.nodes.shape[0], 4], dtype="int64") - 1

            # go along the where_edges array and build up the neighbors list
            # for each element
            for edge,elements in enumerate(where_edges):
                if (elements[0] != -1) and  (elements[1] != -1):
                    self.neighbors[elements[0], edges_element_position[edge,0]] = elements[1]
                    self.neighbors[elements[1], edges_element_position[edge,1]] = elements[0]
        
#------------------------------------------------------------------------------
# FUNCTIONS
#------------------------------------------------------------------------------



    
    
    
    
    
    
    
               