!File _mesher.f90
! Nereides Least Square Spectral Element Flow Solver
! mesher module fortran routines

subroutine get_elements (mesh,mesh_n,mesh_m,node_in_elements,node_n,node_m)
    ! compute for each node of the mesh, which elements contain the node
    implicit none
    ! mesh is an array whose dimensions are N_elements x 4 that contains
    ! in each row the nodes that make up the elements of the mesh
    integer mesh(mesh_n,mesh_m)

    ! node_in_elements is an array whose dimensions are N_nodes x (node_max_repetition+1)
    ! where node_max_repetition is the maximum number of times the nodes are repeated
    ! in different elements. Notice that the last column contains, for each node, 
    ! the number of elements found until the moment,
    ! that contain the given node. This is used to keep track while filling the
    ! node_in_elements array.
    integer node_in_elements(node_n,node_m)

    ! mesh_n and mesh_m are the x and y dimensions of the mesh array
    ! these are passed automatically by python
    integer mesh_n,mesh_m

    ! none_n and node_m are the x and y dimensions of the mesh node_in_elements
    ! these are passed automatically by python
    integer node_n,node_m

    ! node is an index to run over the nodes in a do loop
    integer node

    ! element is an index to run over the elements in a do loop
    integer element

    ! last_column defines the index of the last column of node_in_elements
    ! which is the column that holds the number of elements found for each element
    ! until the moment
    integer last_column

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: mesh
    !f2py intent(in,out) :: node_in_elements
    !f2py intent(hide) :: mesh_n,mesh_m,node_n,node_m
     
    ! define the last_column
    last_column = node_m

    ! run over all nodes
    do node=1, mesh_m
        ! of all elements
        do element=1, mesh_n
            ! assign the value of the current element to the node entry in node_in_elements
            ! in this way registing that the node is at the element
            node_in_elements(mesh(element,node)+1,node_in_elements(mesh(element,node)+1, last_column)+1) = element-1
            ! increment the node_in_how_many entry of the current node, to keep track of in how many
            ! elements the node was found
            node_in_elements(mesh(element,node)+1, last_column) = node_in_elements(mesh(element,node)+1, last_column) + 1
        end do
    end do 
end


subroutine get_edges_of_elements(mesh, mesh_n, mesh_m, edges, edges_n, edges_m,  element_edges, element_edges_n, element_edges_m)
    ! given a mesh, compute for each element the edges that define it
    implicit none
    ! mesh is an array whose dimensions are N_elements x 4 that contains
    ! in each row the nodes that make up the elements of the mesh
    integer mesh(mesh_n,mesh_m)
    
    ! edges is an array whose dimensions are N_edges x 2, where N_edges is
    ! equal to the number of edges of the mesh and is given by 
    !     N_elements + N_nodes -1
    ! Contains in each row the nodes that make up the edge that corresponds to
    ! that row.
    integer edges(edges_n,edges_m)

    ! element_edges is an array whose dimensions are N_elements x 4.
    ! Contains in each row the edges that make up the element that corresponds
    ! to that row. This is what is computed in this function and is the output
    integer element_edges(element_edges_n,element_edges_m)
    
    ! mesh_n and mesh_m are the x and y dimensions of the mesh array
    ! these are passed automatically by python
    integer mesh_n,mesh_m
    
    ! edge_n and edge_m are the x and y dimensions of the edges array
    ! these are passed automatically by python
    integer edges_n,edges_m

    ! element_edges_n and element_edges_m are the x and y dimensions of the
    ! element_edges array, these are passed automatically by python
    integer element_edges_n, element_edges_m

    ! N_elements is the number of elements in the mesh, it is equal to mesh_n
    integer N_elements
    integer edge
    
    ! boundary_element is a do loop iterator that contains the current number of the
    ! boundary_element in the iteration
    integer boundary_element

    ! nodes is a vector whose size is given by mesh_m, that is, the number
    ! of nodes per element is is used in the do loop to store 
    integer nodes(mesh_m)
    
    ! temporary nodes
    integer temp_node_0, temp_node_1
    
    integer jj

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: mesh, edges
    !f2py intent(in,out) :: element_edges
    !f2py intent(hide) :: mesh_n,mesh_m,edges_n,edges_m,element_edges_n,element_edges_m

    N_elements = mesh_n

    if(mesh_m == 4) then
        do boundary_element = 1, N_elements
            nodes = mesh(boundary_element,:)
            ! notice that:
            !     edge 1 is the one with nodes 2 and 3
            !     edge 2 is the one with nodes 3 and 4
            !     edge 3 is the one with nodes 4 and 1
            !     edge 4 is the one with nodes 1 and 2
            ! this is done do compatilibize with triangle elements
            ! where the edge i of the element is the one oposite to node
            ! i of the element.
            
            ! find which edge contains the nodes 2 and 3 of the element
            ! first order the two nodes
            if(nodes(2) < nodes(3)) then
                temp_node_0 = nodes(2)
                temp_node_1 = nodes(3)
            else
                temp_node_0 = nodes(3)
                temp_node_1 = nodes(2)
            endif
            do edge=1, edges_n
                if(edges(edge,2) == temp_node_0) then
                    if(edges(edge,3) == temp_node_1) then
                        element_edges(boundary_element, 1) = edges(edge,1)
                        exit
                    endif
                endif
            end do
            
            ! find which edge contains the nodes 3 and 4 of the element
            ! first order the two nodes
            if(nodes(3) < nodes(4)) then
                temp_node_0 = nodes(3)
                temp_node_1 = nodes(4)
            else
                temp_node_0 = nodes(4)
                temp_node_1 = nodes(3)
            endif
            do edge=1, edges_n
                if(edges(edge,2) == temp_node_0) then
                    if(edges(edge,3) == temp_node_1) then
                        element_edges(boundary_element, 2) = edges(edge,1)
                        exit
                    endif
                endif
            end do
            
            ! find which edge contains the nodes 4 and 1 of the element
            ! first order the two nodes
            if(nodes(4) < nodes(1)) then
                temp_node_0 = nodes(4)
                temp_node_1 = nodes(1)
            else
                temp_node_0 = nodes(1)
                temp_node_1 = nodes(4)
            endif
            do edge=1, edges_n
                if(edges(edge,2) == temp_node_0) then
                    if(edges(edge,3) == temp_node_1) then
                        element_edges(boundary_element, 3) = edges(edge,1)
                        exit
                    endif
                endif
            end do
            
            ! find which edge contains the nodes 1 and 2 of the element
            ! first order the two nodes
            if(nodes(1) < nodes(2)) then
                temp_node_0 = nodes(1)
                temp_node_1 = nodes(2)
            else
                temp_node_0 = nodes(2)
                temp_node_1 = nodes(1)
            endif
            do edge=1, edges_n
                if(edges(edge,2) == temp_node_0) then
                    if(edges(edge,3) == temp_node_1) then
                        element_edges(boundary_element, 4) = edges(edge,1)
                        exit
                    endif
                endif
            end do
        end do
    
    elseif(mesh_m == 3) then
        write(*,*) "It is a triangle!"
    endif
end
