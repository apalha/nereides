import numpy

def boundary_mesh_intersections(line_nodes, myMesh, N_x_elements, N_y_elements,\
                                lower_left_corner, upper_right_corner):
    
    # compute the intersection with the elements
    N_line_nodes = line_nodes.shape[0]

    # loop over the line segments
    for line in range(0,N_line_nodes):
        # compute line length
        line_l = numpy.sqrt(((line_nodes[line,0] - line_nodes[(line+1)%N_line_nodes,0])**2) +\
                 ((line_nodes[line,0] - line_nodes[(line+1)%N_line_nodes,0])**2))
                
        # compute an estimate of the number of intersection points
        N_intersects_estimate = 2*max(int(numpy.ceil(line_l/myMesh.x_spacing)),\
                                    int(numpy.ceil(line_l/myMesh.x_spacing))) + 1
        
        # allocate memory space for a numpy array to hold the element number
        # of the intersected elements
        elements_intersected = numpy.zeros(N_intersects_estimate, dtype="int64")
        
        # allocate memory space for the arrays that hold the start and end coordinates
        # of the intersection with the elements
        start_coordinates_intersect = numpy.zeros([N_intersects_estimate,2], dtype="float64")
        end_coordinates_intersect = numpy.zeros([N_intersects_estimate,2], dtype="float64")
        
        # find in which element the starting point lies
        e_x = int(numpy.floor((line_nodes[line,0] - lower_left_corner[0])/myMesh.x_spacing))
        e_y = int(numpy.floor((line_nodes[line,1] - lower_left_corner[1])/myMesh.y_spacing))
        
        element_start = e_x + N_x_elements*e_y
        
        # find in which element the ending point lies
        e_x = int(numpy.floor((line_nodes[(line+1)%N_line_nodes,0] - lower_left_corner[0])/myMesh.x_spacing))
        e_y = int(numpy.floor((line_nodes[(line+1)%N_line_nodes,1] - lower_left_corner[1])/myMesh.y_spacing))
        
        element_end = e_x + N_x_elements*e_y

        #print "line %d :: from e%d --> e%d" % (line,element_start, element_end)
        
        element = element_start
        
        x0_intersect = line_nodes[line,0]
        y0_intersect = line_nodes[line,1]
        
        not_edge = -1
        
        # initialize counter that keeps track of how many intersections where
        # computed
        intersects_counter = 0
        
        # keep on looping and jumping from element to element
        # until element is element_end, that is the element where the end node
        # of the line lies
        while True:
            # locally store the element_nodes
            element_nodes = myMesh.elements.nodes[element]
            
            # store the element number for the intersection
            elements_intersected[intersects_counter] = element
            
            # store start coodinates of intersection
            start_coordinates_intersect[intersects_counter,:] = numpy.array([x0_intersect,\
                                                                           y0_intersect])
            
            if element != element_end:
                # compute the first intersection with the element edges
                # loop over the element edges
                
                # edge 0
                edge = 0
                if edge != not_edge:
                    # compute parameter value at which intersection occurs
                    t = (myMesh.nodes.coordinates[element_nodes[edge],1] - line_nodes[line,1]) \
                            /(line_nodes[(line+1)%N_line_nodes,1] - line_nodes[line,1])
                    
                    if t >= 0:
                        x_temp_intersect = t*(line_nodes[(line+1)%N_line_nodes,0]-line_nodes[line,0]) + line_nodes[line,0]
                    
                        #    check if in fact is an intersection
                        if myMesh.nodes.coordinates[element_nodes[edge],0] <= x_temp_intersect <= \
                            myMesh.nodes.coordinates[element_nodes[(edge+1)%4],0]:
                            
                            x_intersect = x_temp_intersect
                            y_intersect = t*(line_nodes[(line+1)%N_line_nodes,1]-line_nodes[line,1]) + line_nodes[line,1]
                            
                            # compute the following element where the line lies
                            #print "    e%d :: (%.2f, %.2f) --> (%.2f, %.2f)" % (element, x0_intersect, y0_intersect, x_intersect, y_intersect)        
                            #pylab.scatter([x_intersect], [y_intersect], marker="o",color="r", s=50)
                            
                            # compute the following element
                            # since the edge 0 is the one where intersection occurs
                            # the next element is the neighbour 3 of the current element
                            element = myMesh.elements.neighbors[element,3]
                            intersection_edge = edge
                    
                # edge 1
                edge = 1
                
                if edge != not_edge:
                    # compute parameter value at which intersection occurs
                    t = (myMesh.nodes.coordinates[element_nodes[edge],0] - line_nodes[line,0]) \
                            /(line_nodes[(line+1)%N_line_nodes,0] - line_nodes[line,0])
                    
                    if t >= 0:
                        y_temp_intersect = t*(line_nodes[(line+1)%N_line_nodes,1]-line_nodes[line,1]) + line_nodes[line,1]
                        
                        # check if in fact is an intersection
                        if myMesh.nodes.coordinates[element_nodes[edge],1] <= y_temp_intersect <= \
                            myMesh.nodes.coordinates[element_nodes[(edge+1)%4],1]:
                            
                            x_intersect = t*(line_nodes[(line+1)%N_line_nodes,0]-line_nodes[line,0]) + line_nodes[line,0]
                            y_intersect = y_temp_intersect
                            
                            #print "    e%d :: (%.2f, %.2f) --> (%.2f, %.2f)" % (element, x0_intersect, y0_intersect, x_intersect, y_intersect)
                            #pylab.scatter([x_intersect], [y_intersect], marker="o",color="r", s=50)
                            
                            # compute the following element
                            # since the edge 0 is the one where intersection occurs
                            # the next element is the neighbour 3 of the current element
                            element = myMesh.elements.neighbors[element,0]
                            intersection_edge = edge
                    
                # edge 2
                edge = 2
                
                if edge != not_edge:
                    # compute parameter value at which intersection occurs
                    t = (myMesh.nodes.coordinates[element_nodes[edge],1] - line_nodes[line,1]) \
                            /(line_nodes[(line+1)%N_line_nodes,1] - line_nodes[line,1])
                            
                    if t >= 0:
                        x_temp_intersect = t*(line_nodes[(line+1)%N_line_nodes,0]-line_nodes[line,0]) + line_nodes[line,0]
                        
                        # check if in fact is an intersection
                        if myMesh.nodes.coordinates[element_nodes[(edge+1)%4],0] <= x_temp_intersect <= \
                            myMesh.nodes.coordinates[element_nodes[edge],0]:
                            
                            x_intersect = x_temp_intersect
                            y_intersect = t*(line_nodes[(line+1)%N_line_nodes,1]-line_nodes[line,1]) + line_nodes[line,1]
                            
                            #print "    e%d :: (%.2f, %.2f) --> (%.2f, %.2f)" % (element, x0_intersect, y0_intersect, x_intersect, y_intersect)            
                            #pylab.scatter([x_intersect], [y_intersect], marker="o",color="r", s=50)
                            
                            # compute the following element
                            # since the edge 0 is the one where intersection occurs
                            # the next element is the neighbour 3 of the current element
                            element = myMesh.elements.neighbors[element,1]
                            intersection_edge = edge
                    
                # edge 3
                edge = 3
                
                if edge != not_edge:
                    # compute parameter value at which intersection occurs
                    t = (myMesh.nodes.coordinates[element_nodes[edge],0] - line_nodes[line,0]) \
                            /(line_nodes[(line+1)%N_line_nodes,0] - line_nodes[line,0])
                    
                    if t >= 0:
                        y_temp_intersect = t*(line_nodes[(line+1)%N_line_nodes,1]-line_nodes[line,1]) + line_nodes[line,1]
                        
                        # check if in fact is an intersection
                        if myMesh.nodes.coordinates[element_nodes[(edge+1)%4],1] <= y_temp_intersect <= \
                            myMesh.nodes.coordinates[element_nodes[edge],1]:
                            
                            x_intersect = t*(line_nodes[(line+1)%N_line_nodes,0]-line_nodes[line,0]) + line_nodes[line,0]
                            y_intersect = y_temp_intersect
                            
                            #print "    e%d :: (%.2f, %.2f) --> (%.2f, %.2f)" % (element, x0_intersect, y0_intersect, x_intersect, y_intersect)            
                            #pylab.scatter([x_intersect], [y_intersect], marker="o",color="r", s=50)
                            
                            # compute the following element
                            # since the edge 0 is the one where intersection occurs
                            # the next element is the neighbour 3 of the current element
                            element = myMesh.elements.neighbors[element,2]
                            intersection_edge = edge
                
                not_edge = (intersection_edge + 2)%4
                x0_intersect = x_intersect
                y0_intersect = y_intersect
                
                # store start coodinates of intersection
                end_coordinates_intersect[intersects_counter,:] = numpy.array([x_intersect,\
                                                                           y_intersect])
                
            else:
                x_intersect = line_nodes[(line+1)%N_line_nodes,0]
                y_intersect = line_nodes[(line+1)%N_line_nodes,1]
                #print "    e%d :: (%.2f, %.2f) --> (%.2f, %.2f)" % (element, x0_intersect, y0_intersect, x_intersect, y_intersect)
                #pylab.scatter([x_intersect], [y_intersect], marker="o",color="r", s=50)
                
                # store start coodinates of intersection
                end_coordinates_intersect[intersects_counter,:] = numpy.array([x_intersect,\
                                                                           y_intersect])
                # get out of the while loop since the end of the line has been
                # reached
                break
            
            intersects_counter += 1
            
        try:
            all_elements_intersected = numpy.concatenate([all_elements_intersected, elements_intersected[0:(intersects_counter+1)]])
        except NameError:
            all_elements_intersected = elements_intersected[0:(intersects_counter+1)]
            
        try:
            all_start_coordinates_intersect = numpy.concatenate([all_start_coordinates_intersect, start_coordinates_intersect[0:(intersects_counter+1),:]])
        except NameError:
            all_start_coordinates_intersect = start_coordinates_intersect[0:(intersects_counter+1),:]
        
        try:
            all_end_coordinates_intersect = numpy.concatenate([all_end_coordinates_intersect, end_coordinates_intersect[0:(intersects_counter+1),:]])
        except NameError:
            all_end_coordinates_intersect = end_coordinates_intersect[0:(intersects_counter+1),:]
    
    return all_elements_intersected, all_start_coordinates_intersect, all_end_coordinates_intersect