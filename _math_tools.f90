!File _math_tools.f90
! Nereides Least Square Spectral Element Flow Solver
! math tools module fortran routines

subroutine polyeval_2d(x, mx, nx, y, my, ny, polyX, mpolyX, polyY, mpolyY, result, mresult, nresult)
     ! evaluate the polynomial in (x,y)

    implicit none
    ! comments needed
    real*8 x(mx, nx)

    ! comments needed
    real*8 y(my, ny)
    
    ! comments needed
    real*8 polyX(mpolyX)
    real*8 polyY(mpolyY)

    ! comments needed
    real*8 result(mresult, nresult)

    ! comments needed
    integer mx, nx

    ! comments needed
    integer my, ny

    ! comments needed
    integer mpolyX, mpolyY
    
    ! comments needed
    integer mresult, nresult

    ! dummy variables to loop over the points
    integer i,j,n

    integer px, py

    real*8 tempx, tempy

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: x,y,polyX,polyY
    !f2py intent(in,out) :: result
    !f2py intent(hide) :: mx,nx,my,ny,mpolyX,mpolyY,mresult,nresult

    ! compute the order of the polynomials
    px = mpolyX
    py = mpolyY

    ! the 2d polynomial is the product of the two 1d polynomials in x and y
    ! loop over all (x_i, y_i) points and compute the product
    
    do j=1, mx
        do i=1, nx
            tempx = polyX(1)
            tempy = polyY(1)
            do n=2, px
               tempx = tempx * x(i,j) + polyX(n)
            end do
            do n=2, py
               tempy = tempy * y(i,j) + polyY(n)
            end do
            result(i,j) = tempx * tempy
        end do  
    end do
end


subroutine polyeval_all_2d(x, mx, nx, y, my, ny, polyX, mpolyX, npolyX, polyY, &
                           mpolyY, npolyY, result, mresult, nresult, iresult, jresult)
     ! evaluate the polynomial in (x,y)

    implicit none
    ! comments needed
    real*8 x(mx, nx)

    ! comments needed
    real*8 y(my, ny)
    
    ! comments needed
    real*8 polyX(mpolyX, npolyX)
    real*8 polyY(mpolyY, npolyY)

    ! comments needed
    real*8 result(mresult, nresult, iresult, jresult)

    ! comments needed
    integer mx, nx

    ! comments needed
    integer my, ny

    ! comments needed
    integer mpolyX, npolyX, mpolyY, npolyY
    
    ! comments needed
    integer mresult, nresult, iresult, jresult

    ! dummy variables to loop over the points
    integer i,j,n,phi_i,psi_j

    integer px, py

    real*8 tempx, tempy

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: x,y,polyX,polyY
    !f2py intent(in,out) :: result
    !f2py intent(hide) :: mx,nx,my,ny,mpolyX,npolyX,mpolyY,npolyY,mresult,nresult,iresult,jresult

    ! compute the order of the polynomials
    px = npolyX
    py = npolyY

    ! the 2d polynomial is the product of the two 1d polynomials in x and y
    ! loop over all (x_i, y_i) points and compute the product
    if (nx == 1) then
        do phi_i=1, mpolyX
            do psi_j=1, mpolyY
                do i=1, mx
                    tempx = polyX(phi_i,1)
                    tempy = polyY(psi_j,1)
                    do n=2, px
                        tempx = tempx * x(i,1) + polyX(phi_i,n)
                    end do
                    do n=2, py
                        tempy = tempy * y(i,1) + polyY(psi_j,n)
                    end do
                    result(phi_i,psi_j,i,1) = tempx * tempy
                end do
            end do
        end do
    else
        do phi_i=1, mpolyX
            do psi_j=1, mpolyY
                do j=1, nx
                    do i=1, mx
                        tempx = polyX(phi_i,1)
                        tempy = polyY(psi_j,1)
                        do n=2, px
                            tempx = tempx * x(i,j) + polyX(phi_i,n)
                        end do
                        do n=2, py
                            tempy = tempy * y(i,j) + polyY(psi_j,n)
                        end do
                        result(phi_i,psi_j,i,j) = tempx * tempy
                    end do  
                end do
            end do
        end do
    endif
end


subroutine inner_product(f_0, mf_0, nf_0, f_1, mf_1, nf_1, jacobian, mjacobian, njacobian, weights, mweights, result)
    ! compute the gauss lobatto integration of the product of f_0 with f_1
    ! basically computes the double sum:
    !
    ! \sum_{ij} f_0(u_i,v_j)*f_1(u_i,v_j)*w_i*w_j

    implicit none
    ! comments needed
    real*8 f_0(mf_0, nf_0)

    ! comments needed
    real*8 f_1(mf_1, nf_1)

    ! comments needed
    real*8 jacobian(mjacobian, njacobian)
    
    ! comments needed
    real*8 weights(mweights)

    ! result variable
    real*8 result

    ! comments needed
    integer mf_0, nf_0

    ! comments needed
    integer mf_1, nf_1
    
    ! comments needed
    integer mjacobian, njacobian

    ! comments needed
    integer mweights
    
    ! dummy variables to loop over the points
    integer i,j

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: f_0,f_1,jacobian,weights
    !f2py intent(in,out) :: result
    !f2py intent(hide) :: mf_0, nf_0, mf_1, nf_1, mjacobian, njacobian, mweights

    ! the 2d polynomial is the product of the two 1d polynomials in x and y
    ! loop over all (x_i, y_i) points and compute the product

    if (mjacobian == 1) then
        do j=1, mf_0
            do i=1, nf_0
                result = result + f_0(i,j)*f_1(i,j)*jacobian(1,1)*weights(i)*weights(j)
            end do  
        end do
    else
        do j=1, mf_0
            do i=1, nf_0
                result = result + f_0(i,j)*f_1(i,j)*jacobian(i,j)*weights(i)*weights(j)
            end do  
        end do
    endif
end


subroutine all_basis_inner_products(varphi_ij, ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij, &
                              du_varphi_ij, idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij, &
                              dv_varphi_ij, idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij, &
                              du_t_x, udu_t_x, vdu_t_x, &
                              dv_t_X, udv_t_X, vdv_t_X, &
                              du_t_y, udu_t_y, vdu_t_y, &
                              dv_t_y, udv_t_y, vdv_t_y, &
                              jacobian_0, ujacobian_0, vjacobian_0, &
                              jacobian_1, ujacobian_1, vjacobian_1, &
                              jacobian_2, ujacobian_2, vjacobian_2, &
                              weights_u, uweights_u, weights_v, vweights_v, &
                              result, dresult, iresult, jresult, nresult, mresult)
    ! compute the gauss lobatto integration of all the basis functions between themselves
    ! basically computes the double sum:
    !
    ! \sum_{ij} f_0(u_i,v_j)*f_1(u_i,v_j)*w_i*w_j

    implicit none
    ! comments needed
    real*8 varphi_ij(ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij)

    ! comments needed
    real*8 du_varphi_ij(idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij)

    ! comments needed
    real*8 dv_varphi_ij(idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij)
    
    real*8 du_t_x(udu_t_x, vdu_t_x)
    real*8 dv_t_X(udv_t_X, vdv_t_X)
    real*8 du_t_y(udu_t_y, vdu_t_y)
    real*8 dv_t_y(udv_t_y, vdv_t_y)

    ! comments needed
    real*8 jacobian_0(ujacobian_0, vjacobian_0)
    real*8 jacobian_1(ujacobian_1, vjacobian_1)
    real*8 jacobian_2(ujacobian_2, vjacobian_2)

    real*8 scalling_0, scalling_1, scalling_2, weight_2d
    real*8 dx_varphi_ij_uv, dy_varphi_ij_uv

    ! comments needed
    real*8 weights_u(uweights_u), weights_v(vweights_v)

    ! comments needed
    real*8 result(dresult, iresult, jresult, nresult, mresult)

    ! comments needed
    integer ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij

    ! comments needed
    integer idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij
    
    ! comments needed
    integer idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij

    integer udu_t_x, vdu_t_x, udv_t_X, vdv_t_X, udu_t_y, vdu_t_y, udv_t_y, vdv_t_y

    ! comments needed
    integer ujacobian_0, vjacobian_0
    integer ujacobian_1, vjacobian_1
    integer ujacobian_2, vjacobian_2
    
    integer uweights_u, vweights_v

    integer dresult, iresult, jresult, nresult, mresult
    
    ! dummy variables to loop over the points
    integer i,j,n,m,u,v

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: varphi_ij, du_varphi_ij, dv_varphi_ij
    !f2py intent(in) :: du_t_x, dv_t_X, du_t_y, dv_t_y
    !f2py intent(in) :: jacobian_0, jacobian_1, jacobian_2, weights_u, weights_v
    !f2py intent(in,out) :: result
    !f2py intent(hide) :: ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij
    !f2py intent(hide) :: idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_i
    !f2py intent(hide) :: idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij
    !f2py intent(hide) :: udu_t_x, vdu_t_x, udv_t_X, vdv_t_X, udu_t_y, vdu_t_y, udv_t_y, vdv_t_y
    !f2py intent(hide) :: ujacobian_0, vjacobian_0
    !f2py intent(hide) :: ujacobian_1, vjacobian_1
    !f2py intent(hide) :: ujacobian_2, vjacobian_2
    !f2py intent(hide) :: uweights_u, vweights_v
    !f2py intent(hide) :: dresult, iresult, jresult, nresult, mresult


    ! the 2d polynomial is the product of the two 1d polynomials in x and y
    ! loop over all (x_i, y_i) points and compute the product
    do v=1, vvarphi_ij
        do u=1, uvarphi_ij
            ! this values are computed here to avoid their repetitive computation
            ! when there is no need in doing it, 2x speed up
            weight_2d = weights_u(u)*weights_v(v)
            scalling_0 = jacobian_0(u,v)*weight_2d
            scalling_1 = jacobian_1(u,v)*weight_2d
            scalling_2 = jacobian_2(u,v)*weight_2d
            do i=1, ivarphi_ij
                do j=1, jvarphi_ij
                    dx_varphi_ij_uv = (dv_t_y(u,v)*du_varphi_ij(i,j,u,v) - du_t_y(u,v)*dv_varphi_ij(i,j,u,v))
                    dy_varphi_ij_uv = (du_t_x(u,v)*dv_varphi_ij(i,j,u,v) - dv_t_x(u,v)*du_varphi_ij(i,j,u,v))
                    do n=1, ivarphi_ij
                        do m=1, jvarphi_ij 
                            
                            result(1,i,j,n,m) = result(1,i,j,n,m) + varphi_ij(i,j,u,v)*varphi_ij(n,m,u,v)* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_0(u,v)*weights(u)*weights(v)
                                                 scalling_0
                            result(2,i,j,n,m) = result(2,i,j,n,m) + dx_varphi_ij_uv*varphi_ij(n,m,u,v)* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_1(u,v)*weights(u)*weights(v)
                                                scalling_1
                            result(3,i,j,n,m) = result(3,i,j,n,m) + dy_varphi_ij_uv*varphi_ij(n,m,u,v)* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_1(u,v)*weights(u)*weights(v)
                                                 scalling_1
                            result(4,i,j,n,m) = result(4,i,j,n,m) + dx_varphi_ij_uv*(dv_t_y(u,v)*du_varphi_ij(n,m,u,v) - &
                                                du_t_y(u,v)*dv_varphi_ij(n,m,u,v))* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_2(u,v)*weights(u)*weights(v)
                                                 scalling_2
                            result(5,i,j,n,m) = result(5,i,j,n,m) + dy_varphi_ij_uv*(du_t_x(u,v)*dv_varphi_ij(n,m,u,v) - &
                                                dv_t_x(u,v)*du_varphi_ij(n,m,u,v))* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_2(u,v)*weights(u)*weights(v)
                                                 scalling_2
                            result(6,i,j,n,m) = result(6,i,j,n,m) + dx_varphi_ij_uv*(du_t_x(u,v)*dv_varphi_ij(n,m,u,v) - &
                                                dv_t_x(u,v)*du_varphi_ij(n,m,u,v))* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_2(u,v)*weights(u)*weights(v)
                                                 scalling_2
                        end do
                    end do
                end do
            end do
        end do
    end do
end


subroutine all_other_basis_inner_products(varphi_ij, ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij, &
                              du_varphi_ij, idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij, &
                              dv_varphi_ij, idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij, &
                              other_varphi_ij, iother_varphi_ij, jother_varphi_ij, uother_varphi_ij, vother_varphi_ij, &
                              du_other_varphi_ij, idu_other_varphi_ij, jdu_other_varphi_ij, udu_other_varphi_ij, &
                              vdu_other_varphi_ij, &
                              dv_other_varphi_ij, idv_other_varphi_ij, jdv_other_varphi_ij, udv_other_varphi_ij, &
                              vdv_other_varphi_ij, &
                              du_t_x, udu_t_x, vdu_t_x, &
                              dv_t_X, udv_t_X, vdv_t_X, &
                              du_t_y, udu_t_y, vdu_t_y, &
                              dv_t_y, udv_t_y, vdv_t_y, &
                              jacobian, ujacobian, vjacobian, &
                              weights_u, uweights_u, weights_v, vweights_v, type_du_dv, &
                              result, iresult, jresult, nresult, mresult)
    ! compute the gauss lobatto integration of all the basis functions between themselves
    ! basically computes the double sum:
    !
    ! \sum_{ij} f_0(u_i,v_j)*f_1(u_i,v_j)*w_i*w_j

    implicit none
    ! comments needed
    real*8 varphi_ij(ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij)

    ! comments needed
    real*8 du_varphi_ij(idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij)

    ! comments needed
    real*8 dv_varphi_ij(idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij)

    ! comments needed
    real*8 other_varphi_ij(iother_varphi_ij, jother_varphi_ij, uother_varphi_ij, vother_varphi_ij)

    ! comments needed
    real*8 du_other_varphi_ij(idu_other_varphi_ij, jdu_other_varphi_ij, udu_other_varphi_ij, vdu_other_varphi_ij)

    ! comments needed
    real*8 dv_other_varphi_ij(idv_other_varphi_ij, jdv_other_varphi_ij, udv_other_varphi_ij, vdv_other_varphi_ij)

    
    real*8 du_t_x(udu_t_x, vdu_t_x)
    real*8 dv_t_X(udv_t_X, vdv_t_X)
    real*8 du_t_y(udu_t_y, vdu_t_y)
    real*8 dv_t_y(udv_t_y, vdv_t_y)

    ! comments needed
    real*8 jacobian(ujacobian, vjacobian)

    real*8 scalling_0, scalling_1, scalling_2, weight_2d
    real*8 dx_varphi_ij_uv, dy_varphi_ij_uv, dx_other_varphi_ij_uv, dy_other_varphi_ij_uv

    ! comments needed
    real*8 weights_u(uweights_u), weights_v(vweights_v)
    
    ! comments needed
    integer type_du_dv
    
    ! comments needed
    real*8 result(iresult, jresult, nresult, mresult)

    ! comments needed
    integer ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij

    ! comments needed
    integer idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij
    
    ! comments needed
    integer idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij

    ! comments needed
    integer iother_varphi_ij, jother_varphi_ij, uother_varphi_ij, vother_varphi_ij

    ! comments needed
    integer idu_other_varphi_ij, jdu_other_varphi_ij, udu_other_varphi_ij, vdu_other_varphi_ij
    
    ! comments needed
    integer idv_other_varphi_ij, jdv_other_varphi_ij, udv_other_varphi_ij, vdv_other_varphi_ij

    integer udu_t_x, vdu_t_x, udv_t_X, vdv_t_X, udu_t_y, vdu_t_y, udv_t_y, vdv_t_y

    ! comments needed
    integer ujacobian, vjacobian
    
    integer uweights_u, vweights_v

    integer iresult, jresult, nresult, mresult
    
    ! dummy variables to loop over the points
    integer i,j,n,m,u,v

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: varphi_ij, du_varphi_ij, dv_varphi_ij
    !f2py intent(in) :: other_varphi_ij, du_other_varphi_ij, dv_other_varphi_ij
    !f2py intent(in) :: du_t_x, dv_t_X, du_t_y, dv_t_y
    !f2py intent(in) :: jacobian, weights_u, weights_v
    !f2py intent(in) :: type_du_dv
    !f2py intent(in,out) :: result
    !f2py intent(hide) :: ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij
    !f2py intent(hide) :: idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_i
    !f2py intent(hide) :: idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij
    
    !f2py intent(hide) :: iother_varphi_ij, jother_varphi_ij, uother_varphi_ij, vother_varphi_ij
    !f2py intent(hide) :: idu_other_varphi_ij, jdu_other_varphi_ij, udu_other_varphi_ij, vdu_other_varphi_i
    !f2py intent(hide) :: idv_other_varphi_ij, jdv_other_varphi_ij, udv_other_varphi_ij, vdv_other_arphi_ij
    !f2py intent(hide) :: udu_t_x, vdu_t_x, udv_t_X, vdv_t_X, udu_t_y, vdu_t_y, udv_t_y, vdv_t_y
    !f2py intent(hide) :: ujacobian, vjacobian
    !f2py intent(hide) :: uweights_u, vweights_v
    !f2py intent(hide) :: iresult, jresult, nresult, mresult


    ! this could be done in one group of loops only with the if's inside but that would be slower. Hence
    ! we have chosen to put explicitly all the 9 possible cases and the specific loops
    
    ! the case du_dv = [[0,0], [0,0]] the inner product is <phi_ij, psi_nm> psi_nm is other_varphi_ij
    if (type_du_dv == 0) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_0 = jacobian(u,v)*weight_2d
               do i=1, ivarphi_ij
                   do j=1, jvarphi_ij
                       do n=1, iother_varphi_ij
                           do m=1, jother_varphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + varphi_ij(i,j,u,v)*other_varphi_ij(n,m,u,v)* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_0(u,v)*weights(u)*weights(v)
                                                    scalling_0
                           end do
                       end do
                   end do
               end do
           end do
       end do
    ! the case du_dv = [[1,0], [0,0]] the inner product is <dx_phi_ij, psi_nm>
    else if (type_du_dv == 1) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_1 = jacobian(u,v)*weight_2d
               do i=1, ivarphi_ij
                   do j=1, jvarphi_ij
                       dx_varphi_ij_uv = (dv_t_y(u,v)*du_varphi_ij(i,j,u,v) - du_t_y(u,v)*dv_varphi_ij(i,j,u,v))
                       do n=1, iother_varphi_ij
                           do m=1, jother_varphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + dx_varphi_ij_uv*other_varphi_ij(n,m,u,v)* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_1(u,v)*weights(u)*weights(v)
                                                   scalling_1
                           end do
                       end do
                   end do
               end do
           end do
       end do
    ! the case du_dv = [[0,1], [0,0]] the inner product is <dy_phi_ij, psi_nm>
    else if (type_du_dv == 2) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_1 = jacobian(u,v)*weight_2d
               do i=1, ivarphi_ij
                   do j=1, jvarphi_ij
                       dy_varphi_ij_uv = (du_t_x(u,v)*dv_varphi_ij(i,j,u,v) - dv_t_x(u,v)*du_varphi_ij(i,j,u,v))
                       do n=1, iother_varphi_ij
                           do m=1, jother_varphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + dy_varphi_ij_uv*other_varphi_ij(n,m,u,v)* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_1(u,v)*weights(u)*weights(v)
                                                    scalling_1
                           end do
                       end do
                   end do
               end do
           end do
       end do
    ! the case du_dv = [[0,0], [1,0]] the inner product is <phi_ij, dx_psi_nm>
    else if (type_du_dv == 4) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_1 = jacobian(u,v)*weight_2d
               do n=1, iother_varphi_ij
                   do m=1, jother_varphi_ij
                       dx_other_varphi_ij_uv = (dv_t_y(u,v)*du_other_varphi_ij(n,m,u,v) - du_t_y(u,v)*dv_other_varphi_ij(n,m,u,v))
                       do i=1, ivarphi_ij
                           do j=1, jvarphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + dx_other_varphi_ij_uv*varphi_ij(i,j,u,v)* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_1(u,v)*weights(u)*weights(v)
                                                   scalling_1
                           end do
                       end do
                   end do
               end do
           end do
       end do
    ! the case du_dv = [[1,0], [1,0]] the inner product is <dx_phi_ij, dx_psi_nm>
    else if (type_du_dv == 5) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_2 = jacobian(u,v)*weight_2d
               do i=1, ivarphi_ij
                   do j=1, jvarphi_ij
                       dx_varphi_ij_uv = (dv_t_y(u,v)*du_varphi_ij(i,j,u,v) - du_t_y(u,v)*dv_varphi_ij(i,j,u,v))
                       do n=1, iother_varphi_ij
                           do m=1, jother_varphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + dx_varphi_ij_uv*(dv_t_y(u,v)*du_other_varphi_ij(n,m,u,v) - &
                                                   du_t_y(u,v)*dv_other_varphi_ij(n,m,u,v))* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_2(u,v)*weights(u)*weights(v)
                                                    scalling_2
                           end do
                       end do
                   end do
               end do
           end do
       end do
    ! the case du_dv = [[0,1], [1,0]] the inner product is <dy_phi_ij, dx_psi_nm>
    else if (type_du_dv == 6) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_2 = jacobian(u,v)*weight_2d
               do i=1, ivarphi_ij
                   do j=1, jvarphi_ij
                       dy_varphi_ij_uv = (du_t_x(u,v)*dv_varphi_ij(i,j,u,v) - dv_t_x(u,v)*du_varphi_ij(i,j,u,v))
                       do n=1, iother_varphi_ij
                           do m=1, jother_varphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + dy_varphi_ij_uv*(du_t_x(u,v)*du_other_varphi_ij(n,m,u,v) - &
                                                   dv_t_x(u,v)*dv_other_varphi_ij(n,m,u,v))* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_2(u,v)*weights(u)*weights(v)
                                                    scalling_2
                           end do
                       end do
                   end do
               end do
           end do
       end do
    ! the case du_dv = [[0,0], [0,1]] the inner product is <phi_ij, dy_psi_nm>
    else if (type_du_dv == 8) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_1 = jacobian(u,v)*weight_2d
               do n=1, iother_varphi_ij
                   do m=1, jother_varphi_ij
                       dy_other_varphi_ij_uv = (du_t_x(u,v)*dv_other_varphi_ij(n,m,u,v) - dv_t_x(u,v)*du_other_varphi_ij(n,m,u,v))
                       do i=1, ivarphi_ij
                           do j=1, jvarphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + dy_other_varphi_ij_uv*varphi_ij(i,j,u,v)* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_1(u,v)*weights(u)*weights(v)
                                                    scalling_1
                           end do
                       end do
                   end do
               end do
           end do
       end do    
    ! the case du_dv = [[1,0], [0,1]] the inner product is <dx_phi_ij, dy_psi_nm>
    else if (type_du_dv == 9) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_2 = jacobian(u,v)*weight_2d
               do i=1, ivarphi_ij
                   do j=1, jvarphi_ij
                       dx_varphi_ij_uv = (dv_t_y(u,v)*du_varphi_ij(i,j,u,v) - du_t_y(u,v)*dv_varphi_ij(i,j,u,v))
                       do n=1, iother_varphi_ij
                           do m=1, jother_varphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + dx_varphi_ij_uv*(du_t_x(u,v)*dv_other_varphi_ij(n,m,u,v) - &
                                                   dv_t_x(u,v)*du_other_varphi_ij(n,m,u,v))* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_2(u,v)*weights(u)*weights(v)
                                                    scalling_2
                           end do
                       end do
                   end do
               end do
           end do
       end do
    ! the case du_dv = [[0,1], [0,1]] the inner product is <dy_phi_ij, dy_psi_nm>
    else if (type_du_dv == 10) then
       ! the 2d polynomial is the product of the two 1d polynomials in x and y
       ! loop over all (x_i, y_i) points and compute the product
       do v=1, vvarphi_ij
           do u=1, uvarphi_ij
               ! this values are computed here to avoid their repetitive computation
               ! when there is no need in doing it, 2x speed up
               weight_2d = weights_u(u)*weights_v(v)
               scalling_2 = jacobian(u,v)*weight_2d
               do i=1, ivarphi_ij
                   do j=1, jvarphi_ij
                       dy_varphi_ij_uv = (du_t_x(u,v)*dv_varphi_ij(i,j,u,v) - dv_t_x(u,v)*du_varphi_ij(i,j,u,v))
                       do n=1, ivarphi_ij
                           do m=1, jvarphi_ij 
                               result(i,j,n,m) = result(i,j,n,m) + dy_varphi_ij_uv*(du_t_x(u,v)*dv_varphi_ij(n,m,u,v) - &
                                                   dv_t_x(u,v)*du_varphi_ij(n,m,u,v))* &
                               ! this is the explicit form, the following is the fast one
                               !                    jacobian_2(u,v)*weights(u)*weights(v)
                                                    scalling_2
                           end do
                       end do
                   end do
               end do
           end do
       end do
    end if
end


subroutine all_function_basis_inner_products(evaluated_f, uevaluated_f, vevaluated_f, &
                              varphi_ij, ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij, &
                              du_varphi_ij, idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij, &
                              dv_varphi_ij, idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij, &
                              du_t_x, udu_t_x, vdu_t_x, &
                              dv_t_X, udv_t_X, vdv_t_X, &
                              du_t_y, udu_t_y, vdu_t_y, &
                              dv_t_y, udv_t_y, vdv_t_y, &
                              jacobian_0, ujacobian_0, vjacobian_0, &
                              jacobian_1, ujacobian_1, vjacobian_1, &
                              weights, uweights, &
                              coeffs, ncoeffs, &
                              derivs, nderivs, &
                              result, iresult, jresult)

    ! compute the gauss lobatto integration of all the basis functions between themselves
    ! basically computes the double sum:
    !
    ! \sum_{ij} f_0(u_i,v_j)*f_1(u_i,v_j)*w_i*w_j

    implicit none

    real*8 evaluated_f(uevaluated_f, vevaluated_f)

    ! comments needed
    real*8 varphi_ij(ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij)

    ! comments needed
    real*8 du_varphi_ij(idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij)

    ! comments needed
    real*8 dv_varphi_ij(idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij)
    
    real*8 du_t_x(udu_t_x, vdu_t_x)
    real*8 dv_t_X(udv_t_X, vdv_t_X)
    real*8 du_t_y(udu_t_y, vdu_t_y)
    real*8 dv_t_y(udv_t_y, vdv_t_y)

    ! comments needed
    real*8 jacobian_0(ujacobian_0, vjacobian_0)
    real*8 jacobian_1(ujacobian_1, vjacobian_1)

    real*8 scalling_0, scalling_1, scalling_2, weight_2d
    real*8 dx_varphi_ij_uv, dy_varphi_ij_uv

    ! comments needed
    real*8 weights(uweights)

    real*8 coeffs(ncoeffs)

    integer*2 derivs(nderivs)

    ! comments needed
    real*8 result(iresult, jresult)
    
    integer uevaluated_f, vevaluated_f
    ! comments needed
    integer ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij

    ! comments needed
    integer idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij
    
    ! comments needed
    integer idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij

    integer udu_t_x, vdu_t_x, udv_t_X, vdv_t_X, udu_t_y, vdu_t_y, udv_t_y, vdv_t_y

    ! comments needed
    integer ujacobian_0, vjacobian_0
    integer ujacobian_1, vjacobian_1
    
    integer uweights

    integer ncoeffs
    integer nderivs

    integer iresult, jresult
    
    ! dummy variables to loop over the points
    integer i,j,n,u,v

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: evaluated_f
    !f2py intent(in) :: varphi_ij, du_varphi_ij, dv_varphi_ij
    !f2py intent(in) :: du_t_x, dv_t_X, du_t_y, dv_t_y
    !f2py intent(in) :: jacobian_0, jacobian_1, jacobian_2, weights
    !f2py intent(in) :: coeffs, derivs
    !f2py intent(in,out) :: result
    !f2py intent(hide) :: uevaluated_f, vevaluated_f
    !f2py intent(hide) :: ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij
    !f2py intent(hide) :: idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_i
    !f2py intent(hide) :: idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij
    !f2py intent(hide) :: udu_t_x, vdu_t_x, udv_t_X, vdv_t_X, udu_t_y, vdu_t_y, udv_t_y, vdv_t_y
    !f2py intent(hide) :: ujacobian_0, vjacobian_0
    !f2py intent(hide) :: ujacobian_1, vjacobian_1
    !f2py intent(hide) :: uweights
    !f2py intent(hide) :: ncoeffs, nderivs
    !f2py intent(hide) :: iresult, jresult


    ! the 2d polynomial is the product of the two 1d polynomials in x and y
    ! loop over all (x_i, y_i) points and compute the product
    do v=1, vvarphi_ij
        do u=1, uvarphi_ij
            ! this values are computed here to avoid their repetitive computation
            ! when there is no need in doing it, 2x speed up
            weight_2d = weights(u)*weights(v)
            scalling_0 = jacobian_0(u,v)*weight_2d
            scalling_1 = jacobian_1(u,v)*weight_2d
            do i=1, ivarphi_ij
                do j=1, jvarphi_ij
                    dx_varphi_ij_uv = (dv_t_y(u,v)*du_varphi_ij(i,j,u,v) - du_t_y(u,v)*dv_varphi_ij(i,j,u,v))
                    dy_varphi_ij_uv = (du_t_x(u,v)*dv_varphi_ij(i,j,u,v) - dv_t_x(u,v)*du_varphi_ij(i,j,u,v))
                    do n=1, ncoeffs
                        if (derivs(n) == 0) then
                            result(i,j) = result(i,j) + coeffs(n)*evaluated_f(u,v)*varphi_ij(i,j,u,v)* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_0(u,v)*weights(u)*weights(v)
                                                 scalling_0
                        else if (derivs(n) == 1) then
                            result(i,j) = result(i,j) + coeffs(n)*evaluated_f(u,v)*dx_varphi_ij_uv* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_0(u,v)*weights(u)*weights(v)
                                                 scalling_1
                        else
                            result(i,j) = result(i,j) + coeffs(n)*evaluated_f(u,v)*dy_varphi_ij_uv* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_0(u,v)*weights(u)*weights(v)
                                                 scalling_1
                        end if
                    end do    
                end do
            end do
        end do
    end do
end


subroutine all_function_basis_inner_products_v2(evaluated_f, uevaluated_f, vevaluated_f, &
                              varphi_ij, ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij, &
                              du_varphi_ij, idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij, &
                              dv_varphi_ij, idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij, &
                              du_t_x, udu_t_x, vdu_t_x, &
                              dv_t_X, udv_t_X, vdv_t_X, &
                              du_t_y, udu_t_y, vdu_t_y, &
                              dv_t_y, udv_t_y, vdv_t_y, &
                              jacobian_0, ujacobian_0, vjacobian_0, &
                              jacobian_1, ujacobian_1, vjacobian_1, &
                              weights_x, uweights_x, weights_y, vweights_y, &
                              coeffs, ncoeffs, &
                              derivs, nderivs, &
                              result, iresult, jresult)

    ! compute the gauss lobatto integration of all the basis functions between themselves
    ! basically computes the double sum:
    !
    ! \sum_{ij} f_0(u_i,v_j)*f_1(u_i,v_j)*w_i*w_j

    implicit none

    real*8 evaluated_f(uevaluated_f, vevaluated_f)

    ! comments needed
    real*8 varphi_ij(ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij)

    ! comments needed
    real*8 du_varphi_ij(idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij)

    ! comments needed
    real*8 dv_varphi_ij(idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij)
    
    real*8 du_t_x(udu_t_x, vdu_t_x)
    real*8 dv_t_X(udv_t_X, vdv_t_X)
    real*8 du_t_y(udu_t_y, vdu_t_y)
    real*8 dv_t_y(udv_t_y, vdv_t_y)

    ! comments needed
    real*8 jacobian_0(ujacobian_0, vjacobian_0)
    real*8 jacobian_1(ujacobian_1, vjacobian_1)

    real*8 scalling_0, scalling_1, scalling_2, weight_2d
    real*8 dx_varphi_ij_uv, dy_varphi_ij_uv

    ! comments needed
    real*8 weights_x(uweights_x)
    real*8 weights_y(vweights_y)

    real*8 coeffs(ncoeffs)

    integer*2 derivs(nderivs)

    ! comments needed
    real*8 result(iresult, jresult)
    
    integer uevaluated_f, vevaluated_f
    ! comments needed
    integer ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij

    ! comments needed
    integer idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_ij
    
    ! comments needed
    integer idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij

    integer udu_t_x, vdu_t_x, udv_t_X, vdv_t_X, udu_t_y, vdu_t_y, udv_t_y, vdv_t_y

    ! comments needed
    integer ujacobian_0, vjacobian_0
    integer ujacobian_1, vjacobian_1
    
    integer uweights_x, vweights_y

    integer ncoeffs
    integer nderivs

    integer iresult, jresult
    
    ! dummy variables to loop over the points
    integer i,j,n,u,v

    ! declare to f2py the intentions of the use of the variables
    ! in is as input only
    ! in,out is as input and as output
    ! hide is for the dummy or not so dummy parameters that are passed
    ! automatically by python to the fortran routine
    !f2py intent(in) :: evaluated_f
    !f2py intent(in) :: varphi_ij, du_varphi_ij, dv_varphi_ij
    !f2py intent(in) :: du_t_x, dv_t_X, du_t_y, dv_t_y
    !f2py intent(in) :: jacobian_0, jacobian_1, jacobian_2, weights_x_, weights_y
    !f2py intent(in) :: coeffs, derivs
    !f2py intent(in,out) :: result
    !f2py intent(hide) :: uevaluated_f, vevaluated_f
    !f2py intent(hide) :: ivarphi_ij, jvarphi_ij, uvarphi_ij, vvarphi_ij
    !f2py intent(hide) :: idu_varphi_ij, jdu_varphi_ij, udu_varphi_ij, vdu_varphi_i
    !f2py intent(hide) :: idv_varphi_ij, jdv_varphi_ij, udv_varphi_ij, vdv_varphi_ij
    !f2py intent(hide) :: udu_t_x, vdu_t_x, udv_t_X, vdv_t_X, udu_t_y, vdu_t_y, udv_t_y, vdv_t_y
    !f2py intent(hide) :: ujacobian_0, vjacobian_0
    !f2py intent(hide) :: ujacobian_1, vjacobian_1
    !f2py intent(hide) :: uweights_x, vweights_y
    !f2py intent(hide) :: ncoeffs, nderivs
    !f2py intent(hide) :: iresult, jresult


    ! the 2d polynomial is the product of the two 1d polynomials in x and y
    ! loop over all (x_i, y_i) points and compute the product
    do v=1, vvarphi_ij
        do u=1, uvarphi_ij
            ! this values are computed here to avoid their repetitive computation
            ! when there is no need in doing it, 2x speed up
            weight_2d = weights_x(u)*weights_y(v)
            scalling_0 = jacobian_0(u,v)*weight_2d
            scalling_1 = jacobian_1(u,v)*weight_2d
            do i=1, ivarphi_ij
                do j=1, jvarphi_ij
                    dx_varphi_ij_uv = (dv_t_y(u,v)*du_varphi_ij(i,j,u,v) - du_t_y(u,v)*dv_varphi_ij(i,j,u,v))
                    dy_varphi_ij_uv = (du_t_x(u,v)*dv_varphi_ij(i,j,u,v) - dv_t_x(u,v)*du_varphi_ij(i,j,u,v))
                    do n=1, ncoeffs
                        if (derivs(n) == 0) then
                            result(i,j) = result(i,j) + coeffs(n)*evaluated_f(u,v)*varphi_ij(i,j,u,v)* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_0(u,v)*weights(u)*weights(v)
                                                 scalling_0
                        else if (derivs(n) == 1) then
                            result(i,j) = result(i,j) + coeffs(n)*evaluated_f(u,v)*dx_varphi_ij_uv* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_0(u,v)*weights(u)*weights(v)
                                                 scalling_1
                        else
                            result(i,j) = result(i,j) + coeffs(n)*evaluated_f(u,v)*dy_varphi_ij_uv* &
                            ! this is the explicit form, the following is the fast one
                            !                    jacobian_0(u,v)*weights(u)*weights(v)
                                                 scalling_1
                        end if
                    end do    
                end do
            end do
        end do
    end do
end

