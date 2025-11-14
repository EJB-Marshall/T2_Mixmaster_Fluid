module hdf5_IO_mod
    use real_type_mod, only: dp
    use hdf5
    implicit none

    integer(hid_t) :: file_id, group_id1, group_id2, group_id3
    integer(hid_t) :: dset_ids(13), dspace_id1, dspace_id2
    integer :: nvars, error
    integer :: nt_written = 0 ! Number of timesteps saved to file

    contains 

    !####################################################
    !# Create HDF5 File, Groups etc.
    !####################################################
    subroutine initialise_hdf(filename, Nx, buffer_size, coordinates_x,K)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: Nx, buffer_size
        real(dp), dimension(:), intent(in) :: coordinates_x
        real(dp), intent(in) :: K
        integer(hsize_t), dimension(2) :: dims, maxdims, chunk_dims
        integer(hsize_t), dimension(1) :: dims_time, maxdims_time, chunk_dims_time, grid_dims
        integer(hid_t) :: plist_id, plist_id2, grid_dspace_id, grid_dset_id, K_dspace_id, K_dset_id  ! Property list for chunking
        character(len=120) :: folder


        nt_written = 0 ! Initialise variable

        ! First we initialise the Fortran interface
        call h5open_f(error)

        ! Next, we create the HDF5 file
        folder = '/Users/elliotmarshall/Desktop/T2_Fluid_Fortran/T2_Primitive_Vars/HDF_Files/'

        call h5fcreate_f(trim(folder)//trim(filename), H5F_ACC_TRUNC_F, file_id, error)

        ! Dimension properties for dataspace of dataset
        dims = (/ Nx, 0 /)
        maxdims = (/ int(Nx,8), H5S_UNLIMITED_F /)
        chunk_dims = (/ Nx, buffer_size /)

        ! Create data space for each evolution variable data set
        call h5screate_simple_f(2, dims, dspace_id1, error, maxdims)

        ! Create property list for chunking (required for unlimited dimensions)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
        call h5pset_chunk_f(plist_id, 2, chunk_dims, error)
        call h5pset_deflate_f(plist_id, 6, error) ! Gzip compression

        ! Now, we will create two groups; "Gravitational" and "Fluid"
        ! First create the "Gravitational" group:    
        call h5gcreate_f(file_id, "Gravitational", group_id1, error)

        ! Create datasets inside the "Fluid" group - Notice that the first argument is the *group* id!
        call h5dcreate_f(group_id1, "Sigma_Minus", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(1), error, plist_id)
        call h5dcreate_f(group_id1, "Sigma_Times", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(2), error, plist_id)
        call h5dcreate_f(group_id1, "N_Minus", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(3), error, plist_id)
        call h5dcreate_f(group_id1, "N_Times", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(4), error, plist_id)
        call h5dcreate_f(group_id1, "E11", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(5), error, plist_id)
        call h5dcreate_f(group_id1, "Sigma2", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(6), error, plist_id)
        call h5dcreate_f(group_id1, "Sigma3", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(7), error, plist_id)
        call h5dcreate_f(group_id1, "LambdaTilde", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(8), error, plist_id)

        ! Close the group.
        call h5gclose_f(group_id1, error)

        !Now, repeat for the "Fluid" group
        call h5gcreate_f(file_id, "Fluid", group_id2, error)        
        call h5dcreate_f(group_id2, "tau_log", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(9), error, plist_id)
        call h5dcreate_f(group_id2, "nu1", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(10), error, plist_id)
        call h5dcreate_f(group_id2, "nu2", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(11), error, plist_id)
        call h5dcreate_f(group_id2, "nu3", H5T_NATIVE_DOUBLE, dspace_id1, dset_ids(12), error, plist_id)
        call h5gclose_f(group_id2, error)

        ! Close the 2D dataspace and property list
        call h5sclose_f(dspace_id1, error)
        call h5pclose_f(plist_id, error)


        !For the time data, we can just create a 1D dataset by itself
        dims_time = (/ 0 /)
        maxdims_time = (/ H5S_UNLIMITED_F /)
        chunk_dims_time = (/ buffer_size /)

        ! Create data space for the time data
        call h5screate_simple_f(1, dims_time, dspace_id2, error, maxdims_time)

        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id2, error)
        call h5pset_chunk_f(plist_id2, 1, chunk_dims_time, error)
        call h5pset_deflate_f(plist_id2, 6, error) 

        ! Create the time dataset (Note this is *not* in a group)
        call h5dcreate_f(file_id, "Time", H5T_NATIVE_DOUBLE, dspace_id2, dset_ids(13), error, plist_id2)

        ! Close the dataspace
        call h5sclose_f(dspace_id2, error)

        ! Close property list
        call h5pclose_f(plist_id2, error)

        ! --------------------------------------
        ! Store Grid and Sound Speed Values
        ! --------------------------------------
        ! Create dataspace for grid
        grid_dims = (/ int(Nx,8) /)
        call h5screate_simple_f(1, grid_dims, grid_dspace_id, error)
        
        ! Create grid dataset (no chunking needed - it's fixed size)
        call h5dcreate_f(file_id, "x_coordinates", H5T_NATIVE_DOUBLE, grid_dspace_id, &
                        grid_dset_id, error)
        
        ! Write grid coordinates
        call h5dwrite_f(grid_dset_id, H5T_NATIVE_DOUBLE, coordinates_x(3:Nx+2), grid_dims, error)
        
        ! Close grid resources
        call h5dclose_f(grid_dset_id, error)
        call h5sclose_f(grid_dspace_id, error)

        ! Create dataspace for K
        call h5screate_simple_f(1, (/ 1_8 /), K_dspace_id, error)
        
        ! Create grid dataset (no chunking needed - it's fixed size)
        call h5dcreate_f(file_id, "K", H5T_NATIVE_DOUBLE, K_dspace_id, &
                        K_dset_id, error)
        
        ! Write grid coordinates
        call h5dwrite_f(K_dset_id, H5T_NATIVE_DOUBLE, K, (/ 1_8 /), error)
        
        ! Close datasets
        call h5dclose_f(K_dset_id, error)
        call h5sclose_f(K_dspace_id, error)
        
    end subroutine initialise_hdf

    !####################################################
    !# Write solution buffer to HDF file
    !####################################################
    subroutine write_to_hdf(soln_buffer,time_buffer,buffer_count,Nx)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: soln_buffer
        real(dp), dimension(:), intent(in) :: time_buffer
        integer, intent(in) :: buffer_count, Nx

        ! Local variables
        integer :: i, hstat
        integer(hsize_t), dimension(2) :: new_dims, start_2d, count_2d
        integer(hsize_t), dimension(1) :: new_dims1, start_1d, count_1d
        integer(hid_t) :: memspace_id, filespace_id

        ! -----------------------------
        ! Write Solution Data
        ! -----------------------------

        new_dims  = (/ int(Nx,8), int(nt_written + buffer_count,8) /) ! Append new data after the last set of buffered timesteps
        start_2d  = (/ 0_8, int(nt_written,8) /)         
        count_2d  = (/ int(Nx,8), int(buffer_count,8) /)
        

        ! New total size after appending this buffer
        new_dims = (/ int(Nx,8), int(nt_written + buffer_count,8) /)

        ! Starting position for this write
        start_2d = (/ 0_8, int(nt_written,8) /)

        ! Size of this write
        count_2d = (/ int(Nx,8), int(buffer_count,8) /)

        ! Memory space must match buffer shape
        call h5screate_simple_f(2, count_2d, memspace_id, error)

        do i = 1, 12
            ! Extend dataset
            call h5dset_extent_f(dset_ids(i), new_dims, error)

            ! Select region in file
            call h5dget_space_f(dset_ids(i), filespace_id, error)
            call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, start_2d, count_2d, error)

            ! Write the whole (Nx, buffer_count) block
            call h5dwrite_f(dset_ids(i), H5T_NATIVE_DOUBLE, soln_buffer(:, i, 1:buffer_count), &
                            count_2d, error, memspace_id, filespace_id)

            call h5sclose_f(filespace_id, error)
        end do

        call h5sclose_f(memspace_id, error)
        
        ! -----------------------------
        ! Write Time Data
        ! -----------------------------
        new_dims1 = (/ int(nt_written + buffer_count,8) /)
        start_1d  = (/ int(nt_written,8) /)
        count_1d  = (/ int(buffer_count,8) /)

        call h5dset_extent_f(dset_ids(13), new_dims1, hstat)
        call h5dget_space_f(dset_ids(13), filespace_id, hstat)
        call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, start_1d, count_1d, hstat)

        call h5screate_simple_f(1, count_1d, memspace_id, hstat)
        call h5dwrite_f(dset_ids(13), H5T_NATIVE_DOUBLE, time_buffer(1:buffer_count), &
                        count_1d, hstat, mem_space_id=memspace_id, file_space_id=filespace_id)

        call h5sclose_f(memspace_id, hstat)
        call h5sclose_f(filespace_id, hstat)


        ! Advance the write pointer
        nt_written = nt_written + buffer_count

    end subroutine write_to_hdf


    !####################################################
    !# Close HDF File
    !####################################################
    subroutine close_hdf()
        implicit none
        integer :: i
    
        ! Close all datasets
        do i = 1, 13
            call h5dclose_f(dset_ids(i), error)
        end do
        
        ! Close the file
        call h5fclose_f(file_id, error)
        
        ! Close the Fortran interface
        call h5close_f(error)

    end subroutine close_hdf


end module hdf5_IO_mod