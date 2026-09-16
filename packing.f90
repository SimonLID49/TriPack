module particle_module
    implicit none

    type :: atom
        character(len=2) :: species
        real :: position(3)
    end type atom

    type :: molecule
        real :: com(3)
        real :: nn_distance
        real :: nn_direction(3)
        type(atom), allocatable :: atoms(:)
    end type molecule

contains

    ! IO routines
    ! #####################

    ! Uility
    ! #####################
    subroutine calculate_com(mol)
        type(molecule), intent(inout) :: mol
        integer :: i

        ! Initialize the center of mass
        mol%com = [0.0, 0.0, 0.0]

        ! Calculate the center of mass of the molecule
        do i = 1, size(mol%atoms)
            mol%com = mol%com + mol%atoms(i)%position
        end do
        mol%com = mol%com / real(size(mol%atoms))
    end subroutine calculate_com   

    subroutine rotate_molecule(mol, axis, angle)
        type(molecule), intent(inout) :: mol
        real, intent(in) :: axis(3)
        real, intent(in) :: angle
        real :: angle_rad
        real :: rotation_matrix(3,3)
        real :: c, s, t
        real :: x, y, z

        angle_rad = angle * 3.14159265358979323846 / 180

        ! Normalize the axis
        x = axis(1)
        y = axis(2)
        z = axis(3)
        c = cos(angle_rad)
        s = sin(angle_rad)
        t = 1.0 - c

        call make_rotation_matrix([x, y, z], angle_rad, rotation_matrix)

        call rotate_molecule_bg(mol, rotation_matrix)

    end subroutine rotate_molecule

    subroutine make_rotation_matrix(axis, theta, R)

        real, intent(in)  :: axis(3)
        real, intent(in)  :: theta
        real, intent(out) :: R(3,3)

        real :: ux, uy, uz
        real :: c, s, t

        ux = axis(1)
        uy = axis(2)
        uz = axis(3)

        c = cos(theta)
        s = sin(theta)
        t = 1.0 - c

        R(1,1) = t*ux*ux + c
        R(1,2) = t*ux*uy - s*uz
        R(1,3) = t*ux*uz + s*uy

        R(2,1) = t*ux*uy + s*uz
        R(2,2) = t*uy*uy + c
        R(2,3) = t*uy*uz - s*ux

        R(3,1) = t*ux*uz - s*uy
        R(3,2) = t*uy*uz + s*ux
        R(3,3) = t*uz*uz + c

    end subroutine make_rotation_matrix

    subroutine rotate_molecule_bg(mol, rotation_matrix)
        type(molecule), intent(inout) :: mol
        real, intent(in) :: rotation_matrix(3,3)

        integer :: i
        real :: relative_position(3)

        do i = 1, size(mol%atoms)
            relative_position = mol%atoms(i)%position - mol%com

            mol%atoms(i)%position = mol%com + &
                matmul(rotation_matrix, relative_position)
        end do
    end subroutine rotate_molecule_bg

    subroutine translate_molecule(mol, translation_vector)
        type(molecule), intent(inout) :: mol
        real, intent(in) :: translation_vector(3)
        integer :: i

        ! Translate each atom in the molecule by the translation vector
        do i = 1, size(mol%atoms)
            mol%atoms(i)%position = mol%atoms(i)%position + translation_vector
        end do

        ! Update the center of mass after translation
        call calculate_com(mol)
    end subroutine translate_molecule

    function translate_molecule_pbc(mol, cell_matrix) result(coms_pbc)
        type(molecule), intent(in) :: mol
        real, intent(in) :: cell_matrix(3,3)
        real :: coms_pbc(27,3)

        integer :: i, j, k, idx
        real :: shift(3)

        idx = 0

        do i = -1, 1
            do j = -1, 1
                do k = -1, 1

                    idx = idx + 1

                    shift = matmul(cell_matrix, &
                            [real(i), real(j), real(k)])

                    coms_pbc(idx,:) = mol%com + shift

                end do
            end do
        end do

    end function translate_molecule_pbc

    function init_water(com) result(mol)
        real, intent(in) :: com(3)
        type(molecule) :: mol
        ! Initialize a water molecule and shift its center of mass to the desired position
        allocate(mol%atoms(3))
        mol%atoms(1)%species = 'O'
        mol%atoms(1)%position = [0.0, 0.0, 0.0]
        mol%atoms(2)%species = 'H'
        mol%atoms(2)%position = [0.9572, 0.0, 0.0]
        mol%atoms(3)%species = 'H'
        mol%atoms(3)%position = [-0.2399872, 0.927297, 0.0]
        call calculate_com(mol)

        ! Shift the molecule so that its center of mass is at the desired position
        call translate_molecule(mol, com - mol%com)
    end function init_water

        function init_dme(com) result(mol)
        real, intent(in) :: com(3)
        type(molecule) :: mol
        ! Initialize a dimethyl ether molecule and shift its center of mass to the desired position
        allocate(mol%atoms(9))
        mol%atoms(1)%species = 'C'
        mol%atoms(1)%position = [0.0, 0.0, 0.0]
        mol%atoms(2)%species = 'C'
        mol%atoms(2)%position = [1.414000, 0.000000, 0.000000]
        mol%atoms(3)%species = 'O'
        mol%atoms(3)%position = [-1.414000, 0.000000, 0.000000]
        mol%atoms(4)%species = 'H'
        mol%atoms(4)%position = [1.775000, 0.513000, 0.890000]
        mol%atoms(5)%species = 'H'
        mol%atoms(5)%position = [1.775000, 0.513000, -0.890000]
        mol%atoms(6)%species = 'H'
        mol%atoms(6)%position = [1.775000, -1.026000, 0.000000]
        mol%atoms(7)%species = 'H'
        mol%atoms(7)%position = [-1.775000, 0.513000, 0.890000]
        mol%atoms(8)%species = 'H'
        mol%atoms(8)%position = [-1.775000, 0.513000, -0.890000]
        mol%atoms(9)%species = 'H'
        mol%atoms(9)%position = [-1.775000, -1.026000, 0.000000]
        call calculate_com(mol)

        ! Shift the molecule so that its center of mass is at the desired position
        call translate_molecule(mol, com - mol%com)
    end function init_dme

    function init_h(com) result(mol)
        real, intent(in) :: com(3)
        type(molecule) :: mol
        ! Initialize a hydrogen atom and shift its position to the desired position
        allocate(mol%atoms(1))
        mol%atoms(1)%species = 'H'
        mol%atoms(1)%position = com
        call calculate_com(mol)
    end function init_h

end module particle_module

module system_module
    use particle_module
    implicit none

        type :: bin
        real, allocatable :: atom_positions(:,:)  ! (3, n_atoms_in_bin)
        real :: bin_matrix(3, 3)                  ! edges of this bin: a/na, b/nb, c/nc
        real :: origin(3)                         ! Cartesian corner of the bin
        integer :: n_atoms = 0
        integer :: idx3(3)
    end type bin

    type :: system
        type(molecule), allocatable :: molecules(:)
        type(bin), allocatable :: bins(:)
        real :: cell_tensor(9)
        real :: cell_matrix(3,3)
        real :: safety
        real :: binning_width
        logical :: rdm_rotation, confine
        real, allocatable :: grid_points(:, :)
        integer :: n_bins(3) = [1, 1, 1]          ! bins along a, b, c
    end type system

    type :: bin_neighbor
        integer :: lin_idx      ! linear index into sys%bins
        real    :: shift(3)     ! add this to the neighbor bin's atom
                                 ! positions before computing distances
    end type bin_neighbor
contains

    function cross_product(u, v) result(w)
        real, intent(in) :: u(3), v(3)
        real :: w(3)
        w(1) = u(2)*v(3) - u(3)*v(2)
        w(2) = u(3)*v(1) - u(1)*v(3)
        w(3) = u(1)*v(2) - u(2)*v(1)
    end function cross_product

    function mat3_inverse(m) result(inv)
        real, intent(in) :: m(3,3)
        real :: inv(3,3)
        real :: det

        det =  m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) &
             - m(1,2)*(m(2,1)*m(3,3) - m(2,3)*m(3,1)) &
             + m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))

        if (abs(det) <= tiny(det)) error stop 'Singular cell matrix: cannot invert'

        inv(1,1) =  (m(2,2)*m(3,3) - m(2,3)*m(3,2)) / det
        inv(2,1) = -(m(2,1)*m(3,3) - m(2,3)*m(3,1)) / det
        inv(3,1) =  (m(2,1)*m(3,2) - m(2,2)*m(3,1)) / det
        inv(1,2) = -(m(1,2)*m(3,3) - m(1,3)*m(3,2)) / det
        inv(2,2) =  (m(1,1)*m(3,3) - m(1,3)*m(3,1)) / det
        inv(3,2) = -(m(1,1)*m(3,2) - m(1,2)*m(3,1)) / det
        inv(1,3) =  (m(1,2)*m(2,3) - m(1,3)*m(2,2)) / det
        inv(2,3) = -(m(1,1)*m(2,3) - m(1,3)*m(2,1)) / det
        inv(3,3) =  (m(1,1)*m(2,2) - m(1,2)*m(2,1)) / det
    end function mat3_inverse

    subroutine init_binning(sys)
        ! Sub-divide the (possibly triclinic) cell into bins for faster
        ! overlap checking, keeping all cell angles unchanged: each bin is
        ! a scaled-down copy of the full cell, spanned by a/na, b/nb, c/nc.
        !
        ! Because bins are cut along fractional coordinates, the true
        ! (perpendicular) width of a bin along a given lattice direction is
        ! NOT |a|/na -- it is the face-to-face spacing of the full cell,
        ! V / |b x c|, divided by na. sys%binning_width is treated as the
        ! minimum acceptable perpendicular width.

        type(system), intent(inout) :: sys

        real :: a(3), b(3), c(3)
        real :: volume, h(3)
        integer :: na, nb, nc
        integer :: i, j, k, idx, total_bins

        a = sys%cell_matrix(:,1)
        b = sys%cell_matrix(:,2)
        c = sys%cell_matrix(:,3)

        volume = abs(dot_product(a, cross_product(b, c)))
        if (volume <= tiny(volume)) error stop 'Degenerate cell matrix in init_binning'

        ! Perpendicular face-to-face spacing along each lattice direction
        h(1) = volume / norm2(cross_product(b, c))
        h(2) = volume / norm2(cross_product(a, c))
        h(3) = volume / norm2(cross_product(a, b))

        ! Largest number of equal-width bins per direction such that every
        ! bin's perpendicular width is still >= sys%binning_width
        na = max(1, floor(h(1) / sys%binning_width))
        nb = max(1, floor(h(2) / sys%binning_width))
        nc = max(1, floor(h(3) / sys%binning_width))

        sys%n_bins = [na, nb, nc]

        total_bins = na * nb * nc
        if (allocated(sys%bins)) deallocate(sys%bins)
        allocate(sys%bins(total_bins))

        idx = 0
        do i = 0, na - 1
            do j = 0, nb - 1
                do k = 0, nc - 1
                    idx = idx + 1

                    sys%bins(idx)%bin_matrix(:,1) = a / real(na)
                    sys%bins(idx)%bin_matrix(:,2) = b / real(nb)
                    sys%bins(idx)%bin_matrix(:,3) = c / real(nc)

                    sys%bins(idx)%origin = (real(i)/real(na))*a + &
                                            (real(j)/real(nb))*b + &
                                            (real(k)/real(nc))*c

                    sys%bins(idx)%n_atoms = 0

                    sys%bins(idx)%idx3 = [i, j, k]
                end do
            end do
        end do

        print *, 'init_binning: bins per direction (a,b,c) --> ', na, nb, nc
        print *, 'init_binning: total number of bins        --> ', total_bins
        print *, 'init_binning: perpendicular bin widths    --> ', &
                  h(1)/real(na), h(2)/real(nb), h(3)/real(nc)

    end subroutine init_binning

    subroutine visualize_binning(sys)
        type(system), intent(inout) :: sys
        integer :: i

        ! Visualize the bins by adding a hydrogen atom at each bin's origin
        do i = 1, size(sys%bins)
            call add_molecule(sys, init_h(sys%bins(i)%origin))
        end do
    end subroutine visualize_binning

    function get_bin_from_idx3(sys, idx3) result(bin)
        type(system), intent(in) :: sys
        integer, intent(in) :: idx3(3)
        type(bin) :: bin

        bin = sys%bins(get_bin_linear_index(sys, idx3))
    end function get_bin_from_idx3
    ! Given a Cartesian position, return its (i,j,k) bin index (0-based)
    ! along a, b, c. Wraps into the primary cell assuming PBC.
    function get_bin_index(sys, pos) result(idx3)
        type(system), intent(in) :: sys
        real, intent(in) :: pos(3)
        integer :: idx3(3)
        real :: frac(3)

        frac = matmul(mat3_inverse(sys%cell_matrix), pos)
        frac = frac - floor(frac)

        idx3 = int(frac * real(sys%n_bins))
        idx3 = min(idx3, sys%n_bins - 1)
    end function get_bin_index

    ! Linear index into sys%bins matching the (i outer, j mid, k inner)
    ! loop order used to fill it in init_binning.
    function get_bin_linear_index(sys, idx3) result(lin)
        type(system), intent(in) :: sys
        integer, intent(in) :: idx3(3)
        integer :: lin

        lin = 1 + idx3(3) + sys%n_bins(3)*(idx3(2) + sys%n_bins(2)*idx3(1))
    end function get_bin_linear_index

    subroutine clear_bins(sys)
        type(system), intent(inout) :: sys
        integer :: i

        do i = 1, size(sys%bins)
            sys%bins(i)%n_atoms = 0
            if (allocated(sys%bins(i)%atom_positions)) deallocate(sys%bins(i)%atom_positions)
        end do
    end subroutine clear_bins

    subroutine bin_atoms(sys, mol)
        ! Adds mol's atoms into the existing bins WITHOUT clearing them first.
        ! Call clear_bins(sys) once, up front, before placing any molecules.
        type(system), intent(inout) :: sys
        type(molecule), intent(in) :: mol
        integer :: i, idx3(3), lin_idx, n_old
        real, allocatable :: tmp(:,:)

        do i = 1, size(mol%atoms)
            idx3    = get_bin_index(sys, mol%atoms(i)%position)
            lin_idx = get_bin_linear_index(sys, idx3)

            if (.not. allocated(sys%bins(lin_idx)%atom_positions)) then
                allocate(sys%bins(lin_idx)%atom_positions(3, 1))
                sys%bins(lin_idx)%atom_positions(:, 1) = mol%atoms(i)%position
                sys%bins(lin_idx)%n_atoms = 1
            else
                n_old = sys%bins(lin_idx)%n_atoms

                allocate(tmp(3, n_old + 1))
                tmp(:, 1:n_old)   = sys%bins(lin_idx)%atom_positions
                tmp(:, n_old + 1) = mol%atoms(i)%position

                call move_alloc(tmp, sys%bins(lin_idx)%atom_positions)
                sys%bins(lin_idx)%n_atoms = n_old + 1
            end if
        end do
    end subroutine bin_atoms

    function get_adjacent_bins(sys, current_bin) result(neighbors)
        type(system), intent(in) :: sys
        type(bin), intent(in) :: current_bin
        type(bin_neighbor) :: neighbors(27)

        integer :: i, j, k, n
        integer :: raw(3), wrapped(3)
        real :: shift(3)

        n = 0
        do i = -1, 1
            do j = -1, 1
                do k = -1, 1
                    n = n + 1

                    raw = current_bin%idx3 + [i, j, k]

                    ! modulo (NOT mod) gives the correct floored wrap:
                    ! modulo(-1,3) = 2, whereas mod(-1,3) = -1
                    wrapped = modulo(raw, sys%n_bins)

                    ! Cartesian shift needed whenever raw fell outside
                    ! [0, n_bins-1] along that lattice direction
                    shift = 0.0
                    if (raw(1) < 0) then
                        shift = shift - sys%cell_matrix(:,1)
                    else if (raw(1) >= sys%n_bins(1)) then
                        shift = shift + sys%cell_matrix(:,1)
                    end if
                    if (raw(2) < 0) then
                        shift = shift - sys%cell_matrix(:,2)
                    else if (raw(2) >= sys%n_bins(2)) then
                        shift = shift + sys%cell_matrix(:,2)
                    end if
                    if (raw(3) < 0) then
                        shift = shift - sys%cell_matrix(:,3)
                    else if (raw(3) >= sys%n_bins(3)) then
                        shift = shift + sys%cell_matrix(:,3)
                    end if

                    neighbors(n)%lin_idx = get_bin_linear_index(sys, wrapped)
                    neighbors(n)%shift   = shift
                end do
            end do
        end do

    end function get_adjacent_bins

    subroutine init_system(sys, cell_tensor)
        type(system), intent(out) :: sys
        real, intent(in) :: cell_tensor(9)

        ! Initialize the system with the given cell tensor
        sys%cell_tensor = cell_tensor
        sys%cell_matrix = reshape(cell_tensor, [3, 3])
    end subroutine init_system

    subroutine init_grid(sys, spacing)
        type(system), intent(inout) :: sys
        real, intent(in) :: spacing

        ! Initialize the grid for the system
        real :: a(3), b(3), c(3)
        real :: pos(3)
        real :: fa, fb, fc
        integer :: na, nb, nc
        integer :: i, j, k, n, npoints
        real, allocatable :: grid(:,:)

        a = sys%cell_matrix(:,1)
        b = sys%cell_matrix(:,2)
        c = sys%cell_matrix(:,3)

        ! Number of intervals along each cell vector
        na = max(1, nint(sqrt(sum(a**2)) / spacing))
        nb = max(1, nint(sqrt(sum(b**2)) / spacing))
        nc = max(1, nint(sqrt(sum(c**2)) / spacing))

        ! Number of grid points
        npoints = na * nb * nc

        allocate(sys%grid_points(3,npoints))

        n = 0

        do i = 0, na-1
            do j = 0, nb-1
                do k = 0, nc-1

                    ! Fractional coordinates
                    fa = real(i,8) / real(na,8)
                    fb = real(j,8) / real(nb,8)
                    fc = real(k,8) / real(nc,8)

                    ! Convert fractional -> Cartesian
                    pos = fa*a + fb*b + fc*c

                    n = n + 1
                    sys%grid_points(:,n) = pos

                end do
            end do
        end do

        print *, "Number of grid points -->", npoints

    end subroutine init_grid

    subroutine visualize_grid(sys)
        type(system), intent(inout) :: sys
        integer :: i

        ! Visualie the grid points by adding a hydrogen atom at each grid point
        do i = 1, size(sys%grid_points, 2)
            call add_molecule(sys, init_h(sys%grid_points(:,i)))
        end do
    end subroutine visualize_grid

    subroutine add_molecule(sys, mol)
        type(system), intent(inout) :: sys
        type(molecule), intent(in) :: mol
        type(molecule), allocatable :: molecules(:)
        integer :: n

        ! Add a molecule to the system
        if (.not. allocated(sys%molecules)) then
            allocate(sys%molecules(1))
            sys%molecules(1) = mol
        else
            n = size(sys%molecules)
            allocate(molecules(n + 1))
            molecules(:n) = sys%molecules
            molecules(n + 1) = mol
            call move_alloc(molecules, sys%molecules)
        end if
    end subroutine add_molecule

    subroutine print_system(sys)
        type(system), intent(in) :: sys
        integer :: i, j

        ! Print the details of the system
        do i = 1, size(sys%molecules)
            print *, 'Molecule ', i
            print *, 'Center of Mass: ', sys%molecules(i)%com
            do j = 1, size(sys%molecules(i)%atoms)
                print *, 'Atom ', j, ': ', sys%molecules(i)%atoms(j)%species, &
                         ' Position: ', sys%molecules(i)%atoms(j)%position
            end do
        end do
    end subroutine print_system

    function overlap_sys(sys, mol) result(check)

        type(system), intent(in) :: sys
        type(molecule), intent(in) :: mol

        logical :: check
        integer :: i, s
        real :: shifts(3,27)
        real :: safety2

        ! Precompute the 27 periodic image shifts once (cheap, avoids
        ! recomputing matmul for every molecule)
        call compute_pbc_shifts(sys%cell_matrix, shifts)

        safety2 = sys%safety**2
        check = .false.

        !$omp parallel do default(none) &
        !$omp shared(sys, mol, shifts, safety2) &
        !$omp private(i, s) &
        !$omp reduction(.or.:check)
        do i = 1, size(sys%molecules)
            do s = 1, 27
                if (molecules_overlap_shifted(mol, sys%molecules(i), &
                                            shifts(:,s), safety2)) then
                    check = .true.
                    exit   ! legal here: this "do s" is a plain loop, not the omp worksharing loop
                end if
            end do
        end do
        !$omp end parallel do

    end function overlap_sys

    function check_all_inside_box(sys, mol)
        type(system), intent(in) :: sys
        type(molecule), intent(in) :: mol
        integer :: i
        real :: fractional_coords(3), coords(3)
        real :: inverse_cell(3,3), determinant
        logical :: check_all_inside_box

        determinant = sys%cell_matrix(1,1) * &
            (sys%cell_matrix(2,2)*sys%cell_matrix(3,3) - &
             sys%cell_matrix(2,3)*sys%cell_matrix(3,2)) - &
            sys%cell_matrix(1,2) * &
            (sys%cell_matrix(2,1)*sys%cell_matrix(3,3) - &
             sys%cell_matrix(2,3)*sys%cell_matrix(3,1)) + &
            sys%cell_matrix(1,3) * &
            (sys%cell_matrix(2,1)*sys%cell_matrix(3,2) - &
             sys%cell_matrix(2,2)*sys%cell_matrix(3,1))

        if (abs(determinant) <= tiny(determinant)) then
            check_all_inside_box = .false.
            return
        end if

        inverse_cell(1,1) = (sys%cell_matrix(2,2)*sys%cell_matrix(3,3) - sys%cell_matrix(2,3)*sys%cell_matrix(3,2)) / determinant
        inverse_cell(1,2) = (sys%cell_matrix(1,3)*sys%cell_matrix(3,2) - sys%cell_matrix(1,2)*sys%cell_matrix(3,3)) / determinant
        inverse_cell(1,3) = (sys%cell_matrix(1,2)*sys%cell_matrix(2,3) - sys%cell_matrix(1,3)*sys%cell_matrix(2,2)) / determinant
        inverse_cell(2,1) = (sys%cell_matrix(2,3)*sys%cell_matrix(3,1) - sys%cell_matrix(2,1)*sys%cell_matrix(3,3)) / determinant
        inverse_cell(2,2) = (sys%cell_matrix(1,1)*sys%cell_matrix(3,3) - sys%cell_matrix(1,3)*sys%cell_matrix(3,1)) / determinant
        inverse_cell(2,3) = (sys%cell_matrix(1,3)*sys%cell_matrix(2,1) - sys%cell_matrix(1,1)*sys%cell_matrix(2,3)) / determinant
        inverse_cell(3,1) = (sys%cell_matrix(2,1)*sys%cell_matrix(3,2) - sys%cell_matrix(2,2)*sys%cell_matrix(3,1)) / determinant
        inverse_cell(3,2) = (sys%cell_matrix(1,2)*sys%cell_matrix(3,1) - sys%cell_matrix(1,1)*sys%cell_matrix(3,2)) / determinant
        inverse_cell(3,3) = (sys%cell_matrix(1,1)*sys%cell_matrix(2,2) - sys%cell_matrix(1,2)*sys%cell_matrix(2,1)) / determinant

        do i = 1, size(mol%atoms)
            coords = mol%atoms(i)%position
            fractional_coords = matmul(inverse_cell, coords)
            if (any(fractional_coords < 0.0) .or. any(fractional_coords > 1.0)) then
                check_all_inside_box = .false.
                return
            end if
        end do
        check_all_inside_box = .true.

    end function check_all_inside_box

    subroutine compute_pbc_shifts(cell_matrix, shifts)
        real, intent(in)  :: cell_matrix(3,3)
        real, intent(out) :: shifts(3,27)
        integer :: j, k, l, idx

        idx = 0
        do j = -1, 1
            do k = -1, 1
                do l = -1, 1
                    idx = idx + 1
                    shifts(:,idx) = matmul(cell_matrix, [real(j), real(k), real(l)])
                end do
            end do
        end do
    end subroutine compute_pbc_shifts

    function molecules_overlap_shifted(mol1, mol2, shift, safety2) result(is_overlap)
        ! Returns as soon as ANY atom pair is within sqrt(safety2) of each other.
        ! No sqrt needed, and no need to finish scanning once overlap is found.
        type(molecule), intent(in) :: mol1, mol2
        real, intent(in) :: shift(3)
        real, intent(in) :: safety2
        logical :: is_overlap

        integer :: a_i, a_j
        real :: dr(3)
        real :: distance_squared

        is_overlap = .false. 

        do a_i = 1, size(mol1%atoms)
            do a_j = 1, size(mol2%atoms)

                dr = mol1%atoms(a_i)%position + shift - mol2%atoms(a_j)%position
                distance_squared = sum(dr**2)

                if (distance_squared < safety2) then
                    is_overlap = .true.
                    return   ! legal: ordinary function, exits immediately
                end if

            end do
        end do

    end function molecules_overlap_shifted

    subroutine export_xyz(sys, filename)
        type(system), intent(in) :: sys
        character(len=*), intent(in) :: filename
        integer :: i, j
        integer :: unit
        integer :: total_atoms

        ! Export the system to an XYZ file
        open(newunit=unit, file=filename, status='replace', action='write')

        total_atoms = 0
        do i = 1, size(sys%molecules)
            total_atoms = total_atoms + size(sys%molecules(i)%atoms)
        end do

        write(unit, '(I0)') total_atoms

        write(unit, '(A,9(1X,F0.8),A)') 'Lattice="', &
            sys%cell_matrix(1,1), sys%cell_matrix(2,1), sys%cell_matrix(3,1), &
            sys%cell_matrix(1,2), sys%cell_matrix(2,2), sys%cell_matrix(3,2), &
            sys%cell_matrix(1,3), sys%cell_matrix(2,3), sys%cell_matrix(3,3), &
            '" pbc = "T T T"'

        do i = 1, size(sys%molecules)
            do j = 1, size(sys%molecules(i)%atoms)
                write(unit, '(A,3F13.8)') trim(sys%molecules(i)%atoms(j)%species), &
                     sys%molecules(i)%atoms(j)%position(1), &
                     sys%molecules(i)%atoms(j)%position(2), &
                     sys%molecules(i)%atoms(j)%position(3)
            end do
        end do
        close(unit)
    end subroutine export_xyz

    function overlap_bins(sys, mol) result(check)
        type(system), intent(in) :: sys
        type(molecule), intent(in) :: mol
        logical :: check
        integer :: i, j, s
        real :: shifts(3,27)
        real :: safety2, distance
        integer :: idx3(3), lin_idx
        type(bin_neighbor) :: neighbors(27)

        safety2 = sys%safety**2
        check = .false.

        ! parallelize over atoms in the new molecule
        !$omp parallel do default(none) &
        !$omp shared(sys, mol, shifts, safety2) &
        !$omp private(i, j, s, idx3, lin_idx, neighbors, distance) &
        !$omp reduction(.or.:check)

        do i = 1, size(mol%atoms)
            idx3   = get_bin_index(sys, mol%atoms(i)%position)
            lin_idx = get_bin_linear_index(sys, idx3)

            neighbors = get_adjacent_bins(sys, sys%bins(lin_idx))

            do j = 1, 27
                !print *, "Checking bin ", neighbors(j)%lin_idx, " with size ", sys%bins(neighbors(j)%lin_idx)%n_atoms, " atoms"
                if (sys%bins(neighbors(j)%lin_idx)%n_atoms > 0) then
                    do s = 1, sys%bins(neighbors(j)%lin_idx)%n_atoms
                        distance  = norm2(mol%atoms(i)%position + neighbors(j)%shift - &
                                    sys%bins(neighbors(j)%lin_idx)%atom_positions(:,s))
                        if (distance < sys%safety) then
                            check = .true.
                            return
                        end if
                    end do
                end if
            end do
        end do

        !$omp end parallel do

    end function overlap_bins

    function try_place(sys, mol, rdm_rotation) result(success)
        type(system), intent(inout) :: sys
        type(molecule), intent(inout) :: mol
        logical, intent(in), optional :: rdm_rotation
        real :: rdm_angle, rdm_axis
        logical :: success
        integer :: attempt

        ! Try to place the molecule at every grid point without overlapping existing molecules
        call shuffle_grid_points(sys)

        if (present(rdm_rotation) .and. rdm_rotation) then
            call random_number(rdm_angle)
            rdm_angle = 360.0 * rdm_angle
            ! random axis of rotation x, y, z
            call random_number(rdm_axis)
            rdm_axis = 3.0 * rdm_axis
            if (rdm_axis < 1.0) then
                call rotate_molecule(mol, [1.0, 0.0, 0.0], rdm_angle)
            else if (rdm_axis < 2.0) then
                call rotate_molecule(mol, [0.0, 1.0, 0.0], rdm_angle)
            else
                call rotate_molecule(mol, [0.0, 0.0, 1.0], rdm_angle)
            end if
        end if

        success = .false.
        do attempt = 1, size(sys%grid_points, 2)
            call translate_molecule(mol, sys%grid_points(:,attempt) - mol%com)

            if (sys%confine) then
                if (.not. check_all_inside_box(sys, mol)) cycle
            end if

            if (size(sys%molecules) == 0) then
                call add_molecule(sys, mol)
                call bin_atoms(sys, mol)
                success = .true.
                exit
            end if
            if (.not. overlap_bins(sys, mol)) then
                call add_molecule(sys, mol)
                call bin_atoms(sys, mol)
                success = .true.
                exit
            end if
        end do

    end function try_place

    subroutine shuffle_grid_points(sys)

        type(system), intent(inout) :: sys

        integer :: npoints
        integer :: i, j
        real :: r
        real :: tmp(3)

        npoints = size(sys%grid_points, 2)

        ! Fisher-Yates shuffle
        do i = npoints, 2, -1

            call random_number(r)
            j = 1 + int(r * real(i))

            ! Swap grid points i and j
            tmp = sys%grid_points(:,i)
            sys%grid_points(:,i) = sys%grid_points(:,j)
            sys%grid_points(:,j) = tmp

        end do

    end subroutine shuffle_grid_points

    subroutine initial_guess(sys, mols, n_mols)
        type(system), intent(inout) :: sys
        type(molecule), intent(inout) :: mols(:)
        integer, intent(in) :: n_mols(:)
        integer :: i
        logical :: success
        integer :: j

        ! Try to place n_mol molecules in the system
        ! Evenly distribute the placement attempts across the provided molecules
        do i = 1, size(n_mols)
            do j = 1, n_mols(i)
                success = try_place(sys, mols(i), sys%rdm_rotation)
                if (.not. success) then
                    print *, "Failed to place molecule ", j, " of type ", i
                    exit
                end if
            end do

        end do

    end subroutine initial_guess

end module system_module

module io_module
    use system_module
    implicit none

contains

    subroutine print_logo()
        character(len=512) :: logo(10)
        integer :: i
        
        logo(1) = ' _____    _      _ _       _       ______          _'
        logo(2) = '|_   _|  (_)    | (_)     (_)      | ___ \        | |'
        logo(3) = '  | |_ __ _  ___| |_ _ __  _  ___  | |_/ /_ _  ___| | _____ _ __'
        logo(4) = '  | | ''__| |/ __| | | ''_ \| |/ __| |  __/ _` |/ __| |/ / _ \ ''__|'
        logo(5) = '  | | |  | | (__| | | | | | | (__  | | | (_| | (__|   <  __/ |'
        logo(6) = '  \_/_|  |_|\___|_|_|_| |_|_|\___| \_|  \__,_|\___|_|\_\___|_|'
        logo(7) = ''

        do i = 1, 7
            write(*,'(A)') trim(logo(i))
        end do

    end subroutine print_logo

    function read_xyz(filename) result(mol)
        character(len=*), intent(in) :: filename
        type(molecule) :: mol
        integer :: unit, ios, natoms, i
        character(len=256) :: comment, species
        real :: x, y, z

        ! Read the XYZ file and initialize the molecule
        open(newunit=unit, file=trim(filename), status='old', action='read', &
             iostat=ios)
        if (ios /= 0) error stop 'Could not open XYZ file'

        read(unit, *, iostat=ios) natoms
        if (ios /= 0 .or. natoms < 1) then
            close(unit)
            error stop 'Invalid XYZ atom count'
        end if

        ! The XYZ comment line may contain cell data, which is ignored.
        read(unit, '(A)', iostat=ios) comment
        if (ios /= 0) then
            close(unit)
            error stop 'Invalid XYZ comment line'
        end if

        allocate(mol%atoms(natoms))
        do i = 1, natoms
            read(unit, *, iostat=ios) species, x, y, z
            if (ios /= 0 .or. len_trim(species) == 0) then
                close(unit)
                error stop 'Invalid XYZ atom line'
            end if
            mol%atoms(i)%species = trim(species)
            mol%atoms(i)%position = [x, y, z]
        end do
        close(unit)

        call calculate_com(mol)

    end function read_xyz

    subroutine init_system_from_xyz(sys, filename)
        type(system), intent(inout) :: sys
        type(molecule) :: starting_conf
        character(len=*), intent(in) :: filename
        integer :: unit, ios, natoms
        character(len=4096) :: comment
        real :: cell_tensor(9)
        integer :: lattice_start, lattice_end, i

        ! Read the XYZ file and initialize the system
        open(newunit=unit, file=trim(filename), status='old', action='read', &
             iostat=ios)
        if (ios /= 0) error stop 'Could not open XYZ file'

        read(unit, *, iostat=ios) natoms
        if (ios /= 0 .or. natoms < 1) then
            close(unit)
            error stop 'Invalid XYZ atom count'
        end if

        ! The XYZ comment line may contain cell data.
        read(unit, '(A)', iostat=ios) comment
        if (ios /= 0) then
            close(unit)
            error stop 'Invalid XYZ comment line'
        end if

        ! Extract cell tensor from the comment line if present
        lattice_start = index(comment, 'Lattice="')
        if (lattice_start > 0) then
            lattice_start = lattice_start + len('Lattice="')
            lattice_end = index(comment(lattice_start:), '"')
            if (lattice_end == 0) then
                close(unit)
                error stop 'Invalid Lattice data in XYZ comment line'
            end if
            lattice_end = lattice_start + lattice_end - 2
            read(comment(lattice_start:lattice_end), *, iostat=ios) cell_tensor
            if (ios /= 0) then
                close(unit)
                error stop 'Invalid Lattice data in XYZ comment line'
            end if
            call init_system(sys, cell_tensor)
        else
            error stop 'No lattice information found in XYZ comment line'
        end if

        print *, 'Cell tensor from XYZ:'
        do i = 1, 3
            print '(3F12.4)', sys%cell_matrix(:, i)
        end do

        close(unit)

        ! now read the molecules from the XYZ file and add them to the system

        starting_conf = read_xyz(filename)
        call add_molecule(sys, starting_conf)
        

    end subroutine init_system_from_xyz

    function input_reader() result(sys)
        type(system) :: sys
        real :: cell_tensor(9)
        real :: grid_spacing, binning_width
        integer :: unit, ios, n_mols, i, at_i, at_j
        character(len=256) :: first_line, second_line, third_line, fourth_line, fifth_line
        character(len=256) :: keyword, xyz_filename
        character(len=256), allocatable :: mol_filenames(:)
        integer, allocatable :: mol_counts(:)
        type(molecule), allocatable :: mols(:)

        ! Provide a defined result until the input records are parsed.
        sys%cell_tensor = 0.0
        sys%cell_matrix = 0.0

        ! A bit of output art:
        call print_logo()

        open(newunit=unit, file='simon.inp', status='old', action='read', &
            iostat=ios)

        if (ios /= 0) error stop 'Could not open input file'

        ! Read first line from simon.inp
        read(unit, '(A)', iostat=ios) first_line
        if (ios /= 0) error stop 'Could not read first line'

        read(first_line, *) keyword

        if (trim(keyword) == 'cell') then
            read(first_line, *) keyword, cell_tensor
            print *, 'cell keyword found --> '
            call init_system(sys, cell_tensor)
            do i = 1, 3
                print '(3F12.4)', sys%cell_matrix(:, i)
            end do

        else if (trim(keyword) == 'sysfile') then
            read(first_line, *) keyword, xyz_filename
            print *, 'sysfile keyword found --> ', trim(xyz_filename)
            call init_system_from_xyz(sys, trim(xyz_filename))
        else
            error stop 'Unknown keyword on first line'
        end if

        ! Read second line from simon.inp
        read(unit, '(A)', iostat=ios) second_line
        if (ios /= 0) error stop 'Could not read second line'

        read(second_line, *) keyword

        if (trim(keyword) == 'gridspacing') then
            read(second_line, *) keyword, grid_spacing
            call init_grid(sys, grid_spacing)
            print *, 'gridspacing keyword found -->', grid_spacing

        else
            error stop 'Unknown keyword on second line'
        end if

        ! Read third line from simon.inp
        read(unit, '(A)', iostat=ios) third_line
        if (ios /= 0) error stop 'Could not read third line'

        read(third_line, *) keyword

        if (trim(keyword) == 'safety') then
            read(third_line, *) keyword, sys%safety
            print *, 'safety keyword found --> ', sys%safety
        else
            error stop 'Unknown keyword on third line'
        end if

        ! Read fourth line from input with logical rotation keyword
        read(unit, '(A)', iostat=ios) fourth_line
        if (ios /= 0) error stop 'Could not read fourth line'
        read(fourth_line, *) keyword
        if (trim(keyword) == 'rotate') then
            read(fourth_line, *) keyword, sys%rdm_rotation
            print *, 'rotation keyword found --> ', sys%rdm_rotation
        else
            print *, keyword
            error stop 'Unknown keyword on fourth line'
            
        end if

        ! Read confine information from fifth line:
        read(unit, '(A)', iostat=ios) fifth_line
        if (ios /= 0) error stop 'Could not read fifth line'
        read(fifth_line, *) keyword
        if (trim(keyword) == 'confine') then
            read(fifth_line, *) keyword, sys%confine
            print *, 'confine keyword found --> ', sys%confine
        else
            error stop 'Unknown keyword on fifth line'
        end if

        ! Read molecule information from the input file

        print *, 'molecule keyword found -->'
        read(unit, *) keyword, n_mols
        allocate(mol_filenames(n_mols))
        allocate(mol_counts(n_mols))
        read(unit, *) (mol_filenames(i), i = 1, n_mols)
        read(unit, *) (mol_counts(i), i = 1, n_mols)

        close(unit)

        ! Read each molecule from its XYZ file and add it to mols array
        allocate(mols(n_mols))
        binning_width = sys%safety
        do i = 1, n_mols
            mols(i) = read_xyz(trim(mol_filenames(i)))
            print *, 'Read molecule from file: '
            print *, trim(mol_filenames(i)), &
                     ' with ', size(mols(i)%atoms), ' atoms and count ', mol_counts(i)
        end do

        sys%binning_width = 2.0
        call init_binning(sys)

        call initial_guess(sys, mols, mol_counts)
        call export_xyz(sys, 'output.xyz')
    end function input_reader

end module io_module


program test_packing
    use io_module
    implicit none

    type(molecule) :: mol1, mol2, mol3
    type(molecule) :: mols(2)
    type(bin_neighbor) :: adjacent_bins(27)
    real :: coms(27, 3)
    type(system) :: sys
    integer :: n_mols(2), i
    real :: rdm_pos(3), atom_pos(3)
    logical :: check

    sys = input_reader()

end program test_packing