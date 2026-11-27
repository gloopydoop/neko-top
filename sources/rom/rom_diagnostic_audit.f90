!> Instrumentation for a reproducible ROM/MMA diagnostic audit.
!!
!! This module is deliberately separate from the reduced equations.  It only
!! observes vectors at the MMA input/output boundaries and evaluates ROM-only
!! finite differences.  In particular, it never changes the supplied gradient,
!! the MMA update, the objective weights, or candidate acceptance.
module rom_diagnostic_audit
  use rom_problem, only: rom_problem_t
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use mma, only: mma_t
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use vector_math, only: vector_copy, vector_sub3, vector_absval, vector_glsc2, &
       vector_glmax
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE
  use neko_config, only: NEKO_BCKND_DEVICE
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION, pe_rank, pe_size
  use mpi_f08, only: MPI_Allreduce, MPI_Barrier, MPI_Gather, MPI_Gatherv, &
       MPI_INTEGER, MPI_MIN
  use utils, only: mkdir, neko_error
  implicit none
  private

  public :: rom_diagnostic_audit_t

  !> File-backed ROM diagnostic collector.
  type :: rom_diagnostic_audit_t
     logical :: enabled = .false.
     logical :: initialized = .false.
     character(len=:), allocatable :: directory
     character(len=:), allocatable :: vector_directory
     integer :: n_local = 0
     integer :: n_constraints = 0
     integer :: rank_count = 0
     integer :: trial_rank = 0
     integer :: test_rank = 0
     integer :: enriched_rank = 0
     integer :: n_steps = 0
     real(kind=rp) :: dt = 0.0_rp
     character(len=:), allocatable :: formulation
     character(len=:), allocatable :: gradient_source
     integer :: outer_max_iterations = 0
     integer :: inner_iterations = 0
     real(kind=rp) :: solver_tolerance = 0.0_rp
     real(kind=rp) :: validation_tolerance = 0.0_rp
     real(kind=rp) :: validation_gradient_tolerance = 0.0_rp
     real(kind=rp) :: mma_scale = 1.0_rp
     logical :: mma_auto_scale = .false.
     real(kind=rp) :: mma_asyinit = 0.2_rp
     real(kind=rp) :: mma_asyincr = 1.2_rp
     real(kind=rp) :: mma_asydecr = 0.7_rp
     real(kind=rp) :: mma_move_limit = 0.2_rp
     real(kind=rp) :: mma_xmin = 0.0_rp
     real(kind=rp) :: mma_xmax = 1.0_rp
     integer :: mma_max_iter = 100
     character(len=:), allocatable :: mma_backend_subsolver
     integer :: iteration_unit = -1
     integer :: gradient_unit = -1
     integer :: fd_unit = -1
     integer :: objective_unit = -1
   contains
     procedure, pass(this) :: init => rom_audit_init
     procedure, pass(this) :: free => rom_audit_free
     procedure, pass(this) :: write_evaluation => rom_audit_write_evaluation
     procedure, pass(this) :: write_full_evaluation => rom_audit_write_full_evaluation
     procedure, pass(this) :: write_step => rom_audit_write_step
     procedure, pass(this) :: write_gradient_comparison => &
          rom_audit_write_gradient_comparison
     procedure, pass(this) :: write_objective_comparison => &
          rom_audit_write_objective_comparison
     procedure, pass(this) :: finite_difference => rom_audit_finite_difference
  end type rom_diagnostic_audit_t

contains

  !> Initialise the audit directory, metadata and CSV streams.
  subroutine rom_audit_init(this, directory, rom, mma, outer_max_iterations, &
       inner_iterations, solver_tolerance, validation_tolerance, &
       validation_gradient_tolerance, mma_scale, mma_auto_scale, mma_asyinit, &
       mma_asyincr, mma_asydecr, mma_move_limit, mma_xmin, mma_xmax)
    class(rom_diagnostic_audit_t), intent(inout) :: this
    character(len=*), intent(in) :: directory
    type(rom_problem_t), intent(in) :: rom
    type(mma_t), intent(in) :: mma
    integer, intent(in) :: outer_max_iterations, inner_iterations
    real(kind=rp), intent(in) :: solver_tolerance, validation_tolerance
    real(kind=rp), intent(in) :: validation_gradient_tolerance, mma_scale
    logical, intent(in) :: mma_auto_scale
    real(kind=rp), intent(in) :: mma_asyinit, mma_asyincr, mma_asydecr
    real(kind=rp), intent(in) :: mma_move_limit, mma_xmin, mma_xmax
    integer :: ios, unit
    character(len=:), allocatable :: source_revision, source_dirty

    call this%free()
    this%enabled = .true.
    this%initialized = .true.
    this%directory = trim(directory)
    this%vector_directory = trim(directory)//'/vectors'
    this%n_local = rom%weights%size()
    this%n_constraints = 1
    this%rank_count = pe_size
    this%trial_rank = rom%basis%n_modes
    this%test_rank = rom%test_basis%n_modes
    this%enriched_rank = rom%enriched_adjoint_rank
    this%n_steps = rom%n_steps
    this%dt = rom%dt
    this%formulation = trim(rom%formulation)
    this%gradient_source = trim(rom%sensitivity_source)
    this%outer_max_iterations = outer_max_iterations
    this%inner_iterations = inner_iterations
    this%solver_tolerance = solver_tolerance
    this%validation_tolerance = validation_tolerance
    this%validation_gradient_tolerance = validation_gradient_tolerance
    this%mma_scale = mma_scale
    this%mma_auto_scale = mma_auto_scale
    this%mma_asyinit = mma_asyinit
    this%mma_asyincr = mma_asyincr
    this%mma_asydecr = mma_asydecr
    this%mma_move_limit = mma_move_limit
    this%mma_xmin = mma_xmin
    this%mma_xmax = mma_xmax
    this%mma_max_iter = mma%get_max_iter()
    this%mma_backend_subsolver = mma%get_backend_and_subsolver()

    if (pe_rank == 0) then
       call mkdir(trim(this%directory))
       call mkdir(trim(this%vector_directory))
    end if
    call MPI_Barrier(NEKO_COMM)
    call rom_audit_environment('ROM_AUDIT_SOURCE_REVISION', source_revision)
    call rom_audit_environment('ROM_AUDIT_SOURCE_DIRTY', source_dirty)

    if (pe_rank == 0) then
       open(newunit=unit, file=trim(this%directory)//'/configuration.json', &
            status='replace', action='write', iostat=ios)
       if (ios /= 0) call neko_error('could not open ROM audit configuration.json')
       write(unit, '(A)') '{'
       write(unit, '(A,A,A)') '  "source_revision": "', trim(source_revision), '",'
       write(unit, '(A,A,A)') '  "source_worktree": "', trim(source_dirty), '",'
       write(unit, '(A,I0,A)') '  "mpi_ranks": ', pe_size, ','
       write(unit, '(A,I0,A)') '  "local_design_size": ', this%n_local, ','
       write(unit, '(A,A,A)') '  "formulation": "', trim(rom%formulation), '",'
       write(unit, '(A,A,A)') '  "sensitivity_source": "', &
            trim(rom%sensitivity_source), '",'
       write(unit, '(A,I0,A)') '  "trial_basis_rank": ', rom%basis%n_modes, ','
       write(unit, '(A,I0,A)') '  "test_basis_rank": ', rom%test_basis%n_modes, ','
       write(unit, '(A,I0,A)') '  "enriched_adjoint_rank": ', &
            rom%enriched_adjoint_rank, ','
       write(unit, '(A,ES24.16,A)') '  "time_start": ', rom%time_start, ','
       write(unit, '(A,ES24.16,A)') '  "time_end": ', rom%time_end, ','
       write(unit, '(A,ES24.16,A)') '  "dt": ', rom%dt, ','
       write(unit, '(A,I0,A)') '  "n_steps": ', rom%n_steps, ','
       write(unit, '(A,I0,A)') '  "substeps": ', rom%substeps, ','
       write(unit, '(A,ES24.16,A)') '  "scalar_metric_weight": ', &
            rom%scalar_metric_weight, ','
       write(unit, '(A,ES24.16,A)') '  "scalar_mixing_weight": ', &
            rom%scalar_mixing_weight, ','
       write(unit, '(A,ES24.16,A)') '  "viscous_weight": ', rom%viscous_weight, ','
       write(unit, '(A,ES24.16,A)') '  "brinkman_weight": ', rom%brinkman_weight, ','
       write(unit, '(A,ES24.16,A)') '  "scalar_mixing_start": ', &
            rom%scalar_mixing_start, ','
       write(unit, '(A,ES24.16,A)') '  "viscous_start": ', rom%viscous_start, ','
       write(unit, '(A,ES24.16,A)') '  "brinkman_start": ', rom%brinkman_start, ','
       write(unit, '(A,I0,A)') '  "outer_max_iterations": ', &
            this%outer_max_iterations, ','
       write(unit, '(A,I0,A)') '  "inner_iterations": ', this%inner_iterations, ','
       write(unit, '(A,ES24.16,A)') '  "solver_tolerance": ', &
            this%solver_tolerance, ','
       write(unit, '(A,ES24.16,A)') '  "validation_tolerance": ', &
            this%validation_tolerance, ','
       write(unit, '(A,ES24.16,A)') '  "validation_gradient_tolerance": ', &
            this%validation_gradient_tolerance, ','
       write(unit, '(A,ES24.16,A)') '  "mma_scale": ', this%mma_scale, ','
       if (this%mma_auto_scale) then
          write(unit, '(A)') '  "mma_auto_scale": true,'
       else
          write(unit, '(A)') '  "mma_auto_scale": false,'
       end if
       write(unit, '(A,ES24.16,A)') '  "mma_asyinit": ', this%mma_asyinit, ','
       write(unit, '(A,ES24.16,A)') '  "mma_asyincr": ', this%mma_asyincr, ','
       write(unit, '(A,ES24.16,A)') '  "mma_asydecr": ', this%mma_asydecr, ','
       write(unit, '(A,ES24.16,A)') '  "mma_move_limit": ', this%mma_move_limit, ','
       write(unit, '(A,ES24.16,A)') '  "mma_xmin": ', this%mma_xmin, ','
       write(unit, '(A,ES24.16,A)') '  "mma_xmax": ', this%mma_xmax, ','
       write(unit, '(A,I0,A)') '  "mma_max_iter": ', this%mma_max_iter, ','
       write(unit, '(A,A,A)') '  "mma_backend_subsolver": "', &
            trim(this%mma_backend_subsolver), '",'
       write(unit, '(A,A)') '  "gradient_vector_convention": ', &
            '"selected ROM gradient is copied to MMA unchanged",'
       write(unit, '(A,A)') '  "design_vector_convention": ', &
            '"unfiltered design vector; ROM mapping then avg_B projection",'
       write(unit, '(A,A)') '  "mpi_vector_format": ', &
            '"little-endian float64 concatenated in NEKO_COMM rank order",'
       write(unit, '(A,A)') '  "vector_files": ', &
            '"binary files contain gathered local nodal vectors; shape is recorded by local_design_size and mpi_ranks"'
       write(unit, '(A)') '}'
       close(unit)

    end if
    call rom_audit_write_partition(this)

    if (pe_rank == 0) then
       open(newunit=this%iteration_unit, &
            file=trim(this%directory)//'/iteration_diagnostics.csv', &
            status='replace', action='write', iostat=ios)
       if (ios /= 0) call neko_error('could not open ROM audit iteration CSV')
       write(this%iteration_unit, '(A)') &
            'outer_iteration,inner_iteration,event,design_before,design_after,'// &
            'gradient_before,gradient_after,gradient_exact_before,'// &
            'gradient_projected_before,gradient_exact_after,'// &
            'gradient_projected_after,raw_density_before,raw_density_after,'// &
            'constraint_gradient_before,constraint_gradient_after,'// &
            'rom_objective_before,rom_objective_after,rom_objective_change,'// &
            'rom_scalar_before,rom_scalar_after,rom_viscous_before,'// &
            'rom_viscous_after,rom_brinkman_before,rom_brinkman_after,'// &
            'constraint_before,constraint_after,update_max,update_l2,'// &
            'g_exact_dot_update,g_projected_dot_update,g_mma_dot_update,'// &
            'kkt_max,kkt_l2,gradient_source,formulation,trial_rank,test_rank,'// &
            'enriched_rank,dt,n_steps,selected_gradient_device_error'
       close(this%iteration_unit)
       this%iteration_unit = -1

       open(newunit=this%gradient_unit, &
            file=trim(this%directory)//'/gradient_comparison.csv', &
            status='replace', action='write', iostat=ios)
       if (ios /= 0) call neko_error('could not open ROM audit gradient CSV')
       write(this%gradient_unit, '(A)') &
            'point,outer_iteration,inner_iteration,norm_exact,norm_projected,'// &
            'norm_full,relative_exact_full,relative_projected_full,'// &
            'relative_exact_projected,cosine_exact_full,cosine_projected_full,'// &
            'cosine_exact_projected,dot_exact_projected,dot_exact_full,'// &
            'dot_projected_full'
       close(this%gradient_unit)
       this%gradient_unit = -1

       open(newunit=this%fd_unit, &
            file=trim(this%directory)//'/directional_fd.csv', &
            status='replace', action='write', iostat=ios)
       if (ios /= 0) call neko_error('could not open ROM audit finite-difference CSV')
       write(this%fd_unit, '(A)') &
            'outer_iteration,direction,epsilon,feasible,objective_plus,'// &
            'objective_minus,fd_derivative,exact_dot_direction,'// &
            'projected_dot_direction,fd_minus_exact,fd_minus_projected,'// &
            'exact_localB_dot_direction,projected_localB_dot_direction,'// &
            'fd_minus_exact_localB,fd_minus_projected_localB,'// &
            'direction_file,base_design_file'
       close(this%fd_unit)
       this%fd_unit = -1

       open(newunit=this%objective_unit, &
            file=trim(this%directory)//'/objective_comparison.csv', &
            status='replace', action='write', iostat=ios)
       if (ios /= 0) call neko_error('could not open ROM audit objective CSV')
       write(this%objective_unit, '(A)') &
            'point,outer_iteration,inner_iteration,design_file,rom_total,'// &
            'rom_scalar_mixing,rom_viscous,rom_brinkman,full_order_total,'// &
            'full_scalar_mixing,full_viscous,full_brinkman,rom_minus_full'
       close(this%objective_unit)
    end if
    call MPI_Barrier(NEKO_COMM)
  end subroutine rom_audit_init

  !> Close streams and release the audit path.
  subroutine rom_audit_free(this)
    class(rom_diagnostic_audit_t), intent(inout) :: this
    if (this%iteration_unit > 0) close(this%iteration_unit)
    if (this%gradient_unit > 0) close(this%gradient_unit)
    if (this%fd_unit > 0) close(this%fd_unit)
    if (this%objective_unit > 0) close(this%objective_unit)
    this%iteration_unit = -1
    this%gradient_unit = -1
    this%fd_unit = -1
    this%objective_unit = -1
    this%initialized = .false.
    this%enabled = .false.
    if (allocated(this%directory)) deallocate(this%directory)
    if (allocated(this%vector_directory)) deallocate(this%vector_directory)
    if (allocated(this%mma_backend_subsolver)) deallocate(this%mma_backend_subsolver)
    if (allocated(this%formulation)) deallocate(this%formulation)
    if (allocated(this%gradient_source)) deallocate(this%gradient_source)
  end subroutine rom_audit_free

  !> Write one reduced evaluation and all its full-dimensional vectors.
  subroutine rom_audit_write_evaluation(this, outer, inner, x, raw_selected, &
       gradient_selected, raw_exact, gradient_exact, raw_projected, &
       gradient_projected, constraints, constraint_gradient)
    class(rom_diagnostic_audit_t), intent(inout) :: this
    integer, intent(in) :: outer, inner
    type(vector_t), intent(inout) :: x, raw_selected, gradient_selected
    type(vector_t), intent(inout) :: raw_exact, gradient_exact
    type(vector_t), intent(inout) :: raw_projected, gradient_projected
    real(kind=rp), intent(in) :: constraints(:)
    type(matrix_t), intent(inout) :: constraint_gradient
    character(len=512) :: base, filename
    character(len=64) :: tag

    if (.not. this%enabled) return
    write(tag, '(A,I0,A,I0)') 'outer_', outer, '_inner_', inner
    base = trim(this%vector_directory)//'/'//trim(tag)
    filename = trim(base)//'_design.bin'
    call rom_audit_write_vector(filename, x)
    filename = trim(base)//'_raw_density_selected.bin'
    call rom_audit_write_vector(filename, raw_selected)
    filename = trim(base)//'_gradient_selected.bin'
    call rom_audit_write_vector(filename, gradient_selected)
    filename = trim(base)//'_raw_density_exact.bin'
    call rom_audit_write_vector(filename, raw_exact)
    filename = trim(base)//'_gradient_exact.bin'
    call rom_audit_write_vector(filename, gradient_exact)
    filename = trim(base)//'_raw_density_projected.bin'
    call rom_audit_write_vector(filename, raw_projected)
    filename = trim(base)//'_gradient_projected.bin'
    call rom_audit_write_vector(filename, gradient_projected)
    filename = trim(base)//'_constraints.bin'
    call rom_audit_write_array(filename, constraints)
    filename = trim(base)//'_constraint_gradient.bin'
    call rom_audit_write_matrix(filename, constraint_gradient)
  end subroutine rom_audit_write_evaluation

  !> Write the full-order raw and projected gradients at a matched design.
  subroutine rom_audit_write_full_evaluation(this, outer, inner, x, raw_gradient, &
       gradient, constraints, constraint_gradient)
    class(rom_diagnostic_audit_t), intent(inout) :: this
    integer, intent(in) :: outer, inner
    type(vector_t), intent(inout) :: x, raw_gradient, gradient
    real(kind=rp), intent(in) :: constraints(:)
    type(matrix_t), intent(inout) :: constraint_gradient
    character(len=512) :: base, filename
    character(len=64) :: tag

    if (.not. this%enabled) return
    write(tag, '(A,I0,A,I0,A)') 'outer_', outer, '_inner_', inner, '_fom'
    base = trim(this%vector_directory)//'/'//trim(tag)
    filename = trim(base)//'_design.bin'
    call rom_audit_write_vector(filename, x)
    filename = trim(base)//'_raw_gradient.bin'
    call rom_audit_write_vector(filename, raw_gradient)
    filename = trim(base)//'_gradient.bin'
    call rom_audit_write_vector(filename, gradient)
    filename = trim(base)//'_constraints.bin'
    call rom_audit_write_array(filename, constraints)
    filename = trim(base)//'_constraint_gradient.bin'
    call rom_audit_write_matrix(filename, constraint_gradient)
  end subroutine rom_audit_write_full_evaluation

  !> Record one MMA transition, using gradients at the pre-update design.
  subroutine rom_audit_write_step(this, outer, inner, x_before, x_after, &
       gradient_exact_before, gradient_projected_before, gradient_mma_before, &
       rom_objective_before, rom_objective_after, rom_scalar_before, &
       rom_scalar_after, rom_viscous_before, rom_viscous_after, &
       rom_brinkman_before, rom_brinkman_after, constraints_before, &
       constraints_after, kkt_max, kkt_l2, selected_gradient_device_error)
    class(rom_diagnostic_audit_t), intent(inout) :: this
    integer, intent(in) :: outer, inner
    type(vector_t), intent(inout) :: x_before, x_after
    type(vector_t), intent(inout) :: gradient_exact_before
    type(vector_t), intent(inout) :: gradient_projected_before
    type(vector_t), intent(inout) :: gradient_mma_before
    real(kind=rp), intent(in) :: rom_objective_before, rom_objective_after
    real(kind=rp), intent(in) :: rom_scalar_before, rom_scalar_after
    real(kind=rp), intent(in) :: rom_viscous_before, rom_viscous_after
    real(kind=rp), intent(in) :: rom_brinkman_before, rom_brinkman_after
    real(kind=rp), intent(in) :: constraints_before(:), constraints_after(:)
    real(kind=rp), intent(in) :: kkt_max, kkt_l2, selected_gradient_device_error
    type(vector_t) :: delta, abs_delta
    real(kind=rp) :: update_max, update_l2, dot_exact, dot_projected, dot_mma
    character(len=512) :: filename
    character(len=64) :: tag, tag_before
    integer :: unit, ios

    if (.not. this%enabled) return
    call delta%init(x_after%size())
    call abs_delta%init(x_after%size())
    call vector_sub3(delta, x_after, x_before)
    update_l2 = sqrt(max(0.0_rp, vector_glsc2(delta, delta)))
    call vector_copy(abs_delta, delta)
    call vector_absval(abs_delta)
    update_max = vector_glmax(abs_delta)
    dot_exact = vector_glsc2(gradient_exact_before, delta)
    dot_projected = vector_glsc2(gradient_projected_before, delta)
    dot_mma = vector_glsc2(gradient_mma_before, delta)
    write(tag, '(A,I0,A,I0)') 'outer_', outer, '_inner_', inner
    if (inner > 0) then
       write(tag_before, '(A,I0,A,I0)') 'outer_', outer, '_inner_', inner-1
    else
       tag_before = tag
    end if
    filename = trim(this%vector_directory)//'/'//trim(tag)//'_update.bin'
    call rom_audit_write_vector(filename, delta)

    open(newunit=unit, file=trim(this%directory)//'/iteration_diagnostics.csv', &
         status='old', position='append', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not append ROM audit iteration CSV')
    write(unit, '(*(G0,:,","))') outer, inner, merge('initial','step   ',inner==0), &
         trim(tag_before)//'_design.bin', trim(tag)//'_design.bin', &
         trim(tag_before)//'_gradient_selected.bin', trim(tag)//'_gradient_selected.bin', &
         trim(tag_before)//'_gradient_exact.bin', trim(tag_before)//'_gradient_projected.bin', &
         trim(tag)//'_gradient_exact.bin', trim(tag)//'_gradient_projected.bin', &
         trim(tag_before)//'_raw_density_selected.bin', trim(tag)//'_raw_density_selected.bin', &
         trim(tag_before)//'_constraint_gradient.bin', trim(tag)//'_constraint_gradient.bin', &
         rom_objective_before, rom_objective_after, rom_objective_after-rom_objective_before, &
         rom_scalar_before, rom_scalar_after, rom_viscous_before, rom_viscous_after, &
         rom_brinkman_before, rom_brinkman_after, constraints_before(1), constraints_after(1), &
         update_max, update_l2, dot_exact, dot_projected, dot_mma, kkt_max, kkt_l2, &
         trim(this%gradient_source), trim(this%formulation), this%trial_rank, &
         this%test_rank, this%enriched_rank, this%dt, this%n_steps, &
         selected_gradient_device_error
    close(unit)
    call delta%free(); call abs_delta%free()
    call rom_audit_write_transition_barrier()
  end subroutine rom_audit_write_step

  !> Compare exact-Galerkin, projected-adjoint and full-order gradients.
  subroutine rom_audit_write_gradient_comparison(this, point, outer, inner, &
       gradient_exact, gradient_projected, gradient_full)
    class(rom_diagnostic_audit_t), intent(inout) :: this
    character(len=*), intent(in) :: point
    integer, intent(in) :: outer, inner
    type(vector_t), intent(inout) :: gradient_exact, gradient_projected, gradient_full
    type(vector_t) :: diff_exact_full, diff_projected_full, diff_exact_projected
    real(kind=rp) :: norm_exact, norm_projected, norm_full
    real(kind=rp) :: rel_exact_full, rel_projected_full, rel_exact_projected
    real(kind=rp) :: cosine_exact_full, cosine_projected_full, cosine_exact_projected
    real(kind=rp) :: dot_exact_projected, dot_exact_full, dot_projected_full
    integer :: unit, ios

    if (.not. this%enabled) return
    call diff_exact_full%init(gradient_full%size())
    call diff_projected_full%init(gradient_full%size())
    call diff_exact_projected%init(gradient_full%size())
    norm_exact = sqrt(max(0.0_rp, vector_glsc2(gradient_exact, gradient_exact)))
    norm_projected = sqrt(max(0.0_rp, vector_glsc2(gradient_projected, gradient_projected)))
    norm_full = sqrt(max(0.0_rp, vector_glsc2(gradient_full, gradient_full)))
    call vector_sub3(diff_exact_full, gradient_exact, gradient_full)
    call vector_sub3(diff_projected_full, gradient_projected, gradient_full)
    call vector_sub3(diff_exact_projected, gradient_exact, gradient_projected)
    rel_exact_full = sqrt(max(0.0_rp, vector_glsc2(diff_exact_full, diff_exact_full))) / &
         max(1.0_rp, norm_full)
    rel_projected_full = sqrt(max(0.0_rp, vector_glsc2(diff_projected_full, diff_projected_full))) / &
         max(1.0_rp, norm_full)
    rel_exact_projected = sqrt(max(0.0_rp, vector_glsc2(diff_exact_projected, diff_exact_projected))) / &
         max(1.0_rp, norm_projected)
    dot_exact_projected = vector_glsc2(gradient_exact, gradient_projected)
    dot_exact_full = vector_glsc2(gradient_exact, gradient_full)
    dot_projected_full = vector_glsc2(gradient_projected, gradient_full)
    cosine_exact_full = dot_exact_full / max(sqrt(tiny(1.0_rp)), norm_exact*norm_full)
    cosine_projected_full = dot_projected_full / &
         max(sqrt(tiny(1.0_rp)), norm_projected*norm_full)
    cosine_exact_projected = dot_exact_projected / &
         max(sqrt(tiny(1.0_rp)), norm_exact*norm_projected)
    if (pe_rank == 0) then
       open(newunit=unit, file=trim(this%directory)//'/gradient_comparison.csv', &
            status='old', position='append', action='write', iostat=ios)
       if (ios /= 0) call neko_error('could not append ROM audit gradient CSV')
       write(unit, '(*(G0,:,","))') trim(point), outer, inner, &
            norm_exact, norm_projected, norm_full, rel_exact_full, rel_projected_full, &
            rel_exact_projected, cosine_exact_full, cosine_projected_full, &
            cosine_exact_projected, dot_exact_projected, dot_exact_full, dot_projected_full
       close(unit)
    end if
    call diff_exact_full%free(); call diff_projected_full%free(); &
         call diff_exact_projected%free()
    call rom_audit_write_transition_barrier()
  end subroutine rom_audit_write_gradient_comparison

  !> Record matched ROM/full-order objective values and weighted components.
  subroutine rom_audit_write_objective_comparison(this, point, outer, inner, &
       design_file, rom_total, rom_scalar, rom_viscous, rom_brinkman, full_total, &
       full_scalar, full_viscous, full_brinkman)
    class(rom_diagnostic_audit_t), intent(inout) :: this
    character(len=*), intent(in) :: point, design_file
    integer, intent(in) :: outer, inner
    real(kind=rp), intent(in) :: rom_total, rom_scalar, rom_viscous, rom_brinkman
    real(kind=rp), intent(in) :: full_total, full_scalar, full_viscous, full_brinkman
    integer :: unit, ios
    if (.not. this%enabled .or. pe_rank /= 0) return
    open(newunit=unit, file=trim(this%directory)//'/objective_comparison.csv', &
         status='old', position='append', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not append ROM audit objective CSV')
    write(unit, '(*(G0,:,","))') trim(point), outer, inner, &
         trim(design_file), rom_total, rom_scalar, rom_viscous, rom_brinkman, &
         full_total, full_scalar, full_viscous, full_brinkman, rom_total-full_total
    close(unit)
  end subroutine rom_audit_write_objective_comparison

  !> Run ROM-only central finite differences at the outer reference design.
  subroutine rom_audit_finite_difference(this, outer, rom, design, x_reference, &
       x_final, gradient_exact, gradient_projected, first_update)
    class(rom_diagnostic_audit_t), intent(inout) :: this
    integer, intent(in) :: outer
    type(rom_problem_t), intent(inout) :: rom
    class(design_t), intent(inout) :: design
    type(vector_t), intent(inout) :: x_reference, x_final
    type(vector_t), intent(inout) :: gradient_exact, gradient_projected
    type(vector_t), intent(inout) :: first_update
    type(vector_t) :: d, xplus, xminus, xbase
    type(vector_t) :: density_probe, exact_density_probe
    type(vector_t) :: projected_density_probe
    type(vector_t) :: exact_gradient_probe, projected_gradient_probe
    type(vector_t) :: exact_localB_gradient, projected_localB_gradient
    real(kind=rp) :: eps_max, local_eps, norm_d, eps, obj_plus, obj_minus
    real(kind=rp) :: baseline, fd_value, exact_dot, projected_dot
    real(kind=rp) :: exact_localB_dot, projected_localB_dot
    real(kind=rp) :: eps_list(3)
    real(kind=rp), allocatable :: state(:)
    integer :: direction, level, i, unit, ios
    logical :: feasible
    character(len=32) :: direction_name, exact_density_source
    character(len=512) :: direction_file, base_design_file

    if (.not. this%enabled) return
    call x_reference%copy_from(DEVICE_TO_HOST, sync=.true.)
    call x_final%copy_from(DEVICE_TO_HOST, sync=.true.)
    call gradient_exact%copy_from(DEVICE_TO_HOST, sync=.true.)
    call gradient_projected%copy_from(DEVICE_TO_HOST, sync=.true.)
    call first_update%copy_from(DEVICE_TO_HOST, sync=.true.)
    call d%init(x_reference%size())
    call xplus%init(x_reference%size())
    call xminus%init(x_reference%size())
    call xbase%init(x_reference%size())
    call density_probe%init(x_reference%size())
    call exact_density_probe%init(x_reference%size())
    call projected_density_probe%init(x_reference%size())
    call exact_gradient_probe%init(x_reference%size())
    call projected_gradient_probe%init(x_reference%size())
    call exact_localB_gradient%init(x_reference%size())
    call projected_localB_gradient%init(x_reference%size())

    exact_density_source = 'exact_galerkin'
    if (rom%experimental_algebraic) exact_density_source = 'exact_algebraic'

    do direction = 1, 2
       ! The reference design is a 0/1 corner.  Use the midpoint of the first
       ! bounded MMA move as a common interior probe and recompute gradients
       ! there, so the central FD and both analytic directions share a design.
       xbase%x = x_reference%x + 0.5_rp * first_update%x
       xbase%x = max(this%mma_xmin, min(this%mma_xmax, xbase%x))
       if (NEKO_BCKND_DEVICE .eq. 1) call xbase%copy_from(HOST_TO_DEVICE, sync=.true.)
       call design%update_design(xbase)
       call rom%compute(design, baseline, final_state=state)
       if (allocated(state)) deallocate(state)
       call rom%compute_sensitivity(density_probe)
       call rom%compute_density(exact_density_probe, trim(exact_density_source))
       call rom%compute_density(projected_density_probe, 'projected_adjoint')
       call rom%project_density(design, exact_density_probe, exact_gradient_probe)
       call rom%project_density(design, projected_density_probe, projected_gradient_probe)
       call design%map_backward(exact_density_probe)
       select type (des => design)
       type is (brinkman_design_t)
          call des%get_sensitivity(exact_localB_gradient)
          call des%convert_to_directional_derivative(exact_localB_gradient)
       class default
          call neko_error('ROM local-B diagnostic requires a Brinkman design')
       end select
       call design%map_backward(projected_density_probe)
       select type (des => design)
       type is (brinkman_design_t)
          call des%get_sensitivity(projected_localB_gradient)
          call des%convert_to_directional_derivative(projected_localB_gradient)
       class default
          call neko_error('ROM local-B diagnostic requires a Brinkman design')
       end select
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call xbase%copy_from(DEVICE_TO_HOST, sync=.true.)
          call exact_gradient_probe%copy_from(DEVICE_TO_HOST, sync=.true.)
          call projected_gradient_probe%copy_from(DEVICE_TO_HOST, sync=.true.)
          call exact_localB_gradient%copy_from(DEVICE_TO_HOST, sync=.true.)
          call projected_localB_gradient%copy_from(DEVICE_TO_HOST, sync=.true.)
       end if
       if (direction == 1) then
          direction_name = 'mma_update'
          d%x = first_update%x
       else
          direction_name = 'exact_gradient'
          d%x = exact_gradient_probe%x
       end if
       do i = 1, d%size()
          if (xbase%x(i) <= this%mma_xmin + 1.0e-10_rp .or. &
               xbase%x(i) >= this%mma_xmax - 1.0e-10_rp) d%x(i) = 0.0_rp
       end do
       ! Normalise with the global inner product used by ROM/MMA, not a local
       ! rank sum, so the diagnostic direction is MPI-rank independent.
       if (NEKO_BCKND_DEVICE .eq. 1) call d%copy_from(HOST_TO_DEVICE, sync=.true.)
       norm_d = sqrt(max(0.0_rp, vector_glsc2(d, d)))
       if (norm_d > 0.0_rp) then
          call d%copy_from(DEVICE_TO_HOST, sync=.true.)
          d%x = d%x / norm_d
          local_eps = huge(1.0_rp)
          do i = 1, d%size()
             if (abs(d%x(i)) > tiny(1.0_rp)) local_eps = min(local_eps, &
                  min((xbase%x(i)-this%mma_xmin)/abs(d%x(i)), &
                  (this%mma_xmax-xbase%x(i))/abs(d%x(i))))
          end do
          if (NEKO_BCKND_DEVICE .eq. 1) call d%copy_from(HOST_TO_DEVICE, sync=.true.)
          call MPI_Allreduce(local_eps, eps_max, 1, MPI_REAL_PRECISION, MPI_MIN, &
               NEKO_COMM)
          if (eps_max == huge(1.0_rp)) eps_max = 0.0_rp
       else
          eps_max = 0.0_rp
       end if
       feasible = eps_max > 0.0_rp .and. norm_d > 0.0_rp
       if (feasible) then
          eps_list = [min(1.0e-3_rp, 0.25_rp*eps_max), &
               min(1.0e-4_rp, 0.025_rp*eps_max), &
               min(1.0e-5_rp, 0.0025_rp*eps_max)]
          direction_file = trim(this%vector_directory)//'/outer_'// &
               trim(adjustl(itoa(outer)))//'_'//trim(direction_name)//'.bin'
          base_design_file = trim(this%vector_directory)//'/outer_'// &
               trim(adjustl(itoa(outer)))//'_'//trim(direction_name)//'_base.bin'
          call rom_audit_write_vector(direction_file, d)
          call rom_audit_write_vector(base_design_file, xbase)
          do level = 1, 3
             eps = eps_list(level)
             if (eps <= 0.0_rp) cycle
             xplus%x = xbase%x + eps*d%x
             xminus%x = xbase%x - eps*d%x
             if (NEKO_BCKND_DEVICE .eq. 1) then
                call xplus%copy_from(HOST_TO_DEVICE, sync=.true.)
                call xminus%copy_from(HOST_TO_DEVICE, sync=.true.)
             end if
             call design%update_design(xplus)
             call rom%compute(design, obj_plus, final_state=state)
             if (allocated(state)) deallocate(state)
             call design%update_design(xminus)
             call rom%compute(design, obj_minus, final_state=state)
             if (allocated(state)) deallocate(state)
             fd_value = (obj_plus-obj_minus)/(2.0_rp*eps)
             exact_dot = vector_glsc2(exact_gradient_probe, d)
             projected_dot = vector_glsc2(projected_gradient_probe, d)
             exact_localB_dot = vector_glsc2(exact_localB_gradient, d)
             projected_localB_dot = vector_glsc2(projected_localB_gradient, d)
             if (pe_rank == 0) then
                open(newunit=unit, file=trim(this%directory)//'/directional_fd.csv', &
                     status='old', position='append', action='write', iostat=ios)
                if (ios /= 0) call neko_error('could not append ROM audit FD CSV')
                write(unit, '(*(G0,:,","))') outer, trim(direction_name), eps, feasible, &
                     obj_plus, obj_minus, fd_value, exact_dot, projected_dot, &
                     fd_value-exact_dot, fd_value-projected_dot, &
                     exact_localB_dot, projected_localB_dot, &
                     fd_value-exact_localB_dot, fd_value-projected_localB_dot, &
                     trim(direction_file), trim(base_design_file)
                close(unit)
             end if
             call design%update_design(xbase)
             call rom%compute(design, baseline, final_state=state)
             if (allocated(state)) deallocate(state)
             call rom_audit_write_transition_barrier()
          end do
       else if (pe_rank == 0) then
          open(newunit=unit, file=trim(this%directory)//'/directional_fd.csv', &
               status='old', position='append', action='write', iostat=ios)
          if (ios /= 0) call neko_error('could not append ROM audit FD CSV')
          write(unit, '(*(G0,:,","))') outer, trim(direction_name), 0.0_rp, .false., &
               0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
               0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 'none', 'none'
          close(unit)
       end if
       call design%update_design(xbase)
       call rom_audit_write_transition_barrier()
    end do
    call design%update_design(x_final)
    call rom_audit_write_transition_barrier()
    call d%free(); call xplus%free(); call xminus%free(); call xbase%free()
    call density_probe%free(); call exact_density_probe%free(); &
         call projected_density_probe%free(); call exact_gradient_probe%free(); &
         call projected_gradient_probe%free(); call exact_localB_gradient%free(); &
         call projected_localB_gradient%free()
  end subroutine rom_audit_finite_difference

  !> Gather a vector from all ROM ranks and write a raw float64 file on rank 0.
  subroutine rom_audit_write_vector(filename, values)
    character(len=*), intent(in) :: filename
    type(vector_t), intent(inout) :: values
    call rom_audit_write_array(filename, values%x, values)
  end subroutine rom_audit_write_vector

  !> Array writer with an optional vector whose device mirror must be synced.
  subroutine rom_audit_write_array(filename, local_values, vector)
    character(len=*), intent(in) :: filename
    real(kind=rp), intent(in) :: local_values(:)
    type(vector_t), intent(inout), optional :: vector
    integer, allocatable :: counts(:), displs(:)
    real(kind=rp), allocatable :: gathered(:)
    integer :: local_n, total_n, i, ios, unit
    if (present(vector)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call vector%copy_from(DEVICE_TO_HOST, sync=.true.)
    end if
    local_n = size(local_values)
    allocate(counts(pe_size), displs(pe_size))
    call MPI_Gather(local_n, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, NEKO_COMM)
    total_n = 0
    if (pe_rank == 0) then
       displs(1) = 0
       do i = 1, pe_size
          if (i > 1) displs(i) = displs(i-1) + counts(i-1)
          total_n = total_n + counts(i)
       end do
       allocate(gathered(max(1,total_n)))
    else
       displs = 0
       allocate(gathered(1))
    end if
    call MPI_Gatherv(local_values, local_n, MPI_REAL_PRECISION, gathered, counts, &
         displs, MPI_REAL_PRECISION, 0, NEKO_COMM)
    if (pe_rank == 0) then
       open(newunit=unit, file=trim(filename), status='replace', access='stream', &
            form='unformatted', action='write', iostat=ios)
       if (ios /= 0) call neko_error('could not write ROM audit vector file')
       write(unit) gathered(1:total_n)
       close(unit)
    end if
    deallocate(counts, displs, gathered)
    call MPI_Barrier(NEKO_COMM)
  end subroutine rom_audit_write_array

  !> Write each matrix row as a gathered vector file.
  subroutine rom_audit_write_matrix(filename, values)
    character(len=*), intent(in) :: filename
    type(matrix_t), intent(inout) :: values
    integer :: i, nrows
    character(len=512) :: row_filename
    if (NEKO_BCKND_DEVICE .eq. 1) call values%copy_from(DEVICE_TO_HOST, sync=.true.)
    nrows = values%get_nrows()
    do i = 1, nrows
       write(row_filename, '(A,"_row",I0,".bin")') trim(filename), i
       call rom_audit_write_array(row_filename, values%x(i,:))
    end do
  end subroutine rom_audit_write_matrix

  !> Write the local partition metadata for reproducible vector reconstruction.
  subroutine rom_audit_write_partition(this)
    class(rom_diagnostic_audit_t), intent(in) :: this
    integer :: local_n, ios, unit
    integer, allocatable :: counts(:)
    local_n = this%n_local
    allocate(counts(pe_size))
    call MPI_Gather(local_n, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, NEKO_COMM)
    if (pe_rank == 0) then
       open(newunit=unit, file=trim(this%directory)//'/vector_partition.csv', &
            status='replace', action='write', iostat=ios)
       if (ios /= 0) call neko_error('could not append ROM audit partition CSV')
       write(unit, '(A)') 'rank,local_size,global_offset'
       local_n = 0
       do ios = 1, pe_size
          write(unit, '(I0,",",I0,",",I0)') ios-1, counts(ios), local_n
          local_n = local_n + counts(ios)
       end do
       close(unit)
    end if
    deallocate(counts)
    call MPI_Barrier(NEKO_COMM)
  end subroutine rom_audit_write_partition

  !> Synchronise ranks after a rank-zero-only diagnostic write.
  subroutine rom_audit_write_transition_barrier()
    call MPI_Barrier(NEKO_COMM)
  end subroutine rom_audit_write_transition_barrier

  !> Read an optional environment value for provenance metadata.
  subroutine rom_audit_environment(name, value)
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(out) :: value
    character(len=512) :: buffer
    integer :: length, status
    call get_environment_variable(name, buffer, length, status)
    if (status /= 0 .or. length == 0) then
       value = 'unknown'
    else
       value = buffer(1:length)
    end if
  end subroutine rom_audit_environment

  !> Return an integer as a compact character string.
  function itoa(value) result(string)
    integer, intent(in) :: value
    character(len=32) :: string
    write(string, '(I0)') value
  end function itoa

end module rom_diagnostic_audit
