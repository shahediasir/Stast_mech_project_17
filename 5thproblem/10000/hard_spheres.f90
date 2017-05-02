program hardsphere

   !! Hard sphere Metropolis Monte Carlo program to calculate the radial
   !! distribution function.
   !!
   !! g(r) is written to standard output at the end of the simulation;
   !! all other messages are written to standard error
   !!
   !! Written April 23-25, 2017 by Karl D. Hammond, University of Missouri
   !! Not intended for distribution.

   use, intrinsic :: ISO_Fortran_env

   implicit none

   integer, parameter :: nstep = 10000, nevery = 1000
   double precision, parameter :: sigma = 1.0D0, sigmasq = sigma**2, &
      rho = 0.75D0 ! reduced density
   double precision, parameter :: PI = acos(-1.0D0)
   double precision, parameter :: max_dx = 0.1*sigma, gmax = 2.50D0, &
      gmaxsq = gmax**2, dg = 0.05D0
   double precision, parameter :: k_B = 1.38064852D-23

   integer :: N
   double precision, dimension(:,:), allocatable :: x
   double precision, dimension(:), allocatable :: g
   integer :: i, j, k, m, n_x, n_y, n_z, naccept
   double precision :: V, V0, L_x, L_y, L_z, E, min_a, a_x, a_y, a_z, r
   double precision :: dx, dy, dz, dr(3), drsq
   logical :: accept

   ! The default creates very non-uniformly-distributed particles; this is
   ! an inefficiency that I would hope to correct had I had a bit more time....
   L_x = 10.0
   L_y = L_x ! these should be /multiples/ of a common factor of Lx
   L_z = L_x
   V = L_x * L_y * L_z
   N = rho * V / sigma**3 ! chosen such that rho* is equal to the value given
   allocate(x(3,N))
   allocate(g(int(gmax/dg+0.01)))
   g = 0.0D0
   x = 0.0D0

   ! First, find V0 from geometric concerns (close-packing)
   V0 = N * sigma**3 / sqrt(2.0)

   ! If V < V0, we have overlap => abort!
   if ( V < V0 ) then
      stop 'The specified density exceeds close-packed density; aborting'
   end if

   ! Place particles; assume they're FCC to start (produces non-uniform, but
   ! definitely non-overlapping, starting configurations)
   min_a = sigma * sqrt(2.0)
   n_x = int(L_x / min_a)
   n_y = int(L_y / min_a)
   n_z = int(L_z / min_a)
   a_x = L_x/n_x
   a_y = L_y/n_y
   a_z = L_z/n_z
   m = 1
   creation: do i = 0, n_x - 1
      do j = 0, n_y - 1
         do k = 0, n_z - 1
            ! Each unit cell contains four atoms
            x(:,m) = [i*a_x, j*a_y, k*a_z]
            m = m + 1
            if ( m > N ) exit creation
            x(:,m) = [i*a_x + 0.5D0*a_x, j*a_y + 0.5D0*a_y, k*a_z]
            m = m + 1
            if ( m > N ) exit creation
            x(:,m) = [i*a_x + 0.5D0*a_x, j*a_y, k*a_z + 0.5D0*a_z]
            m = m + 1
            if ( m > N ) exit creation
            x(:,m) = [i*a_x, j*a_y + 0.5D0*a_y, k*a_z + 0.5D0*a_z]
            m = m + 1
            if ( m > N ) exit creation
         end do
      end do
   end do creation

   ! Scale coordinates so that everything is a bit more uniform in all
   ! three spatial dimensions
   dx = maxval(x(1,:)) - minval(x(1,:)) + sigma/sqrt(2.0)
   dy = maxval(x(2,:)) - minval(x(2,:)) + sigma/sqrt(2.0)
   dz = maxval(x(3,:)) - minval(x(3,:)) + sigma/sqrt(2.0)
   x(1,:) = x(1,:) * L_x / dx
   x(2,:) = x(2,:) * L_y / dy
   x(3,:) = x(3,:) * L_z / dz

   ! Initialize random number generator
   call seed_RNG ! uses system clock by default to generate seed

   write (ERROR_UNIT,'(I0,1X,A)') N, 'particles in the simulation'

   !! Metropolis-Hastings algorithm
   do j = 1, nstep
      ! Pick a particle at random
      call random_number(r)
      i = int(r*N) + 1
      ! Move it by a random displacement
      call random_number(r)
      dx = 2.0D0 * max_dx * (r-0.5D0)
      call random_number(r)
      dy = 2.0D0 * max_dx * (r-0.5D0)
      call random_number(r)
      dz = 2.0D0 * max_dx * (r-0.5D0)
      dr = [dx,dy,dz]
      ! Accept or reject move
      accept = .TRUE.
      do k = 1, N
         ! Particle left the box
         if ( any(x(:,i) + dr >= [L_x,L_y,L_z]) .or. &
              any(x(:,i) + dr < 0.0D0) ) then
            accept = .FALSE.
            exit
         end if
         if ( k == i ) cycle ! can't overlap with yourself!
         ! Check for overlapping spheres
         drsq = sum( ( x(:,k) - (x(:,i) + dr) )**2 )
         if ( drsq < sigmasq ) then
            accept = .FALSE.
            exit
         end if
      end do
      if ( accept ) then
         x(:,i) = x(:,i) + dr
         naccept = naccept + 1
      end if
      ! Tally this g(r) value
      do k = 1, N
         do i = k + 1, N ! only look for half of them!
            drsq = sum( ( x(:,k) - (x(:,i) + dr) )**2 )
            if ( drsq > gmaxsq ) cycle
            m = int(sqrt(drsq) / dg) + 1 ! +1 because arrays start at 1 (duh)
            if ( m > size(g) ) write (ERROR_UNIT,'(A,1X,I0)') 'ERROR: m = ', m
            ! we use 2.0 here because we're counting HALF of the particles
            g(m) = g(m) + 2.0D0 / (4.0D0*PI*drsq * N/V * (N-1) * nstep) ! Chandler, page 222
         end do
      end do
      if ( modulo(j,nevery) == 0 ) then
         write (ERROR_UNIT, '(I0,1X,A,1X,I0,1X,A,F0.2,A)') naccept, 'of', j, &
            'moves accepted (', real(naccept)/j * 100.0, '%)'
      end if
   end do

   ! Print g(r)
   do i = 1, size(g)
      write(OUTPUT_UNIT,*) dg*i, g(i)
   end do

   deallocate (x)
   deallocate (g)

contains

   subroutine seed_RNG (seed)
      integer, optional :: seed
      integer :: i, n, clock
      integer, dimension(:), allocatable :: rseed
      call random_seed (size = n)
      allocate (rseed(n))
      if ( present(seed) ) then
         clock = seed
      else
         call system_clock (count=clock)
      end if
      rseed = clock + 37 * [ (i-1, i=1,n) ]
      call random_seed (put = rseed)
   end subroutine seed_RNG

end program hardsphere
