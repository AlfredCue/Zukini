!
! ZDPLASKIN version 2.0a
! (c) 2008, Sergey Pancheshnyi (pancheshnyi@gmail.com)
!
! BOLSIG+
! (c) 2005, Gerjan Hagelaar (gerjan.hagelaar@laplace.univ-tlse.fr)
!
! http://www.zdplaskin.laplace.univ-tlse.fr/
! This software is provided "as is" without warranty and non-commercial use is freely
! granted provided proper reference is made in publications resulting from its use.
! Use of ZDPlasKin in commerical software requires a license.
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Sun Aug 17 18:39:13 2025
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS: E N O 
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE ZDPlasKin
!
!-----------------------------------------------------------------------------------------------------------------------------------
module ZDPlasKin
  use dvode_f90_m, only : vode_opts
  implicit none
  public
!
! config
!
  integer, parameter :: species_max = 53, species_electrons = 53, species_length = 9, reactions_max = 246, reactions_length = 44
  double precision                          :: density(species_max)
  integer                                   :: species_charge(species_max)
  character(species_length)                 :: species_name(species_max)
  character(reactions_length)               :: reaction_sign(reactions_max)
  logical                                   :: lreaction_block(reactions_max)
!
! internal config
!
  double precision, parameter, private      :: vode_atol = 1.00D-10, vode_rtol = 1.00D-05, cfg_rtol = 1.00D-01
  integer, parameter, private               :: vode_neq = species_max + 1
  type (vode_opts), private                 :: vode_options
  integer, private                          :: vode_itask, vode_istate, ifile_unit = 5
  double precision, private                 :: stat_dens(species_max), stat_src(species_max), stat_rrt(reactions_max), stat_time, &
                                               dens_loc(vode_neq,0:3), rrt_loc(reactions_max), tsav = -huge(tsav), &
                                               mach_accur, mach_tiny
  double precision                          :: rrt(reactions_max), mrtm(species_max, reactions_max), ZDPlasKin_cfg(14)
  logical, private                          :: lZDPlasKin_init = .false., lprint, lstat_accum, ldensity_constant, &
                                               density_constant(species_max), lgas_heating
!
! qtplaskin config
!
  logical, private                          :: lqtplaskin, lqtplaskin_first = .true.
  double precision, parameter, private      :: qtplaskin_atol = 1.00D+00, qtplaskin_rtol = 1.00D-02
  character(32), allocatable                :: qtplaskin_user_names(:)
  double precision, allocatable             :: qtplaskin_user_data(:)
!
! physical constants
!
  double precision, parameter, private      :: eV_to_K = 1.16045052d4, q_elem = 1.60217662d-19, k_B = 1.38064852d-23
!
! bolsig+ config
!
  double precision, parameter, private      :: bolsig_rtol = 1.00D-03, bolsig_rtol_half = 3.16D-02, &
                                               bolsig_field_min = 1.00D-01, bolsig_field_max = 1.00D+03, &
                                               bolsig_eecol_frac_def = 1.00D-05
  double precision, private                 :: bolsig_eecol_frac
  integer, parameter, private               :: bolsig_species_max = 21, bolsig_species_length = 6, bolsig_rates_max = 30, &
                                               bolsig_addsect_max = 12 
  character(*), parameter, private          :: bolsigfile = "bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0, &
                                               bolsig_addsect(2,bolsig_addsect_max) 
  logical, private                          :: lbolsig_ignore_gas_temp, lbolsig_Maxwell_EEDF
  double precision, allocatable             :: bolsig_rates(:)
  character(bolsig_species_length), private :: bolsig_species(bolsig_species_max)
  interface
    subroutine ZDPlasKin_bolsig_Init(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_Init
    subroutine ZDPlasKin_bolsig_ReadCollisions(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_ReadCollisions
    subroutine ZDPlasKin_bolsig_GetCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetCollisions
    subroutine ZDPlasKin_bolsig_GetSpeciesName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetSpeciesName
    subroutine ZDPlasKin_bolsig_GetReactionName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetReactionName
    subroutine ZDPlasKin_bolsig_SolveBoltzmann(i,a,j,b)
      integer, intent(in) :: i, j
      double precision, intent(in)  :: a(i)
      double precision, intent(out) :: b(j)
    end subroutine ZDPlasKin_bolsig_SolveBoltzmann
    subroutine ZDPlasKin_bolsig_GetEEDF(i,a,b)
      integer, intent(in) :: i
      double precision, intent(out) :: a,b
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0,&
    0, 0, 1, 1, 1,-1,-1,-1,-1, 1,-1/
  data species_name(1:species_max) &
  /"N2       ","N2(V1)   ","N2(V2)   ","N2(V3)   ","N2(V4)   ","N2(V5)   ","N2(V6)   ","N2(V7)   ","N2(V8)   ","N2(A3)   ",&
   "N2(B3)   ","N2(A`1)  ","N2(C3)   ","N        ","N(2D)    ","N(2P)    ","N^+      ","N2^+     ","N3^+     ","N4^+     ",&
   "O2       ","O2(V1)   ","O2(V2)   ","O2(V3)   ","O2(V4)   ","O2(A1)   ","O2(B1)   ","O2(4.5EV)","O        ","O(1D)    ",&
   "O(1S)    ","O3       ","O^+      ","O2^+     ","O4^+     ","O^-      ","O2^-     ","O3^-     ","O4^-     ","NO       ",&
   "N2O      ","NO2      ","NO3      ","N2O5     ","NO^+     ","N2O^+    ","NO2^+    ","NO^-     ","N2O^-    ","NO2^-    ",&
   "NO3^-    ","O2^+N2   ","E        "/
  data reaction_sign(1:36) &
  /"bolsig:N2->N2(V1RES)                        ","bolsig:N2->N2(V1)                           ",&
   "bolsig:N2->N2(V2)                           ","bolsig:N2->N2(V3)                           ",&
   "bolsig:N2->N2(V4)                           ","bolsig:N2->N2(V5)                           ",&
   "bolsig:N2->N2(V6)                           ","bolsig:N2->N2(V7)                           ",&
   "bolsig:N2->N2(V8)                           ","bolsig:N2(V1)->N2                           ",&
   "bolsig:N2(V2)->N2                           ","bolsig:N2(V3)->N2                           ",&
   "bolsig:N2(V4)->N2                           ","bolsig:N2(V5)->N2                           ",&
   "bolsig:N2(V6)->N2                           ","bolsig:N2(V7)->N2                           ",&
   "bolsig:N2(V8)->N2                           ","N2(V1)+N2=>N2+N2                            ",&
   "N2(V2)+N2=>N2(V1)+N2                        ","N2(V3)+N2=>N2(V2)+N2                        ",&
   "N2(V4)+N2=>N2(V3)+N2                        ","N2(V5)+N2=>N2(V4)+N2                        ",&
   "N2(V6)+N2=>N2(V5)+N2                        ","N2(V7)+N2=>N2(V6)+N2                        ",&
   "N2(V8)+N2=>N2(V7)+N2                        ","N2(V1)+O=>N2+O                              ",&
   "N2(V2)+O=>N2(V1)+O                          ","N2(V3)+O=>N2(V2)+O                          ",&
   "N2(V4)+O=>N2(V3)+O                          ","N2(V5)+O=>N2(V4)+O                          ",&
   "N2(V6)+O=>N2(V5)+O                          ","N2(V7)+O=>N2(V6)+O                          ",&
   "N2(V8)+O=>N2(V7)+O                          ","bolsig:N2->N2(A3)                           ",&
   "bolsig:N2->N2(A3,V5-9)                      ","bolsig:N2->N2(A3,V10-)                      "/
  data reaction_sign(37:72) &
  /"bolsig:N2->N2(B3)                           ","bolsig:N2->N2(W3)                           ",&
   "bolsig:N2->N2(B'3)                          ","bolsig:N2->N2(A'1)                          ",&
   "bolsig:N2->N2(A1)                           ","bolsig:N2->N2(W1)                           ",&
   "bolsig:N2->N2(C3)                           ","bolsig:N2->N2(E3)                           ",&
   "bolsig:N2->N2(SUM)                          ","bolsig:NO->NO^+                             ",&
   "E+NO2^+=>NO+O                               ","E+NO2=>O^-+NO                               ",&
   "E+NO+ANY_NEUTRAL=>NO^-+ANY_NEUTRAL          ","O^-+N=>NO+E                                 ",&
   "O^-+NO=>NO2+E                               ","O2^-+N=>NO2+E                               ",&
   "O3^-+N=>NO+O2+E                             ","N2O^-+N=>NO+N2+E                            ",&
   "NO2^-+N=>NO+NO+E                            ","NO3^-+N=>NO+NO2+E                           ",&
   "NO^-+O=>NO2+E                               ","N2O^-+O=>NO+NO+E                            ",&
   "NO2^-+O=>NO+O2+E                            ","NO3^-+O=>NO+O3+E                            ",&
   "NO^-+N2(A3)=>NO+N2+E                        ","NO2^-+N2(A3)=>NO2+N2+E                      ",&
   "NO^-+N2(B3)=>NO+N2+E                        ","NO2^-+N2(B3)=>NO2+N2+E                      ",&
   "N2(C3)=>N2(B3)                              ","N2(A3)+O=>NO+N(2D)                          ",&
   "N2(A3)+O2=>N2+O+O(1D)                       ","N2(A3)+O2=>N2+O2(B1)                        ",&
   "N2(A3)+N2=>N2+N2                            ","N2(A3)+NO=>N2+NO                            ",&
   "N2(A3)+N2O=>N2+N+NO                         ","N2(A3)+NO2=>N2+O+NO                         "/
  data reaction_sign(73:108) &
  /"N2(A3)+N2(A3)=>N2+N2(B3)                    ","N2(A3)+N2(A3)=>N2+N2(C3)                    ",&
   "N2(B3)+N2=>N2(A3)+N2                        ","N2(B3)+N2=>N2+N2                            ",&
   "N2(B3)+O2=>N2+O+O                           ","N2(B3)+NO=>N2(A3)+NO                        ",&
   "N2(C3)+N2=>N2(A`1)+N2                       ","N2(A`1)+N2=>N2(B3)+N2                       ",&
   "N2(A`1)+O2=>N2+O+O                          ","N2(A`1)+NO=>N2+N+O                          ",&
   "N(2D)+O2=>NO+O                              ","N(2D)+NO=>N2+O                              ",&
   "N(2D)+N2O=>NO+N2                            ","N(2P)+O2=>NO+O                              ",&
   "N(2P)+NO=>N2(A3)+O                          ","O2(A1)+N=>NO+O                              ",&
   "O(1D)+N2=>O+N2                              ","O(1D)+NO=>O2+N                              ",&
   "O(1D)+N2O=>NO+NO                            ","N+NO=>O+N2                                  ",&
   "N+O2=>O+NO                                  ","N+NO2=>O+O+N2                               ",&
   "N+NO2=>O+N2O                                ","N+NO2=>N2+O2                                ",&
   "N+NO2=>NO+NO                                ","O+NO=>N+O2                                  ",&
   "O+NO=>NO2                                   ","O+N2O=>NO+NO                                ",&
   "O+NO2=>NO+O2                                ","O+NO3=>O2+NO2                               ",&
   "NO+O2=>O+NO2                                ","NO+O3=>O2+NO2                               ",&
   "NO+N2O=>N2+NO2                              ","NO+NO3=>NO2+NO2                             ",&
   "O2+NO2=>NO+O3                               ","NO2+NO2=>NO+NO+O2                           "/
  data reaction_sign(109:144) &
  /"NO2+NO2=>NO+NO3                             ","NO2+O3=>O2+NO3                              ",&
   "NO2+NO3=>NO+NO2+O2                          ","NO3+O2=>NO2+O3                              ",&
   "NO3+NO3=>O2+NO2+NO2                         ","NO2+N2=>NO+O+N2                             ",&
   "NO3+N2=>NO2+O+N2                            ","NO3+O2=>NO2+O+O2                            ",&
   "NO3+NO=>NO2+O+NO                            ","NO3+N=>NO2+O+N                              ",&
   "NO3+O=>NO2+O+O                              ","NO3+N2=>NO+O2+N2                            ",&
   "NO3+O2=>NO+O2+O2                            ","NO3+NO=>NO+O2+NO                            ",&
   "NO3+N=>NO+O2+N                              ","NO3+O=>NO+O2+O                              ",&
   "N2O5+ANY_NEUTRAL=>NO2+NO3+ANY_NEUTRAL       ","N+O+N2=>NO+N2                               ",&
   "N+O+O2=>NO+O2                               ","N+O+N=>NO+N                                 ",&
   "N+O+O=>NO+O                                 ","N+O+NO=>NO+NO                               ",&
   "O+O2+N2=>O3+N2                              ","O+NO+N2=>NO2+N2                             ",&
   "O+NO+O2=>NO2+O2                             ","O+NO+NO=>NO2+NO                             ",&
   "O+NO2+N2=>NO3+N2                            ","O+NO2+O2=>NO3+O2                            ",&
   "O+NO2+N=>NO3+N                              ","O+NO2+O=>NO3+O                              ",&
   "O+NO2+NO=>NO3+NO                            ","NO2+NO3+ANY_NEUTRAL=>N2O5+ANY_NEUTRAL       ",&
   "N^++O2=>O^++NO                              ","N^++NO=>NO^++N                              ",&
   "N^++NO=>N2^++O                              ","N^++NO=>O^++N2                              "/
  data reaction_sign(145:180) &
  /"O^++NO=>NO^++O                              ","O^++NO=>O2^++N                              ",&
   "O^++N2O=>NO^++NO                            ","O^++NO2=>NO2^++O                            ",&
   "N2^++NO=>NO^++N2                            ","O2^++N2=>NO^++NO                            ",&
   "O2^++NO=>NO^++O2                            ","O2^++NO2=>NO^++O3                           ",&
   "O2^++NO2=>NO2^++O2                          ","N3^++NO=>NO^++N+N2                          ",&
   "N3^++NO=>N2O^++N2                           ","NO2^++NO=>NO^++NO2                          ",&
   "N2O^++NO=>NO^++N2O                          ","N4^++NO=>NO^++N2+N2                         ",&
   "O4^++NO=>NO^++O2+O2                         ","O^-+NO2=>NO2^-+O                            ",&
   "O^-+N2O=>NO^-+NO                            ","O2^-+NO2=>NO2^-+O2                          ",&
   "O3^-+NO=>NO3^-+O                            ","O3^-+NO=>NO2^-+O2                           ",&
   "O3^-+NO2=>NO2^-+O3                          ","O3^-+NO2=>NO3^-+O2                          ",&
   "NO^-+O2=>O2^-+NO                            ","NO^-+NO2=>NO2^-+NO                          ",&
   "NO2^-+NO2=>NO3^-+NO                         ","NO2^-+NO3=>NO3^-+NO2                        ",&
   "NO2^-+N2O5=>NO3^-+NO2+NO2                   ","NO3^-+NO=>NO2^-+NO2                         ",&
   "O4^-+NO=>NO3^-+O2                           ","O^-+NO+ANY_NEUTRAL=>NO2^-+ANY_NEUTRAL       ",&
   "O^-+NO^+=>O+NO                              ","O^-+NO2^+=>O+NO2                            ",&
   "O2^-+NO^+=>O2+NO                            ","O2^-+NO2^+=>O2+NO2                          ",&
   "O3^-+NO^+=>O3+NO                            ","O3^-+NO2^+=>O3+NO2                          "/
  data reaction_sign(181:216) &
  /"NO^-+N^+=>NO+N                              ","NO^-+N2^+=>NO+N2                            ",&
   "NO^-+O^+=>NO+O                              ","NO^-+O2^+=>NO+O2                            ",&
   "NO^-+NO^+=>NO+NO                            ","NO^-+N2O^+=>NO+N2O                          ",&
   "NO^-+NO2^+=>NO+NO2                          ","N2O^-+NO^+=>N2O+NO                          ",&
   "N2O^-+NO2^+=>N2O+NO2                        ","NO2^-+N^+=>NO2+N                            ",&
   "NO2^-+N2^+=>NO2+N2                          ","NO2^-+O^+=>NO2+O                            ",&
   "NO2^-+O2^+=>NO2+O2                          ","NO2^-+NO^+=>NO2+NO                          ",&
   "NO2^-+N2O^+=>NO2+N2O                        ","NO2^-+NO2^+=>NO2+NO2                        ",&
   "NO3^-+N2^+=>NO3+N2                          ","NO3^-+NO^+=>NO3+NO                          ",&
   "NO3^-+NO2^+=>NO3+NO2                        ","O2^-+O2^+N2=>O2+O2+N2                       ",&
   "NO^-+N2^+=>NO+N+N                           ","NO^-+N3^+=>NO+N+N2                          ",&
   "NO^-+N4^+=>NO+N2+N2                         ","NO^-+O2^+=>NO+O+O                           ",&
   "NO^-+O4^+=>NO+O2+O2                         ","NO^-+NO^+=>NO+N+O                           ",&
   "NO^-+N2O^+=>NO+N2+O                         ","NO^-+NO2^+=>NO+N+O2                         ",&
   "NO^-+O2^+N2=>NO+O2+N2                       ","N2O^-+N2^+=>N2O+N+N                         ",&
   "NO2^-+N2^+=>NO2+N+N                         ","NO2^-+N3^+=>NO2+N+N2                        ",&
   "NO2^-+N4^+=>NO2+N2+N2                       ","NO2^-+O2^+=>NO2+O+O                         ",&
   "NO2^-+O4^+=>NO2+O2+O2                       ","NO2^-+NO^+=>NO2+N+O                         "/
  data reaction_sign(217:246) &
  /"NO2^-+N2O^+=>NO2+N2+O                       ","NO2^-+NO2^+=>NO2+N+O2                       ",&
   "NO2^-+O2^+N2=>NO2+O2+N2                     ","O4^-+NO^+=>O2+O2+NO                         ",&
   "O4^-+NO2^+=>O2+O2+NO2                       ","O^-+NO^++ANY_NEUTRAL=>O+NO+ANY_NEUTRAL      ",&
   "O2^-+NO^++ANY_NEUTRAL=>O2+NO+ANY_NEUTRAL    ","O^-+N^++ANY_NEUTRAL=>NO+ANY_NEUTRAL         ",&
   "O^-+NO^++ANY_NEUTRAL=>NO2+ANY_NEUTRAL       ","O2^-+N^++ANY_NEUTRAL=>NO2+ANY_NEUTRAL       ",&
   "O3^-+NO^++ANY_NEUTRAL=>O3+NO+ANY_NEUTRAL    ","O3^-+NO2^++ANY_NEUTRAL=>O3+NO2+ANY_NEUTRAL  ",&
   "NO^-+N^++ANY_NEUTRAL=>NO+N+ANY_NEUTRAL      ","NO^-+N2^++ANY_NEUTRAL=>NO+N2+ANY_NEUTRAL    ",&
   "NO^-+O^++ANY_NEUTRAL=>NO+O+ANY_NEUTRAL      ","NO^-+O2^++ANY_NEUTRAL=>NO+O2+ANY_NEUTRAL    ",&
   "NO^-+NO^++ANY_NEUTRAL=>NO+NO+ANY_NEUTRAL    ","NO^-+N2O^++ANY_NEUTRAL=>NO+N2O+ANY_NEUTRAL  ",&
   "NO^-+NO2^++ANY_NEUTRAL=>NO+NO2+ANY_NEUTRAL  ","N2O^-+NO^++ANY_NEUTRAL=>N2O+NO+ANY_NEUTRAL  ",&
   "N2O^-+NO2^++ANY_NEUTRAL=>N2O+NO2+ANY_NEUTRAL","NO2^-+N^++ANY_NEUTRAL=>NO2+N+ANY_NEUTRAL    ",&
   "NO2^-+N2^++ANY_NEUTRAL=>NO2+N2+ANY_NEUTRAL  ","NO2^-+O^++ANY_NEUTRAL=>NO2+O+ANY_NEUTRAL    ",&
   "NO2^-+O2^++ANY_NEUTRAL=>NO2+O2+ANY_NEUTRAL  ","NO2^-+NO^++ANY_NEUTRAL=>NO2+NO+ANY_NEUTRAL  ",&
   "NO2^-+N2O^++ANY_NEUTRAL=>NO2+N2O+ANY_NEUTRAL","NO2^-+NO2^++ANY_NEUTRAL=>NO2+NO2+ANY_NEUTRAL",&
   "NO3^-+NO^++ANY_NEUTRAL=>NO3+NO+ANY_NEUTRAL  ","NO3^-+NO2^++ANY_NEUTRAL=>NO3+NO2+ANY_NEUTRAL"/
  data bolsig_species(1:bolsig_species_max) &
  /"N2    ","N2(V1)","N2(V2)","N2(V3)","N2(V4)","N2(V5)","N2(V6)","N2(V7)","N2(V8)","N2(A3)","O2    ","O2(V1)","O2(V2)","O2(V3)",&
   "O2(V4)","O2(A1)","N     ","O     ","NO    ","O3    ","N2O   "/
  data bolsig_addsect(1:2,1:bolsig_addsect_max) &
  / 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9,11,12,11,13,11,14,11,15/
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!
! initialization
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_init()
  implicit none
  character(256) :: string
  integer :: i, j, k
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0a" // ") INIT:"
  if( lZDPlasKin_init ) call ZDPlasKin_stop("   ERROR: the ZDPlasKin library has been initialized")
  write(string,*) species_max
  write(*,"(2x,A)")  "species        ... " // trim(adjustl(string))
  write(string,*) reactions_max
  write(*,"(2x,A)")  "reactions      ... " // trim(adjustl(string))
  if(species_max<=0 .or. reactions_max<=0) call ZDPlasKin_stop("   ERROR: wrong preconfig data")
  write(*,"(2x,A,$)") "BOLSIG+ loader ... " // trim(adjustl(bolsigfile)) // " : "
  call ZDPlasKin_bolsig_Init(bolsigfile)
  do i = 1, bolsig_species_max
    call ZDPlasKin_bolsig_ReadCollisions(trim(bolsig_species(i)))
    j = bolsig_collisions_max
    call ZDPlasKin_bolsig_GetCollisions(k,bolsig_collisions_max)
    if(bolsig_collisions_max <= j) then
      write(*,*)
      call ZDPlasKin_stop("ERROR: wrong file or missing data " // &
                        "(" // trim(adjustl(bolsigfile)) // ": <" // trim(bolsig_species(i)) // ">).")
    endif
  enddo
  if(bolsig_species_max /= k) then
    write(*,*)
    call ZDPlasKin_stop("ERROR: internal error in BOLSIG+ loader")
  endif
  write(string,*) bolsig_species_max
  write(*,"(A,$)") trim(adjustl(string)) // " species & "
  write(string,*) bolsig_collisions_max
  write(*,"(A)")   trim(adjustl(string)) // " collisions"
  write(*,"(2x,A,$)") "species  link  ... "
  j = 0
  do i = 1, bolsig_species_max
    j = j + 1
    k = 1
    do while(k<=species_max .and. bolsig_species_index(i)<=0)
      call ZDPlasKin_bolsig_GetSpeciesName(string,i)
      if(trim(species_name(k)) == trim(string)) then
        bolsig_species_index(i) = k
      else
        k = k + 1
      endif
    enddo
    if(bolsig_species_index(i) <= 0) call ZDPlasKin_stop("cannot find species link for <" // trim(string) // ">")
  enddo
  write(string,*) j
  write(*,"(A)") trim(adjustl(string))
  write(*,"(2x,A,$)") "process  link  ... "
  i = 1
  j = 1
  do while(i<=reactions_max .and. j<=bolsig_rates_max)
    if(reaction_sign(i)(1:7) == "bolsig:") then
      k = 1
      do while(k<=bolsig_collisions_max .and. bolsig_pointer(j)<=0)
        call ZDPlasKin_bolsig_GetReactionName(string,k)
        if(trim(string) == trim(reaction_sign(i)(8:))) then
          bolsig_pointer(j) = k
        else
          k = k + 1
        endif
      enddo
      if(bolsig_pointer(j) <= 0) call ZDPlasKin_stop("cannot find processes link for <" // trim(reaction_sign(i)) // ">")
      j = j + 1
    endif
    i = i + 1
  enddo
  if(j <= bolsig_rates_max) then
    call ZDPlasKin_stop("internal error")
  else
    write(string,*) bolsig_rates_max
    write(*,"(A)") trim(adjustl(string))
  endif
  i = 0
  do while((1.0d0+10.0d0**(i-1)) /= 1.0d0)
    i = i - 1
  enddo
  mach_accur = 10.0d0**i
  mach_tiny  = sqrt( tiny(mach_tiny) )
  lZDPlasKin_init = .true.
  call ZDPlasKin_reset()
  write(*,"(A,/)") "ZDPlasKin INIT DONE"
  return
end subroutine ZDPlasKin_init
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using implicit solver dvode_f90
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep(time,dtime)
  use dvode_f90_m, only : dvode_f90
  implicit none
  double precision, intent(in)    ::  time
  double precision, intent(inout) :: dtime
  double precision, save :: densav(vode_neq) = 0.0d0, cfgsav(3) = 0.0d0
  double precision :: tout
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(time < tsav) vode_istate = 1
  tsav = time
  if(dtime > 0.0d0) then
    vode_itask = 1
    tout = time + dtime
    if(dtime < mach_accur*abs(tout)) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: dtime parameter is too small (subroutine ZDPlasKin_timestep)")
  else
    vode_itask = 2
    tout = ( 1.0d0 + mach_accur ) * time + mach_tiny
  endif
  dens_loc(1:species_max,0) = density(:)
  dens_loc(1:species_max,1) = 0.5d0 * ( density(:) + abs( density(:) ) )
  if(any(dens_loc(1:species_max,1) /= densav(1:species_max))) vode_istate = 1
  densav(1:species_max) = dens_loc(1:species_max,1)
  if(vode_istate /= 1 .and. any( abs(cfgsav(:)-ZDPlasKin_cfg(1:3)) > cfg_rtol*abs(cfgsav(:)+ZDPlasKin_cfg(1:3)) )) vode_istate = 1
  cfgsav(:) = ZDPlasKin_cfg(1:3)
  if( lgas_heating ) then
    if(ZDPlasKin_cfg(1) /= densav(species_max+1)) vode_istate = 1
    densav(species_max+1) = ZDPlasKin_cfg(1)
  endif
  call dvode_f90(ZDPlasKin_fex,vode_neq,densav,tsav,tout,vode_itask,vode_istate,vode_options,j_fcn=ZDPlasKin_jex)
  if(vode_istate < 0) then
    write(*,"(A,1pd11.4)") "Tgas   =", ZDPlasKin_cfg(1)
    write(*,"(A,1pd11.4)") "    EN =", ZDPlasKin_cfg(3)
    write(*,"(A,1pd11.4)") "    Te =", ZDPlasKin_cfg(4)
    call ZDPlasKin_stop("ZDPlasKin ERROR: DVODE solver issued an error (subroutine ZDPlasKin_timestep)")
  endif
  if( lgas_heating ) ZDPlasKin_cfg(1) = densav(species_max+1)
  density(:) = dens_loc(1:species_max,0) - dens_loc(1:species_max,1) + densav(1:species_max)
  if(dtime <= 0.0d0) dtime = tsav - time
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using explicit Euler method
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep_explicit(time,dtime,rtol_loc,atol_loc,switch_implicit)
  implicit none
  double precision, intent(in) ::  time, rtol_loc, atol_loc
  double precision, intent(inout) :: dtime
  double precision, optional, intent(in) :: switch_implicit
  double precision :: time_loc, time_end, dtime_loc, dtime_max
  logical, save :: lwarn = .true.
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(rtol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: rtol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  if(atol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: atol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  tsav     = time
  time_loc = 0.0d0
  time_end = 0.5d0 * ( dtime + abs(dtime) ) + mach_tiny
  do while(time_loc < time_end)
    dens_loc(1:species_max,0) = density(:)
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    dens_loc(:,1) = 0.5d0 * ( dens_loc(:,0) + abs( dens_loc(:,0) ) )
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,1),dens_loc(:,2))
    where(dens_loc(:,2) >= 0.0d0)
      dens_loc(:,3) = + dens_loc(:,2) /    ( rtol_loc * dens_loc(:,1) + atol_loc ) 
    elsewhere
      dens_loc(:,3) = - dens_loc(:,2) / min( rtol_loc * dens_loc(:,1) + atol_loc , dens_loc(:,1) + mach_tiny )
    endwhere
    dtime_loc = 1.0d0 / ( maxval( dens_loc(:,3) ) + mach_tiny )
    if(dtime > 0.0d0) then
      dtime_max = dtime - time_loc
      dtime_loc = min( dtime_loc , dtime_max )
      if( present(switch_implicit) ) then
        if(dtime_loc*switch_implicit < dtime_max) then
          if(lprint .and. lwarn) then
            write(*,"(A,/,A,1pd9.2,A)") "ZDPlasKin INFO: low efficiency of Euler method (subroutine ZDPlasKin_timestep_explicit)", &
                        "                ZDPlasKin_timestep subroutine will be used in similar conditions (", switch_implicit, ")"
            lwarn = .false.
          endif
          time_loc = tsav
          density(:) = dens_loc(1:species_max,0)
          if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0)
          call ZDPlasKin_timestep(time_loc,dtime_max)
          return
        endif
      endif
    else
      dtime = dtime_loc
    endif
    time_loc = time_loc + dtime_loc
    tsav     = time     +  time_loc
    density(:) = dens_loc(1:species_max,0) + dtime_loc * dens_loc(1:species_max,2)
    if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0) + dtime_loc * dens_loc(species_max+1,2)
  enddo
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! update BOLSIG+ solution and get electron parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_bolsig_rates(lbolsig_force)
  implicit none
  logical, optional, intent(in) :: lbolsig_force
  logical :: lforce
  integer :: i, j, k
  integer, save :: bolsig_points_max
  logical, save :: lfirst = .true., leecol = .true.
  double precision :: error, density_loc, cfg_loc(6+bolsig_species_max)
  double precision, save :: low_density_limit = bolsig_rtol, bolsig_mesh_a, bolsig_mesh_b
  double precision, save, allocatable :: bolsig_cfg(:,:), bolsig_reslt(:,:)
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    bolsig_mesh_a = 1.0d0 / log( 1.0d0 + bolsig_rtol )
    bolsig_mesh_b = bolsig_mesh_a * log( bolsig_field_min ) - 0.5d0
    bolsig_points_max = int( bolsig_mesh_a * log( bolsig_field_max ) - bolsig_mesh_b )
    allocate(bolsig_rates(bolsig_collisions_max), &
             bolsig_cfg(6+bolsig_species_max,0:bolsig_points_max), &
             bolsig_reslt(10+bolsig_collisions_max,0:bolsig_points_max),stat=i)
    if(i /= 0) call ZDPlasKin_stop("ZDPlasKin ERROR: memory allocation error (subroutine ZDPlasKin_bolsig_rates)")
    bolsig_cfg(:,:) = 0.0d0
    lfirst = .false.
  endif
  if( present(lbolsig_force) ) then
    lforce = lbolsig_force
  else
    lforce = .false.
  endif
  if(ZDPlasKin_cfg(1) <= 0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
  if(.not. lbolsig_Maxwell_EEDF ) then
    if(ZDPlasKin_cfg(2) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_bolsig_rates)")
    if(ZDPlasKin_cfg(3) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_bolsig_rates)")
    ZDPlasKin_cfg(4) = 0.0d0
  else
    if(ZDPlasKin_cfg(4) <= 0.0d0) then
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELECTRON_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
    elseif(lprint .and. ZDPlasKin_cfg(4) < ZDPlasKin_cfg(1)) then
      write(*,"(A)") "ZDPlasKin INFO: ELECTRON_TEMPERATURE is below GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)"
    endif
    ZDPlasKin_cfg(2:3) = 0.0d0
  endif
  density_loc = 0.5d0 * ( sum(density(bolsig_species_index(:))) + sum(abs(density(bolsig_species_index(:)))) )
  if(density_loc <= mach_tiny) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined densities configured for BOLSIG+ solver " // &
                                                                     "(subroutine ZDPlasKin_bolsig_rates)")
  elseif( lprint ) then
    call ZDPlasKin_get_density_total(ALL_NEUTRAL=error)
    error = abs( 1.0d0 - density_loc / error )
    if(error > low_density_limit) then
      write(*,"(A,1pd9.2)") "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds", error
      low_density_limit = sqrt( low_density_limit )
    endif
  endif
  cfg_loc(1:4) = ZDPlasKin_cfg(1:4)
  cfg_loc(2)   = cfg_loc(2) * 1.0d-6
  cfg_loc(6)   = 0.5d0 * ( density(species_electrons) + abs(density(species_electrons)) )
  cfg_loc(5)   = cfg_loc(6) * 1.0d6
  cfg_loc(6)   = cfg_loc(6) / density_loc
  if(cfg_loc(6) < bolsig_eecol_frac) then
    cfg_loc(6) = 0.0d0
  elseif(lprint .and. leecol) then
    write(*,"(A)") "ZDPlasKin INFO: set electron-electron collisions ON ..."
    leecol = .false.
  endif
  cfg_loc(7:) = 0.5d0 * ( density(bolsig_species_index(:)) + abs(density(bolsig_species_index(:))) ) / density_loc
  if(lbolsig_Maxwell_EEDF .or. bolsig_points_max==0) then
    lforce = .true.
    i = 0
  else
    error = min( max( cfg_loc(3) , bolsig_field_min ) , bolsig_field_max )
    i = int( bolsig_mesh_a * log( error ) - bolsig_mesh_b )
    i = max(0,min(i,bolsig_points_max))
  endif
  if( lforce ) then
    error = 2.0d0
  else
    if(lbolsig_ignore_gas_temp .and. bolsig_cfg(1,i)>0.0d0) then
      cfg_loc(1) = bolsig_cfg(1,i)
      error = 0.0d0
    else
      error = abs( ( cfg_loc(1) - bolsig_cfg(1,i) ) / ( 0.5d0 * ( cfg_loc(1) + bolsig_cfg(1,i) ) + mach_tiny) ) / bolsig_rtol_half
    endif
    if(error <= 1.0d0) then
      error = abs( ( cfg_loc(2) - bolsig_cfg(2,i) ) / ( 0.5d0 * ( cfg_loc(2) + bolsig_cfg(2,i) ) + mach_tiny) ) / bolsig_rtol_half
      if(error <= 1.0d0) then
        error = abs( ( cfg_loc(3) - bolsig_cfg(3,i) ) / ( 0.5d0 * ( cfg_loc(3) + bolsig_cfg(3,i) ) + mach_tiny) ) / bolsig_rtol
        if(error <= 1.0d0) then
          error = abs( ( max(cfg_loc(6),bolsig_eecol_frac) - max(bolsig_cfg(6,i),bolsig_eecol_frac) ) &
           / ( 0.5d0 * ( max(cfg_loc(6),bolsig_eecol_frac) + max(bolsig_cfg(6,i),bolsig_eecol_frac) ) + mach_tiny) ) &
           / bolsig_rtol_half
          if(error <= 1.0d0) error = maxval( abs( cfg_loc(7:) - bolsig_cfg(7:,i) ) ) &
                                     / ( 0.5d0 * maxval( cfg_loc(7:) + bolsig_cfg(7:,i) ) ) / bolsig_rtol
        endif
      endif
    endif
  endif
  if(error > 1.0d0) then
    j = 6 + bolsig_species_max
    k = 10 + bolsig_collisions_max
    bolsig_cfg(:,i) = cfg_loc(:)
    call ZDPlasKin_bolsig_SolveBoltzmann(j,bolsig_cfg(1:j,i),k,bolsig_reslt(1:k,i))
    if(.not. lbolsig_Maxwell_EEDF) then
      bolsig_reslt(2,i) = bolsig_reslt(2, i) * eV_to_K / 1.5d0
    else
      bolsig_reslt(2,i) = cfg_loc(4)
    endif
    bolsig_reslt(3, i) = bolsig_reslt(3, i) * 1.0d-2
    bolsig_reslt(4, i) = bolsig_reslt(4, i) * 1.0d-2 / density_loc
    bolsig_reslt(5, i) = bolsig_reslt(5, i) * 1.0d-2
    bolsig_reslt(6, i) = bolsig_reslt(6, i) * 1.0d-2
    bolsig_reslt(7:,i) = bolsig_reslt(7:,i) * 1.0d6
  endif
  ZDPlasKin_cfg(3:12) = bolsig_reslt(:10,i)
  bolsig_rates(:)     = bolsig_reslt(11:,i)
  return
end subroutine ZDPlasKin_bolsig_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get index of species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_species_index(str,i)
  implicit none
  character(*), intent(in ) :: str
  integer,      intent(out) :: i
  character(species_length) :: string
  integer :: j, istr
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  string = trim(adjustl(str))
  istr   = len_trim(string)
  do i = 1, istr
    j = iachar(string(i:i))
    if(j>=97 .and. j<=122) string(i:i) = achar(j-32)
  enddo
  i = 0
  j = 0
  do while(i==0 .and. j<species_max)
    j = j + 1
    if(string(1:istr) == trim(species_name(j))) i = j
  enddo
  if(i <= 0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: cannot identify species <"//trim(str)//"> (subroutine ZDPlasKin_get_species_index)")
  return
end subroutine ZDPlasKin_get_species_index
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get density for species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in) :: string
  logical, optional, intent(in) :: LDENS_CONST
  double precision, optional, intent(in) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present(DENS) ) density(i) = DENS
  if( present(LDENS_CONST) ) then
    density_constant(i) = LDENS_CONST
    ldensity_constant   = any( density_constant(:) )
  endif
  return
end subroutine ZDPlasKin_set_density
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in ) :: string
  logical, optional, intent(out) :: LDENS_CONST
  double precision, optional, intent(out) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present( DENS)       )  DENS       = density(i)
  if( present(LDENS_CONST) ) LDENS_CONST = density_constant(i)
  return
end subroutine ZDPlasKin_get_density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get total densities
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_density_total(ALL_SPECIES,ALL_NEUTRAL,ALL_ION_POSITIVE,ALL_ION_NEGATIVE,ALL_CHARGE)
  double precision, optional, intent(out) :: ALL_SPECIES, ALL_NEUTRAL, ALL_ION_POSITIVE, ALL_ION_NEGATIVE, ALL_CHARGE
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ALL_SPECIES)      ) ALL_SPECIES      = sum(density(:))
  if( present(ALL_NEUTRAL)      ) ALL_NEUTRAL      = sum(density(:), mask = species_charge(:)==0)
  if( present(ALL_ION_POSITIVE) ) ALL_ION_POSITIVE = sum(density(:), mask = species_charge(:)>0)
  if( present(ALL_ION_NEGATIVE) ) ALL_ION_NEGATIVE = sum(density(:), mask = species_charge(:)<0) - density(species_electrons)
  if( present(ALL_CHARGE)       ) ALL_CHARGE       = sum(density(:) * dble(species_charge(:)))
  return
end subroutine ZDPlasKin_get_density_total
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get species source terms & reaction rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_rates(SOURCE_TERMS,REACTION_RATES,SOURCE_TERMS_MATRIX,MEAN_DENSITY, &
                               MEAN_SOURCE_TERMS,MEAN_REACTION_RATES,MEAN_SOURCE_TERMS_MATRIX)
  double precision, optional, intent(out) :: SOURCE_TERMS(species_max), REACTION_RATES(reactions_max), &
                                             SOURCE_TERMS_MATRIX(species_max,reactions_max), MEAN_DENSITY(species_max), &
                                             MEAN_SOURCE_TERMS(species_max), MEAN_REACTION_RATES(reactions_max), &
                                             MEAN_SOURCE_TERMS_MATRIX(species_max,reactions_max)
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if(present(SOURCE_TERMS) .or. present(REACTION_RATES) .or. present(SOURCE_TERMS_MATRIX)) then
    dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
    if( present(SOURCE_TERMS)                 ) SOURCE_TERMS(:)   = dens_loc(1:species_max,1)
    if( present(REACTION_RATES)               ) REACTION_RATES(:) = rrt(:)
    if( present(SOURCE_TERMS_MATRIX)          ) call ZDPlasKin_reac_source_matrix(rrt(:),SOURCE_TERMS_MATRIX(:,:))
  endif
  if(present(MEAN_DENSITY)        .or. present(MEAN_SOURCE_TERMS) .or. &
     present(MEAN_REACTION_RATES) .or. present(MEAN_SOURCE_TERMS_MATRIX)) then
    if( lstat_accum ) then
      if(stat_time > 0.0d0) then
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = stat_dens(:) / stat_time
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = stat_src(:)  / stat_time
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = stat_rrt(:)  / stat_time
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) then
          call ZDPlasKin_reac_source_matrix(stat_rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
          MEAN_SOURCE_TERMS_MATRIX(:,:)  =  MEAN_SOURCE_TERMS_MATRIX(:,:) / stat_time
        endif
      else
        dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
        if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
        call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = density(:)
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = dens_loc(1:species_max,1)
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = rrt(:)
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) call ZDPlasKin_reac_source_matrix(rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
      endif
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_get_rates)")
    endif
  endif
  return
end subroutine ZDPlasKin_get_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set config
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)
  use dvode_f90_m, only : set_intermediate_opts
  implicit none
  logical, optional, intent(in) :: SILENCE_MODE, STAT_ACCUM, QTPLASKIN_SAVE, BOLSIG_IGNORE_GAS_TEMPERATURE
  double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC
  integer :: i
  logical, save :: lfirst = .true.
  integer, save :: bounded_components(vode_neq)
  double precision :: atol_loc, rtol_loc
  double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    do i = 1, vode_neq
      bounded_components(i) = i
    enddo
    lfirst = .false.
  endif
  if( present(SILENCE_MODE) ) lprint = ( .not. SILENCE_MODE )
  if( present(BOLSIG_EE_FRAC) ) bolsig_eecol_frac = 0.5d0 * ( BOLSIG_EE_FRAC + abs(BOLSIG_EE_FRAC) )
  if( present(BOLSIG_IGNORE_GAS_TEMPERATURE) ) lbolsig_ignore_gas_temp = BOLSIG_IGNORE_GAS_TEMPERATURE
  if( present(STAT_ACCUM) ) then
    if( lprint ) then
      if(lstat_accum .neqv. STAT_ACCUM) then
        if( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition OFF ..."
        endif
      elseif( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: reset statistic acquisition data ..."
      endif
    endif
    stat_dens(:) = 0.0d0
    stat_src(:)  = 0.0d0
    stat_rrt(:)  = 0.0d0
    stat_time    = 0.0d0
    lstat_accum  = STAT_ACCUM
  endif
  if( present(QTPLASKIN_SAVE) ) then
    if( lprint ) then
      if(lqtplaskin .neqv. QTPLASKIN_SAVE) then
        if( QTPLASKIN_SAVE ) then
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format OFF ..."
        endif
      endif
    endif
    lqtplaskin = QTPLASKIN_SAVE
  endif
  if( present(ATOL) ) then
    atol_loc = ATOL
  else
    atol_loc = atol_save
  endif
  if( present(RTOL) ) then
    rtol_loc = RTOL
  else
    rtol_loc = rtol_save
  endif
  if(min(atol_loc,rtol_loc)<0.0d0 .or. max(atol_loc,rtol_loc)<=0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ATOL/RTOL (ZDPlasKin_set_config)")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,"(2(A,1pd9.2),A)") "ZDPlasKin INFO: set accuracy", atol_save, " (absolute) &", rtol_save, " (relative)"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif
  return
end subroutine ZDPlasKin_set_config
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get conditions
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,GAS_HEATING,SPEC_HEAT_RATIO,HEAT_SOURCE,SOFT_RESET)
  implicit none
  double precision, optional, intent(in) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, ELEC_TEMPERATURE, &
                                            SPEC_HEAT_RATIO, HEAT_SOURCE
  logical,          optional, intent(in) :: GAS_HEATING, SOFT_RESET
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(GAS_TEMPERATURE) ) then
    if(GAS_TEMPERATURE <= 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(1) = GAS_TEMPERATURE
  endif
  if( present(REDUCED_FREQUENCY) ) then
    if(REDUCED_FREQUENCY < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(2) = REDUCED_FREQUENCY
  endif
  if( present(REDUCED_FIELD) ) then
    if(REDUCED_FIELD < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(3) = REDUCED_FIELD
  endif
  if( present(SOFT_RESET) ) then
    if( SOFT_RESET ) vode_istate = 1
  endif
  if( present(ELEC_TEMPERATURE) ) then
    if(ELEC_TEMPERATURE < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELEC_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    if(ELEC_TEMPERATURE > 0.0d0) then
      lbolsig_Maxwell_EEDF = .true.
    else
      lbolsig_Maxwell_EEDF = .false.
    endif
    ZDPlasKin_cfg(4) = ELEC_TEMPERATURE
  endif
  if( present(GAS_HEATING) ) then
    if(lgas_heating .neqv. GAS_HEATING) then
      if( GAS_HEATING ) then
        if(present(SPEC_HEAT_RATIO) .or. ZDPlasKin_cfg(13)>0.0d0) then
          if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating ON ..."
        else
          ZDPlasKin_cfg(13) = 2.0d0/3.0d0
          if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set gas heating ON; specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
        endif
      else
        if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating OFF ..."
      endif
      lgas_heating = GAS_HEATING
    endif
  endif
  if( present(SPEC_HEAT_RATIO) ) then
    if(SPEC_HEAT_RATIO > 1.0d0) then
      ZDPlasKin_cfg(13) = SPEC_HEAT_RATIO - 1.0d0
      if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong value of SPEC_HEAT_RATIO (subroutine ZDPlasKin_set_conditions)")
    endif
  endif
  if( present(HEAT_SOURCE) ) then
    ZDPlasKin_cfg(14) = HEAT_SOURCE
    if( lprint ) write(*,"(A,1pd9.2,A)") "ZDPlasKin INFO: set heat source =", ZDPlasKin_cfg(14), " W/cm3"
  endif
end subroutine ZDPlasKin_set_conditions
subroutine ZDPlasKin_get_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,ELEC_DRIFT_VELOCITY,ELEC_DIFF_COEFF,ELEC_MOBILITY_N, &
                                    ELEC_MU_EPS_N,ELEC_DIFF_EPS_N,ELEC_FREQUENCY_N, &
                                    ELEC_POWER_N,ELEC_POWER_ELASTIC_N,ELEC_POWER_INELASTIC_N,ELEC_EEDF)
  implicit none
  double precision, optional, intent(out) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, &
                                             ELEC_TEMPERATURE, ELEC_DRIFT_VELOCITY, ELEC_DIFF_COEFF, ELEC_MOBILITY_N, &
                                             ELEC_MU_EPS_N, ELEC_DIFF_EPS_N, ELEC_FREQUENCY_N, &
                                             ELEC_POWER_N, ELEC_POWER_ELASTIC_N, ELEC_POWER_INELASTIC_N
  double precision, optional, dimension(:,:), intent(out) :: ELEC_EEDF
  integer :: i
  double precision :: x,y
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ELEC_EEDF) ) then
    call ZDPlasKin_bolsig_rates(lbolsig_force=.true.)
  else
    call ZDPlasKin_bolsig_rates()
  endif
  if( present(GAS_TEMPERATURE)        ) GAS_TEMPERATURE        = ZDPlasKin_cfg(1)
  if( present(REDUCED_FREQUENCY)      ) REDUCED_FREQUENCY      = ZDPlasKin_cfg(2)
  if( present(REDUCED_FIELD)          ) REDUCED_FIELD          = ZDPlasKin_cfg(3)
  if( present(ELEC_TEMPERATURE)       ) ELEC_TEMPERATURE       = ZDPlasKin_cfg(4)
  if( present(ELEC_DRIFT_VELOCITY)    ) ELEC_DRIFT_VELOCITY    = ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
  if( present(ELEC_DIFF_COEFF)        ) ELEC_DIFF_COEFF        = ZDPlasKin_cfg(6)
  if( present(ELEC_MOBILITY_N)        ) ELEC_MOBILITY_N        = ZDPlasKin_cfg(5)
  if( present(ELEC_MU_EPS_N)          ) ELEC_MU_EPS_N          = ZDPlasKin_cfg(7)
  if( present(ELEC_DIFF_EPS_N)        ) ELEC_DIFF_EPS_N        = ZDPlasKin_cfg(8)
  if( present(ELEC_FREQUENCY_N)       ) ELEC_FREQUENCY_N       = ZDPlasKin_cfg(9)
  if( present(ELEC_POWER_N)           ) ELEC_POWER_N           = ZDPlasKin_cfg(10)
  if( present(ELEC_POWER_ELASTIC_N)   ) ELEC_POWER_ELASTIC_N   = ZDPlasKin_cfg(11)
  if( present(ELEC_POWER_INELASTIC_N) ) ELEC_POWER_INELASTIC_N = ZDPlasKin_cfg(12)
  if( present(ELEC_EEDF) ) then
    ELEC_EEDF = 0d0
  	 if( size(ELEC_EEDF,dim=1) < 2 ) then
      if(lprint) write(*,"(A)") &
  	     "ZDPlasKin WARNING: insufficient first dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
  	  else
  		y = 1.0d0
  		do i = 1, size(ELEC_EEDF,dim=2)
  		  call ZDPlasKin_bolsig_GetEEDF(i,x,y)
  		  if( x >= 0d0 .and. y > 0d0) then
  			ELEC_EEDF(1,i) = x
  			ELEC_EEDF(2,i) = y
  		  else
  			exit
  		  endif
  		enddo
  		if(lprint .and. y>0d0) write(*,"(A)") &
  		  "ZDPlasKin WARNING: insufficient second dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
     endif
  endif
end subroutine ZDPlasKin_get_conditions
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reset
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reset()
  implicit none
  vode_istate         =  1
  density(:)          =  0.0d0
  ZDPlasKin_cfg(:)    =  0.0d0
  ldensity_constant   = .false.
  density_constant(:) = .false.
  lreaction_block(:)  = .false.
  lprint              = .true.
  lstat_accum         = .false.
  lqtplaskin          = .false.
  lgas_heating        = .false.
  bolsig_eecol_frac       = bolsig_eecol_frac_def
  lbolsig_ignore_gas_temp = .false.
  lbolsig_Maxwell_EEDF    = .false.
  write(*,"(A)") "ZDPlasKin INFO: reset data and configuration"
  call ZDPlasKin_set_config(ATOL=vode_atol,RTOL=vode_rtol)
  return
end subroutine ZDPlasKin_reset
!-----------------------------------------------------------------------------------------------------------------------------------
!
! stop
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_stop(string)
  implicit none
  character(*), intent(in) :: string
  if(string /= "") write(*,"(A)") trim(string)
  write(*,"(A,$)") "PRESS ENTER TO EXIT ... "
  read(*,*)
  stop
end subroutine ZDPlasKin_stop
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data to file
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_file(FILE_SPECIES,FILE_REACTIONS,FILE_SOURCE_MATRIX,FILE_UNIT)
  implicit none
  character(*), optional, intent(in) :: FILE_SPECIES, FILE_REACTIONS, FILE_SOURCE_MATRIX
  integer, optional, intent(in) :: FILE_UNIT
  logical :: lerror
  integer :: i
  if( present(FILE_UNIT) ) ifile_unit = FILE_UNIT
  if( present(FILE_SPECIES) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SPECIES)),action="write",err=100)
    do i = 1, species_max
      write(ifile_unit,111,err=100) i, species_name(i)
    enddo
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SPECIES)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
111 format(i2,1x,A9)
  endif
  if( present(FILE_REACTIONS) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_REACTIONS)),action="write",err=200)
    do i = 1, reactions_max
      write(ifile_unit,211,err=200) i, reaction_sign(i)
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_REACTIONS)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
211 format(i3,1x,A44)
  endif
  if( present(FILE_SOURCE_MATRIX) ) then
    if( lstat_accum ) then
      call ZDPlasKin_reac_source_matrix(stat_rrt(:),mrtm(:,:))
      if(stat_time > 0.0d0) mrtm(:,:) = mrtm(:,:) / stat_time
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_write_file)")
    endif
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SOURCE_MATRIX)),action="write",err=300)
    write(ifile_unit,311,err=300) ( i, i = 1, species_max )
    write(ifile_unit,312,err=300) "N", "reaction", ( trim(species_name(i)), i = 1, species_max )
    do i = 1, reactions_max
      write(ifile_unit,313,err=300) i, reaction_sign(i), mrtm(:,i)
    enddo
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SOURCE_MATRIX)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
311 format(491x,53(1x,i9))
312 format(A3,1x,A44,1x,53(1x,A9))
313 format(i3,1x,A44,1x,53(1x,1pd9.2))
  endif
  return
end subroutine ZDPlasKin_write_file
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data in qtplaskin format
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_qtplaskin(time,LFORCE_WRITE)
  implicit none
  double precision, intent(in) :: time
  logical, optional, intent(in) :: LFORCE_WRITE
  integer, parameter :: idef_data = 5
  character(24), parameter :: qtplaskin_names(idef_data) = (/ "Reduced field [Td]      ", "Gas temperature [K]     ", &
                                  "Electron temperature [K]", "Current density [A/cm2] ", "Power density [W/cm3]   " /)
  double precision, save :: densav(0:species_max,2) = -huge(densav)
  double precision :: rtol, cond(idef_data)
  logical, save :: lfirst = .true.
  logical :: lerror
  integer, save :: iuser_data = 0
  integer :: i
  if( time < densav(0,1) ) lfirst = .true.
  if( lfirst ) then
    call ZDPlasKin_write_file(FILE_SPECIES="qt_species_list.txt",FILE_REACTIONS="qt_reactions_list.txt")
    if( allocated(qtplaskin_user_data) ) then
      iuser_data = size(qtplaskin_user_data)
      iuser_data = min(iuser_data,90)
      if( iuser_data > 0 ) then
        if( allocated(qtplaskin_user_names) ) then
          if( size(qtplaskin_user_names) /= iuser_data ) deallocate(qtplaskin_user_names)
        endif
        if( .not. allocated(qtplaskin_user_names) ) then
          allocate(qtplaskin_user_names(iuser_data))
          do i = 1, iuser_data
            write(qtplaskin_user_names(i),"(A,i2.2)") "user defined #", i
          enddo
        endif
      endif
    endif
    lerror = .true.
    open(ifile_unit,file="qt_conditions_list.txt",action="write",err=100)
    do i = 1, idef_data
      write(ifile_unit,"(i3,1x,A)",err=100) i, trim(adjustl(qtplaskin_names(i)))
    enddo
    if( iuser_data > 0 ) then
      do i = 1, iuser_data
        write(ifile_unit,"(i3,1x,A)",err=100) (i+idef_data), trim(adjustl(qtplaskin_user_names(i)))
      enddo
    endif
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions_list.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    rrt(:) = 1.0d0
    call ZDPlasKin_reac_source_matrix(rrt(:),mrtm(:,:))
    open(ifile_unit,file="qt_matrix.txt",action="write",err=200)
    do i = 1, species_max
      write(ifile_unit,"(246(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,53(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_densities.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_conditions.txt",action="write",err=400)
    write(ifile_unit,"(1x,A12,$)",err=400) "Time_s"
    do i = 1, idef_data + iuser_data
      write(ifile_unit,"(11x,i2.2,$)",err=400) i
    enddo
    write(ifile_unit,*,err=400)
    lerror = .false.
400 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_rates.txt",action="write",err=500)
    write(ifile_unit,"(1x,A12,246(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 6 )
  if( ( time - densav(0,1) ) >= rtol .or. lfirst ) then
    densav(0,2) = time
    densav(1:species_max,2) = density(:)
    where( densav(:,2) < 1.0d-99 ) densav(:,2) = 0.0d0
    if( time > 2.0d0 * densav(0,1) ) then
      rtol = huge(rtol)
    else
      rtol = maxval( abs(densav(1:,1)-densav(1:,2)) / ( abs(densav(1:,1)+densav(1:,2))/2.0d0 + qtplaskin_atol ) )
    endif
    if( rtol > qtplaskin_rtol .or. lfirst ) then
      open(ifile_unit,file="qt_densities.txt",access="append")
      write(ifile_unit,"(1pe15.6,53(1pe13.4))") densav(0,2), densav(1:,2)
      close(ifile_unit)
      open(ifile_unit,file="qt_conditions.txt",access="append")
      cond(1) = ZDPlasKin_cfg(3)
      cond(2) = ZDPlasKin_cfg(1)
      cond(3) = ZDPlasKin_cfg(4)
      cond(4) = q_elem * density(species_electrons) * ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
      call ZDPlasKin_get_density_total(ALL_NEUTRAL=cond(5))
      cond(5) = cond(4) * cond(5) * ZDPlasKin_cfg(3) * 1.0d-17
      where( abs(cond(:)) < 1.0d-99 ) cond(:) = 0.0d0
      write(ifile_unit,"(6(1pe13.4),$)") densav(0,2), cond(:)
      if( iuser_data > 0 ) then
        where( abs(qtplaskin_user_data(1:iuser_data)) < 1.0d-99 ) qtplaskin_user_data(1:iuser_data) = 0.0d0
        write(ifile_unit,"(90(1pe13.4))") qtplaskin_user_data(1:iuser_data)
      else
        write(ifile_unit,*)
      endif
      close(ifile_unit)
      call ZDPlasKin_get_rates(REACTION_RATES=rrt_loc)
      where( abs(rrt_loc(:)) < 1.0d-99 ) rrt_loc(:) = 0.0d0
      open(ifile_unit,file="qt_rates.txt",access="append")
      write(ifile_unit,"(247(1pe13.4))") densav(0,2), rrt_loc(:)
      close(ifile_unit)
      densav(:,1) = densav(:,2)
    endif
  endif
  lfirst = .false.
  lqtplaskin_first = .false.
  return
end subroutine ZDPlasKin_write_qtplaskin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction sensitivity acquisition
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_source_matrix(reac_rate_local,reac_source_local)
  implicit none
  double precision, intent(in)  :: reac_rate_local(reactions_max)
  double precision, intent(out) :: reac_source_local(species_max,reactions_max)
  reac_source_local(:,:) = 0.0d0
  reac_source_local(01,001) = - reac_rate_local(001) 
  reac_source_local(02,001) = + reac_rate_local(001) 
  reac_source_local(01,002) = - reac_rate_local(002) 
  reac_source_local(02,002) = + reac_rate_local(002) 
  reac_source_local(01,003) = - reac_rate_local(003) 
  reac_source_local(03,003) = + reac_rate_local(003) 
  reac_source_local(01,004) = - reac_rate_local(004) 
  reac_source_local(04,004) = + reac_rate_local(004) 
  reac_source_local(01,005) = - reac_rate_local(005) 
  reac_source_local(05,005) = + reac_rate_local(005) 
  reac_source_local(01,006) = - reac_rate_local(006) 
  reac_source_local(06,006) = + reac_rate_local(006) 
  reac_source_local(01,007) = - reac_rate_local(007) 
  reac_source_local(07,007) = + reac_rate_local(007) 
  reac_source_local(01,008) = - reac_rate_local(008) 
  reac_source_local(08,008) = + reac_rate_local(008) 
  reac_source_local(01,009) = - reac_rate_local(009) 
  reac_source_local(09,009) = + reac_rate_local(009) 
  reac_source_local(01,010) = + reac_rate_local(010) 
  reac_source_local(02,010) = - reac_rate_local(010) 
  reac_source_local(01,011) = + reac_rate_local(011) 
  reac_source_local(03,011) = - reac_rate_local(011) 
  reac_source_local(01,012) = + reac_rate_local(012) 
  reac_source_local(04,012) = - reac_rate_local(012) 
  reac_source_local(01,013) = + reac_rate_local(013) 
  reac_source_local(05,013) = - reac_rate_local(013) 
  reac_source_local(01,014) = + reac_rate_local(014) 
  reac_source_local(06,014) = - reac_rate_local(014) 
  reac_source_local(01,015) = + reac_rate_local(015) 
  reac_source_local(07,015) = - reac_rate_local(015) 
  reac_source_local(01,016) = + reac_rate_local(016) 
  reac_source_local(08,016) = - reac_rate_local(016) 
  reac_source_local(01,017) = + reac_rate_local(017) 
  reac_source_local(09,017) = - reac_rate_local(017) 
  reac_source_local(01,018) = + reac_rate_local(018) 
  reac_source_local(02,018) = - reac_rate_local(018) 
  reac_source_local(02,019) = + reac_rate_local(019) 
  reac_source_local(03,019) = - reac_rate_local(019) 
  reac_source_local(03,020) = + reac_rate_local(020) 
  reac_source_local(04,020) = - reac_rate_local(020) 
  reac_source_local(04,021) = + reac_rate_local(021) 
  reac_source_local(05,021) = - reac_rate_local(021) 
  reac_source_local(05,022) = + reac_rate_local(022) 
  reac_source_local(06,022) = - reac_rate_local(022) 
  reac_source_local(06,023) = + reac_rate_local(023) 
  reac_source_local(07,023) = - reac_rate_local(023) 
  reac_source_local(07,024) = + reac_rate_local(024) 
  reac_source_local(08,024) = - reac_rate_local(024) 
  reac_source_local(08,025) = + reac_rate_local(025) 
  reac_source_local(09,025) = - reac_rate_local(025) 
  reac_source_local(01,026) = + reac_rate_local(026) 
  reac_source_local(02,026) = - reac_rate_local(026) 
  reac_source_local(02,027) = + reac_rate_local(027) 
  reac_source_local(03,027) = - reac_rate_local(027) 
  reac_source_local(03,028) = + reac_rate_local(028) 
  reac_source_local(04,028) = - reac_rate_local(028) 
  reac_source_local(04,029) = + reac_rate_local(029) 
  reac_source_local(05,029) = - reac_rate_local(029) 
  reac_source_local(05,030) = + reac_rate_local(030) 
  reac_source_local(06,030) = - reac_rate_local(030) 
  reac_source_local(06,031) = + reac_rate_local(031) 
  reac_source_local(07,031) = - reac_rate_local(031) 
  reac_source_local(07,032) = + reac_rate_local(032) 
  reac_source_local(08,032) = - reac_rate_local(032) 
  reac_source_local(08,033) = + reac_rate_local(033) 
  reac_source_local(09,033) = - reac_rate_local(033) 
  reac_source_local(01,034) = - reac_rate_local(034) 
  reac_source_local(10,034) = + reac_rate_local(034) 
  reac_source_local(01,035) = - reac_rate_local(035) 
  reac_source_local(10,035) = + reac_rate_local(035) 
  reac_source_local(01,036) = - reac_rate_local(036) 
  reac_source_local(10,036) = + reac_rate_local(036) 
  reac_source_local(01,037) = - reac_rate_local(037) 
  reac_source_local(11,037) = + reac_rate_local(037) 
  reac_source_local(01,038) = - reac_rate_local(038) 
  reac_source_local(11,038) = + reac_rate_local(038) 
  reac_source_local(01,039) = - reac_rate_local(039) 
  reac_source_local(11,039) = + reac_rate_local(039) 
  reac_source_local(01,040) = - reac_rate_local(040) 
  reac_source_local(12,040) = + reac_rate_local(040) 
  reac_source_local(01,041) = - reac_rate_local(041) 
  reac_source_local(12,041) = + reac_rate_local(041) 
  reac_source_local(01,042) = - reac_rate_local(042) 
  reac_source_local(12,042) = + reac_rate_local(042) 
  reac_source_local(01,043) = - reac_rate_local(043) 
  reac_source_local(13,043) = + reac_rate_local(043) 
  reac_source_local(01,044) = - reac_rate_local(044) 
  reac_source_local(13,044) = + reac_rate_local(044) 
  reac_source_local(01,045) = - reac_rate_local(045) 
  reac_source_local(14,045) = + reac_rate_local(045) 
  reac_source_local(15,045) = + reac_rate_local(045) 
  reac_source_local(40,046) = - reac_rate_local(046) 
  reac_source_local(45,046) = + reac_rate_local(046) 
  reac_source_local(53,046) = + reac_rate_local(046) 
  reac_source_local(29,047) = + reac_rate_local(047) 
  reac_source_local(40,047) = + reac_rate_local(047) 
  reac_source_local(47,047) = - reac_rate_local(047) 
  reac_source_local(53,047) = - reac_rate_local(047) 
  reac_source_local(36,048) = + reac_rate_local(048) 
  reac_source_local(40,048) = + reac_rate_local(048) 
  reac_source_local(42,048) = - reac_rate_local(048) 
  reac_source_local(53,048) = - reac_rate_local(048) 
  reac_source_local(40,049) = - reac_rate_local(049) 
  reac_source_local(48,049) = + reac_rate_local(049) 
  reac_source_local(53,049) = - reac_rate_local(049) 
  reac_source_local(14,050) = - reac_rate_local(050) 
  reac_source_local(36,050) = - reac_rate_local(050) 
  reac_source_local(40,050) = + reac_rate_local(050) 
  reac_source_local(53,050) = + reac_rate_local(050) 
  reac_source_local(36,051) = - reac_rate_local(051) 
  reac_source_local(40,051) = - reac_rate_local(051) 
  reac_source_local(42,051) = + reac_rate_local(051) 
  reac_source_local(53,051) = + reac_rate_local(051) 
  reac_source_local(14,052) = - reac_rate_local(052) 
  reac_source_local(37,052) = - reac_rate_local(052) 
  reac_source_local(42,052) = + reac_rate_local(052) 
  reac_source_local(53,052) = + reac_rate_local(052) 
  reac_source_local(14,053) = - reac_rate_local(053) 
  reac_source_local(21,053) = + reac_rate_local(053) 
  reac_source_local(38,053) = - reac_rate_local(053) 
  reac_source_local(40,053) = + reac_rate_local(053) 
  reac_source_local(53,053) = + reac_rate_local(053) 
  reac_source_local(01,054) = + reac_rate_local(054) 
  reac_source_local(14,054) = - reac_rate_local(054) 
  reac_source_local(40,054) = + reac_rate_local(054) 
  reac_source_local(49,054) = - reac_rate_local(054) 
  reac_source_local(53,054) = + reac_rate_local(054) 
  reac_source_local(14,055) = - reac_rate_local(055) 
  reac_source_local(40,055) = + reac_rate_local(055) * 2.d0
  reac_source_local(50,055) = - reac_rate_local(055) 
  reac_source_local(53,055) = + reac_rate_local(055) 
  reac_source_local(14,056) = - reac_rate_local(056) 
  reac_source_local(40,056) = + reac_rate_local(056) 
  reac_source_local(42,056) = + reac_rate_local(056) 
  reac_source_local(51,056) = - reac_rate_local(056) 
  reac_source_local(53,056) = + reac_rate_local(056) 
  reac_source_local(29,057) = - reac_rate_local(057) 
  reac_source_local(42,057) = + reac_rate_local(057) 
  reac_source_local(48,057) = - reac_rate_local(057) 
  reac_source_local(53,057) = + reac_rate_local(057) 
  reac_source_local(29,058) = - reac_rate_local(058) 
  reac_source_local(40,058) = + reac_rate_local(058) * 2.d0
  reac_source_local(49,058) = - reac_rate_local(058) 
  reac_source_local(53,058) = + reac_rate_local(058) 
  reac_source_local(21,059) = + reac_rate_local(059) 
  reac_source_local(29,059) = - reac_rate_local(059) 
  reac_source_local(40,059) = + reac_rate_local(059) 
  reac_source_local(50,059) = - reac_rate_local(059) 
  reac_source_local(53,059) = + reac_rate_local(059) 
  reac_source_local(29,060) = - reac_rate_local(060) 
  reac_source_local(32,060) = + reac_rate_local(060) 
  reac_source_local(40,060) = + reac_rate_local(060) 
  reac_source_local(51,060) = - reac_rate_local(060) 
  reac_source_local(53,060) = + reac_rate_local(060) 
  reac_source_local(01,061) = + reac_rate_local(061) 
  reac_source_local(10,061) = - reac_rate_local(061) 
  reac_source_local(40,061) = + reac_rate_local(061) 
  reac_source_local(48,061) = - reac_rate_local(061) 
  reac_source_local(53,061) = + reac_rate_local(061) 
  reac_source_local(01,062) = + reac_rate_local(062) 
  reac_source_local(10,062) = - reac_rate_local(062) 
  reac_source_local(42,062) = + reac_rate_local(062) 
  reac_source_local(50,062) = - reac_rate_local(062) 
  reac_source_local(53,062) = + reac_rate_local(062) 
  reac_source_local(01,063) = + reac_rate_local(063) 
  reac_source_local(11,063) = - reac_rate_local(063) 
  reac_source_local(40,063) = + reac_rate_local(063) 
  reac_source_local(48,063) = - reac_rate_local(063) 
  reac_source_local(53,063) = + reac_rate_local(063) 
  reac_source_local(01,064) = + reac_rate_local(064) 
  reac_source_local(11,064) = - reac_rate_local(064) 
  reac_source_local(42,064) = + reac_rate_local(064) 
  reac_source_local(50,064) = - reac_rate_local(064) 
  reac_source_local(53,064) = + reac_rate_local(064) 
  reac_source_local(11,065) = + reac_rate_local(065) 
  reac_source_local(13,065) = - reac_rate_local(065) 
  reac_source_local(10,066) = - reac_rate_local(066) 
  reac_source_local(15,066) = + reac_rate_local(066) 
  reac_source_local(29,066) = - reac_rate_local(066) 
  reac_source_local(40,066) = + reac_rate_local(066) 
  reac_source_local(01,067) = + reac_rate_local(067) 
  reac_source_local(10,067) = - reac_rate_local(067) 
  reac_source_local(21,067) = - reac_rate_local(067) 
  reac_source_local(29,067) = + reac_rate_local(067) 
  reac_source_local(30,067) = + reac_rate_local(067) 
  reac_source_local(01,068) = + reac_rate_local(068) 
  reac_source_local(10,068) = - reac_rate_local(068) 
  reac_source_local(21,068) = - reac_rate_local(068) 
  reac_source_local(27,068) = + reac_rate_local(068) 
  reac_source_local(01,069) = + reac_rate_local(069) 
  reac_source_local(10,069) = - reac_rate_local(069) 
  reac_source_local(01,070) = + reac_rate_local(070) 
  reac_source_local(10,070) = - reac_rate_local(070) 
  reac_source_local(01,071) = + reac_rate_local(071) 
  reac_source_local(10,071) = - reac_rate_local(071) 
  reac_source_local(14,071) = + reac_rate_local(071) 
  reac_source_local(40,071) = + reac_rate_local(071) 
  reac_source_local(41,071) = - reac_rate_local(071) 
  reac_source_local(01,072) = + reac_rate_local(072) 
  reac_source_local(10,072) = - reac_rate_local(072) 
  reac_source_local(29,072) = + reac_rate_local(072) 
  reac_source_local(40,072) = + reac_rate_local(072) 
  reac_source_local(42,072) = - reac_rate_local(072) 
  reac_source_local(01,073) = + reac_rate_local(073) 
  reac_source_local(10,073) = - reac_rate_local(073) * 2.d0
  reac_source_local(11,073) = + reac_rate_local(073) 
  reac_source_local(01,074) = + reac_rate_local(074) 
  reac_source_local(10,074) = - reac_rate_local(074) * 2.d0
  reac_source_local(13,074) = + reac_rate_local(074) 
  reac_source_local(10,075) = + reac_rate_local(075) 
  reac_source_local(11,075) = - reac_rate_local(075) 
  reac_source_local(01,076) = + reac_rate_local(076) 
  reac_source_local(11,076) = - reac_rate_local(076) 
  reac_source_local(01,077) = + reac_rate_local(077) 
  reac_source_local(11,077) = - reac_rate_local(077) 
  reac_source_local(21,077) = - reac_rate_local(077) 
  reac_source_local(29,077) = + reac_rate_local(077) * 2.d0
  reac_source_local(10,078) = + reac_rate_local(078) 
  reac_source_local(11,078) = - reac_rate_local(078) 
  reac_source_local(12,079) = + reac_rate_local(079) 
  reac_source_local(13,079) = - reac_rate_local(079) 
  reac_source_local(11,080) = + reac_rate_local(080) 
  reac_source_local(12,080) = - reac_rate_local(080) 
  reac_source_local(01,081) = + reac_rate_local(081) 
  reac_source_local(12,081) = - reac_rate_local(081) 
  reac_source_local(21,081) = - reac_rate_local(081) 
  reac_source_local(29,081) = + reac_rate_local(081) * 2.d0
  reac_source_local(01,082) = + reac_rate_local(082) 
  reac_source_local(12,082) = - reac_rate_local(082) 
  reac_source_local(14,082) = + reac_rate_local(082) 
  reac_source_local(29,082) = + reac_rate_local(082) 
  reac_source_local(40,082) = - reac_rate_local(082) 
  reac_source_local(15,083) = - reac_rate_local(083) 
  reac_source_local(21,083) = - reac_rate_local(083) 
  reac_source_local(29,083) = + reac_rate_local(083) 
  reac_source_local(40,083) = + reac_rate_local(083) 
  reac_source_local(01,084) = + reac_rate_local(084) 
  reac_source_local(15,084) = - reac_rate_local(084) 
  reac_source_local(29,084) = + reac_rate_local(084) 
  reac_source_local(40,084) = - reac_rate_local(084) 
  reac_source_local(01,085) = + reac_rate_local(085) 
  reac_source_local(15,085) = - reac_rate_local(085) 
  reac_source_local(40,085) = + reac_rate_local(085) 
  reac_source_local(41,085) = - reac_rate_local(085) 
  reac_source_local(16,086) = - reac_rate_local(086) 
  reac_source_local(21,086) = - reac_rate_local(086) 
  reac_source_local(29,086) = + reac_rate_local(086) 
  reac_source_local(40,086) = + reac_rate_local(086) 
  reac_source_local(10,087) = + reac_rate_local(087) 
  reac_source_local(16,087) = - reac_rate_local(087) 
  reac_source_local(29,087) = + reac_rate_local(087) 
  reac_source_local(40,087) = - reac_rate_local(087) 
  reac_source_local(14,088) = - reac_rate_local(088) 
  reac_source_local(26,088) = - reac_rate_local(088) 
  reac_source_local(29,088) = + reac_rate_local(088) 
  reac_source_local(40,088) = + reac_rate_local(088) 
  reac_source_local(29,089) = + reac_rate_local(089) 
  reac_source_local(30,089) = - reac_rate_local(089) 
  reac_source_local(14,090) = + reac_rate_local(090) 
  reac_source_local(21,090) = + reac_rate_local(090) 
  reac_source_local(30,090) = - reac_rate_local(090) 
  reac_source_local(40,090) = - reac_rate_local(090) 
  reac_source_local(30,091) = - reac_rate_local(091) 
  reac_source_local(40,091) = + reac_rate_local(091) * 2.d0
  reac_source_local(41,091) = - reac_rate_local(091) 
  reac_source_local(01,092) = + reac_rate_local(092) 
  reac_source_local(14,092) = - reac_rate_local(092) 
  reac_source_local(29,092) = + reac_rate_local(092) 
  reac_source_local(40,092) = - reac_rate_local(092) 
  reac_source_local(14,093) = - reac_rate_local(093) 
  reac_source_local(21,093) = - reac_rate_local(093) 
  reac_source_local(29,093) = + reac_rate_local(093) 
  reac_source_local(40,093) = + reac_rate_local(093) 
  reac_source_local(01,094) = + reac_rate_local(094) 
  reac_source_local(14,094) = - reac_rate_local(094) 
  reac_source_local(29,094) = + reac_rate_local(094) * 2.d0
  reac_source_local(42,094) = - reac_rate_local(094) 
  reac_source_local(14,095) = - reac_rate_local(095) 
  reac_source_local(29,095) = + reac_rate_local(095) 
  reac_source_local(41,095) = + reac_rate_local(095) 
  reac_source_local(42,095) = - reac_rate_local(095) 
  reac_source_local(01,096) = + reac_rate_local(096) 
  reac_source_local(14,096) = - reac_rate_local(096) 
  reac_source_local(21,096) = + reac_rate_local(096) 
  reac_source_local(42,096) = - reac_rate_local(096) 
  reac_source_local(14,097) = - reac_rate_local(097) 
  reac_source_local(40,097) = + reac_rate_local(097) * 2.d0
  reac_source_local(42,097) = - reac_rate_local(097) 
  reac_source_local(14,098) = + reac_rate_local(098) 
  reac_source_local(21,098) = + reac_rate_local(098) 
  reac_source_local(29,098) = - reac_rate_local(098) 
  reac_source_local(40,098) = - reac_rate_local(098) 
  reac_source_local(29,099) = - reac_rate_local(099) 
  reac_source_local(40,099) = - reac_rate_local(099) 
  reac_source_local(42,099) = + reac_rate_local(099) 
  reac_source_local(29,100) = - reac_rate_local(100) 
  reac_source_local(40,100) = + reac_rate_local(100) * 2.d0
  reac_source_local(41,100) = - reac_rate_local(100) 
  reac_source_local(21,101) = + reac_rate_local(101) 
  reac_source_local(29,101) = - reac_rate_local(101) 
  reac_source_local(40,101) = + reac_rate_local(101) 
  reac_source_local(42,101) = - reac_rate_local(101) 
  reac_source_local(21,102) = + reac_rate_local(102) 
  reac_source_local(29,102) = - reac_rate_local(102) 
  reac_source_local(42,102) = + reac_rate_local(102) 
  reac_source_local(43,102) = - reac_rate_local(102) 
  reac_source_local(21,103) = - reac_rate_local(103) 
  reac_source_local(29,103) = + reac_rate_local(103) 
  reac_source_local(40,103) = - reac_rate_local(103) 
  reac_source_local(42,103) = + reac_rate_local(103) 
  reac_source_local(21,104) = + reac_rate_local(104) 
  reac_source_local(32,104) = - reac_rate_local(104) 
  reac_source_local(40,104) = - reac_rate_local(104) 
  reac_source_local(42,104) = + reac_rate_local(104) 
  reac_source_local(01,105) = + reac_rate_local(105) 
  reac_source_local(40,105) = - reac_rate_local(105) 
  reac_source_local(41,105) = - reac_rate_local(105) 
  reac_source_local(42,105) = + reac_rate_local(105) 
  reac_source_local(40,106) = - reac_rate_local(106) 
  reac_source_local(42,106) = + reac_rate_local(106) * 2.d0
  reac_source_local(43,106) = - reac_rate_local(106) 
  reac_source_local(21,107) = - reac_rate_local(107) 
  reac_source_local(32,107) = + reac_rate_local(107) 
  reac_source_local(40,107) = + reac_rate_local(107) 
  reac_source_local(42,107) = - reac_rate_local(107) 
  reac_source_local(21,108) = + reac_rate_local(108) 
  reac_source_local(40,108) = + reac_rate_local(108) * 2.d0
  reac_source_local(42,108) = - reac_rate_local(108) * 2.d0
  reac_source_local(40,109) = + reac_rate_local(109) 
  reac_source_local(42,109) = - reac_rate_local(109) * 2.d0
  reac_source_local(43,109) = + reac_rate_local(109) 
  reac_source_local(21,110) = + reac_rate_local(110) 
  reac_source_local(32,110) = - reac_rate_local(110) 
  reac_source_local(42,110) = - reac_rate_local(110) 
  reac_source_local(43,110) = + reac_rate_local(110) 
  reac_source_local(21,111) = + reac_rate_local(111) 
  reac_source_local(40,111) = + reac_rate_local(111) 
  reac_source_local(43,111) = - reac_rate_local(111) 
  reac_source_local(21,112) = - reac_rate_local(112) 
  reac_source_local(32,112) = + reac_rate_local(112) 
  reac_source_local(42,112) = + reac_rate_local(112) 
  reac_source_local(43,112) = - reac_rate_local(112) 
  reac_source_local(21,113) = + reac_rate_local(113) 
  reac_source_local(42,113) = + reac_rate_local(113) * 2.d0
  reac_source_local(43,113) = - reac_rate_local(113) * 2.d0
  reac_source_local(29,114) = + reac_rate_local(114) 
  reac_source_local(40,114) = + reac_rate_local(114) 
  reac_source_local(42,114) = - reac_rate_local(114) 
  reac_source_local(29,115) = + reac_rate_local(115) 
  reac_source_local(42,115) = + reac_rate_local(115) 
  reac_source_local(43,115) = - reac_rate_local(115) 
  reac_source_local(29,116) = + reac_rate_local(116) 
  reac_source_local(42,116) = + reac_rate_local(116) 
  reac_source_local(43,116) = - reac_rate_local(116) 
  reac_source_local(29,117) = + reac_rate_local(117) 
  reac_source_local(42,117) = + reac_rate_local(117) 
  reac_source_local(43,117) = - reac_rate_local(117) 
  reac_source_local(29,118) = + reac_rate_local(118) 
  reac_source_local(42,118) = + reac_rate_local(118) 
  reac_source_local(43,118) = - reac_rate_local(118) 
  reac_source_local(29,119) = + reac_rate_local(119) 
  reac_source_local(42,119) = + reac_rate_local(119) 
  reac_source_local(43,119) = - reac_rate_local(119) 
  reac_source_local(21,120) = + reac_rate_local(120) 
  reac_source_local(40,120) = + reac_rate_local(120) 
  reac_source_local(43,120) = - reac_rate_local(120) 
  reac_source_local(21,121) = + reac_rate_local(121) 
  reac_source_local(40,121) = + reac_rate_local(121) 
  reac_source_local(43,121) = - reac_rate_local(121) 
  reac_source_local(21,122) = + reac_rate_local(122) 
  reac_source_local(40,122) = + reac_rate_local(122) 
  reac_source_local(43,122) = - reac_rate_local(122) 
  reac_source_local(21,123) = + reac_rate_local(123) 
  reac_source_local(40,123) = + reac_rate_local(123) 
  reac_source_local(43,123) = - reac_rate_local(123) 
  reac_source_local(21,124) = + reac_rate_local(124) 
  reac_source_local(40,124) = + reac_rate_local(124) 
  reac_source_local(43,124) = - reac_rate_local(124) 
  reac_source_local(42,125) = + reac_rate_local(125) 
  reac_source_local(43,125) = + reac_rate_local(125) 
  reac_source_local(44,125) = - reac_rate_local(125) 
  reac_source_local(14,126) = - reac_rate_local(126) 
  reac_source_local(29,126) = - reac_rate_local(126) 
  reac_source_local(40,126) = + reac_rate_local(126) 
  reac_source_local(14,127) = - reac_rate_local(127) 
  reac_source_local(29,127) = - reac_rate_local(127) 
  reac_source_local(40,127) = + reac_rate_local(127) 
  reac_source_local(14,128) = - reac_rate_local(128) 
  reac_source_local(29,128) = - reac_rate_local(128) 
  reac_source_local(40,128) = + reac_rate_local(128) 
  reac_source_local(14,129) = - reac_rate_local(129) 
  reac_source_local(29,129) = - reac_rate_local(129) 
  reac_source_local(40,129) = + reac_rate_local(129) 
  reac_source_local(14,130) = - reac_rate_local(130) 
  reac_source_local(29,130) = - reac_rate_local(130) 
  reac_source_local(40,130) = + reac_rate_local(130) 
  reac_source_local(21,131) = - reac_rate_local(131) 
  reac_source_local(29,131) = - reac_rate_local(131) 
  reac_source_local(32,131) = + reac_rate_local(131) 
  reac_source_local(29,132) = - reac_rate_local(132) 
  reac_source_local(40,132) = - reac_rate_local(132) 
  reac_source_local(42,132) = + reac_rate_local(132) 
  reac_source_local(29,133) = - reac_rate_local(133) 
  reac_source_local(40,133) = - reac_rate_local(133) 
  reac_source_local(42,133) = + reac_rate_local(133) 
  reac_source_local(29,134) = - reac_rate_local(134) 
  reac_source_local(40,134) = - reac_rate_local(134) 
  reac_source_local(42,134) = + reac_rate_local(134) 
  reac_source_local(29,135) = - reac_rate_local(135) 
  reac_source_local(42,135) = - reac_rate_local(135) 
  reac_source_local(43,135) = + reac_rate_local(135) 
  reac_source_local(29,136) = - reac_rate_local(136) 
  reac_source_local(42,136) = - reac_rate_local(136) 
  reac_source_local(43,136) = + reac_rate_local(136) 
  reac_source_local(29,137) = - reac_rate_local(137) 
  reac_source_local(42,137) = - reac_rate_local(137) 
  reac_source_local(43,137) = + reac_rate_local(137) 
  reac_source_local(29,138) = - reac_rate_local(138) 
  reac_source_local(42,138) = - reac_rate_local(138) 
  reac_source_local(43,138) = + reac_rate_local(138) 
  reac_source_local(29,139) = - reac_rate_local(139) 
  reac_source_local(42,139) = - reac_rate_local(139) 
  reac_source_local(43,139) = + reac_rate_local(139) 
  reac_source_local(42,140) = - reac_rate_local(140) 
  reac_source_local(43,140) = - reac_rate_local(140) 
  reac_source_local(44,140) = + reac_rate_local(140) 
  reac_source_local(17,141) = - reac_rate_local(141) 
  reac_source_local(21,141) = - reac_rate_local(141) 
  reac_source_local(33,141) = + reac_rate_local(141) 
  reac_source_local(40,141) = + reac_rate_local(141) 
  reac_source_local(14,142) = + reac_rate_local(142) 
  reac_source_local(17,142) = - reac_rate_local(142) 
  reac_source_local(40,142) = - reac_rate_local(142) 
  reac_source_local(45,142) = + reac_rate_local(142) 
  reac_source_local(17,143) = - reac_rate_local(143) 
  reac_source_local(18,143) = + reac_rate_local(143) 
  reac_source_local(29,143) = + reac_rate_local(143) 
  reac_source_local(40,143) = - reac_rate_local(143) 
  reac_source_local(01,144) = + reac_rate_local(144) 
  reac_source_local(17,144) = - reac_rate_local(144) 
  reac_source_local(33,144) = + reac_rate_local(144) 
  reac_source_local(40,144) = - reac_rate_local(144) 
  reac_source_local(29,145) = + reac_rate_local(145) 
  reac_source_local(33,145) = - reac_rate_local(145) 
  reac_source_local(40,145) = - reac_rate_local(145) 
  reac_source_local(45,145) = + reac_rate_local(145) 
  reac_source_local(14,146) = + reac_rate_local(146) 
  reac_source_local(33,146) = - reac_rate_local(146) 
  reac_source_local(34,146) = + reac_rate_local(146) 
  reac_source_local(40,146) = - reac_rate_local(146) 
  reac_source_local(33,147) = - reac_rate_local(147) 
  reac_source_local(40,147) = + reac_rate_local(147) 
  reac_source_local(41,147) = - reac_rate_local(147) 
  reac_source_local(45,147) = + reac_rate_local(147) 
  reac_source_local(29,148) = + reac_rate_local(148) 
  reac_source_local(33,148) = - reac_rate_local(148) 
  reac_source_local(42,148) = - reac_rate_local(148) 
  reac_source_local(47,148) = + reac_rate_local(148) 
  reac_source_local(01,149) = + reac_rate_local(149) 
  reac_source_local(18,149) = - reac_rate_local(149) 
  reac_source_local(40,149) = - reac_rate_local(149) 
  reac_source_local(45,149) = + reac_rate_local(149) 
  reac_source_local(01,150) = - reac_rate_local(150) 
  reac_source_local(34,150) = - reac_rate_local(150) 
  reac_source_local(40,150) = + reac_rate_local(150) 
  reac_source_local(45,150) = + reac_rate_local(150) 
  reac_source_local(21,151) = + reac_rate_local(151) 
  reac_source_local(34,151) = - reac_rate_local(151) 
  reac_source_local(40,151) = - reac_rate_local(151) 
  reac_source_local(45,151) = + reac_rate_local(151) 
  reac_source_local(32,152) = + reac_rate_local(152) 
  reac_source_local(34,152) = - reac_rate_local(152) 
  reac_source_local(42,152) = - reac_rate_local(152) 
  reac_source_local(45,152) = + reac_rate_local(152) 
  reac_source_local(21,153) = + reac_rate_local(153) 
  reac_source_local(34,153) = - reac_rate_local(153) 
  reac_source_local(42,153) = - reac_rate_local(153) 
  reac_source_local(47,153) = + reac_rate_local(153) 
  reac_source_local(01,154) = + reac_rate_local(154) 
  reac_source_local(14,154) = + reac_rate_local(154) 
  reac_source_local(19,154) = - reac_rate_local(154) 
  reac_source_local(40,154) = - reac_rate_local(154) 
  reac_source_local(45,154) = + reac_rate_local(154) 
  reac_source_local(01,155) = + reac_rate_local(155) 
  reac_source_local(19,155) = - reac_rate_local(155) 
  reac_source_local(40,155) = - reac_rate_local(155) 
  reac_source_local(46,155) = + reac_rate_local(155) 
  reac_source_local(40,156) = - reac_rate_local(156) 
  reac_source_local(42,156) = + reac_rate_local(156) 
  reac_source_local(45,156) = + reac_rate_local(156) 
  reac_source_local(47,156) = - reac_rate_local(156) 
  reac_source_local(40,157) = - reac_rate_local(157) 
  reac_source_local(41,157) = + reac_rate_local(157) 
  reac_source_local(45,157) = + reac_rate_local(157) 
  reac_source_local(46,157) = - reac_rate_local(157) 
  reac_source_local(01,158) = + reac_rate_local(158) * 2.d0
  reac_source_local(20,158) = - reac_rate_local(158) 
  reac_source_local(40,158) = - reac_rate_local(158) 
  reac_source_local(45,158) = + reac_rate_local(158) 
  reac_source_local(21,159) = + reac_rate_local(159) * 2.d0
  reac_source_local(35,159) = - reac_rate_local(159) 
  reac_source_local(40,159) = - reac_rate_local(159) 
  reac_source_local(45,159) = + reac_rate_local(159) 
  reac_source_local(29,160) = + reac_rate_local(160) 
  reac_source_local(36,160) = - reac_rate_local(160) 
  reac_source_local(42,160) = - reac_rate_local(160) 
  reac_source_local(50,160) = + reac_rate_local(160) 
  reac_source_local(36,161) = - reac_rate_local(161) 
  reac_source_local(40,161) = + reac_rate_local(161) 
  reac_source_local(41,161) = - reac_rate_local(161) 
  reac_source_local(48,161) = + reac_rate_local(161) 
  reac_source_local(21,162) = + reac_rate_local(162) 
  reac_source_local(37,162) = - reac_rate_local(162) 
  reac_source_local(42,162) = - reac_rate_local(162) 
  reac_source_local(50,162) = + reac_rate_local(162) 
  reac_source_local(29,163) = + reac_rate_local(163) 
  reac_source_local(38,163) = - reac_rate_local(163) 
  reac_source_local(40,163) = - reac_rate_local(163) 
  reac_source_local(51,163) = + reac_rate_local(163) 
  reac_source_local(21,164) = + reac_rate_local(164) 
  reac_source_local(38,164) = - reac_rate_local(164) 
  reac_source_local(40,164) = - reac_rate_local(164) 
  reac_source_local(50,164) = + reac_rate_local(164) 
  reac_source_local(32,165) = + reac_rate_local(165) 
  reac_source_local(38,165) = - reac_rate_local(165) 
  reac_source_local(42,165) = - reac_rate_local(165) 
  reac_source_local(50,165) = + reac_rate_local(165) 
  reac_source_local(21,166) = + reac_rate_local(166) 
  reac_source_local(38,166) = - reac_rate_local(166) 
  reac_source_local(42,166) = - reac_rate_local(166) 
  reac_source_local(51,166) = + reac_rate_local(166) 
  reac_source_local(21,167) = - reac_rate_local(167) 
  reac_source_local(37,167) = + reac_rate_local(167) 
  reac_source_local(40,167) = + reac_rate_local(167) 
  reac_source_local(48,167) = - reac_rate_local(167) 
  reac_source_local(40,168) = + reac_rate_local(168) 
  reac_source_local(42,168) = - reac_rate_local(168) 
  reac_source_local(48,168) = - reac_rate_local(168) 
  reac_source_local(50,168) = + reac_rate_local(168) 
  reac_source_local(40,169) = + reac_rate_local(169) 
  reac_source_local(42,169) = - reac_rate_local(169) 
  reac_source_local(50,169) = - reac_rate_local(169) 
  reac_source_local(51,169) = + reac_rate_local(169) 
  reac_source_local(42,170) = + reac_rate_local(170) 
  reac_source_local(43,170) = - reac_rate_local(170) 
  reac_source_local(50,170) = - reac_rate_local(170) 
  reac_source_local(51,170) = + reac_rate_local(170) 
  reac_source_local(42,171) = + reac_rate_local(171) * 2.d0
  reac_source_local(44,171) = - reac_rate_local(171) 
  reac_source_local(50,171) = - reac_rate_local(171) 
  reac_source_local(51,171) = + reac_rate_local(171) 
  reac_source_local(40,172) = - reac_rate_local(172) 
  reac_source_local(42,172) = + reac_rate_local(172) 
  reac_source_local(50,172) = + reac_rate_local(172) 
  reac_source_local(51,172) = - reac_rate_local(172) 
  reac_source_local(21,173) = + reac_rate_local(173) 
  reac_source_local(39,173) = - reac_rate_local(173) 
  reac_source_local(40,173) = - reac_rate_local(173) 
  reac_source_local(51,173) = + reac_rate_local(173) 
  reac_source_local(36,174) = - reac_rate_local(174) 
  reac_source_local(40,174) = - reac_rate_local(174) 
  reac_source_local(50,174) = + reac_rate_local(174) 
  reac_source_local(29,175) = + reac_rate_local(175) 
  reac_source_local(36,175) = - reac_rate_local(175) 
  reac_source_local(40,175) = + reac_rate_local(175) 
  reac_source_local(45,175) = - reac_rate_local(175) 
  reac_source_local(29,176) = + reac_rate_local(176) 
  reac_source_local(36,176) = - reac_rate_local(176) 
  reac_source_local(42,176) = + reac_rate_local(176) 
  reac_source_local(47,176) = - reac_rate_local(176) 
  reac_source_local(21,177) = + reac_rate_local(177) 
  reac_source_local(37,177) = - reac_rate_local(177) 
  reac_source_local(40,177) = + reac_rate_local(177) 
  reac_source_local(45,177) = - reac_rate_local(177) 
  reac_source_local(21,178) = + reac_rate_local(178) 
  reac_source_local(37,178) = - reac_rate_local(178) 
  reac_source_local(42,178) = + reac_rate_local(178) 
  reac_source_local(47,178) = - reac_rate_local(178) 
  reac_source_local(32,179) = + reac_rate_local(179) 
  reac_source_local(38,179) = - reac_rate_local(179) 
  reac_source_local(40,179) = + reac_rate_local(179) 
  reac_source_local(45,179) = - reac_rate_local(179) 
  reac_source_local(32,180) = + reac_rate_local(180) 
  reac_source_local(38,180) = - reac_rate_local(180) 
  reac_source_local(42,180) = + reac_rate_local(180) 
  reac_source_local(47,180) = - reac_rate_local(180) 
  reac_source_local(14,181) = + reac_rate_local(181) 
  reac_source_local(17,181) = - reac_rate_local(181) 
  reac_source_local(40,181) = + reac_rate_local(181) 
  reac_source_local(48,181) = - reac_rate_local(181) 
  reac_source_local(01,182) = + reac_rate_local(182) 
  reac_source_local(18,182) = - reac_rate_local(182) 
  reac_source_local(40,182) = + reac_rate_local(182) 
  reac_source_local(48,182) = - reac_rate_local(182) 
  reac_source_local(29,183) = + reac_rate_local(183) 
  reac_source_local(33,183) = - reac_rate_local(183) 
  reac_source_local(40,183) = + reac_rate_local(183) 
  reac_source_local(48,183) = - reac_rate_local(183) 
  reac_source_local(21,184) = + reac_rate_local(184) 
  reac_source_local(34,184) = - reac_rate_local(184) 
  reac_source_local(40,184) = + reac_rate_local(184) 
  reac_source_local(48,184) = - reac_rate_local(184) 
  reac_source_local(40,185) = + reac_rate_local(185) * 2.d0
  reac_source_local(45,185) = - reac_rate_local(185) 
  reac_source_local(48,185) = - reac_rate_local(185) 
  reac_source_local(40,186) = + reac_rate_local(186) 
  reac_source_local(41,186) = + reac_rate_local(186) 
  reac_source_local(46,186) = - reac_rate_local(186) 
  reac_source_local(48,186) = - reac_rate_local(186) 
  reac_source_local(40,187) = + reac_rate_local(187) 
  reac_source_local(42,187) = + reac_rate_local(187) 
  reac_source_local(47,187) = - reac_rate_local(187) 
  reac_source_local(48,187) = - reac_rate_local(187) 
  reac_source_local(40,188) = + reac_rate_local(188) 
  reac_source_local(41,188) = + reac_rate_local(188) 
  reac_source_local(45,188) = - reac_rate_local(188) 
  reac_source_local(49,188) = - reac_rate_local(188) 
  reac_source_local(41,189) = + reac_rate_local(189) 
  reac_source_local(42,189) = + reac_rate_local(189) 
  reac_source_local(47,189) = - reac_rate_local(189) 
  reac_source_local(49,189) = - reac_rate_local(189) 
  reac_source_local(14,190) = + reac_rate_local(190) 
  reac_source_local(17,190) = - reac_rate_local(190) 
  reac_source_local(42,190) = + reac_rate_local(190) 
  reac_source_local(50,190) = - reac_rate_local(190) 
  reac_source_local(01,191) = + reac_rate_local(191) 
  reac_source_local(18,191) = - reac_rate_local(191) 
  reac_source_local(42,191) = + reac_rate_local(191) 
  reac_source_local(50,191) = - reac_rate_local(191) 
  reac_source_local(29,192) = + reac_rate_local(192) 
  reac_source_local(33,192) = - reac_rate_local(192) 
  reac_source_local(42,192) = + reac_rate_local(192) 
  reac_source_local(50,192) = - reac_rate_local(192) 
  reac_source_local(21,193) = + reac_rate_local(193) 
  reac_source_local(34,193) = - reac_rate_local(193) 
  reac_source_local(42,193) = + reac_rate_local(193) 
  reac_source_local(50,193) = - reac_rate_local(193) 
  reac_source_local(40,194) = + reac_rate_local(194) 
  reac_source_local(42,194) = + reac_rate_local(194) 
  reac_source_local(45,194) = - reac_rate_local(194) 
  reac_source_local(50,194) = - reac_rate_local(194) 
  reac_source_local(41,195) = + reac_rate_local(195) 
  reac_source_local(42,195) = + reac_rate_local(195) 
  reac_source_local(46,195) = - reac_rate_local(195) 
  reac_source_local(50,195) = - reac_rate_local(195) 
  reac_source_local(42,196) = + reac_rate_local(196) * 2.d0
  reac_source_local(47,196) = - reac_rate_local(196) 
  reac_source_local(50,196) = - reac_rate_local(196) 
  reac_source_local(01,197) = + reac_rate_local(197) 
  reac_source_local(18,197) = - reac_rate_local(197) 
  reac_source_local(43,197) = + reac_rate_local(197) 
  reac_source_local(51,197) = - reac_rate_local(197) 
  reac_source_local(40,198) = + reac_rate_local(198) 
  reac_source_local(43,198) = + reac_rate_local(198) 
  reac_source_local(45,198) = - reac_rate_local(198) 
  reac_source_local(51,198) = - reac_rate_local(198) 
  reac_source_local(42,199) = + reac_rate_local(199) 
  reac_source_local(43,199) = + reac_rate_local(199) 
  reac_source_local(47,199) = - reac_rate_local(199) 
  reac_source_local(51,199) = - reac_rate_local(199) 
  reac_source_local(01,200) = + reac_rate_local(200) 
  reac_source_local(21,200) = + reac_rate_local(200) * 2.d0
  reac_source_local(37,200) = - reac_rate_local(200) 
  reac_source_local(52,200) = - reac_rate_local(200) 
  reac_source_local(14,201) = + reac_rate_local(201) * 2.d0
  reac_source_local(18,201) = - reac_rate_local(201) 
  reac_source_local(40,201) = + reac_rate_local(201) 
  reac_source_local(48,201) = - reac_rate_local(201) 
  reac_source_local(01,202) = + reac_rate_local(202) 
  reac_source_local(14,202) = + reac_rate_local(202) 
  reac_source_local(19,202) = - reac_rate_local(202) 
  reac_source_local(40,202) = + reac_rate_local(202) 
  reac_source_local(48,202) = - reac_rate_local(202) 
  reac_source_local(01,203) = + reac_rate_local(203) * 2.d0
  reac_source_local(20,203) = - reac_rate_local(203) 
  reac_source_local(40,203) = + reac_rate_local(203) 
  reac_source_local(48,203) = - reac_rate_local(203) 
  reac_source_local(29,204) = + reac_rate_local(204) * 2.d0
  reac_source_local(34,204) = - reac_rate_local(204) 
  reac_source_local(40,204) = + reac_rate_local(204) 
  reac_source_local(48,204) = - reac_rate_local(204) 
  reac_source_local(21,205) = + reac_rate_local(205) * 2.d0
  reac_source_local(35,205) = - reac_rate_local(205) 
  reac_source_local(40,205) = + reac_rate_local(205) 
  reac_source_local(48,205) = - reac_rate_local(205) 
  reac_source_local(14,206) = + reac_rate_local(206) 
  reac_source_local(29,206) = + reac_rate_local(206) 
  reac_source_local(40,206) = + reac_rate_local(206) 
  reac_source_local(45,206) = - reac_rate_local(206) 
  reac_source_local(48,206) = - reac_rate_local(206) 
  reac_source_local(01,207) = + reac_rate_local(207) 
  reac_source_local(29,207) = + reac_rate_local(207) 
  reac_source_local(40,207) = + reac_rate_local(207) 
  reac_source_local(46,207) = - reac_rate_local(207) 
  reac_source_local(48,207) = - reac_rate_local(207) 
  reac_source_local(14,208) = + reac_rate_local(208) 
  reac_source_local(21,208) = + reac_rate_local(208) 
  reac_source_local(40,208) = + reac_rate_local(208) 
  reac_source_local(47,208) = - reac_rate_local(208) 
  reac_source_local(48,208) = - reac_rate_local(208) 
  reac_source_local(01,209) = + reac_rate_local(209) 
  reac_source_local(21,209) = + reac_rate_local(209) 
  reac_source_local(40,209) = + reac_rate_local(209) 
  reac_source_local(48,209) = - reac_rate_local(209) 
  reac_source_local(52,209) = - reac_rate_local(209) 
  reac_source_local(14,210) = + reac_rate_local(210) * 2.d0
  reac_source_local(18,210) = - reac_rate_local(210) 
  reac_source_local(41,210) = + reac_rate_local(210) 
  reac_source_local(49,210) = - reac_rate_local(210) 
  reac_source_local(14,211) = + reac_rate_local(211) * 2.d0
  reac_source_local(18,211) = - reac_rate_local(211) 
  reac_source_local(42,211) = + reac_rate_local(211) 
  reac_source_local(50,211) = - reac_rate_local(211) 
  reac_source_local(01,212) = + reac_rate_local(212) 
  reac_source_local(14,212) = + reac_rate_local(212) 
  reac_source_local(19,212) = - reac_rate_local(212) 
  reac_source_local(42,212) = + reac_rate_local(212) 
  reac_source_local(50,212) = - reac_rate_local(212) 
  reac_source_local(01,213) = + reac_rate_local(213) * 2.d0
  reac_source_local(20,213) = - reac_rate_local(213) 
  reac_source_local(42,213) = + reac_rate_local(213) 
  reac_source_local(50,213) = - reac_rate_local(213) 
  reac_source_local(29,214) = + reac_rate_local(214) * 2.d0
  reac_source_local(34,214) = - reac_rate_local(214) 
  reac_source_local(42,214) = + reac_rate_local(214) 
  reac_source_local(50,214) = - reac_rate_local(214) 
  reac_source_local(21,215) = + reac_rate_local(215) * 2.d0
  reac_source_local(35,215) = - reac_rate_local(215) 
  reac_source_local(42,215) = + reac_rate_local(215) 
  reac_source_local(50,215) = - reac_rate_local(215) 
  reac_source_local(14,216) = + reac_rate_local(216) 
  reac_source_local(29,216) = + reac_rate_local(216) 
  reac_source_local(42,216) = + reac_rate_local(216) 
  reac_source_local(45,216) = - reac_rate_local(216) 
  reac_source_local(50,216) = - reac_rate_local(216) 
  reac_source_local(01,217) = + reac_rate_local(217) 
  reac_source_local(29,217) = + reac_rate_local(217) 
  reac_source_local(42,217) = + reac_rate_local(217) 
  reac_source_local(46,217) = - reac_rate_local(217) 
  reac_source_local(50,217) = - reac_rate_local(217) 
  reac_source_local(14,218) = + reac_rate_local(218) 
  reac_source_local(21,218) = + reac_rate_local(218) 
  reac_source_local(42,218) = + reac_rate_local(218) 
  reac_source_local(47,218) = - reac_rate_local(218) 
  reac_source_local(50,218) = - reac_rate_local(218) 
  reac_source_local(01,219) = + reac_rate_local(219) 
  reac_source_local(21,219) = + reac_rate_local(219) 
  reac_source_local(42,219) = + reac_rate_local(219) 
  reac_source_local(50,219) = - reac_rate_local(219) 
  reac_source_local(52,219) = - reac_rate_local(219) 
  reac_source_local(21,220) = + reac_rate_local(220) * 2.d0
  reac_source_local(39,220) = - reac_rate_local(220) 
  reac_source_local(40,220) = + reac_rate_local(220) 
  reac_source_local(45,220) = - reac_rate_local(220) 
  reac_source_local(21,221) = + reac_rate_local(221) * 2.d0
  reac_source_local(39,221) = - reac_rate_local(221) 
  reac_source_local(42,221) = + reac_rate_local(221) 
  reac_source_local(47,221) = - reac_rate_local(221) 
  reac_source_local(29,222) = + reac_rate_local(222) 
  reac_source_local(36,222) = - reac_rate_local(222) 
  reac_source_local(40,222) = + reac_rate_local(222) 
  reac_source_local(45,222) = - reac_rate_local(222) 
  reac_source_local(21,223) = + reac_rate_local(223) 
  reac_source_local(37,223) = - reac_rate_local(223) 
  reac_source_local(40,223) = + reac_rate_local(223) 
  reac_source_local(45,223) = - reac_rate_local(223) 
  reac_source_local(17,224) = - reac_rate_local(224) 
  reac_source_local(36,224) = - reac_rate_local(224) 
  reac_source_local(40,224) = + reac_rate_local(224) 
  reac_source_local(36,225) = - reac_rate_local(225) 
  reac_source_local(42,225) = + reac_rate_local(225) 
  reac_source_local(45,225) = - reac_rate_local(225) 
  reac_source_local(17,226) = - reac_rate_local(226) 
  reac_source_local(37,226) = - reac_rate_local(226) 
  reac_source_local(42,226) = + reac_rate_local(226) 
  reac_source_local(32,227) = + reac_rate_local(227) 
  reac_source_local(38,227) = - reac_rate_local(227) 
  reac_source_local(40,227) = + reac_rate_local(227) 
  reac_source_local(45,227) = - reac_rate_local(227) 
  reac_source_local(32,228) = + reac_rate_local(228) 
  reac_source_local(38,228) = - reac_rate_local(228) 
  reac_source_local(42,228) = + reac_rate_local(228) 
  reac_source_local(47,228) = - reac_rate_local(228) 
  reac_source_local(14,229) = + reac_rate_local(229) 
  reac_source_local(17,229) = - reac_rate_local(229) 
  reac_source_local(40,229) = + reac_rate_local(229) 
  reac_source_local(48,229) = - reac_rate_local(229) 
  reac_source_local(01,230) = + reac_rate_local(230) 
  reac_source_local(18,230) = - reac_rate_local(230) 
  reac_source_local(40,230) = + reac_rate_local(230) 
  reac_source_local(48,230) = - reac_rate_local(230) 
  reac_source_local(29,231) = + reac_rate_local(231) 
  reac_source_local(33,231) = - reac_rate_local(231) 
  reac_source_local(40,231) = + reac_rate_local(231) 
  reac_source_local(48,231) = - reac_rate_local(231) 
  reac_source_local(21,232) = + reac_rate_local(232) 
  reac_source_local(34,232) = - reac_rate_local(232) 
  reac_source_local(40,232) = + reac_rate_local(232) 
  reac_source_local(48,232) = - reac_rate_local(232) 
  reac_source_local(40,233) = + reac_rate_local(233) * 2.d0
  reac_source_local(45,233) = - reac_rate_local(233) 
  reac_source_local(48,233) = - reac_rate_local(233) 
  reac_source_local(40,234) = + reac_rate_local(234) 
  reac_source_local(41,234) = + reac_rate_local(234) 
  reac_source_local(46,234) = - reac_rate_local(234) 
  reac_source_local(48,234) = - reac_rate_local(234) 
  reac_source_local(40,235) = + reac_rate_local(235) 
  reac_source_local(42,235) = + reac_rate_local(235) 
  reac_source_local(47,235) = - reac_rate_local(235) 
  reac_source_local(48,235) = - reac_rate_local(235) 
  reac_source_local(40,236) = + reac_rate_local(236) 
  reac_source_local(41,236) = + reac_rate_local(236) 
  reac_source_local(45,236) = - reac_rate_local(236) 
  reac_source_local(49,236) = - reac_rate_local(236) 
  reac_source_local(41,237) = + reac_rate_local(237) 
  reac_source_local(42,237) = + reac_rate_local(237) 
  reac_source_local(47,237) = - reac_rate_local(237) 
  reac_source_local(49,237) = - reac_rate_local(237) 
  reac_source_local(14,238) = + reac_rate_local(238) 
  reac_source_local(17,238) = - reac_rate_local(238) 
  reac_source_local(42,238) = + reac_rate_local(238) 
  reac_source_local(50,238) = - reac_rate_local(238) 
  reac_source_local(01,239) = + reac_rate_local(239) 
  reac_source_local(18,239) = - reac_rate_local(239) 
  reac_source_local(42,239) = + reac_rate_local(239) 
  reac_source_local(50,239) = - reac_rate_local(239) 
  reac_source_local(29,240) = + reac_rate_local(240) 
  reac_source_local(33,240) = - reac_rate_local(240) 
  reac_source_local(42,240) = + reac_rate_local(240) 
  reac_source_local(50,240) = - reac_rate_local(240) 
  reac_source_local(21,241) = + reac_rate_local(241) 
  reac_source_local(34,241) = - reac_rate_local(241) 
  reac_source_local(42,241) = + reac_rate_local(241) 
  reac_source_local(50,241) = - reac_rate_local(241) 
  reac_source_local(40,242) = + reac_rate_local(242) 
  reac_source_local(42,242) = + reac_rate_local(242) 
  reac_source_local(45,242) = - reac_rate_local(242) 
  reac_source_local(50,242) = - reac_rate_local(242) 
  reac_source_local(41,243) = + reac_rate_local(243) 
  reac_source_local(42,243) = + reac_rate_local(243) 
  reac_source_local(46,243) = - reac_rate_local(243) 
  reac_source_local(50,243) = - reac_rate_local(243) 
  reac_source_local(42,244) = + reac_rate_local(244) * 2.d0
  reac_source_local(47,244) = - reac_rate_local(244) 
  reac_source_local(50,244) = - reac_rate_local(244) 
  reac_source_local(40,245) = + reac_rate_local(245) 
  reac_source_local(43,245) = + reac_rate_local(245) 
  reac_source_local(45,245) = - reac_rate_local(245) 
  reac_source_local(51,245) = - reac_rate_local(245) 
  reac_source_local(42,246) = + reac_rate_local(246) 
  reac_source_local(43,246) = + reac_rate_local(246) 
  reac_source_local(47,246) = - reac_rate_local(246) 
  reac_source_local(51,246) = - reac_rate_local(246) 
  return
end subroutine ZDPlasKin_reac_source_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction source terms
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_fex(neq,t,y,ydot)
  implicit none
  integer,          intent(in)  :: neq
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: ydot(neq)
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(54)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(53) 
  rrt(002) = rrt(002) * density(01) * density(53) 
  rrt(003) = rrt(003) * density(01) * density(53) 
  rrt(004) = rrt(004) * density(01) * density(53) 
  rrt(005) = rrt(005) * density(01) * density(53) 
  rrt(006) = rrt(006) * density(01) * density(53) 
  rrt(007) = rrt(007) * density(01) * density(53) 
  rrt(008) = rrt(008) * density(01) * density(53) 
  rrt(009) = rrt(009) * density(01) * density(53) 
  rrt(010) = rrt(010) * density(02) * density(53) 
  rrt(011) = rrt(011) * density(03) * density(53) 
  rrt(012) = rrt(012) * density(04) * density(53) 
  rrt(013) = rrt(013) * density(05) * density(53) 
  rrt(014) = rrt(014) * density(06) * density(53) 
  rrt(015) = rrt(015) * density(07) * density(53) 
  rrt(016) = rrt(016) * density(08) * density(53) 
  rrt(017) = rrt(017) * density(09) * density(53) 
  rrt(018) = rrt(018) * density(01) * density(02) 
  rrt(019) = rrt(019) * density(01) * density(03) 
  rrt(020) = rrt(020) * density(01) * density(04) 
  rrt(021) = rrt(021) * density(01) * density(05) 
  rrt(022) = rrt(022) * density(01) * density(06) 
  rrt(023) = rrt(023) * density(01) * density(07) 
  rrt(024) = rrt(024) * density(01) * density(08) 
  rrt(025) = rrt(025) * density(01) * density(09) 
  rrt(026) = rrt(026) * density(02) * density(29) 
  rrt(027) = rrt(027) * density(03) * density(29) 
  rrt(028) = rrt(028) * density(04) * density(29) 
  rrt(029) = rrt(029) * density(05) * density(29) 
  rrt(030) = rrt(030) * density(06) * density(29) 
  rrt(031) = rrt(031) * density(07) * density(29) 
  rrt(032) = rrt(032) * density(08) * density(29) 
  rrt(033) = rrt(033) * density(09) * density(29) 
  rrt(034) = rrt(034) * density(01) * density(53) 
  rrt(035) = rrt(035) * density(01) * density(53) 
  rrt(036) = rrt(036) * density(01) * density(53) 
  rrt(037) = rrt(037) * density(01) * density(53) 
  rrt(038) = rrt(038) * density(01) * density(53) 
  rrt(039) = rrt(039) * density(01) * density(53) 
  rrt(040) = rrt(040) * density(01) * density(53) 
  rrt(041) = rrt(041) * density(01) * density(53) 
  rrt(042) = rrt(042) * density(01) * density(53) 
  rrt(043) = rrt(043) * density(01) * density(53) 
  rrt(044) = rrt(044) * density(01) * density(53) 
  rrt(045) = rrt(045) * density(01) * density(53) 
  rrt(046) = rrt(046) * density(40) * density(53) 
  rrt(047) = rrt(047) * density(47) * density(53) 
  rrt(048) = rrt(048) * density(42) * density(53) 
  rrt(049) = rrt(049) * density(40) * density(53) 
  rrt(050) = rrt(050) * density(14) * density(36) 
  rrt(051) = rrt(051) * density(36) * density(40) 
  rrt(052) = rrt(052) * density(14) * density(37) 
  rrt(053) = rrt(053) * density(14) * density(38) 
  rrt(054) = rrt(054) * density(14) * density(49) 
  rrt(055) = rrt(055) * density(14) * density(50) 
  rrt(056) = rrt(056) * density(14) * density(51) 
  rrt(057) = rrt(057) * density(29) * density(48) 
  rrt(058) = rrt(058) * density(29) * density(49) 
  rrt(059) = rrt(059) * density(29) * density(50) 
  rrt(060) = rrt(060) * density(29) * density(51) 
  rrt(061) = rrt(061) * density(10) * density(48) 
  rrt(062) = rrt(062) * density(10) * density(50) 
  rrt(063) = rrt(063) * density(11) * density(48) 
  rrt(064) = rrt(064) * density(11) * density(50) 
  rrt(065) = rrt(065) * density(13) 
  rrt(066) = rrt(066) * density(10) * density(29) 
  rrt(067) = rrt(067) * density(10) * density(21) 
  rrt(068) = rrt(068) * density(10) * density(21) 
  rrt(069) = rrt(069) * density(01) * density(10) 
  rrt(070) = rrt(070) * density(10) * density(40) 
  rrt(071) = rrt(071) * density(10) * density(41) 
  rrt(072) = rrt(072) * density(10) * density(42) 
  rrt(073) = rrt(073) * density(10)**2 
  rrt(074) = rrt(074) * density(10)**2 
  rrt(075) = rrt(075) * density(01) * density(11) 
  rrt(076) = rrt(076) * density(01) * density(11) 
  rrt(077) = rrt(077) * density(11) * density(21) 
  rrt(078) = rrt(078) * density(11) * density(40) 
  rrt(079) = rrt(079) * density(01) * density(13) 
  rrt(080) = rrt(080) * density(01) * density(12) 
  rrt(081) = rrt(081) * density(12) * density(21) 
  rrt(082) = rrt(082) * density(12) * density(40) 
  rrt(083) = rrt(083) * density(15) * density(21) 
  rrt(084) = rrt(084) * density(15) * density(40) 
  rrt(085) = rrt(085) * density(15) * density(41) 
  rrt(086) = rrt(086) * density(16) * density(21) 
  rrt(087) = rrt(087) * density(16) * density(40) 
  rrt(088) = rrt(088) * density(14) * density(26) 
  rrt(089) = rrt(089) * density(01) * density(30) 
  rrt(090) = rrt(090) * density(30) * density(40) 
  rrt(091) = rrt(091) * density(30) * density(41) 
  rrt(092) = rrt(092) * density(14) * density(40) 
  rrt(093) = rrt(093) * density(14) * density(21) 
  rrt(094) = rrt(094) * density(14) * density(42) 
  rrt(095) = rrt(095) * density(14) * density(42) 
  rrt(096) = rrt(096) * density(14) * density(42) 
  rrt(097) = rrt(097) * density(14) * density(42) 
  rrt(098) = rrt(098) * density(29) * density(40) 
  rrt(099) = rrt(099) * density(29) * density(40) 
  rrt(100) = rrt(100) * density(29) * density(41) 
  rrt(101) = rrt(101) * density(29) * density(42) 
  rrt(102) = rrt(102) * density(29) * density(43) 
  rrt(103) = rrt(103) * density(21) * density(40) 
  rrt(104) = rrt(104) * density(32) * density(40) 
  rrt(105) = rrt(105) * density(40) * density(41) 
  rrt(106) = rrt(106) * density(40) * density(43) 
  rrt(107) = rrt(107) * density(21) * density(42) 
  rrt(108) = rrt(108) * density(42)**2 
  rrt(109) = rrt(109) * density(42)**2 
  rrt(110) = rrt(110) * density(32) * density(42) 
  rrt(111) = rrt(111) * density(42) * density(43) 
  rrt(112) = rrt(112) * density(21) * density(43) 
  rrt(113) = rrt(113) * density(43)**2 
  rrt(114) = rrt(114) * density(01) * density(42) 
  rrt(115) = rrt(115) * density(01) * density(43) 
  rrt(116) = rrt(116) * density(21) * density(43) 
  rrt(117) = rrt(117) * density(40) * density(43) 
  rrt(118) = rrt(118) * density(14) * density(43) 
  rrt(119) = rrt(119) * density(29) * density(43) 
  rrt(120) = rrt(120) * density(01) * density(43) 
  rrt(121) = rrt(121) * density(21) * density(43) 
  rrt(122) = rrt(122) * density(40) * density(43) 
  rrt(123) = rrt(123) * density(14) * density(43) 
  rrt(124) = rrt(124) * density(29) * density(43) 
  rrt(125) = rrt(125) * density(44) 
  rrt(126) = rrt(126) * density(01) * density(14) * density(29) 
  rrt(127) = rrt(127) * density(14) * density(21) * density(29) 
  rrt(128) = rrt(128) * density(14)**2 * density(29) 
  rrt(129) = rrt(129) * density(14) * density(29)**2 
  rrt(130) = rrt(130) * density(14) * density(29) * density(40) 
  rrt(131) = rrt(131) * density(01) * density(21) * density(29) 
  rrt(132) = rrt(132) * density(01) * density(29) * density(40) 
  rrt(133) = rrt(133) * density(21) * density(29) * density(40) 
  rrt(134) = rrt(134) * density(29) * density(40)**2 
  rrt(135) = rrt(135) * density(01) * density(29) * density(42) 
  rrt(136) = rrt(136) * density(21) * density(29) * density(42) 
  rrt(137) = rrt(137) * density(14) * density(29) * density(42) 
  rrt(138) = rrt(138) * density(29)**2 * density(42) 
  rrt(139) = rrt(139) * density(29) * density(40) * density(42) 
  rrt(140) = rrt(140) * density(42) * density(43) 
  rrt(141) = rrt(141) * density(17) * density(21) 
  rrt(142) = rrt(142) * density(17) * density(40) 
  rrt(143) = rrt(143) * density(17) * density(40) 
  rrt(144) = rrt(144) * density(17) * density(40) 
  rrt(145) = rrt(145) * density(33) * density(40) 
  rrt(146) = rrt(146) * density(33) * density(40) 
  rrt(147) = rrt(147) * density(33) * density(41) 
  rrt(148) = rrt(148) * density(33) * density(42) 
  rrt(149) = rrt(149) * density(18) * density(40) 
  rrt(150) = rrt(150) * density(01) * density(34) 
  rrt(151) = rrt(151) * density(34) * density(40) 
  rrt(152) = rrt(152) * density(34) * density(42) 
  rrt(153) = rrt(153) * density(34) * density(42) 
  rrt(154) = rrt(154) * density(19) * density(40) 
  rrt(155) = rrt(155) * density(19) * density(40) 
  rrt(156) = rrt(156) * density(40) * density(47) 
  rrt(157) = rrt(157) * density(40) * density(46) 
  rrt(158) = rrt(158) * density(20) * density(40) 
  rrt(159) = rrt(159) * density(35) * density(40) 
  rrt(160) = rrt(160) * density(36) * density(42) 
  rrt(161) = rrt(161) * density(36) * density(41) 
  rrt(162) = rrt(162) * density(37) * density(42) 
  rrt(163) = rrt(163) * density(38) * density(40) 
  rrt(164) = rrt(164) * density(38) * density(40) 
  rrt(165) = rrt(165) * density(38) * density(42) 
  rrt(166) = rrt(166) * density(38) * density(42) 
  rrt(167) = rrt(167) * density(21) * density(48) 
  rrt(168) = rrt(168) * density(42) * density(48) 
  rrt(169) = rrt(169) * density(42) * density(50) 
  rrt(170) = rrt(170) * density(43) * density(50) 
  rrt(171) = rrt(171) * density(44) * density(50) 
  rrt(172) = rrt(172) * density(40) * density(51) 
  rrt(173) = rrt(173) * density(39) * density(40) 
  rrt(174) = rrt(174) * density(36) * density(40) 
  rrt(175) = rrt(175) * density(36) * density(45) 
  rrt(176) = rrt(176) * density(36) * density(47) 
  rrt(177) = rrt(177) * density(37) * density(45) 
  rrt(178) = rrt(178) * density(37) * density(47) 
  rrt(179) = rrt(179) * density(38) * density(45) 
  rrt(180) = rrt(180) * density(38) * density(47) 
  rrt(181) = rrt(181) * density(17) * density(48) 
  rrt(182) = rrt(182) * density(18) * density(48) 
  rrt(183) = rrt(183) * density(33) * density(48) 
  rrt(184) = rrt(184) * density(34) * density(48) 
  rrt(185) = rrt(185) * density(45) * density(48) 
  rrt(186) = rrt(186) * density(46) * density(48) 
  rrt(187) = rrt(187) * density(47) * density(48) 
  rrt(188) = rrt(188) * density(45) * density(49) 
  rrt(189) = rrt(189) * density(47) * density(49) 
  rrt(190) = rrt(190) * density(17) * density(50) 
  rrt(191) = rrt(191) * density(18) * density(50) 
  rrt(192) = rrt(192) * density(33) * density(50) 
  rrt(193) = rrt(193) * density(34) * density(50) 
  rrt(194) = rrt(194) * density(45) * density(50) 
  rrt(195) = rrt(195) * density(46) * density(50) 
  rrt(196) = rrt(196) * density(47) * density(50) 
  rrt(197) = rrt(197) * density(18) * density(51) 
  rrt(198) = rrt(198) * density(45) * density(51) 
  rrt(199) = rrt(199) * density(47) * density(51) 
  rrt(200) = rrt(200) * density(37) * density(52) 
  rrt(201) = rrt(201) * density(18) * density(48) 
  rrt(202) = rrt(202) * density(19) * density(48) 
  rrt(203) = rrt(203) * density(20) * density(48) 
  rrt(204) = rrt(204) * density(34) * density(48) 
  rrt(205) = rrt(205) * density(35) * density(48) 
  rrt(206) = rrt(206) * density(45) * density(48) 
  rrt(207) = rrt(207) * density(46) * density(48) 
  rrt(208) = rrt(208) * density(47) * density(48) 
  rrt(209) = rrt(209) * density(48) * density(52) 
  rrt(210) = rrt(210) * density(18) * density(49) 
  rrt(211) = rrt(211) * density(18) * density(50) 
  rrt(212) = rrt(212) * density(19) * density(50) 
  rrt(213) = rrt(213) * density(20) * density(50) 
  rrt(214) = rrt(214) * density(34) * density(50) 
  rrt(215) = rrt(215) * density(35) * density(50) 
  rrt(216) = rrt(216) * density(45) * density(50) 
  rrt(217) = rrt(217) * density(46) * density(50) 
  rrt(218) = rrt(218) * density(47) * density(50) 
  rrt(219) = rrt(219) * density(50) * density(52) 
  rrt(220) = rrt(220) * density(39) * density(45) 
  rrt(221) = rrt(221) * density(39) * density(47) 
  rrt(222) = rrt(222) * density(36) * density(45) 
  rrt(223) = rrt(223) * density(37) * density(45) 
  rrt(224) = rrt(224) * density(17) * density(36) 
  rrt(225) = rrt(225) * density(36) * density(45) 
  rrt(226) = rrt(226) * density(17) * density(37) 
  rrt(227) = rrt(227) * density(38) * density(45) 
  rrt(228) = rrt(228) * density(38) * density(47) 
  rrt(229) = rrt(229) * density(17) * density(48) 
  rrt(230) = rrt(230) * density(18) * density(48) 
  rrt(231) = rrt(231) * density(33) * density(48) 
  rrt(232) = rrt(232) * density(34) * density(48) 
  rrt(233) = rrt(233) * density(45) * density(48) 
  rrt(234) = rrt(234) * density(46) * density(48) 
  rrt(235) = rrt(235) * density(47) * density(48) 
  rrt(236) = rrt(236) * density(45) * density(49) 
  rrt(237) = rrt(237) * density(47) * density(49) 
  rrt(238) = rrt(238) * density(17) * density(50) 
  rrt(239) = rrt(239) * density(18) * density(50) 
  rrt(240) = rrt(240) * density(33) * density(50) 
  rrt(241) = rrt(241) * density(34) * density(50) 
  rrt(242) = rrt(242) * density(45) * density(50) 
  rrt(243) = rrt(243) * density(46) * density(50) 
  rrt(244) = rrt(244) * density(47) * density(50) 
  rrt(245) = rrt(245) * density(45) * density(51) 
  rrt(246) = rrt(246) * density(47) * density(51) 
  ydot(01) = -rrt(001)-rrt(002)-rrt(003)-rrt(004)-rrt(005)-rrt(006)-rrt(007)-rrt(008)-rrt(009)+rrt(010)+rrt(011)+rrt(012)+rrt(013)&
             +rrt(014)+rrt(015)+rrt(016)+rrt(017)+rrt(018)+rrt(026)-rrt(034)-rrt(035)-rrt(036)-rrt(037)-rrt(038)-rrt(039)-rrt(040)&
             -rrt(041)-rrt(042)-rrt(043)-rrt(044)-rrt(045)+rrt(054)+rrt(061)+rrt(062)+rrt(063)+rrt(064)+rrt(067)+rrt(068)+rrt(069)&
             +rrt(070)+rrt(071)+rrt(072)+rrt(073)+rrt(074)+rrt(076)+rrt(077)+rrt(081)+rrt(082)+rrt(084)+rrt(085)+rrt(092)+rrt(094)&
             +rrt(096)+rrt(105)+rrt(144)+rrt(149)-rrt(150)+rrt(154)+rrt(155)+  2.d0 * rrt(158)+rrt(182)+rrt(191)+rrt(197)+rrt(200)&
             +rrt(202)+  2.d0 * rrt(203)+rrt(207)+rrt(209)+rrt(212)+  2.d0 * rrt(213)+rrt(217)+rrt(219)+rrt(230)+rrt(239) 
  ydot(02) = +rrt(001)+rrt(002)-rrt(010)-rrt(018)+rrt(019)-rrt(026)+rrt(027) 
  ydot(03) = +rrt(003)-rrt(011)-rrt(019)+rrt(020)-rrt(027)+rrt(028) 
  ydot(04) = +rrt(004)-rrt(012)-rrt(020)+rrt(021)-rrt(028)+rrt(029) 
  ydot(05) = +rrt(005)-rrt(013)-rrt(021)+rrt(022)-rrt(029)+rrt(030) 
  ydot(06) = +rrt(006)-rrt(014)-rrt(022)+rrt(023)-rrt(030)+rrt(031) 
  ydot(07) = +rrt(007)-rrt(015)-rrt(023)+rrt(024)-rrt(031)+rrt(032) 
  ydot(08) = +rrt(008)-rrt(016)-rrt(024)+rrt(025)-rrt(032)+rrt(033) 
  ydot(09) = +rrt(009)-rrt(017)-rrt(025)-rrt(033) 
  ydot(10) = +rrt(034)+rrt(035)+rrt(036)-rrt(061)-rrt(062)-rrt(066)-rrt(067)-rrt(068)-rrt(069)-rrt(070)-rrt(071)-rrt(072)&
             -  2.d0 * rrt(073)-  2.d0 * rrt(074)+rrt(075)+rrt(078)+rrt(087) 
  ydot(11) = +rrt(037)+rrt(038)+rrt(039)-rrt(063)-rrt(064)+rrt(065)+rrt(073)-rrt(075)-rrt(076)-rrt(077)-rrt(078)+rrt(080) 
  ydot(12) = +rrt(040)+rrt(041)+rrt(042)+rrt(079)-rrt(080)-rrt(081)-rrt(082) 
  ydot(13) = +rrt(043)+rrt(044)-rrt(065)+rrt(074)-rrt(079) 
  ydot(14) = +rrt(045)-rrt(050)-rrt(052)-rrt(053)-rrt(054)-rrt(055)-rrt(056)+rrt(071)+rrt(082)-rrt(088)+rrt(090)-rrt(092)-rrt(093)&
             -rrt(094)-rrt(095)-rrt(096)-rrt(097)+rrt(098)-rrt(126)-rrt(127)-rrt(128)-rrt(129)-rrt(130)+rrt(142)+rrt(146)+rrt(154)&
             +rrt(181)+rrt(190)+  2.d0 * rrt(201)+rrt(202)+rrt(206)+rrt(208)+  2.d0 * rrt(210)+  2.d0 * rrt(211)+rrt(212)+rrt(216)&
             +rrt(218)+rrt(229)+rrt(238) 
  ydot(15) = +rrt(045)+rrt(066)-rrt(083)-rrt(084)-rrt(085) 
  ydot(16) = -rrt(086)-rrt(087) 
  ydot(17) = -rrt(141)-rrt(142)-rrt(143)-rrt(144)-rrt(181)-rrt(190)-rrt(224)-rrt(226)-rrt(229)-rrt(238) 
  ydot(18) = +rrt(143)-rrt(149)-rrt(182)-rrt(191)-rrt(197)-rrt(201)-rrt(210)-rrt(211)-rrt(230)-rrt(239) 
  ydot(19) = -rrt(154)-rrt(155)-rrt(202)-rrt(212) 
  ydot(20) = -rrt(158)-rrt(203)-rrt(213) 
  ydot(21) = +rrt(053)+rrt(059)-rrt(067)-rrt(068)-rrt(077)-rrt(081)-rrt(083)-rrt(086)+rrt(090)-rrt(093)+rrt(096)+rrt(098)+rrt(101)&
             +rrt(102)-rrt(103)+rrt(104)-rrt(107)+rrt(108)+rrt(110)+rrt(111)-rrt(112)+rrt(113)+rrt(120)+rrt(121)+rrt(122)+rrt(123)&
             +rrt(124)-rrt(131)-rrt(141)+rrt(151)+rrt(153)+  2.d0 * rrt(159)+rrt(162)+rrt(164)+rrt(166)-rrt(167)+rrt(173)+rrt(177)&
             +rrt(178)+rrt(184)+rrt(193)+  2.d0 * rrt(200)+  2.d0 * rrt(205)+rrt(208)+rrt(209)+  2.d0 * rrt(215)+rrt(218)+rrt(219)&
             +  2.d0 * rrt(220)+  2.d0 * rrt(221)+rrt(223)+rrt(232)+rrt(241) 
  ydot(22) = 0.0d0
  ydot(23) = 0.0d0
  ydot(24) = 0.0d0
  ydot(25) = 0.0d0
  ydot(26) = -rrt(088) 
  ydot(27) = +rrt(068) 
  ydot(28) = 0.0d0
  ydot(29) = +rrt(047)-rrt(057)-rrt(058)-rrt(059)-rrt(060)-rrt(066)+rrt(067)+rrt(072)+  2.d0 * rrt(077)+  2.d0 * rrt(081)+rrt(082)&
             +rrt(083)+rrt(084)+rrt(086)+rrt(087)+rrt(088)+rrt(089)+rrt(092)+rrt(093)+  2.d0 * rrt(094)+rrt(095)-rrt(098)-rrt(099)&
             -rrt(100)-rrt(101)-rrt(102)+rrt(103)+rrt(114)+rrt(115)+rrt(116)+rrt(117)+rrt(118)+rrt(119)-rrt(126)-rrt(127)-rrt(128)&
             -rrt(129)-rrt(130)-rrt(131)-rrt(132)-rrt(133)-rrt(134)-rrt(135)-rrt(136)-rrt(137)-rrt(138)-rrt(139)+rrt(143)+rrt(145)&
             +rrt(148)+rrt(160)+rrt(163)+rrt(175)+rrt(176)+rrt(183)+rrt(192)+  2.d0 * rrt(204)+rrt(206)+rrt(207)+  2.d0 * rrt(214)&
             +rrt(216)+rrt(217)+rrt(222)+rrt(231)+rrt(240) 
  ydot(30) = +rrt(067)-rrt(089)-rrt(090)-rrt(091) 
  ydot(31) = 0.0d0
  ydot(32) = +rrt(060)-rrt(104)+rrt(107)-rrt(110)+rrt(112)+rrt(131)+rrt(152)+rrt(165)+rrt(179)+rrt(180)+rrt(227)+rrt(228) 
  ydot(33) = +rrt(141)+rrt(144)-rrt(145)-rrt(146)-rrt(147)-rrt(148)-rrt(183)-rrt(192)-rrt(231)-rrt(240) 
  ydot(34) = +rrt(146)-rrt(150)-rrt(151)-rrt(152)-rrt(153)-rrt(184)-rrt(193)-rrt(204)-rrt(214)-rrt(232)-rrt(241) 
  ydot(35) = -rrt(159)-rrt(205)-rrt(215) 
  ydot(36) = +rrt(048)-rrt(050)-rrt(051)-rrt(160)-rrt(161)-rrt(174)-rrt(175)-rrt(176)-rrt(222)-rrt(224)-rrt(225) 
  ydot(37) = -rrt(052)-rrt(162)+rrt(167)-rrt(177)-rrt(178)-rrt(200)-rrt(223)-rrt(226) 
  ydot(38) = -rrt(053)-rrt(163)-rrt(164)-rrt(165)-rrt(166)-rrt(179)-rrt(180)-rrt(227)-rrt(228) 
  ydot(39) = -rrt(173)-rrt(220)-rrt(221) 
  ydot(40) = -rrt(046)+rrt(047)+rrt(048)-rrt(049)+rrt(050)-rrt(051)+rrt(053)+rrt(054)+  2.d0 * rrt(055)+rrt(056)+  2.d0 * rrt(058)&
             +rrt(059)+rrt(060)+rrt(061)+rrt(063)+rrt(066)+rrt(071)+rrt(072)-rrt(082)+rrt(083)-rrt(084)+rrt(085)+rrt(086)-rrt(087)&
             +rrt(088)-rrt(090)+  2.d0 * rrt(091)-rrt(092)+rrt(093)+  2.d0 * rrt(097)-rrt(098)-rrt(099)+  2.d0 * rrt(100)+rrt(101)&
             -rrt(103)-rrt(104)-rrt(105)-rrt(106)+rrt(107)+  2.d0 * rrt(108)+rrt(109)+rrt(111)+rrt(114)+rrt(120)+rrt(121)+rrt(122)&
             +rrt(123)+rrt(124)+rrt(126)+rrt(127)+rrt(128)+rrt(129)+rrt(130)-rrt(132)-rrt(133)-rrt(134)+rrt(141)-rrt(142)-rrt(143)&
             -rrt(144)-rrt(145)-rrt(146)+rrt(147)-rrt(149)+rrt(150)-rrt(151)-rrt(154)-rrt(155)-rrt(156)-rrt(157)-rrt(158)-rrt(159)&
             +rrt(161)-rrt(163)-rrt(164)+rrt(167)+rrt(168)+rrt(169)-rrt(172)-rrt(173)-rrt(174)+rrt(175)+rrt(177)+rrt(179)+rrt(181)&
             +rrt(182)+rrt(183)+rrt(184)+  2.d0 * rrt(185)+rrt(186)+rrt(187)+rrt(188)+rrt(194)+rrt(198)+rrt(201)+rrt(202)+rrt(203)&
             +rrt(204)+rrt(205)+rrt(206)+rrt(207)+rrt(208)+rrt(209)+rrt(220)+rrt(222)+rrt(223)+rrt(224)+rrt(227)+rrt(229)+rrt(230)&
             +rrt(231)+rrt(232)+  2.d0 * rrt(233)+rrt(234)+rrt(235)+rrt(236)+rrt(242)+rrt(245) 
  ydot(41) = -rrt(071)-rrt(085)-rrt(091)+rrt(095)-rrt(100)-rrt(105)-rrt(147)+rrt(157)-rrt(161)+rrt(186)+rrt(188)+rrt(189)+rrt(195)&
             +rrt(210)+rrt(234)+rrt(236)+rrt(237)+rrt(243) 
  ydot(42) = -rrt(048)+rrt(051)+rrt(052)+rrt(056)+rrt(057)+rrt(062)+rrt(064)-rrt(072)-rrt(094)-rrt(095)-rrt(096)-rrt(097)+rrt(099)&
             -rrt(101)+rrt(102)+rrt(103)+rrt(104)+rrt(105)+  2.d0 * rrt(106)-rrt(107)-  2.d0 * rrt(108)-  2.d0 * rrt(109)-rrt(110)&
             +rrt(112)+  2.d0 * rrt(113)-rrt(114)+rrt(115)+rrt(116)+rrt(117)+rrt(118)+rrt(119)+rrt(125)+rrt(132)+rrt(133)+rrt(134)&
             -rrt(135)-rrt(136)-rrt(137)-rrt(138)-rrt(139)-rrt(140)-rrt(148)-rrt(152)-rrt(153)+rrt(156)-rrt(160)-rrt(162)-rrt(165)&
             -rrt(166)-rrt(168)-rrt(169)+rrt(170)+  2.d0 * rrt(171)+rrt(172)+rrt(176)+rrt(178)+rrt(180)+rrt(187)+rrt(189)+rrt(190)&
             +rrt(191)+rrt(192)+rrt(193)+rrt(194)+rrt(195)+  2.d0 * rrt(196)+rrt(199)+rrt(211)+rrt(212)+rrt(213)+rrt(214)+rrt(215)&
             +rrt(216)+rrt(217)+rrt(218)+rrt(219)+rrt(221)+rrt(225)+rrt(226)+rrt(228)+rrt(235)+rrt(237)+rrt(238)+rrt(239)+rrt(240)&
             +rrt(241)+rrt(242)+rrt(243)+  2.d0 * rrt(244)+rrt(246) 
  ydot(43) = -rrt(102)-rrt(106)+rrt(109)+rrt(110)-rrt(111)-rrt(112)-  2.d0 * rrt(113)-rrt(115)-rrt(116)-rrt(117)-rrt(118)-rrt(119)&
             -rrt(120)-rrt(121)-rrt(122)-rrt(123)-rrt(124)+rrt(125)+rrt(135)+rrt(136)+rrt(137)+rrt(138)+rrt(139)-rrt(140)-rrt(170)&
             +rrt(197)+rrt(198)+rrt(199)+rrt(245)+rrt(246) 
  ydot(44) = -rrt(125)+rrt(140)-rrt(171) 
  ydot(45) = +rrt(046)+rrt(142)+rrt(145)+rrt(147)+rrt(149)+rrt(150)+rrt(151)+rrt(152)+rrt(154)+rrt(156)+rrt(157)+rrt(158)+rrt(159)&
             -rrt(175)-rrt(177)-rrt(179)-rrt(185)-rrt(188)-rrt(194)-rrt(198)-rrt(206)-rrt(216)-rrt(220)-rrt(222)-rrt(223)-rrt(225)&
             -rrt(227)-rrt(233)-rrt(236)-rrt(242)-rrt(245) 
  ydot(46) = +rrt(155)-rrt(157)-rrt(186)-rrt(195)-rrt(207)-rrt(217)-rrt(234)-rrt(243) 
  ydot(47) = -rrt(047)+rrt(148)+rrt(153)-rrt(156)-rrt(176)-rrt(178)-rrt(180)-rrt(187)-rrt(189)-rrt(196)-rrt(199)-rrt(208)-rrt(218)&
             -rrt(221)-rrt(228)-rrt(235)-rrt(237)-rrt(244)-rrt(246) 
  ydot(48) = +rrt(049)-rrt(057)-rrt(061)-rrt(063)+rrt(161)-rrt(167)-rrt(168)-rrt(181)-rrt(182)-rrt(183)-rrt(184)-rrt(185)-rrt(186)&
             -rrt(187)-rrt(201)-rrt(202)-rrt(203)-rrt(204)-rrt(205)-rrt(206)-rrt(207)-rrt(208)-rrt(209)-rrt(229)-rrt(230)-rrt(231)&
             -rrt(232)-rrt(233)-rrt(234)-rrt(235) 
  ydot(49) = -rrt(054)-rrt(058)-rrt(188)-rrt(189)-rrt(210)-rrt(236)-rrt(237) 
  ydot(50) = -rrt(055)-rrt(059)-rrt(062)-rrt(064)+rrt(160)+rrt(162)+rrt(164)+rrt(165)+rrt(168)-rrt(169)-rrt(170)-rrt(171)+rrt(172)&
             +rrt(174)-rrt(190)-rrt(191)-rrt(192)-rrt(193)-rrt(194)-rrt(195)-rrt(196)-rrt(211)-rrt(212)-rrt(213)-rrt(214)-rrt(215)&
             -rrt(216)-rrt(217)-rrt(218)-rrt(219)-rrt(238)-rrt(239)-rrt(240)-rrt(241)-rrt(242)-rrt(243)-rrt(244) 
  ydot(51) = -rrt(056)-rrt(060)+rrt(163)+rrt(166)+rrt(169)+rrt(170)+rrt(171)-rrt(172)+rrt(173)-rrt(197)-rrt(198)-rrt(199)-rrt(245)&
             -rrt(246) 
  ydot(52) = -rrt(200)-rrt(209)-rrt(219) 
  ydot(53) = +rrt(046)-rrt(047)-rrt(048)-rrt(049)+rrt(050)+rrt(051)+rrt(052)+rrt(053)+rrt(054)+rrt(055)+rrt(056)+rrt(057)+rrt(058)&
             +rrt(059)+rrt(060)+rrt(061)+rrt(062)+rrt(063)+rrt(064) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(54) = 0.0d0
  if( lgas_heating ) then
    ydot(54) = ( ZDPlasKin_cfg(14)/k_B + ydot(54) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(54) = ydot(54) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_fex
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction jacobian
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_jex(neq,t,y,ml,mu,pd,nrpd)
  implicit none
  integer,          intent(in)  :: neq, ml, mu, nrpd
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: pd(nrpd,neq)
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(54)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(01,01) = pd(01,01) - rrt(001) * density(53) 
  pd(01,53) = pd(01,53) - rrt(001) * density(01) 
  pd(02,01) = pd(02,01) + rrt(001) * density(53) 
  pd(02,53) = pd(02,53) + rrt(001) * density(01) 
  pd(01,01) = pd(01,01) - rrt(002) * density(53) 
  pd(01,53) = pd(01,53) - rrt(002) * density(01) 
  pd(02,01) = pd(02,01) + rrt(002) * density(53) 
  pd(02,53) = pd(02,53) + rrt(002) * density(01) 
  pd(01,01) = pd(01,01) - rrt(003) * density(53) 
  pd(01,53) = pd(01,53) - rrt(003) * density(01) 
  pd(03,01) = pd(03,01) + rrt(003) * density(53) 
  pd(03,53) = pd(03,53) + rrt(003) * density(01) 
  pd(01,01) = pd(01,01) - rrt(004) * density(53) 
  pd(01,53) = pd(01,53) - rrt(004) * density(01) 
  pd(04,01) = pd(04,01) + rrt(004) * density(53) 
  pd(04,53) = pd(04,53) + rrt(004) * density(01) 
  pd(01,01) = pd(01,01) - rrt(005) * density(53) 
  pd(01,53) = pd(01,53) - rrt(005) * density(01) 
  pd(05,01) = pd(05,01) + rrt(005) * density(53) 
  pd(05,53) = pd(05,53) + rrt(005) * density(01) 
  pd(01,01) = pd(01,01) - rrt(006) * density(53) 
  pd(01,53) = pd(01,53) - rrt(006) * density(01) 
  pd(06,01) = pd(06,01) + rrt(006) * density(53) 
  pd(06,53) = pd(06,53) + rrt(006) * density(01) 
  pd(01,01) = pd(01,01) - rrt(007) * density(53) 
  pd(01,53) = pd(01,53) - rrt(007) * density(01) 
  pd(07,01) = pd(07,01) + rrt(007) * density(53) 
  pd(07,53) = pd(07,53) + rrt(007) * density(01) 
  pd(01,01) = pd(01,01) - rrt(008) * density(53) 
  pd(01,53) = pd(01,53) - rrt(008) * density(01) 
  pd(08,01) = pd(08,01) + rrt(008) * density(53) 
  pd(08,53) = pd(08,53) + rrt(008) * density(01) 
  pd(01,01) = pd(01,01) - rrt(009) * density(53) 
  pd(01,53) = pd(01,53) - rrt(009) * density(01) 
  pd(09,01) = pd(09,01) + rrt(009) * density(53) 
  pd(09,53) = pd(09,53) + rrt(009) * density(01) 
  pd(01,02) = pd(01,02) + rrt(010) * density(53) 
  pd(01,53) = pd(01,53) + rrt(010) * density(02) 
  pd(02,02) = pd(02,02) - rrt(010) * density(53) 
  pd(02,53) = pd(02,53) - rrt(010) * density(02) 
  pd(01,03) = pd(01,03) + rrt(011) * density(53) 
  pd(01,53) = pd(01,53) + rrt(011) * density(03) 
  pd(03,03) = pd(03,03) - rrt(011) * density(53) 
  pd(03,53) = pd(03,53) - rrt(011) * density(03) 
  pd(01,04) = pd(01,04) + rrt(012) * density(53) 
  pd(01,53) = pd(01,53) + rrt(012) * density(04) 
  pd(04,04) = pd(04,04) - rrt(012) * density(53) 
  pd(04,53) = pd(04,53) - rrt(012) * density(04) 
  pd(01,05) = pd(01,05) + rrt(013) * density(53) 
  pd(01,53) = pd(01,53) + rrt(013) * density(05) 
  pd(05,05) = pd(05,05) - rrt(013) * density(53) 
  pd(05,53) = pd(05,53) - rrt(013) * density(05) 
  pd(01,06) = pd(01,06) + rrt(014) * density(53) 
  pd(01,53) = pd(01,53) + rrt(014) * density(06) 
  pd(06,06) = pd(06,06) - rrt(014) * density(53) 
  pd(06,53) = pd(06,53) - rrt(014) * density(06) 
  pd(01,07) = pd(01,07) + rrt(015) * density(53) 
  pd(01,53) = pd(01,53) + rrt(015) * density(07) 
  pd(07,07) = pd(07,07) - rrt(015) * density(53) 
  pd(07,53) = pd(07,53) - rrt(015) * density(07) 
  pd(01,08) = pd(01,08) + rrt(016) * density(53) 
  pd(01,53) = pd(01,53) + rrt(016) * density(08) 
  pd(08,08) = pd(08,08) - rrt(016) * density(53) 
  pd(08,53) = pd(08,53) - rrt(016) * density(08) 
  pd(01,09) = pd(01,09) + rrt(017) * density(53) 
  pd(01,53) = pd(01,53) + rrt(017) * density(09) 
  pd(09,09) = pd(09,09) - rrt(017) * density(53) 
  pd(09,53) = pd(09,53) - rrt(017) * density(09) 
  pd(01,01) = pd(01,01) + rrt(018) * density(02) 
  pd(01,02) = pd(01,02) + rrt(018) * density(01) 
  pd(02,01) = pd(02,01) - rrt(018) * density(02) 
  pd(02,02) = pd(02,02) - rrt(018) * density(01) 
  pd(02,01) = pd(02,01) + rrt(019) * density(03) 
  pd(02,03) = pd(02,03) + rrt(019) * density(01) 
  pd(03,01) = pd(03,01) - rrt(019) * density(03) 
  pd(03,03) = pd(03,03) - rrt(019) * density(01) 
  pd(03,01) = pd(03,01) + rrt(020) * density(04) 
  pd(03,04) = pd(03,04) + rrt(020) * density(01) 
  pd(04,01) = pd(04,01) - rrt(020) * density(04) 
  pd(04,04) = pd(04,04) - rrt(020) * density(01) 
  pd(04,01) = pd(04,01) + rrt(021) * density(05) 
  pd(04,05) = pd(04,05) + rrt(021) * density(01) 
  pd(05,01) = pd(05,01) - rrt(021) * density(05) 
  pd(05,05) = pd(05,05) - rrt(021) * density(01) 
  pd(05,01) = pd(05,01) + rrt(022) * density(06) 
  pd(05,06) = pd(05,06) + rrt(022) * density(01) 
  pd(06,01) = pd(06,01) - rrt(022) * density(06) 
  pd(06,06) = pd(06,06) - rrt(022) * density(01) 
  pd(06,01) = pd(06,01) + rrt(023) * density(07) 
  pd(06,07) = pd(06,07) + rrt(023) * density(01) 
  pd(07,01) = pd(07,01) - rrt(023) * density(07) 
  pd(07,07) = pd(07,07) - rrt(023) * density(01) 
  pd(07,01) = pd(07,01) + rrt(024) * density(08) 
  pd(07,08) = pd(07,08) + rrt(024) * density(01) 
  pd(08,01) = pd(08,01) - rrt(024) * density(08) 
  pd(08,08) = pd(08,08) - rrt(024) * density(01) 
  pd(08,01) = pd(08,01) + rrt(025) * density(09) 
  pd(08,09) = pd(08,09) + rrt(025) * density(01) 
  pd(09,01) = pd(09,01) - rrt(025) * density(09) 
  pd(09,09) = pd(09,09) - rrt(025) * density(01) 
  pd(01,02) = pd(01,02) + rrt(026) * density(29) 
  pd(01,29) = pd(01,29) + rrt(026) * density(02) 
  pd(02,02) = pd(02,02) - rrt(026) * density(29) 
  pd(02,29) = pd(02,29) - rrt(026) * density(02) 
  pd(02,03) = pd(02,03) + rrt(027) * density(29) 
  pd(02,29) = pd(02,29) + rrt(027) * density(03) 
  pd(03,03) = pd(03,03) - rrt(027) * density(29) 
  pd(03,29) = pd(03,29) - rrt(027) * density(03) 
  pd(03,04) = pd(03,04) + rrt(028) * density(29) 
  pd(03,29) = pd(03,29) + rrt(028) * density(04) 
  pd(04,04) = pd(04,04) - rrt(028) * density(29) 
  pd(04,29) = pd(04,29) - rrt(028) * density(04) 
  pd(04,05) = pd(04,05) + rrt(029) * density(29) 
  pd(04,29) = pd(04,29) + rrt(029) * density(05) 
  pd(05,05) = pd(05,05) - rrt(029) * density(29) 
  pd(05,29) = pd(05,29) - rrt(029) * density(05) 
  pd(05,06) = pd(05,06) + rrt(030) * density(29) 
  pd(05,29) = pd(05,29) + rrt(030) * density(06) 
  pd(06,06) = pd(06,06) - rrt(030) * density(29) 
  pd(06,29) = pd(06,29) - rrt(030) * density(06) 
  pd(06,07) = pd(06,07) + rrt(031) * density(29) 
  pd(06,29) = pd(06,29) + rrt(031) * density(07) 
  pd(07,07) = pd(07,07) - rrt(031) * density(29) 
  pd(07,29) = pd(07,29) - rrt(031) * density(07) 
  pd(07,08) = pd(07,08) + rrt(032) * density(29) 
  pd(07,29) = pd(07,29) + rrt(032) * density(08) 
  pd(08,08) = pd(08,08) - rrt(032) * density(29) 
  pd(08,29) = pd(08,29) - rrt(032) * density(08) 
  pd(08,09) = pd(08,09) + rrt(033) * density(29) 
  pd(08,29) = pd(08,29) + rrt(033) * density(09) 
  pd(09,09) = pd(09,09) - rrt(033) * density(29) 
  pd(09,29) = pd(09,29) - rrt(033) * density(09) 
  pd(01,01) = pd(01,01) - rrt(034) * density(53) 
  pd(01,53) = pd(01,53) - rrt(034) * density(01) 
  pd(10,01) = pd(10,01) + rrt(034) * density(53) 
  pd(10,53) = pd(10,53) + rrt(034) * density(01) 
  pd(01,01) = pd(01,01) - rrt(035) * density(53) 
  pd(01,53) = pd(01,53) - rrt(035) * density(01) 
  pd(10,01) = pd(10,01) + rrt(035) * density(53) 
  pd(10,53) = pd(10,53) + rrt(035) * density(01) 
  pd(01,01) = pd(01,01) - rrt(036) * density(53) 
  pd(01,53) = pd(01,53) - rrt(036) * density(01) 
  pd(10,01) = pd(10,01) + rrt(036) * density(53) 
  pd(10,53) = pd(10,53) + rrt(036) * density(01) 
  pd(01,01) = pd(01,01) - rrt(037) * density(53) 
  pd(01,53) = pd(01,53) - rrt(037) * density(01) 
  pd(11,01) = pd(11,01) + rrt(037) * density(53) 
  pd(11,53) = pd(11,53) + rrt(037) * density(01) 
  pd(01,01) = pd(01,01) - rrt(038) * density(53) 
  pd(01,53) = pd(01,53) - rrt(038) * density(01) 
  pd(11,01) = pd(11,01) + rrt(038) * density(53) 
  pd(11,53) = pd(11,53) + rrt(038) * density(01) 
  pd(01,01) = pd(01,01) - rrt(039) * density(53) 
  pd(01,53) = pd(01,53) - rrt(039) * density(01) 
  pd(11,01) = pd(11,01) + rrt(039) * density(53) 
  pd(11,53) = pd(11,53) + rrt(039) * density(01) 
  pd(01,01) = pd(01,01) - rrt(040) * density(53) 
  pd(01,53) = pd(01,53) - rrt(040) * density(01) 
  pd(12,01) = pd(12,01) + rrt(040) * density(53) 
  pd(12,53) = pd(12,53) + rrt(040) * density(01) 
  pd(01,01) = pd(01,01) - rrt(041) * density(53) 
  pd(01,53) = pd(01,53) - rrt(041) * density(01) 
  pd(12,01) = pd(12,01) + rrt(041) * density(53) 
  pd(12,53) = pd(12,53) + rrt(041) * density(01) 
  pd(01,01) = pd(01,01) - rrt(042) * density(53) 
  pd(01,53) = pd(01,53) - rrt(042) * density(01) 
  pd(12,01) = pd(12,01) + rrt(042) * density(53) 
  pd(12,53) = pd(12,53) + rrt(042) * density(01) 
  pd(01,01) = pd(01,01) - rrt(043) * density(53) 
  pd(01,53) = pd(01,53) - rrt(043) * density(01) 
  pd(13,01) = pd(13,01) + rrt(043) * density(53) 
  pd(13,53) = pd(13,53) + rrt(043) * density(01) 
  pd(01,01) = pd(01,01) - rrt(044) * density(53) 
  pd(01,53) = pd(01,53) - rrt(044) * density(01) 
  pd(13,01) = pd(13,01) + rrt(044) * density(53) 
  pd(13,53) = pd(13,53) + rrt(044) * density(01) 
  pd(01,01) = pd(01,01) - rrt(045) * density(53) 
  pd(01,53) = pd(01,53) - rrt(045) * density(01) 
  pd(14,01) = pd(14,01) + rrt(045) * density(53) 
  pd(14,53) = pd(14,53) + rrt(045) * density(01) 
  pd(15,01) = pd(15,01) + rrt(045) * density(53) 
  pd(15,53) = pd(15,53) + rrt(045) * density(01) 
  pd(40,40) = pd(40,40) - rrt(046) * density(53) 
  pd(40,53) = pd(40,53) - rrt(046) * density(40) 
  pd(45,40) = pd(45,40) + rrt(046) * density(53) 
  pd(45,53) = pd(45,53) + rrt(046) * density(40) 
  pd(53,40) = pd(53,40) + rrt(046) * density(53) 
  pd(53,53) = pd(53,53) + rrt(046) * density(40) 
  pd(29,47) = pd(29,47) + rrt(047) * density(53) 
  pd(29,53) = pd(29,53) + rrt(047) * density(47) 
  pd(40,47) = pd(40,47) + rrt(047) * density(53) 
  pd(40,53) = pd(40,53) + rrt(047) * density(47) 
  pd(47,47) = pd(47,47) - rrt(047) * density(53) 
  pd(47,53) = pd(47,53) - rrt(047) * density(47) 
  pd(53,47) = pd(53,47) - rrt(047) * density(53) 
  pd(53,53) = pd(53,53) - rrt(047) * density(47) 
  pd(36,42) = pd(36,42) + rrt(048) * density(53) 
  pd(36,53) = pd(36,53) + rrt(048) * density(42) 
  pd(40,42) = pd(40,42) + rrt(048) * density(53) 
  pd(40,53) = pd(40,53) + rrt(048) * density(42) 
  pd(42,42) = pd(42,42) - rrt(048) * density(53) 
  pd(42,53) = pd(42,53) - rrt(048) * density(42) 
  pd(53,42) = pd(53,42) - rrt(048) * density(53) 
  pd(53,53) = pd(53,53) - rrt(048) * density(42) 
  pd(40,40) = pd(40,40) - rrt(049) * density(53) 
  pd(40,53) = pd(40,53) - rrt(049) * density(40) 
  pd(48,40) = pd(48,40) + rrt(049) * density(53) 
  pd(48,53) = pd(48,53) + rrt(049) * density(40) 
  pd(53,40) = pd(53,40) - rrt(049) * density(53) 
  pd(53,53) = pd(53,53) - rrt(049) * density(40) 
  pd(14,14) = pd(14,14) - rrt(050) * density(36) 
  pd(14,36) = pd(14,36) - rrt(050) * density(14) 
  pd(36,14) = pd(36,14) - rrt(050) * density(36) 
  pd(36,36) = pd(36,36) - rrt(050) * density(14) 
  pd(40,14) = pd(40,14) + rrt(050) * density(36) 
  pd(40,36) = pd(40,36) + rrt(050) * density(14) 
  pd(53,14) = pd(53,14) + rrt(050) * density(36) 
  pd(53,36) = pd(53,36) + rrt(050) * density(14) 
  pd(36,36) = pd(36,36) - rrt(051) * density(40) 
  pd(36,40) = pd(36,40) - rrt(051) * density(36) 
  pd(40,36) = pd(40,36) - rrt(051) * density(40) 
  pd(40,40) = pd(40,40) - rrt(051) * density(36) 
  pd(42,36) = pd(42,36) + rrt(051) * density(40) 
  pd(42,40) = pd(42,40) + rrt(051) * density(36) 
  pd(53,36) = pd(53,36) + rrt(051) * density(40) 
  pd(53,40) = pd(53,40) + rrt(051) * density(36) 
  pd(14,14) = pd(14,14) - rrt(052) * density(37) 
  pd(14,37) = pd(14,37) - rrt(052) * density(14) 
  pd(37,14) = pd(37,14) - rrt(052) * density(37) 
  pd(37,37) = pd(37,37) - rrt(052) * density(14) 
  pd(42,14) = pd(42,14) + rrt(052) * density(37) 
  pd(42,37) = pd(42,37) + rrt(052) * density(14) 
  pd(53,14) = pd(53,14) + rrt(052) * density(37) 
  pd(53,37) = pd(53,37) + rrt(052) * density(14) 
  pd(14,14) = pd(14,14) - rrt(053) * density(38) 
  pd(14,38) = pd(14,38) - rrt(053) * density(14) 
  pd(21,14) = pd(21,14) + rrt(053) * density(38) 
  pd(21,38) = pd(21,38) + rrt(053) * density(14) 
  pd(38,14) = pd(38,14) - rrt(053) * density(38) 
  pd(38,38) = pd(38,38) - rrt(053) * density(14) 
  pd(40,14) = pd(40,14) + rrt(053) * density(38) 
  pd(40,38) = pd(40,38) + rrt(053) * density(14) 
  pd(53,14) = pd(53,14) + rrt(053) * density(38) 
  pd(53,38) = pd(53,38) + rrt(053) * density(14) 
  pd(01,14) = pd(01,14) + rrt(054) * density(49) 
  pd(01,49) = pd(01,49) + rrt(054) * density(14) 
  pd(14,14) = pd(14,14) - rrt(054) * density(49) 
  pd(14,49) = pd(14,49) - rrt(054) * density(14) 
  pd(40,14) = pd(40,14) + rrt(054) * density(49) 
  pd(40,49) = pd(40,49) + rrt(054) * density(14) 
  pd(49,14) = pd(49,14) - rrt(054) * density(49) 
  pd(49,49) = pd(49,49) - rrt(054) * density(14) 
  pd(53,14) = pd(53,14) + rrt(054) * density(49) 
  pd(53,49) = pd(53,49) + rrt(054) * density(14) 
  pd(14,14) = pd(14,14) - rrt(055) * density(50) 
  pd(14,50) = pd(14,50) - rrt(055) * density(14) 
  pd(40,14) = pd(40,14) + rrt(055) * density(50) * 2.0d0
  pd(40,50) = pd(40,50) + rrt(055) * density(14) * 2.0d0
  pd(50,14) = pd(50,14) - rrt(055) * density(50) 
  pd(50,50) = pd(50,50) - rrt(055) * density(14) 
  pd(53,14) = pd(53,14) + rrt(055) * density(50) 
  pd(53,50) = pd(53,50) + rrt(055) * density(14) 
  pd(14,14) = pd(14,14) - rrt(056) * density(51) 
  pd(14,51) = pd(14,51) - rrt(056) * density(14) 
  pd(40,14) = pd(40,14) + rrt(056) * density(51) 
  pd(40,51) = pd(40,51) + rrt(056) * density(14) 
  pd(42,14) = pd(42,14) + rrt(056) * density(51) 
  pd(42,51) = pd(42,51) + rrt(056) * density(14) 
  pd(51,14) = pd(51,14) - rrt(056) * density(51) 
  pd(51,51) = pd(51,51) - rrt(056) * density(14) 
  pd(53,14) = pd(53,14) + rrt(056) * density(51) 
  pd(53,51) = pd(53,51) + rrt(056) * density(14) 
  pd(29,29) = pd(29,29) - rrt(057) * density(48) 
  pd(29,48) = pd(29,48) - rrt(057) * density(29) 
  pd(42,29) = pd(42,29) + rrt(057) * density(48) 
  pd(42,48) = pd(42,48) + rrt(057) * density(29) 
  pd(48,29) = pd(48,29) - rrt(057) * density(48) 
  pd(48,48) = pd(48,48) - rrt(057) * density(29) 
  pd(53,29) = pd(53,29) + rrt(057) * density(48) 
  pd(53,48) = pd(53,48) + rrt(057) * density(29) 
  pd(29,29) = pd(29,29) - rrt(058) * density(49) 
  pd(29,49) = pd(29,49) - rrt(058) * density(29) 
  pd(40,29) = pd(40,29) + rrt(058) * density(49) * 2.0d0
  pd(40,49) = pd(40,49) + rrt(058) * density(29) * 2.0d0
  pd(49,29) = pd(49,29) - rrt(058) * density(49) 
  pd(49,49) = pd(49,49) - rrt(058) * density(29) 
  pd(53,29) = pd(53,29) + rrt(058) * density(49) 
  pd(53,49) = pd(53,49) + rrt(058) * density(29) 
  pd(21,29) = pd(21,29) + rrt(059) * density(50) 
  pd(21,50) = pd(21,50) + rrt(059) * density(29) 
  pd(29,29) = pd(29,29) - rrt(059) * density(50) 
  pd(29,50) = pd(29,50) - rrt(059) * density(29) 
  pd(40,29) = pd(40,29) + rrt(059) * density(50) 
  pd(40,50) = pd(40,50) + rrt(059) * density(29) 
  pd(50,29) = pd(50,29) - rrt(059) * density(50) 
  pd(50,50) = pd(50,50) - rrt(059) * density(29) 
  pd(53,29) = pd(53,29) + rrt(059) * density(50) 
  pd(53,50) = pd(53,50) + rrt(059) * density(29) 
  pd(29,29) = pd(29,29) - rrt(060) * density(51) 
  pd(29,51) = pd(29,51) - rrt(060) * density(29) 
  pd(32,29) = pd(32,29) + rrt(060) * density(51) 
  pd(32,51) = pd(32,51) + rrt(060) * density(29) 
  pd(40,29) = pd(40,29) + rrt(060) * density(51) 
  pd(40,51) = pd(40,51) + rrt(060) * density(29) 
  pd(51,29) = pd(51,29) - rrt(060) * density(51) 
  pd(51,51) = pd(51,51) - rrt(060) * density(29) 
  pd(53,29) = pd(53,29) + rrt(060) * density(51) 
  pd(53,51) = pd(53,51) + rrt(060) * density(29) 
  pd(01,10) = pd(01,10) + rrt(061) * density(48) 
  pd(01,48) = pd(01,48) + rrt(061) * density(10) 
  pd(10,10) = pd(10,10) - rrt(061) * density(48) 
  pd(10,48) = pd(10,48) - rrt(061) * density(10) 
  pd(40,10) = pd(40,10) + rrt(061) * density(48) 
  pd(40,48) = pd(40,48) + rrt(061) * density(10) 
  pd(48,10) = pd(48,10) - rrt(061) * density(48) 
  pd(48,48) = pd(48,48) - rrt(061) * density(10) 
  pd(53,10) = pd(53,10) + rrt(061) * density(48) 
  pd(53,48) = pd(53,48) + rrt(061) * density(10) 
  pd(01,10) = pd(01,10) + rrt(062) * density(50) 
  pd(01,50) = pd(01,50) + rrt(062) * density(10) 
  pd(10,10) = pd(10,10) - rrt(062) * density(50) 
  pd(10,50) = pd(10,50) - rrt(062) * density(10) 
  pd(42,10) = pd(42,10) + rrt(062) * density(50) 
  pd(42,50) = pd(42,50) + rrt(062) * density(10) 
  pd(50,10) = pd(50,10) - rrt(062) * density(50) 
  pd(50,50) = pd(50,50) - rrt(062) * density(10) 
  pd(53,10) = pd(53,10) + rrt(062) * density(50) 
  pd(53,50) = pd(53,50) + rrt(062) * density(10) 
  pd(01,11) = pd(01,11) + rrt(063) * density(48) 
  pd(01,48) = pd(01,48) + rrt(063) * density(11) 
  pd(11,11) = pd(11,11) - rrt(063) * density(48) 
  pd(11,48) = pd(11,48) - rrt(063) * density(11) 
  pd(40,11) = pd(40,11) + rrt(063) * density(48) 
  pd(40,48) = pd(40,48) + rrt(063) * density(11) 
  pd(48,11) = pd(48,11) - rrt(063) * density(48) 
  pd(48,48) = pd(48,48) - rrt(063) * density(11) 
  pd(53,11) = pd(53,11) + rrt(063) * density(48) 
  pd(53,48) = pd(53,48) + rrt(063) * density(11) 
  pd(01,11) = pd(01,11) + rrt(064) * density(50) 
  pd(01,50) = pd(01,50) + rrt(064) * density(11) 
  pd(11,11) = pd(11,11) - rrt(064) * density(50) 
  pd(11,50) = pd(11,50) - rrt(064) * density(11) 
  pd(42,11) = pd(42,11) + rrt(064) * density(50) 
  pd(42,50) = pd(42,50) + rrt(064) * density(11) 
  pd(50,11) = pd(50,11) - rrt(064) * density(50) 
  pd(50,50) = pd(50,50) - rrt(064) * density(11) 
  pd(53,11) = pd(53,11) + rrt(064) * density(50) 
  pd(53,50) = pd(53,50) + rrt(064) * density(11) 
  pd(11,13) = pd(11,13) + rrt(065) 
  pd(13,13) = pd(13,13) - rrt(065) 
  pd(10,10) = pd(10,10) - rrt(066) * density(29) 
  pd(10,29) = pd(10,29) - rrt(066) * density(10) 
  pd(15,10) = pd(15,10) + rrt(066) * density(29) 
  pd(15,29) = pd(15,29) + rrt(066) * density(10) 
  pd(29,10) = pd(29,10) - rrt(066) * density(29) 
  pd(29,29) = pd(29,29) - rrt(066) * density(10) 
  pd(40,10) = pd(40,10) + rrt(066) * density(29) 
  pd(40,29) = pd(40,29) + rrt(066) * density(10) 
  pd(01,10) = pd(01,10) + rrt(067) * density(21) 
  pd(01,21) = pd(01,21) + rrt(067) * density(10) 
  pd(10,10) = pd(10,10) - rrt(067) * density(21) 
  pd(10,21) = pd(10,21) - rrt(067) * density(10) 
  pd(21,10) = pd(21,10) - rrt(067) * density(21) 
  pd(21,21) = pd(21,21) - rrt(067) * density(10) 
  pd(29,10) = pd(29,10) + rrt(067) * density(21) 
  pd(29,21) = pd(29,21) + rrt(067) * density(10) 
  pd(30,10) = pd(30,10) + rrt(067) * density(21) 
  pd(30,21) = pd(30,21) + rrt(067) * density(10) 
  pd(01,10) = pd(01,10) + rrt(068) * density(21) 
  pd(01,21) = pd(01,21) + rrt(068) * density(10) 
  pd(10,10) = pd(10,10) - rrt(068) * density(21) 
  pd(10,21) = pd(10,21) - rrt(068) * density(10) 
  pd(21,10) = pd(21,10) - rrt(068) * density(21) 
  pd(21,21) = pd(21,21) - rrt(068) * density(10) 
  pd(27,10) = pd(27,10) + rrt(068) * density(21) 
  pd(27,21) = pd(27,21) + rrt(068) * density(10) 
  pd(01,01) = pd(01,01) + rrt(069) * density(10) 
  pd(01,10) = pd(01,10) + rrt(069) * density(01) 
  pd(10,01) = pd(10,01) - rrt(069) * density(10) 
  pd(10,10) = pd(10,10) - rrt(069) * density(01) 
  pd(01,10) = pd(01,10) + rrt(070) * density(40) 
  pd(01,40) = pd(01,40) + rrt(070) * density(10) 
  pd(10,10) = pd(10,10) - rrt(070) * density(40) 
  pd(10,40) = pd(10,40) - rrt(070) * density(10) 
  pd(01,10) = pd(01,10) + rrt(071) * density(41) 
  pd(01,41) = pd(01,41) + rrt(071) * density(10) 
  pd(10,10) = pd(10,10) - rrt(071) * density(41) 
  pd(10,41) = pd(10,41) - rrt(071) * density(10) 
  pd(14,10) = pd(14,10) + rrt(071) * density(41) 
  pd(14,41) = pd(14,41) + rrt(071) * density(10) 
  pd(40,10) = pd(40,10) + rrt(071) * density(41) 
  pd(40,41) = pd(40,41) + rrt(071) * density(10) 
  pd(41,10) = pd(41,10) - rrt(071) * density(41) 
  pd(41,41) = pd(41,41) - rrt(071) * density(10) 
  pd(01,10) = pd(01,10) + rrt(072) * density(42) 
  pd(01,42) = pd(01,42) + rrt(072) * density(10) 
  pd(10,10) = pd(10,10) - rrt(072) * density(42) 
  pd(10,42) = pd(10,42) - rrt(072) * density(10) 
  pd(29,10) = pd(29,10) + rrt(072) * density(42) 
  pd(29,42) = pd(29,42) + rrt(072) * density(10) 
  pd(40,10) = pd(40,10) + rrt(072) * density(42) 
  pd(40,42) = pd(40,42) + rrt(072) * density(10) 
  pd(42,10) = pd(42,10) - rrt(072) * density(42) 
  pd(42,42) = pd(42,42) - rrt(072) * density(10) 
  pd(01,10) = pd(01,10) + rrt(073) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(073) * density(10) * 4.0d0
  pd(11,10) = pd(11,10) + rrt(073) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) + rrt(074) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(074) * density(10) * 4.0d0
  pd(13,10) = pd(13,10) + rrt(074) * density(10) * 2.0d0
  pd(10,01) = pd(10,01) + rrt(075) * density(11) 
  pd(10,11) = pd(10,11) + rrt(075) * density(01) 
  pd(11,01) = pd(11,01) - rrt(075) * density(11) 
  pd(11,11) = pd(11,11) - rrt(075) * density(01) 
  pd(01,01) = pd(01,01) + rrt(076) * density(11) 
  pd(01,11) = pd(01,11) + rrt(076) * density(01) 
  pd(11,01) = pd(11,01) - rrt(076) * density(11) 
  pd(11,11) = pd(11,11) - rrt(076) * density(01) 
  pd(01,11) = pd(01,11) + rrt(077) * density(21) 
  pd(01,21) = pd(01,21) + rrt(077) * density(11) 
  pd(11,11) = pd(11,11) - rrt(077) * density(21) 
  pd(11,21) = pd(11,21) - rrt(077) * density(11) 
  pd(21,11) = pd(21,11) - rrt(077) * density(21) 
  pd(21,21) = pd(21,21) - rrt(077) * density(11) 
  pd(29,11) = pd(29,11) + rrt(077) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(077) * density(11) * 2.0d0
  pd(10,11) = pd(10,11) + rrt(078) * density(40) 
  pd(10,40) = pd(10,40) + rrt(078) * density(11) 
  pd(11,11) = pd(11,11) - rrt(078) * density(40) 
  pd(11,40) = pd(11,40) - rrt(078) * density(11) 
  pd(12,01) = pd(12,01) + rrt(079) * density(13) 
  pd(12,13) = pd(12,13) + rrt(079) * density(01) 
  pd(13,01) = pd(13,01) - rrt(079) * density(13) 
  pd(13,13) = pd(13,13) - rrt(079) * density(01) 
  pd(11,01) = pd(11,01) + rrt(080) * density(12) 
  pd(11,12) = pd(11,12) + rrt(080) * density(01) 
  pd(12,01) = pd(12,01) - rrt(080) * density(12) 
  pd(12,12) = pd(12,12) - rrt(080) * density(01) 
  pd(01,12) = pd(01,12) + rrt(081) * density(21) 
  pd(01,21) = pd(01,21) + rrt(081) * density(12) 
  pd(12,12) = pd(12,12) - rrt(081) * density(21) 
  pd(12,21) = pd(12,21) - rrt(081) * density(12) 
  pd(21,12) = pd(21,12) - rrt(081) * density(21) 
  pd(21,21) = pd(21,21) - rrt(081) * density(12) 
  pd(29,12) = pd(29,12) + rrt(081) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(081) * density(12) * 2.0d0
  pd(01,12) = pd(01,12) + rrt(082) * density(40) 
  pd(01,40) = pd(01,40) + rrt(082) * density(12) 
  pd(12,12) = pd(12,12) - rrt(082) * density(40) 
  pd(12,40) = pd(12,40) - rrt(082) * density(12) 
  pd(14,12) = pd(14,12) + rrt(082) * density(40) 
  pd(14,40) = pd(14,40) + rrt(082) * density(12) 
  pd(29,12) = pd(29,12) + rrt(082) * density(40) 
  pd(29,40) = pd(29,40) + rrt(082) * density(12) 
  pd(40,12) = pd(40,12) - rrt(082) * density(40) 
  pd(40,40) = pd(40,40) - rrt(082) * density(12) 
  pd(15,15) = pd(15,15) - rrt(083) * density(21) 
  pd(15,21) = pd(15,21) - rrt(083) * density(15) 
  pd(21,15) = pd(21,15) - rrt(083) * density(21) 
  pd(21,21) = pd(21,21) - rrt(083) * density(15) 
  pd(29,15) = pd(29,15) + rrt(083) * density(21) 
  pd(29,21) = pd(29,21) + rrt(083) * density(15) 
  pd(40,15) = pd(40,15) + rrt(083) * density(21) 
  pd(40,21) = pd(40,21) + rrt(083) * density(15) 
  pd(01,15) = pd(01,15) + rrt(084) * density(40) 
  pd(01,40) = pd(01,40) + rrt(084) * density(15) 
  pd(15,15) = pd(15,15) - rrt(084) * density(40) 
  pd(15,40) = pd(15,40) - rrt(084) * density(15) 
  pd(29,15) = pd(29,15) + rrt(084) * density(40) 
  pd(29,40) = pd(29,40) + rrt(084) * density(15) 
  pd(40,15) = pd(40,15) - rrt(084) * density(40) 
  pd(40,40) = pd(40,40) - rrt(084) * density(15) 
  pd(01,15) = pd(01,15) + rrt(085) * density(41) 
  pd(01,41) = pd(01,41) + rrt(085) * density(15) 
  pd(15,15) = pd(15,15) - rrt(085) * density(41) 
  pd(15,41) = pd(15,41) - rrt(085) * density(15) 
  pd(40,15) = pd(40,15) + rrt(085) * density(41) 
  pd(40,41) = pd(40,41) + rrt(085) * density(15) 
  pd(41,15) = pd(41,15) - rrt(085) * density(41) 
  pd(41,41) = pd(41,41) - rrt(085) * density(15) 
  pd(16,16) = pd(16,16) - rrt(086) * density(21) 
  pd(16,21) = pd(16,21) - rrt(086) * density(16) 
  pd(21,16) = pd(21,16) - rrt(086) * density(21) 
  pd(21,21) = pd(21,21) - rrt(086) * density(16) 
  pd(29,16) = pd(29,16) + rrt(086) * density(21) 
  pd(29,21) = pd(29,21) + rrt(086) * density(16) 
  pd(40,16) = pd(40,16) + rrt(086) * density(21) 
  pd(40,21) = pd(40,21) + rrt(086) * density(16) 
  pd(10,16) = pd(10,16) + rrt(087) * density(40) 
  pd(10,40) = pd(10,40) + rrt(087) * density(16) 
  pd(16,16) = pd(16,16) - rrt(087) * density(40) 
  pd(16,40) = pd(16,40) - rrt(087) * density(16) 
  pd(29,16) = pd(29,16) + rrt(087) * density(40) 
  pd(29,40) = pd(29,40) + rrt(087) * density(16) 
  pd(40,16) = pd(40,16) - rrt(087) * density(40) 
  pd(40,40) = pd(40,40) - rrt(087) * density(16) 
  pd(14,14) = pd(14,14) - rrt(088) * density(26) 
  pd(14,26) = pd(14,26) - rrt(088) * density(14) 
  pd(26,14) = pd(26,14) - rrt(088) * density(26) 
  pd(26,26) = pd(26,26) - rrt(088) * density(14) 
  pd(29,14) = pd(29,14) + rrt(088) * density(26) 
  pd(29,26) = pd(29,26) + rrt(088) * density(14) 
  pd(40,14) = pd(40,14) + rrt(088) * density(26) 
  pd(40,26) = pd(40,26) + rrt(088) * density(14) 
  pd(29,01) = pd(29,01) + rrt(089) * density(30) 
  pd(29,30) = pd(29,30) + rrt(089) * density(01) 
  pd(30,01) = pd(30,01) - rrt(089) * density(30) 
  pd(30,30) = pd(30,30) - rrt(089) * density(01) 
  pd(14,30) = pd(14,30) + rrt(090) * density(40) 
  pd(14,40) = pd(14,40) + rrt(090) * density(30) 
  pd(21,30) = pd(21,30) + rrt(090) * density(40) 
  pd(21,40) = pd(21,40) + rrt(090) * density(30) 
  pd(30,30) = pd(30,30) - rrt(090) * density(40) 
  pd(30,40) = pd(30,40) - rrt(090) * density(30) 
  pd(40,30) = pd(40,30) - rrt(090) * density(40) 
  pd(40,40) = pd(40,40) - rrt(090) * density(30) 
  pd(30,30) = pd(30,30) - rrt(091) * density(41) 
  pd(30,41) = pd(30,41) - rrt(091) * density(30) 
  pd(40,30) = pd(40,30) + rrt(091) * density(41) * 2.0d0
  pd(40,41) = pd(40,41) + rrt(091) * density(30) * 2.0d0
  pd(41,30) = pd(41,30) - rrt(091) * density(41) 
  pd(41,41) = pd(41,41) - rrt(091) * density(30) 
  pd(01,14) = pd(01,14) + rrt(092) * density(40) 
  pd(01,40) = pd(01,40) + rrt(092) * density(14) 
  pd(14,14) = pd(14,14) - rrt(092) * density(40) 
  pd(14,40) = pd(14,40) - rrt(092) * density(14) 
  pd(29,14) = pd(29,14) + rrt(092) * density(40) 
  pd(29,40) = pd(29,40) + rrt(092) * density(14) 
  pd(40,14) = pd(40,14) - rrt(092) * density(40) 
  pd(40,40) = pd(40,40) - rrt(092) * density(14) 
  pd(14,14) = pd(14,14) - rrt(093) * density(21) 
  pd(14,21) = pd(14,21) - rrt(093) * density(14) 
  pd(21,14) = pd(21,14) - rrt(093) * density(21) 
  pd(21,21) = pd(21,21) - rrt(093) * density(14) 
  pd(29,14) = pd(29,14) + rrt(093) * density(21) 
  pd(29,21) = pd(29,21) + rrt(093) * density(14) 
  pd(40,14) = pd(40,14) + rrt(093) * density(21) 
  pd(40,21) = pd(40,21) + rrt(093) * density(14) 
  pd(01,14) = pd(01,14) + rrt(094) * density(42) 
  pd(01,42) = pd(01,42) + rrt(094) * density(14) 
  pd(14,14) = pd(14,14) - rrt(094) * density(42) 
  pd(14,42) = pd(14,42) - rrt(094) * density(14) 
  pd(29,14) = pd(29,14) + rrt(094) * density(42) * 2.0d0
  pd(29,42) = pd(29,42) + rrt(094) * density(14) * 2.0d0
  pd(42,14) = pd(42,14) - rrt(094) * density(42) 
  pd(42,42) = pd(42,42) - rrt(094) * density(14) 
  pd(14,14) = pd(14,14) - rrt(095) * density(42) 
  pd(14,42) = pd(14,42) - rrt(095) * density(14) 
  pd(29,14) = pd(29,14) + rrt(095) * density(42) 
  pd(29,42) = pd(29,42) + rrt(095) * density(14) 
  pd(41,14) = pd(41,14) + rrt(095) * density(42) 
  pd(41,42) = pd(41,42) + rrt(095) * density(14) 
  pd(42,14) = pd(42,14) - rrt(095) * density(42) 
  pd(42,42) = pd(42,42) - rrt(095) * density(14) 
  pd(01,14) = pd(01,14) + rrt(096) * density(42) 
  pd(01,42) = pd(01,42) + rrt(096) * density(14) 
  pd(14,14) = pd(14,14) - rrt(096) * density(42) 
  pd(14,42) = pd(14,42) - rrt(096) * density(14) 
  pd(21,14) = pd(21,14) + rrt(096) * density(42) 
  pd(21,42) = pd(21,42) + rrt(096) * density(14) 
  pd(42,14) = pd(42,14) - rrt(096) * density(42) 
  pd(42,42) = pd(42,42) - rrt(096) * density(14) 
  pd(14,14) = pd(14,14) - rrt(097) * density(42) 
  pd(14,42) = pd(14,42) - rrt(097) * density(14) 
  pd(40,14) = pd(40,14) + rrt(097) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(097) * density(14) * 2.0d0
  pd(42,14) = pd(42,14) - rrt(097) * density(42) 
  pd(42,42) = pd(42,42) - rrt(097) * density(14) 
  pd(14,29) = pd(14,29) + rrt(098) * density(40) 
  pd(14,40) = pd(14,40) + rrt(098) * density(29) 
  pd(21,29) = pd(21,29) + rrt(098) * density(40) 
  pd(21,40) = pd(21,40) + rrt(098) * density(29) 
  pd(29,29) = pd(29,29) - rrt(098) * density(40) 
  pd(29,40) = pd(29,40) - rrt(098) * density(29) 
  pd(40,29) = pd(40,29) - rrt(098) * density(40) 
  pd(40,40) = pd(40,40) - rrt(098) * density(29) 
  pd(29,29) = pd(29,29) - rrt(099) * density(40) 
  pd(29,40) = pd(29,40) - rrt(099) * density(29) 
  pd(40,29) = pd(40,29) - rrt(099) * density(40) 
  pd(40,40) = pd(40,40) - rrt(099) * density(29) 
  pd(42,29) = pd(42,29) + rrt(099) * density(40) 
  pd(42,40) = pd(42,40) + rrt(099) * density(29) 
  pd(29,29) = pd(29,29) - rrt(100) * density(41) 
  pd(29,41) = pd(29,41) - rrt(100) * density(29) 
  pd(40,29) = pd(40,29) + rrt(100) * density(41) * 2.0d0
  pd(40,41) = pd(40,41) + rrt(100) * density(29) * 2.0d0
  pd(41,29) = pd(41,29) - rrt(100) * density(41) 
  pd(41,41) = pd(41,41) - rrt(100) * density(29) 
  pd(21,29) = pd(21,29) + rrt(101) * density(42) 
  pd(21,42) = pd(21,42) + rrt(101) * density(29) 
  pd(29,29) = pd(29,29) - rrt(101) * density(42) 
  pd(29,42) = pd(29,42) - rrt(101) * density(29) 
  pd(40,29) = pd(40,29) + rrt(101) * density(42) 
  pd(40,42) = pd(40,42) + rrt(101) * density(29) 
  pd(42,29) = pd(42,29) - rrt(101) * density(42) 
  pd(42,42) = pd(42,42) - rrt(101) * density(29) 
  pd(21,29) = pd(21,29) + rrt(102) * density(43) 
  pd(21,43) = pd(21,43) + rrt(102) * density(29) 
  pd(29,29) = pd(29,29) - rrt(102) * density(43) 
  pd(29,43) = pd(29,43) - rrt(102) * density(29) 
  pd(42,29) = pd(42,29) + rrt(102) * density(43) 
  pd(42,43) = pd(42,43) + rrt(102) * density(29) 
  pd(43,29) = pd(43,29) - rrt(102) * density(43) 
  pd(43,43) = pd(43,43) - rrt(102) * density(29) 
  pd(21,21) = pd(21,21) - rrt(103) * density(40) 
  pd(21,40) = pd(21,40) - rrt(103) * density(21) 
  pd(29,21) = pd(29,21) + rrt(103) * density(40) 
  pd(29,40) = pd(29,40) + rrt(103) * density(21) 
  pd(40,21) = pd(40,21) - rrt(103) * density(40) 
  pd(40,40) = pd(40,40) - rrt(103) * density(21) 
  pd(42,21) = pd(42,21) + rrt(103) * density(40) 
  pd(42,40) = pd(42,40) + rrt(103) * density(21) 
  pd(21,32) = pd(21,32) + rrt(104) * density(40) 
  pd(21,40) = pd(21,40) + rrt(104) * density(32) 
  pd(32,32) = pd(32,32) - rrt(104) * density(40) 
  pd(32,40) = pd(32,40) - rrt(104) * density(32) 
  pd(40,32) = pd(40,32) - rrt(104) * density(40) 
  pd(40,40) = pd(40,40) - rrt(104) * density(32) 
  pd(42,32) = pd(42,32) + rrt(104) * density(40) 
  pd(42,40) = pd(42,40) + rrt(104) * density(32) 
  pd(01,40) = pd(01,40) + rrt(105) * density(41) 
  pd(01,41) = pd(01,41) + rrt(105) * density(40) 
  pd(40,40) = pd(40,40) - rrt(105) * density(41) 
  pd(40,41) = pd(40,41) - rrt(105) * density(40) 
  pd(41,40) = pd(41,40) - rrt(105) * density(41) 
  pd(41,41) = pd(41,41) - rrt(105) * density(40) 
  pd(42,40) = pd(42,40) + rrt(105) * density(41) 
  pd(42,41) = pd(42,41) + rrt(105) * density(40) 
  pd(40,40) = pd(40,40) - rrt(106) * density(43) 
  pd(40,43) = pd(40,43) - rrt(106) * density(40) 
  pd(42,40) = pd(42,40) + rrt(106) * density(43) * 2.0d0
  pd(42,43) = pd(42,43) + rrt(106) * density(40) * 2.0d0
  pd(43,40) = pd(43,40) - rrt(106) * density(43) 
  pd(43,43) = pd(43,43) - rrt(106) * density(40) 
  pd(21,21) = pd(21,21) - rrt(107) * density(42) 
  pd(21,42) = pd(21,42) - rrt(107) * density(21) 
  pd(32,21) = pd(32,21) + rrt(107) * density(42) 
  pd(32,42) = pd(32,42) + rrt(107) * density(21) 
  pd(40,21) = pd(40,21) + rrt(107) * density(42) 
  pd(40,42) = pd(40,42) + rrt(107) * density(21) 
  pd(42,21) = pd(42,21) - rrt(107) * density(42) 
  pd(42,42) = pd(42,42) - rrt(107) * density(21) 
  pd(21,42) = pd(21,42) + rrt(108) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(108) * density(42) * 4.0d0
  pd(42,42) = pd(42,42) - rrt(108) * density(42) * 4.0d0
  pd(40,42) = pd(40,42) + rrt(109) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(109) * density(42) * 4.0d0
  pd(43,42) = pd(43,42) + rrt(109) * density(42) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(110) * density(42) 
  pd(21,42) = pd(21,42) + rrt(110) * density(32) 
  pd(32,32) = pd(32,32) - rrt(110) * density(42) 
  pd(32,42) = pd(32,42) - rrt(110) * density(32) 
  pd(42,32) = pd(42,32) - rrt(110) * density(42) 
  pd(42,42) = pd(42,42) - rrt(110) * density(32) 
  pd(43,32) = pd(43,32) + rrt(110) * density(42) 
  pd(43,42) = pd(43,42) + rrt(110) * density(32) 
  pd(21,42) = pd(21,42) + rrt(111) * density(43) 
  pd(21,43) = pd(21,43) + rrt(111) * density(42) 
  pd(40,42) = pd(40,42) + rrt(111) * density(43) 
  pd(40,43) = pd(40,43) + rrt(111) * density(42) 
  pd(43,42) = pd(43,42) - rrt(111) * density(43) 
  pd(43,43) = pd(43,43) - rrt(111) * density(42) 
  pd(21,21) = pd(21,21) - rrt(112) * density(43) 
  pd(21,43) = pd(21,43) - rrt(112) * density(21) 
  pd(32,21) = pd(32,21) + rrt(112) * density(43) 
  pd(32,43) = pd(32,43) + rrt(112) * density(21) 
  pd(42,21) = pd(42,21) + rrt(112) * density(43) 
  pd(42,43) = pd(42,43) + rrt(112) * density(21) 
  pd(43,21) = pd(43,21) - rrt(112) * density(43) 
  pd(43,43) = pd(43,43) - rrt(112) * density(21) 
  pd(21,43) = pd(21,43) + rrt(113) * density(43) * 2.0d0
  pd(42,43) = pd(42,43) + rrt(113) * density(43) * 4.0d0
  pd(43,43) = pd(43,43) - rrt(113) * density(43) * 4.0d0
  pd(29,01) = pd(29,01) + rrt(114) * density(42) 
  pd(29,42) = pd(29,42) + rrt(114) * density(01) 
  pd(40,01) = pd(40,01) + rrt(114) * density(42) 
  pd(40,42) = pd(40,42) + rrt(114) * density(01) 
  pd(42,01) = pd(42,01) - rrt(114) * density(42) 
  pd(42,42) = pd(42,42) - rrt(114) * density(01) 
  pd(29,01) = pd(29,01) + rrt(115) * density(43) 
  pd(29,43) = pd(29,43) + rrt(115) * density(01) 
  pd(42,01) = pd(42,01) + rrt(115) * density(43) 
  pd(42,43) = pd(42,43) + rrt(115) * density(01) 
  pd(43,01) = pd(43,01) - rrt(115) * density(43) 
  pd(43,43) = pd(43,43) - rrt(115) * density(01) 
  pd(29,21) = pd(29,21) + rrt(116) * density(43) 
  pd(29,43) = pd(29,43) + rrt(116) * density(21) 
  pd(42,21) = pd(42,21) + rrt(116) * density(43) 
  pd(42,43) = pd(42,43) + rrt(116) * density(21) 
  pd(43,21) = pd(43,21) - rrt(116) * density(43) 
  pd(43,43) = pd(43,43) - rrt(116) * density(21) 
  pd(29,40) = pd(29,40) + rrt(117) * density(43) 
  pd(29,43) = pd(29,43) + rrt(117) * density(40) 
  pd(42,40) = pd(42,40) + rrt(117) * density(43) 
  pd(42,43) = pd(42,43) + rrt(117) * density(40) 
  pd(43,40) = pd(43,40) - rrt(117) * density(43) 
  pd(43,43) = pd(43,43) - rrt(117) * density(40) 
  pd(29,14) = pd(29,14) + rrt(118) * density(43) 
  pd(29,43) = pd(29,43) + rrt(118) * density(14) 
  pd(42,14) = pd(42,14) + rrt(118) * density(43) 
  pd(42,43) = pd(42,43) + rrt(118) * density(14) 
  pd(43,14) = pd(43,14) - rrt(118) * density(43) 
  pd(43,43) = pd(43,43) - rrt(118) * density(14) 
  pd(29,29) = pd(29,29) + rrt(119) * density(43) 
  pd(29,43) = pd(29,43) + rrt(119) * density(29) 
  pd(42,29) = pd(42,29) + rrt(119) * density(43) 
  pd(42,43) = pd(42,43) + rrt(119) * density(29) 
  pd(43,29) = pd(43,29) - rrt(119) * density(43) 
  pd(43,43) = pd(43,43) - rrt(119) * density(29) 
  pd(21,01) = pd(21,01) + rrt(120) * density(43) 
  pd(21,43) = pd(21,43) + rrt(120) * density(01) 
  pd(40,01) = pd(40,01) + rrt(120) * density(43) 
  pd(40,43) = pd(40,43) + rrt(120) * density(01) 
  pd(43,01) = pd(43,01) - rrt(120) * density(43) 
  pd(43,43) = pd(43,43) - rrt(120) * density(01) 
  pd(21,21) = pd(21,21) + rrt(121) * density(43) 
  pd(21,43) = pd(21,43) + rrt(121) * density(21) 
  pd(40,21) = pd(40,21) + rrt(121) * density(43) 
  pd(40,43) = pd(40,43) + rrt(121) * density(21) 
  pd(43,21) = pd(43,21) - rrt(121) * density(43) 
  pd(43,43) = pd(43,43) - rrt(121) * density(21) 
  pd(21,40) = pd(21,40) + rrt(122) * density(43) 
  pd(21,43) = pd(21,43) + rrt(122) * density(40) 
  pd(40,40) = pd(40,40) + rrt(122) * density(43) 
  pd(40,43) = pd(40,43) + rrt(122) * density(40) 
  pd(43,40) = pd(43,40) - rrt(122) * density(43) 
  pd(43,43) = pd(43,43) - rrt(122) * density(40) 
  pd(21,14) = pd(21,14) + rrt(123) * density(43) 
  pd(21,43) = pd(21,43) + rrt(123) * density(14) 
  pd(40,14) = pd(40,14) + rrt(123) * density(43) 
  pd(40,43) = pd(40,43) + rrt(123) * density(14) 
  pd(43,14) = pd(43,14) - rrt(123) * density(43) 
  pd(43,43) = pd(43,43) - rrt(123) * density(14) 
  pd(21,29) = pd(21,29) + rrt(124) * density(43) 
  pd(21,43) = pd(21,43) + rrt(124) * density(29) 
  pd(40,29) = pd(40,29) + rrt(124) * density(43) 
  pd(40,43) = pd(40,43) + rrt(124) * density(29) 
  pd(43,29) = pd(43,29) - rrt(124) * density(43) 
  pd(43,43) = pd(43,43) - rrt(124) * density(29) 
  pd(42,44) = pd(42,44) + rrt(125) 
  pd(43,44) = pd(43,44) + rrt(125) 
  pd(44,44) = pd(44,44) - rrt(125) 
  pd(14,01) = pd(14,01) - rrt(126) * density(14) * density(29) 
  pd(14,14) = pd(14,14) - rrt(126) * density(01) * density(29) 
  pd(14,29) = pd(14,29) - rrt(126) * density(01) * density(14) 
  pd(29,01) = pd(29,01) - rrt(126) * density(14) * density(29) 
  pd(29,14) = pd(29,14) - rrt(126) * density(01) * density(29) 
  pd(29,29) = pd(29,29) - rrt(126) * density(01) * density(14) 
  pd(40,01) = pd(40,01) + rrt(126) * density(14) * density(29) 
  pd(40,14) = pd(40,14) + rrt(126) * density(01) * density(29) 
  pd(40,29) = pd(40,29) + rrt(126) * density(01) * density(14) 
  pd(14,14) = pd(14,14) - rrt(127) * density(21) * density(29) 
  pd(14,21) = pd(14,21) - rrt(127) * density(14) * density(29) 
  pd(14,29) = pd(14,29) - rrt(127) * density(14) * density(21) 
  pd(29,14) = pd(29,14) - rrt(127) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(127) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(127) * density(14) * density(21) 
  pd(40,14) = pd(40,14) + rrt(127) * density(21) * density(29) 
  pd(40,21) = pd(40,21) + rrt(127) * density(14) * density(29) 
  pd(40,29) = pd(40,29) + rrt(127) * density(14) * density(21) 
  pd(14,14) = pd(14,14) - rrt(128) * density(14) * density(29) * 2.0d0
  pd(14,29) = pd(14,29) - rrt(128) * density(14)**2 
  pd(29,14) = pd(29,14) - rrt(128) * density(14) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(128) * density(14)**2 
  pd(40,14) = pd(40,14) + rrt(128) * density(14) * density(29) * 2.0d0
  pd(40,29) = pd(40,29) + rrt(128) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(129) * density(29)**2 
  pd(14,29) = pd(14,29) - rrt(129) * density(14) * density(29) * 2.0d0
  pd(29,14) = pd(29,14) - rrt(129) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(129) * density(14) * density(29) * 2.0d0
  pd(40,14) = pd(40,14) + rrt(129) * density(29)**2 
  pd(40,29) = pd(40,29) + rrt(129) * density(14) * density(29) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(130) * density(29) * density(40) 
  pd(14,29) = pd(14,29) - rrt(130) * density(14) * density(40) 
  pd(14,40) = pd(14,40) - rrt(130) * density(14) * density(29) 
  pd(29,14) = pd(29,14) - rrt(130) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(130) * density(14) * density(40) 
  pd(29,40) = pd(29,40) - rrt(130) * density(14) * density(29) 
  pd(40,14) = pd(40,14) + rrt(130) * density(29) * density(40) 
  pd(40,29) = pd(40,29) + rrt(130) * density(14) * density(40) 
  pd(40,40) = pd(40,40) + rrt(130) * density(14) * density(29) 
  pd(21,01) = pd(21,01) - rrt(131) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(131) * density(01) * density(29) 
  pd(21,29) = pd(21,29) - rrt(131) * density(01) * density(21) 
  pd(29,01) = pd(29,01) - rrt(131) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(131) * density(01) * density(29) 
  pd(29,29) = pd(29,29) - rrt(131) * density(01) * density(21) 
  pd(32,01) = pd(32,01) + rrt(131) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(131) * density(01) * density(29) 
  pd(32,29) = pd(32,29) + rrt(131) * density(01) * density(21) 
  pd(29,01) = pd(29,01) - rrt(132) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(132) * density(01) * density(40) 
  pd(29,40) = pd(29,40) - rrt(132) * density(01) * density(29) 
  pd(40,01) = pd(40,01) - rrt(132) * density(29) * density(40) 
  pd(40,29) = pd(40,29) - rrt(132) * density(01) * density(40) 
  pd(40,40) = pd(40,40) - rrt(132) * density(01) * density(29) 
  pd(42,01) = pd(42,01) + rrt(132) * density(29) * density(40) 
  pd(42,29) = pd(42,29) + rrt(132) * density(01) * density(40) 
  pd(42,40) = pd(42,40) + rrt(132) * density(01) * density(29) 
  pd(29,21) = pd(29,21) - rrt(133) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(133) * density(21) * density(40) 
  pd(29,40) = pd(29,40) - rrt(133) * density(21) * density(29) 
  pd(40,21) = pd(40,21) - rrt(133) * density(29) * density(40) 
  pd(40,29) = pd(40,29) - rrt(133) * density(21) * density(40) 
  pd(40,40) = pd(40,40) - rrt(133) * density(21) * density(29) 
  pd(42,21) = pd(42,21) + rrt(133) * density(29) * density(40) 
  pd(42,29) = pd(42,29) + rrt(133) * density(21) * density(40) 
  pd(42,40) = pd(42,40) + rrt(133) * density(21) * density(29) 
  pd(29,29) = pd(29,29) - rrt(134) * density(40)**2 
  pd(29,40) = pd(29,40) - rrt(134) * density(29) * density(40) * 2.0d0
  pd(40,29) = pd(40,29) - rrt(134) * density(40)**2 
  pd(40,40) = pd(40,40) - rrt(134) * density(29) * density(40) * 2.0d0
  pd(42,29) = pd(42,29) + rrt(134) * density(40)**2 
  pd(42,40) = pd(42,40) + rrt(134) * density(29) * density(40) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(135) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(135) * density(01) * density(42) 
  pd(29,42) = pd(29,42) - rrt(135) * density(01) * density(29) 
  pd(42,01) = pd(42,01) - rrt(135) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(135) * density(01) * density(42) 
  pd(42,42) = pd(42,42) - rrt(135) * density(01) * density(29) 
  pd(43,01) = pd(43,01) + rrt(135) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(135) * density(01) * density(42) 
  pd(43,42) = pd(43,42) + rrt(135) * density(01) * density(29) 
  pd(29,21) = pd(29,21) - rrt(136) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(136) * density(21) * density(42) 
  pd(29,42) = pd(29,42) - rrt(136) * density(21) * density(29) 
  pd(42,21) = pd(42,21) - rrt(136) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(136) * density(21) * density(42) 
  pd(42,42) = pd(42,42) - rrt(136) * density(21) * density(29) 
  pd(43,21) = pd(43,21) + rrt(136) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(136) * density(21) * density(42) 
  pd(43,42) = pd(43,42) + rrt(136) * density(21) * density(29) 
  pd(29,14) = pd(29,14) - rrt(137) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(137) * density(14) * density(42) 
  pd(29,42) = pd(29,42) - rrt(137) * density(14) * density(29) 
  pd(42,14) = pd(42,14) - rrt(137) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(137) * density(14) * density(42) 
  pd(42,42) = pd(42,42) - rrt(137) * density(14) * density(29) 
  pd(43,14) = pd(43,14) + rrt(137) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(137) * density(14) * density(42) 
  pd(43,42) = pd(43,42) + rrt(137) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(138) * density(29) * density(42) * 2.0d0
  pd(29,42) = pd(29,42) - rrt(138) * density(29)**2 
  pd(42,29) = pd(42,29) - rrt(138) * density(29) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(138) * density(29)**2 
  pd(43,29) = pd(43,29) + rrt(138) * density(29) * density(42) * 2.0d0
  pd(43,42) = pd(43,42) + rrt(138) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(139) * density(40) * density(42) 
  pd(29,40) = pd(29,40) - rrt(139) * density(29) * density(42) 
  pd(29,42) = pd(29,42) - rrt(139) * density(29) * density(40) 
  pd(42,29) = pd(42,29) - rrt(139) * density(40) * density(42) 
  pd(42,40) = pd(42,40) - rrt(139) * density(29) * density(42) 
  pd(42,42) = pd(42,42) - rrt(139) * density(29) * density(40) 
  pd(43,29) = pd(43,29) + rrt(139) * density(40) * density(42) 
  pd(43,40) = pd(43,40) + rrt(139) * density(29) * density(42) 
  pd(43,42) = pd(43,42) + rrt(139) * density(29) * density(40) 
  pd(42,42) = pd(42,42) - rrt(140) * density(43) 
  pd(42,43) = pd(42,43) - rrt(140) * density(42) 
  pd(43,42) = pd(43,42) - rrt(140) * density(43) 
  pd(43,43) = pd(43,43) - rrt(140) * density(42) 
  pd(44,42) = pd(44,42) + rrt(140) * density(43) 
  pd(44,43) = pd(44,43) + rrt(140) * density(42) 
  pd(17,17) = pd(17,17) - rrt(141) * density(21) 
  pd(17,21) = pd(17,21) - rrt(141) * density(17) 
  pd(21,17) = pd(21,17) - rrt(141) * density(21) 
  pd(21,21) = pd(21,21) - rrt(141) * density(17) 
  pd(33,17) = pd(33,17) + rrt(141) * density(21) 
  pd(33,21) = pd(33,21) + rrt(141) * density(17) 
  pd(40,17) = pd(40,17) + rrt(141) * density(21) 
  pd(40,21) = pd(40,21) + rrt(141) * density(17) 
  pd(14,17) = pd(14,17) + rrt(142) * density(40) 
  pd(14,40) = pd(14,40) + rrt(142) * density(17) 
  pd(17,17) = pd(17,17) - rrt(142) * density(40) 
  pd(17,40) = pd(17,40) - rrt(142) * density(17) 
  pd(40,17) = pd(40,17) - rrt(142) * density(40) 
  pd(40,40) = pd(40,40) - rrt(142) * density(17) 
  pd(45,17) = pd(45,17) + rrt(142) * density(40) 
  pd(45,40) = pd(45,40) + rrt(142) * density(17) 
  pd(17,17) = pd(17,17) - rrt(143) * density(40) 
  pd(17,40) = pd(17,40) - rrt(143) * density(17) 
  pd(18,17) = pd(18,17) + rrt(143) * density(40) 
  pd(18,40) = pd(18,40) + rrt(143) * density(17) 
  pd(29,17) = pd(29,17) + rrt(143) * density(40) 
  pd(29,40) = pd(29,40) + rrt(143) * density(17) 
  pd(40,17) = pd(40,17) - rrt(143) * density(40) 
  pd(40,40) = pd(40,40) - rrt(143) * density(17) 
  pd(01,17) = pd(01,17) + rrt(144) * density(40) 
  pd(01,40) = pd(01,40) + rrt(144) * density(17) 
  pd(17,17) = pd(17,17) - rrt(144) * density(40) 
  pd(17,40) = pd(17,40) - rrt(144) * density(17) 
  pd(33,17) = pd(33,17) + rrt(144) * density(40) 
  pd(33,40) = pd(33,40) + rrt(144) * density(17) 
  pd(40,17) = pd(40,17) - rrt(144) * density(40) 
  pd(40,40) = pd(40,40) - rrt(144) * density(17) 
  pd(29,33) = pd(29,33) + rrt(145) * density(40) 
  pd(29,40) = pd(29,40) + rrt(145) * density(33) 
  pd(33,33) = pd(33,33) - rrt(145) * density(40) 
  pd(33,40) = pd(33,40) - rrt(145) * density(33) 
  pd(40,33) = pd(40,33) - rrt(145) * density(40) 
  pd(40,40) = pd(40,40) - rrt(145) * density(33) 
  pd(45,33) = pd(45,33) + rrt(145) * density(40) 
  pd(45,40) = pd(45,40) + rrt(145) * density(33) 
  pd(14,33) = pd(14,33) + rrt(146) * density(40) 
  pd(14,40) = pd(14,40) + rrt(146) * density(33) 
  pd(33,33) = pd(33,33) - rrt(146) * density(40) 
  pd(33,40) = pd(33,40) - rrt(146) * density(33) 
  pd(34,33) = pd(34,33) + rrt(146) * density(40) 
  pd(34,40) = pd(34,40) + rrt(146) * density(33) 
  pd(40,33) = pd(40,33) - rrt(146) * density(40) 
  pd(40,40) = pd(40,40) - rrt(146) * density(33) 
  pd(33,33) = pd(33,33) - rrt(147) * density(41) 
  pd(33,41) = pd(33,41) - rrt(147) * density(33) 
  pd(40,33) = pd(40,33) + rrt(147) * density(41) 
  pd(40,41) = pd(40,41) + rrt(147) * density(33) 
  pd(41,33) = pd(41,33) - rrt(147) * density(41) 
  pd(41,41) = pd(41,41) - rrt(147) * density(33) 
  pd(45,33) = pd(45,33) + rrt(147) * density(41) 
  pd(45,41) = pd(45,41) + rrt(147) * density(33) 
  pd(29,33) = pd(29,33) + rrt(148) * density(42) 
  pd(29,42) = pd(29,42) + rrt(148) * density(33) 
  pd(33,33) = pd(33,33) - rrt(148) * density(42) 
  pd(33,42) = pd(33,42) - rrt(148) * density(33) 
  pd(42,33) = pd(42,33) - rrt(148) * density(42) 
  pd(42,42) = pd(42,42) - rrt(148) * density(33) 
  pd(47,33) = pd(47,33) + rrt(148) * density(42) 
  pd(47,42) = pd(47,42) + rrt(148) * density(33) 
  pd(01,18) = pd(01,18) + rrt(149) * density(40) 
  pd(01,40) = pd(01,40) + rrt(149) * density(18) 
  pd(18,18) = pd(18,18) - rrt(149) * density(40) 
  pd(18,40) = pd(18,40) - rrt(149) * density(18) 
  pd(40,18) = pd(40,18) - rrt(149) * density(40) 
  pd(40,40) = pd(40,40) - rrt(149) * density(18) 
  pd(45,18) = pd(45,18) + rrt(149) * density(40) 
  pd(45,40) = pd(45,40) + rrt(149) * density(18) 
  pd(01,01) = pd(01,01) - rrt(150) * density(34) 
  pd(01,34) = pd(01,34) - rrt(150) * density(01) 
  pd(34,01) = pd(34,01) - rrt(150) * density(34) 
  pd(34,34) = pd(34,34) - rrt(150) * density(01) 
  pd(40,01) = pd(40,01) + rrt(150) * density(34) 
  pd(40,34) = pd(40,34) + rrt(150) * density(01) 
  pd(45,01) = pd(45,01) + rrt(150) * density(34) 
  pd(45,34) = pd(45,34) + rrt(150) * density(01) 
  pd(21,34) = pd(21,34) + rrt(151) * density(40) 
  pd(21,40) = pd(21,40) + rrt(151) * density(34) 
  pd(34,34) = pd(34,34) - rrt(151) * density(40) 
  pd(34,40) = pd(34,40) - rrt(151) * density(34) 
  pd(40,34) = pd(40,34) - rrt(151) * density(40) 
  pd(40,40) = pd(40,40) - rrt(151) * density(34) 
  pd(45,34) = pd(45,34) + rrt(151) * density(40) 
  pd(45,40) = pd(45,40) + rrt(151) * density(34) 
  pd(32,34) = pd(32,34) + rrt(152) * density(42) 
  pd(32,42) = pd(32,42) + rrt(152) * density(34) 
  pd(34,34) = pd(34,34) - rrt(152) * density(42) 
  pd(34,42) = pd(34,42) - rrt(152) * density(34) 
  pd(42,34) = pd(42,34) - rrt(152) * density(42) 
  pd(42,42) = pd(42,42) - rrt(152) * density(34) 
  pd(45,34) = pd(45,34) + rrt(152) * density(42) 
  pd(45,42) = pd(45,42) + rrt(152) * density(34) 
  pd(21,34) = pd(21,34) + rrt(153) * density(42) 
  pd(21,42) = pd(21,42) + rrt(153) * density(34) 
  pd(34,34) = pd(34,34) - rrt(153) * density(42) 
  pd(34,42) = pd(34,42) - rrt(153) * density(34) 
  pd(42,34) = pd(42,34) - rrt(153) * density(42) 
  pd(42,42) = pd(42,42) - rrt(153) * density(34) 
  pd(47,34) = pd(47,34) + rrt(153) * density(42) 
  pd(47,42) = pd(47,42) + rrt(153) * density(34) 
  pd(01,19) = pd(01,19) + rrt(154) * density(40) 
  pd(01,40) = pd(01,40) + rrt(154) * density(19) 
  pd(14,19) = pd(14,19) + rrt(154) * density(40) 
  pd(14,40) = pd(14,40) + rrt(154) * density(19) 
  pd(19,19) = pd(19,19) - rrt(154) * density(40) 
  pd(19,40) = pd(19,40) - rrt(154) * density(19) 
  pd(40,19) = pd(40,19) - rrt(154) * density(40) 
  pd(40,40) = pd(40,40) - rrt(154) * density(19) 
  pd(45,19) = pd(45,19) + rrt(154) * density(40) 
  pd(45,40) = pd(45,40) + rrt(154) * density(19) 
  pd(01,19) = pd(01,19) + rrt(155) * density(40) 
  pd(01,40) = pd(01,40) + rrt(155) * density(19) 
  pd(19,19) = pd(19,19) - rrt(155) * density(40) 
  pd(19,40) = pd(19,40) - rrt(155) * density(19) 
  pd(40,19) = pd(40,19) - rrt(155) * density(40) 
  pd(40,40) = pd(40,40) - rrt(155) * density(19) 
  pd(46,19) = pd(46,19) + rrt(155) * density(40) 
  pd(46,40) = pd(46,40) + rrt(155) * density(19) 
  pd(40,40) = pd(40,40) - rrt(156) * density(47) 
  pd(40,47) = pd(40,47) - rrt(156) * density(40) 
  pd(42,40) = pd(42,40) + rrt(156) * density(47) 
  pd(42,47) = pd(42,47) + rrt(156) * density(40) 
  pd(45,40) = pd(45,40) + rrt(156) * density(47) 
  pd(45,47) = pd(45,47) + rrt(156) * density(40) 
  pd(47,40) = pd(47,40) - rrt(156) * density(47) 
  pd(47,47) = pd(47,47) - rrt(156) * density(40) 
  pd(40,40) = pd(40,40) - rrt(157) * density(46) 
  pd(40,46) = pd(40,46) - rrt(157) * density(40) 
  pd(41,40) = pd(41,40) + rrt(157) * density(46) 
  pd(41,46) = pd(41,46) + rrt(157) * density(40) 
  pd(45,40) = pd(45,40) + rrt(157) * density(46) 
  pd(45,46) = pd(45,46) + rrt(157) * density(40) 
  pd(46,40) = pd(46,40) - rrt(157) * density(46) 
  pd(46,46) = pd(46,46) - rrt(157) * density(40) 
  pd(01,20) = pd(01,20) + rrt(158) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(158) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(158) * density(40) 
  pd(20,40) = pd(20,40) - rrt(158) * density(20) 
  pd(40,20) = pd(40,20) - rrt(158) * density(40) 
  pd(40,40) = pd(40,40) - rrt(158) * density(20) 
  pd(45,20) = pd(45,20) + rrt(158) * density(40) 
  pd(45,40) = pd(45,40) + rrt(158) * density(20) 
  pd(21,35) = pd(21,35) + rrt(159) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(159) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(159) * density(40) 
  pd(35,40) = pd(35,40) - rrt(159) * density(35) 
  pd(40,35) = pd(40,35) - rrt(159) * density(40) 
  pd(40,40) = pd(40,40) - rrt(159) * density(35) 
  pd(45,35) = pd(45,35) + rrt(159) * density(40) 
  pd(45,40) = pd(45,40) + rrt(159) * density(35) 
  pd(29,36) = pd(29,36) + rrt(160) * density(42) 
  pd(29,42) = pd(29,42) + rrt(160) * density(36) 
  pd(36,36) = pd(36,36) - rrt(160) * density(42) 
  pd(36,42) = pd(36,42) - rrt(160) * density(36) 
  pd(42,36) = pd(42,36) - rrt(160) * density(42) 
  pd(42,42) = pd(42,42) - rrt(160) * density(36) 
  pd(50,36) = pd(50,36) + rrt(160) * density(42) 
  pd(50,42) = pd(50,42) + rrt(160) * density(36) 
  pd(36,36) = pd(36,36) - rrt(161) * density(41) 
  pd(36,41) = pd(36,41) - rrt(161) * density(36) 
  pd(40,36) = pd(40,36) + rrt(161) * density(41) 
  pd(40,41) = pd(40,41) + rrt(161) * density(36) 
  pd(41,36) = pd(41,36) - rrt(161) * density(41) 
  pd(41,41) = pd(41,41) - rrt(161) * density(36) 
  pd(48,36) = pd(48,36) + rrt(161) * density(41) 
  pd(48,41) = pd(48,41) + rrt(161) * density(36) 
  pd(21,37) = pd(21,37) + rrt(162) * density(42) 
  pd(21,42) = pd(21,42) + rrt(162) * density(37) 
  pd(37,37) = pd(37,37) - rrt(162) * density(42) 
  pd(37,42) = pd(37,42) - rrt(162) * density(37) 
  pd(42,37) = pd(42,37) - rrt(162) * density(42) 
  pd(42,42) = pd(42,42) - rrt(162) * density(37) 
  pd(50,37) = pd(50,37) + rrt(162) * density(42) 
  pd(50,42) = pd(50,42) + rrt(162) * density(37) 
  pd(29,38) = pd(29,38) + rrt(163) * density(40) 
  pd(29,40) = pd(29,40) + rrt(163) * density(38) 
  pd(38,38) = pd(38,38) - rrt(163) * density(40) 
  pd(38,40) = pd(38,40) - rrt(163) * density(38) 
  pd(40,38) = pd(40,38) - rrt(163) * density(40) 
  pd(40,40) = pd(40,40) - rrt(163) * density(38) 
  pd(51,38) = pd(51,38) + rrt(163) * density(40) 
  pd(51,40) = pd(51,40) + rrt(163) * density(38) 
  pd(21,38) = pd(21,38) + rrt(164) * density(40) 
  pd(21,40) = pd(21,40) + rrt(164) * density(38) 
  pd(38,38) = pd(38,38) - rrt(164) * density(40) 
  pd(38,40) = pd(38,40) - rrt(164) * density(38) 
  pd(40,38) = pd(40,38) - rrt(164) * density(40) 
  pd(40,40) = pd(40,40) - rrt(164) * density(38) 
  pd(50,38) = pd(50,38) + rrt(164) * density(40) 
  pd(50,40) = pd(50,40) + rrt(164) * density(38) 
  pd(32,38) = pd(32,38) + rrt(165) * density(42) 
  pd(32,42) = pd(32,42) + rrt(165) * density(38) 
  pd(38,38) = pd(38,38) - rrt(165) * density(42) 
  pd(38,42) = pd(38,42) - rrt(165) * density(38) 
  pd(42,38) = pd(42,38) - rrt(165) * density(42) 
  pd(42,42) = pd(42,42) - rrt(165) * density(38) 
  pd(50,38) = pd(50,38) + rrt(165) * density(42) 
  pd(50,42) = pd(50,42) + rrt(165) * density(38) 
  pd(21,38) = pd(21,38) + rrt(166) * density(42) 
  pd(21,42) = pd(21,42) + rrt(166) * density(38) 
  pd(38,38) = pd(38,38) - rrt(166) * density(42) 
  pd(38,42) = pd(38,42) - rrt(166) * density(38) 
  pd(42,38) = pd(42,38) - rrt(166) * density(42) 
  pd(42,42) = pd(42,42) - rrt(166) * density(38) 
  pd(51,38) = pd(51,38) + rrt(166) * density(42) 
  pd(51,42) = pd(51,42) + rrt(166) * density(38) 
  pd(21,21) = pd(21,21) - rrt(167) * density(48) 
  pd(21,48) = pd(21,48) - rrt(167) * density(21) 
  pd(37,21) = pd(37,21) + rrt(167) * density(48) 
  pd(37,48) = pd(37,48) + rrt(167) * density(21) 
  pd(40,21) = pd(40,21) + rrt(167) * density(48) 
  pd(40,48) = pd(40,48) + rrt(167) * density(21) 
  pd(48,21) = pd(48,21) - rrt(167) * density(48) 
  pd(48,48) = pd(48,48) - rrt(167) * density(21) 
  pd(40,42) = pd(40,42) + rrt(168) * density(48) 
  pd(40,48) = pd(40,48) + rrt(168) * density(42) 
  pd(42,42) = pd(42,42) - rrt(168) * density(48) 
  pd(42,48) = pd(42,48) - rrt(168) * density(42) 
  pd(48,42) = pd(48,42) - rrt(168) * density(48) 
  pd(48,48) = pd(48,48) - rrt(168) * density(42) 
  pd(50,42) = pd(50,42) + rrt(168) * density(48) 
  pd(50,48) = pd(50,48) + rrt(168) * density(42) 
  pd(40,42) = pd(40,42) + rrt(169) * density(50) 
  pd(40,50) = pd(40,50) + rrt(169) * density(42) 
  pd(42,42) = pd(42,42) - rrt(169) * density(50) 
  pd(42,50) = pd(42,50) - rrt(169) * density(42) 
  pd(50,42) = pd(50,42) - rrt(169) * density(50) 
  pd(50,50) = pd(50,50) - rrt(169) * density(42) 
  pd(51,42) = pd(51,42) + rrt(169) * density(50) 
  pd(51,50) = pd(51,50) + rrt(169) * density(42) 
  pd(42,43) = pd(42,43) + rrt(170) * density(50) 
  pd(42,50) = pd(42,50) + rrt(170) * density(43) 
  pd(43,43) = pd(43,43) - rrt(170) * density(50) 
  pd(43,50) = pd(43,50) - rrt(170) * density(43) 
  pd(50,43) = pd(50,43) - rrt(170) * density(50) 
  pd(50,50) = pd(50,50) - rrt(170) * density(43) 
  pd(51,43) = pd(51,43) + rrt(170) * density(50) 
  pd(51,50) = pd(51,50) + rrt(170) * density(43) 
  pd(42,44) = pd(42,44) + rrt(171) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(171) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(171) * density(50) 
  pd(44,50) = pd(44,50) - rrt(171) * density(44) 
  pd(50,44) = pd(50,44) - rrt(171) * density(50) 
  pd(50,50) = pd(50,50) - rrt(171) * density(44) 
  pd(51,44) = pd(51,44) + rrt(171) * density(50) 
  pd(51,50) = pd(51,50) + rrt(171) * density(44) 
  pd(40,40) = pd(40,40) - rrt(172) * density(51) 
  pd(40,51) = pd(40,51) - rrt(172) * density(40) 
  pd(42,40) = pd(42,40) + rrt(172) * density(51) 
  pd(42,51) = pd(42,51) + rrt(172) * density(40) 
  pd(50,40) = pd(50,40) + rrt(172) * density(51) 
  pd(50,51) = pd(50,51) + rrt(172) * density(40) 
  pd(51,40) = pd(51,40) - rrt(172) * density(51) 
  pd(51,51) = pd(51,51) - rrt(172) * density(40) 
  pd(21,39) = pd(21,39) + rrt(173) * density(40) 
  pd(21,40) = pd(21,40) + rrt(173) * density(39) 
  pd(39,39) = pd(39,39) - rrt(173) * density(40) 
  pd(39,40) = pd(39,40) - rrt(173) * density(39) 
  pd(40,39) = pd(40,39) - rrt(173) * density(40) 
  pd(40,40) = pd(40,40) - rrt(173) * density(39) 
  pd(51,39) = pd(51,39) + rrt(173) * density(40) 
  pd(51,40) = pd(51,40) + rrt(173) * density(39) 
  pd(36,36) = pd(36,36) - rrt(174) * density(40) 
  pd(36,40) = pd(36,40) - rrt(174) * density(36) 
  pd(40,36) = pd(40,36) - rrt(174) * density(40) 
  pd(40,40) = pd(40,40) - rrt(174) * density(36) 
  pd(50,36) = pd(50,36) + rrt(174) * density(40) 
  pd(50,40) = pd(50,40) + rrt(174) * density(36) 
  pd(29,36) = pd(29,36) + rrt(175) * density(45) 
  pd(29,45) = pd(29,45) + rrt(175) * density(36) 
  pd(36,36) = pd(36,36) - rrt(175) * density(45) 
  pd(36,45) = pd(36,45) - rrt(175) * density(36) 
  pd(40,36) = pd(40,36) + rrt(175) * density(45) 
  pd(40,45) = pd(40,45) + rrt(175) * density(36) 
  pd(45,36) = pd(45,36) - rrt(175) * density(45) 
  pd(45,45) = pd(45,45) - rrt(175) * density(36) 
  pd(29,36) = pd(29,36) + rrt(176) * density(47) 
  pd(29,47) = pd(29,47) + rrt(176) * density(36) 
  pd(36,36) = pd(36,36) - rrt(176) * density(47) 
  pd(36,47) = pd(36,47) - rrt(176) * density(36) 
  pd(42,36) = pd(42,36) + rrt(176) * density(47) 
  pd(42,47) = pd(42,47) + rrt(176) * density(36) 
  pd(47,36) = pd(47,36) - rrt(176) * density(47) 
  pd(47,47) = pd(47,47) - rrt(176) * density(36) 
  pd(21,37) = pd(21,37) + rrt(177) * density(45) 
  pd(21,45) = pd(21,45) + rrt(177) * density(37) 
  pd(37,37) = pd(37,37) - rrt(177) * density(45) 
  pd(37,45) = pd(37,45) - rrt(177) * density(37) 
  pd(40,37) = pd(40,37) + rrt(177) * density(45) 
  pd(40,45) = pd(40,45) + rrt(177) * density(37) 
  pd(45,37) = pd(45,37) - rrt(177) * density(45) 
  pd(45,45) = pd(45,45) - rrt(177) * density(37) 
  pd(21,37) = pd(21,37) + rrt(178) * density(47) 
  pd(21,47) = pd(21,47) + rrt(178) * density(37) 
  pd(37,37) = pd(37,37) - rrt(178) * density(47) 
  pd(37,47) = pd(37,47) - rrt(178) * density(37) 
  pd(42,37) = pd(42,37) + rrt(178) * density(47) 
  pd(42,47) = pd(42,47) + rrt(178) * density(37) 
  pd(47,37) = pd(47,37) - rrt(178) * density(47) 
  pd(47,47) = pd(47,47) - rrt(178) * density(37) 
  pd(32,38) = pd(32,38) + rrt(179) * density(45) 
  pd(32,45) = pd(32,45) + rrt(179) * density(38) 
  pd(38,38) = pd(38,38) - rrt(179) * density(45) 
  pd(38,45) = pd(38,45) - rrt(179) * density(38) 
  pd(40,38) = pd(40,38) + rrt(179) * density(45) 
  pd(40,45) = pd(40,45) + rrt(179) * density(38) 
  pd(45,38) = pd(45,38) - rrt(179) * density(45) 
  pd(45,45) = pd(45,45) - rrt(179) * density(38) 
  pd(32,38) = pd(32,38) + rrt(180) * density(47) 
  pd(32,47) = pd(32,47) + rrt(180) * density(38) 
  pd(38,38) = pd(38,38) - rrt(180) * density(47) 
  pd(38,47) = pd(38,47) - rrt(180) * density(38) 
  pd(42,38) = pd(42,38) + rrt(180) * density(47) 
  pd(42,47) = pd(42,47) + rrt(180) * density(38) 
  pd(47,38) = pd(47,38) - rrt(180) * density(47) 
  pd(47,47) = pd(47,47) - rrt(180) * density(38) 
  pd(14,17) = pd(14,17) + rrt(181) * density(48) 
  pd(14,48) = pd(14,48) + rrt(181) * density(17) 
  pd(17,17) = pd(17,17) - rrt(181) * density(48) 
  pd(17,48) = pd(17,48) - rrt(181) * density(17) 
  pd(40,17) = pd(40,17) + rrt(181) * density(48) 
  pd(40,48) = pd(40,48) + rrt(181) * density(17) 
  pd(48,17) = pd(48,17) - rrt(181) * density(48) 
  pd(48,48) = pd(48,48) - rrt(181) * density(17) 
  pd(01,18) = pd(01,18) + rrt(182) * density(48) 
  pd(01,48) = pd(01,48) + rrt(182) * density(18) 
  pd(18,18) = pd(18,18) - rrt(182) * density(48) 
  pd(18,48) = pd(18,48) - rrt(182) * density(18) 
  pd(40,18) = pd(40,18) + rrt(182) * density(48) 
  pd(40,48) = pd(40,48) + rrt(182) * density(18) 
  pd(48,18) = pd(48,18) - rrt(182) * density(48) 
  pd(48,48) = pd(48,48) - rrt(182) * density(18) 
  pd(29,33) = pd(29,33) + rrt(183) * density(48) 
  pd(29,48) = pd(29,48) + rrt(183) * density(33) 
  pd(33,33) = pd(33,33) - rrt(183) * density(48) 
  pd(33,48) = pd(33,48) - rrt(183) * density(33) 
  pd(40,33) = pd(40,33) + rrt(183) * density(48) 
  pd(40,48) = pd(40,48) + rrt(183) * density(33) 
  pd(48,33) = pd(48,33) - rrt(183) * density(48) 
  pd(48,48) = pd(48,48) - rrt(183) * density(33) 
  pd(21,34) = pd(21,34) + rrt(184) * density(48) 
  pd(21,48) = pd(21,48) + rrt(184) * density(34) 
  pd(34,34) = pd(34,34) - rrt(184) * density(48) 
  pd(34,48) = pd(34,48) - rrt(184) * density(34) 
  pd(40,34) = pd(40,34) + rrt(184) * density(48) 
  pd(40,48) = pd(40,48) + rrt(184) * density(34) 
  pd(48,34) = pd(48,34) - rrt(184) * density(48) 
  pd(48,48) = pd(48,48) - rrt(184) * density(34) 
  pd(40,45) = pd(40,45) + rrt(185) * density(48) * 2.0d0
  pd(40,48) = pd(40,48) + rrt(185) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(185) * density(48) 
  pd(45,48) = pd(45,48) - rrt(185) * density(45) 
  pd(48,45) = pd(48,45) - rrt(185) * density(48) 
  pd(48,48) = pd(48,48) - rrt(185) * density(45) 
  pd(40,46) = pd(40,46) + rrt(186) * density(48) 
  pd(40,48) = pd(40,48) + rrt(186) * density(46) 
  pd(41,46) = pd(41,46) + rrt(186) * density(48) 
  pd(41,48) = pd(41,48) + rrt(186) * density(46) 
  pd(46,46) = pd(46,46) - rrt(186) * density(48) 
  pd(46,48) = pd(46,48) - rrt(186) * density(46) 
  pd(48,46) = pd(48,46) - rrt(186) * density(48) 
  pd(48,48) = pd(48,48) - rrt(186) * density(46) 
  pd(40,47) = pd(40,47) + rrt(187) * density(48) 
  pd(40,48) = pd(40,48) + rrt(187) * density(47) 
  pd(42,47) = pd(42,47) + rrt(187) * density(48) 
  pd(42,48) = pd(42,48) + rrt(187) * density(47) 
  pd(47,47) = pd(47,47) - rrt(187) * density(48) 
  pd(47,48) = pd(47,48) - rrt(187) * density(47) 
  pd(48,47) = pd(48,47) - rrt(187) * density(48) 
  pd(48,48) = pd(48,48) - rrt(187) * density(47) 
  pd(40,45) = pd(40,45) + rrt(188) * density(49) 
  pd(40,49) = pd(40,49) + rrt(188) * density(45) 
  pd(41,45) = pd(41,45) + rrt(188) * density(49) 
  pd(41,49) = pd(41,49) + rrt(188) * density(45) 
  pd(45,45) = pd(45,45) - rrt(188) * density(49) 
  pd(45,49) = pd(45,49) - rrt(188) * density(45) 
  pd(49,45) = pd(49,45) - rrt(188) * density(49) 
  pd(49,49) = pd(49,49) - rrt(188) * density(45) 
  pd(41,47) = pd(41,47) + rrt(189) * density(49) 
  pd(41,49) = pd(41,49) + rrt(189) * density(47) 
  pd(42,47) = pd(42,47) + rrt(189) * density(49) 
  pd(42,49) = pd(42,49) + rrt(189) * density(47) 
  pd(47,47) = pd(47,47) - rrt(189) * density(49) 
  pd(47,49) = pd(47,49) - rrt(189) * density(47) 
  pd(49,47) = pd(49,47) - rrt(189) * density(49) 
  pd(49,49) = pd(49,49) - rrt(189) * density(47) 
  pd(14,17) = pd(14,17) + rrt(190) * density(50) 
  pd(14,50) = pd(14,50) + rrt(190) * density(17) 
  pd(17,17) = pd(17,17) - rrt(190) * density(50) 
  pd(17,50) = pd(17,50) - rrt(190) * density(17) 
  pd(42,17) = pd(42,17) + rrt(190) * density(50) 
  pd(42,50) = pd(42,50) + rrt(190) * density(17) 
  pd(50,17) = pd(50,17) - rrt(190) * density(50) 
  pd(50,50) = pd(50,50) - rrt(190) * density(17) 
  pd(01,18) = pd(01,18) + rrt(191) * density(50) 
  pd(01,50) = pd(01,50) + rrt(191) * density(18) 
  pd(18,18) = pd(18,18) - rrt(191) * density(50) 
  pd(18,50) = pd(18,50) - rrt(191) * density(18) 
  pd(42,18) = pd(42,18) + rrt(191) * density(50) 
  pd(42,50) = pd(42,50) + rrt(191) * density(18) 
  pd(50,18) = pd(50,18) - rrt(191) * density(50) 
  pd(50,50) = pd(50,50) - rrt(191) * density(18) 
  pd(29,33) = pd(29,33) + rrt(192) * density(50) 
  pd(29,50) = pd(29,50) + rrt(192) * density(33) 
  pd(33,33) = pd(33,33) - rrt(192) * density(50) 
  pd(33,50) = pd(33,50) - rrt(192) * density(33) 
  pd(42,33) = pd(42,33) + rrt(192) * density(50) 
  pd(42,50) = pd(42,50) + rrt(192) * density(33) 
  pd(50,33) = pd(50,33) - rrt(192) * density(50) 
  pd(50,50) = pd(50,50) - rrt(192) * density(33) 
  pd(21,34) = pd(21,34) + rrt(193) * density(50) 
  pd(21,50) = pd(21,50) + rrt(193) * density(34) 
  pd(34,34) = pd(34,34) - rrt(193) * density(50) 
  pd(34,50) = pd(34,50) - rrt(193) * density(34) 
  pd(42,34) = pd(42,34) + rrt(193) * density(50) 
  pd(42,50) = pd(42,50) + rrt(193) * density(34) 
  pd(50,34) = pd(50,34) - rrt(193) * density(50) 
  pd(50,50) = pd(50,50) - rrt(193) * density(34) 
  pd(40,45) = pd(40,45) + rrt(194) * density(50) 
  pd(40,50) = pd(40,50) + rrt(194) * density(45) 
  pd(42,45) = pd(42,45) + rrt(194) * density(50) 
  pd(42,50) = pd(42,50) + rrt(194) * density(45) 
  pd(45,45) = pd(45,45) - rrt(194) * density(50) 
  pd(45,50) = pd(45,50) - rrt(194) * density(45) 
  pd(50,45) = pd(50,45) - rrt(194) * density(50) 
  pd(50,50) = pd(50,50) - rrt(194) * density(45) 
  pd(41,46) = pd(41,46) + rrt(195) * density(50) 
  pd(41,50) = pd(41,50) + rrt(195) * density(46) 
  pd(42,46) = pd(42,46) + rrt(195) * density(50) 
  pd(42,50) = pd(42,50) + rrt(195) * density(46) 
  pd(46,46) = pd(46,46) - rrt(195) * density(50) 
  pd(46,50) = pd(46,50) - rrt(195) * density(46) 
  pd(50,46) = pd(50,46) - rrt(195) * density(50) 
  pd(50,50) = pd(50,50) - rrt(195) * density(46) 
  pd(42,47) = pd(42,47) + rrt(196) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(196) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(196) * density(50) 
  pd(47,50) = pd(47,50) - rrt(196) * density(47) 
  pd(50,47) = pd(50,47) - rrt(196) * density(50) 
  pd(50,50) = pd(50,50) - rrt(196) * density(47) 
  pd(01,18) = pd(01,18) + rrt(197) * density(51) 
  pd(01,51) = pd(01,51) + rrt(197) * density(18) 
  pd(18,18) = pd(18,18) - rrt(197) * density(51) 
  pd(18,51) = pd(18,51) - rrt(197) * density(18) 
  pd(43,18) = pd(43,18) + rrt(197) * density(51) 
  pd(43,51) = pd(43,51) + rrt(197) * density(18) 
  pd(51,18) = pd(51,18) - rrt(197) * density(51) 
  pd(51,51) = pd(51,51) - rrt(197) * density(18) 
  pd(40,45) = pd(40,45) + rrt(198) * density(51) 
  pd(40,51) = pd(40,51) + rrt(198) * density(45) 
  pd(43,45) = pd(43,45) + rrt(198) * density(51) 
  pd(43,51) = pd(43,51) + rrt(198) * density(45) 
  pd(45,45) = pd(45,45) - rrt(198) * density(51) 
  pd(45,51) = pd(45,51) - rrt(198) * density(45) 
  pd(51,45) = pd(51,45) - rrt(198) * density(51) 
  pd(51,51) = pd(51,51) - rrt(198) * density(45) 
  pd(42,47) = pd(42,47) + rrt(199) * density(51) 
  pd(42,51) = pd(42,51) + rrt(199) * density(47) 
  pd(43,47) = pd(43,47) + rrt(199) * density(51) 
  pd(43,51) = pd(43,51) + rrt(199) * density(47) 
  pd(47,47) = pd(47,47) - rrt(199) * density(51) 
  pd(47,51) = pd(47,51) - rrt(199) * density(47) 
  pd(51,47) = pd(51,47) - rrt(199) * density(51) 
  pd(51,51) = pd(51,51) - rrt(199) * density(47) 
  pd(01,37) = pd(01,37) + rrt(200) * density(52) 
  pd(01,52) = pd(01,52) + rrt(200) * density(37) 
  pd(21,37) = pd(21,37) + rrt(200) * density(52) * 2.0d0
  pd(21,52) = pd(21,52) + rrt(200) * density(37) * 2.0d0
  pd(37,37) = pd(37,37) - rrt(200) * density(52) 
  pd(37,52) = pd(37,52) - rrt(200) * density(37) 
  pd(52,37) = pd(52,37) - rrt(200) * density(52) 
  pd(52,52) = pd(52,52) - rrt(200) * density(37) 
  pd(14,18) = pd(14,18) + rrt(201) * density(48) * 2.0d0
  pd(14,48) = pd(14,48) + rrt(201) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(201) * density(48) 
  pd(18,48) = pd(18,48) - rrt(201) * density(18) 
  pd(40,18) = pd(40,18) + rrt(201) * density(48) 
  pd(40,48) = pd(40,48) + rrt(201) * density(18) 
  pd(48,18) = pd(48,18) - rrt(201) * density(48) 
  pd(48,48) = pd(48,48) - rrt(201) * density(18) 
  pd(01,19) = pd(01,19) + rrt(202) * density(48) 
  pd(01,48) = pd(01,48) + rrt(202) * density(19) 
  pd(14,19) = pd(14,19) + rrt(202) * density(48) 
  pd(14,48) = pd(14,48) + rrt(202) * density(19) 
  pd(19,19) = pd(19,19) - rrt(202) * density(48) 
  pd(19,48) = pd(19,48) - rrt(202) * density(19) 
  pd(40,19) = pd(40,19) + rrt(202) * density(48) 
  pd(40,48) = pd(40,48) + rrt(202) * density(19) 
  pd(48,19) = pd(48,19) - rrt(202) * density(48) 
  pd(48,48) = pd(48,48) - rrt(202) * density(19) 
  pd(01,20) = pd(01,20) + rrt(203) * density(48) * 2.0d0
  pd(01,48) = pd(01,48) + rrt(203) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(203) * density(48) 
  pd(20,48) = pd(20,48) - rrt(203) * density(20) 
  pd(40,20) = pd(40,20) + rrt(203) * density(48) 
  pd(40,48) = pd(40,48) + rrt(203) * density(20) 
  pd(48,20) = pd(48,20) - rrt(203) * density(48) 
  pd(48,48) = pd(48,48) - rrt(203) * density(20) 
  pd(29,34) = pd(29,34) + rrt(204) * density(48) * 2.0d0
  pd(29,48) = pd(29,48) + rrt(204) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(204) * density(48) 
  pd(34,48) = pd(34,48) - rrt(204) * density(34) 
  pd(40,34) = pd(40,34) + rrt(204) * density(48) 
  pd(40,48) = pd(40,48) + rrt(204) * density(34) 
  pd(48,34) = pd(48,34) - rrt(204) * density(48) 
  pd(48,48) = pd(48,48) - rrt(204) * density(34) 
  pd(21,35) = pd(21,35) + rrt(205) * density(48) * 2.0d0
  pd(21,48) = pd(21,48) + rrt(205) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(205) * density(48) 
  pd(35,48) = pd(35,48) - rrt(205) * density(35) 
  pd(40,35) = pd(40,35) + rrt(205) * density(48) 
  pd(40,48) = pd(40,48) + rrt(205) * density(35) 
  pd(48,35) = pd(48,35) - rrt(205) * density(48) 
  pd(48,48) = pd(48,48) - rrt(205) * density(35) 
  pd(14,45) = pd(14,45) + rrt(206) * density(48) 
  pd(14,48) = pd(14,48) + rrt(206) * density(45) 
  pd(29,45) = pd(29,45) + rrt(206) * density(48) 
  pd(29,48) = pd(29,48) + rrt(206) * density(45) 
  pd(40,45) = pd(40,45) + rrt(206) * density(48) 
  pd(40,48) = pd(40,48) + rrt(206) * density(45) 
  pd(45,45) = pd(45,45) - rrt(206) * density(48) 
  pd(45,48) = pd(45,48) - rrt(206) * density(45) 
  pd(48,45) = pd(48,45) - rrt(206) * density(48) 
  pd(48,48) = pd(48,48) - rrt(206) * density(45) 
  pd(01,46) = pd(01,46) + rrt(207) * density(48) 
  pd(01,48) = pd(01,48) + rrt(207) * density(46) 
  pd(29,46) = pd(29,46) + rrt(207) * density(48) 
  pd(29,48) = pd(29,48) + rrt(207) * density(46) 
  pd(40,46) = pd(40,46) + rrt(207) * density(48) 
  pd(40,48) = pd(40,48) + rrt(207) * density(46) 
  pd(46,46) = pd(46,46) - rrt(207) * density(48) 
  pd(46,48) = pd(46,48) - rrt(207) * density(46) 
  pd(48,46) = pd(48,46) - rrt(207) * density(48) 
  pd(48,48) = pd(48,48) - rrt(207) * density(46) 
  pd(14,47) = pd(14,47) + rrt(208) * density(48) 
  pd(14,48) = pd(14,48) + rrt(208) * density(47) 
  pd(21,47) = pd(21,47) + rrt(208) * density(48) 
  pd(21,48) = pd(21,48) + rrt(208) * density(47) 
  pd(40,47) = pd(40,47) + rrt(208) * density(48) 
  pd(40,48) = pd(40,48) + rrt(208) * density(47) 
  pd(47,47) = pd(47,47) - rrt(208) * density(48) 
  pd(47,48) = pd(47,48) - rrt(208) * density(47) 
  pd(48,47) = pd(48,47) - rrt(208) * density(48) 
  pd(48,48) = pd(48,48) - rrt(208) * density(47) 
  pd(01,48) = pd(01,48) + rrt(209) * density(52) 
  pd(01,52) = pd(01,52) + rrt(209) * density(48) 
  pd(21,48) = pd(21,48) + rrt(209) * density(52) 
  pd(21,52) = pd(21,52) + rrt(209) * density(48) 
  pd(40,48) = pd(40,48) + rrt(209) * density(52) 
  pd(40,52) = pd(40,52) + rrt(209) * density(48) 
  pd(48,48) = pd(48,48) - rrt(209) * density(52) 
  pd(48,52) = pd(48,52) - rrt(209) * density(48) 
  pd(52,48) = pd(52,48) - rrt(209) * density(52) 
  pd(52,52) = pd(52,52) - rrt(209) * density(48) 
  pd(14,18) = pd(14,18) + rrt(210) * density(49) * 2.0d0
  pd(14,49) = pd(14,49) + rrt(210) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(210) * density(49) 
  pd(18,49) = pd(18,49) - rrt(210) * density(18) 
  pd(41,18) = pd(41,18) + rrt(210) * density(49) 
  pd(41,49) = pd(41,49) + rrt(210) * density(18) 
  pd(49,18) = pd(49,18) - rrt(210) * density(49) 
  pd(49,49) = pd(49,49) - rrt(210) * density(18) 
  pd(14,18) = pd(14,18) + rrt(211) * density(50) * 2.0d0
  pd(14,50) = pd(14,50) + rrt(211) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(211) * density(50) 
  pd(18,50) = pd(18,50) - rrt(211) * density(18) 
  pd(42,18) = pd(42,18) + rrt(211) * density(50) 
  pd(42,50) = pd(42,50) + rrt(211) * density(18) 
  pd(50,18) = pd(50,18) - rrt(211) * density(50) 
  pd(50,50) = pd(50,50) - rrt(211) * density(18) 
  pd(01,19) = pd(01,19) + rrt(212) * density(50) 
  pd(01,50) = pd(01,50) + rrt(212) * density(19) 
  pd(14,19) = pd(14,19) + rrt(212) * density(50) 
  pd(14,50) = pd(14,50) + rrt(212) * density(19) 
  pd(19,19) = pd(19,19) - rrt(212) * density(50) 
  pd(19,50) = pd(19,50) - rrt(212) * density(19) 
  pd(42,19) = pd(42,19) + rrt(212) * density(50) 
  pd(42,50) = pd(42,50) + rrt(212) * density(19) 
  pd(50,19) = pd(50,19) - rrt(212) * density(50) 
  pd(50,50) = pd(50,50) - rrt(212) * density(19) 
  pd(01,20) = pd(01,20) + rrt(213) * density(50) * 2.0d0
  pd(01,50) = pd(01,50) + rrt(213) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(213) * density(50) 
  pd(20,50) = pd(20,50) - rrt(213) * density(20) 
  pd(42,20) = pd(42,20) + rrt(213) * density(50) 
  pd(42,50) = pd(42,50) + rrt(213) * density(20) 
  pd(50,20) = pd(50,20) - rrt(213) * density(50) 
  pd(50,50) = pd(50,50) - rrt(213) * density(20) 
  pd(29,34) = pd(29,34) + rrt(214) * density(50) * 2.0d0
  pd(29,50) = pd(29,50) + rrt(214) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(214) * density(50) 
  pd(34,50) = pd(34,50) - rrt(214) * density(34) 
  pd(42,34) = pd(42,34) + rrt(214) * density(50) 
  pd(42,50) = pd(42,50) + rrt(214) * density(34) 
  pd(50,34) = pd(50,34) - rrt(214) * density(50) 
  pd(50,50) = pd(50,50) - rrt(214) * density(34) 
  pd(21,35) = pd(21,35) + rrt(215) * density(50) * 2.0d0
  pd(21,50) = pd(21,50) + rrt(215) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(215) * density(50) 
  pd(35,50) = pd(35,50) - rrt(215) * density(35) 
  pd(42,35) = pd(42,35) + rrt(215) * density(50) 
  pd(42,50) = pd(42,50) + rrt(215) * density(35) 
  pd(50,35) = pd(50,35) - rrt(215) * density(50) 
  pd(50,50) = pd(50,50) - rrt(215) * density(35) 
  pd(14,45) = pd(14,45) + rrt(216) * density(50) 
  pd(14,50) = pd(14,50) + rrt(216) * density(45) 
  pd(29,45) = pd(29,45) + rrt(216) * density(50) 
  pd(29,50) = pd(29,50) + rrt(216) * density(45) 
  pd(42,45) = pd(42,45) + rrt(216) * density(50) 
  pd(42,50) = pd(42,50) + rrt(216) * density(45) 
  pd(45,45) = pd(45,45) - rrt(216) * density(50) 
  pd(45,50) = pd(45,50) - rrt(216) * density(45) 
  pd(50,45) = pd(50,45) - rrt(216) * density(50) 
  pd(50,50) = pd(50,50) - rrt(216) * density(45) 
  pd(01,46) = pd(01,46) + rrt(217) * density(50) 
  pd(01,50) = pd(01,50) + rrt(217) * density(46) 
  pd(29,46) = pd(29,46) + rrt(217) * density(50) 
  pd(29,50) = pd(29,50) + rrt(217) * density(46) 
  pd(42,46) = pd(42,46) + rrt(217) * density(50) 
  pd(42,50) = pd(42,50) + rrt(217) * density(46) 
  pd(46,46) = pd(46,46) - rrt(217) * density(50) 
  pd(46,50) = pd(46,50) - rrt(217) * density(46) 
  pd(50,46) = pd(50,46) - rrt(217) * density(50) 
  pd(50,50) = pd(50,50) - rrt(217) * density(46) 
  pd(14,47) = pd(14,47) + rrt(218) * density(50) 
  pd(14,50) = pd(14,50) + rrt(218) * density(47) 
  pd(21,47) = pd(21,47) + rrt(218) * density(50) 
  pd(21,50) = pd(21,50) + rrt(218) * density(47) 
  pd(42,47) = pd(42,47) + rrt(218) * density(50) 
  pd(42,50) = pd(42,50) + rrt(218) * density(47) 
  pd(47,47) = pd(47,47) - rrt(218) * density(50) 
  pd(47,50) = pd(47,50) - rrt(218) * density(47) 
  pd(50,47) = pd(50,47) - rrt(218) * density(50) 
  pd(50,50) = pd(50,50) - rrt(218) * density(47) 
  pd(01,50) = pd(01,50) + rrt(219) * density(52) 
  pd(01,52) = pd(01,52) + rrt(219) * density(50) 
  pd(21,50) = pd(21,50) + rrt(219) * density(52) 
  pd(21,52) = pd(21,52) + rrt(219) * density(50) 
  pd(42,50) = pd(42,50) + rrt(219) * density(52) 
  pd(42,52) = pd(42,52) + rrt(219) * density(50) 
  pd(50,50) = pd(50,50) - rrt(219) * density(52) 
  pd(50,52) = pd(50,52) - rrt(219) * density(50) 
  pd(52,50) = pd(52,50) - rrt(219) * density(52) 
  pd(52,52) = pd(52,52) - rrt(219) * density(50) 
  pd(21,39) = pd(21,39) + rrt(220) * density(45) * 2.0d0
  pd(21,45) = pd(21,45) + rrt(220) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(220) * density(45) 
  pd(39,45) = pd(39,45) - rrt(220) * density(39) 
  pd(40,39) = pd(40,39) + rrt(220) * density(45) 
  pd(40,45) = pd(40,45) + rrt(220) * density(39) 
  pd(45,39) = pd(45,39) - rrt(220) * density(45) 
  pd(45,45) = pd(45,45) - rrt(220) * density(39) 
  pd(21,39) = pd(21,39) + rrt(221) * density(47) * 2.0d0
  pd(21,47) = pd(21,47) + rrt(221) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(221) * density(47) 
  pd(39,47) = pd(39,47) - rrt(221) * density(39) 
  pd(42,39) = pd(42,39) + rrt(221) * density(47) 
  pd(42,47) = pd(42,47) + rrt(221) * density(39) 
  pd(47,39) = pd(47,39) - rrt(221) * density(47) 
  pd(47,47) = pd(47,47) - rrt(221) * density(39) 
  pd(29,36) = pd(29,36) + rrt(222) * density(45) 
  pd(29,45) = pd(29,45) + rrt(222) * density(36) 
  pd(36,36) = pd(36,36) - rrt(222) * density(45) 
  pd(36,45) = pd(36,45) - rrt(222) * density(36) 
  pd(40,36) = pd(40,36) + rrt(222) * density(45) 
  pd(40,45) = pd(40,45) + rrt(222) * density(36) 
  pd(45,36) = pd(45,36) - rrt(222) * density(45) 
  pd(45,45) = pd(45,45) - rrt(222) * density(36) 
  pd(21,37) = pd(21,37) + rrt(223) * density(45) 
  pd(21,45) = pd(21,45) + rrt(223) * density(37) 
  pd(37,37) = pd(37,37) - rrt(223) * density(45) 
  pd(37,45) = pd(37,45) - rrt(223) * density(37) 
  pd(40,37) = pd(40,37) + rrt(223) * density(45) 
  pd(40,45) = pd(40,45) + rrt(223) * density(37) 
  pd(45,37) = pd(45,37) - rrt(223) * density(45) 
  pd(45,45) = pd(45,45) - rrt(223) * density(37) 
  pd(17,17) = pd(17,17) - rrt(224) * density(36) 
  pd(17,36) = pd(17,36) - rrt(224) * density(17) 
  pd(36,17) = pd(36,17) - rrt(224) * density(36) 
  pd(36,36) = pd(36,36) - rrt(224) * density(17) 
  pd(40,17) = pd(40,17) + rrt(224) * density(36) 
  pd(40,36) = pd(40,36) + rrt(224) * density(17) 
  pd(36,36) = pd(36,36) - rrt(225) * density(45) 
  pd(36,45) = pd(36,45) - rrt(225) * density(36) 
  pd(42,36) = pd(42,36) + rrt(225) * density(45) 
  pd(42,45) = pd(42,45) + rrt(225) * density(36) 
  pd(45,36) = pd(45,36) - rrt(225) * density(45) 
  pd(45,45) = pd(45,45) - rrt(225) * density(36) 
  pd(17,17) = pd(17,17) - rrt(226) * density(37) 
  pd(17,37) = pd(17,37) - rrt(226) * density(17) 
  pd(37,17) = pd(37,17) - rrt(226) * density(37) 
  pd(37,37) = pd(37,37) - rrt(226) * density(17) 
  pd(42,17) = pd(42,17) + rrt(226) * density(37) 
  pd(42,37) = pd(42,37) + rrt(226) * density(17) 
  pd(32,38) = pd(32,38) + rrt(227) * density(45) 
  pd(32,45) = pd(32,45) + rrt(227) * density(38) 
  pd(38,38) = pd(38,38) - rrt(227) * density(45) 
  pd(38,45) = pd(38,45) - rrt(227) * density(38) 
  pd(40,38) = pd(40,38) + rrt(227) * density(45) 
  pd(40,45) = pd(40,45) + rrt(227) * density(38) 
  pd(45,38) = pd(45,38) - rrt(227) * density(45) 
  pd(45,45) = pd(45,45) - rrt(227) * density(38) 
  pd(32,38) = pd(32,38) + rrt(228) * density(47) 
  pd(32,47) = pd(32,47) + rrt(228) * density(38) 
  pd(38,38) = pd(38,38) - rrt(228) * density(47) 
  pd(38,47) = pd(38,47) - rrt(228) * density(38) 
  pd(42,38) = pd(42,38) + rrt(228) * density(47) 
  pd(42,47) = pd(42,47) + rrt(228) * density(38) 
  pd(47,38) = pd(47,38) - rrt(228) * density(47) 
  pd(47,47) = pd(47,47) - rrt(228) * density(38) 
  pd(14,17) = pd(14,17) + rrt(229) * density(48) 
  pd(14,48) = pd(14,48) + rrt(229) * density(17) 
  pd(17,17) = pd(17,17) - rrt(229) * density(48) 
  pd(17,48) = pd(17,48) - rrt(229) * density(17) 
  pd(40,17) = pd(40,17) + rrt(229) * density(48) 
  pd(40,48) = pd(40,48) + rrt(229) * density(17) 
  pd(48,17) = pd(48,17) - rrt(229) * density(48) 
  pd(48,48) = pd(48,48) - rrt(229) * density(17) 
  pd(01,18) = pd(01,18) + rrt(230) * density(48) 
  pd(01,48) = pd(01,48) + rrt(230) * density(18) 
  pd(18,18) = pd(18,18) - rrt(230) * density(48) 
  pd(18,48) = pd(18,48) - rrt(230) * density(18) 
  pd(40,18) = pd(40,18) + rrt(230) * density(48) 
  pd(40,48) = pd(40,48) + rrt(230) * density(18) 
  pd(48,18) = pd(48,18) - rrt(230) * density(48) 
  pd(48,48) = pd(48,48) - rrt(230) * density(18) 
  pd(29,33) = pd(29,33) + rrt(231) * density(48) 
  pd(29,48) = pd(29,48) + rrt(231) * density(33) 
  pd(33,33) = pd(33,33) - rrt(231) * density(48) 
  pd(33,48) = pd(33,48) - rrt(231) * density(33) 
  pd(40,33) = pd(40,33) + rrt(231) * density(48) 
  pd(40,48) = pd(40,48) + rrt(231) * density(33) 
  pd(48,33) = pd(48,33) - rrt(231) * density(48) 
  pd(48,48) = pd(48,48) - rrt(231) * density(33) 
  pd(21,34) = pd(21,34) + rrt(232) * density(48) 
  pd(21,48) = pd(21,48) + rrt(232) * density(34) 
  pd(34,34) = pd(34,34) - rrt(232) * density(48) 
  pd(34,48) = pd(34,48) - rrt(232) * density(34) 
  pd(40,34) = pd(40,34) + rrt(232) * density(48) 
  pd(40,48) = pd(40,48) + rrt(232) * density(34) 
  pd(48,34) = pd(48,34) - rrt(232) * density(48) 
  pd(48,48) = pd(48,48) - rrt(232) * density(34) 
  pd(40,45) = pd(40,45) + rrt(233) * density(48) * 2.0d0
  pd(40,48) = pd(40,48) + rrt(233) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(233) * density(48) 
  pd(45,48) = pd(45,48) - rrt(233) * density(45) 
  pd(48,45) = pd(48,45) - rrt(233) * density(48) 
  pd(48,48) = pd(48,48) - rrt(233) * density(45) 
  pd(40,46) = pd(40,46) + rrt(234) * density(48) 
  pd(40,48) = pd(40,48) + rrt(234) * density(46) 
  pd(41,46) = pd(41,46) + rrt(234) * density(48) 
  pd(41,48) = pd(41,48) + rrt(234) * density(46) 
  pd(46,46) = pd(46,46) - rrt(234) * density(48) 
  pd(46,48) = pd(46,48) - rrt(234) * density(46) 
  pd(48,46) = pd(48,46) - rrt(234) * density(48) 
  pd(48,48) = pd(48,48) - rrt(234) * density(46) 
  pd(40,47) = pd(40,47) + rrt(235) * density(48) 
  pd(40,48) = pd(40,48) + rrt(235) * density(47) 
  pd(42,47) = pd(42,47) + rrt(235) * density(48) 
  pd(42,48) = pd(42,48) + rrt(235) * density(47) 
  pd(47,47) = pd(47,47) - rrt(235) * density(48) 
  pd(47,48) = pd(47,48) - rrt(235) * density(47) 
  pd(48,47) = pd(48,47) - rrt(235) * density(48) 
  pd(48,48) = pd(48,48) - rrt(235) * density(47) 
  pd(40,45) = pd(40,45) + rrt(236) * density(49) 
  pd(40,49) = pd(40,49) + rrt(236) * density(45) 
  pd(41,45) = pd(41,45) + rrt(236) * density(49) 
  pd(41,49) = pd(41,49) + rrt(236) * density(45) 
  pd(45,45) = pd(45,45) - rrt(236) * density(49) 
  pd(45,49) = pd(45,49) - rrt(236) * density(45) 
  pd(49,45) = pd(49,45) - rrt(236) * density(49) 
  pd(49,49) = pd(49,49) - rrt(236) * density(45) 
  pd(41,47) = pd(41,47) + rrt(237) * density(49) 
  pd(41,49) = pd(41,49) + rrt(237) * density(47) 
  pd(42,47) = pd(42,47) + rrt(237) * density(49) 
  pd(42,49) = pd(42,49) + rrt(237) * density(47) 
  pd(47,47) = pd(47,47) - rrt(237) * density(49) 
  pd(47,49) = pd(47,49) - rrt(237) * density(47) 
  pd(49,47) = pd(49,47) - rrt(237) * density(49) 
  pd(49,49) = pd(49,49) - rrt(237) * density(47) 
  pd(14,17) = pd(14,17) + rrt(238) * density(50) 
  pd(14,50) = pd(14,50) + rrt(238) * density(17) 
  pd(17,17) = pd(17,17) - rrt(238) * density(50) 
  pd(17,50) = pd(17,50) - rrt(238) * density(17) 
  pd(42,17) = pd(42,17) + rrt(238) * density(50) 
  pd(42,50) = pd(42,50) + rrt(238) * density(17) 
  pd(50,17) = pd(50,17) - rrt(238) * density(50) 
  pd(50,50) = pd(50,50) - rrt(238) * density(17) 
  pd(01,18) = pd(01,18) + rrt(239) * density(50) 
  pd(01,50) = pd(01,50) + rrt(239) * density(18) 
  pd(18,18) = pd(18,18) - rrt(239) * density(50) 
  pd(18,50) = pd(18,50) - rrt(239) * density(18) 
  pd(42,18) = pd(42,18) + rrt(239) * density(50) 
  pd(42,50) = pd(42,50) + rrt(239) * density(18) 
  pd(50,18) = pd(50,18) - rrt(239) * density(50) 
  pd(50,50) = pd(50,50) - rrt(239) * density(18) 
  pd(29,33) = pd(29,33) + rrt(240) * density(50) 
  pd(29,50) = pd(29,50) + rrt(240) * density(33) 
  pd(33,33) = pd(33,33) - rrt(240) * density(50) 
  pd(33,50) = pd(33,50) - rrt(240) * density(33) 
  pd(42,33) = pd(42,33) + rrt(240) * density(50) 
  pd(42,50) = pd(42,50) + rrt(240) * density(33) 
  pd(50,33) = pd(50,33) - rrt(240) * density(50) 
  pd(50,50) = pd(50,50) - rrt(240) * density(33) 
  pd(21,34) = pd(21,34) + rrt(241) * density(50) 
  pd(21,50) = pd(21,50) + rrt(241) * density(34) 
  pd(34,34) = pd(34,34) - rrt(241) * density(50) 
  pd(34,50) = pd(34,50) - rrt(241) * density(34) 
  pd(42,34) = pd(42,34) + rrt(241) * density(50) 
  pd(42,50) = pd(42,50) + rrt(241) * density(34) 
  pd(50,34) = pd(50,34) - rrt(241) * density(50) 
  pd(50,50) = pd(50,50) - rrt(241) * density(34) 
  pd(40,45) = pd(40,45) + rrt(242) * density(50) 
  pd(40,50) = pd(40,50) + rrt(242) * density(45) 
  pd(42,45) = pd(42,45) + rrt(242) * density(50) 
  pd(42,50) = pd(42,50) + rrt(242) * density(45) 
  pd(45,45) = pd(45,45) - rrt(242) * density(50) 
  pd(45,50) = pd(45,50) - rrt(242) * density(45) 
  pd(50,45) = pd(50,45) - rrt(242) * density(50) 
  pd(50,50) = pd(50,50) - rrt(242) * density(45) 
  pd(41,46) = pd(41,46) + rrt(243) * density(50) 
  pd(41,50) = pd(41,50) + rrt(243) * density(46) 
  pd(42,46) = pd(42,46) + rrt(243) * density(50) 
  pd(42,50) = pd(42,50) + rrt(243) * density(46) 
  pd(46,46) = pd(46,46) - rrt(243) * density(50) 
  pd(46,50) = pd(46,50) - rrt(243) * density(46) 
  pd(50,46) = pd(50,46) - rrt(243) * density(50) 
  pd(50,50) = pd(50,50) - rrt(243) * density(46) 
  pd(42,47) = pd(42,47) + rrt(244) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(244) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(244) * density(50) 
  pd(47,50) = pd(47,50) - rrt(244) * density(47) 
  pd(50,47) = pd(50,47) - rrt(244) * density(50) 
  pd(50,50) = pd(50,50) - rrt(244) * density(47) 
  pd(40,45) = pd(40,45) + rrt(245) * density(51) 
  pd(40,51) = pd(40,51) + rrt(245) * density(45) 
  pd(43,45) = pd(43,45) + rrt(245) * density(51) 
  pd(43,51) = pd(43,51) + rrt(245) * density(45) 
  pd(45,45) = pd(45,45) - rrt(245) * density(51) 
  pd(45,51) = pd(45,51) - rrt(245) * density(45) 
  pd(51,45) = pd(51,45) - rrt(245) * density(51) 
  pd(51,51) = pd(51,51) - rrt(245) * density(45) 
  pd(42,47) = pd(42,47) + rrt(246) * density(51) 
  pd(42,51) = pd(42,51) + rrt(246) * density(47) 
  pd(43,47) = pd(43,47) + rrt(246) * density(51) 
  pd(43,51) = pd(43,51) + rrt(246) * density(47) 
  pd(47,47) = pd(47,47) - rrt(246) * density(51) 
  pd(47,51) = pd(47,51) - rrt(246) * density(47) 
  pd(51,47) = pd(51,47) - rrt(246) * density(51) 
  pd(51,51) = pd(51,51) - rrt(246) * density(47) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(54,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(54,:) = pd(54,:) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_jex
end module ZDPlasKin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction constant rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_rates(Time)
  use ZDPlasKin, only : ZDPlasKin_bolsig_rates, bolsig_rates, bolsig_pointer, ZDPlasKin_cfg, ZDPlasKin_get_density_total, &
                        lreaction_block, rrt
  implicit none
  double precision, intent(in) :: Time
  double precision :: Tgas
  double precision :: EN
  double precision :: Te
  double precision :: ANY_NEUTRAL
  DOUBLE PRECISION :: DTION, TIONN, TIONN2, TIONN3, TIONN4, TEFFN, TEFFN2, TEFFN3, TEFFN4 ! K
  DOUBLE PRECISION, PARAMETER :: ENERGY_VIBN2 = 0.290D0*11605.0D0, ENERGY_VIBO2 = 0.190D0*11605.0D0 ! K
  DOUBLE PRECISION :: QVIBN2, KVT10_N2N2, KVT01_N2N2, KVT10_N2N, KVT01_N2N, KVT10_N2O, KVT01_N2O ! CM3.S-1
  DOUBLE PRECISION :: QVIBO2, KVT10_O2O2, KVT01_O2O2, KVT10_O2O, KVT01_O2O ! CM3.S-1
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  EN  = ZDPlasKin_cfg(3)
  Te  = ZDPlasKin_cfg(4)
  call ZDPlasKin_get_density_total(ALL_NEUTRAL=ANY_NEUTRAL)
  DTION = 2.0D0 / ( 3.0D0 * 1.3807D-16 ) * 1.6605D-24 * ( 1.0D-17 * EN )**2
  TIONN = TGAS + DTION * 14.0D0 * 8.0D19**2
  TIONN2 = TGAS + DTION * 28.0D0 * 4.1D19**2
  TIONN3 = TGAS + DTION * 42.0D0 * 6.1D19**2
  TIONN4 = TGAS + DTION * 56.0D0 * 7.0D19**2
  TEFFN = ( TIONN + 0.5D0 * TGAS ) / ( 1.0D0 + 0.5D0 )
  TEFFN2 = ( TIONN2 + 1.0D0 * TGAS ) / ( 1.0D0 + 1.0D0 )
  TEFFN3 = ( TIONN3 + 1.5D0 * TGAS ) / ( 1.0D0 + 1.5D0 )
  TEFFN4 = ( TIONN4 + 2.0D0 * TGAS ) / ( 1.0D0 + 2.0D0 )
  QVIBN2 = EXP( - ENERGY_VIBN2 / TGAS )
  KVT10_N2N2 = 7.80D-12 * TGAS * EXP( - 218.0 / TGAS**(1.0/3.0) + 690.0 / TGAS ) / ( 1.0 - QVIBN2 )
  KVT10_N2N = 4.00D-16 * ( TGAS / 300.0D0 )**0.5
  KVT10_N2O = 1.20D-13 * EXP( - 27.6 / TGAS**(1.0/3.0) )
  KVT01_N2N2 = KVT10_N2N2 * QVIBN2
  KVT01_N2N = KVT10_N2N * QVIBN2
  KVT01_N2O = KVT10_N2O * QVIBN2
  QVIBO2 = EXP( - ENERGY_VIBO2 / TGAS )
  KVT10_O2O2 = 1.35D-12 * TGAS * EXP( - 137.9 / TGAS**(1.0/3.0) ) / ( 1.0 - QVIBO2 )
  KVT10_O2O = 4.50D-15 * TGAS
  KVT01_O2O2 = KVT10_O2O2 * QVIBO2
  KVT01_O2O = KVT10_O2O * QVIBO2
  rrt(001) = bolsig_rates(bolsig_pointer(1))
  rrt(002) = bolsig_rates(bolsig_pointer(2))
  rrt(003) = bolsig_rates(bolsig_pointer(3))
  rrt(004) = bolsig_rates(bolsig_pointer(4))
  rrt(005) = bolsig_rates(bolsig_pointer(5))
  rrt(006) = bolsig_rates(bolsig_pointer(6))
  rrt(007) = bolsig_rates(bolsig_pointer(7))
  rrt(008) = bolsig_rates(bolsig_pointer(8))
  rrt(009) = bolsig_rates(bolsig_pointer(9))
  rrt(010) = bolsig_rates(bolsig_pointer(10))
  rrt(011) = bolsig_rates(bolsig_pointer(11))
  rrt(012) = bolsig_rates(bolsig_pointer(12))
  rrt(013) = bolsig_rates(bolsig_pointer(13))
  rrt(014) = bolsig_rates(bolsig_pointer(14))
  rrt(015) = bolsig_rates(bolsig_pointer(15))
  rrt(016) = bolsig_rates(bolsig_pointer(16))
  rrt(017) = bolsig_rates(bolsig_pointer(17))
  rrt(018) = KVT10_N2N2*1.0D0
  rrt(019) = KVT10_N2N2*2.0D0
  rrt(020) = KVT10_N2N2*3.0D0
  rrt(021) = KVT10_N2N2*4.0D0
  rrt(022) = KVT10_N2N2*5.0D0
  rrt(023) = KVT10_N2N2*6.0D0
  rrt(024) = KVT10_N2N2*7.0D0
  rrt(025) = KVT10_N2N2*8.0D0
  rrt(026) = KVT10_N2O*1.0D0
  rrt(027) = KVT10_N2O*2.0D0
  rrt(028) = KVT10_N2O*3.0D0
  rrt(029) = KVT10_N2O*4.0D0
  rrt(030) = KVT10_N2O*5.0D0
  rrt(031) = KVT10_N2O*6.0D0
  rrt(032) = KVT10_N2O*7.0D0
  rrt(033) = KVT10_N2O*8.0D0
  rrt(034) = bolsig_rates(bolsig_pointer(18))
  rrt(035) = bolsig_rates(bolsig_pointer(19))
  rrt(036) = bolsig_rates(bolsig_pointer(20))
  rrt(037) = bolsig_rates(bolsig_pointer(21))
  rrt(038) = bolsig_rates(bolsig_pointer(22))
  rrt(039) = bolsig_rates(bolsig_pointer(23))
  rrt(040) = bolsig_rates(bolsig_pointer(24))
  rrt(041) = bolsig_rates(bolsig_pointer(25))
  rrt(042) = bolsig_rates(bolsig_pointer(26))
  rrt(043) = bolsig_rates(bolsig_pointer(27))
  rrt(044) = bolsig_rates(bolsig_pointer(28))
  rrt(045) = bolsig_rates(bolsig_pointer(29))
  rrt(046) = bolsig_rates(bolsig_pointer(30))
  rrt(047) = 2.0D-7*(300.0D0/TE)**0.5
  rrt(048) = 1.0D-11
  rrt(049) = 8.0D-31*ANY_NEUTRAL
  rrt(050) = 2.6D-10
  rrt(051) = 2.6D-10
  rrt(052) = 5.0D-10
  rrt(053) = 5.0D-10
  rrt(054) = 5.0D-10
  rrt(055) = 5.0D-10
  rrt(056) = 5.0D-10
  rrt(057) = 1.5D-10
  rrt(058) = 1.5D-10
  rrt(059) = 1.5D-10
  rrt(060) = 1.5D-10
  rrt(061) = 2.1D-9
  rrt(062) = 2.1D-9
  rrt(063) = 2.5D-9
  rrt(064) = 2.5D-9
  rrt(065) = 2.45D7
  rrt(066) = 7.0D-12
  rrt(067) = 2.1D-12*(TGAS/300.0D0)**0.55
  rrt(068) = 2.0D-13*(TGAS/300.0D0)**0.55
  rrt(069) = 3.0D-16
  rrt(070) = 6.9D-11
  rrt(071) = 1.0D-11
  rrt(072) = 1.0D-12
  rrt(073) = 3.0D-10
  rrt(074) = 1.5D-10
  rrt(075) = 3.0D-11
  rrt(076) = 2.0D-12
  rrt(077) = 3.0D-10
  rrt(078) = 2.4D-10
  rrt(079) = 1.0D-11
  rrt(080) = 1.9D-13
  rrt(081) = 2.8D-11
  rrt(082) = 3.6D-10
  rrt(083) = 5.2D-12
  rrt(084) = 1.8D-10
  rrt(085) = 3.5D-12
  rrt(086) = 2.6D-12
  rrt(087) = 3.0D-11
  rrt(088) = 2.0D-14*EXP(-600.0D0/TGAS)
  rrt(089) = 2.3D-11
  rrt(090) = 1.7D-10
  rrt(091) = 7.2D-11
  rrt(092) = 1.8D-11*(TGAS/300.0)**0.5
  rrt(093) = 3.2D-12*(TGAS/300.0)*EXP(-3150.0D0/TGAS)
  rrt(094) = 9.1D-13
  rrt(095) = 3.0D-12
  rrt(096) = 7.0D-13
  rrt(097) = 2.3D-12
  rrt(098) = 7.5D-12*(TGAS/300.0)*EXP(-19500.0D0/TGAS)
  rrt(099) = 4.2D-18
  rrt(100) = 1.5D-10*EXP(-14090.0D0/TGAS)
  rrt(101) = 9.1D-12*(TGAS/300.0D0)**0.18
  rrt(102) = 1.0D-11
  rrt(103) = 2.8D-12*EXP(-23400.0D0/TGAS)
  rrt(104) = 2.5D-13*EXP(-765.0D0/TGAS)
  rrt(105) = 4.6D-10*EXP(-25170.0D0/TGAS)
  rrt(106) = 1.7D-11
  rrt(107) = 2.8D-12*EXP(-25400.0D0/TGAS)
  rrt(108) = 3.3D-12*EXP(-13500.0D0/TGAS)
  rrt(109) = 4.5D-10*EXP(-18500.0D0/TGAS)
  rrt(110) = 1.2D-13*EXP(-2450.0D0/TGAS)
  rrt(111) = 2.3D-13*EXP(-1600.0D0/TGAS)
  rrt(112) = 1.5D-12*EXP(-15020.0D0/TGAS)
  rrt(113) = 4.3D-12*EXP(-3850.0D0/TGAS)
  rrt(114) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*1.0D0
  rrt(115) = 3.1D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*1.0D0
  rrt(116) = rrt(115)
  rrt(117) = rrt(115)
  rrt(118) = 3.1D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*10.D0
  rrt(119) = rrt(118)
  rrt(120) = 6.2D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*1.0D0
  rrt(121) = rrt(120)
  rrt(122) = rrt(120)
  rrt(123) = 6.2D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*12.D0
  rrt(124) = rrt(123)
  rrt(125) = 2.1D-11*(300.0D0/TGAS)**4.4*EXP(-11080.0D0/TGAS)*ANY_NEUTRAL
  rrt(126) = 1.0D-32*(300.0D0/TGAS)**0.5
  rrt(127) = rrt(126)
  rrt(128) = 1.8D-31*(300.0D0/TGAS)
  rrt(129) = rrt(128)
  rrt(130) = rrt(128)
  rrt(131) = MAX(5.8D-34*(300.0D0/TGAS)**2.8,5.4D-34*(300.0D0/TGAS)**1.9)
  rrt(132) = 1.2D-31*(300.0D0/TGAS)**1.8*1.0D0
  rrt(133) = 1.2D-31*(300.0D0/TGAS)**1.8*0.78D0
  rrt(134) = rrt(133)
  rrt(135) = 8.9D-32*(300.0D0/TGAS)**2*1.0D0
  rrt(136) = rrt(135)
  rrt(137) = 8.9D-32*(300.0D0/TGAS)**2*13.D0
  rrt(138) = rrt(137)
  rrt(139) = 8.9D-32*(300.0D0/TGAS)**2*2.4D0
  rrt(140) = 3.7D-30*(300.0D0/TGAS)**4.1*ANY_NEUTRAL
  rrt(141) = 2.8D-11
  rrt(142) = 8.0D-10
  rrt(143) = 3.0D-12
  rrt(144) = 1.0D-12
  rrt(145) = 2.4D-11
  rrt(146) = 3.0D-12
  rrt(147) = 2.3D-10
  rrt(148) = 1.6D-9
  rrt(149) = 3.3D-10
  rrt(150) = 1.0D-17
  rrt(151) = 6.3D-10
  rrt(152) = 1.0D-11
  rrt(153) = 6.6D-10
  rrt(154) = 7.0D-11
  rrt(155) = 7.0D-11
  rrt(156) = 2.9D-10
  rrt(157) = 2.9D-10
  rrt(158) = 4.0D-10
  rrt(159) = 1.0D-10
  rrt(160) = 1.2D-9
  rrt(161) = 2.0D-10
  rrt(162) = 7.0D-10
  rrt(163) = 1.0D-11
  rrt(164) = 2.6D-12
  rrt(165) = 7.0D-11
  rrt(166) = 2.0D-11
  rrt(167) = 5.0D-10
  rrt(168) = 7.4D-10
  rrt(169) = 4.0D-12
  rrt(170) = 5.0D-10
  rrt(171) = 7.0D-10
  rrt(172) = 3.0D-15
  rrt(173) = 2.5D-10
  rrt(174) = 1.0D-29*ANY_NEUTRAL
  rrt(175) = 2.0D-7*(300.0D0/TIONN)**0.5
  rrt(176) = rrt(175)
  rrt(177) = rrt(175)
  rrt(178) = rrt(175)
  rrt(179) = rrt(175)
  rrt(180) = rrt(175)
  rrt(181) = rrt(175)
  rrt(182) = rrt(175)
  rrt(183) = rrt(175)
  rrt(184) = rrt(175)
  rrt(185) = rrt(175)
  rrt(186) = rrt(175)
  rrt(187) = rrt(175)
  rrt(188) = rrt(175)
  rrt(189) = rrt(175)
  rrt(190) = rrt(175)
  rrt(191) = rrt(175)
  rrt(192) = rrt(175)
  rrt(193) = rrt(175)
  rrt(194) = rrt(175)
  rrt(195) = rrt(175)
  rrt(196) = rrt(175)
  rrt(197) = rrt(175)
  rrt(198) = rrt(175)
  rrt(199) = rrt(175)
  rrt(200) = 1.0D-7
  rrt(201) = 1.0D-7
  rrt(202) = 1.0D-7
  rrt(203) = 1.0D-7
  rrt(204) = 1.0D-7
  rrt(205) = 1.0D-7
  rrt(206) = 1.0D-7
  rrt(207) = 1.0D-7
  rrt(208) = 1.0D-7
  rrt(209) = 1.0D-7
  rrt(210) = 1.0D-7
  rrt(211) = 1.0D-7
  rrt(212) = 1.0D-7
  rrt(213) = 1.0D-7
  rrt(214) = 1.0D-7
  rrt(215) = 1.0D-7
  rrt(216) = 1.0D-7
  rrt(217) = 1.0D-7
  rrt(218) = 1.0D-7
  rrt(219) = 1.0D-7
  rrt(220) = 1.0D-7
  rrt(221) = 1.0D-7
  rrt(222) = 2.0D-25*(300.0D0/TIONN)**2.5*ANY_NEUTRAL
  rrt(223) = rrt(222)
  rrt(224) = rrt(222)
  rrt(225) = rrt(222)
  rrt(226) = rrt(222)
  rrt(227) = 2.0D-25*(300.0D0/TIONN2)**2.5*ANY_NEUTRAL
  rrt(228) = rrt(227)
  rrt(229) = rrt(227)
  rrt(230) = rrt(227)
  rrt(231) = rrt(227)
  rrt(232) = rrt(227)
  rrt(233) = rrt(227)
  rrt(234) = rrt(227)
  rrt(235) = rrt(227)
  rrt(236) = rrt(227)
  rrt(237) = rrt(227)
  rrt(238) = rrt(227)
  rrt(239) = rrt(227)
  rrt(240) = rrt(227)
  rrt(241) = rrt(227)
  rrt(242) = rrt(227)
  rrt(243) = rrt(227)
  rrt(244) = rrt(227)
  rrt(245) = rrt(227)
  rrt(246) = rrt(227)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
