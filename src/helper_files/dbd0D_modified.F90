!------------------------------------------------------------------------------------------
! BASED ON :
!------------------------------------------------------------------------------------------
! DBD0D TEST CASE
! ZDPLASKIN
!
! June 2008, by S.Pancheshnyi
!

!------------------------------------------------------------------------------------------
! DeNOx by DBD
!------------------------------------------------------------------------------------------
! ZDPLASKIN DBD0D MODEL 
! 
! April 2023, by N.Bente


program DBD0D

  use ZDPlasKin

  implicit none

!------------------------------------------------------------------------------------------
!
! Input Gas Configuration
!
!------------------------------------------------------------------------------------------
!
! gas pressure, temperature and density
!
  real(8), parameter :: no_ppm = 100 ! Input NO concentration, [ppm]
  real(8), parameter :: o2_ppm = 25000 ! Input O2 concentration, [ppm]
  real(8), parameter :: gas_temperature = 300 ! Gas temperature, [K]
  real(8), parameter :: gas_pression = 1 ! Gas pressure, [atm]
  
! Calcul of gas pressure
!
  real(8), parameter :: no_initial_density = ((101325*(no_ppm/1E6)*6.022E23)/(gas_temperature*8.314))/1E6
  real(8), parameter :: o2_initial_density = ((101325*(o2_ppm/1E6)*6.022E23)/(gas_temperature*8.314))/1E6
  real(8), parameter :: n2_initial_density = ((101325*(1-((no_ppm+o2_ppm)/1E6))*6.022E23)/(gas_temperature*8.314))/1E6
  
 
!------------------------------------------------------------------------------------------
!
! DBD Parameters
!
!------------------------------------------------------------------------------------------
!
! 
!
  real(8), parameter :: gap = 1.6E-3 ! DBD gap, [m]
  real(8), parameter :: A = 4.8E-3 ! Current cross section, [m2]
  
  
  integer, parameter          :: ipulses = 1
!
! filename
!
 character(*), parameter     :: file_in = 'data_in.dat'
  
!------------------------------------------------------------------------------------------
!
! local variables
!
  integer :: i, j, idata_max
  double precision :: time, dtime, EN, NE, Vgas, Igas, neutrals_d, e_mobility, no, no2
  double precision, allocatable :: data_in(:,:)

  real(8) :: start_time, end_time

!------------------------------------------------------------------------------------------
!
! write
!

  write(*,'(/,A)') ' DBD0D TEST CASE'
  write(*,*) 'NO Concentration : ', no_ppm, '[ppm]'
  write(*,*) 'NO Initial density : ', no_initial_density, '[cm-3]'
  write(*,*) 'N2 Initial density : ', n2_initial_density, '[cm-3]'


!------------------------------------------------------------------------------------------
!
! load profiles of field and electron density from external file
!

  open(1,file=file_in)
  read(1,*)

  idata_max = 0
  do while( .true. )
    read(1,*,end=11) dtime, EN, NE
    idata_max = idata_max + 1
  enddo

11 write(*,'(x,A,i0,A)') 'found ', idata_max, ' points in ' // file_in // ' file'

  if(idata_max .le. 1) then
    write(*,*) 'wrong or missing data in the file'
    close(1)
    goto 999
  else
    allocate(data_in(idata_max,3))
    rewind(1)
    read(1,*)
    do i = 1, idata_max
      read(1,*) data_in(i,1:3)
    enddo
    close(1)
  endif

  time = data_in(1,1)

!------------------------------------------------------------------------------------------
!
! initialization of ZDPlasKin; initial densities and gas temperature
!
  
  call ZDPlasKin_init()
  call ZDPlasKin_set_density('N2',n2_initial_density)
  call ZDPlasKin_set_density('NO',no_initial_density)
  call ZDPlasKin_set_density('O2',o2_initial_density)
  call ZDPlasKin_set_density('E',LDENS_CONST=.true.) ! fixed electron density
  call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature)
  
  call ZDPlasKin_set_config(QTPLASKIN_SAVE=.false.)


!------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------
!
  

!
! open file for output and write headers
!

!  open(1,file=file_out)
!  write(1,'(99A14)') 'TIME', 'E/N', ( trim(species_name(i)), i = 1, species_max )
!  write(*,'(99A14)') 'TIME'
  
!  open(2,file=file_NO)
!  write(2,'(99A14)') 'TIME', 'NO [ppm]'
  
!  open(3,file=file_NO2)
!  write(3,'(99A14)') 'TIME', 'NO2 [ppm]'
  
!  open(4,file=file_efficiency)
!  write(4,'(99A14)') 'EFFICIENCY [%]'
  
!  open(5,file=file_efficiency_time)
!  write(5,'(99A14)') 'TIME', 'EFFICIENCY [%]'
  
!  open(6,file=file_selectivity)
!  write(6,'(99A14)') 'SELECTIVITY [%]'

!------------------------------------------------------------------------------------------
!
! time integration
!
!------------------------------------------------------------------------------------------

  call CPU_TIME(start_time)

  call ZDPlasKin_get_density_total(ALL_NEUTRAL=neutrals_d)
  write(*,*) 'Initial gas density : ', neutrals_d, '[cm-3]'
  call ZDPlasKin_get_conditions(ELEC_MOBILITY_N=e_mobility)
  write(*,*) 'Initial elec reduced mobility : ', e_mobility
  

  do i = 1, ipulses
  
    write(*,*) '========== ',i,'/',ipulses,' =========='

    do j = 1, idata_max - 1

      dtime = data_in(j+1,1) - data_in(j,1)
      Vgas    = data_in(j,2)
      Igas    = data_in(j,3)
	  
	  IF (Vgas==0) THEN
        EN = 0
		NE = 0
      ELSE
	    EN = ((Vgas)/(gap*1E2*neutrals_d))/(1E-17)
		NE = (Igas)/(A*1E4*1.602E-19*e_mobility*EN*1E-17)
      END IF
	  
      call ZDPlasKin_set_conditions(REDUCED_FIELD=EN)
      call ZDPlasKin_set_density('E',NE)

      call ZDPlasKin_timestep_explicit(time,dtime,1d-1,1d-1)
      time = time + dtime
      
	  call ZDPlasKin_get_density_total(ALL_NEUTRAL=neutrals_d)
	  call ZDPlasKin_get_conditions(ELEC_MOBILITY_N=e_mobility)
	  call ZDPlasKin_get_density('NO',no)
	  call ZDPlasKin_get_density('NO2',no2)
      
      IF( MOD(j,idata_max/100) == 0 .OR. j==1 .OR. j==idata_max - 1 ) THEN
      
            !write(2,'(99(1pe14.5))') time, (no/neutrals_d)*1E6
            !write(3,'(99(1pe14.5))') time, (no2/neutrals_d)*1E6
            !write(5,'(99(1pe14.5))') time, 1-((((no/neutrals_d)*1E6)+((no2/neutrals_d)*1E6))/no_ppm)
            write(*,*) 'Progression :',(j*100)/(idata_max - 1),'%'
            
      
      END IF
      
      IF(MOD(j,idata_max/10000) == 0 .OR. j==1 .OR. j==idata_max - 1) THEN
          call ZDPlasKin_write_qtplaskin(time, .true.)
          !write(*,*) 'Logro la escritura'
      END IF
	  
      !write(1,'(99(1pe14.5))') time, EN, density(:)
	  

    enddo


    !write(*,'(99(1pe14.5))') time
	

  enddo

  
!write(4,'(99(1pe14.5))') 1-((((no/neutrals_d)*1E6)+((no2/neutrals_d)*1E6))/no_ppm)
!write(6,'(99(1pe14.5))') ((no2/neutrals_d)*1E6)/no_ppm

!------------------------------------------------------------------------------------------
!
! close file and end
!
  !close(1)

999 continue
  if( allocated(data_in) ) deallocate(data_in)

  write(*,*) '#SIMULATION ENDED'
  call CPU_TIME(end_time)
  write(*,*) 'Simulation time:',end_time - start_time, 's'
  !write(*,*)' NO [ppm] : ',(no/neutrals_d)*1E6
  !write(*,*)' NO2 [ppm] : ',(no2/neutrals_d)*1E6
  !write(*,*)' deNOx efficiency : ',1-((((no/neutrals_d)*1E6)+((no2/neutrals_d)*1E6))/no_ppm)
  !write(*,*)' NO2 selectivity : ',((no2/neutrals_d)*1E6)/no_ppm
end program DBD0D