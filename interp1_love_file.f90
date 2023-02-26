  subroutine interp1_love_file(yfile)
    use some_functions,only: count_lines,spline3
    implicit none
    character(len=*) :: yfile
    real(kind=8),allocatable :: yDiscrete(:,:),intePara(:)
    real(kind=8),allocatable :: yContinue(:,:)
    integer :: n_line,m_column,ii,jj,max_degree,min_degree,need_terms

    WRITE(*,*) 'in interp1_love_file ',trim(yfile)



    n_line = count_lines(yfile)

    !write(*,*) n_line
    !stop


    allocate(yDiscrete(n_line,9))
    allocate(intePara(n_line))

    open(unit=704,file=yfile,action='read')
    do ii=1,n_line
      read(704,*) yDiscrete(ii,1:9)
      write(*,'(1I10,8(2X,1E20.10E5))') int(yDiscrete(ii,1)), yDiscrete(ii,2:9)
    enddo
    close(unit=704)

    min_degree=int(yDiscrete(1,1))
    max_degree=int(yDiscrete(n_line,1))

    need_terms=max_degree-min_degree+1

    write(*,*) 'need_terms=',need_terms


    allocate(yContinue(need_terms,9))


    ii=0
    do jj=min_degree,max_degree,1
      ii=ii+1
      yContinue(ii,1)=real(jj,kind=8)
    enddo

    do jj=2,9
      ! call spline3( n_line,  yDiscrete(:,1), yDiscrete(:,jj),&
      !         need_terms-1, yContinue(2:,1), yContinue(2:,jj))

      call SPLINE(yDiscrete(:,1), yDiscrete(:,jj),n_line,1.0d40,1.0d40,intePara)
      do ii=1,need_terms
        call SPLINT(yDiscrete(:,1), yDiscrete(:,jj),intePara,n_line,yContinue(ii,1),yContinue(ii,jj))
      enddo

    enddo

    open(unit=72,file='DSLOV33.DAT',action='write')

    do ii=1,need_terms
      write(72,'(1I10,8(2X,1E20.10E5))') int(yContinue(ii,1)), yContinue(ii,2:9)
      write( *,'(1I10,8(2X,1E20.10E5))') int(yContinue(ii,1)), yContinue(ii,2:9)
    enddo
    close(unit=72)

    deallocate(yContinue)
    deallocate(yDiscrete)



  end subroutine interp1_love_file
