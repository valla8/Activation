
	program go


	real a1



	open(14,file='generic.txt',status='old')
	open(15,file='PET-14N-14O.dat',status='unknown')

	do i=1,3000
	read(14,*)a1
	if (a1.ge.0) then
	write(15,*)1000.*a1,','
	else
	write(15,*)0,','
	endif
	enddo
	



	end program go