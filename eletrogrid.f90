!APBS writes out the electrostatic potential in dimensionless units of kb T ec-1 where
!kb is Boltzmann's constant:  1.3806504 × 10−23 J K-1
!T is the temperature of your calculation in K
!ec is the charge of an electron:  1.60217646 × 10-19 C
!As an example, if you ran your calculation at 300 K, then the potential would be written out as multiples of
!
!    kb T (ec)^-1 = (1.3806504 × 10^−23 J K^-1) × (300 K) × (1.60217646 × 10^-19 C)-1
!        = (4.1419512 × 10^-21 J) × (6.241509752 × 10^18 C^-1)
!        = 2.585202 × 10^-2 J C^-1
!        = 25.85202 mV

!Modulo com as variaveis globais.
module pqr
    
  type atomtype ! cria o objeto atom
    real    :: pos(3)  !vamos usar as coordenadas em nanometros 
    real    :: charge
  end type atomtype

  type (atomtype), allocatable :: atom(:)

  real    :: posMax(3)
  real    :: posMin(3)
  real    :: GridSpacing ! Espaco entre os pontos
  real    :: BoxPadding  ! Espaco em volta da molecula
  integer :: natom
  real    :: f ! Fator de conversao eletrica kJ mol^-1 nm e^-2
  parameter (f=138.935485)  !Entao precisa converter as coordenadas para NANOMETROS
!nao usei os dois proximos
!  real      :: KbTe
!  parameter :: (KbTe=25.85202) ! 25.85202 mv ( 2.585202 10^-2 J C^-1 )
end module 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Sets default values if none is provided
subroutine setdefaults
use pqr
implicit none
  natom=0
  GridSpacing=0.05 !nanometros
  BoxPadding=0.5  !nanometros
end subroutine setdefaults




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Le um arquivo .pqr
subroutine readpqr
use pqr
implicit none
integer             ::  i
character(len=30)   ::  line

open(1,file="10ala.pqr")

! Conta quantos atomos existem
do
  read(1,"(a6)",end=10,err=10) line
  if (line(1:6)=="ATOM  ".or.line(1:6)=="HETATM") natom=natom+1
enddo
10  rewind(1)

! Vai para o 1o registro "ATOM"
do
  read(1,"(a6)",end=10,err=10) line
  if (line(1:6)=="ATOM  ".or.line(1:6)=="HETATM") then
    backspace(1)
    exit
  endif
enddo

! Aloca o vetor para guardar os atomos
allocate(atom(natom))

! Le as posicoes e cargas dos atomos
do i=1,natom
!VMD .pqr
  read(1,"(30x,3(f8.3,x),f5.3)") atom(i)%pos(:),atom(i)%charge
enddo

!converte as coordenadas para nanometros
do i=1,natom
  atom(i)%pos(1)=atom(i)%pos(1)/10 
  atom(i)%pos(2)=atom(i)%pos(2)/10 
  atom(i)%pos(3)=atom(i)%pos(3)/10 
enddo

end subroutine readpqr



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Determina os valores minimo e maximo
subroutine minmax
use pqr
implicit none
  posMax(1)=maxval(atom%pos(1))
  posMax(2)=maxval(atom%pos(2))
  posMax(3)=maxval(atom%pos(3))
  posMin(1)=minval(atom%pos(1))
  posMin(2)=minval(atom%pos(2))
  posMin(3)=minval(atom%pos(3))
end subroutine minmax



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Calcula o potencial de Coulomb em uma carga de prova
subroutine calcpot(x,y,z,pot)
use pqr
implicit none
integer           :: i
real,intent(in)   :: x,y,z
real,intent(out)  :: pot
real              :: dist !distace

pot=0.0d0
!calcula a distancia do ponto com a ptn,
do i=1,natom
dist = sqrt (                         &
            ( x-atom(i)%pos(1) )**2 + &
            ( y-atom(i)%pos(2) )**2 + &
            ( z-atom(i)%pos(3) )**2   &
            )

!calcula o potencial
pot = pot + ( f *atom(i)%charge/dist ) ! saida em kJ/mol

enddo
!mudando o valor final para caber no .pqr
pot = pot / 100 ! kJ/mol *10^2 

end subroutine calcpot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Writes the grid
!> O VMD esta lendo normalmente mas o pymol nao :(
!> Talvez seja um problema com o formato de saída.
subroutine writegrid
use pqr
implicit none
real     :: x,y,z    ! Posicao X,Y e Z
integer  :: cx,cy,cz ! Contadores para X,Y e Z
integer  :: nx,ny,nz ! quantos pontos em X,Y e Z
real     :: pot      ! potencial eletrostatico
integer  :: i=0
real,allocatable     :: allpot(:)! potencial eletrostatico

!calcula o numero de pontos em cada direcao
!1) diferenca entre o maximo e minimo
nx = nint (( abs( posMax(1) - posMin(1) ) + BoxPadding*2 ) / GridSpacing )
ny = nint (( abs( posMax(2) - posMin(2) ) + BoxPadding*2 ) / GridSpacing )
nz = nint (( abs( posMax(3) - posMin(3) ) + BoxPadding*2 ) / GridSpacing )

allocate(allpot(nx*ny*nz))

!linhas para o formato .dx
write(*,2000) "object 1 class gridpositions counts ",nx,ny,nz
write(*,2001) "origin ",(posMin(1)-BoxPadding)*10,(posMin(2)-BoxPadding)*10,(posMin(3)-BoxPadding)*10 !em Angstrons
write(*,2001) "delta ",GridSpacing*10,0.0,0.0
write(*,2001) "delta ",0.0,GridSpacing*10,0.0
write(*,2001) "delta ",0.0,0.0,GridSpacing*10
write(*,2000) "object 2 class gridconnections counts ",nx,ny,nz
write(*,2002) "class array type double rank 0 items ",nx*ny*nz," data follows"
2000 format(A,3i5)
2001 format(A,3f8.3)
2002 format(A,i5,A)

!escreve os pontos ! z->y->x (eh a ordem do formato .dx)
  x=posMin(1)-BoxPadding
  do cx=1,nx

    y=posMin(2)-BoxPadding
    do cy=1,ny

      z=posMin(3)-BoxPadding
      do cz=1,nz
      !aqui vc podia colocar uma restricao para ver se o ponto esta muito
      !proximo da proteina. Ex. n escreve nada abaixo de 1.4 Angstrons (raio
      !aproximado de uma molecula de agua.

        call calcpot(x,y,z,pot)  !calcula o potencial eletrostatico antes de escrever
        i=i+1
        allpot(i)=pot
        
        z=z+GridSpacing 
      enddo
 
     y=y+GridSpacing 
    enddo 

  x=x+GridSpacing
  enddo

write(*,"(3(E14.6,x))") allpot
write(*,"(A)") 'attribute "dep" string "positions"'
write(*,"(A)") 'object "regular positions regular connections" class field'
write(*,"(A)") 'component "positions" value 1'
write(*,"(A)") 'component "connections" value 2'
write(*,"(A)") 'component "data" value 3'
!1000   format("ATOM      0  H   GRID    1    ",3f8.3,f8.4)
end subroutine writegrid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program eletrogrid
use pqr
implicit none

  call setdefaults
  call readpqr
  call minmax
  call writegrid

end
