!Obs decidi usar subrotinas pq fica mais facil de integrar com python
!esse programa avalia as integrais uniformes e gaussianas
!limite inferior:l_inf
!limite superior:l_sup
!area: areasn
!ìndices:n
program integral_numerica !programa principal

implicit none !todas as variaveis tem que ser declaradas
real(8)::l_inf,l_sup !limite de integração
real(8)::area,area2,area3,area4 !área calculada método uniforme
real(8)::area5 !área calculada método quadratura
integer:: n
n=6
l_inf=0.1d0
l_sup=10.0d0
call retangular(f,l_inf,l_sup,n,area)
call trapezoidal(f,l_inf,l_sup,n, area2)
call simpson13(f,l_inf,l_sup,n,area3)
call simpson38(f,l_inf,l_sup,n,area5)
write(*,*)"A integral pela quadratura Retangular foi de: ",area
write(*,*)"A integral pela regra do Trapezio foi de: ",area2
write(*,*)"A integral pela regra de Simpson 1foi de: ",area3
write(*,*)"A integral pela regra do Simpson 2 foi de: ",area4
print*,"======================================================================="
call gauss_legendre(f,l_inf,l_sup,n,area5)


write(*,*)"A integral pela quadratura Gauss(legendre): ",area5


!função de todas
contains 
  double precision function f(x)
    double precision, intent(in) :: x
    f = x
  end function f

end program integral_numerica


!este programa irá usar regras de integração númerica do tipo Newton-Cotes (passo uniforme)
!1-Retangular
!2-Trapézio
!3-Simpson 1 e Simpson 2
!Estas integrais são de malha fixa
!São proprias para integrais fechadas
!troquei module por subroutine

!função Retangular-subrotina
subroutine retangular(f,l_inf,l_sup,n,integral)  !é o funcional que recebe a função f(x)
implicit none
!Vamos usar a interface para calular a área da curva
  interface !recebe a função a ser integrada
    double precision function f(x)!função a ser transformada
      double precision, intent(in) :: x !é o argumento da função
    end function f
    end interface
real(8),intent(in)::l_inf,l_sup !limites de integração analitica
real(8),intent(out)::integral !saida
integer::i !contador de integral
real(8)::dx,xi !elementos do somatório
integer, intent(in)::n !quantidade de pontos a serem integradas discretamente. Intent (in) é a entrada
dx=(l_sup-l_inf)/n !é o tamanho do passo: o infinitesimal de integração

integral=0.0
!---------------INICIO DO LOOP DE SOMATÓRIO---------------------------------------
do i= 0, n-1
    xi=l_inf+i*dx
    integral=integral+f(xi)
end do
    integral=integral*dx


end subroutine retangular

!pronta para começar outra integral

!==============================================================================================
!---------------------função trapezoidal-------------------------------------------------
subroutine trapezoidal(f,l_inf,l_sup,n,integral)!é o funcional que recebe a função f(x)
implicit none
!Vamos usar a interface para calular a área da curva
  interface !recebe a função a ser integrada
    double precision function f(x)!função a ser transformada
      double precision, intent(in) :: x !é o argumento da função
    end function f
    end interface
real(8),intent(in)::l_inf,l_sup !limites de integração analitica
real(8),intent(out):: integral
integer::i !contador de integral
real(8)::dx,xi !elementos do somatório
integer, intent(in)::n !quantidade de pontos a serem integradas discretamente
dx=(l_sup-l_inf)/n !é o tamanho do passo: o infinitesimal de integração

integral=0.0
!---------------INICO DO LOOP DE SOMATÓRIO---------------------------------------
do i= 0, n
    xi=l_inf+i*dx
    if(i==0 .or. i==n) then
        integral=integral +f(xi)
    else
        integral=integral+2.0d0*f(xi)
    end if
end do
    integral=(integral)*(dx/2)


end subroutine trapezoidal
!=========================Fim regra do Traopezio================================================


!===========================Início regra de Simpson 1===========================================
subroutine simpson13(f, l_inf, l_sup, n,integral)
implicit none

interface
    real(8) function f(x)
        real(8), intent(in) :: x
    end function f
end interface

real(8), intent(in) :: l_inf, l_sup
real(8), intent(out)::integral
integer, intent(in) :: n
integer :: i
real(8) :: dx, xi

! Simpson exige n par
if (mod(n,2) /= 0) then
    print *, "Erro: n deve ser par para Simpson 1/3"
    stop
end if

dx = (l_sup - l_inf) / n
integral = 0.0d0

do i = 0, n
    xi = l_inf + i*dx

    if (i == 0 .or. i == n) then
        integral = integral + f(xi)

    else if (mod(i,2) == 0) then
        integral = integral + 2.0d0 * f(xi)

    else
        integral = integral + 4.0d0 * f(xi)
    end if

end do

integral = integral * dx / 3.0d0

end subroutine simpson13
!===========================Fim da regra de Simpson 1===============================

!===========================Início regra de simpson 2===================================
subroutine simpson38(f, l_inf, l_sup, n,integral) 
implicit none

interface
    real(8) function f(x)
        real(8), intent(in) :: x
    end function f
end interface

real(8), intent(in) :: l_inf, l_sup
real(8),intent(out)::integral
integer, intent(in) :: n
integer :: i
real(8) :: dx, xi

if (mod(n,3) /= 0) then
    print *, "Erro: n deve ser multiplo de 3 para Simpson 3/8"
    stop
end if

dx = (l_sup - l_inf) / n
integral = 0.0d0

do i = 0, n
    xi = l_inf + i*dx

    if (i == 0 .or. i == n) then
        integral = integral + f(xi)

    else if (mod(i,3) == 0) then
        integral = integral + 2.0d0 * f(xi)

    else
        integral = integral + 3.0d0 * f(xi)
    end if

end do

integral = integral * 3.0d0 * dx / 8.0d0

end subroutine simpson38
!===========================Fim da regra de Simpson 2=======================================
subroutine gauss_legendre(f, l_inf, l_sup, n, integral)   ! Declara a subrotina Gauss-Legendre
  implicit none                                           ! Usa declaração implícita padrão

  ! Interface da função a ser integrada (f(x))
  interface
    real(8) function f(x)                                 ! Define a função f(x) que será passada
      real(8), intent(in) :: x                            ! Argumento de entrada da função
    end function f
  end interface

  ! Argumentos da subrotina
  real(8), intent(in) :: l_inf, l_sup                     ! Limites inferior e superior da integral
  integer, intent(in) :: n                                ! Número de pontos de Gauss-Legendre
  real(8), intent(out) :: integral                        ! Resultado da integral (saída)

  ! Variáveis locais
  integer :: i, j, m                                      ! Contadores de laços
  real(8) :: dx, t                                        ! dx = transformação de intervalo, t = ponto avaliado
  real(8) :: z, z1, p1, p2, p3, pp, pi                    ! Variáveis auxiliares para cálculo das raízes/pesos
  real(8), allocatable :: xg(:), wg(:)                    ! Vetores para abscissas (xg) e pesos (wg)

  ! Inicializa valor de pi
  pi = 4.0d0 * atan(1.0d0)                                ! Calcula π usando arctan

  ! Aloca memória para pontos e pesos
  allocate(xg(n), wg(n))                                  ! Cria vetores de tamanho n

  ! Calcula pontos e pesos de Gauss-Legendre
  m = (n + 1) / 2                                         ! Número de raízes a calcular (metade, por simetria)
  do i = 1, m                                             ! Loop para cada raiz
    z = cos(pi * (dble(i) - 0.25d0) / (dble(n) + 0.5d0))  ! Chute inicial para raiz do polinômio de Legendre
    z1 = z + 1.0d0                                        ! Valor inicial para comparação

    ! Método de Newton-Raphson para refinar raiz
    do while (abs(z - z1) > 1.0d-14)                      ! Itera até convergir
      p1 = 1.0d0                                          ! Polinômio de Legendre inicial
      p2 = 0.0d0                                          ! Valor auxiliar
      do j = 1, n                                         ! Loop para calcular polinômio de Legendre
        p3 = p2                                           ! Atualiza valores anteriores
        p2 = p1
        p1 = ((2*j - 1)*z*p2 - (j - 1)*p3) / j            ! Recorrência de Legendre
      end do
      pp = dble(n) * (z*p1 - p2) / (z*z - 1.0d0)          ! Derivada do polinômio
      z1 = z                                              ! Atualiza valor anterior
      z = z1 - p1/pp                                      ! Passo de Newton-Raphson
    end do

    ! Guarda abscissas simétricas
    xg(i) = -z                                            ! Raiz negativa
    xg(n+1-i) = z                                         ! Raiz positiva

    ! Calcula pesos correspondentes
    wg(i) = 2.0d0 / ((1.0d0 - z*z) * pp * pp)             ! Peso associado à raiz
    wg(n+1-i) = wg(i)                                     ! Peso simétrico
  end do

  ! Integração propriamente dita
  integral = 0.0d0                                        ! Inicializa integral
  dx = (l_sup - l_inf) / 2.0d0                            ! Fator de transformação do intervalo
  do i = 1, n                                             ! Loop sobre todos os pontos
    t = dx * xg(i) + (l_sup + l_inf) / 2.0d0              ! Transforma ponto de [-1,1] para [l_inf,l_sup]
    integral = integral + wg(i) * f(t)                    ! Soma peso * f(t) à integral
  end do
  integral = integral * dx                                ! Multiplica pelo fator dx final

  ! Libera memória
  deallocate(xg, wg)                                      ! Desaloca vetores

end subroutine gauss_legendre                             ! Fim da subrotina
