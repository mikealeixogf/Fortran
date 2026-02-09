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
