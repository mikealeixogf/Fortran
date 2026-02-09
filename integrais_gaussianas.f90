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
