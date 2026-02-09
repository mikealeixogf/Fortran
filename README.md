# Métodos de Integração Numérica em Fortran

Este repositório contém mais de 300 códigos em Fortran
relacionados a métodos de integração numérica, incluindo
integrais gaussianas, uniformes e outras rotinas.

## Estrutura
- `main_integral.f90` → programa principal
- `integrais_gaussianas.f90` → métodos de Gauss
- `integrais_uniformes.f90` → métodos uniformes
- demais arquivos → subrotinas auxiliares

## Como compilar
```bash
gfortran main_integral.f90 -o integral.exe
./integral.exe
