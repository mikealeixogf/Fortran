# Integrais em Fortran

Este projeto contém implementações de métodos numéricos de integração em **Fortran** (Retangular, Trapézio, Simpson 1/3, Simpson 3/8 e Gauss-Legendre).  
O objetivo é didático: mostrar como organizar subrotinas em Fortran e como compilar/rodar no VS Code ou usar no Python via **f2py**.

---

## Estrutura do projeto
- `main_integrais.f90` → programa principal com todas as subrotinas integradas (didático, roda direto no VS Code).
- `integrais_uniformes.f90` → subrotinas de integração por métodos uniformes (Retangular, Trapézio, Simpson).
- `integrais_gaussianas.f90` → subrotina de integração Gauss-Legendre.
- `.vscode/tasks.json` → configuração para compilar e rodar automaticamente no VS Code.

---

## Como compilar no VS Code
1. Instale o compilador **gfortran** (parte do MinGW ou TDM-GCC no Windows).
2. Abra o projeto no VS Code.
3. Pressione **Ctrl+Shift+B** → escolha a task **Compilar**.
4. O executável `main_integral.exe` será gerado.
5. Para rodar:
   ```powershell
   .\main_integral.exe

