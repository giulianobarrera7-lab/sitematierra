# Calculadora de Puesta a Tierra (PAT) — AEA 90364

Calculadora web de predimensionamiento de sistemas de puesta a tierra según la
Reglamentación AEA 90364 (Asociación Electrotécnica Argentina), esquema TT.

No requiere instalación ni build: es HTML/CSS/JS puro. Para usarla, abrí
`index.html` en cualquier navegador, o serví la carpeta con un servidor
estático (por ejemplo `python3 -m http.server`).

## Contenido

- **1. Resistividad y Jabalina** — resistividad del suelo por el método de
  Wenner (con y sin corrección por profundidad de las picas), resistencia de
  una jabalina individual (fórmula de Dwight), jabalinas en paralelo con
  resistencia de acoplamiento mutuo, cantidad de jabalinas necesarias para
  cumplir un objetivo, y distancia mínima de "tierra lejana" (≥10·Re).
- **2. Malla / Electrodo** — dos modos: malla reticulada (fórmulas de Sverak
  e IEEE Std 80, y Laurent como comparación) o cable desnudo enterrado +
  jabalinas (método de las imágenes + Dwight, combinados con acoplamiento
  mutuo). Incluye recomendaciones automáticas de interruptor diferencial,
  sección mínima del conductor de protección (PE) y del conductor de tierra
  enterrado, y código de colores normalizado.
- **3. Mediciones (Telurómetro)** — perfil de resistividad con picas a
  profundidad variable (ρ = 2π·a·R) con estadística de las lecturas, y el
  método del 62% (caída de potencial) para medir la resistencia de una PAT ya
  instalada, con verificación de "curva plana" y de la distancia mínima a la
  pica de corriente.
- **4. Protocolo 62% (multipunto)** — planilla para relevar varios puntos de
  medición de un proyecto, con estado de cumplimiento global.
- **5. Referencias y Normas** — tablas de resistividad de terrenos (AEA
  Anexo 771-C, Tablas 771-C.VIII y 771-C.IX), tabla de distancias de tierra
  lejana según jabalina, fórmulas utilizadas y criterios normativos.

## Fórmulas principales

```
Resistividad (Wenner, simplificada):     ρ = 2·π·a·R
Resistividad (Wenner, con profundidad):  ρ = 4πaR / (1 + 2a/√(a²+4b²) − a/√(a²+b²))
Jabalina (Dwight):                       R₁ = ρ/(2πL)·(ln(4L/a) − 1),  a = Ø/2
Paralelo de n jabalinas:                 Rn = (R₁ + (n−1)·Rm)/n,  Rm = ρ/(2π·s)
Malla (Sverak / IEEE Std 80):            R = ρ·[1/L_T + 1/√(20A)·(1+1/(1+h√(20/A)))]
Cable enterrado (método de imágenes):    R = ρ/(4πL)·[ln(2L/a) + ln(L/h) − 2]
Tierra lejana (AEA):                     d ≥ 10·Re,  Re = L/(ln(4L/a) − 1)
Sección del PE (Tabla 771.18.III):       S≤16→S ; 16<S≤35→16 ; S>35→S/2
Tensión de contacto:                     Uc = R × IΔn ;  Ra máx = Uʟ / IΔn
```

## Criterios normativos aplicados

- Esquema TT — límite general de predimensionamiento: Rpat ≤ 40 Ω (AEA
  90364-7-771, Anexo 771-C).
- Con protección diferencial: Rpat ≤ 10 Ω (preferentemente ≤ 5 Ω).
- Sin protección diferencial adecuada: tensión de contacto indirecto ≤ 24 V.
- Interruptor diferencial de uso general: IΔn ≤ 30 mA (IRAM 2301).
- Toma de tierra de protección ubicada a ≥10·Re de la toma de tierra de
  servicio.
- Conductor de protección (PE): cobre aislado bicolor verde-amarillo,
  sección mínima 2,5 mm².

## Advertencia

Esta herramienta es de **predimensionamiento**. La resistencia de puesta a
tierra real depende de la heterogeneidad del suelo y debe verificarse siempre
en obra con un telurímetro homologado (Norma IRAM 2281, Parte I) antes de dar
por aprobada una instalación.
