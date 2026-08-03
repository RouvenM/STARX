# HLag + Fusion-Penalty für bigtime

Erweiterung des [bigtime](https://github.com/ineswilms/bigtime) R-Pakets
(Nicholson, Wilms, Bien, Matteson 2020, JMLR — HLag-Schätzer) um einen zusätzlichen **Fusions-Strafterm** auf benachbarte
Lag-Koeffizienten.
Der Fusionsterm entspricht strukturell der Group-Fused-Lasso-Strafe
(Bleakley & Vert, 2011), hier entlang der Lag- statt der Zeitachse.

## Verzeichnisstruktur

```
R/
└── hlag_fusion.R                # alle Funktionen, siehe Funktionsreferenz unten
```

## Verwendung

```r
source("R/hlag_fusion.R")

# Y: T x k Matrix von (standardisierten) Zeitreihen
fit <- fit_hlag_fusion_var(
  Y, p = 5,
  lambda1 = 0.05,      # HLag-Stärke
  lambda2 = 0.05,      # Fusions-Stärke
)
```

## bigtime-Original vs. unsere Erweiterung

Funktionen `proxcppelem` und `prox2` entsprechen `bigtime::src/hvar.cpp`, R-Nachbildung 

- `prox_hlag_vec(x, lambda)` — Proximaloperator für eine einzelne
  (Zielreihe i, Quellreihe j)-Kombination. `x` ist ein Vektor der Länge `p`
  (Koeffizienten Lag 1 bis Lag p). Entspricht `proxcppelem`.
- `prox_hlag_row(Phi_i, lambda1)` — wendet `prox_hlag_vec` unabhängig auf
  jede Spalte (Quellreihe j) einer `p x k`-Koeffizientenmatrix an. Entspricht
  `prox2`.

Fusions-Proximaloperator für die Kettenfusion benachbarter Lags (Group-Fused-Lasso entlang der Lag-Achse). Gelöst via ADMM, da wir eine Kettenstruktur im Gegensatz zur genesteten HLag-Struktur haben.
- `build_diff_operator(p)` — baut eine `(p-1) x p`-Matrix `D`, sodass
  `D %*% Phi_i` für jedes Lag-Paar `(l, l+1)` den Unterschied
  `Phi_i^(l+1) - Phi_i^(l)` liefert.
- `prox_fusion_row(V, lambda2, weights, D, ...)` — löst den Proximaloperator  
  der Kettenfusionsstrafe über ein internes ADMM-Verfahren.

Statt eines festen, nur vom Lag-Abstand abhängigen Gewichts wird zuerst eine
schnelle Ridge-Vorabschätzung berechnet, aus deren Lag-zu-Lag-Differenzen
ein individuelles Gewicht je Zielgleichung und Lag-Paar abgeleitet wird.

- `ridge_prelim_row(Zfull, y_i, p, k, ridge_lambda)` — reine
  Ridge-Regression (kein HLag, keine Fusion), nur zur Berechnung der
  adaptiven Gewichte verwendet
- `adaptive_fusion_weights(Phi_prelim, alpha_decay, weight_eps)` —
  berechnet für jedes Lag-Paar `(l, l+1)` ein Gewicht
  `w_l = exp(-alpha) / (||Phi_prelim^(l) - Phi_prelim^(l+1)||_2 + eps)`:
  kleiner Unterschied in der Vorabschätzung → größeres Gewicht → stärkerer
  Zwang zur Fusion.

Verbindung von HLAG und Fusion da die Summe beider Strafterme keinen gemeinsamen
geschlossenen Proximaloperator besitzt, jede der beiden Strafen einzeln
aber effizient lösbar ist.

- `prox_combined_row(V, lambda1, lambda2, weights, D, ...)` — berechnet den
  Proximaloperator von `lambda1 * HLag + lambda2 * Fusion` über Dykstras
  Algorithmus (alternierende Projektionen mit Korrekturtermen). Dies ist im
  Vergleich zu bigtimes Original die einzige Stelle, an der sich der
  Algorithmus tatsächlich unterscheidet: In bigtimes `FistaElem` wird an
  dieser Stelle nur `prox2` (entspricht `prox_hlag_row`) aufgerufen.


bigtimes `FistaElem`/`HVARElemAlgcpp`: Datenaufbau, Schrittweite,
beschleunigter Gradientenabstieg, Zeilen-Entkopplung über Gleichungen. Neu
ist ausschließlich die Gewichts-Weiche (fest vs. adaptiv) und der Aufruf
von `prox_combined_row()` statt eines reinen `prox_hlag_row()`-Aufrufs im
Prox-Schritt.

- `fit_hlag_fusion_var(Y, p, lambda1, lambda2, alpha_decay, adaptive, ...)` —
  zeilenweise beschleunigte proximale Gradientenmethode (FISTA). Gibt ein
  Array der Dimension `p x k x k` zurück (`Phi[lag, Quellreihe j,
  Zielgleichung i]`); das Attribut `weights_used` enthält die tatsächlich
  verwendeten Gewichte je Lag-Paar und Zielgleichung.




## Referenzen

- Nicholson, W. B., Wilms, I., Bien, J., & Matteson, D. S. (2020).
  High-dimensional forecasting via interpretable vector autoregression.
  *Journal of Machine Learning Research*, 21(166), 1-52.
- Jenatton, R., Mairal, J., Obozinski, G., & Bach, F. (2011).
  Proximal methods for hierarchical sparse coding.
  *Journal of Machine Learning Research*, 12, 2297-2334.
- Bleakley, K., & Vert, J. P. (2011). The group fused lasso for multiple
  change-point detection. *arXiv:1106.4199*.
