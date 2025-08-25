# Elektrofit Refactor & Architektur Leitfaden

Stand: 2025-08-19 (aktualisiert nach IO-Modularisierung, Step6 Paritäts-Fix + Outlier-Feature (EXPERIMENTELL), Net-Charge Normalisierung)

Dieses Dokument dient als lebende Referenz für die laufende Struktur‑Konsolidierung und Architektur‑Verbesserungen. Ziel: Klare Schichten, minimale Kopplung, eindeutige Verantwortlichkeiten, konsistente Namensgebung und wartbare Public API.

---
## Leitprinzipien
1. Single Source of Truth: Konfiguration & Layering nur an einer Stelle (compose/build Snapshot).
2. Dünne Workflows: Step-* Dateien orchestrieren nur, keine Geschäftslogik.
3. Domain vor Technik: Fachlogik (Symmetrie, Ladungen, Konformer) in `domain/*` statt in generischen Sammelordnern.
4. Infrastruktur kapseln: Logging, Scratch, Config, externe Tool-Adapter getrennt in `infra/` bzw. `adapters/`.
5. Öffentliche vs. Interne API: Explizit via `__all__` + Namenskonvention (`_internal`).
6. Predictable Imports: Pipelines importieren Services, nicht Einzel-Funktionsmodule querbeet.
7. Testbarkeit: Domain-Funktionen ohne Seiteneffekte (kein global logging, kein cwd change) testbar.
8. Idempotenz & Reentrance: Wiederholter Aufruf eines Steps verändert nur notwendige Artefakte.
9. Minimaler Log-Lärm: Keine redundanten Merge-Logs pro Konformer.
10. Optionalität: Plotting & Visualisierung als optionale Extras (lazy import / optional dependency).

---
## Ziel-Schichten (Soll-Struktur)
```
electrofit/
  api/                # stabile Fassaden (optional)
  cli/                # Argument-Parsing + Dispatcher
  pipeline/steps/     # step1.py ... stepN.py (nur Orchestrierung)
  domain/
    symmetry/
    charges/
    conformers/
    prep/
  infra/
    config_snapshot.py
    logging.py
    scratch_manager.py
  adapters/
    gromacs.py
    gaussian.py (potenziell)
    acpype.py (potenziell)
  io/
    mol2.py
    resp.py
    files.py
  viz/ (oder plotting/)
  internal/ (falls nötig für Migrationsphase)
```

---
## Aktueller Ist-Zustand (Kurz)
- `workflows/` enthält Schritte + Snapshot-Builder (gemischt: Orchestrierung & Infrastruktur).
- `core/` mischt charges/symmetry/prep.
- `external/` nur GROMACS Adapter (Name zu generisch).
- `logging.py` top-level; `scratch/manager.py` eigenes Mini-Paket.
- `utils_curly_brace.py` unklare Funktion / Name.

---
## Phasenplan
### Phase 1 – Reorganisation (Low Risk)
- [x] Verschiebe `workflows/snapshot.py` → `infra/config_snapshot.py` (Export: `compose_snapshot` (Alias für `build_snapshot_with_layers`), `CONFIG_ARG_HELP`).
- [x] Umbenenne / entferne `external/` → `adapters/` (alter Namespace gelöscht, kein Shim nötig).
- [x] Verschiebe `logging.py` → `infra/logging.py`.
- [x] Verschiebe `scratch/manager.py` → `infra/scratch_manager.py`; entferne Paket `scratch`.
- [x] Verschiebe `utils_curly_brace.py` → `viz/curly_brace.py` (finaler Funktionsname `draw_curly_brace`).
- [x] Migriere frühere `plotting/helpers.py` Funktionen nach `viz/helpers.py` und entferne Altdatei.
- [x] Füge Deprecation-Shim in `plotting/__init__.py` hinzu (re-exportiert `viz`, gibt `DeprecationWarning`).
- [ ] Entferne Shim nach mindestens einem Minor Release (Release-Note ankündigen).
- [ ] Bereinige veraltete Einträge in `rewrite_imports.py` (alte Pfade `electrofit.plotting.helpers`).

### Phase 2 – Domain-Isolation
Zwischenstand aktualisiert: `equiv_groups` und `write_symmetry` vollständig in `domain.symmetry`; `process_conformer` und `process_initial` domänisiert; Legacy-Shims vorhanden mit DeprecationWarning & Forwarding für Tests. Step6 Aggregationslogik bereinigt (acpype Parität in allen Pfaden, experimenteller Outlier-IQR Filter + Net Charge Normalisierung) – Vorbereitung auf spätere Auslagerung Plotting.

- [x] Verschiebe `core/equiv_groups.py` → `domain/symmetry/equiv_groups.py` (Funktions-API `generate_equivalence_groups`).
- [x] CLI Wrapper `electrofit-create-equiv-groups` hinzugefügt.
- [x] Extrahiere gemeinsame RESP‑Symmetrie Hilfsfunktion (jetzt `domain/symmetry/resp_constraints.py`, genutzt von prep & charges).
- [x] Verschiebe `core/process_conform.py` → `domain/charges/process_conformer.py` (erste Zerlegung erfolgt; weitergehende Feingliederung offen).
- [x] Verschiebe `core/process_initial_structure.py` → `domain/prep/process_initial.py` (Injection-Punkt `run_acpype` + Shim Warnung).
- [ ] Ergänze Services: `services/conformer_sampling.py`, `services/config_service.py` (noch offen).

Design Notes für Zerlegung (Vorab):
1. Keine `print` Statements in Domain – nur Logger.
2. Scratch Handling via dünne Orchestrator-Funktion, Pure Functions erhalten explizite Pfade.
3. Externe Tools (gaussian, respgen, bondtype, espgen) über Adapter‐Layer kapseln (künftige Phase, jetzt nur leichte Isolierung durch Wrapper‑Funktionen).
4. Fehlerklassen vereinheitlichen (`ChargeComputationError`, `SymmetryDerivationError`).

### Phase 3 – API & Services
- [ ] `pipeline/steps/stepX.py` einführen; vorhandene `stepX_*.py` umbenennen / verschieben.
- [ ] Öffentlicher Einstieg: `electrofit.api` definiert stabile Funktionen (z.B. `run_step(step:int, project:Path, config:Path|None)`).
- [ ] Interne Module prefix `_`; unveröffentlichte Objekte entfernen aus `__all__`.
- [ ] AggregationResult Dataclass für Step6 (ersetzt implizite Rückgabewerte) + Plot-Service entkoppelt.

### Phase 4 – Bereinigung & Deprecation
- [ ] Entferne Legacy-Dateien / Aliase (`run_process_*`).
- [ ] Entferne `merge_into_snapshot` aus Public Surface (nur intern verfügbar).
- [ ] Einheitliche Schritt-Namen: `step1.py`, `step2.py` ...
- [ ] Dokumentiere neue Modulkarte in README / docs/architecture.md.
- [ ] Entferne Deprecation-Shim (`plotting` → `viz`).
- [ ] Entferne verbliebene Alt-Import-Rewrites.
- [ ] Entferne Re-Exports der Legacy binary_mapping Helfer; ersetze durch direkte Nutzung oder entferne Bedarf.
- [ ] DeprecationWarnings für verbleibende Übergangsfunktionen in `io/files.py` aufnehmen (Zeitstempel + Entferndatum).

### Phase 5 – Optional Enhancements
- [ ] Optionales Extra: `pip install electrofit[viz]` für plotting.
- [ ] Typsicherheit: mypy-Konfiguration + strictere Typen in infra & domain.
- [ ] Structured Logging (JSON) optional für CI.
- [ ] Performance Messpunkte (Timer-Dekoratoren) in Services.

---
## Migrationsstrategie
1. Verschieben mit Re-Exports + DeprecationWarnings.
2. CI Test-Suite sichern (jede Phase grün bevor nächste startet).
3. Release Notizen: „0.2.x – interne Reorganisation, API unverändert / 0.3.0 – neue API-Fassade“.
4. Nach zwei Minor-Releases: Entfernen der Deprecation-Re-Exports.

---
## Naming Guidelines
- Dateien: snake_case, prägnant (keine überflüssigen Präfixe wie `step3_start_sim` → `step3.py`).
- Funktionen: Verben für Aktionen (`compose_snapshot`, `sample_conformers`), Substantive für reine Modelle.
- Keine Doppelungen in Pfad + Modul (kein `symmetry/symmetry.py` → besser `symmetry/core.py` oder `symmetry/ops.py`).
- Interne Implementierungen: `_helper.py` oder `_function` Präfix.

---
## Logging Guidelines
- Workflows: Nur High-Level Fortschritt.
- Domain/Services: Debug-fähige Detail Logs via Logger `electrofit.<sub>`.
- Konfig-Layering Log-Kategorien: `[config][override]`, `[config][fill]`, `[config][warn]`.

---
## Offene Fragen / ToDo klären
- Benötigen wir einen stabilen Python API-Modus außerhalb CLI? (Batch/Library Use)
- Sollen Konformer-Snapshots langfristig wegfallen (nur Root-Snapshot)?
- Einführen eines Settings-Objekts statt direkte weitergabe von Path/Sting?

---
## Aktueller Status (Fortschritt markieren)
| Phase | Aufgabe | Status |
|-------|---------|--------|
| 1 | snapshot -> infra | done |
| 1 | external -> adapters (Namespace entfernt) | done |
| 1 | logging -> infra | done |
| 1 | scratch -> infra | done |
| 1 | utils_curly_brace -> viz/curly_brace.py | done |
| 1 | plotting/helpers -> viz/helpers Migration | done |
| 1 | Deprecation Shim plotting/__init__ | done |
| 1 | rewrite_imports alte Einträge entfernen | done |
| 2 | domain symmetry | done (write_symmetry + equiv_groups + CLI) |
| 2 | domain charges | partial (process_conformer extrahiert, weitere Zerlegung offen) |
| 2 | domain prep | partial (process_initial extrahiert, RESP-Symmetrie jetzt shared) |
| 2 | step6 aggregation parity + experimental outlier & net charge normalization | done |
| 2 | io modularisierung (ff, mol2_ops, equiv, resp_edit, legacy) | done |
| 3 | pipeline/steps Struktur | in-progress (step0-6 migrated) |
| 3 | api facade | open |
| 3 | step6 AggregationResult + Plot-Service Entkopplung | open |
| 4 | remove legacy aliases | open |
| 4 | hide merge_into_snapshot | partial (intern markiert) |
| 4 | Entferne plotting Shim | in-progress (Dateisystem-Cleanup: leeres 'plotting/' Verzeichnis + Bytecode entfernen) |
| 4 | Entferne alte Import-Rewrites | partial (Hauptmapping bereinigt) |
| 4 | Entferne binary_mapping Re-Exports | open |
| 4 | DeprecationWarnings io/files Übergangsfunktionen | open |

(Beim Fortschritt jeweils "open | in-progress | done" einsetzen.)

---
## Pflege
Dieses Dokument wird bei jeder abgeschlossenen Umsetzungs-Phase aktualisiert. Bitte Pull Requests mit Änderungen an Architektur / Struktur verlinken und Status-Tabelle pflegen.

---
## Ergänzende durchgeführte Arbeiten
- Dependency-Inventar Tool (`tools/dependency_inventory.py`) hinzugefügt: erzeugt `dependency_graph.json` & `dependency_report.md` (AST-basierte Importanalyse) – unterstützt sichere Refactors.
- Timestamps des Inventars auf `datetime.now(timezone.utc)` (TZ-aware) umgestellt.
- Tests nach jeder Migrationsetappe ausgeführt (alle grün; nur externe DeprecationWarnings von SWIG-Typen).
// Removed legacy `core.process_conform` shim (fully deleted; domain API only).

## Nächste kurzfristige Low-Risk Schritte (Update 2025-08-19, nach Migration Step5 Orchestrator)

1. Feingranulare Aufteilung `process_conformer` (Gaussian Phase, RESP Phase) zur Vorbereitung Adapter-Layer.
2. Services-Skelett: `services/config_service.py` & `services/conformer_sampling.py` (Sampling bereits teilweise extrahiert, API-Fassade fehlt).
3. Deprecation-Zeitplan dokumentieren (Entfernung Shims earliest nach 2 Minor Releases) + README Abschnitt „Migration“.
4. Cleanup verwaiste Alt-Workflows (Entfernung step5_* Altcode nach Stabilisierung /2 Releases).
5. Migration Steps 6–8 analog (Durchschnittsladungen, final setup, final sim) mit Domain-Extraktion.
6. Einführung ENV-Schalter für erweiterte Diagnostik (ELECTROFIT_DISABLE_ADV_DIAG) – implementiert für Step5; ggf. globalisieren.
7. Feingranulare Zerlegung `process_conformer` (Stufen: prepare_inputs, generate_esp, prepare_resp, write_symmetry_doc) inkl. Timing-Decorator; Transitional-Helper entfernt.

