# Electrofit Paketdokumentation

Diese Datei bietet einen vollständigen Überblick über Aufbau, Datenflüsse, wichtige Module/Funktionen sowie eine priorisierte Roadmap für die verbleibende REFAKTORIERUNG. Ziel: Vereinheitlichung der neuen Domänenschichten gegenüber Legacy-Workflows, Reduktion von Duplikation und Vorbereitung auf eine stabile öffentliche API.

---
 
## 1. Zweck & Gesamtprozess

Elektrostatische Ladungsableitung (Gaussian/ESP/RESP oder BCC) und GROMACS‑Simulationsvorbereitung in reproduzierbaren, schrittweisen Pipeline‑Stages (Step0–Step8) mit konfigurationsgetriebenen „Snapshots“ (`electrofit.toml`) pro Arbeitsverzeichnis.

Hauptpipeline (vereinheitlicht über `electrofit pipeline steps`):


| Step | Zweck (Kurz) | Ergebnis-Artefakte | Hinweise |
|------|--------------|--------------------|----------|
| 0 | Kopie Eingabe-Moleküle nach `process/` | Arbeitsverzeichnisse | stabil |
| 1 | Initiale Vorbereitung (Gaussian Input, Antechamber/RESP oder BCC) | `run_gau_create_gmx_in/` + Snapshot | stabil |
| 2 | Vorbereitung Simulationslauf (Topologie + MDP) | `run_gmx_simulation/` + `run.json` | stabil |
| 3 | Erzeugung Produktions-Inputs (solvatisiert, Ionisierung etc.) | GROMACS Input-Dateien | stabil |
| 4 | Trajektorien-Sampling → Konformations-Extraktion | `extracted_conforms/` | stabil |
| 5 | Per-Konformer Gaussian/ESP/RESP oder BCC (parallelisiert) | RESP/BCC Resultate, Symmetrie-Dokumente | refaktorisiert |
| 6 | Aggregation & Mittelung der Ensemble-Ladungen | `results/` (avg MOL2, Histogramm-Manifest) | Outlier-Filter EXPERIMENTELL |
| 7 | Finaler Simulations-Setup | `run_final_gmx_simulation/` | geplant |
| 8 | Start finaler Produktion | Laufende Simulation | geplant |


---

## 2. Schichten & Verantwortlichkeiten

| Layer | Rolle | Reinheit / Seiteneffekte | Hauptmodule |
|-------|-------|-------------------------|-------------|
| domain | Fachlogik (Berechnung, Selektion, Aggregation) | möglichst rein (I/O injiziert) | `domain/aggregation`, `domain/charges`, `domain/sampling`, `domain/symmetry`, `domain/final_sim` |
| adapters | Abstraktion externer Tools (Gaussian, RESP, GROMACS) | gezielte Subprocess-/Datei‑Seiteneffekte | `adapters/*.py` |
| infra | Querschnitt: Logging, Snapshots, Scratch Mgmt, Decisions | kontrollierte Seiteneffekte | `infra/*.py` |
| config | Laden / Mergen TOML → Dataclasses | rein | `config/loader.py`, `config/merger.py` |
| io | Datei-Parsing, Formatkonvertierung, Sammelhilfen | viele Seiteneffekte (Zersetzung nötig) | `io/files.py`, `io/mol2.py` |
| viz | Diagnostische Plots + Manifestvalidierung | Seiteneffekt (Dateiannahmen) | `viz/histograms.py`, `viz/hist_manifest.py` |
| pipeline.steps | Dünne Orchestratoren (CLI-Fokus) | Imperativ | `pipeline/steps/step*.py` |
| cli | Legacy / Kommandos / Flags | Imperativ | `cli/*` |
| workflows (removed) | Ehemalige Legacy-Orchestrierung | entfernt | (ImportError) |

Architektur-Ziel: Domain ↔ Adapter Schnitt klar; IO-Helfer modularisieren (kein „God Module“). Pipeline‑Steps bleiben dünn und rufen Domain/APIs auf.

---

## 3. Zentrale Konzepte

1. Snapshot Layering: `compose_snapshot` vereinigt (Projekt Defaults → Molekül Input → Prozess-Kontext → Upstream Snapshot → Optionale Override). Ergebnis ist deterministisch und auditierbar.
2. DecisionModel (`infra.decisions`): Strukturierte, maschinenlesbare Laufentscheidungen (Sampling, Protokollpfad) für Logging & Tests.
3. Histogramm-Manifest (Step6): JSON mit Expectation/Creation/Reason Feldern zur robusten Testbarkeit statt fragiler Platzhalter-PDFs.
4. Scratch Management: Temporäre Rechenverzeichnisse (Gaussian etc.) kontrolliert finalisieren (`scratch_manager`).
5. Symmetriebehandlung: Generierung äquivalenter Atomgruppen, optionales Anpassen von RESP Eingaben, Gruppendurchschnitt in Aggregation.

---

## 4. Modul- & Funktionsindex (Kurzzusammenfassung)

Nur wichtige / öffentliche oder stark genutzte Funktionen; Hilfsfunktionen im Code näher dokumentieren (In-Code Docstrings ausbaubar).

### domain

- `aggregation.average_charges.process_molecule_average_charges` – Kernfluss Step6: Einlesen RESP Resultate, Outlier-IQR Filter, Gruppenmittel, MOL2‑Update, Histogramm/Manifest.
- `charges.process_conformer.process_conformer` – Orchestriert pro Konformer: Vorbereitung → ESP → RESP/BCC → Symmetrie-Dokument.
- `charges.conformer_batch.process_conformer_dir` – Batch-Wrapper + Logging für einen Konformerpfad.
- `sampling.frames.select_frame_indices` – Frame-Selektion (linear, random, maxmin).
- `sampling.conformer_io.prepare_conformer_directory` – Strukturiert Ausgabeverzeichnis & Snapshot.
- `symmetry.equiv_groups.generate_equivalence_groups` – Äquivalenzgruppen aus RESP Eingaben ableiten.
- `symmetry.resp_constraints.apply_and_optionally_modify` – Anwenden/Modifizieren Symmetrie in RESP Inputs.
- `final_sim.prepare_final_sim_directory` / `launch_final_sim_run` – Finaler Simulationssetup & Start.

### adapters

- `gaussian.run_with_result` – Ausführung inkl. strukturierter Rückgabe & optionalem Cache.
- `resp.run_two_stage_with_result` – Zweistufige RESP-Berechnung mit Symmetrieapplikation.
- `gromacs.set_up_production` – Produktions-Workflow (Solvatisierung, Ionisierung, MDP‑Sequenz). (Derzeit noch umfangreicher Imperativ; Kandidat für Feingranularisierung.)

### infra

- `config_snapshot.compose_snapshot` – Layering-Engine.
- `scratch_manager.finalize_scratch_directory` – Sicheres Verschieben/Committen temporärer Ergebnisse.
- `decisions.build_*_decision` – Builder je Pipelinephase für nachvollziehbares Logging.
- `logging.setup_logging` / `log_run_header` – Einheitliche Log-Header.

### io

Modularisiert (Aufspaltung des früheren God-Moduls `files.py`).

- `ff.py` – Forcefield / Topologie Includes & POSRES Handling.
- `mol2_ops.py` – MOL2↔PDB Konvertierungen & Parsing (inkl. Ladungsextraktion).
- `equiv.py` – Erzeugung & Persistenz von Äquivalenzgruppen.
- `resp_edit.py` – Modifikation von RESP Input bzgl. Symmetrie.
- `files.py` – Übergangscontainer: generische File/Path Utilities (Deprecation geplant, keine parsing Logik mehr hinzufügen!).
- `legacy/binary_mapping.py` – Isolierte Legacy-Hilfen (nur temporär re-exportiert).
- `mol2.update_mol2_charges` – Validiertes Umschreiben von MOL2 Partial Charges.

### viz

- `histograms.*` – Erstellung Diagnostik-PDFs (per Atom / Gruppen overlay).
- `hist_manifest.validate_hist_manifest` – Schema/Integritätscheck Manifest.

### pipeline.steps

- Schritt-orientierte schlanke Orchestratoren (`step0`…`step8`). Halten sich an Muster: Verzeichnisse entdecken → Snapshot layern → Domain-Funktion aufrufen → Konsolenstatus.

### cli (aktuell)

CLI bleibt dünne Shell um `pipeline.steps`. Legacy `workflows.*` wurde vollständig entfernt (harte Migration; ImportError mit Hinweis auf neue Pfade).

---

## 5. Aktuelle Schwachstellen & Risiken

| Thema | Beobachtung | Risiko | Priorität |
|-------|-------------|--------|----------|
| God Module `io/files.py` | Vermischt Parsing, Ausgabe, Pfad- & Businesslogik | Erschwerte Testbarkeit & Erweiterung | Hoch |
| Legacy Workflows parallel zu Pipeline | Doppelter Codepfad | Divergierende Semantik / Wartungskosten | Hoch |
| Direkte Nutzung von `cli.run_commands` in Domain/adapters | Vermischt CLI-Hilfen mit Kernlogik | Geringere Portabilität/API Sauberkeit | Mittel |
| Fehlende Typhints/Contracts in einigen Adapterfunktionen | Implizite Annahmen | Spätere Integrationsfehler | Mittel |
| Unvollständige Fehlerklassifizierung (meist bool + msg) | Grobe Fehlerbehandlung | Fehlende granulare Recovery | Mittel |
| Plot- & I/O-Seiteneffekte tief in Aggregationslogik | Geringere Reinheit / Testbarkeit | Schwerer Unit-Test Outlier-Pfad | Niedrig-Mittel |
| GROMACS Adapter monolithisch | Erschwert alternative Backends / Unit Tests | Erweiterungshemmnis | Mittel |

---

## 6. Refactor Roadmap (Phasen)

### Phase 1 (Schnelle Verbesserungen)

1. Modulaufsplittung `io/files.py` (mechanisch, keine Semantikänderung) in klar benannte Submodule. *Akzeptanzkriterium:* Import-Auflösungen & Tests unverändert grün.
2. Typhints & Rückgabedataclasses für Kernfunktionen mit `(ok, msg)` Pattern → Einführung `Result[T]` (dataclass: `ok: bool`, `value: T | None`, `error: str | None`, `reason: Enum`).
3. Adapter-Aufrufe entkoppeln von `cli.run_commands`: extrahiere reine Subprozess-Helfer (`infra.process` oder `adapters.exec`).
4. Ergänze minimalen Architektur-Readme Abschnitt (Kurzlink auf dieses Dokument) im Paket-Root.

### Phase 2 (Strukturelle Konsolidierung)

1. Auslagerung Plot-Erzeugung aus `average_charges.process_molecule_average_charges` in Service (`viz.service.histograms.generate(...)`), Domain-Funktion gibt nur Sammel-Objekte zurück.
2. Zerlegung `adapters.gromacs.set_up_production` in kleinere klare Schritte (Box Anlegung, Solvatation, Ionen, Minimierung, Equilibration, Produktion). Einführung eines Kommandosequenz-Descriptors (Liste von `GmxStage` dataclasses) zur besseren Testbarkeit (Trockenlauf).

### Phase 3 (API Verfestigung & Erweiterbarkeit)

1. Einheitliches Fehler & Logging Model: Exceptions mit Kontext (`ElectrofitError` Hierarchie) + optional strukturiertes JSON Logging.
2. Public API Oberfläche (`electrofit.api`): stabile Funktionen für Skriptintegration (z.B. `run_step5(project: Path, **opts)`).
3. Konfigurationsschema validieren (lightweight Validator) vor Snapshot Persistenz; klare Fehlermeldungen.

### Phase 4 (Qualität & Zukunftssicherheit)

1. Parallelisierung vereinheitlichen: Zentraler Executor/Worker Layer (Progress / Error Aggregation) statt ad-hoc `ProcessPoolExecutor`.
2. Caching-Strategie standardisieren (Gaussian Cache, RESP Intermediates) → einheitliches `cache/` Layout + Prüfsummen.
3. Optional: Plugin-System (Entry Points) für alternative Ladungs-Methoden (AM1-BCC Varianten, XTB, etc.).

---

## 7. Migrations-/Deprecation-Strategie

| Schritt | Maßnahme | Status |
|---------|----------|--------|
| A | Workflows entfernt (ImportError) | erledigt |
| B | Öffentliche API definieren (`electrofit.api`) | offen |
| C | Konsistentes Fehlermodell & Validator | geplant |
| D | SemVer strikt ab API Freeze | geplant |

---

## 8. Teststrategie Ausbau

| Ebene | Ergänzung | Tools |
|-------|-----------|-------|
| Unit | Neue Submodule aus `io` isoliert testen (reines Parsing / keine Filesysteme durch Strings + `tempfile`) | pytest |
| Integration | Szenarien Step4→Step6 (mit/ohne Outlier, Symmetrie an/aus) | vorhandene e2e Tests erweitern |
| Property | Sampling Indizes (Random, MaxMin) – Invarianten: keine Duplikate, Bereichsgrenzen | hypothesis (optional) |
| Performance (Smoke) | Gaussian Cache Treffer / Miss Ratio Logging validate | Timing Harness |
| CLI Snapshot | Golden Files für `electrofit.toml` nach Compose (Hash) | pytest + fixture |

---

## 9. Qualitätsmetriken / Definition of Done

| Bereich | Ziel |
|---------|-----|
| Testabdeckung Domain | >90% Linien / >95% Zweige kritische Logik (Sampling, Aggregation) |
| Zyklomatische Komplexität einzelner Funktionen | <15 (außer orchestrierende CLI) |
| `io/files.py` Größe nach Aufspaltung | <400 Zeilen pro Submodul, thematisch fokussiert |
| Lint (ruff/mypy optional) | Keine neuen Warnungen in PRs |
| Start-zu-Ende Pipeline (Demo Projekt) | < X Minuten (Dokumentierter Referenzlauf) |

---

## 10. Offene Fragen (zur Klärung vor Phase 2)

1. Soll ein generischer Subprozess-Wrapper (Retry, Timeout, Logging) eingeführt werden? (Empfehlung: Ja.)
2. Benötigt die Histogramm-Pipeline zukünftig maschinenlesbare Statistik-Ausgabe (JSON) für CI‑Regeln? (Empfehlung: hinzufügen in Phase 2.)
3. Sollen alternative Ladungsmodelle (z.B. XTB Mulliken) roadmapfähig gemacht werden? (Falls ja → Adapter-Schnitt design früh stabilisieren.)

---

## 11. Kurz-Fazit

Die aktuelle Struktur ist funktional und modularisiert viele ehemals monolithische Bereiche. Aufspaltung von `io/files.py` abgeschlossen; weitere Schritte: Deprecation verbleibender Re-Exports, Entfernung Legacy-Workflows, Feingranularisierung GROMACS Adapter. Histogramm-Manifest ist Muster für reproduzierbare Diagnostik; Outlier-Filter (EXPERIMENTELL) + Net Charge Normalisierung hinzugefügt.

---

## 12. Nächste konkrete Schritte (operativ)

1. (PR 1) Aufspaltung `io/files.py` (mechanisch) + Update Imports + Tests.
2. (PR 2) Adapter-Subprozess-Wrapper + Entfernen direkter `cli.run_commands` Aufrufe aus Domain.
3. (PR 3) `average_charges` entkoppeln: Plot-Service auslagern; Funktion gibt strukturiertes `AggregationResult` zurück.
4. (PR 4) Public API Draft + Stabilitäts-Hinweis (Einführung `electrofit.api`).
5. (PR 5) Fehler-/Logging-Hierarchie + Validator Einbindung.

---

## 13. Experimentelle Features

| Feature | Beschreibung | Aktivierung | Stabilität |
|---------|--------------|-------------|------------|
| Step6 Outlier Removal | IQR-basierter Ladungs-Outlier Filter vor Gruppenmittelung; schreibt bereinigte Artefakte (`cleaned_adjusted_charges.json`). | Konfig `remove_outlier = true` | EXPERIMENTELL (Algorithmus & Dateinamen können sich ändern) |
| Net Charge Normalization | Normalisiert nach Outlier-Entfernung Summenladung auf Sollwert (proportional, Fallback gleichmäßig). | Implizit bei aktivem Outlier-Filter | Stabil (intern optimierbar) |

Für produktive Runs: Outlier-Filter nur nach Sichtung Histogramme aktivieren.

---

*Stand: Automatisch generiert (Audit Session, aktualisiert nach IO-Modularisierung & Step6 Paritäts-Verbesserungen). Bitte bei Strukturänderungen synchron halten.*
